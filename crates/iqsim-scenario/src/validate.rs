//! Semantic checks that the schema alone cannot express.
//!
//! - The format version matches this build.
//! - Entity names (platforms, emitters, receivers, antenna ports, sinks,
//!   resources) are unique across the scenario; mount names within a platform and
//!   channel names within a receiver.
//! - Every reference names an entity of the right type, and a mount reference
//!   names a mount on the referenced platform.
//! - Frequencies that must be positive are. Each channel or antenna-port response has
//!   0 < bandwidth < stopband ≤ complex sample rate (the stopband defaults to the
//!   rate), attenuation 20 to 200 dB and group delay 0 to 1 ms.
//! - Clock and LO `fixed_offset` models have an `offset`, and `ideal` ones do not.
//! - Render duration and block length are positive.
//! - A `world.geoid` path is not empty, and `world.pressure_offset` is within
//!   ±200 hPa.
//! - A receiver has at least one channel.
//! - An emitter has an `iq_file` table exactly when `waveform = "iq_file"`; its start
//!   is not negative, its path is not empty, and its source response has
//!   0 < bandwidth < stopband (≤ the stated sample rate) and 20 to 200 dB. Checks that
//!   need the recording itself happen when the scenario is resolved.
//!
//! - A platform trajectory has the fields its kind needs and no others, each maneuver
//!   has one rate and one end condition, limits are positive, and track points are in
//!   strictly increasing time order (contracts §11.4, §13). Plans that need geometry
//!   (a waypoint inside the turn circle, an unreachable altitude) are checked when the
//!   scenario is resolved.
//!
//! Library references (antennas, waveforms, behaviors) are not resolved yet.

use std::collections::BTreeMap;
use std::ops::Range;

use iqsim_units::{Frequency, SampleRate};
use toml::Spanned;

use crate::diagnostic::Diagnostic;

/// Largest |world.pressure_offset|, Pa (`iqsim_geo::atmosphere::MAX_PRESSURE_OFFSET_PA`).
const MAX_PRESSURE_OFFSET_PA: f64 = 20_000.0;
use crate::ident::Ident;
use crate::model::{
    ClockModel, ClockModelSpec, LoModel, LoModelSpec, ManeuverKind, ManeuverSpec, Platform,
    ResponseDetail, SCENARIO_SCHEMA_VERSION, Scenario, Trajectory, TrajectoryKind,
};

/// Every problem found, in source order.
pub fn validate(scenario: &Scenario) -> Vec<Diagnostic> {
    let mut v = Validator::default();
    v.run(scenario);
    v.diagnostics.sort_by_key(|d| {
        d.span
            .as_ref()
            .map_or((usize::MAX, 0), |s| (s.start, s.end))
    });
    v.diagnostics
}

#[derive(Default)]
struct Validator {
    diagnostics: Vec<Diagnostic>,
}

/// Names declared in one namespace, with their kind and where they were defined.
#[derive(Default)]
struct Names<'a> {
    defined: BTreeMap<&'a str, (&'static str, Range<usize>)>,
}

impl<'a> Names<'a> {
    fn declare(&mut self, v: &mut Validator, name: &'a Ident, kind: &'static str, scope: &str) {
        if let Some((first_kind, first)) = self.defined.get(name.as_str()) {
            v.diagnostics.push(
                Diagnostic::new(
                    name.span(),
                    format!("duplicate name \"{}\"{scope}", name.as_str()),
                )
                .with_note(
                    first.clone(),
                    format!(
                        "\"{}\" is first defined here, as {}",
                        name.as_str(),
                        with_article(first_kind)
                    ),
                ),
            );
        } else {
            self.defined.insert(name.as_str(), (kind, name.span()));
        }
    }

    fn kind_of(&self, name: &str) -> Option<&'static str> {
        self.defined.get(name).map(|(kind, _)| *kind)
    }
}

/// "a platform", "an emitter".
fn with_article(kind: &str) -> String {
    let article = if kind.starts_with(['a', 'e', 'i', 'o', 'u']) {
        "an"
    } else {
        "a"
    };
    format!("{article} {kind}")
}

impl Validator {
    fn error(&mut self, span: Range<usize>, message: String) {
        self.diagnostics.push(Diagnostic::new(span, message));
    }

    /// Checks that `reference` names an entity of `kind`.
    fn reference(&mut self, names: &Names<'_>, reference: &Ident, kind: &'static str, user: &str) {
        match names.kind_of(reference.as_str()) {
            Some(found) if found == kind => {}
            Some(found) => self.error(
                reference.span(),
                format!(
                    "{user} refers to \"{}\", which is {}, not {}",
                    reference.as_str(),
                    with_article(found),
                    with_article(kind)
                ),
            ),
            None => self.error(
                reference.span(),
                format!("{user} refers to unknown {kind} \"{}\"", reference.as_str()),
            ),
        }
    }

    fn positive(&mut self, value: &Spanned<Frequency>, what: &str) {
        if !value.get_ref().is_positive() {
            self.error(value.span(), format!("{what} must be positive"));
        }
    }

    /// 0 < bandwidth < stopband ≤ complex sample rate, compared exactly
    /// (`width_µHz · den` against `10⁶ · num`), plus the attenuation and delay ranges.
    fn response(
        &mut self,
        bandwidth: &Spanned<Frequency>,
        rate: &Spanned<SampleRate>,
        detail: &ResponseDetail,
        owner: &str,
    ) {
        self.positive(bandwidth, "the bandwidth");
        let rate = rate.get_ref().get();
        let scaled = |f: Frequency| i128::from(f.as_micro_hz()) * i128::from(rate.den());
        let rate_uhz = 1_000_000 * i128::from(rate.num());
        let shown = iqsim_units::SampleRate::new(rate);
        if scaled(*bandwidth.get_ref()) > rate_uhz {
            self.error(
                bandwidth.span(),
                format!(
                    "{owner}: bandwidth {} exceeds the complex sample rate {shown}",
                    bandwidth.get_ref()
                ),
            );
            return;
        }
        match &detail.stopband {
            Some(stop) => {
                if scaled(*stop.get_ref()) > rate_uhz {
                    self.error(
                        stop.span(),
                        format!(
                            "{owner}: stopband {} exceeds the complex sample rate {shown}",
                            stop.get_ref()
                        ),
                    );
                } else if stop.get_ref() <= bandwidth.get_ref() {
                    self.error(
                        stop.span(),
                        format!(
                            "{owner}: stopband {} must be wider than the bandwidth {}",
                            stop.get_ref(),
                            bandwidth.get_ref()
                        ),
                    );
                }
            }
            None => {
                if scaled(*bandwidth.get_ref()) == rate_uhz {
                    self.error(
                        bandwidth.span(),
                        format!(
                            "{owner}: bandwidth {} leaves no transition band below the \
                             default stopband (the sample rate {shown}); narrow it or set \
                             response.stopband",
                            bandwidth.get_ref()
                        ),
                    );
                }
            }
        }
        if let Some(a) = &detail.attenuation
            && !(20.0..=200.0).contains(&a.get_ref().db())
        {
            self.error(
                a.span(),
                format!("{owner}: stopband attenuation must be 20 to 200 dB"),
            );
        }
        if let Some(d) = &detail.group_delay {
            let ps = d.get_ref().get().as_ps();
            if !(0..=1_000_000_000).contains(&ps) {
                self.error(d.span(), format!("{owner}: group delay must be 0 to 1 ms"));
            }
        }
    }

    fn iq_file(&mut self, e: &crate::model::Emitter, user: &str) {
        let wants = e.waveform.as_str() == "iq_file";
        let Some(q) = &e.iq_file else {
            if wants {
                self.error(
                    e.waveform.span(),
                    format!("{user}: waveform \"iq_file\" needs an iq_file table"),
                );
            }
            return;
        };
        if !wants {
            self.error(
                q.path.span(),
                format!(
                    "{user}: an iq_file table needs waveform = \"iq_file\" (it is \"{}\")",
                    e.waveform.as_str()
                ),
            );
        }
        if q.path.get_ref().trim().is_empty() {
            self.error(q.path.span(), format!("{user}: iq_file.path is empty"));
        }
        if let Some(st) = &q.start
            && st.get_ref().get().as_ps() < 0
        {
            self.error(
                st.span(),
                format!("{user}: iq_file.start cannot be negative"),
            );
        }
        let r = &q.response;
        let width_ok = |v: &mut Self, f: &Spanned<Frequency>, what: &str| {
            v.positive(f, what);
            if let Some(rate) = &q.sample_rate {
                let rate = rate.get_ref().get();
                let scaled = i128::from(f.get_ref().as_micro_hz()) * i128::from(rate.den());
                if scaled > 1_000_000 * i128::from(rate.num()) {
                    v.error(
                        f.span(),
                        format!("{user}: iq_file {what} exceeds the native sample rate"),
                    );
                }
            }
        };
        if let Some(b) = &r.bandwidth {
            width_ok(self, b, "response bandwidth");
        }
        if let Some(st) = &r.stopband {
            width_ok(self, st, "response stopband");
        }
        if let (Some(b), Some(st)) = (&r.bandwidth, &r.stopband)
            && st.get_ref() <= b.get_ref()
        {
            self.error(
                st.span(),
                format!("{user}: iq_file response stopband must be wider than its bandwidth"),
            );
        }
        if let Some(a) = &r.attenuation
            && !(20.0..=200.0).contains(&a.get_ref().db())
        {
            self.error(
                a.span(),
                format!("{user}: iq_file response attenuation must be 20 to 200 dB"),
            );
        }
    }

    fn clock_and_lo(
        &mut self,
        clock: &ClockModelSpec,
        lo: &LoModelSpec,
        at: Range<usize>,
        owner: &str,
    ) {
        match (clock.model, &clock.offset) {
            (ClockModel::Ideal, Some(o)) => self.error(
                o.span(),
                format!("{owner}: an ideal clock has no offset; use model = \"fixed_offset\""),
            ),
            (ClockModel::FixedOffset, None) => self.error(
                at.clone(),
                format!("{owner}: clock model \"fixed_offset\" needs an offset"),
            ),
            _ => {}
        }
        match (lo.model, &lo.offset) {
            (LoModel::Ideal, Some(o)) => self.error(
                o.span(),
                format!("{owner}: an ideal LO has no offset; use model = \"fixed_offset\""),
            ),
            (LoModel::FixedOffset, None) => self.error(
                at,
                format!("{owner}: LO model \"fixed_offset\" needs an offset"),
            ),
            _ => {}
        }
    }

    fn trajectory(&mut self, p: &Platform) {
        let user = format!("platform \"{}\"", p.name.as_str());
        let t: &Trajectory = &p.trajectory;
        let kind = t.kind;
        let kind_name = match kind {
            TrajectoryKind::Static => "static",
            TrajectoryKind::Maneuvers => "maneuvers",
            TrajectoryKind::Route => "route",
            TrajectoryKind::Track => "track",
            TrajectoryKind::Ground => "ground",
        };
        let name_span = p.name.span();
        let track = kind == TrajectoryKind::Track;
        match (&p.position, track) {
            (None, false) => self.error(name_span.clone(), format!("{user} needs a position")),
            (Some(_), true) => self.error(
                name_span.clone(),
                format!("{user}: a track trajectory starts at its first point; remove `position`"),
            ),
            _ => {}
        }
        let airborne = matches!(
            kind,
            TrajectoryKind::Maneuvers | TrajectoryKind::Route | TrajectoryKind::Track
        );
        let moving = kind != TrajectoryKind::Static;
        let flies = matches!(
            kind,
            TrajectoryKind::Maneuvers | TrajectoryKind::Route | TrajectoryKind::Ground
        );
        let not_for = |v: &mut Self, span: Option<Range<usize>>, field: &str, ok: bool| {
            if let Some(span) = span
                && !ok
            {
                v.error(
                    span,
                    format!("{user}: `{field}` does not apply to a {kind_name} trajectory"),
                );
            }
        };
        not_for(self, t.start.as_ref().map(Spanned::span), "start", flies);
        not_for(
            self,
            t.heading.as_ref().map(Spanned::span),
            "heading",
            flies,
        );
        not_for(self, t.speed.as_ref().map(Spanned::span), "speed", flies);
        not_for(
            self,
            t.attitude.as_ref().map(Spanned::span),
            "attitude",
            !moving,
        );
        not_for(
            self,
            t.angle_of_attack.as_ref().map(Spanned::span),
            "angle_of_attack",
            airborne,
        );
        not_for(self, t.wind.as_ref().map(Spanned::span), "wind", airborne);
        not_for(
            self,
            t.smoothing.as_ref().map(Spanned::span),
            "smoothing",
            track,
        );
        not_for(
            self,
            t.maneuvers.first().map(Spanned::span),
            "maneuvers",
            matches!(kind, TrajectoryKind::Maneuvers | TrajectoryKind::Ground),
        );
        not_for(
            self,
            t.waypoints.first().map(Spanned::span),
            "waypoints",
            matches!(kind, TrajectoryKind::Route | TrajectoryKind::Ground),
        );
        not_for(self, t.points.first().map(Spanned::span), "points", track);
        let limits = &t.limits;
        let limit_spans = [
            limits.bank.as_ref().map(Spanned::span),
            limits.roll_rate.as_ref().map(Spanned::span),
            limits.accel.as_ref().map(Spanned::span),
            limits.jerk.as_ref().map(Spanned::span),
            limits.vertical_accel.as_ref().map(Spanned::span),
            limits.vertical_jerk.as_ref().map(Spanned::span),
            limits.climb_rate.as_ref().map(Spanned::span),
        ];
        if !moving {
            for span in limit_spans.iter().flatten() {
                self.error(
                    span.clone(),
                    format!("{user}: a static trajectory has no limits"),
                );
            }
        }

        let need = |v: &mut Self, present: bool, field: &str| {
            if !present {
                v.error(
                    name_span.clone(),
                    format!("{user}: a {kind_name} trajectory needs `{field}`"),
                );
            }
        };
        match kind {
            TrajectoryKind::Static | TrajectoryKind::Track => {}
            TrajectoryKind::Maneuvers => {
                need(self, t.heading.is_some(), "heading");
                need(self, t.speed.is_some(), "speed");
                need(self, !t.maneuvers.is_empty(), "maneuvers");
            }
            TrajectoryKind::Route => {
                need(self, t.speed.is_some(), "speed");
                need(self, !t.waypoints.is_empty(), "waypoints");
            }
            TrajectoryKind::Ground => {
                need(
                    self,
                    t.heading.is_some() || !t.waypoints.is_empty(),
                    "heading",
                );
                need(self, t.speed.is_some(), "speed");
                if t.maneuvers.is_empty() == t.waypoints.is_empty() {
                    self.error(
                        name_span.clone(),
                        format!(
                            "{user}: a ground trajectory needs either `maneuvers` or `waypoints`"
                        ),
                    );
                }
                if let Some(pos) = &p.position
                    && pos.altitude().reference() != iqsim_units::AltitudeRef::Agl
                {
                    self.error(
                        name_span.clone(),
                        format!("{user}: a ground platform's altitude must be AGL (its height above the ground)"),
                    );
                }
                for w in &t.waypoints {
                    if w.get_ref().alt.get_ref().reference() != iqsim_units::AltitudeRef::Agl {
                        self.error(
                            w.get_ref().alt.span(),
                            format!("{user}: ground waypoint altitudes must be AGL"),
                        );
                    }
                }
            }
        }
        if let Some(s) = &t.speed
            && s.get_ref().meters_per_second() < 1.0
        {
            self.error(
                s.span(),
                format!("{user}: the speed must be at least 1 m/s"),
            );
        }
        if let Some(a) = &t.angle_of_attack
            && a.get_ref().degrees().abs() >= 45.0
        {
            self.error(
                a.span(),
                format!("{user}: the angle of attack must be within ±45°"),
            );
        }
        if let Some(w) = &t.wind
            && w.get_ref().speed.meters_per_second() < 0.0
        {
            self.error(
                w.span(),
                format!("{user}: the wind speed cannot be negative"),
            );
        }
        if let Some(s) = &t.smoothing
            && s.get_ref().get().as_ps() < 0
        {
            self.error(
                s.span(),
                format!("{user}: the smoothing time constant cannot be negative"),
            );
        }
        let positive = |v: &mut Self, x: Option<(f64, Range<usize>)>, what: &str| {
            if let Some((x, span)) = x
                && (x.is_nan() || x <= 0.0)
            {
                v.error(span, format!("{user}: {what} must be positive"));
            }
        };
        if let Some(b) = &limits.bank {
            let d = b.get_ref().degrees();
            if !(d > 0.0 && d < 80.0) {
                self.error(
                    b.span(),
                    format!("{user}: limits.bank must be between 0° and 80°"),
                );
            }
        }
        positive(
            self,
            limits
                .roll_rate
                .as_ref()
                .map(|x| (x.get_ref().degrees_per_second(), x.span())),
            "limits.roll_rate",
        );
        positive(
            self,
            limits
                .accel
                .as_ref()
                .map(|x| (x.get_ref().meters_per_second_squared(), x.span())),
            "limits.accel",
        );
        positive(
            self,
            limits
                .jerk
                .as_ref()
                .map(|x| (x.get_ref().meters_per_second_cubed(), x.span())),
            "limits.jerk",
        );
        positive(
            self,
            limits
                .vertical_accel
                .as_ref()
                .map(|x| (x.get_ref().meters_per_second_squared(), x.span())),
            "limits.vertical_accel",
        );
        positive(
            self,
            limits
                .vertical_jerk
                .as_ref()
                .map(|x| (x.get_ref().meters_per_second_cubed(), x.span())),
            "limits.vertical_jerk",
        );
        positive(
            self,
            limits
                .climb_rate
                .as_ref()
                .map(|x| (x.get_ref().meters_per_second(), x.span())),
            "limits.climb_rate",
        );

        for (i, m) in t.maneuvers.iter().enumerate() {
            let owner = format!("{user}, maneuver {}", i + 1);
            self.maneuver(m, &owner, kind == TrajectoryKind::Ground);
        }
        for w in &t.waypoints {
            if let Some(s) = &w.get_ref().speed
                && s.get_ref().meters_per_second() < 1.0
            {
                self.error(
                    s.span(),
                    format!("{user}: a waypoint speed must be at least 1 m/s"),
                );
            }
        }
        if track && t.points.len() < 2 {
            self.error(
                name_span,
                format!("{user}: a track needs at least two points"),
            );
        }
        for w in t.points.windows(2) {
            let (a, b) = (&w[0].get_ref().time, &w[1].get_ref().time);
            if b.get_ref().get() <= a.get_ref().get() {
                self.error(b.span(), format!("{user}: track point times must increase"));
            }
        }
    }

    fn maneuver(&mut self, m: &Spanned<ManeuverSpec>, owner: &str, ground: bool) {
        let span = m.span();
        let m = m.get_ref();
        let fields: [(&str, Option<Range<usize>>); 9] = [
            ("duration", m.duration.as_ref().map(Spanned::span)),
            ("rate", m.rate.as_ref().map(Spanned::span)),
            ("bank", m.bank.as_ref().map(Spanned::span)),
            ("heading", m.heading.as_ref().map(Spanned::span)),
            ("by", m.by.as_ref().map(Spanned::span)),
            ("vertical_rate", m.vertical_rate.as_ref().map(Spanned::span)),
            ("angle", m.angle.as_ref().map(Spanned::span)),
            ("altitude", m.altitude.as_ref().map(Spanned::span)),
            ("speed", m.speed.as_ref().map(Spanned::span)),
        ];
        let (name, allowed): (&str, &[&str]) = match m.kind {
            ManeuverKind::Straight => ("straight", &["duration"]),
            ManeuverKind::Turn => ("turn", &["duration", "rate", "bank", "heading", "by"]),
            ManeuverKind::Climb => ("climb", &["duration", "vertical_rate", "angle", "altitude"]),
            ManeuverKind::Speed => ("speed", &["speed"]),
        };
        for (field, s) in &fields {
            if let Some(s) = s
                && !allowed.contains(field)
            {
                self.error(
                    s.clone(),
                    format!("{owner}: `{field}` does not apply to a {name} maneuver"),
                );
            }
        }
        if m.direction.is_some() && m.kind != ManeuverKind::Turn {
            self.error(
                span.clone(),
                format!("{owner}: `direction` only applies to a turn"),
            );
        }
        let count = |names: &[&str]| {
            fields
                .iter()
                .filter(|(f, s)| s.is_some() && names.contains(f))
                .count()
        };
        let exactly_one = |v: &mut Self, names: &[&str]| {
            if count(names) != 1 {
                v.error(
                    span.clone(),
                    format!(
                        "{owner}: a {name} maneuver needs exactly one of {}",
                        names
                            .iter()
                            .map(|n| format!("`{n}`"))
                            .collect::<Vec<_>>()
                            .join(", ")
                    ),
                );
            }
        };
        let positive_duration = |v: &mut Self| {
            if let Some(d) = &m.duration
                && d.get_ref().get().as_ps() <= 0
            {
                v.error(d.span(), format!("{owner}: the duration must be positive"));
            }
        };
        match m.kind {
            ManeuverKind::Straight => {
                exactly_one(self, &["duration"]);
                positive_duration(self);
            }
            ManeuverKind::Turn => {
                exactly_one(self, &["rate", "bank"]);
                exactly_one(self, &["duration", "heading", "by"]);
                positive_duration(self);
                if (m.duration.is_some() || m.by.is_some()) && m.direction.is_none() {
                    self.error(
                        span.clone(),
                        format!("{owner}: a turn by duration or by angle needs a `direction`"),
                    );
                }
                if let Some(r) = &m.rate
                    && r.get_ref().degrees_per_second() <= 0.0
                {
                    self.error(r.span(), format!("{owner}: the turn rate must be positive"));
                }
                if let Some(b) = &m.bank {
                    let d = b.get_ref().degrees();
                    if !(d > 0.0 && d < 80.0) {
                        self.error(
                            b.span(),
                            format!("{owner}: the bank must be between 0° and 80°"),
                        );
                    }
                }
                if let Some(b) = &m.by
                    && b.get_ref().degrees() <= 0.0
                {
                    self.error(
                        b.span(),
                        format!("{owner}: `by` must be positive (use `direction`)"),
                    );
                }
            }
            ManeuverKind::Climb => {
                if ground {
                    self.error(
                        span.clone(),
                        format!("{owner}: a ground platform cannot climb"),
                    );
                }
                exactly_one(self, &["vertical_rate", "angle"]);
                exactly_one(self, &["altitude", "duration"]);
                positive_duration(self);
                if let Some(r) = &m.vertical_rate
                    && r.get_ref().meters_per_second() == 0.0
                {
                    self.error(
                        r.span(),
                        format!("{owner}: the vertical rate cannot be zero"),
                    );
                }
                if let Some(a) = &m.angle {
                    let d = a.get_ref().degrees().abs();
                    if !(d > 0.0 && d < 90.0) {
                        self.error(
                            a.span(),
                            format!("{owner}: the climb angle must be within ±90° and not zero"),
                        );
                    }
                }
            }
            ManeuverKind::Speed => {
                exactly_one(self, &["speed"]);
                if let Some(s) = &m.speed
                    && s.get_ref().meters_per_second() < 1.0
                {
                    self.error(
                        s.span(),
                        format!("{owner}: the speed must be at least 1 m/s"),
                    );
                }
            }
        }
    }

    fn run(&mut self, s: &Scenario) {
        if *s.schema_version.get_ref() != SCENARIO_SCHEMA_VERSION {
            self.error(
                s.schema_version.span(),
                format!(
                    "schema_version {} is not supported; this build reads version {SCENARIO_SCHEMA_VERSION}",
                    s.schema_version.get_ref()
                ),
            );
        }
        if s.metadata.name.get_ref().trim().is_empty() {
            self.error(
                s.metadata.name.span(),
                "metadata.name cannot be empty".to_owned(),
            );
        }

        if let Some(g) = &s.world.geoid
            && g.path.get_ref().trim().is_empty()
        {
            self.error(g.path.span(), "world.geoid.path is empty".to_owned());
        }
        if let Some(p) = &s.world.pressure_offset
            && p.get_ref().pascals().abs() > MAX_PRESSURE_OFFSET_PA
        {
            self.error(
                p.span(),
                "world.pressure_offset must be within ±200 hPa".to_owned(),
            );
        }

        if let Some(render) = &s.render {
            if render.duration.get_ref().get().as_ps() <= 0 {
                self.error(
                    render.duration.span(),
                    "render.duration must be positive".to_owned(),
                );
            }
            if let Some(block) = &render.block
                && block.get_ref().get().as_ps() <= 0
            {
                self.error(block.span(), "render.block must be positive".to_owned());
            }
        }

        let mut names = Names::default();
        for p in &s.platforms {
            names.declare(self, &p.name, "platform", "");
        }
        for e in &s.emitters {
            names.declare(self, &e.name, "emitter", "");
        }
        for r in &s.receivers {
            names.declare(self, &r.name, "receiver", "");
        }
        for a in &s.antenna_ports {
            names.declare(self, &a.name, "antenna port", "");
        }
        for k in &s.sinks {
            names.declare(self, &k.name, "sink", "");
        }
        for r in &s.resources {
            names.declare(self, &r.name, "resource", "");
        }

        for p in &s.platforms {
            self.trajectory(p);
        }

        let mut mounts: BTreeMap<&str, Names<'_>> = BTreeMap::new();
        for p in &s.platforms {
            let scope = format!(" among the mounts of platform \"{}\"", p.name.as_str());
            let mut local = Names::default();
            for m in &p.mounts {
                local.declare(self, &m.name, "mount", &scope);
            }
            // A duplicate platform keeps the first platform's mounts.
            mounts.entry(p.name.as_str()).or_insert(local);
        }

        let attachment = |v: &mut Self, platform: &Ident, mount: Option<&Ident>, user: &str| {
            v.reference(&names, platform, "platform", user);
            if let (Some(mount), Some(platform_mounts)) = (mount, mounts.get(platform.as_str()))
                && platform_mounts.kind_of(mount.as_str()).is_none()
            {
                v.error(
                    mount.span(),
                    format!(
                        "{user} refers to mount \"{}\", which platform \"{}\" does not have",
                        mount.as_str(),
                        platform.as_str()
                    ),
                );
            }
        };

        for e in &s.emitters {
            let user = format!("emitter \"{}\"", e.name.as_str());
            attachment(self, &e.platform, e.mount.as_ref(), &user);
            self.positive(&e.frequency, "the emitter frequency");
            self.iq_file(e, &user);
        }
        for r in &s.receivers {
            let user = format!("receiver \"{}\"", r.name.as_str());
            attachment(self, &r.platform, r.mount.as_ref(), &user);
            if r.channels.is_empty() {
                self.error(
                    r.name.span(),
                    format!("{user} has no channels; add at least one [[receivers.channels]]"),
                );
            }
            let scope = format!(" among the channels of receiver \"{}\"", r.name.as_str());
            let mut channels = Names::default();
            for c in &r.channels {
                channels.declare(self, &c.name, "channel", &scope);
                let owner = format!("channel \"{}\" of {user}", c.name.as_str());
                self.positive(&c.center_frequency, "the centre frequency");
                self.response(&c.bandwidth, &c.sample_rate, &c.response, &owner);
                self.clock_and_lo(&c.clock, &c.lo, c.name.span(), &owner);
                for sink in &c.sinks {
                    self.reference(&names, sink, "sink", &owner);
                }
            }
        }
        for a in &s.antenna_ports {
            let user = format!("antenna port \"{}\"", a.name.as_str());
            attachment(self, &a.platform, a.mount.as_ref(), &user);
            self.positive(&a.center_frequency, "the centre frequency");
            self.response(&a.bandwidth, &a.sample_rate, &a.response, &user);
            for sink in &a.sinks {
                self.reference(&names, sink, "sink", &user);
            }
        }
        for k in &s.sinks {
            let user = format!("sink \"{}\"", k.name.as_str());
            for resource in &k.resources {
                self.reference(&names, resource, "resource", &user);
            }
        }
    }
}
