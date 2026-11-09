//! Semantic checks that the schema alone cannot express.
//!
//! - The format version matches this build.
//! - Entity names (platforms, emitters, receivers, antenna ports, sinks,
//!   resources) are unique across the scenario; mount names within a platform and
//!   channel names within a receiver.
//! - Every reference names an entity of the right type, and a mount reference
//!   names a mount on the referenced platform.
//! - Frequencies that must be positive are, and each channel or antenna-port
//!   bandwidth is at most its complex sample rate.
//! - A receiver has at least one channel.
//!
//! Library references (antennas, waveforms, behaviors) are not resolved yet.

use std::collections::BTreeMap;
use std::ops::Range;

use iqsim_units::{Frequency, SampleRate};
use toml::Spanned;

use crate::diagnostic::Diagnostic;
use crate::ident::Ident;
use crate::model::{SCENARIO_SCHEMA_VERSION, Scenario};

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

    /// Bandwidth must be positive and at most the complex sample rate.
    fn bandwidth(
        &mut self,
        bandwidth: &Spanned<Frequency>,
        rate: &Spanned<SampleRate>,
        owner: &str,
    ) {
        self.positive(bandwidth, "the bandwidth");
        let rate = rate.get_ref().get();
        // bandwidth_µHz ≤ 10⁶ · num/den  ⇔  bandwidth_µHz · den ≤ 10⁶ · num (exact).
        let lhs = i128::from(bandwidth.get_ref().as_micro_hz()) * i128::from(rate.den());
        let rhs = 1_000_000 * i128::from(rate.num());
        if lhs > rhs {
            self.error(
                bandwidth.span(),
                format!(
                    "{owner}: bandwidth {} exceeds the complex sample rate {}",
                    bandwidth.get_ref(),
                    iqsim_units::SampleRate::new(rate)
                ),
            );
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
                self.bandwidth(&c.bandwidth, &c.sample_rate, &owner);
                for sink in &c.sinks {
                    self.reference(&names, sink, "sink", &owner);
                }
            }
        }
        for a in &s.antenna_ports {
            let user = format!("antenna port \"{}\"", a.name.as_str());
            attachment(self, &a.platform, a.mount.as_ref(), &user);
            self.positive(&a.center_frequency, "the centre frequency");
            self.bandwidth(&a.bandwidth, &a.sample_rate, &user);
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
