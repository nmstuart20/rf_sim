//! Scalar motion profiles: turn rate, speed and vertical rate (contracts §13.4).
//!
//! A [`Channel`] is a value `x(s)` over segment-local seconds `s`, made of
//! [`Phase`]s of constant second derivative: `x = x₀ + ẋ₀τ + ½Jτ²`. Its integral is
//! exact too (`X(s) = ∫₀ˢ x`), which gives the turn angle ψ from the turn rate ω and
//! the height change from the vertical rate ḣ. After the last phase the value holds.
//!
//! - [`scurve`]: from `(x₀, ẋ₀)` to `x₁` with `ẋ = 0` at the end, `|ẋ| ≤ A`,
//!   `|ẍ| = J` during ramps (an S-curve; ẋ is continuous). Speed and vertical rate use
//!   it, so acceleration is continuous and jerk bounded.
//! - [`ramp`]: x changes linearly at rate R (ẋ steps). The turn rate uses it, because
//!   only ω itself must be continuous for C2 position.
//! - `move_by` (crate-internal): rise to a peak, hold, return to 0, with a prescribed integral (turn
//!   by an angle, climb by a height). The peak is found by bisection (fixed iteration
//!   count) and the hold time in closed form.

use crate::GeoError;

/// One phase: constant second derivative `jerk` for `duration` seconds.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Phase {
    /// Local start time, seconds.
    pub start: f64,
    /// Duration, seconds (≥ 0).
    pub duration: f64,
    /// Value at the start.
    pub x0: f64,
    /// First derivative at the start.
    pub d0: f64,
    /// Constant second derivative.
    pub jerk: f64,
    /// Integral of the channel from local time 0 to the phase start.
    pub integral0: f64,
}

impl Phase {
    fn at(&self, tau: f64) -> (f64, f64, f64, f64) {
        let x = self.x0 + tau * (self.d0 + tau * (0.5 * self.jerk));
        let d = self.d0 + tau * self.jerk;
        let i = self.integral0 + tau * (self.x0 + tau * (0.5 * self.d0 + tau * (self.jerk / 6.0)));
        (x, d, self.jerk, i)
    }

    fn end(&self) -> (f64, f64, f64) {
        let (x, d, _, i) = self.at(self.duration);
        (x, d, i)
    }
}

/// A value over local time with its derivatives and integral.
#[derive(Clone, Debug, PartialEq)]
pub struct Channel {
    initial: f64,
    phases: Vec<Phase>,
}

/// A planned phase before placement: (duration, second derivative).
pub type Step = (f64, f64);

impl Channel {
    /// A constant channel.
    pub fn constant(x: f64) -> Self {
        Self {
            initial: x,
            phases: Vec::new(),
        }
    }

    /// Steps placed one after another from `(x0, d0)` at local time 0. Steps of zero
    /// duration are dropped.
    pub fn from_steps(x0: f64, d0: f64, steps: &[Step]) -> Self {
        let mut phases: Vec<Phase> = Vec::new();
        let (mut start, mut x, mut d, mut integral) = (0.0, x0, d0, 0.0);
        for &(duration, jerk) in steps {
            if duration <= 0.0 {
                continue;
            }
            let p = Phase {
                start,
                duration,
                x0: x,
                d0: d,
                jerk,
                integral0: integral,
            };
            (x, d, integral) = p.end();
            start += duration;
            phases.push(p);
        }
        Self {
            initial: x0,
            phases,
        }
    }

    /// Local time at which the last phase ends (0 for a constant).
    pub fn duration(&self) -> f64 {
        self.phases.last().map_or(0.0, |p| p.start + p.duration)
    }

    /// Phase boundaries strictly inside (0, duration], local seconds.
    pub fn breakpoints(&self) -> impl Iterator<Item = f64> + '_ {
        self.phases.iter().map(|p| p.start + p.duration)
    }

    /// `(x, ẋ, ẍ, ∫x)` at local time `s ≥ 0`. After the last phase, x holds with zero
    /// derivatives (and the integral grows linearly).
    pub fn at(&self, s: f64) -> (f64, f64, f64, f64) {
        let i = self.phases.partition_point(|p| p.start <= s);
        if i == 0 {
            return (self.initial, 0.0, 0.0, self.initial * s);
        }
        let p = &self.phases[i - 1];
        let tau = s - p.start;
        if tau <= p.duration {
            return p.at(tau);
        }
        let (x, _, integral) = p.end();
        (x, 0.0, 0.0, integral + x * (tau - p.duration))
    }

    /// The largest |x| over the channel (x is linear or quadratic in each phase, so
    /// its extremes are at phase ends or at a vertex).
    pub fn max_abs(&self) -> f64 {
        let mut m = self.initial.abs();
        for p in &self.phases {
            let (x_end, _, _) = p.end();
            m = m.max(p.x0.abs()).max(x_end.abs());
            if p.jerk != 0.0 {
                let tau = -p.d0 / p.jerk;
                if tau > 0.0 && tau < p.duration {
                    m = m.max(p.at(tau).0.abs());
                }
            }
        }
        m
    }

    /// The value after the last phase.
    pub fn final_value(&self) -> f64 {
        self.phases.last().map_or(self.initial, |p| p.end().0)
    }

    /// The integral from 0 to the end of the last phase.
    pub fn total_integral(&self) -> f64 {
        self.phases.last().map_or(0.0, |p| p.end().2)
    }
}

fn check_limits(what: &str, limits: &[f64]) -> Result<(), GeoError> {
    if limits.iter().all(|l| l.is_finite() && *l > 0.0) {
        Ok(())
    } else {
        Err(GeoError::Trajectory(format!(
            "{what}: rate, acceleration and jerk limits must be positive"
        )))
    }
}

/// Steps of an S-curve from `(x0, d0)` to `x1` with zero end derivative, derivative
/// bounded by `max_d` (raised to |d0| if that is larger) and second derivative ±`jerk`.
///
/// The stopping value `x₀ + d₀|d₀|/(2J)` decides the direction s; then the derivative
/// ramps to `s·p`, holds, and ramps to 0. With `b₀ = s·d₀` and `D = s(x₁ − x₀)`, a
/// triangular profile has `p² = J·D + b₀²/2`; if `p > A` it is clamped to A and the
/// hold lasts `(D − (2A² − b₀²)/(2J))/A`.
pub fn scurve(x0: f64, d0: f64, x1: f64, max_d: f64, jerk: f64) -> Result<Vec<Step>, GeoError> {
    check_limits("S-curve", &[max_d, jerk])?;
    if !(x0.is_finite() && d0.is_finite() && x1.is_finite()) {
        return Err(GeoError::NonFinite("profile endpoint"));
    }
    let max_d = max_d.max(d0.abs());
    let stop = x0 + d0 * d0.abs() / (2.0 * jerk);
    let gap = x1 - stop;
    let scale = x0.abs().max(x1.abs()).max(1.0);
    if gap.abs() <= 1.0e-15 * scale {
        // Only bring the derivative to zero.
        return Ok(vec![(d0.abs() / jerk, -d0.signum() * jerk)]);
    }
    let s = gap.signum();
    let b0 = s * d0;
    let dist = s * (x1 - x0);
    let tri = jerk * dist + 0.5 * (b0 * b0);
    let (peak, hold) = if tri <= max_d * max_d {
        (libm::sqrt(tri.max(0.0)), 0.0)
    } else {
        let hold = (dist - (2.0 * (max_d * max_d) - b0 * b0) / (2.0 * jerk)) / max_d;
        (max_d, hold.max(0.0))
    };
    Ok(vec![
        ((peak - b0).max(0.0) / jerk, s * jerk),
        (hold, 0.0),
        (peak / jerk, -s * jerk),
    ])
}

/// One linear ramp of x from `x0` to `x1` at `rate` (a step of the first derivative).
pub fn ramp(x0: f64, x1: f64, rate: f64) -> Result<Channel, GeoError> {
    check_limits("ramp", &[rate])?;
    let d = x1 - x0;
    Ok(Channel::from_steps(x0, d.signum() * rate, &[(d.abs() / rate, 0.0)]).with_ramp_end(x1))
}

impl Channel {
    /// Pins the final value of a ramp channel exactly (the last phase's slope is
    /// re-derived so it ends at `x1`).
    fn with_ramp_end(mut self, x1: f64) -> Self {
        if let Some(p) = self.phases.last_mut()
            && p.duration > 0.0
        {
            p.d0 = (x1 - p.x0) / p.duration;
        }
        self
    }
}

/// The shape of a rise-hold-return profile, as steps for a given signed peak.
pub(crate) trait Shape {
    /// Steps from the start state to `peak` (derivative zero at the end for S-curves).
    fn rise(&self, peak: f64) -> Result<Channel, GeoError>;
    /// Steps from `peak` back to 0, starting where `rise` ended.
    fn fall(&self, peak: f64) -> Result<Channel, GeoError>;
}

/// Linear ramps at `rate` (the turn rate ω, whose integral is the turn angle).
pub(crate) struct RampShape {
    /// Starting value.
    pub x0: f64,
    /// Ramp rate.
    pub rate: f64,
}

impl Shape for RampShape {
    fn rise(&self, peak: f64) -> Result<Channel, GeoError> {
        ramp(self.x0, peak, self.rate)
    }
    fn fall(&self, peak: f64) -> Result<Channel, GeoError> {
        ramp(peak, 0.0, self.rate)
    }
}

/// S-curves (the vertical rate ḣ, whose integral is the height change).
pub(crate) struct ScurveShape {
    /// Starting value.
    pub x0: f64,
    /// Starting derivative.
    pub d0: f64,
    /// Derivative limit.
    pub max_d: f64,
    /// Second-derivative limit.
    pub jerk: f64,
}

impl Shape for ScurveShape {
    fn rise(&self, peak: f64) -> Result<Channel, GeoError> {
        let steps = scurve(self.x0, self.d0, peak, self.max_d, self.jerk)?;
        Ok(Channel::from_steps(self.x0, self.d0, &steps))
    }
    fn fall(&self, peak: f64) -> Result<Channel, GeoError> {
        let steps = scurve(peak, 0.0, 0.0, self.max_d, self.jerk)?;
        Ok(Channel::from_steps(peak, 0.0, &steps))
    }
}

/// Concatenates rise, a hold of `hold` seconds at the rise's final value, and fall.
fn join(rise: &Channel, hold: f64, fall: &Channel) -> Channel {
    let mut out = rise.clone();
    let (mut start, mut integral) = (rise.duration(), rise.total_integral());
    let x = rise.final_value();
    if hold > 0.0 {
        out.phases.push(Phase {
            start,
            duration: hold,
            x0: x,
            d0: 0.0,
            jerk: 0.0,
            integral0: integral,
        });
        start += hold;
        integral += x * hold;
    }
    for p in &fall.phases {
        out.phases.push(Phase {
            start: start + p.start,
            integral0: integral + p.integral0,
            ..*p
        });
    }
    out
}

/// Plans rise → hold → fall with total integral `delta` and |peak| ≤ `peak_max`.
///
/// If the integral at `±peak_max` (towards `delta`) with no hold falls short, the peak
/// is `±peak_max` and the hold makes up the rest. Otherwise the peak is a root of
/// `integral(peak) = delta`. The integral need not be monotone in the peak (a start
/// with non-zero rate and acceleration has interior extrema), so the peak range is
/// scanned from `±peak_max` towards `∓peak_max` in [`MOVE_BY_SCAN`] steps. The first
/// sign change is then refined by 200 bisection steps, which picks the root with the
/// largest peak towards `delta` (the fastest). No root is an error (the start state
/// overshoots `delta` whatever the peak).
pub(crate) fn move_by(shape: &dyn Shape, delta: f64, peak_max: f64) -> Result<Channel, GeoError> {
    check_limits("profile", &[peak_max])?;
    if !delta.is_finite() {
        return Err(GeoError::NonFinite("profile target"));
    }
    let area = |peak: f64| -> Result<(Channel, Channel, f64), GeoError> {
        let r = shape.rise(peak)?;
        let f = shape.fall(peak)?;
        let a = r.total_integral() + f.total_integral();
        Ok((r, f, a))
    };
    let s = if delta >= 0.0 { 1.0 } else { -1.0 };
    let (r, f, a_max) = area(s * peak_max)?;
    if s * a_max <= s * delta {
        let hold = (delta - a_max) / (s * peak_max);
        return Ok(join(&r, hold, &f));
    }
    // s·(area − delta) > 0 at the top; find where it first becomes ≤ 0.
    let mut hi = s * peak_max;
    let mut lo = None;
    for i in 1..=MOVE_BY_SCAN {
        let c = s * peak_max * (1.0 - 2.0 * f64::from(i) / f64::from(MOVE_BY_SCAN));
        let (_, _, a) = area(c)?;
        if s * (a - delta) <= 0.0 {
            lo = Some(c);
            break;
        }
        hi = c;
    }
    let Some(mut lo) = lo else {
        return Err(GeoError::Trajectory(format!(
            "a change of {delta} cannot be met without overshooting from the current state"
        )));
    };
    // Invariant: s·(area(lo) − delta) ≤ 0 < s·(area(hi) − delta).
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if mid == lo || mid == hi {
            break;
        }
        let (_, _, a) = area(mid)?;
        if s * (a - delta) <= 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let (r, f, _) = area(0.5 * (lo + hi))?;
    Ok(join(&r, 0.0, &f))
}

/// Scan steps across the peak range in [`move_by`].
pub(crate) const MOVE_BY_SCAN: u32 = 64;

/// Rise to `peak`, hold, fall, lasting `total` seconds; an error if the transitions
/// alone are longer.
pub(crate) fn hold_for(
    shape: &dyn Shape,
    peak: f64,
    total: f64,
    what: &str,
) -> Result<Channel, GeoError> {
    let r = shape.rise(peak)?;
    let f = shape.fall(peak)?;
    let hold = total - r.duration() - f.duration();
    if hold < 0.0 {
        return Err(GeoError::Trajectory(format!(
            "{what}: the duration ({total} s) is shorter than its transitions ({} s)",
            r.duration() + f.duration()
        )));
    }
    Ok(join(&r, hold, &f))
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    fn check_scurve(x0: f64, d0: f64, x1: f64, a: f64, j: f64) {
        let steps = scurve(x0, d0, x1, a, j).unwrap();
        let c = Channel::from_steps(x0, d0, &steps);
        let end = c.final_value();
        let scale = x0.abs().max(x1.abs()).max(d0 * d0 / j).max(1.0);
        assert!((end - x1).abs() <= 1e-9 * scale, "end {end} vs {x1}");
        let (_, d_end, _, _) = c.at(c.duration());
        assert!(
            d_end.abs() <= 1e-9 * a.max(d0.abs()),
            "end derivative {d_end}"
        );
        // Derivative continuous at every breakpoint and bounded.
        let lim = a.max(d0.abs()) * (1.0 + 1e-9);
        for b in c.breakpoints() {
            let l = c.at(b - 1e-9);
            let r = c.at(b + 1e-9);
            assert!((l.0 - r.0).abs() < 1e-6 * scale);
            assert!((l.1 - r.1).abs() < 1e-6 * lim);
            assert!(l.1.abs() <= lim && r.1.abs() <= lim);
        }
    }

    #[test]
    fn scurve_cases() {
        check_scurve(100.0, 0.0, 150.0, 1.0, 0.5); // trapezoidal
        check_scurve(100.0, 0.0, 100.5, 1.0, 0.5); // triangular
        check_scurve(0.0, 0.8, 0.1, 1.0, 0.5); // overshoot: reverses
        check_scurve(0.0, -0.8, 5.0, 1.0, 0.5); // starts the wrong way
        check_scurve(3.0, 0.0, 3.0, 1.0, 0.5); // nothing to do
        check_scurve(0.0, 2.0, 50.0, 1.0, 0.5); // starts above the limit
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(2048))]
        #[test]
        fn scurve_properties(
            x0 in -300.0..300.0_f64,
            d0 in -2.0..2.0_f64,
            x1 in -300.0..300.0_f64,
            a in 0.05..3.0_f64,
            j in 0.05..2.0_f64,
        ) {
            check_scurve(x0, d0, x1, a, j);
        }

        #[test]
        fn move_by_meets_its_integral(
            d0 in -0.02..0.02_f64,
            delta in -6.0..6.0_f64,
            rate in 0.005..0.05_f64,
        ) {
            let peak = 0.0524;
            let shape = RampShape { x0: d0, rate };
            if let Ok(c) = move_by(&shape, delta, peak) {
                prop_assert!((c.total_integral() - delta).abs() < 1e-9);
                prop_assert!(c.final_value().abs() < 1e-15);
                prop_assert!(c.phases.iter().all(|p| p.x0.abs() <= peak.max(d0.abs()) + 1e-12));
            } else {
                // Reversing at the full peak always reaches these targets.
                prop_assert!(false, "move_by failed");
            }
        }
    }

    #[test]
    fn move_by_finds_roots_of_a_non_monotone_integral() {
        // Found by review: from ḣ = 0.4 m/s, ḧ = −0.5 m/s², a +0.05 m change capped at
        // 0.1 m/s is reachable (peak ≈ −0.0964 m/s) though the integral at −0.1 m/s
        // exceeds the target.
        let shape = ScurveShape {
            x0: 0.4,
            d0: -0.5,
            max_d: 0.5,
            jerk: 0.25,
        };
        let c = move_by(&shape, 0.05, 0.1).unwrap();
        assert!(
            (c.total_integral() - 0.05).abs() < 1e-9,
            "{}",
            c.total_integral()
        );
        assert!(c.final_value().abs() < 1e-12);
    }

    #[test]
    fn climb_by_height() {
        let shape = ScurveShape {
            x0: 0.0,
            d0: 0.0,
            max_d: 0.5,
            jerk: 0.25,
        };
        for dh in [600.0, -600.0, 3.0, 0.01] {
            let c = move_by(&shape, dh, 7.62).unwrap();
            assert!(
                (c.total_integral() - dh).abs() < 1e-8 * dh.abs().max(1.0),
                "{dh}"
            );
            assert!(c.final_value().abs() < 1e-12);
            let (_, _, _, i) = c.at(c.duration() + 10.0);
            assert!((i - dh).abs() < 1e-8 * dh.abs().max(1.0));
        }
        let t = hold_for(&shape, 5.0, 100.0, "climb").unwrap();
        assert!((t.duration() - 100.0).abs() < 1e-12);
        assert!(hold_for(&shape, 5.0, 5.0, "climb").is_err());
    }

    #[test]
    fn ramp_ends_exactly() {
        let r = ramp(0.01, -0.03, 0.007).unwrap();
        assert_eq!(r.final_value(), -0.03);
        assert!((r.duration() - 0.04 / 0.007).abs() < 1e-12);
        assert_eq!(ramp(0.2, 0.2, 1.0).unwrap().duration(), 0.0);
    }
}
