//! One-path light time in the pseudo-inertial ECI frame ([`crate::eci`]).
//!
//! Each endpoint is rotated at its own instant, including the Earth-rotation term
//! in its inertial velocity. Time offsets use signed Q64.64 picoseconds; absolute
//! truth time is never multiplied by a floating-point light time.
//!
//! The fixed-point iteration stops after computing τₙ₊₁ when
//! `|τₙ₊₁ − τₙ| ≤ max(1e-16 s, 4·ulp(τₙ₊₁))`, within ten steps. For a
//! contraction factor β ≈ |v_tx|/c (|v_rx|/c for emission solves), the remaining
//! error in τₙ₊₁ is bounded by β/(1−β) times the last step, apart from f64
//! rounding. LEO has β ≲ 3e-5 and usually takes three to five steps. We evaluate
//! the moving endpoint once more at the converged offset, and report its actual
//! distance and that distance divided by c. Endpoint time and the reported delay
//! consequently agree to the fixed-point residual and rounding, rather than
//! carrying the preceding iterate's state.

use iqsim_core::AcqInstant;

use crate::eci::{EarthOrientation, Eci};
use crate::geodetic::SPEED_OF_LIGHT_M_S;
use crate::{GeoError, State, Trajectory, Vec3};

/// Maximum number of fixed-point updates in a light-time solve.
pub const LIGHT_TIME_MAX_ITERATIONS: u32 = 10;

/// Seconds as signed Q64.64 picoseconds, rounding ties to even.
///
/// First rounds `s * 1e12` to f64 picoseconds, then scales exactly by 2⁶⁴ and
/// rounds to an integer. Non-finite inputs and magnitudes above 10⁶ s are errors.
pub fn seconds_to_q64_ps(s: f64) -> Result<i128, GeoError> {
    if !s.is_finite() {
        return Err(GeoError::NonFinite("light-time offset in seconds"));
    }
    if s.abs() > 1.0e6 {
        return Err(GeoError::LightTime("offset exceeds 1e6 seconds".into()));
    }
    let x = s * 1.0e12;
    let scaled = x * 18_446_744_073_709_551_616.0;
    // |scaled| <= 1e18 * 2^64 < 2^127, so the conversion cannot saturate.
    Ok(scaled.round_ties_even() as i128)
}

/// Add seconds to an acquisition instant, preserving its fractional picoseconds.
///
/// Uses checked Q64.64 addition; leaving the i64 whole-picosecond range is an
/// error. The offset has the bounds and rounding of [`seconds_to_q64_ps`].
pub fn instant_offset(t: AcqInstant, s: f64) -> Result<AcqInstant, GeoError> {
    let offset = seconds_to_q64_ps(s)?;
    let sum = t
        .to_q64_64()
        .checked_add(offset)
        .ok_or_else(|| GeoError::LightTime("acquisition instant overflow".into()))?;
    // Every i128 is a valid signed Q64.64 instant, with an i64 integer part.
    Ok(AcqInstant::from_q64_64(sum))
}

/// One endpoint of a propagation path, evaluated at its own instant.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PathEnd {
    /// Acquisition instant of this endpoint.
    pub instant: AcqInstant,
    /// Earth-fixed trajectory state at this instant.
    pub state: State,
    /// Position rotated into ECI at this instant, metres.
    pub position_eci: Eci,
    /// Inertial velocity, including Earth rotation, in ECI axes, m/s.
    pub velocity_eci: Vec3,
}

impl PathEnd {
    fn at(trajectory: &dyn Trajectory, instant: AcqInstant) -> Result<Self, GeoError> {
        let state = trajectory.state(instant)?;
        let (position_eci, velocity_eci) =
            EarthOrientation::at(instant).state_to_eci(state.position, state.velocity);
        Ok(Self {
            instant,
            state,
            position_eci,
            velocity_eci,
        })
    }
}

/// Solved propagation delay and endpoint states.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LightTime {
    /// Light time τ in seconds, equal to `range_m / c`.
    pub tau_s: f64,
    /// dτ/dt_rx, dimensionless; received carrier Doppler factor is `1 − tau_rate`.
    ///
    /// Differentiating `cτ = |r_rx(t) − r_tx(t − τ)|` gives
    /// `cτ' = û·(v_rx − v_tx(1 − τ'))`, hence
    /// `τ' = û·(v_rx − v_tx) / (c − û·v_tx)`, with û pointing TX → RX.
    pub tau_rate: f64,
    /// Distance between the returned ECI positions, metres.
    pub range_m: f64,
    /// Receiver at the receive instant.
    pub rx: PathEnd,
    /// Transmitter at the emission instant.
    pub tx: PathEnd,
    /// Number of fixed-point updates (excluding initialization and final evaluation).
    pub iterations: u32,
}

fn distance(rx: &PathEnd, tx: &PathEnd) -> Result<f64, GeoError> {
    let range = rx.position_eci.vec().distance(tx.position_eci.vec());
    if !range.is_finite() {
        return Err(GeoError::NonFinite("light-time range"));
    }
    if range == 0.0 {
        return Err(GeoError::LightTime("coincident path endpoints".into()));
    }
    Ok(range)
}

// x is positive and finite, as checked by distance().
fn ulp(x: f64) -> f64 {
    f64::from_bits(x.to_bits() + 1) - x
}

fn finish(rx: PathEnd, tx: PathEnd, iterations: u32) -> Result<LightTime, GeoError> {
    let range_m = distance(&rx, &tx)?;
    let unit = (rx.position_eci.vec() - tx.position_eci.vec()) * (1.0 / range_m);
    let denominator = SPEED_OF_LIGHT_M_S - unit.dot(tx.velocity_eci);
    let tau_rate = unit.dot(rx.velocity_eci - tx.velocity_eci) / denominator;
    if denominator <= 0.0 || !tau_rate.is_finite() {
        return Err(GeoError::LightTime("invalid light-time derivative".into()));
    }
    Ok(LightTime {
        tau_s: range_m / SPEED_OF_LIGHT_M_S,
        tau_rate,
        range_m,
        rx,
        tx,
        iterations,
    })
}

/// `t` moved into `trajectory`'s span, if it lies outside.
fn clamp_into(trajectory: &dyn Trajectory, t: AcqInstant) -> AcqInstant {
    let span = trajectory.span();
    let mut t = t;
    if let Some(s) = span.start {
        t = t.max(AcqInstant::from_truth_time(s));
    }
    if let Some(e) = span.end {
        t = t.min(AcqInstant::from_truth_time(e));
    }
    t
}

fn solve(
    rx: &dyn Trajectory,
    tx: &dyn Trajectory,
    instant: AcqInstant,
    receive: bool,
) -> Result<LightTime, GeoError> {
    // The moving endpoint starts at the fixed instant clamped into its span (review
    // fix: a transmitter knot before the receiver's start, or an emission before the
    // transmitter's end, made the first guess fail although the solution is valid).
    let (rx_at, tx_at) = if receive {
        (instant, clamp_into(tx, instant))
    } else {
        (clamp_into(rx, instant), instant)
    };
    let mut rx_end = PathEnd::at(rx, rx_at)?;
    let mut tx_end = PathEnd::at(tx, tx_at)?;
    let mut tau = distance(&rx_end, &tx_end)? / SPEED_OF_LIGHT_M_S;
    for iterations in 1..=LIGHT_TIME_MAX_ITERATIONS {
        if receive {
            tx_end = PathEnd::at(tx, instant_offset(instant, -tau)?)?;
        } else {
            rx_end = PathEnd::at(rx, instant_offset(instant, tau)?)?;
        }
        let next = distance(&rx_end, &tx_end)? / SPEED_OF_LIGHT_M_S;
        if (next - tau).abs() <= 1.0e-16_f64.max(4.0 * ulp(next)) {
            if receive {
                tx_end = PathEnd::at(tx, instant_offset(instant, -next)?)?;
            } else {
                rx_end = PathEnd::at(rx, instant_offset(instant, next)?)?;
            }
            return finish(rx_end, tx_end, iterations);
        }
        tau = next;
    }
    Err(GeoError::NoConvergence {
        what: "light time",
        iterations: LIGHT_TIME_MAX_ITERATIONS,
    })
}

/// Solve `|r_rx(t_rx) − r_tx(t_rx − τ)| = cτ` in ECI at a fixed receive instant.
///
/// Uses the stopping rule in the module docs. Trajectory errors propagate unchanged;
/// coincidence, invalid offsets and nonconvergence are errors.
pub fn solve_receive(
    rx: &dyn Trajectory,
    tx: &dyn Trajectory,
    t_rx: AcqInstant,
) -> Result<LightTime, GeoError> {
    solve(rx, tx, t_rx, true)
}

/// Solve `|r_rx(t_tx + τ) − r_tx(t_tx)| = cτ` in ECI at a fixed emission instant.
///
/// Maps an emission knot to its receive time, using the module's stopping rule.
/// The returned [`LightTime::tau_rate`] remains dτ/dt_rx at the solved reception.
/// Trajectory errors propagate unchanged.
pub fn solve_emit(
    rx: &dyn Trajectory,
    tx: &dyn Trajectory,
    t_tx: AcqInstant,
) -> Result<LightTime, GeoError> {
    solve(rx, tx, t_tx, false)
}
