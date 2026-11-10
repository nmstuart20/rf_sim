//! Light-time path geometry and per-block polynomial fitting (prompt 9, contracts §14).
//!
//! For one path (an emitter's transmitter trajectory into a receiver trajectory) and
//! one block of a render target, [`fit_path`] produces the fixed-point polynomials the
//! renderer evaluates: delay, retarded waveform coordinate, carrier phase and
//! amplitude, as one or more [`Segment`]s that tile the block.
//!
//! # Exact model at a sample
//!
//! Sample k of the block is acquired at `t = ClockMap.instant(n0 + k)` (Q64.64 ps, the
//! fraction kept). The light time τ solves `|r_rx(t) − r_tx(t − τ)| = c·τ` in ECI
//! ([`iqsim_geo::light_time::solve_receive`]), with its rate τ̇ = dτ/dt. Then:
//!
//! - **delay** `τ·F_s` target samples;
//! - **phase** `φ = f_off·t − f_e·τ` cycles, plus for CW the response fold
//!   `−f_inst·τ_g` with the instantaneous offset `f_inst = f_off − f_e·τ̇`. The exact
//!   value is computed with modular integers only: [`OffsetFrequency::phase_at`] for
//!   `f_off·t`, and [`carrier_phase`] (an exact integer product of the µHz carrier and
//!   τ as Q64.64 ps) for `f_e·τ`. No absolute time or large cycle count ever meets a
//!   floating-point product;
//! - **waveform coordinate** `F·(t − t_ref) − F·τ` (F the coordinate rate: the target
//!   rate for pulses, the native rate for continuous sources), as an exact Q32.32
//!   integer: the interval in exact rational arithmetic minus `F·τ` rounded once;
//! - **amplitude** `√(P·R₀)·λ/(4π·c·τ)` volts (Friis with unit gains at the light-time
//!   range), times `|H(f_inst)|` for CW.
//!
//! # Fit
//!
//! An interval of n samples is interpolated by a cubic in the segment-local index j
//! through 4 **Chebyshev–Lobatto nodes** `u = 0, ¼, ¾, 1` (`j = round(u·(n − 1))`,
//! duplicates dropped for short intervals, so the degree is ≤ min(3, n − 1)). The
//! nodes are the extrema of the degree-3 Chebyshev polynomial, which include the
//! interval ends, so neighbouring segments and blocks join at exactly evaluated
//! values. Each quantity is fitted as its *variation* from the first node, then:
//!
//! - delay: D0 = τ₀·F_s, D1..D3 from the variation;
//! - phase: D0 is the exact modular phase at j = 0, D1..D3 from the variation in f64;
//! - coordinate: the exact Q32.32 value at j = 0 is kept separately as an i128 (the
//!   caller rebases it), D1..D3 from the variation;
//! - amplitude: an f32 cubic in `x = f32(j)·f32(1/n)` (the kernel's own x).
//!
//! # Checks and refinement (contracts §14.3)
//!
//! Every fit is compared with the exact model at 15 points: the ends, the inner nodes,
//! the three extrema of the nodal polynomial `ω(u) = u(u − ¼)(u − ¾)(u − 1)` (u = ½
//! and ½ ± √(5/32)), where the interpolation error of a smooth function peaks, and the
//! eight odd sixteenths. The fixed-point polynomials are evaluated exactly as the
//! kernel does. An interval is accepted when every residual is at most a quarter of
//! its budget ([`FIT_MARGIN`]): delay 10⁻³ samples, phase 10⁻² rad, coordinate 10⁻³
//! target samples, amplitude 0.01 dB. Point checks are not a proof between points;
//! the margin is the conservative part: on an interval without known discontinuities
//! the error is ≈ f⁽⁴⁾/24·ω(u)·T⁴, so it stays inside the budget unless f⁽⁴⁾ varies by
//! more than 4× across the interval. Otherwise the interval is split in half and both
//! halves are fitted again, down to [`MAX_REFINE_DEPTH`] halvings; beyond that the fit
//! fails ([`FitError::Budget`]).
//!
//! Before any refinement, the block is split at the receiver's trajectory knots and at
//! the transmitter's knots mapped from emission to receive time
//! ([`iqsim_geo::light_time::solve_emit`]): a jerk discontinuity is never inside a
//! fitted interval. Every split is logged ([`Split`]).

use std::fmt;
use std::sync::Arc;

use iqsim_core::{
    AcqInstant, ClockMap, DelayPoly, NumericError, OffsetFrequency, Phase, PhasePoly, SampleIndex,
    SampleRate, TapResponse, TimeError, TruthTime,
};
use iqsim_geo::light_time::{self, LightTime};
use iqsim_geo::{GeoError, Trajectory};

use iqsim_worklist::AmplitudePoly;

/// Delay budget, target samples (contracts §10).
pub const DELAY_BUDGET_SAMPLES: f64 = 1e-3;
/// Carrier-phase budget, radians (contracts §10).
pub const PHASE_BUDGET_RAD: f64 = 1e-2;
/// Amplitude budget above the null floor: 0.01 dB as a relative error,
/// 10^(0.01/20) − 1 rounded down.
pub const AMPLITUDE_BUDGET_REL: f64 = 1.1519e-3;
/// A fit is accepted when every residual is at most budget / `FIT_MARGIN`.
pub const FIT_MARGIN: f64 = 4.0;
/// Most halvings of one knot-free interval before the fit fails.
pub const MAX_REFINE_DEPTH: u32 = 12;

/// Fit nodes in u ∈ [0, 1]: the degree-3 Chebyshev–Lobatto points.
pub const NODES: [f64; 4] = [0.0, 0.25, 0.75, 1.0];

/// Check points in u ∈ [0, 1] (module docs).
pub const CHECKS: [f64; 15] = [
    0.0,
    0.0625,
    0.104_715_292_478_522_3,
    0.1875,
    0.25,
    0.3125,
    0.4375,
    0.5,
    0.5625,
    0.6875,
    0.75,
    0.8125,
    0.895_284_707_521_477_7,
    0.9375,
    1.0,
];

const TWO_POW_32: f64 = 4_294_967_296.0;
const TWO_POW_64: f64 = 18_446_744_073_709_551_616.0;

/// The two ends of a path.
#[derive(Clone, Debug)]
pub struct PathGeometry {
    /// Receiver trajectory (receive instants).
    pub rx: Arc<dyn Trajectory>,
    /// Transmitter trajectory (emission instants; normally a
    /// [`iqsim_geo::trajectory::Lookback`] so light time can look before its start).
    pub tx: Arc<dyn Trajectory>,
}

impl PathGeometry {
    /// The light-time solution for the sample acquired at `t`.
    pub fn at(&self, t: AcqInstant) -> Result<LightTime, GeoError> {
        light_time::solve_receive(&*self.rx, &*self.tx, t)
    }
}

/// How a path's waveform coordinate is measured.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoordModel {
    /// Coordinate rate F (target rate for pulses, native rate for continuous sources).
    pub rate: SampleRate,
    /// Origin `t_ref` of the source segment (emission reference or sample-0 time).
    pub origin: TruthTime,
}

impl CoordModel {
    fn rate_hz(&self) -> f64 {
        self.rate.num() as f64 / self.rate.den() as f64
    }
}

/// Everything the fit needs about one path into one target.
#[derive(Clone, Debug)]
pub struct PathModel<'a> {
    /// The two trajectories.
    pub geometry: &'a PathGeometry,
    /// Emitter carrier f_e, µHz.
    pub carrier_uhz: i64,
    /// f_off = f_e − f_c, exact.
    pub offset: OffsetFrequency,
    /// f_off, Hz.
    pub offset_hz: f64,
    /// Target sample rate.
    pub target_rate: SampleRate,
    /// The waveform coordinate, if the emission has one (pulses, continuous).
    pub coord: Option<CoordModel>,
    /// f_off, µHz (the nominal response offset).
    pub offset_uhz: i64,
    /// For CW: the target response, folded at the instantaneous offset with its group
    /// delay applied exactly (the geometry is evaluated at `t − τ_g`).
    pub fold: Option<&'a TapResponse>,
    /// Pulses and continuous sources: the target response their tables and kernels are
    /// built from at a segment's Doppler-shifted offset ([`Segment::offset_uhz`]), and
    /// the largest |f_inst − offset| that keeps the response error within a quarter of
    /// the budgets ([`doppler_tolerance_hz`]). `None` for CW.
    pub filtered: Option<(&'a TapResponse, f64)>,
    /// Amplitude × range: `√(P·R₀)·λ/(4π)`, V·m.
    pub volts_m: f64,
}

/// Fit residuals (maximum over the check points).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Residuals {
    /// Delay, target samples.
    pub delay_samples: f64,
    /// Carrier phase, radians.
    pub phase_rad: f64,
    /// Waveform coordinate, target samples (0 without a coordinate).
    pub coord_samples: f64,
    /// Amplitude: relative above the null floor; below it the error over
    /// 10⁻⁵·V_ref/[`AMPLITUDE_BUDGET_REL`] (so the acceptance is the contracts §2
    /// absolute bound there).
    pub amplitude_rel: f64,
    /// Response mismatch from Doppler: for CW the largest change of |H| between the
    /// interval's instantaneous offsets and its key, relative (as `amplitude_rel`); for
    /// pulses and continuous sources the largest |f_inst − key| over the budget's
    /// tolerance (accepted at ≤ 1).
    pub response: f64,
}

impl Residuals {
    /// The componentwise maximum.
    pub fn max(self, o: Self) -> Self {
        Self {
            delay_samples: self.delay_samples.max(o.delay_samples),
            phase_rad: self.phase_rad.max(o.phase_rad),
            coord_samples: self.coord_samples.max(o.coord_samples),
            amplitude_rel: self.amplitude_rel.max(o.amplitude_rel),
            response: self.response.max(o.response),
        }
    }

    /// Whether every residual is within budget / [`FIT_MARGIN`].
    pub fn accepted(&self) -> bool {
        self.delay_samples <= DELAY_BUDGET_SAMPLES / FIT_MARGIN
            && self.phase_rad <= PHASE_BUDGET_RAD / FIT_MARGIN
            && self.coord_samples <= DELAY_BUDGET_SAMPLES / FIT_MARGIN
            && self.amplitude_rel <= AMPLITUDE_BUDGET_REL / FIT_MARGIN
            && self.response <= 1.0
    }
}

impl fmt::Display for Residuals {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "delay {:.3e} samples, phase {:.3e} rad, coordinate {:.3e} samples, amplitude {:.3e}, response {:.3e}",
            self.delay_samples,
            self.phase_rad,
            self.coord_samples,
            self.amplitude_rel,
            self.response
        )
    }
}

/// A fitted waveform coordinate: the exact Q32.32 value at the segment start, and the
/// variation polynomial (D0 = 0).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Coord {
    /// Exact coordinate at j = 0, Q32.32 bits in the coordinate's units.
    pub start: i128,
    /// `w(j) − w(0)`, Q32.96 registers with D0 = 0.
    pub variation: DelayPoly,
}

impl Coord {
    /// The coordinate at segment index j, Q32.32 bits (as the kernel computes it).
    pub fn at(&self, j: u64) -> Result<i128, NumericError> {
        Ok(self.start + i128::from(self.variation.eval(j)?.to_bits()))
    }
}

/// One fitted interval of a block.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Segment {
    /// Block-local index of the first sample.
    pub first: u64,
    /// Samples.
    pub len: u64,
    /// Delay, target samples.
    pub delay: DelayPoly,
    /// Carrier phase.
    pub phase: PhasePoly,
    /// Waveform coordinate, if the path has one.
    pub coord: Option<Coord>,
    /// Amplitude, volts.
    pub amplitude: AmplitudePoly,
    /// Offset the segment's edge tables or resampling kernel are built at, µHz: f_off
    /// plus the interval's Doppler shift rounded to the tolerance grid (exactly f_off
    /// when there is no Doppler). Unused for CW.
    pub offset_uhz: i64,
    /// Worst residuals at the check points.
    pub residuals: Residuals,
}

/// Why a block's path was split.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SplitReason {
    /// A receiver trajectory knot.
    ReceiverKnot,
    /// A transmitter trajectory knot, mapped to receive time.
    TransmitterKnot,
    /// The fit missed its budget; the interval was halved.
    Refine,
}

impl SplitReason {
    /// Name used in logs and the truth record.
    pub const fn name(self) -> &'static str {
        match self {
            Self::ReceiverKnot => "rx_knot",
            Self::TransmitterKnot => "tx_knot",
            Self::Refine => "refine",
        }
    }
}

/// One split of one block's path.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Split {
    /// Block-local index where the later interval starts.
    pub at: u64,
    /// Why.
    pub reason: SplitReason,
    /// For a refinement: the residuals of the interval that was halved.
    pub residuals: Option<Residuals>,
    /// For a refinement: that interval's first sample and length.
    pub interval: Option<(u64, u64)>,
}

/// A fit failure.
#[derive(Clone, Debug, PartialEq)]
pub enum FitError {
    /// Geometry (light time, trajectory) failed.
    Geo(GeoError),
    /// A fixed-point conversion failed.
    Numeric(NumericError),
    /// A clock or index operation failed.
    Time(TimeError),
    /// The delay leaves the Q32.32 range (with the contracts §1.3 margin).
    DelayRange {
        /// Delay, target samples.
        samples: f64,
    },
    /// The budget could not be met at the refinement limit.
    Budget {
        /// Block-local first sample of the failing interval.
        first: u64,
        /// Its length.
        len: u64,
        /// Its residuals.
        residuals: Residuals,
    },
}

impl fmt::Display for FitError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Geo(e) => e.fmt(f),
            Self::Numeric(e) => e.fmt(f),
            Self::Time(e) => e.fmt(f),
            Self::DelayRange { samples } => {
                write!(
                    f,
                    "path delay {samples} samples is outside the Q32.32 range"
                )
            }
            Self::Budget {
                first,
                len,
                residuals,
            } => write!(
                f,
                "the geometry fit misses its budget after {MAX_REFINE_DEPTH} halvings \
                 (samples {first}..{} of the block: {residuals})",
                first + len
            ),
        }
    }
}

impl std::error::Error for FitError {}

impl From<GeoError> for FitError {
    fn from(e: GeoError) -> Self {
        Self::Geo(e)
    }
}
impl From<NumericError> for FitError {
    fn from(e: NumericError) -> Self {
        Self::Numeric(e)
    }
}
impl From<TimeError> for FitError {
    fn from(e: TimeError) -> Self {
        Self::Time(e)
    }
}

/// `f_e·τ` mod 1 cycle, in 2⁻⁶⁴ cycle units, floored: `⌊f_uhz·d / 10¹⁸⌋ mod 2⁶⁴` with d
/// = τ as Q64.64 ps ([`light_time::seconds_to_q64_ps`]). Exact for the given τ: a
/// 64 × 128-bit product in three 64-bit limbs, then long division by 10¹⁸ < 2⁶⁰.
pub fn carrier_phase(carrier_uhz: i64, tau_s: f64) -> Result<Phase, FitError> {
    let d = light_time::seconds_to_q64_ps(tau_s)?;
    if carrier_uhz < 0 || d < 0 {
        return Err(FitError::Numeric(NumericError::OutOfRange));
    }
    let f = carrier_uhz as u128;
    let d = d as u128;
    let (d_hi, d_lo) = (d >> 64, d & u128::from(u64::MAX));
    let lo = f * d_lo;
    let hi = f * d_hi;
    let mid = (lo >> 64) + (hi & u128::from(u64::MAX));
    let limbs = [
        (hi >> 64) + (mid >> 64),
        mid & u128::from(u64::MAX),
        lo & u128::from(u64::MAX),
    ];
    const D: u128 = 1_000_000_000_000_000_000;
    let mut r: u128 = 0;
    let mut q: u128 = 0;
    for limb in limbs {
        let cur = (r << 64) | limb;
        q = cur / D;
        r = cur % D;
    }
    Ok(Phase::from_bits(q as u64))
}

/// Q32.32 samples at `rate` in the interval `ps + frac·2⁻⁶⁴` picoseconds:
/// ⌊dt·num/(den·10¹²)·2³²⌋ to within one unit (the low 32 bits of the fractional
/// product are dropped).
///
/// The whole picoseconds come separately (review fix: two instants ±100 days apart
/// overflow a single Q64.64 i128). |ps| < 2⁶⁵, so `ps·num` < 2¹¹³, and the remainder
/// pieces stay below 2¹⁰⁵: no intermediate leaves i128.
pub(crate) fn interval_to_q32(ps: i128, frac: u64, rate: SampleRate) -> i128 {
    let frac = u128::from(frac);
    let num = i128::from(rate.num());
    let d = i128::from(rate.den()) * 1_000_000_000_000;
    let a = ps * num;
    let q = a.div_euclid(d);
    let r = a.rem_euclid(d);
    let f = ((frac * rate.num() as u128) >> 32) as i128;
    (q << 32) + ((r << 32) + f) / d
}

/// The largest |f_inst − key| (Hz) for which building a pulse's tables or a continuous
/// source's kernel at `key` instead of the instantaneous offset keeps the response
/// error within a quarter of the budgets: |ΔH| ≤ max|H′|·δ with max|H′| = 1/(√(2π)·σ)
/// for the Gaussian transition (≤ ¼ of the 0.01 dB amplitude budget relative to the
/// passband), and the group-delay phase 2π·δ·τ_g ≤ ¼ of the phase budget.
pub fn doppler_tolerance_hz(h: &TapResponse) -> f64 {
    let amp =
        AMPLITUDE_BUDGET_REL / FIT_MARGIN * libm::sqrt(2.0 * core::f64::consts::PI) * h.sigma_hz();
    let tg = h.group_delay_s();
    if tg > 0.0 {
        amp.min(PHASE_BUDGET_RAD / FIT_MARGIN / (2.0 * core::f64::consts::PI * tg))
    } else {
        amp
    }
}

/// The exact model at one sample.
#[derive(Clone, Copy, Debug)]
struct Eval {
    /// Instant the geometry is evaluated at (the acquisition instant, or for CW with a
    /// group delay `t − τ_g`), Q64.64 ps.
    t: i128,
    tau: f64,
    /// Instantaneous offset f_off − f_e·τ̇, Hz.
    f_inst: f64,
    /// Exact phase, including the CW fold.
    phase: Phase,
    /// Exact coordinate, Q32.32 bits (0 without a coordinate).
    coord: i128,
    /// Amplitude, volts.
    amplitude: f64,
    /// V_ref: the amplitude without the response (Friis at the light-time range).
    v_ref: f64,
}

impl PathModel<'_> {
    fn target_rate_hz(&self) -> f64 {
        self.target_rate.num() as f64 / self.target_rate.den() as f64
    }

    fn carrier_hz(&self) -> f64 {
        self.carrier_uhz as f64 * 1e-6
    }

    fn eval(&self, clock: &dyn ClockMap, n: SampleIndex) -> Result<Eval, FitError> {
        let acq = clock.instant(n)?;
        // A linear-phase response delays an in-band tone's whole trajectory by τ_g
        // (review fix: the first-order fold −f_inst·τ_g omitted ½τ_g²·φ̈, 0.03 rad for a
        // 31 GHz LEO pass with τ_g = 1 ms). So CW geometry is evaluated at t − τ_g.
        let t = match self.fold {
            Some(h) if h.group_delay_s() > 0.0 => {
                light_time::instant_offset(acq, -h.group_delay_s())?
            }
            _ => acq,
        };
        let lt = self.geometry.at(t)?;
        let tau = lt.tau_s;
        let f_inst = self.offset_hz - self.carrier_hz() * lt.tau_rate;
        let gain = self.fold.map_or(1.0, |h| h.magnitude(f_inst));
        let phase = self
            .offset
            .phase_at(t)
            .wrapping_sub(carrier_phase(self.carrier_uhz, tau)?);
        let coord = match self.coord {
            Some(c) => {
                let ps = i128::from(t.ps().as_ps()) - i128::from(c.origin.as_ps());
                let shift = iqsim_core::f64_to_i64_ties_even(tau * c.rate_hz() * TWO_POW_32)?;
                interval_to_q32(ps, t.frac(), c.rate) - i128::from(shift)
            }
            None => 0,
        };
        let v_ref = self.volts_m / lt.range_m;
        let amplitude = v_ref * gain;
        Ok(Eval {
            t: t.to_q64_64(),
            tau,
            f_inst,
            phase,
            coord,
            amplitude,
            v_ref,
        })
    }

    /// The exact waveform coordinate (Q32.32 bits) at channel sample `n`.
    pub fn exact_coord(&self, clock: &dyn ClockMap, n: SampleIndex) -> Result<i128, FitError> {
        Ok(self.eval(clock, n)?.coord)
    }
}

/// min and max of the power-basis cubic `a` over x ∈ [0, x_max]: the ends and the
/// derivative's real roots inside.
fn cubic_range(a: [f64; 4], x_max: f64) -> (f64, f64) {
    let p = |x: f64| a[0] + x * (a[1] + x * (a[2] + x * a[3]));
    let mut xs = vec![0.0, x_max];
    // p′ = a1 + 2a2·x + 3a3·x².
    let (qa, qb, qc) = (3.0 * a[3], 2.0 * a[2], a[1]);
    if qa != 0.0 {
        let disc = qb * qb - 4.0 * qa * qc;
        if disc >= 0.0 {
            let r = libm::sqrt(disc);
            xs.push((-qb - r) / (2.0 * qa));
            xs.push((-qb + r) / (2.0 * qa));
        }
    } else if qb != 0.0 {
        xs.push(-qc / qb);
    }
    xs.iter()
        .filter(|x| x.is_finite() && (0.0..=x_max).contains(*x))
        .map(|&x| p(x))
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), v| {
            (lo.min(v), hi.max(v))
        })
}

/// Power-basis coefficients `a` (in x − xs\[0\]) of the polynomial through the points
/// (Newton divided differences, then expansion). 1 to 4 points with distinct xs.
fn interpolate(xs: &[f64], ys: &[f64]) -> [f64; 4] {
    let m = xs.len();
    let mut c = [0.0_f64; 4];
    c[..m].copy_from_slice(&ys[..m]);
    for level in 1..m {
        for i in (level..m).rev() {
            c[i] = (c[i] - c[i - 1]) / (xs[i] - xs[i - level]);
        }
    }
    // p(x) = c0 + c1·x + c2·x(x − x1) + c3·x(x − x1)(x − x2), x relative to xs[0].
    let x1 = if m > 1 { xs[1] - xs[0] } else { 0.0 };
    let x2 = if m > 2 { xs[2] - xs[0] } else { 0.0 };
    [
        c[0],
        c[1] - c[2] * x1 + c[3] * (x1 * x2),
        c[2] - c[3] * (x1 + x2),
        c[3],
    ]
}

/// Forward differences at j = 0 of the power-basis cubic `a` with unit step.
fn forward_differences(a: [f64; 4]) -> [f64; 4] {
    [
        a[0],
        a[1] + a[2] + a[3],
        2.0 * a[2] + 6.0 * a[3],
        6.0 * a[3],
    ]
}

/// Local indices for fractions `u` of an n-sample interval, deduplicated, ascending.
fn positions(us: &[f64], n: u64) -> Vec<u64> {
    let span = (n - 1) as f64;
    let mut v: Vec<u64> = us
        .iter()
        .map(|u| (libm::round(u * span) as u64).min(n - 1))
        .collect();
    v.sort_unstable();
    v.dedup();
    v
}

/// The signed difference of two phases, radians.
fn phase_error_rad(a: Phase, b: Phase) -> f64 {
    let d = a.wrapping_sub(b).to_bits() as i64;
    (d as f64 / TWO_POW_64 * 2.0 * core::f64::consts::PI).abs()
}

/// Fits one interval `[first, first + n)` of the block starting at sample `n0`.
fn fit_interval(
    model: &PathModel<'_>,
    clock: &dyn ClockMap,
    n0: SampleIndex,
    first: u64,
    n: u64,
) -> Result<Segment, FitError> {
    let at = |j: u64| -> Result<Eval, FitError> { model.eval(clock, n0.checked_add(first + j)?) };
    let nodes = positions(&NODES, n);
    let evals = nodes
        .iter()
        .map(|&j| at(j))
        .collect::<Result<Vec<_>, _>>()?;
    let e0 = evals[0];
    let f_s = model.target_rate_hz();
    let f_e = model.carrier_hz();
    let xs: Vec<f64> = nodes.iter().map(|&j| j as f64).collect();
    let dt = |e: &Eval| (e.t - e0.t) as f64 / TWO_POW_64 / 1e12;

    // Delay. The fitted cubic must stay inside ±(2³¹ − 1) samples over the whole
    // interval, interior extrema included (review fix: node checks alone missed an
    // interior crossing that wraps in Q32.32): the contracts §1.3 margin with room to
    // spare.
    let ys: Vec<f64> = evals.iter().map(|e| (e.tau - e0.tau) * f_s).collect();
    let mut a = interpolate(&xs, &ys);
    a[0] = e0.tau * f_s;
    let (lo, hi) = cubic_range(a, (n - 1) as f64);
    for v in [lo, hi] {
        if !v.is_finite() || v.abs() >= 2_147_483_647.0 {
            return Err(FitError::DelayRange { samples: v });
        }
    }
    let mut d = forward_differences(a);
    d[0] = e0.tau * f_s;
    let delay = DelayPoly::from_f64_differences(d, n)?;

    // Phase: the variation in f64, D0 exact.
    let ys: Vec<f64> = evals
        .iter()
        .map(|e| model.offset_hz * dt(e) - f_e * (e.tau - e0.tau))
        .collect();
    let mut d = forward_differences(interpolate(&xs, &ys));
    d[0] = 0.0;
    let mut regs = PhasePoly::from_f64_differences(d, n)?.registers();
    regs[0] = iqsim_core::wide::U128::new(e0.phase.to_bits(), 0);
    let phase = PhasePoly::from_registers(regs, n)?;

    // Coordinate: the variation; the exact start is kept apart.
    let coord = match model.coord {
        Some(c) => {
            let f = c.rate_hz();
            let ys: Vec<f64> = evals
                .iter()
                .map(|e| f * (dt(e) - (e.tau - e0.tau)))
                .collect();
            let mut d = forward_differences(interpolate(&xs, &ys));
            d[0] = 0.0;
            Some(Coord {
                start: e0.coord,
                variation: DelayPoly::from_f64_differences(d, n)?,
            })
        }
        None => None,
    };

    // Amplitude: f32 cubic in the kernel's x.
    let inv_len = (1.0 / n as f64) as f32;
    let amplitude = if evals.iter().all(|e| e.amplitude == e0.amplitude) {
        AmplitudePoly::constant(e0.amplitude as f32, n)
    } else {
        let xa: Vec<f64> = nodes
            .iter()
            .map(|&j| f64::from((j as f32) * inv_len))
            .collect();
        let ya: Vec<f64> = evals.iter().map(|e| e.amplitude).collect();
        let a = interpolate(&xa, &ya);
        AmplitudePoly {
            c: a.map(|v| v as f32),
            inv_len,
        }
    };

    // The response offset: pulses and continuous sources are filtered at the
    // interval's Doppler-shifted offset, on a grid of half the tolerance so blocks of a
    // slowly changing Doppler share tables and kernels (review fix: they were filtered
    // at the nominal offset whatever the Doppler).
    let (f_lo, f_hi) = evals
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), e| {
            (l.min(e.f_inst), h.max(e.f_inst))
        });
    let offset_uhz = match model.filtered {
        Some((_, tol)) => {
            let grid = libm::floor(tol * 0.5e6).max(1.0);
            let doppler_uhz = (0.5 * (f_lo + f_hi) - model.offset_hz) * 1e6;
            let steps = libm::round(doppler_uhz / grid);
            model.offset_uhz + iqsim_core::f64_to_i64_ties_even(steps * grid)?
        }
        None => model.offset_uhz,
    };
    let key_hz = offset_uhz as f64 * 1e-6;

    // Residuals against the exact model.
    let mut residuals = Residuals::default();
    let mut prev_f: Option<f64> = None;
    for j in positions(&CHECKS, n) {
        let e = at(j)?;
        let mut r = Residuals {
            delay_samples: (delay.eval(j)?.to_samples_f64() - e.tau * f_s).abs(),
            phase_rad: phase_error_rad(phase.eval(j)?, e.phase),
            ..Residuals::default()
        };
        if let (Some(c), Some(m)) = (&coord, model.coord) {
            let units = (c.at(j)? - e.coord) as f64 / TWO_POW_32;
            r.coord_samples = (units * f_s / m.rate_hz()).abs();
        }
        // Relative above the null floor; below it the absolute bound 10⁻⁵·V_ref
        // (contracts §2), expressed as the relative budget of an amplitude at
        // 10⁻⁵/budget·V_ref (review fix: a deep-stopband value that underflows f32 is
        // a 100 % relative error but no error at all against V_ref).
        let scale = e.amplitude.max(1e-5 / AMPLITUDE_BUDGET_REL * e.v_ref);
        if scale > 0.0 {
            r.amplitude_rel = (f64::from(amplitude.eval(j)) - e.amplitude).abs() / scale;
        }
        match (model.filtered, model.fold) {
            (Some((_, tol)), _) => r.response = (e.f_inst - key_hz).abs() / tol,
            // CW: |H(f_inst)| must be resolved by the checks: f_inst may move at most
            // σ/8 between adjacent check points (review fix: a sweep across a sharp
            // transition could vanish between checks).
            (None, Some(h)) => {
                if let Some(p) = prev_f {
                    r.response = (e.f_inst - p).abs() / (0.125 * h.sigma_hz());
                }
                prev_f = Some(e.f_inst);
            }
            (None, None) => {}
        }
        residuals = residuals.max(r);
    }
    Ok(Segment {
        first,
        len: n,
        delay,
        phase,
        coord,
        amplitude,
        offset_uhz,
        residuals,
    })
}

/// Fits `[first, first + n)`, halving on a missed budget (module docs).
#[allow(clippy::too_many_arguments)]
fn fit_refined(
    model: &PathModel<'_>,
    clock: &dyn ClockMap,
    n0: SampleIndex,
    first: u64,
    n: u64,
    depth: u32,
    out: &mut Vec<Segment>,
    splits: &mut Vec<Split>,
) -> Result<(), FitError> {
    let seg = fit_interval(model, clock, n0, first, n)?;
    if seg.residuals.accepted() {
        out.push(seg);
        return Ok(());
    }
    if depth >= MAX_REFINE_DEPTH || n < 2 {
        return Err(FitError::Budget {
            first,
            len: n,
            residuals: seg.residuals,
        });
    }
    let half = n / 2;
    splits.push(Split {
        at: first + half,
        reason: SplitReason::Refine,
        residuals: Some(seg.residuals),
        interval: Some((first, n)),
    });
    fit_refined(model, clock, n0, first, half, depth + 1, out, splits)?;
    fit_refined(
        model,
        clock,
        n0,
        first + half,
        n - half,
        depth + 1,
        out,
        splits,
    )
}

/// The fitted segments of one path over one block, in sample order, and its splits.
#[derive(Clone, Debug, PartialEq)]
pub struct PathFit {
    /// Segments tiling `[0, len)`.
    pub segments: Vec<Segment>,
    /// Every split, knots first, then refinements in the order they happened.
    pub splits: Vec<Split>,
}

impl PathFit {
    /// Worst residuals over every segment.
    pub fn residuals(&self) -> Residuals {
        self.segments
            .iter()
            .fold(Residuals::default(), |r, s| r.max(s.residuals))
    }
}

/// Fits one path over the block of `len` samples starting at channel sample `n0`
/// (module docs).
pub fn fit_path(
    model: &PathModel<'_>,
    clock: &dyn ClockMap,
    n0: SampleIndex,
    len: u64,
) -> Result<PathFit, FitError> {
    let mut cuts: Vec<(u64, SplitReason)> = Vec::new();
    let local = |t: AcqInstant| -> Result<Option<u64>, FitError> {
        let idx = clock.index_at_or_after(t)?.get();
        Ok(idx.checked_sub(n0.get()).filter(|&k| k > 0 && k < len))
    };
    let last = n0.checked_add(len - 1)?;
    let (t_lo, t_hi) = (clock.instant(n0)?, clock.instant(last)?);
    let geo = model.geometry;
    let hi_ps =
        |t: AcqInstant| TruthTime::from_ps(t.ps().as_ps().saturating_add(i64::from(t.frac() > 0)));
    for k in geo.rx.knots(t_lo.ps(), hi_ps(t_hi)) {
        if let Some(c) = local(AcqInstant::from_truth_time(k))? {
            cuts.push((c, SplitReason::ReceiverKnot));
        }
    }
    let (e_lo, e_hi) = (geo.at(t_lo)?.tx.instant, geo.at(t_hi)?.tx.instant);
    for k in geo.tx.knots(e_lo.ps(), hi_ps(e_hi)) {
        let lt = light_time::solve_emit(&*geo.rx, &*geo.tx, AcqInstant::from_truth_time(k))?;
        if let Some(c) = local(lt.rx.instant)? {
            cuts.push((c, SplitReason::TransmitterKnot));
        }
    }
    cuts.sort_by_key(|&(c, r)| (c, r as u8));
    cuts.dedup_by_key(|c| c.0);
    let mut splits: Vec<Split> = cuts
        .iter()
        .map(|&(at, reason)| Split {
            at,
            reason,
            residuals: None,
            interval: None,
        })
        .collect();
    let mut bounds: Vec<u64> = vec![0];
    bounds.extend(cuts.iter().map(|c| c.0));
    bounds.push(len);
    let mut segments = Vec::new();
    for w in bounds.windows(2) {
        fit_refined(
            model,
            clock,
            n0,
            w[0],
            w[1] - w[0],
            0,
            &mut segments,
            &mut splits,
        )?;
    }
    Ok(PathFit { segments, splits })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn carrier_phase_is_the_exact_fraction() {
        // 1 GHz × 1.5 ns = 1.5 cycles → ½ cycle.
        let p = carrier_phase(1_000_000_000_000_000, 1.5e-9).unwrap();
        let half = 1_u64 << 63;
        assert!(p.to_bits().abs_diff(half) < 1 << 20, "{:x}", p.to_bits());
        // A huge product: 31 GHz over 1.3 s (4·10¹⁰ cycles) keeps its fraction: compare
        // with an independent decomposition f·τ = f·τ_ps·10⁻¹² in exact decimal steps.
        let tau = 1.3;
        let p = carrier_phase(31_000_000_000_000_000, tau).unwrap();
        let d = light_time::seconds_to_q64_ps(tau).unwrap();
        // 31·10⁹ Hz · d·2⁻⁶⁴·10⁻¹² s = 31·10⁻³·d·2⁻⁶⁴ cycles; frac in 2⁻⁶⁴ units is
        // (31·d / 1000) mod 2⁶⁴.
        let want = ((31 * d as u128) / 1000) as u64;
        assert_eq!(p.to_bits(), want);
        assert_eq!(carrier_phase(5, 0.0).unwrap(), Phase::ZERO);
        assert!(carrier_phase(-1, 1.0).is_err());
    }

    #[test]
    fn intervals_far_apart_do_not_overflow() {
        // Review P1: +100 days against an origin at −100 days overflowed a single
        // Q64.64 i128 difference. At 1 S/s that is 17 280 000 samples.
        let day_ps: i128 = 86_400_000_000_000_000;
        let rate = SampleRate::from_hz(1).unwrap();
        let w = interval_to_q32(200 * day_ps, 1 << 63, rate);
        assert_eq!(w >> 32, 17_280_000);
        // The fraction: ½ ps at 1 S/s is 2⁻¹·10⁻¹² samples, below one Q32.32 unit.
        assert_eq!(w & 0xFFFF_FFFF, 0);
        // A fraction that is visible: ½ ps at 1 GS/s = 5·10⁻⁴ samples.
        let g = SampleRate::from_hz(1_000_000_000).unwrap();
        let w = interval_to_q32(0, 1 << 63, g);
        assert_eq!(w, (0.0005_f64 * 4_294_967_296.0) as i128);
        assert_eq!(interval_to_q32(-1, 0, g), -(4_294_967_296_i128 / 1000) - 1);
    }

    #[test]
    fn cubic_range_finds_interior_extrema() {
        // Review P1: 2³¹ + 0.001 − 100(u − 0.4)² over u ∈ [0, 1] (x = u·(n − 1)) peaks
        // inside, above every node value.
        let n = 4_000_000.0_f64;
        let s = n - 1.0;
        let a = [2_147_483_648.001 - 16.0, 80.0 / s, -100.0 / (s * s), 0.0];
        let (lo, hi) = cubic_range(a, s);
        assert!((hi - 2_147_483_648.001).abs() < 1e-3, "{hi}");
        assert!((lo - (2_147_483_648.001 - 36.0)).abs() < 1e-3, "{lo}");
        // A cubic with both critical points inside.
        let (lo, hi) = cubic_range([0.0, -3.0, 0.0, 1.0], 3.0);
        assert_eq!((lo, hi), (-2.0, 18.0));
        assert_eq!(cubic_range([5.0, 0.0, 0.0, 0.0], 0.0), (5.0, 5.0));
    }

    #[test]
    fn interpolation_reproduces_cubics_and_forward_differences_match() {
        let p = |x: f64| 3.0 - 0.5 * x + 1e-3 * x * x + 2e-7 * x * x * x;
        for xs in [
            vec![0.0],
            vec![0.0, 7.0],
            vec![0.0, 2.0, 9.0],
            vec![0.0, 250.0, 750.0, 999.0],
        ] {
            let ys: Vec<f64> = xs.iter().map(|&x| p(x)).collect();
            let a = interpolate(&xs, &ys);
            let q = |x: f64| a[0] + x * (a[1] + x * (a[2] + x * a[3]));
            for &x in &xs {
                assert!((q(x) - p(x)).abs() < 1e-9);
            }
            if xs.len() == 4 {
                let d = forward_differences(a);
                assert!((d[1] - (p(1.0) - p(0.0))).abs() < 1e-12);
                assert!((d[2] - (p(2.0) - 2.0 * p(1.0) + p(0.0))).abs() < 1e-12);
                assert!((d[3] - 6.0 * 2e-7).abs() < 1e-15);
            }
        }
        assert_eq!(positions(&NODES, 1), vec![0]);
        assert_eq!(positions(&NODES, 2), vec![0, 1]);
        assert_eq!(positions(&NODES, 1001), vec![0, 250, 750, 1000]);
    }
}
