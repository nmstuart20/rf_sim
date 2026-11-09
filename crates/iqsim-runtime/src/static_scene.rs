//! Temporary static-geometry work-list builder (prompt 5).
//!
//! It builds [`BlockWorkList`]s for fixed emitters and one fixed receiver. It stands in
//! for the phase 3 geometry engine and keeps that engine's interfaces: every
//! acquisition instant comes from a [`ClockMap`], every response effect comes from
//! the target's [`TapResponse`], and work lists are built per tap.
//!
//! Model, per emitter → receiver link:
//!
//! - Isotropic antennas and free-space loss. The voltage at the target's reference
//!   plane is `|v| = √(P_t·R₀)·λ/(4πR)` (Friis with unit gains, 50 Ω RMS).
//! - Constant delay `τ = R/c`. Straight line in ECEF, with no Sagnac or light-time
//!   iteration. Phase 3 replaces this.
//! - Carrier phase `φ(k) = f_off·t_acq(k) − f_c·τ`. The block-start term uses the exact
//!   [`OffsetFrequency::phase_at`]. The slope is f_off times the clock map's sample
//!   spacing over the block.
//! - Waveform coordinate `w(k) = F_s·(t_acq(k) − τ − t_emission)` for pulse trains.
//! - CW: |H(f_off)| folds into the amplitude and −f_off·τ_g into the phase. Pulses
//!   carry the response in their edge tables.
//! - Continuous sources (IQ files, contracts §5.6): the combined kernel K (source
//!   reconstruction φ ⊛ the target response shifted to the carrier offset) is built
//!   once per (source, target, offset). Per block, the native coordinate
//!   `w(k) = F_n·(t_acq(k) − start) − F_n·τ` is computed exactly at k = 0 (Q64.64
//!   instant × rational rate in i128, then τ rounded once to Q32.32), one stage-1
//!   window is filtered per group of links that share a source, kernel and start, and
//!   every path's coordinate is rebased to its window's z\[0\], so coordinates stay
//!   small for any run length or loop count.
//!
//! Links with |f_c·τ| above 2⁴⁰ cycles are rejected, because the f64 product would no
//! longer resolve the fractional cycle within the phase budget.
//!
//! The clock map must be linear over each block, which the ideal map is to 2⁻⁶⁴ ps.
//! The builder checks the midpoint residual and rejects a non-linear map, because a
//! linear fit is all a static builder does.

use std::collections::BTreeMap;
use std::fmt;
use std::sync::{Arc, Mutex};

use iqsim_core::amplitude::rms_volts;
use iqsim_core::wide::U128;
use iqsim_core::{
    AmplitudeError, BlockSpec, ClockMap, DelayPoly, Headroom, NumericError, OffsetFrequency, Phase,
    PhasePoly, ResponseError, ResponseSpec, SampleIndex, SampleRate, TapKind, TapResponse,
    TargetScale, TimeError, TruthTime,
};
use iqsim_geo::{Ecef, SPEED_OF_LIGHT_M_S};
use iqsim_render::FARROW_LEBESGUE;
use iqsim_waveform::continuous::source_response;
use iqsim_waveform::{
    ContinuousError, ContinuousSource, EdgeError, EdgeShape, EdgeTableKey, EndPolicy,
    ResampleKernel, ResampleKey, SourceTimeline, build_edge_table,
};
use iqsim_worklist::{
    AmplitudePoly, BlockWorkList, COORD_NATIVE_SAMPLES, COORD_TARGET_SAMPLES, ContinuousEmission,
    CwEmission, EMISSION_CONTINUOUS, EMISSION_CW, EMISSION_PULSE_TRAIN, EdgeTableArena,
    EmissionRef, FARROW_AFTER, FARROW_BEFORE, MAX_COORD_BITS, MAX_COORD_SAMPLES, MAX_WINDOW_LEN,
    PathSegment, PolyRegisters, Pulse, PulseTrain, SourceWindowArena, TargetDesc, WorkListError,
    target_desc,
};

/// A continuous source shared by every target and block that reads it. Reads are pure
/// functions of the sample index, so sharing one reader (and its cache) cannot change
/// any value.
#[derive(Clone)]
pub struct SharedSource(pub Arc<Mutex<Box<dyn ContinuousSource + Send>>>);

impl SharedSource {
    /// Wraps a source.
    pub fn new(source: impl ContinuousSource + Send + 'static) -> Self {
        Self(Arc::new(Mutex::new(Box::new(source))))
    }
}

impl fmt::Debug for SharedSource {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("SharedSource(..)")
    }
}

impl PartialEq for SharedSource {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }
}

/// A continuous emission: native samples played from `start` on the emitter carrier.
#[derive(Clone, Debug, PartialEq)]
pub struct ContinuousSpec {
    /// The samples.
    pub source: SharedSource,
    /// Length and end policy.
    pub timeline: SourceTimeline,
    /// Native rate.
    pub rate: SampleRate,
    /// Truth time at which sample 0 is emitted.
    pub start: TruthTime,
    /// Source reconstruction response φ.
    pub response: ResponseSpec,
    /// Level factor applied to the normalized samples (1 for full scale, 1/RMS for RMS
    /// normalization).
    pub scale: f64,
    /// max |x| of the normalized samples (for the pre-filter cull bound).
    pub peak: f64,
}

/// One pulse, in emitted time relative to its train's emission time.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PulseSpec {
    /// Start of the rising edge, ps.
    pub start_ps: i64,
    /// Rise duration (linear ramp), ps.
    pub rise_ps: i64,
    /// Start of the falling edge, ps (≥ `start_ps`).
    pub fall_start_ps: i64,
    /// Fall duration (linear ramp), ps.
    pub fall_ps: i64,
    /// Relative amplitude.
    pub amplitude: f32,
}

/// What an emitter transmits.
#[derive(Clone, Debug, PartialEq)]
pub enum StaticWaveform {
    /// An unmodulated carrier, always on.
    Cw,
    /// A continuous sampled source (an IQ file).
    Continuous(ContinuousSpec),
    /// Pulses relative to one emission time.
    Pulses {
        /// Emission reference time (the waveform-coordinate origin).
        emission: TruthTime,
        /// Pulses, sorted by start.
        pulses: Vec<PulseSpec>,
    },
}

/// A fixed emitter.
#[derive(Clone, Debug, PartialEq)]
pub struct StaticEmitter {
    /// Position.
    pub position: Ecef,
    /// Carrier frequency, µHz.
    pub carrier_uhz: i64,
    /// Transmit power (pulse peak or CW), W.
    pub power_w: f64,
    /// Waveform.
    pub waveform: StaticWaveform,
}

/// A fixed receiver and its emitters.
#[derive(Clone, Debug, PartialEq)]
pub struct StaticScene {
    /// Receiver position.
    pub receiver: Ecef,
    /// Emitters, in a fixed order.
    pub emitters: Vec<StaticEmitter>,
}

/// One render target: a tap of a receiver channel.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StaticTarget {
    /// Antenna port (B_out) or virtual receiver (H_rx).
    pub kind: TapKind,
    /// Channel centre, µHz.
    pub center_uhz: i64,
    /// Output sample rate.
    pub rate: SampleRate,
    /// The tap's pre-sampling response.
    pub response: ResponseSpec,
    /// Full-scale voltage of this target's accumulator scale.
    pub full_scale_volts: f64,
    /// Paths whose peak contribution is below this many LSB are culled. The lossless
    /// default is [`iqsim_core::amplitude::LOSSLESS_CULL_LSB`].
    pub cull_threshold_lsb: f64,
}

/// A failure while building a scene or a block.
#[derive(Clone, Debug, PartialEq)]
pub enum SceneError {
    /// The tap response is invalid.
    Response(ResponseError),
    /// The scale or the headroom proof failed.
    Amplitude(AmplitudeError),
    /// An edge table could not be built.
    Edge(EdgeError),
    /// A work list was malformed.
    WorkList(WorkListError),
    /// A fixed-point conversion failed.
    Numeric(NumericError),
    /// A clock or index operation failed.
    Time(TimeError),
    /// The clock map is not linear over the block.
    NonLinearClock,
    /// A waveform coordinate or pulse offset is outside ±2³⁰ samples.
    CoordinateRange,
    /// An emitter definition is invalid (power, frequency, pulse order).
    Emitter {
        /// Emitter index.
        index: usize,
    },
    /// A continuous source could not be prepared or filtered.
    Continuous(ContinuousError),
    /// A block's stage-1 window would exceed `MAX_WINDOW_LEN` points (use shorter
    /// blocks).
    WindowTooLong {
        /// Emitter index.
        index: usize,
        /// Points the window would need.
        points: u64,
    },
    /// A `stop` source does not cover the block: the receiver would hear past its last
    /// sample.
    SourceExhausted {
        /// Emitter index.
        index: usize,
        /// Block number.
        block: u64,
    },
}

impl fmt::Display for SceneError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Response(e) => e.fmt(f),
            Self::Amplitude(e) => e.fmt(f),
            Self::Edge(e) => e.fmt(f),
            Self::WorkList(e) => e.fmt(f),
            Self::Numeric(e) => e.fmt(f),
            Self::Time(e) => e.fmt(f),
            Self::NonLinearClock => {
                f.write_str("static builder needs a clock map linear over each block")
            }
            Self::CoordinateRange => {
                f.write_str("waveform coordinate outside ±2^30 samples (rebase the source)")
            }
            Self::Emitter { index } => write!(f, "emitter {index} is invalid"),
            Self::Continuous(e) => e.fmt(f),
            Self::WindowTooLong { index, points } => write!(
                f,
                "emitter {index}: a block needs a {points}-point source window, above the \
                 {MAX_WINDOW_LEN}-point limit; use shorter blocks"
            ),
            Self::SourceExhausted { index, block } => write!(
                f,
                "emitter {index}: its IQ file (end = \"stop\") ends before block {block}; \
                 use a longer recording, a shorter run, or end = \"zero_pad\" or \"loop\""
            ),
        }
    }
}

impl std::error::Error for SceneError {}

macro_rules! from_err {
    ($($t:ty => $v:ident),*) => {$(
        impl From<$t> for SceneError {
            fn from(e: $t) -> Self { Self::$v(e) }
        }
    )*};
}
from_err!(ResponseError => Response, AmplitudeError => Amplitude, EdgeError => Edge,
    WorkListError => WorkList, NumericError => Numeric, TimeError => Time,
    ContinuousError => Continuous);

/// A continuous link's prepared state.
struct ContinuousLink {
    spec: ContinuousSpec,
    /// Index of the shared source among the builder's distinct sources.
    source_id: usize,
    kernel: Arc<ResampleKernel>,
    rate_hz: f64,
    /// F_n·τ, Q32.32 native samples.
    tau_q: i64,
}

enum LinkWave {
    Cw,
    Continuous(Box<ContinuousLink>),
    Pulses {
        emission: TruthTime,
        pulses: Vec<Pulse>,
    },
}

struct Link {
    tau_s: f64,
    volts: f64,
    carrier_hz: f64,
    offset_hz: f64,
    offset: OffsetFrequency,
    wave: LinkWave,
}

/// Relative slack on the Farrow envelope bound for its f32 evaluation (a dozen
/// roundings of at most 2⁻²⁴ each, with margin).
const FARROW_ROUNDING_SLACK: f64 = 1e-5;

/// A continuous path's stage-1 result for one block.
struct PreparedPath {
    window: u32,
    coord: PolyRegisters,
    peak_env: f64,
}

/// Builds work lists for one target of a static scene.
pub struct StaticSceneBuilder<'c> {
    response: TapResponse,
    scale: TargetScale,
    desc: TargetDesc,
    rate_hz: f64,
    cull_threshold_lsb: f64,
    clock: &'c dyn ClockMap,
    blocks: BlockSpec,
    tables: Arc<EdgeTableArena>,
    links: Vec<Link>,
    /// Samples the run keeps in total, if known: a final partial block is checked
    /// against the `stop` policy only up to it.
    sample_limit: Option<u64>,
}

/// Largest |f_c·τ| in cycles the static builder accepts: f64 then resolves the
/// fractional cycle to 2⁴⁰⁻⁵³ = 2⁻¹³ cycle (7.7·10⁻⁴ rad), inside the 0.01 rad budget.
const MAX_PROPAGATION_CYCLES: f64 = 1_099_511_627_776.0;

/// Q32.32 samples at `rate` in a Q64.64 ps interval `dt`: ⌊dt·num/(den·10¹²·2⁶⁴)·2³²⌋
/// to within one unit (the low 32 bits of the fractional product are dropped).
///
/// dt = ps·2⁶⁴ + frac. `ps·num` < 2¹¹² and the remainder pieces stay below 2¹⁰⁵, so no
/// intermediate leaves i128.
fn q64_interval_to_q32(dt: i128, rate: SampleRate) -> i128 {
    let ps = dt >> 64;
    let frac = (dt & ((1_i128 << 64) - 1)) as u128;
    let num = i128::from(rate.num());
    let d = i128::from(rate.den()) * 1_000_000_000_000;
    let a = ps * num;
    let q = a.div_euclid(d);
    let r = a.rem_euclid(d);
    let f = ((frac * rate.num() as u128) >> 32) as i128;
    (q << 32) + ((r << 32) + f) / d
}

/// Q32.32 target samples for `ps` picoseconds at `rate`: ⌊ps·num·2³²/(den·10¹²)⌋.
///
/// The whole-sample quotient and the remainder are split before the shift, so no
/// intermediate exceeds i128: |ps·num| < 2¹¹¹, and the remainder times 2³² < 2¹⁰⁴.
fn ps_to_coord(ps: i64, rate: SampleRate) -> Result<i64, SceneError> {
    let n = i128::from(ps) * i128::from(rate.num());
    let d = i128::from(rate.den()) * 1_000_000_000_000;
    let whole = n.div_euclid(d);
    let rem = n.rem_euclid(d);
    if whole.abs() >= i128::from(MAX_COORD_SAMPLES) {
        return Err(SceneError::CoordinateRange);
    }
    let frac = (rem << 32) / d;
    // |whole| < 2³⁰ and 0 ≤ frac < 2³², so this fits i64 and stays inside ±2⁶².
    Ok(((whole << 32) + frac) as i64)
}

impl<'c> StaticSceneBuilder<'c> {
    /// Validates the target, computes every link and builds every edge table the
    /// scene needs (before any block renders).
    pub fn new(
        scene: &StaticScene,
        target: &StaticTarget,
        clock: &'c dyn ClockMap,
        blocks: BlockSpec,
    ) -> Result<Self, SceneError> {
        let response =
            TapResponse::new(target.kind, target.response, target.rate, target.center_uhz)?;
        let scale = TargetScale::new(target.full_scale_volts)?;
        if !(target.cull_threshold_lsb.is_finite() && target.cull_threshold_lsb >= 0.0) {
            return Err(SceneError::Amplitude(AmplitudeError::Bound));
        }
        let rate = target.rate;
        let mut arena = EdgeTableArena::default();
        let mut index: BTreeMap<EdgeTableKey, u32> = BTreeMap::new();
        let mut links = Vec::with_capacity(scene.emitters.len());
        let mut kernels: BTreeMap<ResampleKey, Arc<ResampleKernel>> = BTreeMap::new();
        let mut sources: Vec<SharedSource> = Vec::new();
        for (i, e) in scene.emitters.iter().enumerate() {
            if !(e.power_w.is_finite() && e.power_w >= 0.0) || e.carrier_uhz <= 0 {
                return Err(SceneError::Emitter { index: i });
            }
            let range = scene.receiver.distance(e.position);
            let carrier_hz = e.carrier_uhz as f64 * 1e-6;
            let lambda = SPEED_OF_LIGHT_M_S / carrier_hz;
            let volts = if range > 0.0 {
                rms_volts(e.power_w) * lambda / (4.0 * core::f64::consts::PI * range)
            } else {
                return Err(SceneError::Emitter { index: i });
            };
            let offset_uhz = e
                .carrier_uhz
                .checked_sub(target.center_uhz)
                .ok_or(SceneError::Emitter { index: i })?;
            let offset_hz = offset_uhz as f64 * 1e-6;
            let tau_s = range / SPEED_OF_LIGHT_M_S;
            let wave = match &e.waveform {
                StaticWaveform::Cw => LinkWave::Cw,
                StaticWaveform::Continuous(spec) => {
                    if !(spec.scale.is_finite() && spec.scale >= 0.0)
                        || !(spec.peak.is_finite() && spec.peak >= 0.0)
                    {
                        return Err(SceneError::Emitter { index: i });
                    }
                    let phi = source_response(spec.rate, spec.response)?;
                    let key = ResampleKey {
                        source: phi.id(),
                        target: response.id(),
                        offset_uhz,
                    };
                    let kernel = match kernels.get(&key) {
                        Some(k) => Arc::clone(k),
                        None => {
                            let k = Arc::new(ResampleKernel::build(&phi, &response, offset_uhz)?);
                            kernels.insert(key, Arc::clone(&k));
                            k
                        }
                    };
                    let source_id = match sources.iter().position(|s| *s == spec.source) {
                        Some(p) => p,
                        None => {
                            sources.push(spec.source.clone());
                            sources.len() - 1
                        }
                    };
                    let rate_hz = spec.rate.num() as f64 / spec.rate.den() as f64;
                    let tau_q =
                        iqsim_core::f64_to_i64_ties_even(tau_s * rate_hz * 4_294_967_296.0)?;
                    LinkWave::Continuous(Box::new(ContinuousLink {
                        spec: spec.clone(),
                        source_id,
                        kernel,
                        rate_hz,
                        tau_q,
                    }))
                }
                StaticWaveform::Pulses { emission, pulses } => {
                    let mut built = Vec::with_capacity(pulses.len());
                    for (j, p) in pulses.iter().enumerate() {
                        if p.fall_start_ps < p.start_ps
                            || !p.amplitude.is_finite()
                            || p.amplitude.abs() > iqsim_worklist::MAX_PULSE_AMPLITUDE
                            || (j > 0 && pulses[j - 1].start_ps > p.start_ps)
                        {
                            return Err(SceneError::Emitter { index: i });
                        }
                        let mut table = |rise_ps: i64| -> Result<u32, SceneError> {
                            let key = EdgeTableKey {
                                response: response.id(),
                                shape: EdgeShape::Linear { rise_ps },
                                offset_uhz,
                            };
                            if let Some(&t) = index.get(&key) {
                                return Ok(t);
                            }
                            let t = build_edge_table(&response, key.shape, offset_uhz)?;
                            let id = arena.push(t.u_lo, t.log2_res, &t.points)?;
                            index.insert(key, id);
                            Ok(id)
                        };
                        let rise_table = table(p.rise_ps)?;
                        let fall_table = table(p.fall_ps)?;
                        built.push(Pulse {
                            start: ps_to_coord(p.start_ps, rate)?,
                            fall_start: ps_to_coord(p.fall_start_ps, rate)?,
                            amplitude: p.amplitude,
                            rise_table,
                            fall_table,
                            _pad: 0,
                        });
                    }
                    LinkWave::Pulses {
                        emission: *emission,
                        pulses: built,
                    }
                }
            };
            if carrier_hz * tau_s > MAX_PROPAGATION_CYCLES {
                return Err(SceneError::Emitter { index: i });
            }
            links.push(Link {
                tau_s,
                volts,
                carrier_hz,
                offset_hz,
                offset: OffsetFrequency::from_uhz(offset_uhz)?,
                wave,
            });
        }
        Ok(Self {
            desc: target_desc(target.kind, response.id(), scale.lsb_per_volt()),
            response,
            scale,
            rate_hz: rate.num() as f64 / rate.den() as f64,
            cull_threshold_lsb: target.cull_threshold_lsb,
            clock,
            blocks,
            tables: Arc::new(arena),
            links,
            sample_limit: None,
        })
    }

    /// Declares how many samples the run keeps, so the `stop` check ignores the
    /// rendered-but-discarded tail of the final block.
    pub fn with_sample_limit(mut self, total: u64) -> Self {
        self.sample_limit = Some(total);
        self
    }

    /// The target's response.
    pub fn response(&self) -> &TapResponse {
        &self.response
    }

    /// The target's scale.
    pub fn scale(&self) -> &TargetScale {
        &self.scale
    }

    /// Builds the work list of block `b`: culls paths below the threshold, proves
    /// headroom and validates the result.
    pub fn build_block(&self, b: u64) -> Result<BlockWorkList, SceneError> {
        let block = self.blocks.block(b)?;
        let n0 = block.start();
        let len = block.len();
        let t0 = self.clock.instant(n0)?.to_q64_64();
        let q64 = 18_446_744_073_709_551_616.0_f64;
        // Sample spacing over the block, ps, from the clock map.
        let spacing_ps = if len > 1 {
            let t_end = self.clock.instant(n0.checked_add(len - 1)?)?.to_q64_64();
            let mid_k = (len - 1) / 2;
            let t_mid = self.clock.instant(n0.checked_add(mid_k)?)?.to_q64_64();
            let spacing = (t_end - t0) as f64 / q64 / (len - 1) as f64;
            let predicted = spacing * mid_k as f64;
            let actual = (t_mid - t0) as f64 / q64;
            if (actual - predicted).abs() > 1e-6 * spacing {
                return Err(SceneError::NonLinearClock);
            }
            spacing
        } else {
            1e12 / self.rate_hz
        };
        let spacing_s = spacing_ps * 1e-12;
        let samples_per_ps = self.rate_hz * 1e-12;

        let (windows, prepared) = self.continuous_windows(b, t0, spacing_s, len)?;
        let mut work = BlockWorkList {
            target: self.desc,
            block: b,
            first_sample: n0.get(),
            len,
            paths: Vec::new(),
            cw: Vec::new(),
            trains: Vec::new(),
            pulses: Vec::new(),
            continuous: Vec::new(),
            tables: Arc::clone(&self.tables),
            windows: Arc::new(windows),
        };
        let mut headroom = Headroom::new();
        let zero = PolyRegisters::from_poly(&DelayPoly::from_f64_differences([0.0; 4], len)?);
        for (source, link) in self.links.iter().enumerate() {
            let tau_samples = link.tau_s * self.rate_hz;
            let delay = DelayPoly::from_f64_differences([tau_samples, 0.0, 0.0, 0.0], len)?;
            let mut phase_extra = -(link.carrier_hz * link.tau_s);
            let mut window = None;
            let (amp, peak_env, coord, emission) = match &link.wave {
                LinkWave::Cw => {
                    let gain = self.response.magnitude(link.offset_hz);
                    phase_extra += self.response.phase_cycles(link.offset_hz);
                    let amp = link.volts * gain;
                    (amp, 1.0, zero, None)
                }
                LinkWave::Continuous(c) => {
                    // K carries the response, including its group delay, and stage 1
                    // applies the level factor in f64 (review fix: a large 1/RMS factor
                    // combined with tiny samples underflowed in f32), so the amplitude
                    // is Friis only.
                    let Some(p) = prepared.get(&source) else {
                        continue;
                    };
                    debug_assert!(c.spec.scale.is_finite());
                    window = Some(p.window);
                    (link.volts, p.peak_env, p.coord, None)
                }
                LinkWave::Pulses { emission, pulses } => {
                    let dt = t0 - i128::from(emission.as_ps()) * (1_i128 << 64);
                    let w0 = (dt as f64 / q64) * samples_per_ps - tau_samples;
                    let slope = spacing_ps * samples_per_ps;
                    let w_end = w0 + slope * (len - 1) as f64;
                    let limit = MAX_COORD_SAMPLES as f64;
                    // Clamped to ±2³⁰ samples (exact in Q32.32) only for the activity test.
                    let to_q = |x: f64| {
                        iqsim_core::f64_to_i64_ties_even(x.clamp(-limit, limit) * 4_294_967_296.0)
                    };
                    let (lo, hi) = (to_q(w0)?, to_q(w_end)?);
                    let tables = &*self.tables;
                    let active: Vec<&Pulse> = pulses
                        .iter()
                        .filter(|p| {
                            let (on, off) = tables.pulse_extent(p);
                            on <= hi && off >= lo
                        })
                        .collect();
                    if active.is_empty() {
                        continue;
                    }
                    if !(w0.abs() < limit && w_end.abs() < limit) {
                        return Err(SceneError::CoordinateRange);
                    }
                    // Peak envelope: the largest sum of A·(peak_rise + peak_fall) over
                    // pulses whose extents overlap one point (a maximum is reached at
                    // some pulse's start).
                    let extent = |p: &Pulse| {
                        let (on, off) = tables.pulse_extent(p);
                        let pk = f64::from(tables.headers[p.rise_table as usize].peak)
                            + f64::from(tables.headers[p.fall_table as usize].peak);
                        (on, off, f64::from(p.amplitude.abs()) * pk)
                    };
                    let mut peak = 0.0_f64;
                    for a in &active {
                        let (on_a, _, _) = extent(a);
                        let s: f64 = active
                            .iter()
                            .map(|p| extent(p))
                            .filter(|&(on, off, _)| on <= on_a && off >= on_a)
                            .map(|(_, _, v)| v)
                            .sum();
                        peak = peak.max(s);
                    }
                    let coord = DelayPoly::from_f64_differences([w0, slope, 0.0, 0.0], len)?;
                    (
                        link.volts,
                        peak,
                        PolyRegisters::from_poly(&coord),
                        Some(active),
                    )
                }
            };
            let peak_lsb = self.scale.peak_lsb(amp, peak_env)?;
            if peak_lsb < self.cull_threshold_lsb {
                continue;
            }
            let amplitude = AmplitudePoly::constant(amp as f32, len);
            // The envelope bound is taken as at least 1, so the proof also covers the
            // kernel's intermediate m = a·lsb_per_volt (contracts §2), not only the
            // final product.
            headroom.add(self.scale.peak_lsb(amplitude.bound(), peak_env.max(1.0))?)?;

            // Phase: exact offset phase at the block start, plus the static terms.
            let p0 = link
                .offset
                .phase_at(iqsim_core::AcqInstant::from_q64_64(t0))
                .wrapping_add(Phase::from_cycles_f64(phase_extra)?);
            let d1 = link.offset_hz * spacing_s;
            let slope_poly = PhasePoly::from_f64_differences([0.0, d1, 0.0, 0.0], len)?;
            let mut regs = slope_poly.registers();
            regs[0] = U128::new(p0.to_bits(), 0);
            let phase = PhasePoly::from_registers(regs, len)?;

            let emission = match emission {
                None if window.is_some() => {
                    work.continuous.push(ContinuousEmission {
                        window: window.unwrap_or_default(),
                        _pad: 0,
                    });
                    EmissionRef {
                        kind: EMISSION_CONTINUOUS,
                        index: (work.continuous.len() - 1) as u32,
                    }
                }
                None => {
                    work.cw.push(CwEmission {
                        envelope: [1.0, 0.0],
                    });
                    EmissionRef {
                        kind: EMISSION_CW,
                        index: (work.cw.len() - 1) as u32,
                    }
                }
                Some(active) => {
                    let first = work.pulses.len() as u32;
                    work.pulses.extend(active.iter().map(|p| **p));
                    work.trains.push(PulseTrain {
                        first_pulse: first,
                        pulse_count: active.len() as u32,
                    });
                    EmissionRef {
                        kind: EMISSION_PULSE_TRAIN,
                        index: (work.trains.len() - 1) as u32,
                    }
                }
            };
            work.paths.push(PathSegment {
                delay: PolyRegisters::from_poly(&delay),
                coord,
                phase: PolyRegisters::from_poly(&phase),
                block: b,
                amplitude,
                target: 0,
                coord_units: if window.is_some() {
                    COORD_NATIVE_SAMPLES
                } else {
                    COORD_TARGET_SAMPLES
                },
                coord_source: source as u32,
                emission,
            });
        }
        work.validate()?;
        Ok(work)
    }

    /// Stage 1 for block `b`: the windows of every continuous link that is not culled,
    /// and per link index its window, rebased coordinate and envelope bound.
    ///
    /// Links that share a source, kernel and start share one window covering all of
    /// their coordinates.
    fn continuous_windows(
        &self,
        b: u64,
        t0: i128,
        spacing_s: f64,
        len: u64,
    ) -> Result<(SourceWindowArena, BTreeMap<usize, PreparedPath>), SceneError> {
        struct Pending {
            link: usize,
            /// Absolute native coordinate at k = 0, Q32.32 bits.
            w0: i128,
            /// Absolute native coordinate at k = len − 1, Q32.32 bits.
            w1: i128,
            coord: DelayPoly,
        }
        // Grouped by everything the window depends on: source, kernel, start, timeline
        // (review fix: the end policy was missing) and level factor.
        type Key = (usize, ResampleKey, i64, u64, EndPolicy, u64);
        let mut groups: BTreeMap<Key, Vec<Pending>> = BTreeMap::new();
        // Samples of this block the run keeps (review fix: `stop` checked the discarded
        // tail of a final partial block).
        let keep = match self.sample_limit {
            Some(total) => total.saturating_sub(b.saturating_mul(len)).min(len),
            None => len,
        };
        for (i, link) in self.links.iter().enumerate() {
            let LinkWave::Continuous(c) = &link.wave else {
                continue;
            };
            let dt = t0 - i128::from(c.spec.start.as_ps()) * (1_i128 << 64);
            let w0 = q64_interval_to_q32(dt, c.spec.rate) - i128::from(c.tau_q);
            let slope = c.rate_hz * spacing_s;
            // The slope-only polynomial; its D0 is set once the window origin is known.
            let coord = DelayPoly::from_f64_differences([0.0, slope, 0.0, 0.0], len)?;
            let w1 = w0 + i128::from(coord.eval(len - 1)?.to_bits());
            let w_kept = w0 + i128::from(coord.eval(keep.max(1) - 1)?.to_bits());
            if c.spec.timeline.exhausted(w_kept) {
                return Err(SceneError::SourceExhausted { index: i, block: b });
            }
            // Lossless pre-filter cull: |z| ≤ l1_max·peak|x|, so nothing in this block
            // can reach the threshold.
            let bound =
                FARROW_LEBESGUE * c.kernel.l1_max() * c.spec.peak * (1.0 + FARROW_ROUNDING_SLACK);
            if self.scale.peak_lsb(link.volts * c.spec.scale, bound)? < self.cull_threshold_lsb {
                continue;
            }
            let t = c.spec.timeline;
            let key = (
                c.source_id,
                c.kernel.key(),
                c.spec.start.as_ps(),
                t.len,
                t.policy,
                c.spec.scale.to_bits(),
            );
            groups.entry(key).or_default().push(Pending {
                link: i,
                w0,
                w1,
                coord,
            });
        }
        let mut arena = SourceWindowArena::default();
        let mut out = BTreeMap::new();
        for pending in groups.values() {
            let LinkWave::Continuous(c) = &self.links[pending[0].link].wave else {
                unreachable!("grouped links are continuous");
            };
            let log2_up = c.kernel.log2_up();
            let shift = 32 - log2_up;
            let lo = pending.iter().map(|p| p.w0).min().unwrap_or_default();
            let hi = pending.iter().map(|p| p.w1).max().unwrap_or_default();
            // z cells, with the Farrow taps and one cell of margin on each side.
            let m_lo = (lo >> shift) - FARROW_BEFORE as i128 - 1;
            let m_hi = (hi >> shift) + FARROW_AFTER as i128 + 1;
            let m_lo = i64::try_from(m_lo).map_err(|_| SceneError::CoordinateRange)?;
            let z_len = usize::try_from(m_hi - i128::from(m_lo))
                .map_err(|_| SceneError::CoordinateRange)?;
            // Review fix: reject an oversized window before allocating or reading it.
            if z_len > MAX_WINDOW_LEN as usize {
                return Err(SceneError::WindowTooLong {
                    index: pending[0].link,
                    points: z_len as u64,
                });
            }
            let (first, last) = c.kernel.source_span(m_lo, z_len);
            if c.spec.timeline.all_zero(first, last) {
                continue;
            }
            let z = {
                let mut src = c
                    .spec
                    .source
                    .0
                    .lock()
                    .map_err(|_| ContinuousError::Source("source lock poisoned".to_owned()))?;
                c.kernel
                    .filter_window(&c.spec.timeline, src.as_mut(), m_lo, z_len, c.spec.scale)?
            };
            let window = arena.push(log2_up, &z)?;
            let peak = f64::from(arena.headers[window as usize].peak);
            let origin = i128::from(m_lo) << shift;
            for p in pending {
                let local = p.w0 - origin;
                if !(0..i128::from(MAX_COORD_BITS)).contains(&local) {
                    return Err(SceneError::CoordinateRange);
                }
                let mut regs = p.coord.registers();
                regs[0] = U128::new(local as u64, 0);
                let coord = DelayPoly::from_registers(regs, len)?;
                out.insert(
                    p.link,
                    PreparedPath {
                        window,
                        coord: PolyRegisters::from_poly(&coord),
                        peak_env: FARROW_LEBESGUE * peak * (1.0 + FARROW_ROUNDING_SLACK),
                    },
                );
            }
        }
        Ok((arena, out))
    }

    /// Channel sample index of block `b`'s first sample.
    pub fn block_start(&self, b: u64) -> Result<SampleIndex, SceneError> {
        Ok(self.blocks.block(b)?.start())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn coordinate_conversion_never_overflows_its_intermediates() {
        // Regression (prompt 5 review): rate (2⁴⁸ − 1)/2³² and 10¹⁶ ps is ≈ 6.55·10⁸
        // samples, in range, but ps·num·2³² needs 134 bits.
        let rate = SampleRate::from_ratio((1 << 48) - 1, 1 << 32).unwrap();
        let ps = 10_000_000_000_000_000_i64;
        let q = ps_to_coord(ps, rate).unwrap();
        // Exact reference: ⌊ps·num·2³²/(den·10¹²)⌋ by big-integer long division in u128
        // pieces: ps·num < 2¹¹¹, so split as (whole, rem) and shift the remainder.
        let n = u128::from(ps as u64) * u128::from(rate.num());
        let d = u128::from(rate.den()) * 1_000_000_000_000;
        let expect = ((n / d) << 32) + ((n % d) << 32) / d;
        assert_eq!(q as u128, expect);
        assert_eq!(q >> 32, 655_359_999);
        // −1 ps = −1.6·10⁻⁵ samples = −68 719.48 Q32.32 units, floored.
        assert_eq!(
            ps_to_coord(-1, SampleRate::from_hz(16_000_000).unwrap()).unwrap(),
            -68_720
        );
        assert_eq!(
            ps_to_coord(i64::MAX, rate),
            Err(SceneError::CoordinateRange)
        );
        assert_eq!(
            ps_to_coord(i64::MIN, rate),
            Err(SceneError::CoordinateRange)
        );
    }
}
