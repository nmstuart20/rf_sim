//! Work-list builder over the light-time geometry engine (prompt 9; prompt 5 built
//! its static predecessor).
//!
//! It builds [`BlockWorkList`]s for emitters and one receiver, each on its own
//! [`Trajectory`]. Every acquisition instant comes from a [`ClockMap`], every response
//! effect from the target's [`TapResponse`], and work lists are built per tap.
//!
//! Model, per emitter → receiver path ([`crate::geometry`] has the details):
//!
//! - Light time in ECI at every evaluated receive instant, with the transmitter at the
//!   retarded instant (Earth rotation, Sagnac and both motions included). Platform
//!   origins stand in for antenna phase centres until prompt 10, and antennas are
//!   isotropic: the voltage at the target's reference plane is
//!   `|v| = √(P_t·R₀)·λ/(4π·c·τ)` (Friis with unit gains, 50 Ω RMS).
//! - Delay, phase, waveform coordinate and amplitude polynomials are fitted per block
//!   on Chebyshev–Lobatto nodes and checked against the exact model, split at
//!   trajectory knots and refined until they meet the budgets ([`fit_path`]). A block
//!   can therefore hold several [`PathSegment`]s of one path.
//! - Carrier phase `φ(k) = f_off·t_acq(k) − f_e·τ(k)`, with exact modular integers at
//!   every segment start.
//! - CW: |H(f_inst)| folds into the amplitude and −f_inst·τ_g into the phase, at the
//!   instantaneous (Doppler-shifted) offset. Pulses carry the response in their edge
//!   tables, continuous sources in their resampling kernel.
//! - Continuous sources (IQ files, contracts §5.6): the combined kernel K (source
//!   reconstruction φ ⊛ the target response shifted to the carrier offset) is built
//!   once per (source, target, offset). Per block, one stage-1 window is filtered per
//!   group of paths that share a source, kernel and start, covering all of their
//!   segments' native coordinates, and every segment's coordinate is rebased to its
//!   window's z\[0\], so coordinates stay small for any run length or loop count.
//!
//! A transmitter trajectory with a start is wrapped in a 1 s [`Lookback`]: light time
//! looks before the start, where the start state is extrapolated (C2). The receiver is
//! evaluated only at its own acquisition instants.
//!
//! [`StaticScene`] is a convenience for fixed ECEF positions; it converts into a
//! [`Scene`] of [`Static`] trajectories.

use std::collections::BTreeMap;
use std::fmt;
use std::sync::{Arc, Mutex};

use iqsim_core::amplitude::rms_volts;
use iqsim_core::wide::U128;
use iqsim_core::{
    AmplitudeError, BlockSpec, ClockMap, DelayPoly, Headroom, NumericError, OffsetFrequency,
    ResponseError, ResponseSpec, SampleIndex, SampleRate, TapKind, TapResponse, TargetScale,
    TimeError, TruthDuration, TruthTime,
};
use iqsim_geo::trajectory::Lookback;
use iqsim_geo::{Ecef, Quat, SPEED_OF_LIGHT_M_S, Static, Trajectory};
use iqsim_render::FARROW_LEBESGUE;
use iqsim_waveform::continuous::source_response;
use iqsim_waveform::{
    BuiltEdgeTable, ContinuousError, ContinuousSource, EdgeError, EdgeShape, EdgeTableKey,
    EndPolicy, ResampleKernel, ResampleKey, SourceTimeline, build_edge_table,
};
use iqsim_worklist::{
    BlockWorkList, COORD_NATIVE_SAMPLES, COORD_TARGET_SAMPLES, ContinuousEmission, CwEmission,
    EMISSION_CONTINUOUS, EMISSION_CW, EMISSION_PULSE_TRAIN, EdgeTableArena, EmissionRef,
    FARROW_AFTER, FARROW_BEFORE, MAX_COORD_BITS, MAX_COORD_SAMPLES, MAX_WINDOW_LEN, PathSegment,
    PolyRegisters, Pulse, PulseTrain, SourceWindowArena, TargetDesc, WorkListError, target_desc,
};

use crate::geometry::{
    CoordModel, FitError, PathFit, PathGeometry, PathModel, Segment, doppler_tolerance_hz, fit_path,
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

/// A fixed emitter (see [`StaticScene`]).
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

/// A fixed receiver and fixed emitters, at exact ECEF positions. It converts into a
/// [`Scene`] of [`Static`] trajectories (identity attitude).
#[derive(Clone, Debug, PartialEq)]
pub struct StaticScene {
    /// Receiver position.
    pub receiver: Ecef,
    /// Emitters, in a fixed order.
    pub emitters: Vec<StaticEmitter>,
}

/// An emitter on a trajectory.
#[derive(Clone, Debug)]
pub struct SceneEmitter {
    /// The transmitter's motion (platform origin).
    pub trajectory: Arc<dyn Trajectory>,
    /// Carrier frequency, µHz.
    pub carrier_uhz: i64,
    /// Transmit power (pulse peak or CW), W.
    pub power_w: f64,
    /// Waveform.
    pub waveform: StaticWaveform,
}

/// A receiver and its emitters, each on a trajectory.
#[derive(Clone, Debug)]
pub struct Scene {
    /// The receiver's motion (platform origin).
    pub receiver: Arc<dyn Trajectory>,
    /// Emitters, in a fixed order.
    pub emitters: Vec<SceneEmitter>,
}

fn fixed(p: Ecef) -> Arc<dyn Trajectory> {
    Arc::new(Static::from_ecef(p, Quat::IDENTITY))
}

impl From<&StaticScene> for Scene {
    fn from(s: &StaticScene) -> Self {
        Self {
            receiver: fixed(s.receiver),
            emitters: s
                .emitters
                .iter()
                .map(|e| SceneEmitter {
                    trajectory: fixed(e.position),
                    carrier_uhz: e.carrier_uhz,
                    power_w: e.power_w,
                    waveform: e.waveform.clone(),
                })
                .collect(),
        }
    }
}

impl From<&Scene> for Scene {
    fn from(s: &Scene) -> Self {
        s.clone()
    }
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
    /// A trajectory could not be prepared.
    Geo(iqsim_geo::GeoError),
    /// A path's geometry fit failed (light time, delay range or budget).
    Fit {
        /// Emitter index.
        index: usize,
        /// Block number.
        block: u64,
        /// What failed.
        error: FitError,
    },
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
            Self::Geo(e) => e.fmt(f),
            Self::Fit {
                index,
                block,
                error,
            } => write!(f, "emitter {index}, block {block}: {error}"),
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
    ContinuousError => Continuous, iqsim_geo::GeoError => Geo);

/// A continuous link's prepared state.
struct ContinuousLink {
    spec: ContinuousSpec,
    /// Index of the shared source among the builder's distinct sources.
    source_id: usize,
    /// The source reconstruction response φ.
    phi: TapResponse,
}

/// One pulse in target coordinates, with its edge shapes (tables are built per
/// segment offset).
#[derive(Clone, Copy, Debug)]
struct PulseShape {
    start: i64,
    fall_start: i64,
    amplitude: f32,
    rise_ps: i64,
    fall_ps: i64,
}

enum LinkWave {
    Cw,
    Continuous(Box<ContinuousLink>),
    Pulses {
        emission: TruthTime,
        pulses: Vec<PulseShape>,
    },
}

struct Link {
    geometry: PathGeometry,
    carrier_uhz: i64,
    /// Amplitude × range, V·m.
    volts_m: f64,
    offset_uhz: i64,
    offset_hz: f64,
    offset: OffsetFrequency,
    wave: LinkWave,
}

/// Relative slack on the Farrow envelope bound for its f32 evaluation (a dozen
/// roundings of at most 2⁻²⁴ each, with margin).
const FARROW_ROUNDING_SLACK: f64 = 1e-5;

/// Prepared continuous segments of one block, by (link, segment index).
type PreparedWindows = BTreeMap<(usize, usize), PreparedPath>;

/// A continuous segment's stage-1 window for one block.
struct PreparedPath {
    window: u32,
    /// Native coordinate of the window's z\[0\], Q32.32 bits.
    origin: i128,
    peak_env: f64,
}

/// What one block's geometry fit did for one path.
#[derive(Clone, Debug, PartialEq)]
pub struct PathLog {
    /// Emitter index.
    pub link: usize,
    /// The fit: segments (with residuals) and splits.
    pub fit: PathFit,
}

/// The geometry-fit log of one block (contracts §14.4).
#[derive(Clone, Debug, PartialEq)]
pub struct BlockLog {
    /// Block number.
    pub block: u64,
    /// One entry per emitter, in emitter order (culled paths included).
    pub paths: Vec<PathLog>,
}

/// Builds work lists for one target of a [`Scene`].
pub struct SceneBuilder<'c> {
    response: TapResponse,
    scale: TargetScale,
    desc: TargetDesc,
    rate: SampleRate,
    cull_threshold_lsb: f64,
    clock: &'c dyn ClockMap,
    blocks: BlockSpec,
    /// Doppler tolerance of the target response (pulses and continuous sources).
    doppler_tol_hz: f64,
    /// Edge tables built so far, by key (built when a segment first needs them).
    tables: Mutex<BTreeMap<EdgeTableKey, Arc<BuiltEdgeTable>>>,
    /// Resampling kernels built so far, by key.
    kernels: Mutex<BTreeMap<ResampleKey, Arc<ResampleKernel>>>,
    links: Vec<Link>,
    /// Samples the run keeps in total, if known: a final partial block is checked
    /// against the `stop` policy only up to it.
    sample_limit: Option<u64>,
}

/// The prompt 5 name of [`SceneBuilder`].
pub type StaticSceneBuilder<'c> = SceneBuilder<'c>;

/// Look-back allowed before a transmitter trajectory's start (prompt 9 decision).
const LOOKBACK: TruthDuration = iqsim_geo::MAX_LOOKBACK;

/// The coordinate of `seg` at its index j, Q32.32 bits.
fn coord_of(seg: &Segment, j: u64) -> Result<i128, SceneError> {
    seg.coord
        .as_ref()
        .ok_or(SceneError::CoordinateRange)?
        .at(j)
        .map_err(SceneError::Numeric)
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

/// Coordinate registers: `start` (Q32.32 bits, |start| < 2⁶²) as D0, the variation's
/// D1..D3.
fn coord_registers(
    start: i128,
    variation: &DelayPoly,
    len: u64,
) -> Result<PolyRegisters, SceneError> {
    if !(-i128::from(MAX_COORD_BITS)..i128::from(MAX_COORD_BITS)).contains(&start) {
        return Err(SceneError::CoordinateRange);
    }
    let mut regs = variation.registers();
    regs[0] = U128::new(start as i64 as u64, 0);
    Ok(PolyRegisters::from_poly(&DelayPoly::from_registers(
        regs, len,
    )?))
}

impl<'c> SceneBuilder<'c> {
    /// Validates the target and prepares every link. Edge tables and resampling
    /// kernels at each emitter's nominal offset are built here (so a bad pulse or
    /// source fails at once); tables at Doppler-shifted offsets are built when a
    /// segment first needs them, before the block that uses them renders.
    pub fn new<S: Into<Scene>>(
        scene: S,
        target: &StaticTarget,
        clock: &'c dyn ClockMap,
        blocks: BlockSpec,
    ) -> Result<Self, SceneError> {
        let scene: Scene = scene.into();
        let response =
            TapResponse::new(target.kind, target.response, target.rate, target.center_uhz)?;
        let scale = TargetScale::new(target.full_scale_volts)?;
        if !(target.cull_threshold_lsb.is_finite() && target.cull_threshold_lsb >= 0.0) {
            return Err(SceneError::Amplitude(AmplitudeError::Bound));
        }
        let rate = target.rate;
        let mut links = Vec::with_capacity(scene.emitters.len());
        let mut sources: Vec<SharedSource> = Vec::new();
        for (i, e) in scene.emitters.iter().enumerate() {
            if !(e.power_w.is_finite() && e.power_w >= 0.0) || e.carrier_uhz <= 0 {
                return Err(SceneError::Emitter { index: i });
            }
            let carrier_hz = e.carrier_uhz as f64 * 1e-6;
            let lambda = SPEED_OF_LIGHT_M_S / carrier_hz;
            let volts_m = rms_volts(e.power_w) * lambda / (4.0 * core::f64::consts::PI);
            let offset_uhz = e
                .carrier_uhz
                .checked_sub(target.center_uhz)
                .ok_or(SceneError::Emitter { index: i })?;
            let offset_hz = offset_uhz as f64 * 1e-6;
            let wave = match &e.waveform {
                StaticWaveform::Cw => LinkWave::Cw,
                StaticWaveform::Continuous(spec) => {
                    if !(spec.scale.is_finite() && spec.scale >= 0.0)
                        || !(spec.peak.is_finite() && spec.peak >= 0.0)
                    {
                        return Err(SceneError::Emitter { index: i });
                    }
                    let phi = source_response(spec.rate, spec.response)?;
                    let source_id = match sources.iter().position(|s| *s == spec.source) {
                        Some(p) => p,
                        None => {
                            sources.push(spec.source.clone());
                            sources.len() - 1
                        }
                    };
                    LinkWave::Continuous(Box::new(ContinuousLink {
                        spec: spec.clone(),
                        source_id,
                        phi,
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
                        built.push(PulseShape {
                            start: ps_to_coord(p.start_ps, rate)?,
                            fall_start: ps_to_coord(p.fall_start_ps, rate)?,
                            amplitude: p.amplitude,
                            rise_ps: p.rise_ps,
                            fall_ps: p.fall_ps,
                        });
                    }
                    LinkWave::Pulses {
                        emission: *emission,
                        pulses: built,
                    }
                }
            };
            let tx: Arc<dyn Trajectory> =
                Arc::new(Lookback::new(Arc::clone(&e.trajectory), LOOKBACK)?);
            links.push(Link {
                geometry: PathGeometry {
                    rx: Arc::clone(&scene.receiver),
                    tx,
                },
                carrier_uhz: e.carrier_uhz,
                volts_m,
                offset_uhz,
                offset_hz,
                offset: OffsetFrequency::from_uhz(offset_uhz)?,
                wave,
            });
        }
        let builder = Self {
            desc: target_desc(target.kind, response.id(), scale.lsb_per_volt()),
            doppler_tol_hz: doppler_tolerance_hz(&response),
            response,
            scale,
            rate,
            cull_threshold_lsb: target.cull_threshold_lsb,
            clock,
            blocks,
            tables: Mutex::new(BTreeMap::new()),
            kernels: Mutex::new(BTreeMap::new()),
            links,
            sample_limit: None,
        };
        for link in &builder.links {
            match &link.wave {
                LinkWave::Cw => {}
                LinkWave::Continuous(c) => {
                    builder.kernel(&c.phi, link.offset_uhz)?;
                }
                LinkWave::Pulses { pulses, .. } => {
                    for p in pulses {
                        builder.table(p.rise_ps, link.offset_uhz)?;
                        builder.table(p.fall_ps, link.offset_uhz)?;
                    }
                }
            }
        }
        Ok(builder)
    }

    /// The edge table for a linear edge of `rise_ps` at `offset_uhz`, built once.
    fn table(
        &self,
        rise_ps: i64,
        offset_uhz: i64,
    ) -> Result<(EdgeTableKey, Arc<BuiltEdgeTable>), SceneError> {
        let key = EdgeTableKey {
            response: self.response.id(),
            shape: EdgeShape::Linear { rise_ps },
            offset_uhz,
        };
        let mut cache = self
            .tables
            .lock()
            .map_err(|_| ContinuousError::Source("table cache lock poisoned".to_owned()))?;
        if let Some(t) = cache.get(&key) {
            return Ok((key, Arc::clone(t)));
        }
        let t = Arc::new(build_edge_table(&self.response, key.shape, offset_uhz)?);
        cache.insert(key, Arc::clone(&t));
        Ok((key, t))
    }

    /// The resampling kernel of φ into this target at `offset_uhz`, built once.
    fn kernel(
        &self,
        phi: &TapResponse,
        offset_uhz: i64,
    ) -> Result<Arc<ResampleKernel>, SceneError> {
        let key = ResampleKey {
            source: phi.id(),
            target: self.response.id(),
            offset_uhz,
        };
        let mut cache = self
            .kernels
            .lock()
            .map_err(|_| ContinuousError::Source("kernel cache lock poisoned".to_owned()))?;
        if let Some(k) = cache.get(&key) {
            return Ok(Arc::clone(k));
        }
        let k = Arc::new(ResampleKernel::build(phi, &self.response, offset_uhz)?);
        cache.insert(key, Arc::clone(&k));
        Ok(k)
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

    /// The geometry of emitter `index`'s path (receiver and look-back transmitter).
    pub fn path_geometry(&self, index: usize) -> Option<&PathGeometry> {
        self.links.get(index).map(|l| &l.geometry)
    }

    fn model<'a>(&'a self, link: &'a Link) -> PathModel<'a> {
        let coord = match &link.wave {
            LinkWave::Cw => None,
            LinkWave::Continuous(c) => Some(CoordModel {
                rate: c.spec.rate,
                origin: c.spec.start,
            }),
            LinkWave::Pulses { emission, .. } => Some(CoordModel {
                rate: self.rate,
                origin: *emission,
            }),
        };
        let cw = matches!(link.wave, LinkWave::Cw);
        PathModel {
            geometry: &link.geometry,
            carrier_uhz: link.carrier_uhz,
            offset: link.offset,
            offset_hz: link.offset_hz,
            offset_uhz: link.offset_uhz,
            target_rate: self.rate,
            coord,
            fold: cw.then_some(&self.response),
            filtered: (!cw).then_some((&self.response, self.doppler_tol_hz)),
            volts_m: link.volts_m,
        }
    }

    /// Builds the work list of block `b`: fits every path, culls segments below the
    /// threshold, proves headroom and validates the result.
    pub fn build_block(&self, b: u64) -> Result<BlockWorkList, SceneError> {
        self.build_block_logged(b).map(|(w, _)| w)
    }

    /// [`build_block`](Self::build_block), also returning the block's fit log.
    pub fn build_block_logged(&self, b: u64) -> Result<(BlockWorkList, BlockLog), SceneError> {
        let block = self.blocks.block(b)?;
        let n0 = block.start();
        let len = block.len();
        let mut fits = Vec::with_capacity(self.links.len());
        for (i, link) in self.links.iter().enumerate() {
            let fit = fit_path(&self.model(link), self.clock, n0, len).map_err(|error| {
                SceneError::Fit {
                    index: i,
                    block: b,
                    error,
                }
            })?;
            fits.push(fit);
        }

        let (windows, prepared) = self.continuous_windows(b, n0, len, &fits)?;
        let mut arena = EdgeTableArena::default();
        let mut arena_index: BTreeMap<EdgeTableKey, u32> = BTreeMap::new();
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
            tables: Arc::new(EdgeTableArena::default()),
            windows: Arc::new(windows),
        };
        let mut headroom = Headroom::new();
        let mut cw_emission: Option<EmissionRef> = None;
        for (source, link) in self.links.iter().enumerate() {
            let fit = &fits[source];
            let mut path_bound = 0.0_f64;
            let mut any = false;
            if matches!(link.wave, LinkWave::Cw) {
                cw_emission = None;
            }
            for (si, seg) in fit.segments.iter().enumerate() {
                // Per segment: the envelope bound, the emission, the window origin.
                let (peak_env, emission, origin, native) = match &link.wave {
                    LinkWave::Cw => (1.0, None, 0_i128, false),
                    LinkWave::Continuous(_) => {
                        let Some(p) = prepared.get(&(source, si)) else {
                            continue;
                        };
                        (p.peak_env, Some(Err(p.window)), p.origin, true)
                    }
                    LinkWave::Pulses { pulses, .. } => {
                        // The coordinate is non-decreasing, so the segment's ends bound
                        // it.
                        let lo = coord_of(seg, 0)?;
                        let hi = coord_of(seg, seg.len - 1)?;
                        let limit = i128::from(MAX_COORD_BITS);
                        // Clamped to ±2³⁰ samples only for the activity test.
                        let (lo_c, hi_c) = (
                            lo.clamp(-limit, limit) as i64,
                            hi.clamp(-limit, limit) as i64,
                        );
                        // Activity from the nominal-offset tables (built at
                        // construction). A table's extent is its response support,
                        // the same at any offset up to one cell of resolution, so one
                        // sample of margin keeps the test conservative; an extra pulse
                        // only adds zeros.
                        let mut active = Vec::new();
                        for p in pulses {
                            let (_, rt) = self.table(p.rise_ps, link.offset_uhz)?;
                            let (_, ft) = self.table(p.fall_ps, link.offset_uhz)?;
                            let end = |t: &BuiltEdgeTable| {
                                t.u_lo + ((t.points.len() as i64 - 1) << (32 - t.log2_res))
                            };
                            let margin = 1_i64 << 32;
                            let on = (p.start + rt.u_lo).min(p.fall_start + ft.u_lo) - margin;
                            let off = (p.start + end(&rt)).max(p.fall_start + end(&ft)) + margin;
                            if on <= hi_c && off >= lo_c {
                                let (rk, rt) = self.table(p.rise_ps, seg.offset_uhz)?;
                                let (fk, ft) = self.table(p.fall_ps, seg.offset_uhz)?;
                                active.push((p, rk, rt, fk, ft));
                            }
                        }
                        if active.is_empty() {
                            continue;
                        }
                        if !(lo.abs() < limit && hi.abs() < limit) {
                            return Err(SceneError::CoordinateRange);
                        }
                        let mut built = Vec::with_capacity(active.len());
                        for (p, rk, rt, fk, ft) in active {
                            let mut id =
                                |k: EdgeTableKey, t: &BuiltEdgeTable| -> Result<u32, SceneError> {
                                    if let Some(&i) = arena_index.get(&k) {
                                        return Ok(i);
                                    }
                                    let i = arena.push(t.u_lo, t.log2_res, &t.points)?;
                                    arena_index.insert(k, i);
                                    Ok(i)
                                };
                            let rise_table = id(rk, &rt)?;
                            let fall_table = id(fk, &ft)?;
                            built.push(Pulse {
                                start: p.start,
                                fall_start: p.fall_start,
                                amplitude: p.amplitude,
                                rise_table,
                                fall_table,
                                _pad: 0,
                            });
                        }
                        // Peak envelope: the largest sum of A·(peak_rise + peak_fall)
                        // over pulses whose extents overlap one point (a maximum is
                        // reached at some pulse's start).
                        let extent = |p: &Pulse| {
                            let (on, off) = arena.pulse_extent(p);
                            let pk = f64::from(arena.headers[p.rise_table as usize].peak)
                                + f64::from(arena.headers[p.fall_table as usize].peak);
                            (on, off, f64::from(p.amplitude.abs()) * pk)
                        };
                        let mut peak = 0.0_f64;
                        for a in &built {
                            let (on_a, _, _) = extent(a);
                            let s: f64 = built
                                .iter()
                                .map(extent)
                                .filter(|&(on, off, _)| on <= on_a && off >= on_a)
                                .map(|(_, _, v)| v)
                                .sum();
                            peak = peak.max(s);
                        }
                        (peak, Some(Ok(built)), 0, false)
                    }
                };
                // Cull segment by segment (lossless at the default threshold: every
                // component of a culled segment rounds to 0), and bound the path by its
                // largest segment (segments never overlap in time).
                let bound = seg.amplitude.bound();
                if self.scale.peak_lsb(bound, peak_env)? < self.cull_threshold_lsb {
                    continue;
                }
                // The envelope bound is taken as at least 1, so the proof also covers
                // the kernel's intermediate m = a·lsb_per_volt (contracts §2).
                path_bound = path_bound.max(self.scale.peak_lsb(bound, peak_env.max(1.0))?);
                any = true;
                let emission = match emission {
                    Some(Err(window)) => {
                        work.continuous.push(ContinuousEmission { window, _pad: 0 });
                        EmissionRef {
                            kind: EMISSION_CONTINUOUS,
                            index: (work.continuous.len() - 1) as u32,
                        }
                    }
                    Some(Ok(pulses)) => {
                        let first = work.pulses.len() as u32;
                        let count = pulses.len() as u32;
                        work.pulses.extend(pulses);
                        work.trains.push(PulseTrain {
                            first_pulse: first,
                            pulse_count: count,
                        });
                        EmissionRef {
                            kind: EMISSION_PULSE_TRAIN,
                            index: (work.trains.len() - 1) as u32,
                        }
                    }
                    None => *cw_emission.get_or_insert_with(|| {
                        work.cw.push(CwEmission {
                            envelope: [1.0, 0.0],
                        });
                        EmissionRef {
                            kind: EMISSION_CW,
                            index: (work.cw.len() - 1) as u32,
                        }
                    }),
                };
                let coord = match &seg.coord {
                    Some(c) => coord_registers(c.start - origin, &c.variation, seg.len)?,
                    None => PolyRegisters::from_poly(&DelayPoly::from_f64_differences(
                        [0.0; 4], seg.len,
                    )?),
                };
                work.paths.push(PathSegment {
                    delay: PolyRegisters::from_poly(&seg.delay),
                    coord,
                    phase: PolyRegisters::from_poly(&seg.phase),
                    block: b,
                    amplitude: seg.amplitude,
                    target: 0,
                    coord_units: if native {
                        COORD_NATIVE_SAMPLES
                    } else {
                        COORD_TARGET_SAMPLES
                    },
                    coord_source: source as u32,
                    emission,
                    first: seg.first as u32,
                    len: seg.len as u32,
                });
            }
            if any {
                headroom.add(path_bound)?;
            }
        }
        work.tables = Arc::new(arena);
        work.validate()?;
        let log = BlockLog {
            block: b,
            paths: fits
                .into_iter()
                .enumerate()
                .map(|(link, fit)| PathLog { link, fit })
                .collect(),
        };
        Ok((work, log))
    }

    /// Stage 1 for block `b`: the windows of every continuous segment that is not
    /// culled, keyed by (link, segment index), with origin and envelope bound.
    ///
    /// Segments that share a source, kernel (so a Doppler-shifted offset), start,
    /// timeline and level factor share one window covering all of their coordinates.
    fn continuous_windows(
        &self,
        b: u64,
        n0: SampleIndex,
        len: u64,
        fits: &[PathFit],
    ) -> Result<(SourceWindowArena, PreparedWindows), SceneError> {
        struct Pending {
            key: (usize, usize),
            /// Absolute native coordinate at the segment start, Q32.32 bits.
            w0: i128,
            /// Absolute native coordinate at the segment end, Q32.32 bits.
            w1: i128,
            kernel: Arc<ResampleKernel>,
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
            // `stop`: the exact retarded coordinate of the last kept sample (review
            // fix: the fitted one can fall short by up to the fit residual).
            if keep > 0 {
                let w = self
                    .model(link)
                    .exact_coord(self.clock, n0.checked_add(keep - 1)?)
                    .map_err(|error| SceneError::Fit {
                        index: i,
                        block: b,
                        error,
                    })?;
                if c.spec.timeline.exhausted(w) {
                    return Err(SceneError::SourceExhausted { index: i, block: b });
                }
            }
            for (si, seg) in fits[i].segments.iter().enumerate() {
                let kernel = self.kernel(&c.phi, seg.offset_uhz)?;
                // Lossless pre-filter cull: |z| ≤ l1_max·peak|x|, so nothing in this
                // segment can reach the threshold.
                let bound =
                    FARROW_LEBESGUE * kernel.l1_max() * c.spec.peak * (1.0 + FARROW_ROUNDING_SLACK);
                if self
                    .scale
                    .peak_lsb(seg.amplitude.bound() * c.spec.scale, bound)?
                    < self.cull_threshold_lsb
                {
                    continue;
                }
                let t = c.spec.timeline;
                let key = (
                    c.source_id,
                    kernel.key(),
                    c.spec.start.as_ps(),
                    t.len,
                    t.policy,
                    c.spec.scale.to_bits(),
                );
                groups.entry(key).or_default().push(Pending {
                    key: (i, si),
                    w0: coord_of(seg, 0)?,
                    w1: coord_of(seg, seg.len - 1)?,
                    kernel,
                });
            }
        }
        let mut arena = SourceWindowArena::default();
        let mut out = BTreeMap::new();
        for pending in groups.values() {
            let LinkWave::Continuous(c) = &self.links[pending[0].key.0].wave else {
                unreachable!("grouped links are continuous");
            };
            let kernel = &pending[0].kernel;
            let log2_up = kernel.log2_up();
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
                    index: pending[0].key.0,
                    points: z_len as u64,
                });
            }
            let (first, last) = kernel.source_span(m_lo, z_len);
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
                kernel.filter_window(&c.spec.timeline, src.as_mut(), m_lo, z_len, c.spec.scale)?
            };
            let window = arena.push(log2_up, &z)?;
            let peak = f64::from(arena.headers[window as usize].peak);
            let origin = i128::from(m_lo) << shift;
            for p in pending {
                let local = p.w0 - origin;
                if !(0..i128::from(MAX_COORD_BITS)).contains(&local) {
                    return Err(SceneError::CoordinateRange);
                }
                out.insert(
                    p.key,
                    PreparedPath {
                        window,
                        origin,
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
