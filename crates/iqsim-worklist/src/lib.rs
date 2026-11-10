//! Flat block work lists shared by all render backends.
//! Per the plan section "Architecture", this crate owns pointer-free path,
//! emission, and immutable-table references for each block and output tap.
//! Delay, phase, and retarded waveform coordinates preserve exact render timing.
//!
//! A [`BlockWorkList`] is everything one render target needs for one block. It is
//! the only render input (contracts §5.2). Its element types are `#[repr(C)]` plain data
//! with explicit padding and no pointers: every reference is a `u32` index into a flat
//! array of the same work list, or into its [`EdgeTableArena`]. A GPU backend uploads
//! the slices unchanged.
//!
//! # Coordinates and units
//!
//! - **Local sample index** k ∈ [0, len) within the block. A [`PathSegment`] covers
//!   the samples `[first, first + len)` of its block (prompt 9: the geometry engine
//!   splits a path at trajectory knots and wherever its fit needs refining), and every
//!   polynomial of the segment is in the **segment-local index** j = k − first over
//!   its own `len` samples. A path whose fit needs no split is one segment covering the
//!   whole block.
//! - **Delay** ([`PathSegment::delay`]): Q32.32 output samples (`DelayPoly` registers).
//!   It carries the path delay for truth and diagnostics. The renderer reads time from
//!   the waveform coordinate instead.
//! - **Waveform coordinate** ([`PathSegment::coord`]): the retarded time
//!   `w(k) = F_s·(t_acq(k) − τ(k) − t_ref)` in Q32.32 *target* samples
//!   ([`COORD_TARGET_SAMPLES`]). `t_acq` is the ClockMap acquisition instant, so clock
//!   error moves pulse timing, and `t_ref` is the origin of the source segment
//!   ([`PathSegment::coord_source`]): for a pulse train, its emission reference time.
//!   **Rebasing:** a builder must keep |w| < 2³⁰ samples over the block
//!   ([`MAX_COORD_SAMPLES`]), choosing a later origin for long-running sources, so no
//!   fixed-point coordinate can overflow in a long run. Pulse offsets obey the same
//!   bound, so `w − offset` always fits in i64.
//!   A continuous emission ([`EMISSION_CONTINUOUS`]) uses [`COORD_NATIVE_SAMPLES`]
//!   instead: Q32.32 *native source samples* relative to z\[0\] of its stage-1 window
//!   ([`SourceWindowHeader`]), so the window origin is the rebasing origin.
//! - **Phase** ([`PathSegment::phase`]): the total carrier rotation in cycles,
//!   `f_off·t_acq(k) − f_c·τ(k)` plus any folded response phase (`PhasePoly`
//!   registers, wrapping by design).
//! - **Amplitude** ([`AmplitudePoly`]): volts at the target's reference plane, an f32
//!   cubic in x = j·(1/len) (segment index and length), evaluated by Horner from scratch at every sample so the
//!   value never depends on where a chunk starts.
//!
//! Dependencies: iqsim-core.

#![forbid(unsafe_code)]

use std::fmt;
use std::sync::Arc;

use iqsim_core::wide::U128;
use iqsim_core::{BlockPoly, NumericError, PolyOutput, TapKind, TapResponseId};
use sha2::{Digest, Sha256};

/// Coordinate units code: Q32.32 target samples relative to the source origin.
pub const COORD_TARGET_SAMPLES: u32 = 1;
/// Coordinate units code: Q32.32 native-rate source samples relative to z\[0\] of the
/// path's source window (continuous emissions, contracts §5.6).
pub const COORD_NATIVE_SAMPLES: u32 = 2;

/// Bound on |waveform coordinate| and |pulse offset|, in samples (2³⁰).
pub const MAX_COORD_SAMPLES: i64 = 1 << 30;
/// [`MAX_COORD_SAMPLES`] in Q32.32 bits (2⁶²).
pub const MAX_COORD_BITS: i64 = MAX_COORD_SAMPLES << 32;

/// Emission kind code: [`CwEmission`].
pub const EMISSION_CW: u32 = 1;
/// Emission kind code: [`PulseTrain`].
pub const EMISSION_PULSE_TRAIN: u32 = 2;
/// Emission kind code: [`ContinuousEmission`].
pub const EMISSION_CONTINUOUS: u32 = 3;

/// Farrow taps before the interpolation cell (nodes −2, −1).
pub const FARROW_BEFORE: u64 = 2;
/// Farrow taps from the cell start on (nodes 0..3).
pub const FARROW_AFTER: u64 = 4;
/// Largest stage-1 oversampling exponent: U = 2^log2_up ≤ 32 grid points per native
/// sample.
pub const MAX_WINDOW_LOG2_UP: u32 = 5;
/// Largest source window, z samples (2²⁶), so `w·U` can never leave Q32.32.
pub const MAX_WINDOW_LEN: u32 = 1 << 26;

/// Largest accepted |pulse amplitude|. With the headroom proof's amplitude-only term
/// (contracts §2) it keeps every f32 intermediate of the kernel finite.
pub const MAX_PULSE_AMPLITUDE: f32 = 1e6;

/// Largest edge-table resolution exponent: 2¹² points per target sample.
pub const MAX_TABLE_LOG2_RES: u32 = 12;

/// Raw Newton-form registers D0..D3 of a block polynomial as `[lo, hi]` word pairs.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct PolyRegisters {
    /// `[D0.lo, D0.hi, D1.lo, D1.hi, D2.lo, D2.hi, D3.lo, D3.hi]`.
    pub words: [u64; 8],
}

impl PolyRegisters {
    /// The registers of `p`.
    pub fn from_poly<T: PolyOutput>(p: &BlockPoly<T>) -> Self {
        let d = p.registers();
        let mut words = [0; 8];
        for (j, r) in d.iter().enumerate() {
            words[2 * j] = r.lo;
            words[2 * j + 1] = r.hi;
        }
        Self { words }
    }

    /// The polynomial over `len` samples with these registers.
    pub fn to_poly<T: PolyOutput>(&self, len: u64) -> Result<BlockPoly<T>, NumericError> {
        let w = &self.words;
        let d = [0, 1, 2, 3].map(|j| U128::new(w[2 * j + 1], w[2 * j]));
        BlockPoly::from_registers(d, len)
    }
}

/// Path amplitude: `c0 + x·(c1 + x·(c2 + x·c3))` volts with `x = f32(k)·inv_len`.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct AmplitudePoly {
    /// Power-basis coefficients in x, volts.
    pub c: [f32; 4],
    /// f32(1/len).
    pub inv_len: f32,
}

impl AmplitudePoly {
    /// A constant amplitude over a block of `len` samples.
    pub fn constant(volts: f32, len: u64) -> Self {
        Self {
            c: [volts, 0.0, 0.0, 0.0],
            inv_len: (1.0 / len as f64) as f32,
        }
    }

    /// Value at local index k (exact f32 for k < 2²⁴, which covers every block).
    /// Every operation is separate, in this order, on every backend.
    #[inline]
    pub fn eval(&self, k: u64) -> f32 {
        let x = (k as f32) * self.inv_len;
        let t = self.c[3] * x;
        let t = self.c[2] + t;
        let t = t * x;
        let t = self.c[1] + t;
        let t = t * x;
        self.c[0] + t
    }

    /// Bound on |a| for x ∈ [0, 1]: Σ|c_j|, widened by 2⁻²⁰ relative for the f32
    /// evaluation.
    pub fn bound(&self) -> f64 {
        let s: f64 = self.c.iter().map(|c| f64::from(c.abs())).sum();
        s * (1.0 + 1.0 / 1_048_576.0)
    }
}

/// A reference from a path to its emission.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct EmissionRef {
    /// [`EMISSION_CW`], [`EMISSION_PULSE_TRAIN`] or [`EMISSION_CONTINUOUS`].
    pub kind: u32,
    /// Index into [`BlockWorkList::cw`], [`BlockWorkList::trains`] or
    /// [`BlockWorkList::continuous`].
    pub index: u32,
}

/// One propagation path of one emission into one render target, over the samples
/// `[first, first + len)` of one block. Its polynomials are in the segment-local index
/// j = k − first, j ∈ [0, len).
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PathSegment {
    /// Delay polynomial registers (Q32.96 samples).
    pub delay: PolyRegisters,
    /// Retarded waveform-coordinate registers (Q32.96 in [`Self::coord_units`]).
    pub coord: PolyRegisters,
    /// Carrier-phase registers (0.128 cycles).
    pub phase: PolyRegisters,
    /// Block number this segment belongs to.
    pub block: u64,
    /// Amplitude, volts.
    pub amplitude: AmplitudePoly,
    /// Render target index (always 0 within one work list: one list per target).
    pub target: u32,
    /// [`COORD_TARGET_SAMPLES`], or [`COORD_NATIVE_SAMPLES`] for a continuous emission.
    pub coord_units: u32,
    /// Source segment whose origin `t_ref` the coordinate is relative to.
    pub coord_source: u32,
    /// The emission.
    pub emission: EmissionRef,
    /// Block-local index of the segment's first sample.
    pub first: u32,
    /// Samples in the segment (≥ 1; `first + len` ≤ the block length).
    pub len: u32,
}

/// A continuous-wave emission: a constant complex envelope (normally 1 + 0j). Any tap
/// response is already folded into the path amplitude and phase.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CwEmission {
    /// Complex envelope (re, im), unitless.
    pub envelope: [f32; 2],
}

/// A continuous emission: the path's envelope is the stage-1 window `window`
/// interpolated by the 6-tap Farrow kernel at the path's native coordinate.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct ContinuousEmission {
    /// Index into [`SourceWindowArena::headers`].
    pub window: u32,
    /// Explicit padding (zero).
    pub _pad: u32,
}

/// Header of one stage-1 window: the source filtered to one target's response, on a
/// grid of U = 2^log2_up points per native sample. Point i is z at native coordinate
/// i/U relative to the window origin.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SourceWindowHeader {
    /// First point in [`SourceWindowArena::samples`].
    pub first: u32,
    /// Number of points (≥ [`FARROW_BEFORE`] + [`FARROW_AFTER`], ≤ [`MAX_WINDOW_LEN`]).
    pub len: u32,
    /// log2 of the oversampling U, ≤ [`MAX_WINDOW_LOG2_UP`].
    pub log2_up: u32,
    /// Explicit padding (zero).
    pub _pad: u32,
    /// max |point|, the envelope bound before the Farrow kernel's Lebesgue constant.
    pub peak: f32,
    /// Explicit padding (zero).
    pub _pad2: u32,
}

/// Immutable stage-1 windows of one block, shared by every path that references them.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SourceWindowArena {
    /// Window headers.
    pub headers: Vec<SourceWindowHeader>,
    /// Complex points (re, im) of every window, concatenated.
    pub samples: Vec<[f32; 2]>,
}

impl SourceWindowArena {
    /// Appends a window and returns its index.
    pub fn push(&mut self, log2_up: u32, samples: &[[f32; 2]]) -> Result<u32, WorkListError> {
        let len = u32::try_from(samples.len()).map_err(|_| WorkListError::Window)?;
        if log2_up > MAX_WINDOW_LOG2_UP
            || u64::from(len) < FARROW_BEFORE + FARROW_AFTER
            || len > MAX_WINDOW_LEN
            || !samples.iter().all(|p| p[0].is_finite() && p[1].is_finite())
        {
            return Err(WorkListError::Window);
        }
        // In f64, so neither tiny nor huge components underflow or overflow, then
        // rounded up to the next f32 so the stored peak is a bound (review fix).
        let peak64 = samples
            .iter()
            .map(|p| {
                let (r, i) = (f64::from(p[0]), f64::from(p[1]));
                (r * r + i * i).sqrt()
            })
            .fold(0.0_f64, f64::max);
        let mut peak = peak64 as f32;
        if f64::from(peak) < peak64 {
            peak = f32::from_bits(peak.to_bits() + 1);
        }
        if !peak.is_finite() {
            return Err(WorkListError::Window);
        }
        let first = u32::try_from(self.samples.len()).map_err(|_| WorkListError::Window)?;
        let index = u32::try_from(self.headers.len()).map_err(|_| WorkListError::Window)?;
        self.samples.extend_from_slice(samples);
        self.headers.push(SourceWindowHeader {
            first,
            len,
            log2_up,
            _pad: 0,
            peak,
            _pad2: 0,
        });
        Ok(index)
    }
}

/// A run of pulses sharing one emission reference time.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct PulseTrain {
    /// First index in [`BlockWorkList::pulses`].
    pub first_pulse: u32,
    /// Number of pulses, sorted by `start`.
    pub pulse_count: u32,
}

/// One pulse: `amplitude·[G_rise(w − start) − G_fall(w − fall_start)]`.
///
/// Both edge tables rise from 0 to the same plateau, so the pulse is zero before
/// `start + rise.u_lo` and exactly zero again after `fall_start + fall.u_hi`.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Pulse {
    /// Rise-edge start, Q32.32 target samples relative to the emission reference time.
    pub start: i64,
    /// Fall-edge start, Q32.32 target samples relative to the emission reference time.
    pub fall_start: i64,
    /// Relative pulse amplitude (unitless).
    pub amplitude: f32,
    /// Rise edge table, index into the arena.
    pub rise_table: u32,
    /// Fall edge table, index into the arena.
    pub fall_table: u32,
    /// Explicit padding (zero).
    pub _pad: u32,
}

/// Header of one band-limited edge table in an [`EdgeTableArena`].
///
/// Point i is at coordinate `u_lo + i·2^(−log2_res)` target samples. Below `u_lo` the
/// edge is exactly 0. From the last point on it is exactly `plateau` (the last stored
/// point equals `plateau`). In between it is linearly interpolated.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct EdgeTableHeader {
    /// Coordinate of point 0, Q32.32 target samples.
    pub u_lo: i64,
    /// First point in [`EdgeTableArena::points`].
    pub first: u32,
    /// Number of points (≥ 2).
    pub len: u32,
    /// Points per target sample = 2^log2_res, ≤ 2^[`MAX_TABLE_LOG2_RES`].
    pub log2_res: u32,
    /// Explicit padding (zero).
    pub _pad: u32,
    /// Final value (re, im).
    pub plateau: [f32; 2],
    /// max |point|, which also bounds every interpolated value.
    pub peak: f32,
    /// Explicit padding (zero).
    pub _pad2: u32,
}

/// Immutable edge tables shared by every block that references them.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct EdgeTableArena {
    /// Table headers.
    pub headers: Vec<EdgeTableHeader>,
    /// Complex points (re, im) of every table, concatenated.
    pub points: Vec<[f32; 2]>,
}

impl EdgeTableArena {
    /// Appends a table and returns its index. `points` must hold at least two points,
    /// and its last point must equal `plateau`.
    pub fn push(
        &mut self,
        u_lo: i64,
        log2_res: u32,
        points: &[[f32; 2]],
    ) -> Result<u32, WorkListError> {
        let (Some(&last), true) = (points.last(), points.len() >= 2) else {
            return Err(WorkListError::Table);
        };
        if log2_res > MAX_TABLE_LOG2_RES || u_lo.abs() >= MAX_COORD_BITS {
            return Err(WorkListError::Table);
        }
        let peak = points
            .iter()
            .map(|p| (p[0] * p[0] + p[1] * p[1]).sqrt())
            .fold(0.0_f32, f32::max);
        let first = u32::try_from(self.points.len()).map_err(|_| WorkListError::Table)?;
        let len = u32::try_from(points.len()).map_err(|_| WorkListError::Table)?;
        let index = u32::try_from(self.headers.len()).map_err(|_| WorkListError::Table)?;
        self.points.extend_from_slice(points);
        self.headers.push(EdgeTableHeader {
            u_lo,
            first,
            len,
            log2_res,
            _pad: 0,
            plateau: last,
            peak,
            _pad2: 0,
        });
        Ok(index)
    }

    /// Coordinate (Q32.32 target samples) of a table's last point.
    pub fn u_hi(&self, table: u32) -> i64 {
        let h = &self.headers[table as usize];
        h.u_lo + ((i64::from(h.len) - 1) << (32 - h.log2_res))
    }

    /// Coordinates `[on, off]` outside which a pulse is exactly zero: from the earlier
    /// table start to the later table end of its two shifted edges. A long rise can
    /// end after the fall's table, so `off` is the maximum of both ends.
    pub fn pulse_extent(&self, p: &Pulse) -> (i64, i64) {
        let r = &self.headers[p.rise_table as usize];
        let f = &self.headers[p.fall_table as usize];
        let on = (p.start + r.u_lo).min(p.fall_start + f.u_lo);
        let off = (p.start + self.u_hi(p.rise_table)).max(p.fall_start + self.u_hi(p.fall_table));
        (on, off)
    }
}

/// Description of the render target a work list is for.
#[repr(C)]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TargetDesc {
    /// Response identity (tap kind, response, rate, centre).
    pub response: [u8; 32],
    /// [`TapKind`] code.
    pub tap_kind: u32,
    /// The target's f32 volts → LSB factor.
    pub lsb_per_volt: f32,
}

/// Everything one render target needs for one block.
#[derive(Clone, Debug, PartialEq)]
pub struct BlockWorkList {
    /// The target.
    pub target: TargetDesc,
    /// Block number.
    pub block: u64,
    /// Channel sample index of local index 0.
    pub first_sample: u64,
    /// Block length (samples).
    pub len: u64,
    /// Paths, in builder order (order cannot change the integer sum).
    pub paths: Vec<PathSegment>,
    /// CW emissions.
    pub cw: Vec<CwEmission>,
    /// Pulse trains.
    pub trains: Vec<PulseTrain>,
    /// Pulses of every train.
    pub pulses: Vec<Pulse>,
    /// Continuous emissions.
    pub continuous: Vec<ContinuousEmission>,
    /// Edge tables (shared, immutable).
    pub tables: Arc<EdgeTableArena>,
    /// Stage-1 source windows of this block (shared, immutable).
    pub windows: Arc<SourceWindowArena>,
}

/// A malformed work list or table.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum WorkListError {
    /// A path's block, target, units, polynomial length or emission is invalid.
    Path {
        /// Index of the offending path.
        index: usize,
    },
    /// A pulse train or pulse is out of range.
    Pulse {
        /// Index of the offending pulse train.
        train: usize,
    },
    /// An edge table header is malformed.
    Table,
    /// A source window header is malformed.
    Window,
}

impl fmt::Display for WorkListError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Path { index } => write!(f, "work list path {index} is malformed"),
            Self::Pulse { train } => write!(f, "work list pulse train {train} is malformed"),
            Self::Table => f.write_str("edge table is malformed"),
            Self::Window => f.write_str("source window is malformed"),
        }
    }
}

impl std::error::Error for WorkListError {}

impl BlockWorkList {
    /// Checks every index, code, length and coordinate bound a backend relies on.
    /// Backends call this before rendering.
    pub fn validate(&self) -> Result<(), WorkListError> {
        let tables = &self.tables;
        for h in &tables.headers {
            let end = u64::from(h.first) + u64::from(h.len);
            if h.len < 2
                || end > tables.points.len() as u64
                || h.log2_res > MAX_TABLE_LOG2_RES
                || h.u_lo.abs() >= MAX_COORD_BITS
                || tables.points[end as usize - 1] != h.plateau
            {
                return Err(WorkListError::Table);
            }
        }
        let windows = &self.windows;
        for h in &windows.headers {
            let end = u64::from(h.first) + u64::from(h.len);
            if u64::from(h.len) < FARROW_BEFORE + FARROW_AFTER
                || h.len > MAX_WINDOW_LEN
                || h.log2_up > MAX_WINDOW_LOG2_UP
                || end > windows.samples.len() as u64
            {
                return Err(WorkListError::Window);
            }
        }
        if self.len == 0 || self.len > iqsim_core::MAX_BLOCK_LEN {
            return Err(WorkListError::Path { index: 0 });
        }
        for (index, p) in self.paths.iter().enumerate() {
            let bad = || WorkListError::Path { index };
            let i = p.emission.index as usize;
            // The window a continuous path reads, if any.
            let mut window = None;
            let emission_ok = match p.emission.kind {
                EMISSION_CW => i < self.cw.len() && p.coord_units == COORD_TARGET_SAMPLES,
                EMISSION_PULSE_TRAIN => {
                    i < self.trains.len() && p.coord_units == COORD_TARGET_SAMPLES
                }
                EMISSION_CONTINUOUS => {
                    window = self
                        .continuous
                        .get(i)
                        .and_then(|c| windows.headers.get(c.window as usize));
                    window.is_some() && p.coord_units == COORD_NATIVE_SAMPLES
                }
                _ => false,
            };
            if p.block != self.block
                || p.target != 0
                || !emission_ok
                || !p.amplitude.c.iter().all(|c| c.is_finite())
                || !p.amplitude.inv_len.is_finite()
            {
                return Err(bad());
            }
            let seg_len = u64::from(p.len);
            if p.len == 0 || u64::from(p.first) + seg_len > self.len {
                return Err(bad());
            }
            let coord = p
                .coord
                .to_poly::<iqsim_core::Delay>(seg_len)
                .map_err(|_| bad())?;
            p.delay
                .to_poly::<iqsim_core::Delay>(seg_len)
                .map_err(|_| bad())?;
            p.phase
                .to_poly::<iqsim_core::Phase>(seg_len)
                .map_err(|_| bad())?;
            // The coordinate must be non-decreasing over the segment, and its exact
            // value must not wrap mod 2¹²⁸ anywhere in it (contracts §5.2). Then the
            // endpoints bound every interior value (review fix, prompt 6a: a huge slope
            // could wrap back into range).
            if !coordinate_is_monotone(&coord.registers(), seg_len) {
                return Err(bad());
            }
            for k in [0, seg_len - 1] {
                let w = coord.eval(k).map_err(|_| bad())?.to_bits();
                if w.abs() >= MAX_COORD_BITS {
                    return Err(bad());
                }
                // Every Farrow tap of a continuous path lies inside its window. The
                // coordinate is non-decreasing, so the endpoints bound it.
                if let Some(h) = window {
                    let cell = w >> (32 - h.log2_up);
                    if cell < FARROW_BEFORE as i64 || cell + FARROW_AFTER as i64 > i64::from(h.len)
                    {
                        return Err(bad());
                    }
                }
            }
        }
        for (train, t) in self.trains.iter().enumerate() {
            let bad = WorkListError::Pulse { train };
            let end = u64::from(t.first_pulse) + u64::from(t.pulse_count);
            if end > self.pulses.len() as u64 {
                return Err(bad);
            }
            let pulses = &self.pulses[t.first_pulse as usize..end as usize];
            for (i, p) in pulses.iter().enumerate() {
                let n = tables.headers.len();
                if p.start.abs() >= MAX_COORD_BITS
                    || p.fall_start.abs() >= MAX_COORD_BITS
                    || p.fall_start < p.start
                    || p.rise_table as usize >= n
                    || p.fall_table as usize >= n
                    || tables.headers[p.rise_table as usize].plateau
                        != tables.headers[p.fall_table as usize].plateau
                    || !p.amplitude.is_finite()
                    || p.amplitude.abs() > MAX_PULSE_AMPLITUDE
                    || (i > 0 && pulses[i - 1].start > p.start)
                {
                    return Err(bad);
                }
            }
        }
        Ok(())
    }

    /// SHA-256 of a canonical little-endian encoding of the whole work list,
    /// including the referenced tables. Used for tap-isolation and cache checks.
    pub fn content_hash(&self) -> [u8; 32] {
        let mut h = Sha256::new();
        h.update(b"iqsim.worklist.v2\0");
        h.update(self.target.response);
        h.update(self.target.tap_kind.to_le_bytes());
        h.update(self.target.lsb_per_volt.to_le_bytes());
        for v in [self.block, self.first_sample, self.len] {
            h.update(v.to_le_bytes());
        }
        h.update((self.paths.len() as u64).to_le_bytes());
        for p in &self.paths {
            for r in [&p.delay, &p.coord, &p.phase] {
                for w in r.words {
                    h.update(w.to_le_bytes());
                }
            }
            h.update(p.block.to_le_bytes());
            for c in p.amplitude.c {
                h.update(c.to_le_bytes());
            }
            h.update(p.amplitude.inv_len.to_le_bytes());
            for v in [
                p.target,
                p.coord_units,
                p.coord_source,
                p.emission.kind,
                p.emission.index,
                p.first,
                p.len,
            ] {
                h.update(v.to_le_bytes());
            }
        }
        h.update((self.cw.len() as u64).to_le_bytes());
        for c in &self.cw {
            h.update(c.envelope[0].to_le_bytes());
            h.update(c.envelope[1].to_le_bytes());
        }
        h.update((self.trains.len() as u64).to_le_bytes());
        for t in &self.trains {
            h.update(t.first_pulse.to_le_bytes());
            h.update(t.pulse_count.to_le_bytes());
        }
        h.update((self.pulses.len() as u64).to_le_bytes());
        for p in &self.pulses {
            h.update(p.start.to_le_bytes());
            h.update(p.fall_start.to_le_bytes());
            h.update(p.amplitude.to_le_bytes());
            h.update(p.rise_table.to_le_bytes());
            h.update(p.fall_table.to_le_bytes());
        }
        h.update((self.continuous.len() as u64).to_le_bytes());
        for c in &self.continuous {
            h.update(c.window.to_le_bytes());
        }
        h.update((self.windows.headers.len() as u64).to_le_bytes());
        for w in &self.windows.headers {
            for v in [w.first, w.len, w.log2_up] {
                h.update(v.to_le_bytes());
            }
        }
        for p in &self.windows.samples {
            h.update(p[0].to_le_bytes());
            h.update(p[1].to_le_bytes());
        }
        h.update((self.tables.headers.len() as u64).to_le_bytes());
        for t in &self.tables.headers {
            h.update(t.u_lo.to_le_bytes());
            for v in [t.first, t.len, t.log2_res] {
                h.update(v.to_le_bytes());
            }
        }
        for p in &self.tables.points {
            h.update(p[0].to_le_bytes());
            h.update(p[1].to_le_bytes());
        }
        h.finalize().into()
    }
}

/// Whether the coordinate polynomial with registers `d` (Q32.96, read as signed) is
/// non-decreasing over local indices `[0, len)` without wrapping.
///
/// Its exact value `D0 + j·D1 + C(j,2)·D2 + C(j,3)·D3` must fit i128 at both ends, and
/// the step `Δ1(j) = D1 + j·D2 + C(j,2)·D3`, a quadratic in j, must be ≥ 0 for every
/// j ∈ [0, len − 2]. The quadratic's minimum over the integers is at an end, or, when
/// D3 > 0, at one of the two integers around its vertex `j* = ½ − D2/D3`. Every
/// product and sum is checked, and an overflow rejects the polynomial. A value that
/// is monotone and wrap-free between in-range endpoints never wraps in between.
pub fn coordinate_is_monotone(d: &[U128; 4], len: u64) -> bool {
    fn exact(d: &[i128; 4], j: i128) -> Option<i128> {
        let c2 = j.checked_mul(j - 1)? / 2;
        let c3 = c2.checked_mul(j - 2)? / 3;
        d[0].checked_add(j.checked_mul(d[1])?)?
            .checked_add(c2.checked_mul(d[2])?)?
            .checked_add(c3.checked_mul(d[3])?)
    }
    fn step(d: &[i128; 4], j: i128) -> Option<i128> {
        let c2 = j.checked_mul(j - 1)? / 2;
        d[1].checked_add(j.checked_mul(d[2])?)?
            .checked_add(c2.checked_mul(d[3])?)
    }
    if len == 0 {
        return false;
    }
    let d = d.map(|r| ((u128::from(r.hi) << 64) | u128::from(r.lo)) as i128);
    let last = i128::from(len - 1);
    if exact(&d, 0).is_none() || exact(&d, last).is_none() {
        return false;
    }
    if len == 1 {
        return true;
    }
    let top = last - 1;
    let mut candidates = vec![0, top];
    if d[3] > 0 {
        // ⌊j*⌋ = ⌊(D3 − 2·D2) / (2·D3)⌋.
        let Some(num) = d[2].checked_mul(2).and_then(|t| d[3].checked_sub(t)) else {
            return false;
        };
        let Some(den) = d[3].checked_mul(2) else {
            return false;
        };
        let jc = num.div_euclid(den);
        for j in [jc, jc.saturating_add(1)] {
            candidates.push(j.clamp(0, top));
        }
    }
    candidates
        .into_iter()
        .all(|j| step(&d, j).is_some_and(|v| v >= 0))
}

/// The [`TargetDesc`] for a response and scale.
pub fn target_desc(kind: TapKind, response: TapResponseId, lsb_per_volt: f32) -> TargetDesc {
    TargetDesc {
        response: response.0,
        tap_kind: kind.code(),
        lsb_per_volt,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::mem::{align_of, size_of};
    use iqsim_core::{DelayPoly, PhasePoly};

    #[test]
    fn layouts_are_padding_free_and_stable() {
        // Sizes equal the sum of the declared fields, so there is no implicit padding
        // a GPU struct could disagree on.
        assert_eq!(size_of::<PolyRegisters>(), 64);
        assert_eq!(size_of::<AmplitudePoly>(), 20);
        assert_eq!(size_of::<EmissionRef>(), 8);
        assert_eq!(
            size_of::<PathSegment>(),
            3 * 64 + 8 + 20 + 3 * 4 + 8 + 2 * 4
        );
        assert_eq!(align_of::<PathSegment>(), 8);
        assert_eq!(size_of::<CwEmission>(), 8);
        assert_eq!(size_of::<PulseTrain>(), 8);
        assert_eq!(size_of::<Pulse>(), 8 + 8 + 4 * 4);
        assert_eq!(size_of::<EdgeTableHeader>(), 8 + 4 * 4 + 8 + 4 + 4);
        assert_eq!(size_of::<TargetDesc>(), 32 + 4 + 4);
        assert_eq!(size_of::<ContinuousEmission>(), 8);
        assert_eq!(size_of::<SourceWindowHeader>(), 6 * 4);
    }

    #[test]
    fn registers_round_trip() {
        let p = DelayPoly::from_f64_differences([123.456, -1e-3, 2e-9, -3e-15], 1000).unwrap();
        let r = PolyRegisters::from_poly(&p);
        assert_eq!(r.to_poly::<iqsim_core::Delay>(1000).unwrap(), p);
        let q = PhasePoly::from_f64_differences([0.25, 0.125, 1e-6, 0.0], 77).unwrap();
        assert_eq!(
            PolyRegisters::from_poly(&q)
                .to_poly::<iqsim_core::Phase>(77)
                .unwrap(),
            q
        );
    }

    #[test]
    fn amplitude_is_chunk_independent_and_bounded() {
        let a = AmplitudePoly {
            c: [1.0, -0.5, 0.25, 0.125],
            inv_len: (1.0 / 1000.0_f64) as f32,
        };
        assert_eq!(a.eval(0), 1.0);
        for k in 0..1000 {
            assert!(f64::from(a.eval(k).abs()) <= a.bound());
        }
        assert_eq!(AmplitudePoly::constant(2.0, 10).eval(9), 2.0);
    }

    #[test]
    fn pulse_extent_covers_a_rise_that_outlasts_the_fall() {
        // Rise table spans [0, 10] samples, fall table [0, 1].
        let mut a = EdgeTableArena::default();
        let mut pts = vec![[0.0_f32, 0.0_f32]; 11];
        pts[10] = [1.0, 0.0];
        let rise = a.push(0, 0, &pts).unwrap();
        let fall = a.push(0, 0, &[[0.0, 0.0], [1.0, 0.0]]).unwrap();
        let p = Pulse {
            start: 0,
            fall_start: 1 << 32,
            amplitude: 1.0,
            rise_table: rise,
            fall_table: fall,
            _pad: 0,
        };
        assert_eq!(a.pulse_extent(&p), (0, 10 << 32));
    }

    #[test]
    fn arena_and_validation() {
        let mut arena = EdgeTableArena::default();
        assert_eq!(arena.push(0, 4, &[[0.0, 0.0]]), Err(WorkListError::Table));
        let t = arena
            .push(-(5 << 32), 4, &[[0.0, 0.0], [0.5, 0.0], [1.0, 0.0]])
            .unwrap();
        assert_eq!(arena.headers[0].peak, 1.0);
        assert_eq!(arena.u_hi(t), -(5 << 32) + (2 << 28));
        let len = 100;
        let coord = DelayPoly::from_f64_differences([-50.0, 1.0, 0.0, 0.0], len).unwrap();
        let zero = DelayPoly::from_f64_differences([0.0; 4], len).unwrap();
        let phase = PhasePoly::from_f64_differences([0.0; 4], len).unwrap();
        let mut wl = BlockWorkList {
            target: TargetDesc {
                response: [0; 32],
                tap_kind: 2,
                lsb_per_volt: 1.0,
            },
            block: 3,
            first_sample: 300,
            len,
            paths: vec![PathSegment {
                delay: PolyRegisters::from_poly(&zero),
                coord: PolyRegisters::from_poly(&coord),
                phase: PolyRegisters::from_poly(&phase),
                block: 3,
                amplitude: AmplitudePoly::constant(1.0, len),
                target: 0,
                coord_units: COORD_TARGET_SAMPLES,
                coord_source: 0,
                emission: EmissionRef {
                    kind: EMISSION_PULSE_TRAIN,
                    index: 0,
                },
                first: 0,
                len: len as u32,
            }],
            cw: vec![],
            trains: vec![PulseTrain {
                first_pulse: 0,
                pulse_count: 1,
            }],
            pulses: vec![Pulse {
                start: 0,
                fall_start: 10 << 32,
                amplitude: 1.0,
                rise_table: 0,
                fall_table: 0,
                _pad: 0,
            }],
            continuous: vec![],
            tables: Arc::new(arena),
            windows: Arc::new(SourceWindowArena::default()),
        };
        wl.validate().unwrap();
        let h = wl.content_hash();
        wl.pulses[0].fall_table = 1;
        assert_eq!(wl.validate(), Err(WorkListError::Pulse { train: 0 }));
        wl.pulses[0].fall_table = 0;
        wl.paths[0].block = 4;
        assert_eq!(wl.validate(), Err(WorkListError::Path { index: 0 }));
        wl.paths[0].block = 3;
        // A coordinate that turns back is rejected before any endpoint check: endpoint
        // checks alone miss an interior wrap such as w(k) = 2²⁹k(4 − k) over 5 samples.
        let curved = DelayPoly::from_f64_differences(
            [0.0, 3.0 * 536_870_912.0, -2.0 * 536_870_912.0, 0.0],
            len,
        )
        .unwrap();
        let saved = wl.paths[0].coord;
        wl.paths[0].coord = PolyRegisters::from_poly(&curved);
        assert_eq!(wl.validate(), Err(WorkListError::Path { index: 0 }));
        wl.paths[0].coord = saved;
        // Segment ranges: empty or past the block end are rejected; a sub-range is fine.
        wl.paths[0].len = 0;
        assert_eq!(wl.validate(), Err(WorkListError::Path { index: 0 }));
        wl.paths[0].first = 1;
        wl.paths[0].len = len as u32;
        assert_eq!(wl.validate(), Err(WorkListError::Path { index: 0 }));
        wl.paths[0].len = len as u32 - 1;
        // The registers were built for `len` samples; over `len − 1` they still are a
        // valid, monotone polynomial.
        wl.validate().unwrap();
        wl.paths[0].first = 0;
        wl.paths[0].len = len as u32;
        wl.pulses[0].amplitude = 2e6;
        assert_eq!(wl.validate(), Err(WorkListError::Pulse { train: 0 }));
        wl.pulses[0].amplitude = 1.0;
        wl.validate().unwrap();
        wl.paths[0].emission.kind = EMISSION_CW;
        assert_eq!(wl.validate(), Err(WorkListError::Path { index: 0 }));
        assert_ne!(wl.content_hash(), h);
    }

    #[test]
    fn monotone_cubic_coordinates_are_accepted_and_interior_dips_are_not() {
        let regs =
            |d: [f64; 4], len: u64| DelayPoly::from_f64_differences(d, len).unwrap().registers();
        // Linear, and a gently curving cubic (geometry fits look like this).
        assert!(coordinate_is_monotone(
            &regs([5.0, 1.0, 0.0, 0.0], 1000),
            1000
        ));
        assert!(coordinate_is_monotone(
            &regs([5.0, 1.0, 1e-7, -1e-10], 1000),
            1000
        ));
        // Step 1 − 0.002·j + 1e-6·C(j,2): its vertex is near j = 2001, where the step is
        // ≈ −1.0. Both ends are increasing, so only the vertex test catches it.
        let dip = regs([0.0, 1.0, -0.002, 1e-6], 4000);
        assert!(!coordinate_is_monotone(&dip, 4000));
        // Over the first 100 samples the same polynomial is still increasing.
        assert!(coordinate_is_monotone(&dip, 100));
        // Decreasing at the start (concave), and a value that wraps i128.
        assert!(!coordinate_is_monotone(
            &regs([0.0, -1e-9, 0.0, 0.0], 10),
            10
        ));
        let huge = [U128::ZERO, U128::new(1 << 62, 0), U128::ZERO, U128::ZERO];
        assert!(!coordinate_is_monotone(&huge, 4));
        assert!(coordinate_is_monotone(&huge, 2));
        assert!(coordinate_is_monotone(&regs([1.0, -5.0, 0.0, 0.0], 1), 1));
        assert!(!coordinate_is_monotone(&regs([0.0; 4], 1), 0));
    }

    #[test]
    fn continuous_paths_must_keep_every_farrow_tap_in_their_window() {
        let len = 100;
        let mut windows = SourceWindowArena::default();
        assert_eq!(
            windows.push(2, &[[0.0, 0.0]; 5]),
            Err(WorkListError::Window)
        );
        assert_eq!(
            windows.push(6, &[[0.0, 0.0]; 64]),
            Err(WorkListError::Window)
        );
        assert_eq!(
            windows.push(2, &[[f32::NAN, 0.0]; 64]),
            Err(WorkListError::Window)
        );
        // 64 points at U = 4: native coordinates 0 .. 16.
        let win = windows.push(2, &[[0.5, -0.5]; 64]).unwrap();
        // The stored peak is the f64 magnitude rounded up to an f32 bound.
        let p = windows.headers[0].peak;
        assert!(f64::from(p) >= 0.5_f64.sqrt() && f64::from(p) - 0.5_f64.sqrt() < 1e-7);
        // Tiny values keep a non-zero peak (review fix: f32 squares underflowed).
        let mut tiny = SourceWindowArena::default();
        tiny.push(0, &[[1e-30, -1e-30]; 8]).unwrap();
        assert!(tiny.headers[0].peak >= 1.4e-30);
        // Coordinate 0.5 + 0.1·k native samples: cells 2 .. 41, taps 0 .. 44.
        let path = |w0: f64| {
            let coord = DelayPoly::from_f64_differences([w0, 0.1, 0.0, 0.0], len).unwrap();
            let zero = DelayPoly::from_f64_differences([0.0; 4], len).unwrap();
            let phase = PhasePoly::from_f64_differences([0.0; 4], len).unwrap();
            PathSegment {
                delay: PolyRegisters::from_poly(&zero),
                coord: PolyRegisters::from_poly(&coord),
                phase: PolyRegisters::from_poly(&phase),
                block: 0,
                amplitude: AmplitudePoly::constant(1.0, len),
                target: 0,
                coord_units: COORD_NATIVE_SAMPLES,
                coord_source: 0,
                emission: EmissionRef {
                    kind: EMISSION_CONTINUOUS,
                    index: 0,
                },
                first: 0,
                len: len as u32,
            }
        };
        let mut wl = BlockWorkList {
            target: TargetDesc {
                response: [0; 32],
                tap_kind: 2,
                lsb_per_volt: 1.0,
            },
            block: 0,
            first_sample: 0,
            len,
            paths: vec![path(0.5)],
            cw: vec![],
            trains: vec![],
            pulses: vec![],
            continuous: vec![ContinuousEmission {
                window: win,
                _pad: 0,
            }],
            tables: Arc::new(EdgeTableArena::default()),
            windows: Arc::new(windows),
        };
        wl.validate().unwrap();
        // Too early: cell 1 needs tap −1.
        wl.paths[0] = path(0.25);
        assert_eq!(wl.validate(), Err(WorkListError::Path { index: 0 }));
        // Too late: 6.2 + 9.9 = 16.1 native samples runs past the window.
        wl.paths[0] = path(6.2);
        assert_eq!(wl.validate(), Err(WorkListError::Path { index: 0 }));
        // Just inside: last cell 60, taps up to 63.
        wl.paths[0] = path(5.1);
        wl.validate().unwrap();
        // Target-sample units are wrong for a continuous emission.
        wl.paths[0].coord_units = COORD_TARGET_SAMPLES;
        assert_eq!(wl.validate(), Err(WorkListError::Path { index: 0 }));
        wl.paths[0].coord_units = COORD_NATIVE_SAMPLES;
        // Review regression: a slope that wraps the 128-bit register over the block
        // ends back inside the window but is rejected.
        let mut regs = DelayPoly::from_f64_differences([0.5, 0.0, 0.0, 0.0], len)
            .unwrap()
            .registers();
        // D1 = 2¹²⁸/99 + ε: 99 steps wrap once and land just past the start.
        regs[1] = U128::new(0x0295_FAD4_0A57_EB50, 0x295F_AD40_A57E_B504);
        let wrapped = DelayPoly::from_registers(regs, len).unwrap();
        let saved = wl.paths[0].coord;
        wl.paths[0].coord = PolyRegisters::from_poly(&wrapped);
        assert_eq!(wl.validate(), Err(WorkListError::Path { index: 0 }));
        wl.paths[0].coord = saved;
        wl.validate().unwrap();
        wl.continuous[0].window = 1;
        assert_eq!(wl.validate(), Err(WorkListError::Path { index: 0 }));
    }
}
