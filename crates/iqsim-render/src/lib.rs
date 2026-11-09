//! Render backend interface and CPU reference backend.
//! Per the plan section "Compute backends and scaling", this crate consumes flat
//! block work lists with stateless kernels and checked i64 accumulation.
//! CPU output supplies the bit-exact reference for every alternate backend.
//!
//! # Per-sample kernel (contracts §5.2)
//!
//! For each path and local sample k, in this order, with every f32 operation separate:
//!
//! 1. `a = amplitude.eval(k)` (volts, Horner from scratch).
//! 2. Envelope `(er, ei)`: the CW envelope, or for a pulse train
//!    `Σ_i A_i·(G_rise(w − start_i) − G_fall(w − fall_start_i))` over the train's pulses
//!    in order, with w the waveform coordinate at k. Each G is a linear interpolation
//!    of its edge table ([`edge_value`]).
//! 3. `(s, c) = sincos(phase(k))`.
//! 4. `pr = er·c − ei·s`, `pi = er·s + ei·c`.
//! 5. `m = a·lsb_per_volt`, `I = pr·m`, `Q = pi·m`.
//! 6. Each component rounds ties-to-even to i64 and is added with a checked add.
//!
//! The polynomials jump to the chunk start and then step, and jump and step agree bit
//! for bit. Every other quantity is computed from k or w alone. Every per-sample
//! contribution is therefore independent of chunking, and integer addition makes the sum
//! independent of path order and thread count.
//!
//! Dependencies: iqsim-core, iqsim-worklist, iqsim-waveform.

#![forbid(unsafe_code)]

use std::fmt;

use iqsim_core::{Delay, NumericError, Phase, f32_to_i64_ties_even, sincos};
use iqsim_worklist::{
    BlockWorkList, EMISSION_CW, EMISSION_PULSE_TRAIN, EdgeTableArena, PathSegment, Pulse,
    WorkListError,
};

/// A render failure. Every variant is a hard error: no partial block is valid.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum RenderError {
    /// The work list failed validation.
    WorkList(WorkListError),
    /// The accumulator slice is not `2 · len` long.
    AccumulatorLength {
        /// Expected length (I/Q interleaved).
        expected: usize,
        /// Supplied length.
        actual: usize,
    },
    /// A contribution was not finite or does not fit i64.
    Contribution {
        /// Local sample index.
        k: u64,
    },
    /// An accumulator add overflowed i64 (the headroom proof was violated).
    Overflow {
        /// Local sample index.
        k: u64,
    },
    /// A polynomial operation failed.
    Numeric(NumericError),
}

impl fmt::Display for RenderError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::WorkList(e) => write!(f, "invalid work list: {e}"),
            Self::AccumulatorLength { expected, actual } => {
                write!(f, "accumulator holds {actual} values, expected {expected}")
            }
            Self::Contribution { k } => {
                write!(f, "contribution at sample {k} is not finite or exceeds i64")
            }
            Self::Overflow { k } => write!(f, "accumulator overflow at sample {k}"),
            Self::Numeric(e) => write!(f, "numeric error: {e}"),
        }
    }
}

impl std::error::Error for RenderError {}

impl From<NumericError> for RenderError {
    fn from(e: NumericError) -> Self {
        Self::Numeric(e)
    }
}

/// A render backend: adds a block's contributions into a target's accumulators.
pub trait Backend {
    /// Adds every path of `work` into `acc`, I/Q interleaved (`acc[2k]` = I,
    /// `acc[2k+1]` = Q, `len = 2·work.len`). The caller owns zeroing and reduction.
    fn render(&self, work: &BlockWorkList, acc: &mut [i64]) -> Result<(), RenderError>;
}

/// Value of an edge table at coordinate `u` (Q32.32 target samples, relative to the
/// edge start): 0 below the table, the plateau from the last point on, linear
/// interpolation in between. The interpolation weight is the top ≤ 24 bits of the
/// fractional cell position, so its conversion to f32 is exact.
#[inline]
pub fn edge_value(tables: &EdgeTableArena, table: u32, u: i64) -> [f32; 2] {
    let h = &tables.headers[table as usize];
    // |u| < 2⁶³ and |u_lo| < 2⁶², so the difference can leave i64 only far outside the
    // table, where saturation gives the same answer (0 below, plateau above).
    let d = u.saturating_sub(h.u_lo);
    if d < 0 {
        return [0.0, 0.0];
    }
    let frac_bits = 32 - h.log2_res;
    let pos = (d >> frac_bits) as u64;
    if pos >= u64::from(h.len) - 1 {
        return h.plateau;
    }
    let frac = (d as u64) & ((1u64 << frac_bits) - 1);
    let w = if frac_bits > 24 {
        ((frac >> (frac_bits - 24)) as f32) * (1.0 / 16_777_216.0)
    } else {
        (frac as f32) * (1.0 / (1u64 << frac_bits) as f32)
    };
    let i = h.first as usize + pos as usize;
    let g0 = tables.points[i];
    let g1 = tables.points[i + 1];
    let dr = g1[0] - g0[0];
    let di = g1[1] - g0[1];
    [g0[0] + dr * w, g0[1] + di * w]
}

/// The multithreaded CPU golden-reference backend.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct CpuBackend {
    threads: usize,
    chunk_len: Option<u64>,
}

impl CpuBackend {
    /// A backend using `threads` worker threads (at least 1). The block is split into
    /// `threads` contiguous chunks.
    pub fn new(threads: usize) -> Self {
        Self {
            threads: threads.max(1),
            chunk_len: None,
        }
    }

    /// Uses chunks of `chunk_len` samples (at least 1), dealt round-robin to the
    /// threads. Output does not depend on it; tests use odd sizes to show that.
    pub fn with_chunk_len(mut self, chunk_len: u64) -> Self {
        self.chunk_len = Some(chunk_len.max(1));
        self
    }
}

/// One path's state prepared for a chunk.
struct PathChunk<'a> {
    path: &'a PathSegment,
    /// Pulses of the train that can be non-zero inside the chunk (empty for CW).
    pulses: Vec<&'a Pulse>,
    cw: [f32; 2],
    is_cw: bool,
}

fn render_chunk(
    work: &BlockWorkList,
    k0: u64,
    acc: &mut [i64],
    lsb_per_volt: f32,
) -> Result<(), RenderError> {
    let n = (acc.len() / 2) as u64;
    let k1 = k0 + n;
    let tables = &*work.tables;
    for path in &work.paths {
        let coord = path.coord.to_poly::<Delay>(work.len)?;
        let phase = path.phase.to_poly::<Phase>(work.len)?;
        let mut pc = PathChunk {
            path,
            pulses: Vec::new(),
            cw: [0.0, 0.0],
            is_cw: path.emission.kind == EMISSION_CW,
        };
        if pc.is_cw {
            pc.cw = work.cw[path.emission.index as usize].envelope;
        } else {
            debug_assert_eq!(path.emission.kind, EMISSION_PULSE_TRAIN);
            // The coordinate is non-decreasing over a block (builder contract), so the
            // chunk spans [w(k0), w(k1 − 1)].
            let w_lo = coord.eval(k0)?.to_bits();
            let w_hi = coord.eval(k1 - 1)?.to_bits();
            let train = work.trains[path.emission.index as usize];
            let first = train.first_pulse as usize;
            for p in &work.pulses[first..first + train.pulse_count as usize] {
                let (on, off) = tables.pulse_extent(p);
                if on <= w_hi && off >= w_lo {
                    pc.pulses.push(p);
                }
            }
            if pc.pulses.is_empty() {
                continue;
            }
        }
        let mut coord_c = coord.cursor(k0)?;
        let mut phase_c = phase.cursor(k0)?;
        for i in 0..n as usize {
            let k = k0 + i as u64;
            let a = pc.path.amplitude.eval(k);
            let env = if pc.is_cw {
                pc.cw
            } else {
                let w = coord_c.value().to_bits();
                let mut e = [0.0_f32, 0.0_f32];
                for p in &pc.pulses {
                    let r = edge_value(tables, p.rise_table, w - p.start);
                    let f = edge_value(tables, p.fall_table, w - p.fall_start);
                    let dr = r[0] - f[0];
                    let di = r[1] - f[1];
                    e[0] += p.amplitude * dr;
                    e[1] += p.amplitude * di;
                }
                e
            };
            let (s, c) = sincos(phase_c.value());
            let pr = env[0] * c - env[1] * s;
            let pi = env[0] * s + env[1] * c;
            let m = a * lsb_per_volt;
            let vi = f32_to_i64_ties_even(pr * m).map_err(|_| RenderError::Contribution { k })?;
            let vq = f32_to_i64_ties_even(pi * m).map_err(|_| RenderError::Contribution { k })?;
            let slot = &mut acc[2 * i..2 * i + 2];
            slot[0] = slot[0].checked_add(vi).ok_or(RenderError::Overflow { k })?;
            slot[1] = slot[1].checked_add(vq).ok_or(RenderError::Overflow { k })?;
            coord_c.step();
            phase_c.step();
        }
    }
    Ok(())
}

impl Backend for CpuBackend {
    fn render(&self, work: &BlockWorkList, acc: &mut [i64]) -> Result<(), RenderError> {
        work.validate().map_err(RenderError::WorkList)?;
        let expected = 2 * work.len as usize;
        if acc.len() != expected {
            return Err(RenderError::AccumulatorLength {
                expected,
                actual: acc.len(),
            });
        }
        let chunk = self
            .chunk_len
            .unwrap_or_else(|| work.len.div_ceil(self.threads as u64))
            .min(work.len);
        let lsb = work.target.lsb_per_volt;
        let mut chunks: Vec<(u64, &mut [i64])> = acc
            .chunks_mut(2 * chunk as usize)
            .enumerate()
            .map(|(i, c)| (i as u64 * chunk, c))
            .collect();
        let threads = self.threads.min(chunks.len());
        if threads <= 1 {
            for (k0, c) in chunks {
                render_chunk(work, k0, c, lsb)?;
            }
            return Ok(());
        }
        // Deal chunks round-robin; report the error of the lowest failing chunk so the
        // outcome does not depend on scheduling.
        let mut buckets: Vec<Vec<(u64, &mut [i64])>> = (0..threads).map(|_| Vec::new()).collect();
        for (i, c) in chunks.drain(..).enumerate() {
            buckets[i % threads].push(c);
        }
        let results: Vec<Vec<(u64, Result<(), RenderError>)>> = std::thread::scope(|s| {
            let handles: Vec<_> = buckets
                .into_iter()
                .map(|bucket| {
                    s.spawn(move || {
                        bucket
                            .into_iter()
                            .map(|(k0, c)| (k0, render_chunk(work, k0, c, lsb)))
                            .collect::<Vec<_>>()
                    })
                })
                .collect();
            handles
                .into_iter()
                .map(|h| h.join().expect("render worker panicked"))
                .collect()
        });
        let mut errors: Vec<(u64, RenderError)> = results
            .into_iter()
            .flatten()
            .filter_map(|(k0, r)| r.err().map(|e| (k0, e)))
            .collect();
        errors.sort_by_key(|(k0, _)| *k0);
        match errors.first() {
            Some((_, e)) => Err(*e),
            None => Ok(()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn edge_value_regions_and_interpolation() {
        let mut t = EdgeTableArena::default();
        // Points at u = 2, 2.25, 2.5 (log2_res 2).
        let id = t
            .push(2 << 32, 2, &[[0.0, 0.0], [1.0, -1.0], [2.0, 0.5]])
            .unwrap();
        assert_eq!(edge_value(&t, id, (2 << 32) - 1), [0.0, 0.0]);
        assert_eq!(edge_value(&t, id, 2 << 32), [0.0, 0.0]);
        assert_eq!(edge_value(&t, id, (2 << 32) + (1 << 29)), [0.5, -0.5]);
        assert_eq!(edge_value(&t, id, (2 << 32) + (1 << 30)), [1.0, -1.0]);
        assert_eq!(edge_value(&t, id, (2 << 32) + (1 << 31)), [2.0, 0.5]);
        assert_eq!(edge_value(&t, id, 100 << 32), [2.0, 0.5]);
    }
}
