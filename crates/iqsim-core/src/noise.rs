//! Complex Gaussian noise from Philox4x32-10 by Box-Muller (contracts §1.5).
//!
//! `noise(seed, entity, stream, n)` is circularly symmetric complex Gaussian with
//! **unit complex power**: `E[I] = E[Q] = 0`, `E[I²] = E[Q²] = ½`, `E[I·Q] = 0` and
//! `E[|z|²] = 1`. Scale it by the RMS voltage of the wanted noise power (contracts §2).
//!
//! Sample `n` uses one Philox block, at counter (n, [`NOISE_SUB_INDEX`]):
//!
//! 1. Words 0 and 1 form `x = w1·2³² + w0`. The radius uniform is `u = x·2⁻⁶⁴`
//!    truncated to 24 significant bits, built as f32 from the bit pattern with integer
//!    operations only (exact). `x = 0` maps to `u = 2⁻⁶⁵`, so `u ∈ [2⁻⁶⁵, 1 − 2⁻²⁴]`.
//! 2. `r = sqrt(−ln u)`, with the in-house [`ln`](crate::ln) and IEEE sqrt. Since `−ln u`
//!    is Exp(1), `r²` has unit mean, and `r ≤ sqrt(65·ln 2) ≈ 6.71` (9.5 σ per
//!    component).
//! 3. Word 2 is the angle, as the top 32 bits of a [`Phase`]. `(s, c) = sincos(θ)`.
//! 4. The output is `(r·c, r·s)`.
//!
//! Word 3 is unused. Every step is integer or a single IEEE f32 operation, so the CUDA
//! port gives the same bits. The quantization of `u` (24 bits) and θ (2⁻³² cycle) and
//! the `ln`/`sincos` errors bias the variance by under 10⁻⁶ relative.

use crate::fixed::Phase;
use crate::log::ln_positive;
use crate::philox::{PhiloxStream, StreamId};
use crate::timebase::SampleIndex;
use crate::trig::sincos;

/// The Philox sub-index the noise draw uses. Noise takes one block per sample.
pub const NOISE_SUB_INDEX: u64 = 0;

/// The radius uniform: `x·2⁻⁶⁴` truncated to 24 significant bits, or 2⁻⁶⁵ for `x = 0`.
#[inline]
pub(crate) const fn radius_uniform(x: u64) -> f32 {
    if x == 0 {
        // 2⁻⁶⁵: biased exponent 127 − 65.
        return f32::from_bits(62 << 23);
    }
    let lz = x.leading_zeros();
    // The leading one is bit p = 63 − lz, so x·2⁻⁶⁴ = 1.m · 2^(p − 64), with biased
    // exponent p − 64 + 127 = 126 − lz ∈ [63, 126].
    let m24 = ((x << lz) >> 40) as u32;
    f32::from_bits(((126 - lz) << 23) | (m24 & 0x007F_FFFF))
}

/// The Box-Muller transform of one Philox block (module docs, steps 1–4).
#[inline]
pub(crate) fn gaussian_from_block(w: [u32; 4]) -> (f32, f32) {
    let x = ((w[1] as u64) << 32) | w[0] as u64;
    let neg_ln = -ln_positive(radius_uniform(x));
    let r = neg_ln.sqrt();
    let (s, c) = sincos(Phase::from_bits((w[2] as u64) << 32));
    (r * c, r * s)
}

/// A keyed noise stream for one (seed, entity, stream). The key is mixed once, here.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct NoiseStream {
    philox: PhiloxStream,
}

impl NoiseStream {
    /// The noise stream for (scenario seed, entity id, stream id).
    pub const fn new(seed: u64, entity: u32, stream: StreamId) -> Self {
        Self {
            philox: PhiloxStream::new(seed, entity, stream),
        }
    }

    /// Unit-power complex Gaussian sample `(I, Q)` at sample index `n`.
    #[inline]
    pub fn sample(self, n: SampleIndex) -> (f32, f32) {
        gaussian_from_block(self.philox.block(n.get(), NOISE_SUB_INDEX))
    }
}

/// Unit-power complex Gaussian sample `(I, Q)` for (seed, entity, stream) at sample
/// index `n`. Equal to `NoiseStream::new(seed, entity, stream).sample(n)`. Use the
/// stream form in loops, so the key is mixed once.
pub fn noise(seed: u64, entity: u32, stream: StreamId, n: SampleIndex) -> (f32, f32) {
    NoiseStream::new(seed, entity, stream).sample(n)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 2ᵉ as f64, from the bit pattern (normal range only).
    fn pow2(e: i64) -> f64 {
        f64::from_bits(((1023 + e) as u64) << 52)
    }

    /// The exact value the radius uniform must take (f64 holds it exactly).
    fn want_uniform(x: u64) -> f64 {
        if x == 0 {
            return pow2(-65);
        }
        let p = 63 - x.leading_zeros();
        let keep = if p >= 23 {
            (x >> (p - 23)) << (p - 23)
        } else {
            x
        };
        // keep has ≤ 24 significant bits, so the conversion and the scaling are exact.
        keep as f64 * pow2(-64)
    }

    #[test]
    fn radius_uniform_is_truncated_value() {
        let mut xs = vec![
            0,
            1,
            2,
            3,
            (1 << 23) - 1,
            1 << 23,
            (1 << 24) + 1,
            u64::MAX,
            u64::MAX - 1,
        ];
        xs.extend((0..64).map(|p| 1u64 << p));
        xs.extend((1..64).map(|p| (1u64 << p) - 1));
        let mut z = 0x1234_5678u64;
        for _ in 0..100_000 {
            z = crate::philox::mix64(z.wrapping_add(0x9E37_79B9_7F4A_7C15));
            xs.push(z >> (z % 64));
        }
        for x in xs {
            let u = radius_uniform(x);
            assert_eq!(u as f64, want_uniform(x), "x = {x:#x}");
            assert!(u > 0.0 && u < 1.0);
        }
        assert_eq!(radius_uniform(u64::MAX) as f64, 1.0 - pow2(-24));
    }

    #[test]
    fn radius_uniform_is_monotone() {
        let mut prev = radius_uniform(0);
        for p in 0..64 {
            for x in [
                1u64 << p,
                (1u64 << p) | ((1u64 << p) >> 1),
                ((1u64 << p) - 1) | (1u64 << p),
            ] {
                let u = radius_uniform(x);
                assert!(u >= prev, "x = {x:#x}");
                prev = u;
            }
        }
    }

    #[test]
    fn output_bounds_at_extremes() {
        // x = 0 gives the largest radius; x = MAX the smallest.
        let (i, q) = gaussian_from_block([0, 0, 0, 0]);
        let r_max = ((i as f64) * (i as f64) + (q as f64) * (q as f64)).sqrt();
        assert!(r_max > 6.70 && r_max < 6.72, "{r_max}");
        let (i, q) = gaussian_from_block([u32::MAX, u32::MAX, 0, 0]);
        assert!(i > 0.0 && i < 3e-4 && q == 0.0, "{i} {q}");
        // Angle a quarter cycle: I ≈ 0, Q = r.
        let (i, q) = gaussian_from_block([0, 1 << 31, 1 << 30, 0]);
        assert!(i.abs() < 1e-6 && q > 0.8 && q < 0.9, "{i} {q}");
    }

    #[test]
    fn free_function_matches_stream() {
        let id = StreamId::receiver_thermal_noise(2);
        let s = NoiseStream::new(99, 4, id);
        for n in [0, 1, 1 << 32, u64::MAX] {
            let n = SampleIndex::new(n);
            assert_eq!(noise(99, 4, id, n), s.sample(n));
        }
    }
}
