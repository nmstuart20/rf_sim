//! Deterministic sin/cos of a [`Phase`] (contracts §1.5).
//!
//! **Method.** The top [`SIN_TABLE_BITS`] = 10 phase bits select a table angle
//! θ = 2π·i/1024. The next 24 bits give the remainder δ ∈ [0, 2π/1024) as an f32
//! (exactly representable integer times a constant). The remaining 30 low bits are
//! dropped, an angle error below 2π·2⁻³⁴ ≈ 3.7·10⁻¹⁰ rad. Then
//!
//! ```text
//! sin δ ≈ sd = δ − δ·(δ²·(1/6))      (truncation < δ⁵/120 < 7.4·10⁻¹⁴)
//! 1 − cos δ ≈ h = δ²·0.5             (truncation < δ⁴/24  < 5.9·10⁻¹¹)
//! sin(θ+δ) = s + (c·sd − s·h)
//! cos(θ+δ) = c − (s·sd + c·h)
//! ```
//!
//! where s = sin θ and c = cos θ come from one 1024-entry f32 sine table (cos θ is the
//! entry a quarter cycle ahead). Every f32 operation is written separately, in source
//! order, with no fused multiply-add. The large table term is added last, so only the
//! table rounding (≤ 2⁻²⁵) and the final add (≤ ½ ulp) contribute materially. Measured
//! maximum absolute error over a 2²⁴-point sweep is recorded in `docs/progress.md`;
//! the contract is ≤ 3·10⁻⁷.
//!
//! No subnormal value can arise: the smallest nonzero δ is ≈ 3.7·10⁻¹⁰, so δ³ ≈ 5·10⁻²⁹,
//! and every nonzero table entry is ≥ sin(2π/1024) ≈ 6.1·10⁻³.
//!
//! **Table construction** ([`SIN_TABLE`], built at compile time, integer only, no
//! floating-point library):
//!
//! 1. For j = 0..=256, x_j = (π/2)·j/256 in unsigned Q2.62, from the constant
//!    ⌊(π/2)·2⁶²⌋ = `0x6487ED5110B4611A`, multiplied by j and shifted right by 8.
//! 2. sin x_j by the Taylor series Σ (−1)ⁿ x²ⁿ⁺¹/(2n+1)!, with each term formed from the
//!    previous one as ⌊⌊term·x²/2⁶²⌋ / ((2n)(2n+1))⌋ until it is zero. Total error is a
//!    few units of 2⁻⁶².
//! 3. Round each Q2.62 value to the nearest f32, ties to even, by explicit bit
//!    manipulation. A unit test checks that no value is within 2⁸ units of a rounding
//!    tie, so the few-unit error above cannot change the result.
//! 4. Fill the table by quarter-wave symmetry: entry 256q + j is sin x_j, sin x_{256−j},
//!    −sin x_j or −sin x_{256−j} for quadrants q = 0..3. Zero is always +0.0.
//!
//! The table bytes have a golden hash (`core/sincos_table`).

use crate::fixed::Phase;

/// log₂ of the table length.
pub const SIN_TABLE_BITS: u32 = 10;

const TABLE_LEN: usize = 1 << SIN_TABLE_BITS;
const QUARTER: usize = TABLE_LEN / 4;

/// Bits of the phase remainder that are used after the table index.
const REMAINDER_BITS: u32 = 24;

/// ⌊(π/2)·2⁶²⌋.
const HALF_PI_Q62: u64 = 0x6487_ED51_10B4_611A;

/// 2π rounded to f32.
const TWO_PI_F32: f32 = 6.283_185_5;

/// Radians per remainder unit: 2π·2⁻³⁴. Dividing by a power of two is exact.
const DELTA_SCALE: f32 = TWO_PI_F32 / ((1u64 << (SIN_TABLE_BITS + REMAINDER_BITS)) as f32);

/// 1/6 rounded to f32.
const INV6: f32 = 1.0 / 6.0;

/// sin((π/2)·j/256) in Q2.62 by the integer Taylor series (step 2 of the module docs).
const fn sin_quarter_q62(j: u64) -> u64 {
    let x = ((HALF_PI_Q62 as u128) * (j as u128)) >> 8;
    let x2 = (x * x) >> 62;
    let mut sum = x;
    let mut term = x;
    let mut n: u128 = 1;
    loop {
        term = ((term * x2) >> 62) / ((2 * n) * (2 * n + 1));
        if term == 0 {
            break;
        }
        if n % 2 == 1 {
            sum -= term;
        } else {
            sum += term;
        }
        n += 1;
    }
    // sum ≤ 1 + a few units < 2⁶³, so it fits in u64.
    sum as u64
}

/// Nearest f32 to the Q2.62 value `v` (v < 2⁶³), ties to even, with the distance of
/// the discarded bits from the tie point (for the construction test).
const fn q62_to_f32(v: u64) -> (f32, u64) {
    if v == 0 {
        return (0.0, u64::MAX);
    }
    let p = 63 - v.leading_zeros(); // index of the leading one bit
    if p <= 23 {
        // Fewer than 24 significant bits: exact. Unreachable for table values (≥ 2⁵⁴).
        let m = v << (23 - p);
        let e = p as i32 - 62;
        let bits = (((e + 127) as u32) << 23) | ((m as u32) & 0x7F_FFFF);
        return (f32::from_bits(bits), u64::MAX);
    }
    let shift = p - 23;
    let mut m = v >> shift;
    let rem = v & ((1u64 << shift) - 1);
    let half = 1u64 << (shift - 1);
    if rem > half || (rem == half && (m & 1) == 1) {
        m += 1;
    }
    let mut e = p as i32 - 62;
    if m == 1 << 24 {
        m >>= 1;
        e += 1;
    }
    let bits = (((e + 127) as u32) << 23) | ((m as u32) & 0x7F_FFFF);
    let tie_distance = rem.abs_diff(half);
    (f32::from_bits(bits), tie_distance)
}

/// Negate, keeping zero as +0.0.
const fn neg_nonzero(x: f32) -> f32 {
    let b = x.to_bits();
    if b == 0 {
        x
    } else {
        f32::from_bits(b ^ 0x8000_0000)
    }
}

const fn build_sin_table() -> [f32; TABLE_LEN] {
    let mut quarter = [0.0f32; QUARTER + 1];
    let mut j = 0;
    while j <= QUARTER {
        quarter[j] = q62_to_f32(sin_quarter_q62(j as u64)).0;
        j += 1;
    }
    let mut table = [0.0f32; TABLE_LEN];
    let mut i = 0;
    while i < TABLE_LEN {
        let q = i / QUARTER;
        let j = i % QUARTER;
        table[i] = match q {
            0 => quarter[j],
            1 => quarter[QUARTER - j],
            2 => neg_nonzero(quarter[j]),
            _ => neg_nonzero(quarter[QUARTER - j]),
        };
        i += 1;
    }
    table
}

/// sin(2π·i/1024) rounded to f32, for i = 0..1024 (construction in the module docs).
/// The CUDA backend uploads these exact bytes.
pub const SIN_TABLE: [f32; TABLE_LEN] = build_sin_table();

/// `(sin, cos)` of the angle 2π·phase/2⁶⁴, as f32.
///
/// Deterministic: the same bits on every IEEE 754 platform with round-to-nearest-even
/// and no FMA contraction. Maximum absolute error ≤ 3·10⁻⁷ (contracts §1.5).
pub fn sincos(phase: Phase) -> (f32, f32) {
    let b = phase.to_bits();
    let i = (b >> (64 - SIN_TABLE_BITS)) as usize;
    // 24 remainder bits: an integer < 2²⁴, so the u32 → f32 conversion is exact.
    let r = ((b >> (64 - SIN_TABLE_BITS - REMAINDER_BITS)) & ((1 << REMAINDER_BITS) - 1)) as u32;
    let s = SIN_TABLE[i];
    let c = SIN_TABLE[(i + QUARTER) & (TABLE_LEN - 1)];

    let d = (r as f32) * DELTA_SCALE;
    let d2 = d * d;
    let d3_6 = d * (d2 * INV6);
    let sd = d - d3_6;
    let h = d2 * 0.5;

    let cs = c * sd;
    let sh = s * h;
    let sin = s + (cs - sh);
    let ss = s * sd;
    let ch = c * h;
    let cos = c - (ss + ch);
    (sin, cos)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn table_constants() {
        assert_eq!(SIN_TABLE[0].to_bits(), 0);
        assert_eq!(SIN_TABLE[256], 1.0);
        assert_eq!(SIN_TABLE[512].to_bits(), 0);
        assert_eq!(SIN_TABLE[768], -1.0);
        for j in 1..256 {
            assert_eq!(SIN_TABLE[j], SIN_TABLE[512 - j]);
            assert_eq!(SIN_TABLE[j], -SIN_TABLE[512 + j]);
            assert_eq!(SIN_TABLE[j], -SIN_TABLE[1024 - j]);
        }
    }

    #[test]
    fn no_table_value_near_a_rounding_tie() {
        for j in 1..=QUARTER as u64 {
            let (_, dist) = q62_to_f32(sin_quarter_q62(j));
            assert!(dist > 1 << 8, "j = {j}: {dist} units from a tie");
        }
    }

    // Reference oracle only: std f64 sin is accurate to ~1 ulp, far finer than the
    // f32 result being checked, and it never feeds simulation output.
    #[allow(clippy::disallowed_methods)]
    #[test]
    fn table_is_correctly_rounded() {
        for (i, &v) in SIN_TABLE.iter().enumerate() {
            let x = (i as f64) * (2.0 * std::f64::consts::PI / 1024.0);
            let want = x.sin() as f32;
            // sin(π) and sin(2π) in f64 are ~1e-16, not 0; the table holds exact zeros.
            if i % 512 == 0 {
                assert_eq!(v, 0.0);
            } else {
                assert_eq!(v, want, "entry {i}");
            }
        }
    }

    #[test]
    fn rounding_helper_edges() {
        assert_eq!(q62_to_f32(1 << 62).0, 1.0);
        assert_eq!(q62_to_f32(1 << 61).0, 0.5);
        // 1 + 2⁻²⁴ is a tie between 1 and 1 + 2⁻²³: rounds to even (1.0).
        assert_eq!(q62_to_f32((1 << 62) + (1 << 38)).0, 1.0);
        // 1 + 3·2⁻²⁴ is a tie that rounds up to 1 + 2⁻²².
        assert_eq!(
            q62_to_f32((1 << 62) + 3 * (1 << 38)).0,
            1.0 + f32::EPSILON * 2.0
        );
        // Just below 2 rounds up to 2.0 (mantissa carry into the exponent).
        assert_eq!(q62_to_f32((1 << 63) - 1).0, 2.0);
        assert_eq!(q62_to_f32(3).0, 3.0 / (1u64 << 62) as f32);
    }

    #[test]
    fn exact_quadrant_points() {
        let q = 1u64 << 62;
        assert_eq!(sincos(Phase::from_bits(0)), (0.0, 1.0));
        assert_eq!(sincos(Phase::from_bits(q)), (1.0, 0.0));
        assert_eq!(sincos(Phase::from_bits(2 * q)), (0.0, -1.0));
        assert_eq!(sincos(Phase::from_bits(3 * q)), (-1.0, 0.0));
    }
}
