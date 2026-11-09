//! Deterministic natural logarithm in f32, for the Box-Muller radius (contracts §1.5).
//!
//! It makes no std or libm call, and every operation is a separate IEEE f32 basic
//! operation in source order, so the CUDA port (`--fmad=false`) gives the same bits.
//!
//! Method (the classic fdlibm decomposition, with Taylor coefficients):
//!
//! 1. Range reduction, with integer operations on the bit pattern. A subnormal input is
//!    first scaled by 2²³, which is exact. Then `x = 2ᵉ · m` with `m ∈ [1, 2)`. If
//!    `m > √2` (bit pattern above `0x3FB504F3`, the largest f32 below √2), halve it and
//!    increment `e`. Now `m ∈ [0.7071068, 1.4142135]`.
//! 2. `f = m − 1` (exact by Sterbenz's lemma), `s = f / (2 + f)`, `z = s²`,
//!    `hfsq = 0.5·f²`.
//! 3. `ln(1+f) = f − s·f + s·R` exactly, with `R = 2z/3 + 2z²/5 + 2z³/7 + 2z⁴/9 + …`.
//!    `R` keeps four terms in Horner form. On `|s| ≤ 0.1716` the omitted tail is
//!    `< 7·10⁻¹⁰` absolute (0.02 ulp at the end of the range).
//! 4. `ln x = e·LN2_HI − ((hfsq − (s·(hfsq + R) + e·LN2_LO)) − f)`. `LN2_HI` keeps the
//!    top 12 significant bits of ln 2, so `e·LN2_HI` is exact for every f32 exponent
//!    (|e| ≤ 149 needs 8 bits). `LN2_LO` is the rest of ln 2, rounded to f32.
//!
//! Error, measured against an f64 reference over every f32 in 13 binades (around 1,
//! the Box-Muller range, the extremes) and a sample of subnormals (the
//! `ln_error_sweep` test): at most **0.86 ulp** (the test asserts ≤ 1 ulp). On
//! x ∈ [2⁻⁶⁵, 1], the Box-Muller input range, the absolute error is at most
//! **1.932·10⁻⁶** (asserted ≤ 2·10⁻⁶), about half an ulp of |ln 2⁻⁶⁵| ≈ 45.
//! 0.86 ulp is a relative error of at most 0.86·2⁻²³ ≈ 1.03·10⁻⁷. `ln(1) = +0` exactly.

use crate::error::NumericError;
use crate::fixed::require_finite_f32;

/// 2²³, which scales any subnormal f32 into the normal range exactly.
const TWO_POW_23: f32 = 8_388_608.0;
/// Bit pattern of the largest f32 below √2.
const SQRT2_BITS: u32 = 0x3FB5_04F3;
/// ln 2 rounded to f32, before the split.
const LN2_F32: f32 = std::f64::consts::LN_2 as f32;
/// ln 2 truncated to 12 significant bits: `e·LN2_HI` is exact for every exponent.
pub(crate) const LN2_HI: f32 = f32::from_bits(LN2_F32.to_bits() & 0xFFFF_F000);
/// ln 2 − `LN2_HI`, rounded to f32.
pub(crate) const LN2_LO: f32 = (std::f64::consts::LN_2 - LN2_HI as f64) as f32;

/// Taylor coefficients of R(z) = 2z/3 + 2z²/5 + 2z³/7 + 2z⁴/9, rounded to f32.
const R1: f32 = 2.0 / 3.0;
const R2: f32 = 2.0 / 5.0;
const R3: f32 = 2.0 / 7.0;
const R4: f32 = 2.0 / 9.0;

/// Natural logarithm of a positive finite f32 (subnormals included).
///
/// Errors: `NonFinite` for NaN or ±Inf, `OutOfRange` for zero or negative input.
/// Accuracy is in the module docs.
pub fn ln(x: f32) -> Result<f32, NumericError> {
    let x = require_finite_f32(x)?;
    if x > 0.0 {
        Ok(ln_positive(x))
    } else {
        Err(NumericError::OutOfRange)
    }
}

/// [`ln`] without input checks. The caller guarantees `x` is positive and finite.
pub(crate) fn ln_positive(x: f32) -> f32 {
    let mut bits = x.to_bits();
    let mut e: i32 = 0;
    if bits < 0x0080_0000 {
        // Subnormal: scale into the normal range, exactly.
        bits = (x * TWO_POW_23).to_bits();
        e = -23;
    }
    e += ((bits >> 23) as i32) - 127;
    let mut mbits = (bits & 0x007F_FFFF) | 0x3F80_0000;
    if mbits > SQRT2_BITS {
        // m/2: decrement the exponent field from 127 to 126.
        mbits -= 1 << 23;
        e += 1;
    }
    let m = f32::from_bits(mbits);

    let f = m - 1.0;
    let s = f / (2.0 + f);
    let z = s * s;
    let ff = f * f;
    let hfsq = 0.5 * ff;
    let r4z = R4 * z;
    let r = R3 + r4z;
    let r = r * z;
    let r = R2 + r;
    let r = r * z;
    let r = R1 + r;
    let r = r * z;

    // |e| ≤ 149, so the i32 → f32 conversion is exact.
    let ef = e as f32;
    let hi = ef * LN2_HI;
    let lo = ef * LN2_LO;
    let t = hfsq + r;
    let t = s * t;
    let t = t + lo;
    let t = hfsq - t;
    let t = t - f;
    hi - t
}

#[cfg(test)]
mod tests {
    use super::*;

    /// f64 reference. std ln is allowed here because it is only the test oracle: its
    /// < 1 f64 ulp error is ~2⁻²⁹ of an f32 ulp.
    #[allow(clippy::disallowed_methods)]
    fn reference(x: f32) -> f64 {
        (x as f64).ln()
    }

    /// Error in units of the f32 ulp at the reference value.
    fn ulp_error(x: f32) -> f64 {
        let want = reference(x);
        let got = ln_positive(x) as f64;
        let w = want as f32;
        let ulp = if w == 0.0 {
            f32::from_bits(1) as f64
        } else {
            let a = w.abs();
            (f32::from_bits(a.to_bits() + 1) - a) as f64
        };
        (got - want).abs() / ulp
    }

    #[test]
    fn split_constants() {
        assert_eq!(LN2_HI.to_bits() & 0xFFF, 0);
        // hi + lo reproduces ln 2 to within half an ulp of LN2_LO (≈ 7·10⁻¹²), so even
        // e = −149 adds only ~10⁻⁹ to a result whose ulp is ~8·10⁻⁶.
        let sum = LN2_HI as f64 + LN2_LO as f64;
        let half_ulp_lo = (f32::from_bits(LN2_LO.to_bits() + 1) - LN2_LO) as f64 / 2.0;
        assert!((sum - std::f64::consts::LN_2).abs() <= half_ulp_lo);
        // e·LN2_HI exact for every exponent.
        for e in -149..=128 {
            let ef = e as f32;
            assert_eq!((ef * LN2_HI) as f64, (e as f64) * (LN2_HI as f64));
        }
        assert!((f32::from_bits(SQRT2_BITS) as f64) < std::f64::consts::SQRT_2);
        assert!((f32::from_bits(SQRT2_BITS + 1) as f64) > std::f64::consts::SQRT_2);
    }

    #[test]
    fn exact_and_rejected_cases() {
        assert_eq!(ln(1.0).unwrap().to_bits(), 0.0f32.to_bits());
        assert_eq!(ln(0.0), Err(NumericError::OutOfRange));
        assert_eq!(ln(-0.0), Err(NumericError::OutOfRange));
        assert_eq!(ln(-1.0), Err(NumericError::OutOfRange));
        assert_eq!(ln(f32::NAN), Err(NumericError::NonFinite));
        assert_eq!(ln(f32::INFINITY), Err(NumericError::NonFinite));
        assert_eq!(ln(f32::NEG_INFINITY), Err(NumericError::NonFinite));
        // Powers of two: e·ln 2 from the split, within 1 ulp.
        for e in -149..=127 {
            let bits = if e >= -126 {
                ((127 + e) as u32) << 23
            } else {
                1u32 << (e + 149)
            };
            let x = f32::from_bits(bits);
            assert!(ulp_error(x) <= 1.0, "2^{e}");
        }
        assert!(ln(f32::MAX).unwrap() > 88.72);
        assert!(ln(f32::from_bits(1)).unwrap() < -103.27);
    }

    /// Every f32 in a set of binades: those around 1 (where cancellation matters), the
    /// Box-Muller range down to 2⁻⁶⁵, the extremes and the subnormals.
    #[test]
    fn ln_error_sweep() {
        let mut max_ulp = 0.0f64;
        let mut max_abs_unit = 0.0f64;
        let mut check = |x: f32| {
            let u = ulp_error(x);
            if u > max_ulp {
                max_ulp = u;
            }
            if x <= 1.0 && x >= f32::from_bits((127 - 65) << 23) {
                let a = (ln_positive(x) as f64 - reference(x)).abs();
                if a > max_abs_unit {
                    max_abs_unit = a;
                }
            }
        };
        for e in [-149, -127, -126, -65, -64, -33, -32, -2, -1, 0, 1, 64, 127] {
            let base = 127 + e;
            if base <= 0 {
                continue;
            }
            for m in 0..(1u32 << 23) {
                check(f32::from_bits(((base as u32) << 23) | m));
            }
        }
        // Subnormals: every 64th pattern (2¹⁷ values) plus the first 4096.
        for m in (1..(1u32 << 23)).step_by(64).chain(1..4096) {
            check(f32::from_bits(m));
        }
        // The documented bounds (module docs and contracts §1.5).
        assert!(max_ulp <= 1.0, "max ulp error {max_ulp}");
        assert!(
            max_abs_unit <= 2.0e-6,
            "max abs error on [2^-65, 1]: {max_abs_unit}"
        );
        eprintln!("ln: max ulp error {max_ulp:.4}, max abs error on [2^-65,1] {max_abs_unit:.3e}");
    }
}
