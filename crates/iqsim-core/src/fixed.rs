//! Fixed-point delay and phase outputs, and the float → integer conversion rules.
//!
//! - [`Delay`]: path delay in output samples, Q32.32 in an i64 (contracts §1.1).
//!   Range-checked; it never wraps.
//! - [`Phase`]: carrier phase as a u64 fraction of a cycle (2⁶⁴ = one cycle). It wraps
//!   by definition.
//!
//! Float → integer conversion follows contracts §1.4: reject NaN and ±Inf, round to
//! nearest with ties to even (`round_ties_even`), then convert with an explicit range
//! check. Rust's `as` is used only after that check, where it is exact. Scaling by a
//! power of two 2ⁿ (n ≥ 0) is an exact IEEE multiply: scaling up never discards
//! significand bits, even when the input or the product is subnormal (for example
//! 2⁻¹⁰⁷⁴ · 2³² = 2⁻¹⁰⁴²). A product that overflows to ±Inf is rejected, and a
//! subnormal product rounds to zero.

use crate::NumericError;
use crate::wide::U128;

/// Fractional bits of the [`Delay`] output (Q32.32).
pub const DELAY_FRAC_BITS: u32 = 32;

/// Extra fractional bits a block-polynomial register carries beyond its output.
pub const POLY_EXTRA_BITS: u32 = 64;

/// The f64 value 2ⁿ, built from its bit pattern. Exact for n ≤ 1023.
const fn pow2(n: u32) -> f64 {
    f64::from_bits(((1023 + n) as u64) << 52)
}

/// `x` unchanged if it is finite, else [`NumericError::NonFinite`].
pub fn require_finite_f64(x: f64) -> Result<f64, NumericError> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(NumericError::NonFinite)
    }
}

/// `x` unchanged if it is finite, else [`NumericError::NonFinite`].
pub fn require_finite_f32(x: f32) -> Result<f32, NumericError> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(NumericError::NonFinite)
    }
}

/// Round `x` to the nearest integer, ties to even, and convert to i64.
///
/// [`NumericError::NonFinite`] for NaN/±Inf and [`NumericError::OutOfRange`] if the
/// rounded value is outside `[−2⁶³, 2⁶³)`. This is the CPU equivalent of CUDA
/// `__double2ll_rn` plus the range check.
pub fn f64_to_i64_ties_even(x: f64) -> Result<i64, NumericError> {
    let r = require_finite_f64(x)?.round_ties_even();
    if r >= -pow2(63) && r < pow2(63) {
        // r is an integer in range, so the cast is exact.
        Ok(r as i64)
    } else {
        Err(NumericError::OutOfRange)
    }
}

/// [`f64_to_i64_ties_even`] for an f32 input. The f32 → f64 widening is exact, so the
/// result equals rounding the f32 directly (CUDA `__float2ll_rn`).
pub fn f32_to_i64_ties_even(x: f32) -> Result<i64, NumericError> {
    f64_to_i64_ties_even(f64::from(x))
}

/// `round_ties_even(x · 2^frac_bits)` as a signed 128-bit register, checked to lie in
/// `[−2¹²⁷, 2¹²⁷)`.
pub(crate) fn f64_to_fixed_i128(x: f64, frac_bits: u32) -> Result<U128, NumericError> {
    debug_assert!(frac_bits <= 128);
    let y = require_finite_f64(x)? * pow2(frac_bits);
    let r = require_finite_f64(y)?.round_ties_even();
    if r >= -pow2(127) && r < pow2(127) {
        // r is an integer in range, so the cast is exact.
        Ok(U128::from_i128(r as i128))
    } else {
        Err(NumericError::OutOfRange)
    }
}

/// `round_ties_even(x · 2^frac_bits) mod 2^frac_bits`, for a quantity measured in
/// whole cycles that wraps (phase). `frac_bits` is 64 or 128; the result occupies the
/// low `frac_bits` bits.
///
/// The integer cycles are removed first with `x − trunc(x)`, which is exact in IEEE
/// arithmetic, so no precision is lost for large `|x|`. Rounding commutes with adding
/// whole cycles, so this equals rounding `x · 2^frac_bits` and then reducing.
pub(crate) fn f64_cycles_to_fixed_wrapping(x: f64, frac_bits: u32) -> Result<U128, NumericError> {
    debug_assert!(frac_bits == 64 || frac_bits == 128);
    let x = require_finite_f64(x)?;
    let f = x - x.trunc(); // exact, in (−1, 1)
    let r = (f * pow2(frac_bits)).round_ties_even(); // |r| ≤ 2^frac_bits
    let mag = r.abs();
    let one_cycle = pow2(frac_bits);
    let m = if mag >= one_cycle {
        // |r| = 2^frac_bits: a whole cycle, which is 0 mod one cycle.
        U128::ZERO
    } else {
        // mag is an integer < 2^frac_bits ≤ 2¹²⁸, so the cast is exact.
        U128::from_u128(mag as u128)
    };
    let v = if r < 0.0 { m.wrapping_neg() } else { m };
    Ok(if frac_bits == 64 {
        U128::new(0, v.lo)
    } else {
        v
    })
}

/// Path delay in output samples, Q32.32 fixed point in an i64.
///
/// Range ±2³¹ samples, resolution 2⁻³² sample. It is range-checked and never wraps.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct Delay(i64);

impl Delay {
    /// The delay whose Q32.32 bit pattern is `bits`.
    pub const fn from_bits(bits: i64) -> Self {
        Self(bits)
    }

    /// The Q32.32 bit pattern.
    pub const fn to_bits(self) -> i64 {
        self.0
    }

    /// Nearest Q32.32 delay to `samples`, ties to even. Errors for non-finite input or
    /// a value outside `[−2³¹, 2³¹)` samples after rounding.
    pub fn from_samples_f64(samples: f64) -> Result<Self, NumericError> {
        let y = require_finite_f64(samples)? * pow2(DELAY_FRAC_BITS);
        f64_to_i64_ties_even(y).map(Self)
    }

    /// The delay in samples as f64. Diagnostics only: it rounds when the value needs
    /// more than 53 significant bits, and the render path never uses it.
    pub fn to_samples_f64(self) -> f64 {
        // i64 → f64 rounds to nearest (an integer → float conversion, not float → int).
        (self.0 as f64) / pow2(DELAY_FRAC_BITS)
    }
}

/// Carrier phase as a fraction of one cycle: the value is `bits / 2⁶⁴` cycles.
///
/// Arithmetic wraps mod one cycle by definition (contracts §1.1).
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct Phase(u64);

impl Phase {
    /// Zero phase.
    pub const ZERO: Self = Self(0);

    /// The phase `bits / 2⁶⁴` cycles.
    pub const fn from_bits(bits: u64) -> Self {
        Self(bits)
    }

    /// The u64 cycle fraction.
    pub const fn to_bits(self) -> u64 {
        self.0
    }

    /// `self + rhs` mod one cycle.
    #[must_use]
    pub const fn wrapping_add(self, rhs: Self) -> Self {
        Self(self.0.wrapping_add(rhs.0))
    }

    /// `self − rhs` mod one cycle.
    #[must_use]
    pub const fn wrapping_sub(self, rhs: Self) -> Self {
        Self(self.0.wrapping_sub(rhs.0))
    }

    /// Nearest phase to `cycles` (any finite value), ties to even, reduced mod one
    /// cycle without losing precision for large `|cycles|`.
    pub fn from_cycles_f64(cycles: f64) -> Result<Self, NumericError> {
        f64_cycles_to_fixed_wrapping(cycles, 64).map(|v| Self(v.lo))
    }

    /// The phase in cycles, in `[0, 1)`, as f64. Diagnostics only (rounds to 53 bits).
    pub fn to_cycles_f64(self) -> f64 {
        // u64 → f64 rounds to nearest; a value just below 2⁶⁴ rounds to 1.0, so clamp
        // to keep the documented half-open range.
        let c = (self.0 as f64) / pow2(64);
        if c >= 1.0 { 0.0 } else { c }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ties_to_even() {
        let cases = [
            (0.5, 0),
            (1.5, 2),
            (2.5, 2),
            (-0.5, 0),
            (-1.5, -2),
            (-2.5, -2),
            (0.499_999_999_999_999_94, 0),
            (1.0 - f64::EPSILON / 2.0, 1),
        ];
        for (x, want) in cases {
            assert_eq!(f64_to_i64_ties_even(x), Ok(want), "{x}");
        }
        assert_eq!(f32_to_i64_ties_even(2.5), Ok(2));
        assert_eq!(f32_to_i64_ties_even(3.5), Ok(4));
    }

    #[test]
    fn i64_range_boundaries() {
        assert_eq!(f64_to_i64_ties_even(-pow2(63)), Ok(i64::MIN));
        assert_eq!(
            f64_to_i64_ties_even(pow2(63)),
            Err(NumericError::OutOfRange)
        );
        // Largest f64 below 2⁶³ is 2⁶³ − 1024.
        let below = f64::from_bits(pow2(63).to_bits() - 1);
        assert_eq!(f64_to_i64_ties_even(below), Ok(i64::MAX - 1023));
        assert_eq!(f64_to_i64_ties_even(-1e300), Err(NumericError::OutOfRange));
    }

    #[test]
    fn non_finite_rejected() {
        for x in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert_eq!(f64_to_i64_ties_even(x), Err(NumericError::NonFinite));
            assert_eq!(Delay::from_samples_f64(x), Err(NumericError::NonFinite));
            assert_eq!(Phase::from_cycles_f64(x), Err(NumericError::NonFinite));
            assert_eq!(f64_to_fixed_i128(x, 96), Err(NumericError::NonFinite));
        }
        for x in [f32::NAN, f32::INFINITY, f32::NEG_INFINITY] {
            assert_eq!(f32_to_i64_ties_even(x), Err(NumericError::NonFinite));
            assert_eq!(require_finite_f32(x), Err(NumericError::NonFinite));
        }
        // Finite input whose scaled value overflows to Inf.
        assert_eq!(
            f64_to_fixed_i128(f64::MAX, 96),
            Err(NumericError::NonFinite)
        );
    }

    #[test]
    fn delay_rounding_and_range() {
        let lsb = pow2(DELAY_FRAC_BITS).recip();
        assert_eq!(Delay::from_samples_f64(0.5 * lsb).unwrap().to_bits(), 0);
        assert_eq!(Delay::from_samples_f64(1.5 * lsb).unwrap().to_bits(), 2);
        assert_eq!(Delay::from_samples_f64(-1.5 * lsb).unwrap().to_bits(), -2);
        assert_eq!(Delay::from_samples_f64(-3.25).unwrap().to_bits(), -13 << 30);
        assert_eq!(
            Delay::from_samples_f64(-pow2(31)).unwrap().to_bits(),
            i64::MIN
        );
        assert_eq!(
            Delay::from_samples_f64(pow2(31)),
            Err(NumericError::OutOfRange)
        );
        assert_eq!(Delay::from_bits(-13 << 30).to_samples_f64(), -3.25);
    }

    #[test]
    fn phase_wraps_from_cycles() {
        assert_eq!(Phase::from_cycles_f64(1.25).unwrap().to_bits(), 1 << 62);
        assert_eq!(Phase::from_cycles_f64(-0.25).unwrap().to_bits(), 3 << 62);
        assert_eq!(Phase::from_cycles_f64(-1.0).unwrap().to_bits(), 0);
        assert_eq!(Phase::from_cycles_f64(1e30).unwrap().to_bits(), 0);
        // 2⁴⁰ + 0.5 cycles: the half cycle survives the large integer part.
        assert_eq!(
            Phase::from_cycles_f64(pow2(40) + 0.5).unwrap().to_bits(),
            1 << 63
        );
        // A tie at the 2⁻⁶⁴ grid rounds to even.
        let half_lsb = pow2(65).recip();
        assert_eq!(Phase::from_cycles_f64(half_lsb).unwrap().to_bits(), 0);
        assert_eq!(Phase::from_cycles_f64(3.0 * half_lsb).unwrap().to_bits(), 2);
        // −½ LSB rounds to 0 (even), not to 2⁶⁴ − 1.
        assert_eq!(Phase::from_cycles_f64(-half_lsb).unwrap().to_bits(), 0);
        // −1.5 LSB rounds to −2 ≡ 2⁶⁴ − 2.
        assert_eq!(
            Phase::from_cycles_f64(-3.0 * half_lsb).unwrap().to_bits(),
            u64::MAX - 1
        );
        assert_eq!(Phase::from_bits(u64::MAX).to_cycles_f64(), 0.0);
        assert_eq!(Phase::from_bits(1 << 63).to_cycles_f64(), 0.5);
    }

    #[test]
    fn wrapping_cycles_128() {
        let v = f64_cycles_to_fixed_wrapping(-0.25, 128).unwrap();
        assert_eq!(v, U128::new(3 << 62, 0));
        // 1 − 2⁻⁵³ cycles is exact at 128 bits: 2¹²⁸ − 2⁷⁵.
        let v = f64_cycles_to_fixed_wrapping(1.0 - f64::EPSILON / 2.0, 128).unwrap();
        assert_eq!(v, U128::new(u64::MAX << 11, 0));
        // −2⁻¹³⁰ cycles is −¼ LSB: rounds to zero.
        let v = f64_cycles_to_fixed_wrapping(-pow2(130).recip(), 128).unwrap();
        assert_eq!(v, U128::ZERO);
    }

    #[test]
    fn fixed_i128_range() {
        assert_eq!(
            f64_to_fixed_i128(-pow2(31), 96).unwrap(),
            U128::new(1 << 63, 0)
        );
        assert_eq!(
            f64_to_fixed_i128(pow2(31), 96),
            Err(NumericError::OutOfRange)
        );
        assert_eq!(f64_to_fixed_i128(-1.0, 96).unwrap().to_i128(), -(1 << 96));
    }
}
