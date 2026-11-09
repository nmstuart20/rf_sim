//! Carrier-offset phase from exact time (contracts §1.6).
//!
//! φ(t) = frac(Δf · t), computed in integer arithmetic from the exact Q64.64 ps
//! acquisition instant. Nothing forms a large floating-point Δf·t and drops its integer
//! cycles, so precision near 100 days equals precision near t = 0.
//!
//! **Representation.** Δf is an i64 in units of 2⁻³² Hz ([`OffsetFrequency`]): range
//! ±2³¹ Hz ≈ ±2.1 GHz, resolution 2.3·10⁻¹⁰ Hz. With t = T·2⁻⁶⁴ ps (T an i128), the phase
//! in units of 2⁻⁶⁴ cycle is
//!
//! ```text
//! φ = ⌊ F·T / D ⌋ mod 2⁶⁴,   D = 2³²·10¹² = 2⁴⁴·5¹²
//! ```
//!
//! floored toward −∞, so the error against the exact phase is in [0, 1) unit
//! (< 2⁻⁶⁴ cycle) for every t.
//!
//! **Algorithm** (64-bit primitives only, so a CUDA port can mirror it):
//!
//! 1. A = |F|·|T|, a 192-bit product from two `u64 × u64 → u128` products.
//! 2. Divide by 2⁴⁴ with a three-word shift, remembering whether the dropped bits were 0.
//! 3. Divide by 5¹² = 244 140 625 with schoolbook division over 32-bit digits: the
//!    running remainder is < 5¹² < 2²⁸, so each partial dividend is < 2⁶⁰ and fits a u64.
//!    Keep the low 64 quotient bits and the remainder.
//! 4. For a negative product, ⌊−A/D⌋ = −⌈A/D⌉: add 1 if the division was inexact, then
//!    negate mod 2⁶⁴.

use crate::fixed::{Phase, f64_to_i64_ties_even, require_finite_f64};
use crate::wide::{U128, add_with_carry};
use crate::{AcqInstant, NumericError};

/// 5¹², the odd factor of D = 2³²·10¹².
const FIVE_POW_12: u64 = 244_140_625;

/// A carrier frequency offset Δf in units of 2⁻³² Hz.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct OffsetFrequency(i64);

impl OffsetFrequency {
    /// Fractional bits: one unit is 2⁻³² Hz.
    pub const FRAC_BITS: u32 = 32;

    /// The offset `bits · 2⁻³²` Hz.
    pub const fn from_bits(bits: i64) -> Self {
        Self(bits)
    }

    /// The raw value in units of 2⁻³² Hz.
    pub const fn to_bits(self) -> i64 {
        self.0
    }

    /// An exact whole-hertz offset. [`NumericError::OutOfRange`] unless
    /// `−2³¹ ≤ hz < 2³¹`.
    pub fn from_hz(hz: i64) -> Result<Self, NumericError> {
        hz.checked_mul(1 << Self::FRAC_BITS)
            .map(Self)
            .ok_or(NumericError::OutOfRange)
    }

    /// The nearest representable offset to `hz`, ties to even. Errors for a non-finite
    /// value or one outside `[−2³¹, 2³¹)` Hz after rounding.
    pub fn from_hz_f64(hz: f64) -> Result<Self, NumericError> {
        // Multiplying by 2³² is exact (see `fixed`).
        let scaled = require_finite_f64(hz)? * 4_294_967_296.0;
        f64_to_i64_ties_even(scaled).map(Self)
    }

    /// Offset phase φ = frac(Δf·t) at the acquisition instant `t`, floored to the
    /// 2⁻⁶⁴ cycle grid. Exact for every i64 frequency and every instant.
    pub fn phase_at(self, t: AcqInstant) -> Phase {
        let tq = t.to_q64_64();
        let negative = (self.0 < 0) != (tq < 0);
        let f = self.0.unsigned_abs();
        let tm = U128::from_u128(tq.unsigned_abs());

        // 1. A = f · |T| as three 64-bit words (l2, l1, l0). A < 2¹⁹⁰.
        let p0 = U128::mul_wide_u64(f, tm.lo);
        let p1 = U128::mul_wide_u64(f, tm.hi);
        let l0 = p0.lo;
        let (l1, carry) = add_with_carry(p0.hi, p1.lo, 0);
        let l2 = p1.hi + carry; // no overflow: A < 2¹⁹⁰

        // 2. B = A >> 44.
        let low_bits_zero = (l0 & ((1 << 44) - 1)) == 0;
        let b0 = (l0 >> 44) | (l1 << 20);
        let b1 = (l1 >> 44) | (l2 << 20);
        let b2 = l2 >> 44;

        // 3. B / 5¹² over 32-bit digits, most significant first.
        let digits = [
            b2 >> 32,
            b2 & 0xFFFF_FFFF,
            b1 >> 32,
            b1 & 0xFFFF_FFFF,
            b0 >> 32,
            b0 & 0xFFFF_FFFF,
        ];
        let mut rem: u64 = 0;
        let mut q_hi: u64 = 0; // second-lowest quotient digit
        let mut q_lo: u64 = 0; // lowest quotient digit
        for digit in digits {
            let cur = (rem << 32) | digit; // < 5¹²·2³² < 2⁶⁰
            q_hi = q_lo;
            q_lo = cur / FIVE_POW_12; // < 2³²
            rem = cur % FIVE_POW_12;
        }
        let q = (q_hi << 32) | q_lo;

        // 4. Floor for a negative product.
        let bits = if negative {
            let inexact = !low_bits_zero || rem != 0;
            q.wrapping_add(inexact as u64).wrapping_neg()
        } else {
            q
        };
        Phase::from_bits(bits)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{PS_PER_SECOND, TruthTime};

    fn at_ps(ps: i64) -> AcqInstant {
        AcqInstant::from_truth_time(TruthTime::from_ps(ps))
    }

    const DAY_PS: i64 = 86_400 * PS_PER_SECOND;

    #[test]
    fn closed_form_values() {
        let one_hz = OffsetFrequency::from_hz(1).unwrap();
        let khz = OffsetFrequency::from_hz(1000).unwrap();
        let q = PS_PER_SECOND / 4;
        assert_eq!(one_hz.phase_at(at_ps(q)).to_bits(), 1 << 62);
        assert_eq!(one_hz.phase_at(at_ps(100 * DAY_PS + q)).to_bits(), 1 << 62);
        assert_eq!(one_hz.phase_at(at_ps(-q)).to_bits(), 3 << 62);
        let minus = OffsetFrequency::from_hz(-1).unwrap();
        assert_eq!(minus.phase_at(at_ps(q)).to_bits(), 3 << 62);
        assert_eq!(minus.phase_at(at_ps(-q)).to_bits(), 1 << 62);
        // 1 kHz at 100 days + 0.5 ms is exactly half a cycle.
        let half_ms = PS_PER_SECOND / 2000;
        assert_eq!(
            khz.phase_at(at_ps(100 * DAY_PS + half_ms)).to_bits(),
            1 << 63
        );
        assert_eq!(khz.phase_at(at_ps(0)).to_bits(), 0);
    }

    #[test]
    fn floor_for_tiny_negative_phase() {
        // 1 Hz for one 2⁻⁶⁴ ps unit before the epoch: a tiny negative phase floors to
        // 2⁶⁴ − 1, not to 0.
        let one_hz = OffsetFrequency::from_hz(1).unwrap();
        let t = AcqInstant::from_q64_64(-1);
        assert_eq!(one_hz.phase_at(t).to_bits(), u64::MAX);
        assert_eq!(one_hz.phase_at(AcqInstant::from_q64_64(1)).to_bits(), 0);
    }

    #[test]
    fn extremes_do_not_overflow() {
        for f in [i64::MIN, i64::MAX, -1, 1] {
            for t in [i128::MIN, i128::MAX, -1, 0, 1] {
                let _ = OffsetFrequency::from_bits(f).phase_at(AcqInstant::from_q64_64(t));
            }
        }
        // Δf = −2³¹ Hz over 1 ps is exactly −2³¹·10⁻¹² cycles.
        let f = OffsetFrequency::from_bits(i64::MIN);
        let got = f.phase_at(at_ps(1)).to_bits();
        // ⌊−2⁶³·2⁶⁴ / (2⁴⁴·5¹²)⌋ = −⌈2⁸³/5¹²⌉.
        let want = (1u128 << 83).div_ceil(u128::from(FIVE_POW_12)) as u64;
        assert_eq!(got, want.wrapping_neg());
    }

    #[test]
    fn frequency_constructors() {
        assert_eq!(OffsetFrequency::from_hz(1).unwrap().to_bits(), 1 << 32);
        assert_eq!(
            OffsetFrequency::from_hz(1 << 31),
            Err(NumericError::OutOfRange)
        );
        assert_eq!(
            OffsetFrequency::from_hz(-(1 << 31)).unwrap().to_bits(),
            i64::MIN
        );
        assert_eq!(
            OffsetFrequency::from_hz_f64(0.5).unwrap().to_bits(),
            1 << 31
        );
        // Ties to even at half a unit.
        let half_unit = 1.0 / 8_589_934_592.0;
        assert_eq!(
            OffsetFrequency::from_hz_f64(half_unit).unwrap().to_bits(),
            0
        );
        assert_eq!(
            OffsetFrequency::from_hz_f64(3.0 * half_unit)
                .unwrap()
                .to_bits(),
            2
        );
        assert_eq!(
            OffsetFrequency::from_hz_f64(f64::NAN),
            Err(NumericError::NonFinite)
        );
        assert_eq!(
            OffsetFrequency::from_hz_f64(3e9),
            Err(NumericError::OutOfRange)
        );
    }
}
