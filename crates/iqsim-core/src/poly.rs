//! Block polynomials for delay and carrier phase in Newton forward-difference form
//! (contracts §1.3).
//!
//! A [`BlockPoly`] holds four 128-bit registers D0..D3 for a block of `len` samples.
//! The value at local sample index k is
//!
//! ```text
//! p(k) = D0 + C(k,1)·D1 + C(k,2)·D2 + C(k,3)·D3      (mod 2¹²⁸)
//! ```
//!
//! A thread can jump straight to any k with that formula ([`BlockPoly::jump`]), then
//! advance with three 128-bit adds per sample ([`PolyCursor::step`]). Because both use
//! the same wrapping integer arithmetic and C(k,j) is exact, they agree bit for bit at
//! every k. The registers carry [`POLY_EXTRA_BITS`] fractional bits beyond the output:
//!
//! - [`DelayPoly`]: registers are Q32.96 samples; output is [`Delay`] (Q32.32).
//! - [`PhasePoly`]: registers are 0.128 cycle fractions; output is [`Phase`] (0.64).
//!
//! Output extraction keeps the high 64 bits (floor toward −∞). It does not round.
//!
//! **Maximum block length.** The jump needs C(k,3) as an exact u64, so
//! k ≤ 4 801 280, the largest k with C(k,3) < 2⁶⁴, and [`MAX_BLOCK_LEN`] = 4 801 281.
//! That is 4.8 ms at 1 GS/s and 19 ms at 250 MS/s.
//!
//! **Coefficient quantization.** If each Dj is the nearest register value to the
//! exact difference, |error_j| ≤ ½ register unit (2⁻⁶⁵ output LSB). The value error at
//! k is at most ½·(1 + C(k,1) + C(k,2) + C(k,3)) register units, which for every
//! k < [`MAX_BLOCK_LEN`] is below 2⁶³ + 2⁴³ units, i.e. ½ + 2⁻²¹ output LSB. With the
//! floor at extraction, the output is within 1.5 + 2⁻²¹ LSB (so under 2 LSB) of the
//! exact polynomial: < 4.7·10⁻¹⁰ samples of delay (budget 10⁻³) and < 6.9·10⁻¹⁹ rad of
//! phase (budget 10⁻²). For a 1 ms block at 500 MS/s, C(k,3) < 2⁵⁵ and the register
//! error is below 2⁻⁹ LSB.
//!
//! **Delay range margin.** The registers wrap, so evaluation cannot detect a delay that
//! leaves the Q32.32 range. Quantization can move the value by up to 2 LSB, so a true
//! delay that stays inside the range but within 2 LSB of −2³¹ can still wrap (to near
//! `i64::MAX`). The geometry engine (phase 3) must therefore keep the true delay over
//! the whole block within [−2³¹ + m, 2³¹ − m) samples, m = [`DELAY_RANGE_MARGIN_LSB`]
//! LSB. A delay polynomial outside that margin is a fit error.

use std::marker::PhantomData;

use crate::NumericError;
use crate::fixed::{
    DELAY_FRAC_BITS, Delay, POLY_EXTRA_BITS, Phase, f64_cycles_to_fixed_wrapping, f64_to_fixed_i128,
};
use crate::wide::U128;

/// Margin, in Q32.32 LSB, the true delay must keep from both ends of the Q32.32 range
/// so that coefficient quantization (< 2 LSB) cannot wrap the output.
pub const DELAY_RANGE_MARGIN_LSB: i64 = 2;

/// Largest supported block-polynomial length. Local indices are `0..MAX_BLOCK_LEN`,
/// and C(4 801 280, 3) is the largest C(k,3) below 2⁶⁴.
pub const MAX_BLOCK_LEN: u64 = 4_801_281;

/// Exact binomial coefficients `(C(k,1), C(k,2), C(k,3))` in u64 arithmetic.
///
/// Requires `k < MAX_BLOCK_LEN`; for larger k, C(k,3) overflows (a panic with overflow
/// checks). No intermediate exceeds the final C(k,3): k(k−1) < 2⁴⁵, and C(k,3) is
/// formed as C(k,2)·(k−2)/3 with the division by 3 applied to whichever factor it
/// divides exactly (one of k, k−1, k−2 is a multiple of 3, and if it is not k−2, then 3
/// divides C(k,2) because gcd(2,3) = 1).
pub const fn binomials(k: u64) -> (u64, u64, u64) {
    let c1 = k;
    let c2 = if k < 2 { 0 } else { k * (k - 1) / 2 };
    let c3 = if k < 3 {
        0
    } else {
        let m = k - 2;
        if m.is_multiple_of(3) {
            c2 * (m / 3)
        } else {
            (c2 / 3) * m
        }
    };
    (c1, c2, c3)
}

/// An output format of a block polynomial: how coefficients are built from f64 and how
/// the register's high word becomes the output.
pub trait PolyOutput: Copy {
    /// The output for a register value (high 64 bits, floor).
    fn from_register(r: U128) -> Self;

    /// The register nearest to `x` output units (ties to even), or an error.
    fn coefficient_from_f64(x: f64) -> Result<U128, NumericError>;
}

impl PolyOutput for Delay {
    /// Arithmetic shift right by 64: the high word, read as signed.
    fn from_register(r: U128) -> Self {
        Delay::from_bits(r.hi as i64)
    }

    /// `x` samples (or samples per sampleʲ for Dj) as Q32.96, range-checked to
    /// `[−2³¹, 2³¹)`.
    fn coefficient_from_f64(x: f64) -> Result<U128, NumericError> {
        f64_to_fixed_i128(x, DELAY_FRAC_BITS + POLY_EXTRA_BITS)
    }
}

impl PolyOutput for Phase {
    /// The high word: a u64 cycle fraction.
    fn from_register(r: U128) -> Self {
        Phase::from_bits(r.hi)
    }

    /// `x` cycles (any finite value) as a 0.128 cycle fraction, reduced mod one cycle.
    /// Reducing Dj mod one cycle is exact because C(k,j) is an integer.
    fn coefficient_from_f64(x: f64) -> Result<U128, NumericError> {
        f64_cycles_to_fixed_wrapping(x, 128)
    }
}

/// A cubic (or lower) polynomial over one block, in Newton forward-difference form.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct BlockPoly<T: PolyOutput> {
    d: [U128; 4],
    len: u64,
    _output: PhantomData<T>,
}

/// Delay polynomial: Q32.96 registers, [`Delay`] output.
pub type DelayPoly = BlockPoly<Delay>;

/// Carrier-phase polynomial: 0.128 cycle registers, [`Phase`] output.
pub type PhasePoly = BlockPoly<Phase>;

impl<T: PolyOutput> BlockPoly<T> {
    /// The polynomial with registers `d = [D0, D1, D2, D3]` over `len` samples.
    /// [`NumericError::BlockLength`] unless `1 ≤ len ≤ MAX_BLOCK_LEN`.
    pub fn from_registers(d: [U128; 4], len: u64) -> Result<Self, NumericError> {
        if len == 0 || len > MAX_BLOCK_LEN {
            return Err(NumericError::BlockLength { len });
        }
        Ok(Self {
            d,
            len,
            _output: PhantomData,
        })
    }

    /// The polynomial whose forward differences at k = 0 are `d` (in output units:
    /// samples for delay, cycles for phase), each rounded to the nearest register.
    pub fn from_f64_differences(d: [f64; 4], len: u64) -> Result<Self, NumericError> {
        let r = [
            T::coefficient_from_f64(d[0])?,
            T::coefficient_from_f64(d[1])?,
            T::coefficient_from_f64(d[2])?,
            T::coefficient_from_f64(d[3])?,
        ];
        Self::from_registers(r, len)
    }

    /// The registers `[D0, D1, D2, D3]`.
    pub const fn registers(&self) -> [U128; 4] {
        self.d
    }

    /// Block length in samples.
    pub const fn block_len(&self) -> u64 {
        self.len
    }

    fn check_index(&self, k: u64) -> Result<(), NumericError> {
        if k < self.len {
            Ok(())
        } else {
            Err(NumericError::IndexOutsideBlock { k, len: self.len })
        }
    }

    /// Register value at local index `k`: D0 + C(k,1)·D1 + C(k,2)·D2 + C(k,3)·D3.
    pub fn jump(&self, k: u64) -> Result<U128, NumericError> {
        self.check_index(k)?;
        let (c1, c2, c3) = binomials(k);
        let [d0, d1, d2, d3] = self.d;
        Ok(d0
            .wrapping_add(d1.wrapping_mul_u64(c1))
            .wrapping_add(d2.wrapping_mul_u64(c2))
            .wrapping_add(d3.wrapping_mul_u64(c3)))
    }

    /// Output at local index `k` (high word of [`jump`](Self::jump)).
    pub fn eval(&self, k: u64) -> Result<T, NumericError> {
        self.jump(k).map(T::from_register)
    }

    /// A cursor positioned at local index `k`, for stepping through the block.
    ///
    /// Its difference registers are Δ1(k) = D1 + C(k,1)·D2 + C(k,2)·D3 and
    /// Δ2(k) = D2 + C(k,1)·D3.
    pub fn cursor(&self, k: u64) -> Result<PolyCursor<T>, NumericError> {
        let p = self.jump(k)?;
        let (c1, c2, _) = binomials(k);
        let [_, d1, d2, d3] = self.d;
        Ok(PolyCursor {
            p,
            d1: d1
                .wrapping_add(d2.wrapping_mul_u64(c1))
                .wrapping_add(d3.wrapping_mul_u64(c2)),
            d2: d2.wrapping_add(d3.wrapping_mul_u64(c1)),
            d3,
            _output: PhantomData,
        })
    }
}

/// Stepping state of a [`BlockPoly`]: the value and its first, second and third
/// forward differences at the current index.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PolyCursor<T: PolyOutput> {
    p: U128,
    d1: U128,
    d2: U128,
    d3: U128,
    _output: PhantomData<T>,
}

impl<T: PolyOutput> PolyCursor<T> {
    /// Register value at the current index.
    pub const fn register(&self) -> U128 {
        self.p
    }

    /// Output at the current index.
    pub fn value(&self) -> T {
        T::from_register(self.p)
    }

    /// Advance one sample: three wrapping 128-bit adds, in this order.
    ///
    /// The cursor does not track the block end. Stepping past it continues the same
    /// polynomial, which callers must not use as block output.
    pub fn step(&mut self) {
        self.p = self.p.wrapping_add(self.d1);
        self.d1 = self.d1.wrapping_add(self.d2);
        self.d2 = self.d2.wrapping_add(self.d3);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn binom_ref(k: u64, j: u32) -> u128 {
        let k = u128::from(k);
        match j {
            1 => k,
            2 => k * k.saturating_sub(1) / 2,
            3 => k * k.saturating_sub(1) * k.saturating_sub(2) / 6,
            _ => unreachable!(),
        }
    }

    #[test]
    fn binomials_exact_small_and_at_limit() {
        for k in (0..2000).chain(MAX_BLOCK_LEN - 2000..MAX_BLOCK_LEN) {
            let (c1, c2, c3) = binomials(k);
            assert_eq!(u128::from(c1), binom_ref(k, 1));
            assert_eq!(u128::from(c2), binom_ref(k, 2));
            assert_eq!(u128::from(c3), binom_ref(k, 3), "k = {k}");
        }
    }

    #[test]
    fn max_block_len_is_tight() {
        let last = MAX_BLOCK_LEN - 1;
        assert!(binom_ref(last, 3) <= u128::from(u64::MAX));
        assert!(binom_ref(MAX_BLOCK_LEN, 3) > u128::from(u64::MAX));
    }

    #[test]
    fn quantization_bound_below_one_output_lsb() {
        // ½·(1 + C(k,1) + C(k,2) + C(k,3)) register units, at the worst k.
        let k = MAX_BLOCK_LEN - 1;
        let sum = 1 + binom_ref(k, 1) + binom_ref(k, 2) + binom_ref(k, 3);
        let half = sum / 2 + 1;
        assert!(half < (1u128 << 63) + (1u128 << 43));
        assert!(
            half < 1u128 << 64,
            "must be below one output LSB (2^64 units)"
        );
        // 1 ms at 500 MS/s.
        let k = 499_999;
        let sum = 1 + binom_ref(k, 1) + binom_ref(k, 2) + binom_ref(k, 3);
        assert!(sum / 2 < 1u128 << 55);
    }

    #[test]
    fn length_and_index_checks() {
        let d = [U128::ZERO; 4];
        assert_eq!(
            DelayPoly::from_registers(d, 0),
            Err(NumericError::BlockLength { len: 0 })
        );
        assert_eq!(
            DelayPoly::from_registers(d, MAX_BLOCK_LEN + 1),
            Err(NumericError::BlockLength {
                len: MAX_BLOCK_LEN + 1
            })
        );
        let p = PhasePoly::from_registers(d, MAX_BLOCK_LEN).unwrap();
        assert!(p.eval(MAX_BLOCK_LEN - 1).is_ok());
        assert_eq!(
            p.eval(MAX_BLOCK_LEN),
            Err(NumericError::IndexOutsideBlock {
                k: MAX_BLOCK_LEN,
                len: MAX_BLOCK_LEN
            })
        );
    }

    #[test]
    fn delay_extraction_floors() {
        // D0 = −2⁻⁹⁶ samples (one register unit below zero) → floor gives −1 output LSB.
        let p =
            DelayPoly::from_registers([U128::from_i128(-1), U128::ZERO, U128::ZERO, U128::ZERO], 1)
                .unwrap();
        assert_eq!(p.eval(0).unwrap().to_bits(), -1);
    }
}
