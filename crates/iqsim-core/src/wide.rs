//! 128-bit wrapping integer arithmetic from explicit 64-bit primitives.
//!
//! The render path's block polynomials (contracts §1.3) are 128-bit registers that
//! wrap mod 2¹²⁸ by design. This module builds them from three 64-bit primitives,
//! [`mul_hi_u64`], wrapping 64-bit multiply and [`add_with_carry`], instead of native
//! `u128`/`i128`, so the CUDA C port (prompt 20) can mirror every line with
//! `__umul64hi`, plain `*` and a carry compare.
//!
//! Signed values are two's complement: the same bit patterns and the same wrapping
//! operations serve signed and unsigned interpretations. Conversions to and from
//! native `i128`/`u128` exist for host-side construction and tests only.
//!
//! Nothing here checks for overflow: every operation is wrapping on purpose. Use it
//! only for quantities whose definition is modular (polynomial registers, phase).

/// High 64 bits of the full 128-bit product `a · b`.
///
/// Built from four 32 × 32 → 64-bit partial products, so it needs no wider type.
/// Matches CUDA `__umul64hi`.
pub const fn mul_hi_u64(a: u64, b: u64) -> u64 {
    const MASK: u64 = 0xFFFF_FFFF;
    let a_lo = a & MASK;
    let a_hi = a >> 32;
    let b_lo = b & MASK;
    let b_hi = b >> 32;

    // Each partial product is < 2⁶⁴, so none of these multiplies can overflow.
    let ll = a_lo * b_lo;
    let lh = a_lo * b_hi;
    let hl = a_hi * b_lo;
    let hh = a_hi * b_hi;

    // Middle column: (ll >> 32) + low(lh) + low(hl) < 3 · 2³², no overflow.
    let mid = (ll >> 32) + (lh & MASK) + (hl & MASK);
    // hh + high(lh) + high(hl) + carry(mid) is the true high word, which is < 2⁶⁴.
    hh + (lh >> 32) + (hl >> 32) + (mid >> 32)
}

/// `a + b + carry_in` as `(sum mod 2⁶⁴, carry_out)`. `carry_in` and `carry_out` are 0 or 1.
pub const fn add_with_carry(a: u64, b: u64, carry_in: u64) -> (u64, u64) {
    let s = a.wrapping_add(b);
    let c1 = (s < a) as u64;
    let s2 = s.wrapping_add(carry_in);
    let c2 = (s2 < s) as u64;
    // At most one of c1 and c2 is set, since a + b + 1 ≤ 2¹²⁹ − 1.
    (s2, c1 | c2)
}

/// `a − b − borrow_in` as `(difference mod 2⁶⁴, borrow_out)`. Borrows are 0 or 1.
pub const fn sub_with_borrow(a: u64, b: u64, borrow_in: u64) -> (u64, u64) {
    let d = a.wrapping_sub(b);
    let b1 = (a < b) as u64;
    let d2 = d.wrapping_sub(borrow_in);
    let b2 = (d < borrow_in) as u64;
    (d2, b1 | b2)
}

/// A 128-bit register as two 64-bit words. All arithmetic wraps mod 2¹²⁸.
///
/// Read as unsigned or as two's-complement signed, whichever the caller defines.
/// The derived ordering is unsigned and is not used by any numeric code path.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct U128 {
    /// High 64 bits (bit 63 is the sign bit when read as signed).
    pub hi: u64,
    /// Low 64 bits.
    pub lo: u64,
}

impl U128 {
    /// Zero.
    pub const ZERO: Self = Self { hi: 0, lo: 0 };

    /// The register `hi · 2⁶⁴ + lo`.
    pub const fn new(hi: u64, lo: u64) -> Self {
        Self { hi, lo }
    }

    /// `self + rhs` mod 2¹²⁸.
    pub const fn wrapping_add(self, rhs: Self) -> Self {
        let (lo, c) = add_with_carry(self.lo, rhs.lo, 0);
        let (hi, _) = add_with_carry(self.hi, rhs.hi, c);
        Self { hi, lo }
    }

    /// `self − rhs` mod 2¹²⁸.
    pub const fn wrapping_sub(self, rhs: Self) -> Self {
        let (lo, b) = sub_with_borrow(self.lo, rhs.lo, 0);
        let (hi, _) = sub_with_borrow(self.hi, rhs.hi, b);
        Self { hi, lo }
    }

    /// `−self` mod 2¹²⁸.
    pub const fn wrapping_neg(self) -> Self {
        Self::ZERO.wrapping_sub(self)
    }

    /// `self · m` mod 2¹²⁸ for an unsigned 64-bit multiplier.
    ///
    /// Correct for a signed (two's-complement) `self` too, because multiplication
    /// mod 2¹²⁸ does not depend on the interpretation.
    pub const fn wrapping_mul_u64(self, m: u64) -> Self {
        // (hi·2⁶⁴ + lo)·m = lo·m + (hi·m)·2⁶⁴; the second term keeps only its low word.
        let lo = self.lo.wrapping_mul(m);
        let carry = mul_hi_u64(self.lo, m);
        let hi = self.hi.wrapping_mul(m).wrapping_add(carry);
        Self { hi, lo }
    }

    /// The full product `a · b` of two u64 values (never wraps).
    pub const fn mul_wide_u64(a: u64, b: u64) -> Self {
        Self {
            hi: mul_hi_u64(a, b),
            lo: a.wrapping_mul(b),
        }
    }

    /// Two's-complement bits of a native `i128` (host-side construction and tests).
    pub const fn from_i128(v: i128) -> Self {
        Self::from_u128(v as u128)
    }

    /// The register as a native `i128` (two's complement, host side and tests).
    pub const fn to_i128(self) -> i128 {
        self.to_u128() as i128
    }

    /// Bits of a native `u128` (host-side construction and tests).
    pub const fn from_u128(v: u128) -> Self {
        // Bit extraction: the shift leaves < 2⁶⁴ and the truncating cast keeps the low word.
        Self {
            hi: (v >> 64) as u64,
            lo: v as u64,
        }
    }

    /// The register as a native `u128` (host side and tests).
    pub const fn to_u128(self) -> u128 {
        ((self.hi as u128) << 64) | (self.lo as u128)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mul_hi_edges() {
        assert_eq!(mul_hi_u64(u64::MAX, u64::MAX), u64::MAX - 1);
        assert_eq!(mul_hi_u64(1 << 32, 1 << 32), 1);
        assert_eq!(mul_hi_u64(u64::MAX, 1), 0);
        assert_eq!(mul_hi_u64(0, u64::MAX), 0);
    }

    #[test]
    fn carries_and_borrows() {
        assert_eq!(add_with_carry(u64::MAX, 1, 0), (0, 1));
        assert_eq!(add_with_carry(u64::MAX, 0, 1), (0, 1));
        assert_eq!(add_with_carry(u64::MAX, u64::MAX, 1), (u64::MAX, 1));
        assert_eq!(sub_with_borrow(0, 1, 0), (u64::MAX, 1));
        assert_eq!(sub_with_borrow(0, 0, 1), (u64::MAX, 1));
        assert_eq!(sub_with_borrow(0, u64::MAX, 1), (0, 1));
        assert_eq!(
            U128::new(0, u64::MAX).wrapping_add(U128::new(0, 1)),
            U128::new(1, 0)
        );
        assert_eq!(U128::ZERO.wrapping_neg(), U128::ZERO);
        assert_eq!(
            U128::new(0, 1).wrapping_neg(),
            U128::new(u64::MAX, u64::MAX)
        );
    }
}
