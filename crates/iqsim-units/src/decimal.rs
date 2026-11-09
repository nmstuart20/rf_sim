//! Exact decimal numbers parsed from text.
//!
//! A quantity's number is kept as `±mantissa · 10^exp10` with a u128 mantissa, so
//! exact quantities (frequency, sample rate, duration, data rate) convert to their
//! integer representation with no binary rounding. Float quantities convert once,
//! through Rust's correctly rounded `str → f64` parser, which gives the same bits on
//! every platform.

/// Most significant digits a mantissa may hold (10³⁸ < 2¹²⁸).
pub(crate) const MAX_DIGITS: usize = 38;
/// Largest accepted written exponent magnitude.
const MAX_EXPONENT: i32 = 10_000;

/// `±mantissa · 10^exp10`, normalized so the mantissa has no trailing zeros
/// (zero is `+0 · 10⁰`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Decimal {
    negative: bool,
    mantissa: u128,
    exp10: i32,
}

/// Why a decimal could not be scaled to an integer.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ScaleError {
    /// Nonzero digits remain below the requested resolution.
    Inexact,
    /// The scaled value does not fit.
    Overflow,
}

/// Length in bytes of the longest prefix of `s` that is a number:
/// `[+-]? (digits ('.' digits*)? | '.' digits) ([eE] [+-]? digits)?`.
pub(crate) fn number_len(s: &str) -> usize {
    let b = s.as_bytes();
    let mut i = 0;
    if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
        i += 1;
    }
    let int_start = i;
    while i < b.len() && b[i].is_ascii_digit() {
        i += 1;
    }
    let mut digits = i - int_start;
    if i < b.len() && b[i] == b'.' {
        let frac_start = i + 1;
        let mut j = frac_start;
        while j < b.len() && b[j].is_ascii_digit() {
            j += 1;
        }
        if digits > 0 || j > frac_start {
            digits += j - frac_start;
            i = j;
        }
    }
    if digits == 0 {
        return 0;
    }
    if i < b.len() && (b[i] == b'e' || b[i] == b'E') {
        let mut j = i + 1;
        if j < b.len() && (b[j] == b'+' || b[j] == b'-') {
            j += 1;
        }
        let exp_start = j;
        while j < b.len() && b[j].is_ascii_digit() {
            j += 1;
        }
        if j > exp_start {
            i = j;
        }
    }
    i
}

/// Splits `s` (already trimmed) into its leading number and the rest, with the
/// rest's leading whitespace removed.
pub(crate) fn split_number(s: &str) -> (&str, &str) {
    let n = number_len(s);
    (&s[..n], s[n..].trim_start())
}

const fn pow10(k: u32) -> Option<u128> {
    10u128.checked_pow(k)
}

impl Decimal {
    /// The value zero.
    pub(crate) const ZERO: Self = Self {
        negative: false,
        mantissa: 0,
        exp10: 0,
    };

    /// Parses a complete number. Returns `None` unless all of `text` is a number with
    /// at most [`MAX_DIGITS`] significant digits and a written exponent of at most
    /// 10 000 in magnitude.
    pub(crate) fn parse(text: &str) -> Option<Self> {
        if text.is_empty() || number_len(text) != text.len() {
            return None;
        }
        let (negative, body) = match text.as_bytes()[0] {
            b'-' => (true, &text[1..]),
            b'+' => (false, &text[1..]),
            _ => (false, text),
        };
        let (significand, exponent) = match body.find(['e', 'E']) {
            Some(at) => (&body[..at], body[at + 1..].parse::<i32>().ok()?),
            None => (body, 0),
        };
        if !(-MAX_EXPONENT..=MAX_EXPONENT).contains(&exponent) {
            return None;
        }
        let (int_part, frac_part) = significand.split_once('.').unwrap_or((significand, ""));
        let frac_len = i32::try_from(frac_part.len()).ok()?;
        let digits = format!("{int_part}{frac_part}");
        let digits = digits.trim_start_matches('0');
        let trimmed = digits.trim_end_matches('0');
        if trimmed.is_empty() {
            return Some(Self::ZERO);
        }
        if trimmed.len() > MAX_DIGITS {
            return None;
        }
        let trailing = i32::try_from(digits.len() - trimmed.len()).ok()?;
        Some(Self {
            negative,
            mantissa: trimmed.parse().ok()?,
            exp10: exponent - frac_len + trailing,
        })
    }

    /// A non-negative integer as a decimal.
    pub(crate) fn from_u128(value: u128) -> Self {
        Self::normalized(false, value, 0)
    }

    fn normalized(negative: bool, mut mantissa: u128, mut exp10: i32) -> Self {
        if mantissa == 0 {
            return Self::ZERO;
        }
        while mantissa.is_multiple_of(10) {
            mantissa /= 10;
            exp10 += 1;
        }
        Self {
            negative,
            mantissa,
            exp10,
        }
    }

    /// Whether the value is below zero.
    pub(crate) fn is_negative(self) -> bool {
        self.negative
    }

    /// Whether the value is zero.
    pub(crate) fn is_zero(self) -> bool {
        self.mantissa == 0
    }

    /// Whether the value is a whole number.
    pub(crate) fn is_integer(self) -> bool {
        self.exp10 >= 0 || self.is_zero()
    }

    /// The value times `k`, exactly. `None` if the mantissa would overflow.
    pub(crate) fn mul_u128(self, k: u128) -> Option<Self> {
        Some(Self::normalized(
            self.negative,
            self.mantissa.checked_mul(k)?,
            self.exp10,
        ))
    }

    /// The value times `10^k`, exactly. `None` if the exponent overflows.
    pub(crate) fn shifted(self, k: i32) -> Option<Self> {
        if self.is_zero() {
            return Some(self);
        }
        Some(Self {
            exp10: self.exp10.checked_add(k)?,
            ..self
        })
    }

    /// Whether the value is below the integer `n`, exactly.
    pub(crate) fn less_than(self, n: u128) -> bool {
        if self.negative {
            return true;
        }
        if self.is_zero() {
            return n > 0;
        }
        let shift = self.exp10.unsigned_abs();
        if self.exp10 >= 0 {
            // mantissa · 10^e < n; an overflowing left side is certainly not below n.
            pow10(shift)
                .and_then(|p| self.mantissa.checked_mul(p))
                .is_some_and(|lhs| lhs < n)
        } else {
            // mantissa < n · 10^-e; an overflowing right side is certainly above.
            pow10(shift)
                .and_then(|p| n.checked_mul(p))
                .is_none_or(|rhs| self.mantissa < rhs)
        }
    }

    /// `self + other` for non-negative operands, exactly. `None` on overflow.
    pub(crate) fn add_nonnegative(self, other: Self) -> Option<Self> {
        debug_assert!(!self.negative && !other.negative);
        let exp10 = self.exp10.min(other.exp10);
        let align = |d: Self| -> Option<u128> {
            let shift = u32::try_from(d.exp10 - exp10).ok()?;
            d.mantissa.checked_mul(pow10(shift)?)
        };
        let lhs = if self.is_zero() { 0 } else { align(self)? };
        let rhs = if other.is_zero() { 0 } else { align(other)? };
        Some(Self::normalized(false, lhs.checked_add(rhs)?, exp10))
    }

    /// `value · 10^scale` as an integer, exactly.
    pub(crate) fn scaled_i128(self, scale: i32) -> Result<i128, ScaleError> {
        if self.is_zero() {
            return Ok(0);
        }
        let e = self.exp10.checked_add(scale).ok_or(ScaleError::Overflow)?;
        let magnitude = if e >= 0 {
            let factor = u32::try_from(e)
                .ok()
                .and_then(pow10)
                .ok_or(ScaleError::Overflow)?;
            self.mantissa
                .checked_mul(factor)
                .ok_or(ScaleError::Overflow)?
        } else {
            // A normalized mantissa has no trailing zeros, so a negative exponent
            // always leaves a fraction.
            return Err(ScaleError::Inexact);
        };
        let magnitude = i128::try_from(magnitude).map_err(|_| ScaleError::Overflow)?;
        Ok(if self.negative { -magnitude } else { magnitude })
    }

    /// `value · 10^scale` as a reduced-by-powers-of-ten ratio `num / 10^k` of
    /// non-negative integers. `None` if negative or too large.
    pub(crate) fn ratio(self, scale: i32) -> Option<(u128, u128)> {
        if self.negative {
            return None;
        }
        if self.is_zero() {
            return Some((0, 1));
        }
        let e = self.exp10.checked_add(scale)?;
        if e >= 0 {
            let factor = pow10(u32::try_from(e).ok()?)?;
            Some((self.mantissa.checked_mul(factor)?, 1))
        } else {
            Some((self.mantissa, pow10(e.unsigned_abs())?))
        }
    }

    /// The nearest f64 (ties to even), through the correctly rounded std parser.
    /// May be infinite or zero for extreme exponents; callers range-check.
    pub(crate) fn to_f64(self) -> f64 {
        let sign = if self.negative { "-" } else { "" };
        format!("{sign}{}e{}", self.mantissa, self.exp10)
            .parse()
            .expect("a formatted decimal always parses as f64")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn d(s: &str) -> Decimal {
        Decimal::parse(s).unwrap()
    }

    #[test]
    fn number_prefixes() {
        for (s, n) in [
            ("1090 MHz", 4),
            ("1090MHz", 4),
            ("-1.5e3 m", 6),
            ("+.5us", 3),
            ("5.", 2),
            ("1e", 1),
            ("1e+", 1),
            ("1E-9s", 4),
            (".", 0),
            ("-", 0),
            ("MHz", 0),
            ("FL350", 0),
            ("", 0),
        ] {
            assert_eq!(number_len(s), n, "{s}");
        }
    }

    #[test]
    fn parses_and_normalizes() {
        assert_eq!(d("1090"), d("1.09e3"));
        assert_eq!(d("61.44"), d("6144e-2"));
        assert_eq!(d("0.000"), Decimal::ZERO);
        assert_eq!(d("-0"), Decimal::ZERO);
        assert_eq!(d("1.50").scaled_i128(1), Ok(15));
        assert_eq!(d("-2.5").scaled_i128(1), Ok(-25));
        assert_eq!(d("2.5").scaled_i128(0), Err(ScaleError::Inexact));
        assert_eq!(d("1e40").scaled_i128(0), Err(ScaleError::Overflow));
        assert_eq!(d("1e-10000").scaled_i128(0), Err(ScaleError::Inexact));
        assert_eq!(d("61.44").ratio(6), Some((61_440_000, 1)));
        assert_eq!(d("1.5").ratio(-1), Some((15, 100)));
        assert_eq!(d("-1").ratio(0), None);
        assert!(Decimal::parse("1.2.3").is_none());
        assert!(Decimal::parse("1e10001").is_none());
        assert!(Decimal::parse(&"1".repeat(39)).is_none());
        assert!(Decimal::parse(&format!("{}000", "9".repeat(38))).is_some());
        assert!(Decimal::parse(&format!("0.000{}", "7".repeat(38))).is_some());
    }

    #[test]
    fn arithmetic_is_exact() {
        assert_eq!(d("0.3048").mul_u128(100).unwrap(), d("30.48"));
        assert_eq!(
            d("56.7").add_nonnegative(Decimal::from_u128(3600)).unwrap(),
            d("3656.7")
        );
        assert_eq!(Decimal::ZERO.add_nonnegative(d("1.25")).unwrap(), d("1.25"));
        assert!(d("-12.5").is_negative());
        assert!(!d("-0.0").is_negative());
        assert!(d("59.999").less_than(60));
        assert!(!d("60").less_than(60));
        assert!(!d("6e1").less_than(60));
        assert!(!d("60.0001").less_than(60));
        assert!(d("0").less_than(60));
        assert!(!d("1e300").less_than(60));
        assert!(d("1e-300").less_than(60));
    }

    #[test]
    fn f64_conversion_is_correctly_rounded() {
        assert_eq!(d("0.1").to_f64(), 0.1);
        assert_eq!(d("-1090.5e6").to_f64(), -1_090_500_000.0);
        assert_eq!(d("1e400").to_f64(), f64::INFINITY);
        assert_eq!(d("1e-400").to_f64(), 0.0);
    }
}
