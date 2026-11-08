//! Exact rational sample rates.

use crate::TimeError;

/// Largest accepted reduced denominator: 2³².
pub const MAX_RATE_DEN: u64 = 1 << 32;
/// Largest accepted reduced numerator: 2⁴⁸.
pub const MAX_RATE_NUM: u64 = 1 << 48;
/// Largest accepted rate in samples per second: 1 GS/s.
pub const MAX_RATE_HZ: u64 = 1_000_000_000;

/// A sample rate `num / den` samples per second, stored reduced (gcd = 1).
///
/// # Bounds (contracts §10, decision 1)
///
/// After reduction a rate is accepted only if `den <= 2³²`, `num <= 2⁴⁸` and
/// `num / den <= 10⁹`. These bounds make every intermediate in
/// [`ChannelTimebase`](crate::ChannelTimebase) provably fit in u128; the argument is
/// written next to each conversion. Any integer rate from 1 S/s to 1 GS/s is accepted,
/// as is any rational with a denominator up to 2³² (for example 10⁹/3).
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub struct SampleRate {
    num: u64,
    den: u64,
}

const fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let r = a % b;
        a = b;
        b = r;
    }
    a
}

impl SampleRate {
    /// A rate of `hz` samples per second.
    pub fn from_hz(hz: u64) -> Result<Self, TimeError> {
        Self::from_ratio(hz, 1)
    }

    /// A rate of `num / den` samples per second. The ratio is reduced before the
    /// bounds are checked, so `122880000/2` is accepted as `61440000/1`.
    pub fn from_ratio(num: u64, den: u64) -> Result<Self, TimeError> {
        if num == 0 || den == 0 {
            return Err(TimeError::ZeroRate);
        }
        let g = gcd(num, den);
        let (num, den) = (num / g, den / g);
        // num <= MAX_RATE_HZ * den, evaluated in u128 so it cannot overflow.
        let within_rate = u128::from(num) <= u128::from(MAX_RATE_HZ) * u128::from(den);
        if den > MAX_RATE_DEN || num > MAX_RATE_NUM || !within_rate {
            return Err(TimeError::RateOutOfBounds { num, den });
        }
        Ok(Self { num, den })
    }

    /// Reduced numerator (samples).
    pub const fn num(self) -> u64 {
        self.num
    }

    /// Reduced denominator (seconds).
    pub const fn den(self) -> u64 {
        self.den
    }

    /// Whether the rate is a whole number of samples per second.
    pub const fn is_integer(self) -> bool {
        self.den == 1
    }
}

/// Orders by numeric value. Rates are stored reduced, so equal values have equal
/// fields and this agrees with the derived `Eq` and `Hash`.
impl Ord for SampleRate {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // Cross-multiplied in u128: each product is below 2^48 · 2^32 = 2^80.
        let lhs = u128::from(self.num) * u128::from(other.den);
        let rhs = u128::from(other.num) * u128::from(self.den);
        lhs.cmp(&rhs)
    }
}

impl PartialOrd for SampleRate {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reduces_and_rejects() {
        let r = SampleRate::from_ratio(122_880_000, 2).unwrap();
        assert_eq!((r.num(), r.den()), (61_440_000, 1));
        assert!(r.is_integer());
        assert_eq!(SampleRate::from_hz(0), Err(TimeError::ZeroRate));
        assert_eq!(SampleRate::from_ratio(1, 0), Err(TimeError::ZeroRate));
        assert!(SampleRate::from_hz(MAX_RATE_HZ).is_ok());
        assert!(matches!(
            SampleRate::from_hz(MAX_RATE_HZ + 1),
            Err(TimeError::RateOutOfBounds { .. })
        ));
        // Reduces to within bounds even though the inputs are large.
        assert!(SampleRate::from_ratio(u64::MAX - 1, u64::MAX - 1).is_ok());
        assert!(matches!(
            SampleRate::from_ratio(1, MAX_RATE_DEN + 1),
            Err(TimeError::RateOutOfBounds { .. })
        ));
        assert!(SampleRate::from_ratio(1_000_000_000, 3).is_ok());
        assert!(SampleRate::from_ratio(1, MAX_RATE_DEN).is_ok());
    }

    #[test]
    fn orders_by_value() {
        let half = SampleRate::from_ratio(1, 2).unwrap();
        let one = SampleRate::from_hz(1).unwrap();
        let third_ghz = SampleRate::from_ratio(1_000_000_000, 3).unwrap();
        assert!(half < one);
        assert!(third_ghz < SampleRate::from_hz(333_333_334).unwrap());
        assert!(third_ghz > SampleRate::from_hz(333_333_333).unwrap());
        assert_eq!(
            one.cmp(&SampleRate::from_ratio(7, 7).unwrap()),
            std::cmp::Ordering::Equal
        );
    }

    #[test]
    fn numerator_bound_is_enforced() {
        // num > 2^48 with a large enough den to stay under 1 GS/s.
        let num = MAX_RATE_NUM + 1; // odd, so coprime with a power of two
        let den = 1 << 20;
        assert!(matches!(
            SampleRate::from_ratio(num, den),
            Err(TimeError::RateOutOfBounds { .. })
        ));
    }
}
