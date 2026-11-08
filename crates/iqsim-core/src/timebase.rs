//! Sample indices and the exact rational map between sample index and truth time.
//!
//! # Overflow argument (contracts §1.2)
//!
//! Let `M = den · 10¹²`. With the [`SampleRate`] bounds `den <= 2³²`, `num <= 2⁴⁸` and
//! `10¹² < 2⁴⁰`:
//!
//! - **Index → time.** `n · den < 2⁶⁴ · 2³² = 2⁹⁶`, so the first product always fits.
//!   The second product `n · den · 10¹²` is checked. For every run of at most 100 days
//!   at at most 1 GS/s, `n <= 8.64·10¹⁵ < 2⁵³`, so it is below `2⁵³⁺³²⁺⁴⁰ = 2¹²⁵` and never
//!   overflows. When it does overflow (only for far larger `n`), the product is at least
//!   2¹²⁸, so the quotient is at least `2¹²⁸ / 2⁴⁸ = 2⁸⁰`, which cannot be a valid i64
//!   picosecond time anyway: the error is [`TimeError::TimeOverflow`] either way, so the
//!   check never rejects a representable result.
//! - **Time → index.** `t - epoch` is between `-2⁶⁴` and `2⁶⁴` and is computed in i128.
//!   For non-negative spans `d < 2⁶⁴`, `d · num < 2¹¹²` and `M < 2⁷²`, so both fit in u128.
//!   The quotient is at most `d · 10⁹ / 10¹² < 2⁵⁴`, so it always fits in u64.
//! - **Remainder → binary fraction.** The remainder `r < num <= 2⁴⁸`, so `r · 2⁶⁴ < 2¹¹²`.

use crate::{AcqInstant, SampleRate, TimeError, TruthTime};

const PS_PER_SECOND_U128: u128 = 1_000_000_000_000;

/// A sample index on one receiver channel, counted from that channel's epoch.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct SampleIndex(u64);

impl SampleIndex {
    /// Index zero (the channel epoch).
    pub const ZERO: Self = Self(0);

    /// The index `n`.
    pub const fn new(n: u64) -> Self {
        Self(n)
    }

    /// The raw index.
    pub const fn get(self) -> u64 {
        self.0
    }

    /// `self + k`, or [`TimeError::IndexOverflow`].
    pub fn checked_add(self, k: u64) -> Result<Self, TimeError> {
        self.0
            .checked_add(k)
            .map(Self)
            .ok_or(TimeError::IndexOverflow)
    }

    /// `self - k`, or [`TimeError::IndexOverflow`] below zero.
    pub fn checked_sub(self, k: u64) -> Result<Self, TimeError> {
        self.0
            .checked_sub(k)
            .map(Self)
            .ok_or(TimeError::IndexOverflow)
    }
}

/// The exact truth time of a sample: `floor` picoseconds plus `rem / rem_den` of a
/// picosecond, with `0 <= rem < rem_den`.
///
/// Only [`ChannelTimebase::index_to_time`] constructs it, so the invariant always holds.
///
/// `rem_den` is the reduced rate numerator, so the value is exact for every index. It is
/// [`exact`](Self::is_exact) as an integer picosecond time when `rem == 0`.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ExactTime {
    floor: TruthTime,
    rem: u64,
    rem_den: u64,
}

impl ExactTime {
    /// Floor of the exact time, in picoseconds.
    pub const fn floor(self) -> TruthTime {
        self.floor
    }

    /// Remainder numerator, `0 <= rem < rem_den`.
    pub const fn rem(self) -> u64 {
        self.rem
    }

    /// Remainder denominator: the reduced sample-rate numerator, so nonzero.
    pub const fn rem_den(self) -> u64 {
        self.rem_den
    }

    /// Whether the time is an exact integer number of picoseconds.
    pub const fn is_exact(self) -> bool {
        self.rem == 0
    }

    /// The time as an [`AcqInstant`]: the remainder is converted to a 2⁻⁶⁴ ps fraction,
    /// rounded toward −∞. The conversion error is below 2⁻⁶⁴ ps (≈ 5.4·10⁻²⁰ ps).
    pub fn to_acq_instant(self) -> AcqInstant {
        // rem < rem_den <= 2^48, so rem << 64 < 2^112: no overflow. The quotient is < 2^64.
        let frac = (u128::from(self.rem) << 64) / u128::from(self.rem_den);
        let frac = u64::try_from(frac).expect("rem < rem_den implies frac < 2^64");
        AcqInstant::new(self.floor, frac)
    }
}

/// The ideal timebase of one channel: its epoch in truth time and its exact sample rate.
///
/// Sample `n` is acquired at exactly `epoch + n · den / num` seconds. This is the
/// reference map behind [`IdealClockMap`](crate::IdealClockMap); perturbed clocks
/// (prompt 15) do not use it directly.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub struct ChannelTimebase {
    epoch: TruthTime,
    rate: SampleRate,
}

impl ChannelTimebase {
    /// A timebase whose sample 0 is at `epoch`.
    pub const fn new(epoch: TruthTime, rate: SampleRate) -> Self {
        Self { epoch, rate }
    }

    /// Truth time of sample 0.
    pub const fn epoch(self) -> TruthTime {
        self.epoch
    }

    /// The sample rate.
    pub const fn rate(self) -> SampleRate {
        self.rate
    }

    /// `M = den · 10¹²`, the picoseconds in `den` seconds. Below 2⁷².
    fn ps_per_num_samples(self) -> u128 {
        u128::from(self.rate.den()) * PS_PER_SECOND_U128
    }

    /// Exact truth time of sample `n`: floor in ps plus the exact remainder.
    ///
    /// Returns [`TimeError::TimeOverflow`] if the floor does not fit in i64 ps.
    pub fn index_to_time(self, n: SampleIndex) -> Result<ExactTime, TimeError> {
        let num = u128::from(self.rate.num());
        // n * den < 2^96: cannot overflow. The 10^12 product is checked (module docs).
        let product = (u128::from(n.get()) * u128::from(self.rate.den()))
            .checked_mul(PS_PER_SECOND_U128)
            .ok_or(TimeError::TimeOverflow)?;
        let q = product / num;
        let rem = u64::try_from(product % num).expect("remainder < num <= 2^48");
        let q = i128::from(u64::try_from(q).map_err(|_| TimeError::TimeOverflow)?);
        let floor = i64::try_from(i128::from(self.epoch.as_ps()) + q)
            .map_err(|_| TimeError::TimeOverflow)?;
        Ok(ExactTime {
            floor: TruthTime::from_ps(floor),
            rem,
            rem_den: self.rate.num(),
        })
    }

    /// Largest `n` whose exact time is `<= t` (floor). [`TimeError::BeforeEpoch`] if
    /// `t` is before sample 0.
    pub fn index_at_or_before(self, t: TruthTime) -> Result<SampleIndex, TimeError> {
        let d = self.span_from_epoch(t)?;
        // d < 2^64, num <= 2^48: d * num < 2^112.
        let n = d * u128::from(self.rate.num()) / self.ps_per_num_samples();
        u64::try_from(n)
            .map(SampleIndex)
            .map_err(|_| TimeError::IndexOverflow)
    }

    /// Smallest `n` whose exact time is `>= t` (ceiling). Times at or before the epoch
    /// map to index 0.
    pub fn index_at_or_after(self, t: TruthTime) -> Result<SampleIndex, TimeError> {
        let d = match self.span_from_epoch(t) {
            Ok(d) => d,
            Err(TimeError::BeforeEpoch) => return Ok(SampleIndex::ZERO),
            Err(e) => return Err(e),
        };
        let m = self.ps_per_num_samples();
        let n = (d * u128::from(self.rate.num())).div_ceil(m);
        u64::try_from(n)
            .map(SampleIndex)
            .map_err(|_| TimeError::IndexOverflow)
    }

    /// `t - epoch` in ps as a non-negative u128 (always below 2⁶⁴).
    fn span_from_epoch(self, t: TruthTime) -> Result<u128, TimeError> {
        let d = i128::from(t.as_ps()) - i128::from(self.epoch.as_ps());
        u128::try_from(d).map_err(|_| TimeError::BeforeEpoch)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tb(epoch_ps: i64, num: u64, den: u64) -> ChannelTimebase {
        ChannelTimebase::new(
            TruthTime::from_ps(epoch_ps),
            SampleRate::from_ratio(num, den).unwrap(),
        )
    }

    #[test]
    fn remainder_is_exact() {
        // 61.44 MS/s: one sample is 16276 + 1/24 ps.
        let t = tb(0, 61_440_000, 1)
            .index_to_time(SampleIndex::new(1))
            .unwrap();
        assert_eq!(t.floor.as_ps(), 16_276);
        assert_eq!(t.rem * 24, t.rem_den);
        assert!(!t.is_exact());
        // The fraction is floor(2^64 / 24).
        assert_eq!(t.to_acq_instant().frac(), u64::MAX / 24);
    }

    #[test]
    fn epoch_offset_and_overflow() {
        let base = tb(-5, 250_000_000, 1);
        assert_eq!(
            base.index_to_time(SampleIndex::ZERO).unwrap().floor.as_ps(),
            -5
        );
        assert_eq!(
            base.index_at_or_before(TruthTime::from_ps(-6)),
            Err(TimeError::BeforeEpoch)
        );
        assert_eq!(
            base.index_at_or_after(TruthTime::from_ps(-6)),
            Ok(SampleIndex::ZERO)
        );
        assert_eq!(
            base.index_to_time(SampleIndex::new(u64::MAX)),
            Err(TimeError::TimeOverflow)
        );
        // Negative epoch lets q exceed i64::MAX - |epoch| only up to the i64 limit.
        let neg = tb(i64::MIN, 1_000_000_000, 1);
        let last = neg
            .index_at_or_before(TruthTime::MAX)
            .expect("a full 2^64 ps span fits");
        assert!(neg.index_to_time(last).is_ok());
        assert_eq!(
            neg.index_to_time(last.checked_add(1).unwrap()),
            Err(TimeError::TimeOverflow)
        );
    }

    #[test]
    fn sample_index_checked() {
        assert_eq!(
            SampleIndex::new(u64::MAX).checked_add(1),
            Err(TimeError::IndexOverflow)
        );
        assert_eq!(
            SampleIndex::ZERO.checked_sub(1),
            Err(TimeError::IndexOverflow)
        );
    }
}
