//! Truth time: i64 picoseconds from the scenario epoch, with checked arithmetic.

use crate::TimeError;

/// Picoseconds per second.
pub const PS_PER_SECOND: i64 = 1_000_000_000_000;
/// Picoseconds per millisecond.
pub const PS_PER_MILLISECOND: i64 = 1_000_000_000;
/// Picoseconds per microsecond.
pub const PS_PER_MICROSECOND: i64 = 1_000_000;
/// Picoseconds per nanosecond.
pub const PS_PER_NANOSECOND: i64 = 1_000;

/// An instant of truth time: signed picoseconds from the scenario epoch.
///
/// The positive range is about 106.7 days. All arithmetic is checked: overflow returns
/// [`TimeError::TimeOverflow`] and never wraps (contracts §1.1).
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct TruthTime(i64);

/// A signed span of truth time in picoseconds.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct TruthDuration(i64);

impl TruthTime {
    /// The scenario epoch (0 ps).
    pub const EPOCH: Self = Self(0);
    /// The largest representable instant.
    pub const MAX: Self = Self(i64::MAX);
    /// The smallest representable instant.
    pub const MIN: Self = Self(i64::MIN);

    /// Instant `ps` picoseconds after the scenario epoch.
    pub const fn from_ps(ps: i64) -> Self {
        Self(ps)
    }

    /// Picoseconds after the scenario epoch.
    pub const fn as_ps(self) -> i64 {
        self.0
    }

    /// Instant `s` whole seconds after the epoch.
    pub fn from_secs(s: i64) -> Result<Self, TimeError> {
        TruthDuration::from_secs(s).map(|d| Self(d.0))
    }

    /// `self + d`, or [`TimeError::TimeOverflow`].
    pub fn checked_add(self, d: TruthDuration) -> Result<Self, TimeError> {
        self.0
            .checked_add(d.0)
            .map(Self)
            .ok_or(TimeError::TimeOverflow)
    }

    /// `self - d`, or [`TimeError::TimeOverflow`].
    pub fn checked_sub(self, d: TruthDuration) -> Result<Self, TimeError> {
        self.0
            .checked_sub(d.0)
            .map(Self)
            .ok_or(TimeError::TimeOverflow)
    }

    /// The span `self - earlier`, or [`TimeError::TimeOverflow`] if it does not fit in i64.
    pub fn checked_since(self, earlier: TruthTime) -> Result<TruthDuration, TimeError> {
        self.0
            .checked_sub(earlier.0)
            .map(TruthDuration)
            .ok_or(TimeError::TimeOverflow)
    }
}

impl TruthDuration {
    /// Zero span.
    pub const ZERO: Self = Self(0);

    /// A span of `ps` picoseconds.
    pub const fn from_ps(ps: i64) -> Self {
        Self(ps)
    }

    /// The span in picoseconds.
    pub const fn as_ps(self) -> i64 {
        self.0
    }

    /// A span of `ns` nanoseconds, checked.
    pub fn from_ns(ns: i64) -> Result<Self, TimeError> {
        Self::scaled(ns, PS_PER_NANOSECOND)
    }

    /// A span of `us` microseconds, checked.
    pub fn from_us(us: i64) -> Result<Self, TimeError> {
        Self::scaled(us, PS_PER_MICROSECOND)
    }

    /// A span of `ms` milliseconds, checked.
    pub fn from_ms(ms: i64) -> Result<Self, TimeError> {
        Self::scaled(ms, PS_PER_MILLISECOND)
    }

    /// A span of `s` seconds, checked.
    pub fn from_secs(s: i64) -> Result<Self, TimeError> {
        Self::scaled(s, PS_PER_SECOND)
    }

    fn scaled(value: i64, unit: i64) -> Result<Self, TimeError> {
        value
            .checked_mul(unit)
            .map(Self)
            .ok_or(TimeError::TimeOverflow)
    }

    /// `self + other`, checked.
    pub fn checked_add(self, other: Self) -> Result<Self, TimeError> {
        self.0
            .checked_add(other.0)
            .map(Self)
            .ok_or(TimeError::TimeOverflow)
    }

    /// `self - other`, checked.
    pub fn checked_sub(self, other: Self) -> Result<Self, TimeError> {
        self.0
            .checked_sub(other.0)
            .map(Self)
            .ok_or(TimeError::TimeOverflow)
    }

    /// `self * k`, checked.
    pub fn checked_mul(self, k: i64) -> Result<Self, TimeError> {
        self.0
            .checked_mul(k)
            .map(Self)
            .ok_or(TimeError::TimeOverflow)
    }

    /// `-self`, checked (`i64::MIN` has no negation).
    pub fn checked_neg(self) -> Result<Self, TimeError> {
        self.0
            .checked_neg()
            .map(Self)
            .ok_or(TimeError::TimeOverflow)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn overflow_is_an_error() {
        assert_eq!(
            TruthTime::MAX.checked_add(TruthDuration::from_ps(1)),
            Err(TimeError::TimeOverflow)
        );
        assert_eq!(
            TruthTime::MIN.checked_sub(TruthDuration::from_ps(1)),
            Err(TimeError::TimeOverflow)
        );
        assert_eq!(
            TruthTime::MAX.checked_since(TruthTime::from_ps(-1)),
            Err(TimeError::TimeOverflow)
        );
        assert_eq!(
            TruthDuration::from_ps(i64::MIN).checked_neg(),
            Err(TimeError::TimeOverflow)
        );
        assert_eq!(
            TruthDuration::from_secs(9_300_000),
            Err(TimeError::TimeOverflow)
        );
    }

    #[test]
    fn hundred_days_fits() {
        let d = TruthDuration::from_secs(100 * 86_400).unwrap();
        assert_eq!(d.as_ps(), 8_640_000_000_000_000_000);
        let t = TruthTime::EPOCH.checked_add(d).unwrap();
        assert_eq!(t.checked_since(TruthTime::EPOCH), Ok(d));
    }

    #[test]
    fn unit_constructors() {
        assert_eq!(TruthDuration::from_ns(3).unwrap().as_ps(), 3_000);
        assert_eq!(TruthDuration::from_us(3).unwrap().as_ps(), 3_000_000);
        assert_eq!(TruthDuration::from_ms(-3).unwrap().as_ps(), -3_000_000_000);
        assert_eq!(TruthTime::from_secs(2).unwrap().as_ps(), 2 * PS_PER_SECOND);
    }
}
