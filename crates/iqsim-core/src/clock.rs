//! Sampling clocks ([`ClockMap`]) and timestamp clocks ([`TimestampClock`]).
//!
//! These are distinct models (contracts §3). A [`ClockMap`] says *when* a sample was
//! actually acquired; geometry, delay and carrier phase are evaluated at that instant.
//! A [`TimestampClock`] only says what time the receiver *reports* for an instant. It
//! never moves sampling instants, and its output type ([`ReportedTime`]) cannot be used
//! where an [`AcqInstant`] is expected. The LO model is separate again and arrives with
//! the receiver chain.

use std::fmt;

use sha2::{Digest, Sha256};

use crate::{ChannelTimebase, SampleIndex, TimeError, TruthTime};

/// An acquisition instant in truth time: `ps + frac / 2⁶⁴` picoseconds from the
/// scenario epoch (contracts §10, decision 2).
///
/// `ps` is the floor and `frac` the sub-picosecond fraction in units of 2⁻⁶⁴ ps, so
/// the value is a signed Q64.64 picosecond count. Clock maps round the exact instant
/// toward −∞ onto this grid; the error is below 2⁻⁶⁴ ps. The derived ordering is the
/// numeric ordering.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct AcqInstant {
    ps: TruthTime,
    frac: u64,
}

impl AcqInstant {
    /// The instant `ps + frac · 2⁻⁶⁴` picoseconds.
    pub const fn new(ps: TruthTime, frac: u64) -> Self {
        Self { ps, frac }
    }

    /// The integer-picosecond instant `t`.
    pub const fn from_truth_time(t: TruthTime) -> Self {
        Self { ps: t, frac: 0 }
    }

    /// Floor of the instant, in picoseconds.
    pub const fn ps(self) -> TruthTime {
        self.ps
    }

    /// Sub-picosecond fraction, in units of 2⁻⁶⁴ ps.
    pub const fn frac(self) -> u64 {
        self.frac
    }

    /// The instant as a signed Q64.64 picosecond count. Exact, and every i128 value is
    /// a valid instant, so differences of instants are exact i128 subtractions (checked
    /// by the caller).
    pub const fn to_q64_64(self) -> i128 {
        ((self.ps.as_ps() as i128) << 64) | (self.frac as i128)
    }

    /// The inverse of [`to_q64_64`](Self::to_q64_64).
    pub const fn from_q64_64(v: i128) -> Self {
        // An arithmetic shift of any i128 by 64 fits in i64; the low 64 bits are the
        // fraction. Both casts are exact bit extractions.
        Self {
            ps: TruthTime::from_ps((v >> 64) as i64),
            frac: v as u64,
        }
    }
}

/// Stable identity of a clock model and its parameters: SHA-256 of a canonical,
/// versioned byte encoding. Used in cache keys and run manifests.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct ClockMapId(pub [u8; 32]);

impl fmt::Display for ClockMapId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for b in self.0 {
            write!(f, "{b:02x}")?;
        }
        Ok(())
    }
}

/// Map from a channel's sample index to its actual acquisition instant (contracts §3.1).
///
/// Implementations must be deterministic and strictly increasing in the index, so the
/// inverses are well defined. Render code always goes through this trait, never through
/// a nominal-rate shortcut.
pub trait ClockMap {
    /// Actual acquisition instant of sample `n`. Strictly increasing in `n`.
    ///
    /// Returns [`TimeError::TimeOverflow`] if the instant is not representable.
    fn instant(&self, n: SampleIndex) -> Result<AcqInstant, TimeError>;

    /// Smallest `n` with `instant(n) >= t`.
    ///
    /// Returns [`TimeError::TimeOverflow`] or [`TimeError::IndexOverflow`] if no
    /// representable sample is at or after `t`.
    fn index_at_or_after(&self, t: AcqInstant) -> Result<SampleIndex, TimeError>;

    /// Largest `n` with `instant(n) <= t`, or `None` if `t` is before sample 0.
    fn index_at_or_before(&self, t: AcqInstant) -> Result<Option<SampleIndex>, TimeError>;

    /// Stable identity of the model and its parameters.
    fn id(&self) -> ClockMapId;
}

/// The ideal sampling clock: the exact rational map of a [`ChannelTimebase`], with each
/// instant rounded toward −∞ to the 2⁻⁶⁴ ps grid of [`AcqInstant`].
///
/// Rounding preserves strict monotonicity: samples are at least 1000 ps apart (rate
/// at most 1 GS/s), far more than one grid step.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub struct IdealClockMap {
    timebase: ChannelTimebase,
}

impl IdealClockMap {
    /// Domain tag at the start of the identity encoding. Bump the version if the
    /// encoding or the model's rounding ever changes.
    const ID_TAG: &'static [u8] = b"iqsim.clockmap.ideal.v1\0";

    /// The ideal clock of `timebase`.
    pub const fn new(timebase: ChannelTimebase) -> Self {
        Self { timebase }
    }

    /// The underlying timebase.
    pub const fn timebase(self) -> ChannelTimebase {
        self.timebase
    }

    /// First-guess index for `t`: the floor index of `t`'s integer picoseconds, whose
    /// instant is `<= t`. `None` if `t` is before the epoch.
    fn floor_guess(&self, t: AcqInstant) -> Result<Option<SampleIndex>, TimeError> {
        match self.timebase.index_at_or_before(t.ps()) {
            Ok(n) => Ok(Some(n)),
            Err(TimeError::BeforeEpoch) => Ok(None),
            Err(e) => Err(e),
        }
    }
}

/// Returns `Some(instant(n + 1))`, or `None` if that sample is not representable.
/// An unrepresentable sample lies beyond every representable instant.
fn next_instant<C: ClockMap + ?Sized>(
    c: &C,
    n: SampleIndex,
) -> Result<Option<(SampleIndex, AcqInstant)>, TimeError> {
    let Ok(next) = n.checked_add(1) else {
        return Ok(None);
    };
    match c.instant(next) {
        Ok(i) => Ok(Some((next, i))),
        Err(TimeError::TimeOverflow) => Ok(None),
        Err(e) => Err(e),
    }
}

impl ClockMap for IdealClockMap {
    fn instant(&self, n: SampleIndex) -> Result<AcqInstant, TimeError> {
        Ok(self.timebase.index_to_time(n)?.to_acq_instant())
    }

    fn index_at_or_after(&self, t: AcqInstant) -> Result<SampleIndex, TimeError> {
        let Some(mut n) = self.floor_guess(t)? else {
            // t < epoch = instant(0).
            return Ok(SampleIndex::ZERO);
        };
        // instant(n) <= t.ps <= t. Samples are >= 1000 ps apart and t < t.ps + 1, so the
        // answer is n (if equal), n + 1 or, at most, n + 2.
        if self.instant(n)? >= t {
            return Ok(n);
        }
        loop {
            match next_instant(self, n)? {
                Some((next, i)) if i >= t => return Ok(next),
                Some((next, _)) => n = next,
                None => return Err(TimeError::TimeOverflow),
            }
        }
    }

    fn index_at_or_before(&self, t: AcqInstant) -> Result<Option<SampleIndex>, TimeError> {
        let Some(mut n) = self.floor_guess(t)? else {
            return Ok(None);
        };
        // instant(n) <= t. Advance while the next sample is still <= t (at most once).
        while let Some((next, i)) = next_instant(self, n)? {
            if i > t {
                break;
            }
            n = next;
        }
        Ok(Some(n))
    }

    fn id(&self) -> ClockMapId {
        let rate = self.timebase.rate();
        let mut h = Sha256::new();
        h.update(Self::ID_TAG);
        h.update(self.timebase.epoch().as_ps().to_le_bytes());
        h.update(rate.num().to_le_bytes());
        h.update(rate.den().to_le_bytes());
        ClockMapId(h.finalize().into())
    }
}

/// A time as reported by a receiver's timestamp clock: `ps + frac / 2⁶⁴` picoseconds
/// on the receiver's own time scale.
///
/// This is a label for consumers (stream headers, VITA-49, SigMF). It is deliberately
/// a different type from [`AcqInstant`] so it cannot drive geometry or phase.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct ReportedTime {
    /// Floor of the reported time, in picoseconds.
    pub ps: i64,
    /// Sub-picosecond fraction, in units of 2⁻⁶⁴ ps.
    pub frac: u64,
}

/// Labels acquisition instants with the time the receiver reports (contracts §3.2).
///
/// It has its own offset and drift (prompt 15) and never changes when samples are taken.
pub trait TimestampClock {
    /// The reported time for a sample acquired at `t`.
    fn label(&self, t: AcqInstant) -> Result<ReportedTime, TimeError>;

    /// Stable identity of the model and its parameters.
    fn id(&self) -> ClockMapId;
}

/// A perfect timestamp clock: reports truth time unchanged.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Hash)]
pub struct IdealTimestampClock;

impl TimestampClock for IdealTimestampClock {
    fn label(&self, t: AcqInstant) -> Result<ReportedTime, TimeError> {
        Ok(ReportedTime {
            ps: t.ps().as_ps(),
            frac: t.frac(),
        })
    }

    fn id(&self) -> ClockMapId {
        ClockMapId(Sha256::digest(b"iqsim.timestamp.ideal.v1\0").into())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::SampleRate;

    fn map(epoch: i64, num: u64, den: u64) -> IdealClockMap {
        IdealClockMap::new(ChannelTimebase::new(
            TruthTime::from_ps(epoch),
            SampleRate::from_ratio(num, den).unwrap(),
        ))
    }

    #[test]
    fn q64_64_round_trip() {
        for v in [
            AcqInstant::new(TruthTime::MIN, 0),
            AcqInstant::new(TruthTime::from_ps(-1), u64::MAX),
            AcqInstant::new(TruthTime::from_ps(7), 3),
            AcqInstant::new(TruthTime::MAX, u64::MAX),
        ] {
            assert_eq!(AcqInstant::from_q64_64(v.to_q64_64()), v);
        }
        let a = AcqInstant::new(TruthTime::from_ps(-1), u64::MAX);
        let b = AcqInstant::new(TruthTime::from_ps(0), 0);
        assert!(a < b);
        assert_eq!(b.to_q64_64() - a.to_q64_64(), 1);
    }

    #[test]
    fn inverse_at_sub_ps_boundaries() {
        // 61.44 MS/s: instant(1) = 16276 ps + floor(2^64/24) · 2^-64 ps.
        let m = map(0, 61_440_000, 1);
        let i1 = m.instant(SampleIndex::new(1)).unwrap();
        let just_before = AcqInstant::from_q64_64(i1.to_q64_64() - 1);
        let just_after = AcqInstant::from_q64_64(i1.to_q64_64() + 1);
        assert_eq!(m.index_at_or_after(i1).unwrap().get(), 1);
        assert_eq!(m.index_at_or_after(just_before).unwrap().get(), 1);
        assert_eq!(m.index_at_or_after(just_after).unwrap().get(), 2);
        assert_eq!(m.index_at_or_before(i1).unwrap().unwrap().get(), 1);
        assert_eq!(m.index_at_or_before(just_before).unwrap().unwrap().get(), 0);
        assert_eq!(m.index_at_or_before(just_after).unwrap().unwrap().get(), 1);
        let before_epoch = AcqInstant::new(TruthTime::from_ps(-1), u64::MAX);
        assert_eq!(m.index_at_or_before(before_epoch).unwrap(), None);
        assert_eq!(m.index_at_or_after(before_epoch).unwrap().get(), 0);
    }

    #[test]
    fn inverse_at_end_of_time() {
        let m = map(0, 250_000_000, 1);
        let end = AcqInstant::new(TruthTime::MAX, u64::MAX);
        let last = m.index_at_or_before(end).unwrap().unwrap();
        assert!(m.instant(last).is_ok());
        assert_eq!(
            m.instant(last.checked_add(1).unwrap()),
            Err(TimeError::TimeOverflow)
        );
        assert_eq!(m.index_at_or_after(end), Err(TimeError::TimeOverflow));
    }

    #[test]
    fn ids_are_stable_and_distinct() {
        let a = map(0, 61_440_000, 1);
        let b = map(0, 122_880_000, 2);
        let c = map(1, 61_440_000, 1);
        assert_eq!(a.id(), b.id());
        assert_ne!(a.id(), c.id());
        assert_ne!(a.id(), IdealTimestampClock.id());
        assert_eq!(a.id().to_string().len(), 64);
    }

    #[test]
    fn ideal_timestamp_is_identity() {
        let t = AcqInstant::new(TruthTime::from_ps(42), 9);
        let r = IdealTimestampClock.label(t).unwrap();
        assert_eq!((r.ps, r.frac), (42, 9));
    }
}
