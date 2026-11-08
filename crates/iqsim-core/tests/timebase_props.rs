//! Property and edge-case tests for truth time, sample rates, the ideal clock map and
//! block tiling (prompt 1).

use iqsim_core::{
    AcqInstant, Block, BlockSpec, ChannelTimebase, ClockMap, IdealClockMap, MAX_RATE_DEN,
    MAX_RATE_HZ, MAX_RATE_NUM, PS_PER_SECOND, SampleIndex, SampleRate, TimeError, TruthTime,
};
use proptest::prelude::*;

/// Cases per property when `PROPTEST_CASES` is not set.
const DEFAULT_CASES: u32 = 2048;

/// Proptest settings: `PROPTEST_CASES` overrides the case count, and every other
/// setting (including `PROPTEST_RNG_SEED`) comes from proptest's environment defaults.
/// Failing cases are saved under `proptest-regressions/`; commit those files.
fn config() -> ProptestConfig {
    let mut config = ProptestConfig::default();
    if std::env::var_os("PROPTEST_CASES").is_none() {
        config.cases = DEFAULT_CASES;
    }
    config
}

const HUNDRED_DAYS_PS: i64 = 100 * 86_400 * PS_PER_SECOND;

/// Any accepted rate, including large numerators and denominators.
fn rate() -> impl Strategy<Value = SampleRate> {
    (1..=MAX_RATE_DEN)
        .prop_flat_map(|den| {
            let max_num = u64::try_from(
                (u128::from(MAX_RATE_HZ) * u128::from(den)).min(u128::from(MAX_RATE_NUM)),
            )
            .unwrap();
            (1..=max_num, Just(den))
        })
        .prop_map(|(num, den)| SampleRate::from_ratio(num, den).unwrap())
}

/// Rates that sample-rate users actually configure, plus the two named edge rates.
fn practical_rate() -> impl Strategy<Value = SampleRate> {
    prop_oneof![
        Just(SampleRate::from_hz(61_440_000).unwrap()),
        Just(SampleRate::from_hz(250_000_000).unwrap()),
        Just(SampleRate::from_hz(MAX_RATE_HZ).unwrap()),
        Just(SampleRate::from_ratio(1_000_000_000, 3).unwrap()),
        (1..=MAX_RATE_HZ).prop_map(|hz| SampleRate::from_hz(hz).unwrap()),
        rate(),
    ]
}

/// Epochs that leave 100 days of positive range after them.
fn epoch() -> impl Strategy<Value = TruthTime> {
    (i64::MIN..=i64::MAX - HUNDRED_DAYS_PS).prop_map(TruthTime::from_ps)
}

/// Last sample index within 100 days of the epoch at `r`.
fn hundred_day_limit(r: SampleRate) -> u64 {
    let n = u128::try_from(HUNDRED_DAYS_PS).unwrap() * u128::from(r.num())
        / (u128::from(r.den()) * 1_000_000_000_000);
    u64::try_from(n).unwrap()
}

fn timebase_and_index() -> impl Strategy<Value = (ChannelTimebase, SampleIndex)> {
    (practical_rate(), epoch()).prop_flat_map(|(r, e)| {
        let tb = ChannelTimebase::new(e, r);
        (
            Just(tb),
            (0..=hundred_day_limit(r)).prop_map(SampleIndex::new),
        )
    })
}

/// Exact `n · den · 10¹²` for checking, or `None` if it exceeds u128.
fn exact_product(r: SampleRate, n: SampleIndex) -> Option<u128> {
    (u128::from(n.get()) * u128::from(r.den())).checked_mul(1_000_000_000_000)
}

proptest! {
    #![proptest_config(config())]

    /// index → time is exact: floor and remainder reconstruct n · den · 10¹² / num.
    #[test]
    fn index_to_time_is_exact((tb, n) in timebase_and_index()) {
        let r = tb.rate();
        let t = tb.index_to_time(n).unwrap();
        let span = u128::try_from(
            i128::from(t.floor().as_ps()) - i128::from(tb.epoch().as_ps()),
        ).unwrap();
        prop_assert!(t.rem() < t.rem_den());
        prop_assert_eq!(t.rem_den(), r.num());
        prop_assert_eq!(
            span * u128::from(r.num()) + u128::from(t.rem()),
            exact_product(r, n).unwrap()
        );
    }

    /// Round trips through the floor/ceiling inverses of the exact map.
    #[test]
    fn timebase_round_trip((tb, n) in timebase_and_index()) {
        let t = tb.index_to_time(n).unwrap();
        prop_assert_eq!(tb.index_at_or_after(t.floor()).unwrap(), n);
        let floor_index = tb.index_at_or_before(t.floor()).unwrap();
        if t.is_exact() {
            prop_assert_eq!(floor_index, n);
        } else {
            // Samples are >= 1000 ps apart, so the previous sample is <= floor.
            prop_assert_eq!(floor_index.checked_add(1).unwrap(), n);
        }
    }

    /// The clock map's inverses return exactly the sample whose instant was asked for,
    /// and the sub-ps fraction is the floor of the exact remainder.
    #[test]
    fn clock_map_round_trip((tb, n) in timebase_and_index()) {
        let m = IdealClockMap::new(tb);
        let i = m.instant(n).unwrap();
        let exact = tb.index_to_time(n).unwrap();
        prop_assert_eq!(i.ps(), exact.floor());
        prop_assert_eq!(
            u128::from(i.frac()),
            (u128::from(exact.rem()) << 64) / u128::from(exact.rem_den())
        );
        prop_assert_eq!(m.index_at_or_after(i).unwrap(), n);
        prop_assert_eq!(m.index_at_or_before(i).unwrap(), Some(n));
    }

    /// For any instant t, the inverses bracket t correctly.
    #[test]
    fn clock_map_inverse_brackets(
        tb in (practical_rate(), epoch()).prop_map(|(r, e)| ChannelTimebase::new(e, r)),
        offset in -1_000_000_000_i64..=HUNDRED_DAYS_PS,
        frac in any::<u64>(),
    ) {
        let m = IdealClockMap::new(tb);
        let t = AcqInstant::new(
            TruthTime::from_ps(tb.epoch().as_ps().saturating_add(offset)),
            frac,
        );
        // At very low rates the sample after t may not be representable; then
        // index_at_or_after must report TimeOverflow rather than an index.
        let after = match m.index_at_or_after(t) {
            Ok(after) => {
                prop_assert!(m.instant(after).unwrap() >= t);
                if let Ok(prev) = after.checked_sub(1) {
                    prop_assert!(m.instant(prev).unwrap() < t);
                }
                Some(after)
            }
            Err(err) => {
                prop_assert_eq!(err, TimeError::TimeOverflow);
                None
            }
        };
        match m.index_at_or_before(t).unwrap() {
            Some(before) => {
                prop_assert!(m.instant(before).unwrap() <= t);
                let next = before.checked_add(1).unwrap();
                match m.instant(next) {
                    Ok(i) => {
                        prop_assert!(i > t);
                        let a = after.unwrap();
                        prop_assert!(a == before || a == next);
                    }
                    Err(err) => {
                        prop_assert_eq!(err, TimeError::TimeOverflow);
                        prop_assert!(after.is_none() || after == Some(before));
                    }
                }
            }
            None => {
                prop_assert!(t < m.instant(SampleIndex::ZERO).unwrap());
                prop_assert_eq!(after, Some(SampleIndex::ZERO));
            }
        }
    }

    /// Instants are strictly increasing; the time → index maps are non-decreasing.
    #[test]
    fn monotonicity(
        (tb, a) in timebase_and_index(),
        delta in 1_u64..=1_000_000,
        dt in 0_i64..=1_000_000_000,
    ) {
        let m = IdealClockMap::new(tb);
        let b = a.checked_add(delta).unwrap();
        if let (Ok(ia), Ok(ib)) = (m.instant(a), m.instant(b)) {
            prop_assert!(ia < ib);
        }
        let t0 = tb.index_to_time(a).unwrap().floor();
        let Ok(t1) = t0.checked_add(iqsim_core::TruthDuration::from_ps(dt)) else {
            return Ok(());
        };
        prop_assert!(tb.index_at_or_before(t0).unwrap() <= tb.index_at_or_before(t1).unwrap());
        prop_assert!(tb.index_at_or_after(t0).unwrap() <= tb.index_at_or_after(t1).unwrap());
    }

    /// Any u64 index either converts exactly or overflows only when the true time does
    /// not fit in i64 picoseconds. Never panics.
    #[test]
    fn overflow_is_never_spurious(r in rate(), e in any::<i64>(), n in any::<u64>()) {
        let tb = ChannelTimebase::new(TruthTime::from_ps(e), r);
        let n = SampleIndex::new(n);
        // Independent oracle: every step is checked, so it cannot itself overflow.
        let fits = exact_product(r, n)
            .and_then(|p| i128::try_from(p / u128::from(r.num())).ok())
            .and_then(|q| q.checked_add(i128::from(e)))
            .is_some_and(|t| i64::try_from(t).is_ok());
        match tb.index_to_time(n) {
            Ok(_) => prop_assert!(fits),
            Err(err) => {
                prop_assert_eq!(err, TimeError::TimeOverflow);
                prop_assert!(!fits);
            }
        }
    }

    /// Block tiling: pieces are non-empty, in order, contiguous, inside one grid block,
    /// and cover the span exactly.
    #[test]
    fn blocks_tile_without_gaps_or_overlaps(
        len in 1_u64..=5_000_000,
        start in any::<u64>(),
        count in 0_u64..=20_000_000,
    ) {
        let spec = BlockSpec::new(len).unwrap();
        let end = start.saturating_add(count);
        let span = Block::new(SampleIndex::new(start), SampleIndex::new(end)).unwrap();
        let mut cursor = start;
        let mut last_number = None;
        for (number, piece) in spec.tile(span) {
            prop_assert!(!piece.is_empty());
            prop_assert!(piece.len() <= len);
            prop_assert_eq!(piece.start().get(), cursor);
            prop_assert_eq!(spec.block_number_of(piece.start()), number);
            prop_assert_eq!(spec.block_number_of(SampleIndex::new(piece.end().get() - 1)), number);
            if let Some(prev) = last_number {
                prop_assert_eq!(number, prev + 1);
            }
            if piece.start().get() != start && piece.end().get() != end {
                prop_assert_eq!(spec.block(number).unwrap(), piece);
            }
            last_number = Some(number);
            cursor = piece.end().get();
        }
        prop_assert_eq!(cursor, end);
    }
}

#[test]
fn rate_bounds_cover_every_accepted_rate() {
    // The corner rates of the accepted region all convert at 100 days without overflow.
    for (num, den) in [
        (MAX_RATE_HZ, 1),
        (MAX_RATE_NUM, MAX_RATE_NUM / MAX_RATE_HZ + 1),
        (MAX_RATE_NUM - 1, MAX_RATE_DEN - 1),
        (1, MAX_RATE_DEN),
        (1_000_000_000, 3),
    ] {
        let r = SampleRate::from_ratio(num, den).unwrap();
        let tb = ChannelTimebase::new(TruthTime::EPOCH, r);
        let n = SampleIndex::new(hundred_day_limit(r));
        let t = tb.index_to_time(n).unwrap();
        assert!(t.floor().as_ps() <= HUNDRED_DAYS_PS);
        assert_eq!(tb.index_at_or_after(t.floor()).unwrap(), n);
    }
}

#[test]
fn edge_cases_61_44_msps() {
    let r = SampleRate::from_hz(61_440_000).unwrap();
    let tb = ChannelTimebase::new(TruthTime::EPOCH, r);
    let m = IdealClockMap::new(tb);
    // Period is 16276 + 1/24 ps, so every 24th sample is exact.
    let t24 = tb.index_to_time(SampleIndex::new(24)).unwrap();
    assert!(t24.is_exact());
    assert_eq!(t24.floor().as_ps(), 390_625);
    let t23 = tb.index_to_time(SampleIndex::new(23)).unwrap();
    assert_eq!(
        (t23.floor().as_ps(), t23.rem() * 24 / t23.rem_den()),
        (374_348, 23)
    );
    // 1 ms is exactly 61 440 samples.
    let ms = tb.index_to_time(SampleIndex::new(61_440)).unwrap();
    assert!(ms.is_exact());
    assert_eq!(ms.floor().as_ps(), 1_000_000_000);
    // 100 days.
    let n100 = SampleIndex::new(61_440_000 * 86_400 * 100);
    let t100 = tb.index_to_time(n100).unwrap();
    assert!(t100.is_exact());
    assert_eq!(t100.floor().as_ps(), HUNDRED_DAYS_PS);
    assert_eq!(m.index_at_or_after(m.instant(n100).unwrap()).unwrap(), n100);
    // One sample after 100 days has the 1/24 ps remainder.
    let t_next = m.instant(n100.checked_add(1).unwrap()).unwrap();
    assert_eq!(t_next.ps().as_ps(), HUNDRED_DAYS_PS + 16_276);
    assert_eq!(t_next.frac(), u64::MAX / 24);
    // The last representable sample, and the first that overflows.
    let last = tb.index_at_or_before(TruthTime::MAX).unwrap();
    assert!(tb.index_to_time(last).is_ok());
    assert_eq!(
        tb.index_to_time(last.checked_add(1).unwrap()),
        Err(TimeError::TimeOverflow)
    );
}

#[test]
fn edge_cases_250_msps() {
    let r = SampleRate::from_hz(250_000_000).unwrap();
    let tb = ChannelTimebase::new(TruthTime::from_ps(-123), r);
    let m = IdealClockMap::new(tb);
    // Every sample is an exact 4000 ps step.
    for n in [0_u64, 1, 249_999, 250_000, 250_000 * 86_400 * 100] {
        let t = tb.index_to_time(SampleIndex::new(n)).unwrap();
        assert!(t.is_exact());
        assert_eq!(t.floor().as_ps(), -123 + 4_000 * i64::try_from(n).unwrap());
    }
    // Between samples, floor and ceiling pick neighbours.
    let mid = TruthTime::from_ps(-123 + 4_000 * 7 + 1);
    assert_eq!(tb.index_at_or_before(mid).unwrap().get(), 7);
    assert_eq!(tb.index_at_or_after(mid).unwrap().get(), 8);
    let at = AcqInstant::from_truth_time(TruthTime::from_ps(-123 + 4_000 * 7));
    assert_eq!(m.index_at_or_before(at).unwrap().unwrap().get(), 7);
    assert_eq!(m.index_at_or_after(at).unwrap().get(), 7);
    // 1 ms blocks are 250 000 samples and line up with exact millisecond times.
    let spec = BlockSpec::new(250_000).unwrap();
    let b = spec.block(10).unwrap();
    assert_eq!(
        tb.index_to_time(b.start()).unwrap().floor().as_ps(),
        -123 + 10 * 1_000_000_000
    );
}

/// Regression cases from review: the slowest accepted rate, where a few samples
/// already exceed the i64 time range.
#[test]
fn slowest_rate_overflow_cases() {
    let r = SampleRate::from_ratio(1, MAX_RATE_DEN).unwrap();
    let tb = ChannelTimebase::new(TruthTime::from_ps(1), r);
    assert_eq!(
        tb.index_to_time(SampleIndex::new(1 << 56)),
        Err(TimeError::TimeOverflow)
    );
    let m = IdealClockMap::new(ChannelTimebase::new(TruthTime::EPOCH, r));
    // Sample 1 is 2^32 s ≈ 4.3e24 ps away, beyond the i64 range.
    let t = AcqInstant::from_truth_time(TruthTime::from_ps(1));
    assert_eq!(m.index_at_or_after(t), Err(TimeError::TimeOverflow));
    assert_eq!(m.index_at_or_before(t), Ok(Some(SampleIndex::ZERO)));
    assert_eq!(m.instant(SampleIndex::new(1)), Err(TimeError::TimeOverflow));
}
