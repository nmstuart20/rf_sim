//! Picosecond UTC labels for receiver-reported times, without floating point.

use iqsim_core::ReportedTime;

use crate::IoError;

const PS_PER_SECOND: i128 = 1_000_000_000_000;

/// Format `epoch + reported` as `YYYY-MM-DDTHH:MM:SS.ffffffffffffZ`.
///
/// The epoch is POSIX seconds plus nanoseconds. The calendar is proleptic
/// Gregorian, without leap seconds. Sub-picosecond fractions are floored toward
/// negative infinity, including for times before the epoch.
///
/// # Errors
///
/// Rejects nanoseconds at or above one billion and years outside `0000..=9999`.
pub fn format_datetime(
    epoch_unix_seconds: i64,
    epoch_nanos: u32,
    reported: ReportedTime,
) -> Result<String, IoError> {
    if epoch_nanos >= 1_000_000_000 {
        return Err(IoError::invalid("epoch_nanos must be below 1000000000"));
    }
    // The epoch is an integer number of ps. ReportedTime.ps is already its
    // floor, even when negative, so dropping the unsigned fraction is exact
    // flooring. i64 seconds times 10^12 comfortably fits in i128.
    let ps = i128::from(epoch_unix_seconds) * PS_PER_SECOND
        + i128::from(epoch_nanos) * 1_000
        + i128::from(reported.ps);
    let seconds = ps.div_euclid(PS_PER_SECOND);
    let fraction = ps.rem_euclid(PS_PER_SECOND);
    let days = seconds.div_euclid(86_400);
    let time_of_day = seconds.rem_euclid(86_400);
    let (year, month, day) = civil_from_days(days);
    if !(0..=9999).contains(&year) {
        return Err(IoError::invalid("UTC year must be in 0000..=9999"));
    }
    let hour = time_of_day / 3_600;
    let minute = time_of_day % 3_600 / 60;
    let second = time_of_day % 60;
    Ok(format!(
        "{year:04}-{month:02}-{day:02}T{hour:02}:{minute:02}:{second:02}.{fraction:012}Z"
    ))
}

// Howard Hinnant's civil_from_days, with Euclidean division for negative eras.
// Day zero is 1970-01-01. The intermediate year may be outside the output range.
fn civil_from_days(days: i128) -> (i128, i128, i128) {
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1_460 + doe / 36_524 - doe / 146_096) / 365;
    let mut year = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = mp + if mp < 10 { 3 } else { -9 };
    year += i128::from(month <= 2);
    (year, month, day)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn label(seconds: i64, nanos: u32, ps: i64, frac: u64) -> String {
        format_datetime(seconds, nanos, ReportedTime { ps, frac }).unwrap()
    }

    #[test]
    fn known_instants_and_negative_crossings() {
        assert_eq!(label(0, 0, 0, 0), "1970-01-01T00:00:00.000000000000Z");
        assert_eq!(
            label(1_767_225_600, 0, 0, 0),
            "2026-01-01T00:00:00.000000000000Z"
        );
        assert_eq!(
            label(1_709_164_800, 0, 0, 0),
            "2024-02-29T00:00:00.000000000000Z"
        );
        assert_eq!(label(0, 0, -1, 0), "1969-12-31T23:59:59.999999999999Z");
        assert_eq!(
            label(1_767_225_600, 0, -1, u64::MAX),
            "2025-12-31T23:59:59.999999999999Z"
        );
        assert_eq!(
            label(1_709_251_200, 0, -1_000_000_000_001, 0),
            "2024-02-29T23:59:58.999999999999Z"
        );
    }

    #[test]
    fn fractions_are_floored_and_nanos_carry() {
        assert_eq!(
            label(0, 123_456_789, 123, u64::MAX),
            "1970-01-01T00:00:00.123456789123Z"
        );
        assert_eq!(
            label(0, 999_999_999, 1_001, 1 << 63),
            "1970-01-01T00:00:01.000000000001Z"
        );
        assert_eq!(
            label(-1, 999_999_999, -1, 1 << 63),
            "1969-12-31T23:59:59.999999998999Z"
        );
    }

    // Independent calendar oracle: count leap years before January 1 and add
    // explicit month lengths, rather than invert the March-based era algorithm.
    fn days_from_civil(year: i64, month: usize, day: i64) -> i64 {
        fn before_year(y: i64) -> i64 {
            365 * y + (y + 3).div_euclid(4) - (y + 99).div_euclid(100) + (y + 399).div_euclid(400)
        }
        before_year(year) - before_year(1970)
            + month_lengths(year)[..month - 1].iter().sum::<i64>()
            + day
            - 1
    }

    fn month_lengths(year: i64) -> [i64; 12] {
        let leap = year % 4 == 0 && (year % 100 != 0 || year % 400 == 0);
        [
            31,
            if leap { 29 } else { 28 },
            31,
            30,
            31,
            30,
            31,
            31,
            30,
            31,
            30,
            31,
        ]
    }

    #[test]
    fn calendar_agrees_with_independent_days_from_civil() {
        for year in (0..=9999)
            .step_by(37)
            .chain([0, 1, 100, 400, 1900, 2000, 9999])
        {
            for (month_index, last_day) in month_lengths(year).into_iter().enumerate() {
                let month = month_index + 1;
                for day in [1, last_day] {
                    let seconds = days_from_civil(year, month, day) * 86_400 + 45_296;
                    assert_eq!(
                        label(seconds, 987_654_321, 123, u64::MAX),
                        format!("{year:04}-{month:02}-{day:02}T12:34:56.987654321123Z")
                    );
                }
            }
        }
    }

    #[test]
    fn range_errors_and_calendar_endpoints() {
        let first = days_from_civil(0, 1, 1) * 86_400;
        let after_last = days_from_civil(10_000, 1, 1) * 86_400;
        assert_eq!(label(first, 0, 0, 0), "0000-01-01T00:00:00.000000000000Z");
        assert_eq!(
            label(after_last, 0, -1, u64::MAX),
            "9999-12-31T23:59:59.999999999999Z"
        );
        for seconds in [i64::MIN, first - 1, after_last, i64::MAX] {
            assert!(format_datetime(seconds, 0, ReportedTime::default()).is_err());
        }
        for nanos in [1_000_000_000, u32::MAX] {
            assert!(format_datetime(0, nanos, ReportedTime::default()).is_err());
        }
        // Extreme reported values still do not overflow intermediate arithmetic.
        assert!(
            format_datetime(
                0,
                0,
                ReportedTime {
                    ps: i64::MIN,
                    frac: 0
                }
            )
            .is_ok()
        );
        assert!(
            format_datetime(
                0,
                0,
                ReportedTime {
                    ps: i64::MAX,
                    frac: u64::MAX
                }
            )
            .is_ok()
        );
    }
}
