//! The scenario epoch: a UTC calendar instant.

use std::borrow::Cow;
use std::fmt;

use schemars::JsonSchema;
use serde::{Deserialize, Deserializer};
use toml::value::{Datetime, Offset};

/// The UTC instant of truth time zero (`TruthTime::EPOCH`).
///
/// Written as a TOML offset date-time in UTC, `epoch = 2026-01-01T00:00:00Z` (or
/// `+00:00`), or the same text as a string. Fractional seconds are kept to the
/// nanosecond. Leap seconds (`:60`) are rejected; the epoch is a label and the
/// simulation itself never steps across a leap second.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Epoch {
    /// Whole seconds since 1970-01-01T00:00:00Z, not counting leap seconds.
    unix_seconds: i64,
    /// Nanoseconds within the second.
    nanosecond: u32,
}

/// Days from 1970-01-01 to the proleptic Gregorian date (Hinnant's algorithm).
fn days_from_civil(year: i64, month: u32, day: u32) -> i64 {
    let y = if month <= 2 { year - 1 } else { year };
    let era = y.div_euclid(400);
    let yoe = y - era * 400;
    let mp = i64::from((month + 9) % 12);
    let doy = (153 * mp + 2) / 5 + i64::from(day) - 1;
    let doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    era * 146_097 + doe - 719_468
}

fn days_in_month(year: u16, month: u8) -> u8 {
    match month {
        1 | 3 | 5 | 7 | 8 | 10 | 12 => 31,
        4 | 6 | 9 | 11 => 30,
        _ if year.is_multiple_of(4) && (!year.is_multiple_of(100) || year.is_multiple_of(400)) => {
            29
        }
        _ => 28,
    }
}

impl Epoch {
    /// Converts a TOML date-time, which must be a full UTC date and time.
    pub fn from_datetime(datetime: &Datetime) -> Result<Self, String> {
        let (Some(date), Some(time)) = (datetime.date, datetime.time) else {
            return Err(format!(
                "the epoch \"{datetime}\" needs a date, a time and a UTC offset, e.g. 2026-01-01T00:00:00Z"
            ));
        };
        match datetime.offset {
            Some(Offset::Z) | Some(Offset::Custom { minutes: 0 }) => {}
            Some(Offset::Custom { .. }) => {
                return Err(format!(
                    "the epoch \"{datetime}\" must be in UTC: write it with a Z offset, e.g. 2026-01-01T00:00:00Z"
                ));
            }
            None => {
                return Err(format!(
                    "the epoch \"{datetime}\" has no UTC offset; add Z, e.g. 2026-01-01T00:00:00Z"
                ));
            }
        }
        if !(1..=12).contains(&date.month)
            || date.day == 0
            || date.day > days_in_month(date.year, date.month)
        {
            return Err(format!("the epoch date in \"{datetime}\" does not exist"));
        }
        let second = time.second.unwrap_or(0);
        if second == 60 {
            return Err(format!(
                "the epoch \"{datetime}\" is a leap second, which is not supported as an epoch"
            ));
        }
        if time.hour > 23 || time.minute > 59 || second > 59 {
            return Err(format!("the epoch time in \"{datetime}\" is out of range"));
        }
        let days = days_from_civil(
            i64::from(date.year),
            u32::from(date.month),
            u32::from(date.day),
        );
        let seconds = i64::from(time.hour) * 3600 + i64::from(time.minute) * 60 + i64::from(second);
        Ok(Self {
            unix_seconds: days * 86_400 + seconds,
            nanosecond: time.nanosecond.unwrap_or(0),
        })
    }

    /// Whole seconds since 1970-01-01T00:00:00Z (POSIX time, no leap seconds).
    pub fn unix_seconds(self) -> i64 {
        self.unix_seconds
    }

    /// Nanoseconds within the second.
    pub fn nanosecond(self) -> u32 {
        self.nanosecond
    }
}

impl fmt::Display for Epoch {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let days = self.unix_seconds.div_euclid(86_400);
        let secs = self.unix_seconds.rem_euclid(86_400);
        // Inverse of days_from_civil (Hinnant).
        let z = days + 719_468;
        let era = z.div_euclid(146_097);
        let doe = z - era * 146_097;
        let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
        let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
        let mp = (5 * doy + 2) / 153;
        let day = doy - (153 * mp + 2) / 5 + 1;
        let month = if mp < 10 { mp + 3 } else { mp - 9 };
        let year = yoe + era * 400 + i64::from(month <= 2);
        write!(
            f,
            "{year:04}-{month:02}-{day:02}T{:02}:{:02}:{:02}",
            secs / 3600,
            secs / 60 % 60,
            secs % 60
        )?;
        if self.nanosecond != 0 {
            let frac = format!("{:09}", self.nanosecond);
            write!(f, ".{}", frac.trim_end_matches('0'))?;
        }
        f.write_str("Z")
    }
}

impl<'de> Deserialize<'de> for Epoch {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        use serde::de::Error;
        let datetime = match toml::Value::deserialize(d)? {
            toml::Value::Datetime(datetime) => datetime,
            toml::Value::String(text) => text.parse::<Datetime>().map_err(|_| {
                D::Error::custom(format!(
                    "\"{text}\" is not a date-time; write the epoch as 2026-01-01T00:00:00Z"
                ))
            })?,
            other => {
                return Err(D::Error::custom(format!(
                    "expected the epoch as a UTC date-time such as 2026-01-01T00:00:00Z, got a {}",
                    other.type_str()
                )));
            }
        };
        Self::from_datetime(&datetime).map_err(D::Error::custom)
    }
}

impl JsonSchema for Epoch {
    fn schema_name() -> Cow<'static, str> {
        "Epoch".into()
    }

    fn json_schema(_: &mut schemars::SchemaGenerator) -> schemars::Schema {
        schemars::json_schema!({
            "type": "string",
            "format": "date-time",
            "description": "UTC instant of truth time zero, as a TOML date-time with a Z (or +00:00) offset, e.g. 2026-01-01T00:00:00Z. Leap seconds are not accepted.",
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn epoch(s: &str) -> Result<Epoch, String> {
        Epoch::from_datetime(&s.parse::<Datetime>().unwrap())
    }

    #[test]
    fn converts_utc_instants() {
        assert_eq!(epoch("1970-01-01T00:00:00Z").unwrap().unix_seconds(), 0);
        let e = epoch("2026-01-01T00:00:00Z").unwrap();
        assert_eq!(e.unix_seconds(), 1_767_225_600);
        assert_eq!(e.to_string(), "2026-01-01T00:00:00Z");
        let e = epoch("2000-02-29T23:59:59.250+00:00").unwrap();
        assert_eq!(
            (e.unix_seconds(), e.nanosecond()),
            (951_868_799, 250_000_000)
        );
        assert_eq!(e.to_string(), "2000-02-29T23:59:59.25Z");
        assert_eq!(epoch("1969-12-31T23:59:59Z").unwrap().unix_seconds(), -1);
        assert_eq!(
            epoch("1900-03-01T00:00:00Z").unwrap().to_string(),
            "1900-03-01T00:00:00Z"
        );
    }

    #[test]
    fn rejects_non_utc_and_incomplete_values() {
        assert!(
            epoch("2026-01-01T00:00:00+01:00")
                .unwrap_err()
                .contains("must be in UTC")
        );
        assert!(
            epoch("2026-01-01T00:00:00")
                .unwrap_err()
                .contains("no UTC offset")
        );
        assert!(
            epoch("2026-01-01")
                .unwrap_err()
                .contains("needs a date, a time")
        );
        assert!(
            epoch("2026-12-31T23:59:60Z")
                .unwrap_err()
                .contains("leap second")
        );
        // toml's parser rejects impossible dates already; the check here is a backstop
        // for values built in code.
        let mut feb29 = "2024-02-29T00:00:00Z".parse::<Datetime>().unwrap();
        assert!(Epoch::from_datetime(&feb29).is_ok());
        feb29.date.as_mut().unwrap().year = 2025;
        assert!(
            Epoch::from_datetime(&feb29)
                .unwrap_err()
                .contains("does not exist")
        );
    }

    #[test]
    fn civil_round_trip_over_four_centuries() {
        for days in (-146_097 * 2..146_097 * 2).step_by(97) {
            let e = Epoch {
                unix_seconds: days * 86_400 + 3723,
                nanosecond: 0,
            };
            assert_eq!(epoch(&e.to_string()), Ok(e));
        }
    }
}
