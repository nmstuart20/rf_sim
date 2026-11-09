//! Angles (deg, rad, DMS) and geodetic latitude/longitude.

use std::borrow::Cow;
use std::fmt;

use serde::de::{self, Visitor};

use crate::decimal::{Decimal, number_len};
use crate::error::{QuantityKind, UnitError};
use crate::text::{TextQuantity, finite, number_and_unit, text_quantity_impls};

static ANGLE: QuantityKind = QuantityKind {
    name: "angle",
    units: &["deg", "°", "rad", "DMS (12°34'56.7\" or 12d34m56.7s)"],
    example: "12.5 deg",
};

/// A signed angle in degrees.
///
/// Written as `"12.5 deg"`, `"12.5°"`, `"0.2 rad"`, or degrees-minutes-seconds:
/// `12°34'56.7"` or the ASCII form `12d34m56.7s` (minutes and seconds optional,
/// each below 60; only the last component may have a fraction).
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct Angle {
    degrees: f64,
}

impl Angle {
    /// An angle of `degrees`. `None` unless finite.
    pub fn from_degrees(degrees: f64) -> Option<Self> {
        degrees.is_finite().then_some(Self { degrees })
    }

    /// The angle in degrees.
    pub fn degrees(self) -> f64 {
        self.degrees
    }

    /// The angle in radians.
    pub fn radians(self) -> f64 {
        self.degrees.to_radians()
    }
}

/// Whether a unit starts with a DMS marker (so a DMS error explains it best).
fn starts_dms(unit: &str) -> bool {
    unit.starts_with(['°', 'º', 'd', '\'', '′', '"', '″'])
}

impl TextQuantity for Angle {
    const KIND: &'static QuantityKind = &ANGLE;

    fn parse_text(text: &str) -> Result<Self, UnitError> {
        let kind = Self::KIND;
        let (value, unit) = number_and_unit(kind, text)?;
        let degrees = match unit {
            "deg" => finite(kind, text, value.to_f64())?,
            "rad" => finite(kind, text, value.to_f64())?.to_degrees(),
            _ if starts_dms(unit) => parse_dms(kind, text, text.trim())?,
            _ => {
                return Err(UnitError::UnknownUnit {
                    kind,
                    unit: unit.to_owned(),
                });
            }
        };
        Ok(Self {
            degrees: finite(kind, text, degrees)?,
        })
    }
}

impl fmt::Display for Angle {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} deg", self.degrees)
    }
}

text_quantity_impls!(Angle, None);

const MARKERS: [&[&str]; 3] = [&["°", "º", "d"], &["'", "′", "m"], &["\"", "″", "''", "s"]];

/// Parses `[+-]D° [M' [S"]]` (or `D d M m S s`) into degrees.
fn parse_dms(kind: &'static QuantityKind, text: &str, body: &str) -> Result<f64, UnitError> {
    let invalid = |reason: &str| UnitError::Invalid {
        kind,
        text: text.to_owned(),
        reason: reason.to_owned(),
    };
    let (negative, mut rest) = match body.strip_prefix('-') {
        Some(rest) => (true, rest),
        None => (false, body.strip_prefix('+').unwrap_or(body)),
    };
    let mut parts: Vec<(usize, Decimal)> = Vec::new();
    while !rest.trim_start().is_empty() {
        rest = rest.trim_start();
        if rest.starts_with(['+', '-']) {
            return Err(invalid("only the degrees may carry a sign"));
        }
        let n = number_len(rest);
        let value = Decimal::parse(&rest[..n])
            .ok_or_else(|| invalid("expected a number in degrees-minutes-seconds"))?;
        rest = rest[n..].trim_start();
        let previous = parts.last().map_or(-1, |&(level, _)| level as isize);
        let found = MARKERS.iter().enumerate().find_map(|(level, markers)| {
            markers
                .iter()
                .find(|m| rest.starts_with(**m))
                .map(|m| (level, m.len()))
        });
        let Some((level, len)) = found else {
            return Err(invalid(
                "each component needs a marker: degrees ° or d, minutes ' or m, seconds \" or s",
            ));
        };
        if level as isize <= previous || (parts.is_empty() && level != 0) {
            return Err(invalid(
                "components must be degrees, then minutes, then seconds",
            ));
        }
        rest = &rest[len..];
        parts.push((level, value));
    }
    let Some(&(last_level, _)) = parts.last() else {
        return Err(invalid("expected degrees"));
    };
    let factor = |level: usize| -> u128 {
        match level {
            0 => 3600,
            1 => 60,
            _ => 1,
        }
    };
    for &(level, value) in &parts {
        if level != last_level && !value.is_integer() {
            return Err(invalid("only the last component may have a fraction"));
        }
        if level > 0 && !value.less_than(60) {
            return Err(invalid("minutes and seconds must be below 60"));
        }
    }
    // Total arcseconds, exactly when the aligned sum fits in 38 digits (one rounding
    // to f64), otherwise as an f64 sum (e.g. `1d0m1e-100s`).
    let exact = parts
        .iter()
        .try_fold(Decimal::ZERO, |sum, &(level, value)| {
            sum.add_nonnegative(value.mul_u128(factor(level))?)
        });
    let arcseconds = match exact {
        Some(total) => total.to_f64(),
        None => parts
            .iter()
            .map(|&(level, value)| value.to_f64() * factor(level) as f64)
            .fold(0.0, |sum, term| sum + term),
    };
    let degrees = finite(kind, text, arcseconds)? / 3600.0;
    Ok(if negative { -degrees } else { degrees })
}

// ---------------------------------------------------------------- Latitude / Longitude

static LATITUDE: QuantityKind = QuantityKind {
    name: "latitude",
    units: &["decimal degrees", "DMS", "N/S suffix"],
    example: "38°53'23.3\"N",
};

static LONGITUDE: QuantityKind = QuantityKind {
    name: "longitude",
    units: &["decimal degrees", "DMS", "E/W suffix"],
    example: "77.0365 W",
};

fn parse_geodetic(
    kind: &'static QuantityKind,
    text: &str,
    positive: char,
    negative: char,
    limit: f64,
) -> Result<f64, UnitError> {
    let invalid = |reason: String| UnitError::Invalid {
        kind,
        text: text.to_owned(),
        reason,
    };
    let trimmed = text.trim();
    if trimmed.is_empty() {
        return Err(UnitError::Empty { kind });
    }
    // Upper case only: a lower-case `s` is the seconds marker in `38d30m36s`.
    let hemisphere = |c: char| (c == positive || c == negative).then_some(c);
    let (body, sign) = if let Some(h) = trimmed.chars().last().and_then(hemisphere) {
        (trimmed[..trimmed.len() - 1].trim_end(), Some(h))
    } else if let Some(h) = trimmed.chars().next().and_then(hemisphere) {
        (trimmed[1..].trim_start(), Some(h))
    } else {
        (trimmed, None)
    };
    if sign.is_some() && body.starts_with(['+', '-']) {
        return Err(invalid(format!(
            "use either a sign or a {positive}/{negative} hemisphere letter, not both"
        )));
    }
    let degrees = if let Some(value) = Decimal::parse(body) {
        finite(kind, text, value.to_f64())?
    } else {
        let (value, unit) = number_and_unit(kind, body)?;
        match unit {
            "deg" => finite(kind, text, value.to_f64())?,
            _ if starts_dms(unit) => parse_dms(kind, text, body)?,
            _ => {
                return Err(invalid(format!(
                    "expected decimal degrees or degrees-minutes-seconds with an optional {positive}/{negative}, e.g. \"{}\"",
                    kind.example
                )));
            }
        }
    };
    let degrees = if sign == Some(negative) {
        -degrees
    } else {
        degrees
    };
    if !(-limit..=limit).contains(&degrees) {
        return Err(UnitError::OutOfRange {
            kind,
            text: text.to_owned(),
            limit: format!("it must be within ±{limit}°"),
        });
    }
    Ok(degrees)
}

/// Geodetic latitude in degrees, within [-90, 90].
///
/// Written as a number (decimal degrees, north positive) or a string: `"38.8898"`,
/// `"38.8898 N"`, `"S33.85"` (hemisphere letters are upper case), `"38°53'23.3\"N"`, `"38d53m23.3s N"`, `"-33.85 deg"`.
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct Latitude(f64);

/// Geodetic longitude in degrees, within [-180, 180].
///
/// Written like [`Latitude`], east positive, with E/W hemisphere letters.
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct Longitude(f64);

macro_rules! geodetic {
    ($ty:ident, $kind:ident, $pos:literal, $neg:literal, $limit:literal) => {
        impl $ty {
            /// Checked constructor from degrees.
            pub fn from_degrees(degrees: f64) -> Option<Self> {
                (-$limit..=$limit).contains(&degrees).then_some(Self(degrees))
            }

            /// The value in degrees.
            pub fn degrees(self) -> f64 {
                self.0
            }

            fn from_number(number: f64) -> Result<Self, UnitError> {
                Self::from_degrees(number).ok_or_else(|| UnitError::OutOfRange {
                    kind: &$kind,
                    text: number.to_string(),
                    limit: format!("it must be within ±{}°", $limit),
                })
            }
        }

        impl ::std::str::FromStr for $ty {
            type Err = UnitError;
            fn from_str(text: &str) -> Result<Self, UnitError> {
                parse_geodetic(&$kind, text, $pos, $neg, $limit).map(Self)
            }
        }

        impl<'de> ::serde::Deserialize<'de> for $ty {
            fn deserialize<D: ::serde::Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
                struct V;
                impl Visitor<'_> for V {
                    type Value = $ty;
                    fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                        write!(f, "a {} in degrees, as a number or a string such as \"{}\"", $kind.name, $kind.example)
                    }
                    fn visit_str<E: de::Error>(self, v: &str) -> Result<$ty, E> {
                        v.parse().map_err(E::custom)
                    }
                    fn visit_f64<E: de::Error>(self, v: f64) -> Result<$ty, E> {
                        $ty::from_number(v).map_err(E::custom)
                    }
                    fn visit_i64<E: de::Error>(self, v: i64) -> Result<$ty, E> {
                        $ty::from_number(v as f64).map_err(E::custom)
                    }
                    fn visit_u64<E: de::Error>(self, v: u64) -> Result<$ty, E> {
                        $ty::from_number(v as f64).map_err(E::custom)
                    }
                }
                d.deserialize_any(V)
            }
        }

        impl ::schemars::JsonSchema for $ty {
            fn schema_name() -> Cow<'static, str> {
                stringify!($ty).into()
            }

            fn json_schema(_: &mut ::schemars::SchemaGenerator) -> ::schemars::Schema {
                ::schemars::json_schema!({
                    "description": format!(
                        "{} in degrees: a number ({} positive) or a string in decimal degrees or degrees-minutes-seconds with an optional {}/{} letter, e.g. \"{}\".",
                        $kind.name, $pos, $pos, $neg, $kind.example
                    ),
                    "anyOf": [
                        { "type": "number", "minimum": -$limit, "maximum": $limit },
                        { "type": "string" }
                    ]
                })
            }
        }
    };
}

geodetic!(Latitude, LATITUDE, 'N', 'S', 90.0);
geodetic!(Longitude, LONGITUDE, 'E', 'W', 180.0);

#[cfg(test)]
mod tests {
    use super::*;

    fn deg(s: &str) -> Result<f64, UnitError> {
        s.parse::<Angle>().map(Angle::degrees)
    }

    #[test]
    fn angle_forms() {
        assert_eq!(deg("12.5 deg"), Ok(12.5));
        assert_eq!(deg("12.5°"), Ok(12.5));
        assert_eq!(deg("-90deg"), Ok(-90.0));
        assert_eq!(deg("3.141592653589793 rad"), Ok(180.0));
        assert_eq!(deg("12°30'"), Ok(12.5));
        assert_eq!(deg("12°30'36\""), Ok(12.51));
        assert_eq!(deg("12d30m36s"), Ok(12.51));
        assert_eq!(deg("12d 30m 36s"), Ok(12.51));
        assert_eq!(deg("-12°30'"), Ok(-12.5));
        assert_eq!(deg("1°0'0.36\""), Ok(1.0001));
        assert_eq!(deg("12°59.5'"), Ok(12.0 + 59.5 / 60.0));
        assert_eq!(deg("0°0'59.999\""), Ok(59.999 / 3600.0));
        // Found by review: tiny fractions no longer overflow the exact sum.
        assert_eq!(deg("1d0m1e-100s"), Ok(1.0));
        for bad in [
            "12°60'",
            "12°30'60\"",
            "12.5°30'",
            "12°-30'",
            "12'",
            "12°30\"5'",
            "12°30",
            "12°°",
        ] {
            assert!(
                matches!(deg(bad), Err(UnitError::Invalid { .. })),
                "{bad}: {:?}",
                deg(bad)
            );
        }
        assert!(matches!(deg("12"), Err(UnitError::MissingUnit { .. })));
        assert!(matches!(deg("12 grad"), Err(UnitError::UnknownUnit { .. })));
    }

    #[test]
    fn latitude_and_longitude_forms() {
        let lat = |s: &str| s.parse::<Latitude>().map(Latitude::degrees);
        let lon = |s: &str| s.parse::<Longitude>().map(Longitude::degrees);
        assert_eq!(lat("38.5"), Ok(38.5));
        assert_eq!(lat("38.5 N"), Ok(38.5));
        assert_eq!(lat("38.5S"), Ok(-38.5));
        assert_eq!(lat("S 38.5"), Ok(-38.5));
        assert_eq!(lat("-38.5 deg"), Ok(-38.5));
        assert_eq!(lat("38°30'S"), Ok(-38.5));
        assert_eq!(lat("38d30m S"), Ok(-38.5));
        assert_eq!(lat("38d30m36s"), Ok(38.51));
        assert_eq!(lat("38d30m36sS"), Ok(-38.51));
        assert_eq!(lon("77°02'24\"W"), Ok(-77.04));
        assert_eq!(lon("180 E"), Ok(180.0));
        assert!(matches!(lat("90.1"), Err(UnitError::OutOfRange { .. })));
        assert!(matches!(lon("-180.5"), Err(UnitError::OutOfRange { .. })));
        assert!(matches!(lat("-38.5 S"), Err(UnitError::Invalid { .. })));
        assert!(matches!(lat("38.5 E"), Err(UnitError::Invalid { .. })));
        assert!(matches!(lon("12 rad"), Err(UnitError::Invalid { .. })));
        assert!(matches!(lat(""), Err(UnitError::Empty { .. })));
    }
}
