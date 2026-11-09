//! Altitude with an explicit vertical reference, and geodetic positions.

use std::borrow::Cow;
use std::fmt;

use schemars::JsonSchema;
use serde::Deserialize;

use crate::angle::{Latitude, Longitude};
use crate::decimal::Decimal;
use crate::error::{QuantityKind, UnitError};
use crate::float::length_meters;
use crate::mgrs::Mgrs;
use crate::text::{NUMBER_PATTERN, TextQuantity, finite, number_and_unit, text_quantity_impls};

// ---------------------------------------------------------------- Altitude

static ALTITUDE: QuantityKind = QuantityKind {
    name: "altitude",
    units: &["m", "km", "ft", "nmi", "FL"],
    example: "35000 ft MSL",
};

/// The vertical datum an [`Altitude`] is measured from.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum AltitudeRef {
    /// Height above the WGS84 ellipsoid (`ellipsoid` or `HAE`).
    Ellipsoid,
    /// Height above mean sea level, i.e. the geoid (`MSL`).
    Msl,
    /// Height above the terrain below (`AGL`).
    Agl,
    /// Pressure altitude in the ISA standard atmosphere (1013.25 hPa datum), as
    /// written with a flight level (`FL350`). It is not a geometric height.
    PressureStd,
}

impl AltitudeRef {
    /// The keyword used in scenario text.
    pub const fn keyword(self) -> &'static str {
        match self {
            Self::Ellipsoid => "ellipsoid",
            Self::Msl => "MSL",
            Self::Agl => "AGL",
            Self::PressureStd => "FL",
        }
    }
}

/// An altitude in metres together with its reference.
///
/// Written as a number, a length unit (m, km, ft, nmi) and a required reference:
/// `"35000 ft MSL"`, `"120 m AGL"`, `"10 m ellipsoid"` (or `HAE`). A flight level
/// `"FL350"` is a standard pressure altitude of 35 000 ft and takes no reference.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Altitude {
    meters: f64,
    reference: AltitudeRef,
}

impl Altitude {
    /// An altitude of `meters` above `reference`. `None` unless finite.
    pub fn new(meters: f64, reference: AltitudeRef) -> Option<Self> {
        meters.is_finite().then_some(Self { meters, reference })
    }

    /// The altitude in metres (pressure altitude for [`AltitudeRef::PressureStd`]).
    pub fn meters(self) -> f64 {
        self.meters
    }

    /// The vertical reference.
    pub fn reference(self) -> AltitudeRef {
        self.reference
    }
}

impl TextQuantity for Altitude {
    const KIND: &'static QuantityKind = &ALTITUDE;

    fn parse_text(text: &str) -> Result<Self, UnitError> {
        let kind = Self::KIND;
        let invalid = |reason: &str| UnitError::Invalid {
            kind,
            text: text.to_owned(),
            reason: reason.to_owned(),
        };
        let trimmed = text.trim();
        if let Some(level) = trimmed.strip_prefix("FL") {
            let level = level.trim_start();
            if level.is_empty() || level.len() > 3 || !level.bytes().all(|b| b.is_ascii_digit()) {
                return Err(invalid(
                    "a flight level is FL and up to three digits of hundreds of feet, e.g. FL350, with no unit or reference",
                ));
            }
            let hundreds: u128 = level.parse().expect("ASCII digits");
            // FL n = n · 100 ft = n · 30.48 m, exact before the one rounding.
            let meters = Decimal::from_u128(hundreds)
                .mul_u128(3048)
                .and_then(|d| d.shifted(-2))
                .expect("at most 999 · 3048")
                .to_f64();
            return Ok(Self {
                meters,
                reference: AltitudeRef::PressureStd,
            });
        }
        let (value, rest) = number_and_unit(kind, text)?;
        let mut words = rest.split_whitespace();
        let unit = words
            .next()
            .expect("number_and_unit leaves a non-empty unit");
        let meters =
            length_meters(kind, text, value, unit).ok_or_else(|| UnitError::UnknownUnit {
                kind,
                unit: unit.to_owned(),
            })??;
        let reference = match (words.next(), words.next()) {
            (Some("MSL"), None) => AltitudeRef::Msl,
            (Some("AGL"), None) => AltitudeRef::Agl,
            (Some("ellipsoid" | "HAE"), None) => AltitudeRef::Ellipsoid,
            (None, _) => {
                return Err(invalid(
                    "an altitude needs a reference after the unit: MSL, AGL or ellipsoid (HAE), e.g. \"35000 ft MSL\"",
                ));
            }
            _ => {
                return Err(invalid(
                    "the reference must be one of MSL, AGL, ellipsoid or HAE (case-sensitive)",
                ));
            }
        };
        Ok(Self {
            meters: finite(kind, text, meters)?,
            reference,
        })
    }
}

/// Parseable for every value parsed from text: `FLnnn` for a flight level,
/// `"<m> m <reference>"` otherwise. A pressure altitude that is not a whole flight
/// level (only constructible in code) prints as `"<m> m pressure altitude"`.
impl fmt::Display for Altitude {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.reference == AltitudeRef::PressureStd {
            let level = libm::round(self.meters / 30.48);
            if (0.0..=999.0).contains(&level)
                && format!("FL{level:03}").parse::<Self>().ok() == Some(*self)
            {
                return write!(f, "FL{level:03}");
            }
            return write!(f, "{} m pressure altitude", self.meters);
        }
        write!(f, "{} m {}", self.meters, self.reference.keyword())
    }
}

text_quantity_impls!(
    Altitude,
    Some(format!(
        r"^\s*(?:FL\s*[0-9]{{1,3}}|{NUMBER_PATTERN}\s*(?:m|km|ft|nmi)\s+(?:MSL|AGL|HAE|ellipsoid))\s*$"
    ))
);

// ---------------------------------------------------------------- MGRS serde

impl<'de> Deserialize<'de> for Mgrs {
    fn deserialize<D: serde::Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let text = String::deserialize(d)?;
        Mgrs::parse(&text).map_err(serde::de::Error::custom)
    }
}

impl JsonSchema for Mgrs {
    fn schema_name() -> Cow<'static, str> {
        "Mgrs".into()
    }

    fn json_schema(_: &mut schemars::SchemaGenerator) -> schemars::Schema {
        schemars::json_schema!({
            "type": "string",
            "description": "An MGRS reference on WGS84 (zone, band, 100 km square, 0–5 digits each of easting and northing), e.g. \"18SUJ2337106519\". It denotes the south-west corner of the grid square. Polar (UPS) bands A, B, Y, Z are not supported.",
            "pattern": r"^\s*[0-9]{1,2}\s*[C-HJ-NP-Xc-hj-np-x]\s*[A-HJ-NP-Za-hj-np-z]{2}\s*(?:[0-9]\s*){0,10}$",
            "examples": ["18SUJ2337106519"],
        })
    }
}

// ---------------------------------------------------------------- Position

/// A position as written in a scenario: `lat` and `lon`, or `mgrs`, plus `alt`.
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct RawPosition {
    /// Geodetic latitude (WGS84).
    lat: Option<Latitude>,
    /// Geodetic longitude (WGS84).
    lon: Option<Longitude>,
    /// MGRS reference (instead of lat/lon); the south-west corner of its square.
    mgrs: Option<Mgrs>,
    /// Altitude with its reference, e.g. "120 m AGL".
    alt: Altitude,
}

/// A geodetic WGS84 position with an altitude.
///
/// In TOML: `{ lat = "38°53'23.3\"N", lon = -77.0365, alt = "100 m MSL" }` or
/// `{ mgrs = "18SUJ2337106519", alt = "FL350" }`. An MGRS position is the south-west
/// corner of the grid square (NGA convention); the reference is kept for display.
#[derive(Clone, Debug, PartialEq, Deserialize)]
#[serde(try_from = "RawPosition")]
pub struct Position {
    latitude: Latitude,
    longitude: Longitude,
    altitude: Altitude,
    mgrs: Option<Mgrs>,
}

impl Position {
    /// A position from latitude, longitude and altitude.
    pub fn new(latitude: Latitude, longitude: Longitude, altitude: Altitude) -> Self {
        Self {
            latitude,
            longitude,
            altitude,
            mgrs: None,
        }
    }

    /// A position at the south-west corner of an MGRS square.
    pub fn from_mgrs(mgrs: Mgrs, altitude: Altitude) -> Self {
        let (lat, lon) = mgrs.sw_corner();
        Self {
            latitude: Latitude::from_degrees(lat).expect("UTM inverse gives a valid latitude"),
            longitude: Longitude::from_degrees(lon).expect("UTM inverse gives [-180, 180)"),
            altitude,
            mgrs: Some(mgrs),
        }
    }

    /// Geodetic latitude.
    pub fn latitude(&self) -> Latitude {
        self.latitude
    }

    /// Geodetic longitude.
    pub fn longitude(&self) -> Longitude {
        self.longitude
    }

    /// Altitude and its reference.
    pub fn altitude(&self) -> Altitude {
        self.altitude
    }

    /// The MGRS reference this position was written as, if any.
    pub fn mgrs(&self) -> Option<&Mgrs> {
        self.mgrs.as_ref()
    }
}

/// Why a position table is incomplete or ambiguous.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PositionError(&'static str);

impl fmt::Display for PositionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.0)
    }
}

impl std::error::Error for PositionError {}

impl TryFrom<RawPosition> for Position {
    type Error = PositionError;

    fn try_from(raw: RawPosition) -> Result<Self, PositionError> {
        match (raw.lat, raw.lon, raw.mgrs) {
            (Some(lat), Some(lon), None) => Ok(Self::new(lat, lon, raw.alt)),
            (None, None, Some(mgrs)) => Ok(Self::from_mgrs(mgrs, raw.alt)),
            (_, _, Some(_)) => Err(PositionError(
                "a position takes either `lat` and `lon`, or `mgrs`, not both",
            )),
            (Some(_), None, None) | (None, Some(_), None) => Err(PositionError(
                "a position needs both `lat` and `lon` (or `mgrs` instead)",
            )),
            (None, None, None) => Err(PositionError(
                "a position needs `lat` and `lon`, or `mgrs`, plus `alt`",
            )),
        }
    }
}

impl JsonSchema for Position {
    fn schema_name() -> Cow<'static, str> {
        "Position".into()
    }

    fn json_schema(generator: &mut schemars::SchemaGenerator) -> schemars::Schema {
        let mut schema = RawPosition::json_schema(generator);
        schema.insert(
            "description".to_owned(),
            "A WGS84 position: `lat` and `lon`, or `mgrs`, plus `alt`.".into(),
        );
        schema.insert(
            "oneOf".to_owned(),
            serde_json::json!([
                { "required": ["lat", "lon"], "not": { "required": ["mgrs"] } },
                {
                    "required": ["mgrs"],
                    "not": { "anyOf": [{ "required": ["lat"] }, { "required": ["lon"] }] }
                }
            ]),
        );
        schema
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn alt(s: &str) -> Result<(f64, AltitudeRef), UnitError> {
        s.parse::<Altitude>().map(|a| (a.meters(), a.reference()))
    }

    #[test]
    fn altitude_forms() {
        assert_eq!(alt("35000 ft MSL"), Ok((10_668.0, AltitudeRef::Msl)));
        assert_eq!(alt("120 m AGL"), Ok((120.0, AltitudeRef::Agl)));
        assert_eq!(alt("-20 m ellipsoid"), Ok((-20.0, AltitudeRef::Ellipsoid)));
        assert_eq!(alt("10m HAE"), Ok((10.0, AltitudeRef::Ellipsoid)));
        assert_eq!(alt("FL350"), Ok((10_668.0, AltitudeRef::PressureStd)));
        assert_eq!(alt("FL 050"), Ok((1524.0, AltitudeRef::PressureStd)));
        assert!(matches!(alt("FL350 MSL"), Err(UnitError::Invalid { .. })));
        assert!(matches!(alt("FL3500"), Err(UnitError::Invalid { .. })));
        assert!(matches!(alt("1000 ft"), Err(UnitError::Invalid { .. })));
        assert!(matches!(alt("1000 ft msl"), Err(UnitError::Invalid { .. })));
        assert!(matches!(
            alt("1000 ft MSL AGL"),
            Err(UnitError::Invalid { .. })
        ));
        assert!(matches!(
            alt("1000 yd MSL"),
            Err(UnitError::UnknownUnit { .. })
        ));
        assert!(matches!(alt("1000"), Err(UnitError::MissingUnit { .. })));
        // Found by review: Display parses back, flight levels included.
        for text in [
            "FL350",
            "FL005",
            "35000 ft MSL",
            "-20 m ellipsoid",
            "0.3 km AGL",
        ] {
            let a: Altitude = text.parse().unwrap();
            assert_eq!(a.to_string().parse::<Altitude>(), Ok(a), "{text}");
        }
        assert_eq!("FL350".parse::<Altitude>().unwrap().to_string(), "FL350");
        let odd = Altitude::new(100.0, AltitudeRef::PressureStd).unwrap();
        assert_eq!(odd.to_string(), "100 m pressure altitude");
        let message = alt("1000 ft").unwrap_err().to_string();
        assert!(message.contains("MSL, AGL or ellipsoid"), "{message}");
    }
}
