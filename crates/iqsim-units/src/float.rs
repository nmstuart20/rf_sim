//! Control-plane quantities kept as f64 in SI units: power, gain, distance, speed.
//!
//! Each converts from its exact decimal text with as few roundings as possible:
//! decimal unit factors (km, mW, ft = 0.3048 m, nmi = 1852 m) are applied exactly
//! before the single correctly rounded conversion to f64. dB conversions use `libm`.

use std::fmt;

use crate::decimal::Decimal;
use crate::error::{QuantityKind, UnitError};
use crate::text::{TextQuantity, finite, number_and_unit, text_quantity_impls, unit_pattern};

fn unknown(kind: &'static QuantityKind, unit: &str) -> UnitError {
    UnitError::UnknownUnit {
        kind,
        unit: unit.to_owned(),
    }
}

fn overflow(kind: &'static QuantityKind, text: &str) -> UnitError {
    UnitError::OutOfRange {
        kind,
        text: text.to_owned(),
        limit: "the value overflows a 64-bit float".to_owned(),
    }
}

/// `value · factor · 10^shift` exactly, then rounded once to f64. If the exact
/// product needs more than 38 digits, it rounds twice instead (decimal to f64, then
/// the multiply).
fn exact_f64(
    kind: &'static QuantityKind,
    text: &str,
    value: Decimal,
    factor: u128,
    shift: i32,
) -> Result<f64, UnitError> {
    let shifted = value.shifted(shift).ok_or_else(|| overflow(kind, text))?;
    let result = match shifted.mul_u128(factor) {
        Some(exact) => exact.to_f64(),
        None => shifted.to_f64() * factor as f64,
    };
    finite(kind, text, result)
}

// ---------------------------------------------------------------- Power

static POWER: QuantityKind = QuantityKind {
    name: "power",
    units: &["W", "mW", "dBm", "dBW"],
    example: "30 dBm",
};

/// A non-negative power in watts.
///
/// Units: W, mW, dBm, dBW. dB values convert with `libm::exp10`.
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct Power {
    watts: f64,
}

impl Power {
    /// A power of `watts`. `None` unless finite and non-negative.
    pub fn from_watts(watts: f64) -> Option<Self> {
        (watts.is_finite() && watts >= 0.0).then_some(Self { watts })
    }

    /// The power in watts.
    pub fn watts(self) -> f64 {
        self.watts
    }

    /// The power in dBW (−∞ for zero).
    pub fn dbw(self) -> f64 {
        10.0 * libm::log10(self.watts)
    }

    /// The power in dBm (−∞ for zero).
    pub fn dbm(self) -> f64 {
        self.dbw() + 30.0
    }
}

impl TextQuantity for Power {
    const KIND: &'static QuantityKind = &POWER;

    fn parse_text(text: &str) -> Result<Self, UnitError> {
        let kind = Self::KIND;
        let (value, unit) = number_and_unit(kind, text)?;
        let watts = match unit {
            "W" => exact_f64(kind, text, value, 1, 0)?,
            "mW" => exact_f64(kind, text, value, 1, -3)?,
            "dBW" => libm::exp10(finite(kind, text, value.to_f64())? / 10.0),
            "dBm" => libm::exp10((finite(kind, text, value.to_f64())? - 30.0) / 10.0),
            _ => return Err(unknown(kind, unit)),
        };
        if value.is_negative() && matches!(unit, "W" | "mW") {
            return Err(UnitError::OutOfRange {
                kind,
                text: text.to_owned(),
                limit: "a power in W or mW cannot be negative (use dBm or dBW for small powers)"
                    .to_owned(),
            });
        }
        Ok(Self {
            watts: finite(kind, text, watts)?,
        })
    }
}

impl fmt::Display for Power {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} W", self.watts)
    }
}

text_quantity_impls!(Power, Some(unit_pattern(POWER.units)));

// ---------------------------------------------------------------- Gain

static GAIN: QuantityKind = QuantityKind {
    name: "gain",
    units: &["dB", "dBi"],
    example: "3 dBi",
};

/// A signed gain in decibels. Units: dB, dBi (both stored as dB).
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct Gain {
    db: f64,
}

impl Gain {
    /// A gain of `db` decibels. `None` unless finite.
    pub fn from_db(db: f64) -> Option<Self> {
        db.is_finite().then_some(Self { db })
    }

    /// The gain in dB.
    pub fn db(self) -> f64 {
        self.db
    }

    /// The linear power ratio, 10^(dB/10).
    pub fn power_ratio(self) -> f64 {
        libm::exp10(self.db / 10.0)
    }
}

impl TextQuantity for Gain {
    const KIND: &'static QuantityKind = &GAIN;

    fn parse_text(text: &str) -> Result<Self, UnitError> {
        let (value, unit) = number_and_unit(Self::KIND, text)?;
        match unit {
            "dB" | "dBi" => Ok(Self {
                db: finite(Self::KIND, text, value.to_f64())?,
            }),
            _ => Err(unknown(Self::KIND, unit)),
        }
    }
}

impl fmt::Display for Gain {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} dB", self.db)
    }
}

text_quantity_impls!(Gain, Some(unit_pattern(GAIN.units)));

// ---------------------------------------------------------------- Distance

static DISTANCE: QuantityKind = QuantityKind {
    name: "distance",
    units: &["m", "km", "ft", "nmi"],
    example: "20 nmi",
};

/// A signed distance in metres. Units: m, km, ft (0.3048 m), nmi (1852 m).
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct Distance {
    meters: f64,
}

impl Distance {
    /// A distance of `meters`. `None` unless finite.
    pub fn from_meters(meters: f64) -> Option<Self> {
        meters.is_finite().then_some(Self { meters })
    }

    /// The distance in metres.
    pub fn meters(self) -> f64 {
        self.meters
    }
}

/// Metres from a decimal in `unit`, if it is a length unit.
pub(crate) fn length_meters(
    kind: &'static QuantityKind,
    text: &str,
    value: Decimal,
    unit: &str,
) -> Option<Result<f64, UnitError>> {
    Some(match unit {
        "m" => exact_f64(kind, text, value, 1, 0),
        "km" => exact_f64(kind, text, value, 1, 3),
        "ft" => exact_f64(kind, text, value, 3048, -4),
        "nmi" => exact_f64(kind, text, value, 1852, 0),
        _ => return None,
    })
}

impl TextQuantity for Distance {
    const KIND: &'static QuantityKind = &DISTANCE;

    fn parse_text(text: &str) -> Result<Self, UnitError> {
        let (value, unit) = number_and_unit(Self::KIND, text)?;
        let meters = length_meters(Self::KIND, text, value, unit)
            .ok_or_else(|| unknown(Self::KIND, unit))??;
        Ok(Self { meters })
    }
}

impl fmt::Display for Distance {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} m", self.meters)
    }
}

text_quantity_impls!(Distance, Some(unit_pattern(DISTANCE.units)));

// ---------------------------------------------------------------- Speed

static SPEED: QuantityKind = QuantityKind {
    name: "speed",
    units: &["m/s", "km/h", "kt"],
    example: "250 kt",
};

/// A signed speed in metres per second. Units: m/s, km/h, kt (1852/3600 m/s).
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct Speed {
    meters_per_second: f64,
}

impl Speed {
    /// A speed of `meters_per_second`. `None` unless finite.
    pub fn from_meters_per_second(meters_per_second: f64) -> Option<Self> {
        meters_per_second
            .is_finite()
            .then_some(Self { meters_per_second })
    }

    /// The speed in metres per second.
    pub fn meters_per_second(self) -> f64 {
        self.meters_per_second
    }
}

impl TextQuantity for Speed {
    const KIND: &'static QuantityKind = &SPEED;

    fn parse_text(text: &str) -> Result<Self, UnitError> {
        let kind = Self::KIND;
        let (value, unit) = number_and_unit(kind, text)?;
        // Exact numerator, then one division by an exactly representable integer.
        let meters_per_second = match unit {
            "m/s" => exact_f64(kind, text, value, 1, 0)?,
            "km/h" => exact_f64(kind, text, value, 1000, 0)? / 3600.0,
            "kt" => exact_f64(kind, text, value, 1852, 0)? / 3600.0,
            _ => return Err(unknown(kind, unit)),
        };
        Ok(Self { meters_per_second })
    }
}

impl fmt::Display for Speed {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} m/s", self.meters_per_second)
    }
}

text_quantity_impls!(Speed, Some(unit_pattern(SPEED.units)));

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn power_units() {
        let w = |s: &str| s.parse::<Power>().map(Power::watts);
        assert_eq!(w("2.5 W"), Ok(2.5));
        assert_eq!(w("0.1 mW"), Ok(0.0001));
        assert_eq!(w("0 dBW"), Ok(1.0));
        assert_eq!(w("30 dBm"), Ok(1.0));
        assert_eq!(w("20 dBW"), Ok(100.0));
        assert!((w("-3 dBm").unwrap() - 5.011_872_336_272_722e-4).abs() < 1e-18);
        assert!(matches!(w("-1 W"), Err(UnitError::OutOfRange { .. })));
        assert!(matches!(w("1e400 W"), Err(UnitError::OutOfRange { .. })));
        assert!(matches!(w("1e400 dBm"), Err(UnitError::OutOfRange { .. })));
        assert!(matches!(w("10 dbm"), Err(UnitError::UnknownUnit { .. })));
        let p: Power = "36 dBm".parse().unwrap();
        assert!((p.dbm() - 36.0).abs() < 1e-12);
    }

    #[test]
    fn gain_distance_speed_units() {
        assert_eq!("-3 dBi".parse::<Gain>().unwrap().db(), -3.0);
        assert_eq!("10 dB".parse::<Gain>().unwrap().power_ratio(), 10.0);
        let m = |s: &str| s.parse::<Distance>().map(Distance::meters);
        assert_eq!(m("1 ft"), Ok(0.3048));
        assert_eq!(m("35000 ft"), Ok(10_668.0));
        assert_eq!(m("2.5 km"), Ok(2500.0));
        assert_eq!(m("20 nmi"), Ok(37_040.0));
        // Found by review: a 38-digit mantissa times 3048 no longer overflows.
        let long = m("1.2345678901234567890123456789012345678 ft").unwrap();
        assert!((long - 1.234_567_890_123_456_7 * 0.3048).abs() < 1e-15);
        assert!(matches!(m("3 yd"), Err(UnitError::UnknownUnit { .. })));
        let v = |s: &str| s.parse::<Speed>().map(Speed::meters_per_second);
        assert_eq!(v("250 kt"), Ok(463_000.0 / 3600.0));
        assert_eq!(v("36 km/h"), Ok(10.0));
        assert_eq!(v("-5 m/s"), Ok(-5.0));
        assert!(matches!(v("250 knots"), Err(UnitError::UnknownUnit { .. })));
    }
}
