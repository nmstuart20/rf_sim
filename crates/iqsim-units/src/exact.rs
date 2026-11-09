//! Quantities kept exactly as integers: frequency, sample rate, duration, data rate.
//!
//! These values reach the render path (carrier frequency, sample clock, event
//! times) or are summed for capacity checks, so the decimal text converts to the
//! integer representation with no binary rounding. A value finer than the
//! representation's resolution is an error, never silently rounded.

use std::fmt;

use iqsim_core::{TimeError, TruthDuration};

use crate::decimal::ScaleError;
use crate::error::{QuantityKind, UnitError};
use crate::text::{TextQuantity, number_and_unit, text_quantity_impls, unit_pattern};

fn unknown(kind: &'static QuantityKind, unit: &str) -> UnitError {
    UnitError::UnknownUnit {
        kind,
        unit: unit.to_owned(),
    }
}

/// Formats `value / 10^digits` exactly, without trailing fractional zeros.
fn fmt_fixed(f: &mut fmt::Formatter<'_>, value: i128, digits: u32) -> fmt::Result {
    let scale = 10i128.pow(digits);
    let sign = if value < 0 { "-" } else { "" };
    let magnitude = value.unsigned_abs();
    let whole = magnitude / scale.unsigned_abs();
    let frac = magnitude % scale.unsigned_abs();
    if frac == 0 {
        write!(f, "{sign}{whole}")
    } else {
        let width = usize::try_from(digits).expect("small");
        let frac = format!("{frac:0width$}");
        write!(f, "{sign}{whole}.{}", frac.trim_end_matches('0'))
    }
}

// ---------------------------------------------------------------- Frequency

static FREQUENCY: QuantityKind = QuantityKind {
    name: "frequency",
    units: &["Hz", "kHz", "MHz", "GHz"],
    example: "1090 MHz",
};

/// A frequency, exact to 1 µHz, as a signed i64 count of micro-hertz (±9.2 THz).
///
/// Units: Hz, kHz, MHz, GHz. Negative values are allowed here (offsets); fields
/// that need a positive frequency check it.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Frequency {
    micro_hz: i64,
}

impl Frequency {
    /// A frequency of `micro_hz` micro-hertz.
    pub const fn from_micro_hz(micro_hz: i64) -> Self {
        Self { micro_hz }
    }

    /// A whole number of hertz, or `None` if it exceeds the range.
    pub const fn from_hz(hz: i64) -> Option<Self> {
        match hz.checked_mul(1_000_000) {
            Some(micro_hz) => Some(Self { micro_hz }),
            None => None,
        }
    }

    /// The exact value in micro-hertz.
    pub const fn as_micro_hz(self) -> i64 {
        self.micro_hz
    }

    /// The nearest f64 in hertz (for control-plane arithmetic only).
    pub fn as_hz_f64(self) -> f64 {
        self.micro_hz as f64 / 1e6
    }

    /// Whether the frequency is above zero.
    pub const fn is_positive(self) -> bool {
        self.micro_hz > 0
    }
}

impl TextQuantity for Frequency {
    const KIND: &'static QuantityKind = &FREQUENCY;

    fn parse_text(text: &str) -> Result<Self, UnitError> {
        let (value, unit) = number_and_unit(Self::KIND, text)?;
        let scale = match unit {
            "Hz" => 6,
            "kHz" => 9,
            "MHz" => 12,
            "GHz" => 15,
            _ => return Err(unknown(Self::KIND, unit)),
        };
        match value.scaled_i128(scale).map(i64::try_from) {
            Ok(Ok(micro_hz)) => Ok(Self { micro_hz }),
            Err(ScaleError::Inexact) => Err(UnitError::Inexact {
                kind: Self::KIND,
                text: text.to_owned(),
                resolution: "1 µHz",
            }),
            _ => Err(UnitError::OutOfRange {
                kind: Self::KIND,
                text: text.to_owned(),
                limit: "the magnitude must be below 9.2 THz".to_owned(),
            }),
        }
    }
}

impl fmt::Display for Frequency {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_fixed(f, i128::from(self.micro_hz), 6)?;
        f.write_str(" Hz")
    }
}

text_quantity_impls!(Frequency, Some(unit_pattern(FREQUENCY.units)));

// ---------------------------------------------------------------- SampleRate

static SAMPLE_RATE: QuantityKind = QuantityKind {
    name: "sample rate",
    units: &["S/s", "kS/s", "MS/s", "GS/s"],
    example: "61.44 MS/s",
};

/// A sample rate, exact, as the core's reduced rational [`iqsim_core::SampleRate`].
///
/// Units: S/s, kS/s, MS/s, GS/s. `"61.44 MS/s"` is exactly 61 440 000 S/s. The
/// core's bounds apply (denominator ≤ 2³², numerator ≤ 2⁴⁸, at most 1 GS/s).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SampleRate(iqsim_core::SampleRate);

impl SampleRate {
    /// Wraps a core sample rate.
    pub const fn new(rate: iqsim_core::SampleRate) -> Self {
        Self(rate)
    }

    /// The core sample rate.
    pub const fn get(self) -> iqsim_core::SampleRate {
        self.0
    }
}

const fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        let r = a % b;
        a = b;
        b = r;
    }
    a
}

impl TextQuantity for SampleRate {
    const KIND: &'static QuantityKind = &SAMPLE_RATE;

    fn parse_text(text: &str) -> Result<Self, UnitError> {
        let (value, unit) = number_and_unit(Self::KIND, text)?;
        let scale = match unit {
            "S/s" => 0,
            "kS/s" => 3,
            "MS/s" => 6,
            "GS/s" => 9,
            _ => return Err(unknown(Self::KIND, unit)),
        };
        let out_of_range = |limit: &str| UnitError::OutOfRange {
            kind: Self::KIND,
            text: text.to_owned(),
            limit: limit.to_owned(),
        };
        if value.is_negative() || value.is_zero() {
            return Err(out_of_range("a sample rate must be positive"));
        }
        let bounds = "rates go up to 1 GS/s, with an exact denominator up to 2^32 (at most 9 decimal places in S/s)";
        let (num, den) = value.ratio(scale).ok_or_else(|| out_of_range(bounds))?;
        let g = gcd(num, den);
        let (num, den) = (num / g, den / g);
        let (Ok(num), Ok(den)) = (u64::try_from(num), u64::try_from(den)) else {
            return Err(out_of_range(bounds));
        };
        match iqsim_core::SampleRate::from_ratio(num, den) {
            Ok(rate) => Ok(Self(rate)),
            Err(TimeError::RateOutOfBounds { .. }) => Err(out_of_range(bounds)),
            Err(_) => Err(out_of_range("a sample rate must be positive")),
        }
    }
}

/// The `k` with `den | 10^k` (smallest), if `den` has no prime factors but 2 and 5.
fn decimal_places(den: u64) -> Option<u32> {
    let (mut d, mut twos, mut fives) = (den, 0, 0);
    while d.is_multiple_of(2) {
        d /= 2;
        twos += 1;
    }
    while d.is_multiple_of(5) {
        d /= 5;
        fives += 1;
    }
    (d == 1).then_some(twos.max(fives))
}

/// Exact decimal S/s when the denominator allows it (always true for rates parsed
/// from text, so those round-trip); otherwise `num/den S/s`, which does not parse.
impl fmt::Display for SampleRate {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let rate = self.0;
        let scaled = decimal_places(rate.den()).and_then(|k| {
            let factor = 10u128.checked_pow(k)? / u128::from(rate.den());
            let value = i128::try_from(u128::from(rate.num()).checked_mul(factor)?).ok()?;
            Some((value, k))
        });
        match scaled {
            Some((value, k)) => {
                fmt_fixed(f, value, k)?;
                f.write_str(" S/s")
            }
            None => write!(f, "{}/{} S/s", rate.num(), rate.den()),
        }
    }
}

text_quantity_impls!(SampleRate, Some(unit_pattern(SAMPLE_RATE.units)));

// ---------------------------------------------------------------- Duration

static DURATION: QuantityKind = QuantityKind {
    name: "duration",
    units: &["s", "ms", "us", "µs", "ns", "ps"],
    example: "20.3 us",
};

/// A signed duration, exact to 1 ps, as the core's [`TruthDuration`] (±106 days).
///
/// Units: s, ms, us (or µs), ns, ps.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Duration(TruthDuration);

impl Duration {
    /// Wraps a core duration.
    pub const fn new(duration: TruthDuration) -> Self {
        Self(duration)
    }

    /// The core duration.
    pub const fn get(self) -> TruthDuration {
        self.0
    }
}

impl TextQuantity for Duration {
    const KIND: &'static QuantityKind = &DURATION;

    fn parse_text(text: &str) -> Result<Self, UnitError> {
        let (value, unit) = number_and_unit(Self::KIND, text)?;
        let scale = match unit {
            "s" => 12,
            "ms" => 9,
            "us" | "µs" | "μs" => 6,
            "ns" => 3,
            "ps" => 0,
            _ => return Err(unknown(Self::KIND, unit)),
        };
        match value.scaled_i128(scale).map(i64::try_from) {
            Ok(Ok(ps)) => Ok(Self(TruthDuration::from_ps(ps))),
            Err(ScaleError::Inexact) => Err(UnitError::Inexact {
                kind: Self::KIND,
                text: text.to_owned(),
                resolution: "1 ps",
            }),
            _ => Err(UnitError::OutOfRange {
                kind: Self::KIND,
                text: text.to_owned(),
                limit: "the magnitude must be below about 106 days (2^63 ps)".to_owned(),
            }),
        }
    }
}

impl fmt::Display for Duration {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_fixed(f, i128::from(self.0.as_ps()), 12)?;
        f.write_str(" s")
    }
}

text_quantity_impls!(
    Duration,
    Some(unit_pattern(&["s", "ms", "us", "µs", "μs", "ns", "ps"]))
);

// ---------------------------------------------------------------- DataRate

static DATA_RATE: QuantityKind = QuantityKind {
    name: "data rate",
    units: &["b/s", "kb/s", "Mb/s", "Gb/s", "B/s", "kB/s", "MB/s", "GB/s"],
    example: "1 Gb/s",
};

/// A non-negative data rate, exact, in whole bits per second (u64).
///
/// Units (SI, powers of 1000): b/s, kb/s, Mb/s, Gb/s for bits and B/s, kB/s, MB/s,
/// GB/s for bytes (8 bits).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DataRate {
    bits_per_second: u64,
}

impl DataRate {
    /// A rate of `bits_per_second`.
    pub const fn from_bits_per_second(bits_per_second: u64) -> Self {
        Self { bits_per_second }
    }

    /// The rate in bits per second.
    pub const fn bits_per_second(self) -> u64 {
        self.bits_per_second
    }
}

impl TextQuantity for DataRate {
    const KIND: &'static QuantityKind = &DATA_RATE;

    fn parse_text(text: &str) -> Result<Self, UnitError> {
        let (value, unit) = number_and_unit(Self::KIND, text)?;
        let (scale, bits) = match unit {
            "b/s" => (0, 1),
            "kb/s" => (3, 1),
            "Mb/s" => (6, 1),
            "Gb/s" => (9, 1),
            "B/s" => (0, 8),
            "kB/s" => (3, 8),
            "MB/s" => (6, 8),
            "GB/s" => (9, 8),
            _ => return Err(unknown(Self::KIND, unit)),
        };
        let out_of_range = |limit: &str| UnitError::OutOfRange {
            kind: Self::KIND,
            text: text.to_owned(),
            limit: limit.to_owned(),
        };
        if value.is_negative() {
            return Err(out_of_range("a data rate cannot be negative"));
        }
        let scaled = value
            .mul_u128(bits)
            .ok_or_else(|| out_of_range("too large"))?
            .scaled_i128(scale);
        match scaled.map(u64::try_from) {
            Ok(Ok(bits_per_second)) => Ok(Self { bits_per_second }),
            Err(ScaleError::Inexact) => Err(UnitError::Inexact {
                kind: Self::KIND,
                text: text.to_owned(),
                resolution: "1 b/s",
            }),
            _ => Err(out_of_range("the rate must fit in 64 bits of b/s")),
        }
    }
}

impl fmt::Display for DataRate {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} b/s", self.bits_per_second)
    }
}

text_quantity_impls!(DataRate, Some(unit_pattern(DATA_RATE.units)));

#[cfg(test)]
mod tests {
    use super::*;

    fn freq(s: &str) -> Result<Frequency, UnitError> {
        s.parse()
    }

    #[test]
    fn frequency_is_exact() {
        assert_eq!(
            freq("1090 MHz").unwrap().as_micro_hz(),
            1_090_000_000_000_000
        );
        assert_eq!(freq("1090MHz").unwrap(), freq(" 1.09 GHz ").unwrap());
        assert_eq!(freq("0.000001 Hz").unwrap().as_micro_hz(), 1);
        assert_eq!(freq("-2.5 kHz").unwrap().as_micro_hz(), -2_500_000_000);
        assert_eq!(
            freq("1030.000000001 MHz").unwrap().as_micro_hz(),
            1_030_000_000_001_000
        );
        assert!(matches!(
            freq("0.0000001 Hz"),
            Err(UnitError::Inexact { .. })
        ));
        assert!(matches!(
            freq("10000 GHz"),
            Err(UnitError::OutOfRange { .. })
        ));
        assert!(matches!(freq("1090"), Err(UnitError::MissingUnit { .. })));
        assert!(matches!(
            freq("1090 mhz"),
            Err(UnitError::UnknownUnit { .. })
        ));
        assert!(matches!(freq("10 dBm"), Err(UnitError::UnknownUnit { .. })));
        assert!(matches!(freq("MHz"), Err(UnitError::InvalidNumber { .. })));
        assert!(matches!(freq("  "), Err(UnitError::Empty { .. })));
        assert_eq!(freq("1090.5 MHz").unwrap().to_string(), "1090500000 Hz");
        assert_eq!(freq("-0.25 Hz").unwrap().to_string(), "-0.25 Hz");
    }

    #[test]
    fn sample_rate_is_exact_rational() {
        let r: SampleRate = "61.44 MS/s".parse().unwrap();
        assert_eq!((r.get().num(), r.get().den()), (61_440_000, 1));
        let r: SampleRate = "2.5 S/s".parse().unwrap();
        assert_eq!((r.get().num(), r.get().den()), (5, 2));
        let r: SampleRate = "0.001 S/s".parse().unwrap();
        assert_eq!((r.get().num(), r.get().den()), (1, 1000));
        assert_eq!(r.to_string(), "0.001 S/s");
        assert_eq!(
            "2.5 S/s".parse::<SampleRate>().unwrap().to_string(),
            "2.5 S/s"
        );
        let third = SampleRate::new(iqsim_core::SampleRate::from_ratio(1_000_000_000, 3).unwrap());
        assert_eq!(third.to_string(), "1000000000/3 S/s");
        assert!("1 GS/s".parse::<SampleRate>().is_ok());
        for bad in ["1.000000001 GS/s", "0 S/s", "-1 MS/s", "0.0000000001 S/s"] {
            assert!(
                matches!(bad.parse::<SampleRate>(), Err(UnitError::OutOfRange { .. })),
                "{bad}"
            );
        }
        assert!(matches!(
            "10 MHz".parse::<SampleRate>(),
            Err(UnitError::UnknownUnit { .. })
        ));
    }

    #[test]
    fn duration_is_exact_picoseconds() {
        let d = |s: &str| s.parse::<Duration>().map(|d| d.get().as_ps());
        assert_eq!(d("20.3 us"), Ok(20_300_000));
        assert_eq!(d("20.3 µs"), d("20.3us"));
        assert_eq!(d("1.5 s"), Ok(1_500_000_000_000));
        assert_eq!(d("-3 ns"), Ok(-3000));
        assert_eq!(d("1 ps"), Ok(1));
        assert!(matches!(d("0.5 ps"), Err(UnitError::Inexact { .. })));
        assert!(matches!(d("107 days"), Err(UnitError::UnknownUnit { .. })));
        assert!(matches!(d("9300000 s"), Err(UnitError::OutOfRange { .. })));
        assert_eq!(
            "2.5 ms".parse::<Duration>().unwrap().to_string(),
            "0.0025 s"
        );
    }

    #[test]
    fn data_rate_is_exact_bits() {
        let r = |s: &str| s.parse::<DataRate>().map(DataRate::bits_per_second);
        assert_eq!(r("1 Gb/s"), Ok(1_000_000_000));
        assert_eq!(r("1.5 B/s"), Ok(12));
        assert_eq!(r("1.25 GB/s"), Ok(10_000_000_000));
        assert!(matches!(r("0.5 b/s"), Err(UnitError::Inexact { .. })));
        assert!(matches!(r("-1 Mb/s"), Err(UnitError::OutOfRange { .. })));
        assert!(matches!(r("1 Gbps"), Err(UnitError::UnknownUnit { .. })));
    }
}
