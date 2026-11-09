//! Quantities through TOML, and properties of the exact parsers.

#![forbid(unsafe_code)]

use iqsim_units::{Altitude, Angle, Duration, Frequency, Latitude, Position, Power, SampleRate};
use proptest::prelude::*;
use proptest::test_runner::FileFailurePersistence;
use serde::Deserialize;

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct Sample {
    f: Option<Frequency>,
    p: Option<Power>,
    lat: Option<Latitude>,
    pos: Option<Position>,
}

fn load(text: &str) -> Result<Sample, String> {
    toml::from_str(text).map_err(|e| e.message().to_owned())
}

#[test]
fn bare_numbers_are_rejected_with_a_suggestion() {
    let message = load("f = 1090").unwrap_err();
    assert!(message.contains("frequency needs a unit"), "{message}");
    assert!(message.contains("bare number 1090"), "{message}");
    assert!(message.contains("\"1090 MHz\""), "{message}");
    assert!(message.contains("Hz, kHz, MHz, GHz"), "{message}");
    let message = load("p = 30.5").unwrap_err();
    assert!(
        message.contains("power needs a unit: got the bare number 30.5"),
        "{message}"
    );
    assert!(
        load("f = true")
            .unwrap_err()
            .contains("frequency string such as")
    );
}

#[test]
fn toml_values_parse() {
    let s = load("f = \"1090 MHz\"\np = \"0 dBW\"\nlat = -33.5").unwrap();
    assert_eq!(s.f.unwrap().as_micro_hz(), 1_090_000_000_000_000);
    assert_eq!(s.p.unwrap().watts(), 1.0);
    assert_eq!(s.lat.unwrap().degrees(), -33.5);
    let s = load("pos = { lat = \"38°53'23.3\\\"N\", lon = -77.0365, alt = \"FL350\" }").unwrap();
    let pos = s.pos.unwrap();
    assert!((pos.latitude().degrees() - (38.0 + 53.0 / 60.0 + 23.3 / 3600.0)).abs() < 1e-12);
    assert_eq!(pos.longitude().degrees(), -77.0365);
    assert!(pos.mgrs().is_none());
}

#[test]
fn position_table_errors() {
    for (text, fragment) in [
        (
            "pos = { lat = 1, alt = \"1 m MSL\" }",
            "both `lat` and `lon`",
        ),
        (
            "pos = { mgrs = \"18SUJ\", lat = 1, lon = 2, alt = \"1 m MSL\" }",
            "not both",
        ),
        ("pos = { alt = \"1 m MSL\" }", "`lat` and `lon`, or `mgrs`"),
        ("pos = { lat = 1, lon = 2 }", "missing field `alt`"),
        (
            "pos = { lat = 1, lon = 2, alt = \"1 m MSL\", height = 3 }",
            "unknown field `height`",
        ),
        (
            "pos = { lat = 91, lon = 2, alt = \"1 m MSL\" }",
            "within ±90",
        ),
        (
            "pos = { mgrs = \"18ZUJ\", alt = \"1 m MSL\" }",
            "polar (UPS)",
        ),
    ] {
        let message = load(text).unwrap_err();
        assert!(message.contains(fragment), "{text}: {message}");
    }
}

/// Cases per property unless `PROPTEST_CASES` is set (AGENTS.md).
const DEFAULT_CASES: u32 = 2048;

/// Proptest settings. Integration tests have no lib.rs beside them, so failures
/// are persisted to an explicit path (relative to the crate directory).
fn config() -> ProptestConfig {
    let mut config = ProptestConfig {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            "proptest-regressions/units_props.txt",
        ))),
        ..ProptestConfig::default()
    };
    if std::env::var_os("PROPTEST_CASES").is_none() {
        config.cases = DEFAULT_CASES;
    }
    config
}

proptest! {
    #![proptest_config(config())]

    /// Display then parse returns the identical frequency.
    #[test]
    fn frequency_display_round_trips(micro_hz in any::<i64>()) {
        let f = Frequency::from_micro_hz(micro_hz);
        prop_assert_eq!(f.to_string().parse::<Frequency>(), Ok(f));
    }

    /// A decimal in MHz with up to 12 fractional digits converts exactly.
    #[test]
    fn frequency_decimal_is_exact(whole in 0i64..9_000_000, frac in 0i64..1_000_000_000_000, digits in 0usize..=12) {
        let frac = frac % 10i64.pow(u32::try_from(digits).unwrap());
        let text = if digits == 0 {
            format!("{whole} MHz")
        } else {
            format!("{whole}.{frac:0digits$}MHz")
        };
        let expected = whole * 1_000_000_000_000 + frac * 10i64.pow(12 - u32::try_from(digits).unwrap());
        prop_assert_eq!(text.parse::<Frequency>().map(Frequency::as_micro_hz), Ok(expected));
    }

    /// Durations round-trip through Display at picosecond resolution.
    #[test]
    fn duration_display_round_trips(ps in any::<i64>()) {
        let d = Duration::new(iqsim_core::TruthDuration::from_ps(ps));
        prop_assert_eq!(d.to_string().parse::<Duration>(), Ok(d));
    }

    /// Integer rates in S/s, kS/s and MS/s agree with the core constructor.
    #[test]
    fn sample_rate_matches_core(hz in 1u64..=1_000_000_000) {
        let core = iqsim_core::SampleRate::from_hz(hz).unwrap();
        prop_assert_eq!(format!("{hz} S/s").parse::<SampleRate>().map(SampleRate::get), Ok(core));
        let khz = format!("{}.{:03}kS/s", hz / 1000, hz % 1000);
        prop_assert_eq!(khz.parse::<SampleRate>().map(SampleRate::get), Ok(core));
    }

    /// Any accepted rate with a decimal denominator (2^a·5^b) prints as an exact
    /// decimal that parses back to the same rate.
    #[test]
    fn sample_rate_display_round_trips(num in 1u64..(1 << 48), twos in 0u32..=32, fives in 0u32..=13) {
        let den = 2u64.pow(twos).checked_mul(5u64.pow(fives));
        if let Some(Ok(core)) = den.map(|den| iqsim_core::SampleRate::from_ratio(num, den)) {
            let rate = SampleRate::new(core);
            let text = rate.to_string();
            // num · 10^k / den ≤ 2^48 · 5^32 < 2^127, so the decimal form always fits.
            prop_assert!(!text.trim_end_matches(" S/s").contains('/'), "{}", text);
            prop_assert_eq!(text.parse::<SampleRate>(), Ok(rate));
        }
    }

    /// DMS equals d + m/60 + s/3600 (to f64 rounding) and the ASCII form agrees.
    #[test]
    fn dms_matches_decimal(d in 0u32..360, m in 0u32..60, ms in 0u32..60_000) {
        let s = f64::from(ms) / 1000.0;
        let symbols = format!("{d}°{m}'{}.{:03}\"", ms / 1000, ms % 1000);
        let ascii = format!("{d}d{m}m{}.{:03}s", ms / 1000, ms % 1000);
        let a = symbols.parse::<Angle>().unwrap().degrees();
        prop_assert_eq!(ascii.parse::<Angle>().map(Angle::degrees), Ok(a));
        let reference = f64::from(d) + f64::from(m) / 60.0 + s / 3600.0;
        prop_assert!((a - reference).abs() <= 1e-12, "{} vs {}", a, reference);
    }

    /// Feet and flight levels agree: FLn = n·100 ft MSL in metres.
    #[test]
    fn flight_level_is_hundreds_of_feet(n in 0u32..1000) {
        let fl = format!("FL{n:03}").parse::<Altitude>().unwrap().meters();
        let ft = format!("{} ft MSL", n * 100).parse::<Altitude>().unwrap().meters();
        prop_assert_eq!(fl, ft);
    }

    /// Arbitrary text never panics any parser.
    #[test]
    fn parsers_never_panic(text in "\\PC{0,24}") {
        let _ = text.parse::<Frequency>();
        let _ = text.parse::<SampleRate>();
        let _ = text.parse::<Duration>();
        let _ = text.parse::<Power>();
        let _ = text.parse::<Angle>();
        let _ = text.parse::<Altitude>();
        let _ = text.parse::<Latitude>();
        let _ = iqsim_units::Mgrs::parse(&text);
    }
}
