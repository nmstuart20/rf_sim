//! Semantic validation beyond what the schema expresses.

#![forbid(unsafe_code)]

use iqsim_scenario::{Diagnostic, line_col, lint};

const HEAD: &str = "schema_version = 1\nepoch = 2026-01-01T00:00:00Z\nseed = 7\nrun_mode = \"realtime\"\n[metadata]\nname = \"t\"\n";

fn errors(body: &str) -> Vec<(usize, usize, String)> {
    let source = format!("{HEAD}{body}");
    lint(&source)
        .into_iter()
        .map(|d: Diagnostic| {
            let (line, col) = d.location(&source).unwrap();
            (line, col, d.message)
        })
        .collect()
}

#[test]
fn header_only_scenario_is_clean() {
    assert!(errors("").is_empty());
}

#[test]
fn scoped_duplicates_and_empty_channels() {
    let body = r#"
[[platforms]]
name = "p"
position = { lat = 0, lon = 0, alt = "0 m ellipsoid" }
mounts = [{ name = "m", antenna = "iso" }, { name = "m", antenna = "dipole" }]

[[platforms]]
name = "q"
position = { lat = 1, lon = 1, alt = "0 m ellipsoid" }
mounts = [{ name = "m", antenna = "iso" }]

[[receivers]]
name = "r"
platform = "p"
channels = []

[[receivers]]
name = "r2"
platform = "q"

[[receivers.channels]]
name = "c"
center_frequency = "1030 MHz"
sample_rate = "20 MS/s"
bandwidth = "18 MHz"

[[receivers.channels]]
name = "c"
center_frequency = "0 Hz"
sample_rate = "20 MS/s"
bandwidth = "-1 MHz"
"#;
    let found = errors(body);
    let messages: Vec<_> = found.iter().map(|(_, _, m)| m.as_str()).collect();
    assert_eq!(
        messages,
        [
            "duplicate name \"m\" among the mounts of platform \"p\"",
            "receiver \"r\" has no channels; add at least one [[receivers.channels]]",
            "duplicate name \"c\" among the channels of receiver \"r2\"",
            "the centre frequency must be positive",
            "the bandwidth must be positive",
        ],
        "{found:#?}"
    );
    // Mount "m" on platform "q" is not a duplicate of p's mount "m"; errors are in
    // source order.
    assert!(
        found
            .windows(2)
            .all(|w| (w[0].0, w[0].1) <= (w[1].0, w[1].1))
    );
}

/// The response needs 0 < bandwidth < stopband ≤ rate, compared exactly: a bandwidth
/// equal to a fractional rate is caught (no transition band below the default
/// stopband), and one a micro-hertz above it exceeds the rate. (Prompt 4 allowed
/// bandwidth = rate; prompt 6 needs a transition band to build H_rx.)
#[test]
fn bandwidth_is_compared_exactly_against_a_fractional_sample_rate() {
    let body = r#"
[[platforms]]
name = "p"
position = { lat = 0, lon = 0, alt = "0 m MSL" }
[[antenna_ports]]
name = "port"
platform = "p"
center_frequency = "1090 MHz"
sample_rate = "2.5 S/s"
bandwidth = "2.5 Hz"
[[antenna_ports]]
name = "port2"
platform = "p"
center_frequency = "1090 MHz"
sample_rate = "2.5 S/s"
bandwidth = "2.500001 Hz"
[[antenna_ports]]
name = "port3"
platform = "p"
center_frequency = "1090 MHz"
sample_rate = "2.5 S/s"
bandwidth = "2.499999 Hz"
[[antenna_ports]]
name = "port4"
platform = "p"
center_frequency = "1090 MHz"
sample_rate = "2.5 S/s"
bandwidth = "2 Hz"
response = { stopband = "2 Hz" }
"#;
    let found = errors(body);
    assert_eq!(found.len(), 3, "{found:#?}");
    assert!(found[0].2.contains(
        "antenna port \"port\": bandwidth 2.5 Hz leaves no transition band below the default stopband"
    ));
    assert!(found[1].2.contains(
        "antenna port \"port2\": bandwidth 2.500001 Hz exceeds the complex sample rate 2.5 S/s"
    ));
    assert!(
        found[2]
            .2
            .contains("stopband 2 Hz must be wider than the bandwidth 2 Hz")
    );
}

#[test]
fn reference_kinds_and_mounts() {
    let body = r#"
[[platforms]]
name = "p"
position = { lat = 0, lon = 0, alt = "0 m MSL" }
[[emitters]]
name = "e"
platform = "e"
mount = "nose"
frequency = "1090 MHz"
power = "1 W"
behavior = "lib/always_on"
waveform = "cw"
[[sinks]]
name = "s"
kind = "vita49"
resources = ["p", "nic9"]
"#;
    let messages: Vec<_> = errors(body).into_iter().map(|(_, _, m)| m).collect();
    assert_eq!(
        messages,
        [
            "emitter \"e\" refers to \"e\", which is an emitter, not a platform",
            "sink \"s\" refers to \"p\", which is a platform, not a resource",
            "sink \"s\" refers to unknown resource \"nic9\"",
        ]
    );
}

#[test]
fn syntax_errors_are_located() {
    let source = "schema_version = 1\nepoch = 2026-01-01T00:00:00Z\nseed = -1\n";
    let found = lint(source);
    assert_eq!(found.len(), 1);
    assert_eq!(found[0].location(source), Some((3, 8)));
    assert_eq!(line_col(source, 0), (1, 1));
}

#[test]
fn iq_file_emitters() {
    let platform =
        "[[platforms]]\nname = \"p\"\nposition = { lat = 0, lon = 0, alt = \"0 m ellipsoid\" }\n";
    let emitter = |name: &str, waveform: &str, iq: &str| {
        format!(
            "[[emitters]]\nname = \"{name}\"\nplatform = \"p\"\nfrequency = \"1 GHz\"\npower = \"0 dBm\"\nbehavior = \"always_on\"\nwaveform = \"{waveform}\"\n{iq}\n"
        )
    };
    let clean = format!(
        "{platform}{}",
        emitter(
            "a",
            "iq_file",
            r#"iq_file = { path = "rec", start = "1 ms", end = "loop", normalization = "rms", sample_rate = "2 MS/s", response = { bandwidth = "1.5 MHz", stopband = "2 MHz", attenuation = "90 dB" } }"#
        )
    );
    assert!(errors(&clean).is_empty(), "{:?}", errors(&clean));
    let body = format!(
        "{platform}{}{}{}{}",
        emitter("missing", "iq_file", ""),
        emitter(
            "stray",
            "cw",
            r#"iq_file = { path = "rec", end = "stop", normalization = "full_scale" }"#
        ),
        emitter(
            "bad",
            "iq_file",
            r#"iq_file = { path = " ", start = "-1 ms", end = "zero_pad", normalization = "rms", sample_rate = "2 MS/s", response = { bandwidth = "3 MHz", stopband = "1 MHz", attenuation = "10 dB" } }"#
        ),
        emitter(
            "ok",
            "iq_file",
            r#"iq_file = { path = "x", end = "loop", normalization = "rms" }"#
        ),
    );
    let messages: Vec<String> = errors(&body).into_iter().map(|(_, _, m)| m).collect();
    for needle in [
        "emitter \"missing\": waveform \"iq_file\" needs an iq_file table",
        "emitter \"stray\": an iq_file table needs waveform = \"iq_file\"",
        "emitter \"bad\": iq_file.path is empty",
        "emitter \"bad\": iq_file.start cannot be negative",
        "emitter \"bad\": iq_file response bandwidth exceeds the native sample rate",
        "emitter \"bad\": iq_file response stopband must be wider than its bandwidth",
        "emitter \"bad\": iq_file response attenuation must be 20 to 200 dB",
    ] {
        assert!(
            messages.iter().any(|m| m.contains(needle)),
            "missing {needle:?} in {messages:?}"
        );
    }
    assert_eq!(messages.len(), 7, "{messages:?}");
    // Unknown end policy or normalization is a parse error.
    let typo = format!(
        "{platform}{}",
        emitter(
            "t",
            "iq_file",
            r#"iq_file = { path = "x", end = "forever", normalization = "rms" }"#
        )
    );
    assert!(iqsim_scenario::parse(&format!("{HEAD}{typo}")).is_err());
}

#[test]
fn examples_lint_clean() {
    for name in ["minimal.toml", "iq_file.toml"] {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../examples")
            .join(name);
        let source = std::fs::read_to_string(&path).unwrap();
        let found = lint(&source);
        assert!(found.is_empty(), "{name}: {:?}", found);
    }
}

#[test]
fn world_table() {
    assert!(
        errors("[world]\ngeoid = { path = \"egm96-5.pgm\" }\npressure_offset = \"-35.5 hPa\"\n")
            .is_empty()
    );
    assert!(errors("[world]\npressure_offset = \"29.92 inHg\"\n").len() == 1);
    let e = errors("[world]\ngeoid = { path = \" \" }\npressure_offset = \"201 hPa\"\n");
    let m: Vec<&str> = e.iter().map(|(_, _, m)| m.as_str()).collect();
    assert_eq!(
        m,
        [
            "world.geoid.path is empty",
            "world.pressure_offset must be within ±200 hPa"
        ]
    );
    // Line 8 of the full source is `geoid = …`; the span starts at the path string.
    assert_eq!((e[0].0, e[0].1), (8, 18));
    // Units are checked like every other quantity, and unknown fields are rejected.
    assert!(!errors("[world]\npressure_offset = \"5 hpa\"\n").is_empty());
    assert!(!errors("[world]\ngeoid_path = \"x\"\n").is_empty());
}
