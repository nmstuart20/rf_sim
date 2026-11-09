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
bandwidth = "20 MHz"

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

#[test]
fn bandwidth_may_equal_a_fractional_sample_rate() {
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
"#;
    let found = errors(body);
    assert_eq!(found.len(), 1, "{found:#?}");
    assert!(found[0].2.contains(
        "antenna port \"port2\": bandwidth 2.500001 Hz exceeds the complex sample rate 2.5 S/s"
    ));
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
