//! The checked-in example scenario parses and validates clean.

#![forbid(unsafe_code)]

use iqsim_scenario::{RunMode, SinkKind, lint, parse};

const MINIMAL: &str = include_str!("../../../examples/minimal.toml");

#[test]
fn minimal_example_is_clean() {
    let diagnostics = lint(MINIMAL);
    assert!(
        diagnostics.is_empty(),
        "{}",
        diagnostics
            .iter()
            .map(|d| d.render("examples/minimal.toml", MINIMAL))
            .collect::<String>()
    );
    let s = parse(MINIMAL).unwrap();
    assert_eq!(s.epoch.unix_seconds(), 1_767_225_600);
    assert_eq!(s.run_mode, RunMode::Offline);
    assert_eq!(
        s.emitters[0].frequency.get_ref().as_micro_hz(),
        1_090_000_000_000_000
    );
    assert!((s.emitters[0].power.watts() - 1.0).abs() < 1e-15);
    let channel = &s.receivers[0].channels[0];
    assert_eq!(channel.sample_rate.get_ref().get().num(), 10_000_000);
    assert_eq!(s.sinks[0].kind, SinkKind::Sigmf);
    assert_eq!(
        s.resources[0].capacity.unwrap().bits_per_second(),
        3_200_000_000
    );
    let site = s.platforms[1].position.as_ref().unwrap();
    assert_eq!(site.mgrs().unwrap().zone(), 18);
}
