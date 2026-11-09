//! Verify the CLI through its executable interface.

#![forbid(unsafe_code)]

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn sim(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_sim"))
        .args(args)
        .output()
        .expect("run sim")
}

fn workspace_file(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(relative)
}

#[test]
fn version_flags_report_the_package_version() {
    for flag in ["--version", "-V"] {
        let output = sim(&[flag]);
        assert!(output.status.success());
        assert!(output.stderr.is_empty());
        assert_eq!(
            String::from_utf8(output.stdout).unwrap(),
            concat!("sim ", env!("CARGO_PKG_VERSION"), "\n")
        );
    }
}

#[test]
fn help_and_no_arguments_list_commands() {
    for args in [vec![], vec!["--help"], vec!["-h"]] {
        let output = sim(&args);
        assert!(output.status.success());
        assert!(output.stderr.is_empty());
        let stdout = String::from_utf8(output.stdout).unwrap();
        assert!(stdout.starts_with("Usage: sim <COMMAND>\n"));
        let line = |command: &str| {
            stdout
                .lines()
                .find(|line| line.split_whitespace().next() == Some(command))
                .unwrap()
                .to_owned()
        };
        for command in ["lint", "run", "replay"] {
            assert!(!line(command).contains("not yet implemented"));
        }
        assert!(line("dry-run").contains("not yet implemented"));
    }
}

#[test]
fn unsupported_commands_and_bad_arguments_exit_two() {
    for args in [
        vec!["lint"],
        vec!["lint", "a.toml", "b.toml"],
        vec!["dry-run"],
        vec!["run"],
        vec!["run", "x.toml", "--out", "d"],
        vec!["run", "--offline", "x.toml"],
        vec!["run", "--offline", "x.toml", "--out", "d", "--threads", "0"],
        vec!["replay"],
        vec!["replay", "a.json", "b.json"],
        vec!["unknown"],
        vec!["--unknown"],
        vec!["--version", "extra"],
    ] {
        let output = sim(&args);
        assert_eq!(output.status.code(), Some(2), "{args:?}");
        assert!(output.stdout.is_empty());
        let stderr = String::from_utf8(output.stderr).unwrap();
        assert!(stderr.contains(args[0]));
        assert!(stderr.contains("sim --help"));
    }
}

#[test]
fn minimal_example_lints_clean() {
    let path = workspace_file("examples/minimal.toml");
    let path = path.to_str().unwrap();
    let output = sim(&["lint", path]);
    assert_eq!(
        output.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output.stderr.is_empty());
    assert_eq!(
        String::from_utf8(output.stdout).unwrap(),
        format!("{path}: ok\n")
    );
}

#[test]
fn unreadable_file_exits_one() {
    let output = sim(&["lint", "/nonexistent/scenario.toml"]);
    assert_eq!(output.status.code(), Some(1));
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.starts_with("/nonexistent/scenario.toml: error: cannot read the scenario"));
}

/// Each `tests/lint/*.toml` starts with `# expect: LINE:COL: FRAGMENT`; the first
/// diagnostic must be at that location and contain the fragment.
#[test]
fn negative_fixtures_report_file_line_column() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/lint");
    let mut fixtures: Vec<_> = fs::read_dir(&dir)
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .filter(|path| path.extension().is_some_and(|ext| ext == "toml"))
        .collect();
    fixtures.sort();
    assert!(fixtures.len() >= 15, "found {}", fixtures.len());
    for fixture in fixtures {
        let text = fs::read_to_string(&fixture).unwrap();
        let expect = text
            .lines()
            .next()
            .and_then(|line| line.strip_prefix("# expect: "))
            .unwrap_or_else(|| panic!("{} has no expect header", fixture.display()));
        let (location, fragment) = expect.split_once(": ").unwrap();
        let path = fixture.to_str().unwrap();
        let output = sim(&["lint", path]);
        let stderr = String::from_utf8(output.stderr).unwrap();
        assert_eq!(output.status.code(), Some(1), "{path}\n{stderr}");
        assert!(output.stdout.is_empty());
        let first = stderr.lines().next().unwrap();
        let prefix = format!("{path}:{location}: error: ");
        assert!(
            first.starts_with(&prefix) && first.contains(fragment),
            "{path}: expected {location}: …{fragment}…\ngot:\n{stderr}"
        );
        // The caret line sits under the reported column.
        let caret = stderr.lines().nth(2).unwrap();
        let column: usize = location.split(':').nth(1).unwrap().parse().unwrap();
        let bar = caret.find('|').unwrap();
        assert_eq!(
            caret[bar + 2..].chars().take_while(|&c| c == ' ').count() + 1,
            column,
            "{stderr}"
        );
        assert!(stderr.trim_end().ends_with("error") || stderr.trim_end().ends_with("errors"));
    }
}

#[test]
fn run_offline_then_replay_succeeds_and_rejects_reuse_of_the_output() {
    let dir = std::env::temp_dir().join(format!("iqsim-cli-run-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    let scenario = workspace_file("examples/minimal.toml");
    let out = dir.join("out");
    let run = sim(&[
        "run",
        "--offline",
        scenario.to_str().unwrap(),
        "--out",
        out.to_str().unwrap(),
        "--threads",
        "2",
    ]);
    let stdout = String::from_utf8(run.stdout).unwrap();
    assert_eq!(
        run.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&run.stderr)
    );
    assert!(stdout.contains("consumer/recording/rx.ch0.sigmf-data: 100000 samples"));
    assert!(stdout.contains("consumer/port_recording/port.sigmf-data: 100000 samples"));
    let manifest = out.join("run/manifest.json");
    let replay = sim(&["replay", manifest.to_str().unwrap()]);
    assert_eq!(
        replay.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&replay.stderr)
    );
    assert!(
        String::from_utf8(replay.stdout)
            .unwrap()
            .contains("replay ok, 2 outputs match")
    );
    // A second run into the same directory is refused.
    let again = sim(&[
        "run",
        "--offline",
        scenario.to_str().unwrap(),
        "--out",
        out.to_str().unwrap(),
    ]);
    assert_eq!(again.status.code(), Some(1));
    assert!(
        String::from_utf8(again.stderr)
            .unwrap()
            .contains("not empty")
    );
    let missing = sim(&["replay", dir.join("nope.json").to_str().unwrap()]);
    assert_eq!(missing.status.code(), Some(1));
    fs::remove_dir_all(&dir).unwrap();
}
