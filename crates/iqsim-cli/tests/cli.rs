//! Verify the bootstrap CLI through its executable interface.

#![forbid(unsafe_code)]

use std::process::{Command, Output};

fn sim(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_sim"))
        .args(args)
        .output()
        .expect("run sim")
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
fn help_and_no_arguments_list_planned_commands() {
    for args in [vec![], vec!["--help"], vec!["-h"]] {
        let output = sim(&args);
        assert!(output.status.success());
        assert!(output.stderr.is_empty());
        let stdout = String::from_utf8(output.stdout).unwrap();
        assert!(stdout.starts_with("Usage: sim <COMMAND>\n"));
        for command in ["lint", "dry-run", "run", "replay"] {
            let line = stdout
                .lines()
                .find(|line| line.split_whitespace().next() == Some(command))
                .unwrap();
            assert!(line.contains("not yet implemented"));
        }
    }
}

#[test]
fn unsupported_commands_and_extra_arguments_exit_two() {
    for args in [
        vec!["lint"],
        vec!["dry-run"],
        vec!["run"],
        vec!["replay"],
        vec!["unknown"],
        vec!["--unknown"],
        vec!["--version", "extra"],
    ] {
        let output = sim(&args);
        assert_eq!(output.status.code(), Some(2));
        assert!(output.stdout.is_empty());
        let stderr = String::from_utf8(output.stderr).unwrap();
        assert!(stderr.contains(args[0]));
        assert!(stderr.contains("sim --help"));
    }
}
