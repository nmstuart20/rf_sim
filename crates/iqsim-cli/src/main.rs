//! Command-line entry point for the simulator.
//! Per the plan section "Interfaces" (CLI), this binary will expose scenario
//! linting, capacity dry runs, execution, and verified replay.
//! `sim lint` is implemented; the other commands are not yet.
//!
//! Dependencies: iqsim-runtime, iqsim-control, iqsim-hil, iqsim-scenario, iqsim-contracts.

#![forbid(unsafe_code)]

use std::env;
use std::ffi::OsString;
use std::fs;
use std::process::ExitCode;

const USAGE: &str = "Usage: sim <COMMAND>

Commands:
  lint <FILE>  Validate a scenario file
  dry-run      Check scenario capacity (not yet implemented)
  run          Execute a scenario (not yet implemented)
  replay       Replay a verified run (not yet implemented)

Options:
  -h, --help       Print help
  -V, --version    Print version

Exit status: 0 success, 1 lint errors or an unreadable file, 2 usage error.";

/// `sim lint <FILE>`: diagnostics go to stderr as `file:line:col: error: ...`.
fn lint(path: &OsString) -> ExitCode {
    let shown = path.to_string_lossy();
    let source = match fs::read_to_string(path) {
        Ok(source) => source,
        Err(error) => {
            eprintln!("{shown}: error: cannot read the scenario: {error}");
            return ExitCode::FAILURE;
        }
    };
    let diagnostics = iqsim_scenario::lint(&source);
    if diagnostics.is_empty() {
        println!("{shown}: ok");
        return ExitCode::SUCCESS;
    }
    for diagnostic in &diagnostics {
        eprint!("{}", diagnostic.render(&shown, &source));
    }
    let count = diagnostics.len();
    let plural = if count == 1 { "" } else { "s" };
    eprintln!("{shown}: {count} error{plural}");
    ExitCode::FAILURE
}

fn usage_error(message: &str) -> ExitCode {
    eprintln!("sim: {message}\nRun `sim --help` for usage.");
    ExitCode::from(2)
}

fn main() -> ExitCode {
    let args: Vec<_> = env::args_os().skip(1).collect();
    match args.as_slice() {
        [] => println!("{USAGE}"),
        [arg] if arg == "--help" || arg == "-h" => println!("{USAGE}"),
        [arg] if arg == "--version" || arg == "-V" => {
            println!("sim {}", env!("CARGO_PKG_VERSION"));
        }
        [command, file] if command == "lint" => return lint(file),
        [command] if command == "lint" => return usage_error("lint: missing <FILE>"),
        [command, ..] if command == "lint" => return usage_error("lint: takes exactly one <FILE>"),
        _ => {
            return usage_error(&format!(
                "unsupported command or arguments: {}",
                args[0].to_string_lossy()
            ));
        }
    }
    ExitCode::SUCCESS
}
