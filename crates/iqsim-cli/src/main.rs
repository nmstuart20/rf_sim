//! Command-line entry point for the simulator.
//! Per the plan section "Interfaces" (CLI), this binary will expose scenario
//! linting, capacity dry runs, execution, and verified replay.
//! This bootstrap exposes help and version reporting only.
//!
//! Dependencies: iqsim-runtime, iqsim-control, iqsim-hil, iqsim-scenario, iqsim-contracts.

#![forbid(unsafe_code)]

use std::env;
use std::process::ExitCode;

const USAGE: &str = "Usage: sim <COMMAND>

Commands:
  lint      Validate a scenario (not yet implemented)
  dry-run   Check scenario capacity (not yet implemented)
  run       Execute a scenario (not yet implemented)
  replay    Replay a verified run (not yet implemented)

Options:
  -h, --help       Print help
  -V, --version    Print version";

fn main() -> ExitCode {
    let args: Vec<_> = env::args_os().skip(1).collect();
    match args.as_slice() {
        [] => println!("{USAGE}"),
        [arg] if arg == "--help" || arg == "-h" => println!("{USAGE}"),
        [arg] if arg == "--version" || arg == "-V" => {
            println!("sim {}", env!("CARGO_PKG_VERSION"));
        }
        _ => {
            eprintln!(
                "sim: unsupported command or arguments: {}\nRun `sim --help` for usage.",
                args[0].to_string_lossy()
            );
            return ExitCode::from(2);
        }
    }
    ExitCode::SUCCESS
}
