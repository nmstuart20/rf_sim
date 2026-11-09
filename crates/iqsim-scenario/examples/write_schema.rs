//! Writes the scenario JSON Schema to the path given as the only argument.
//! Run through `cargo xtask schema`, which writes `schema/scenario.schema.json`.

#![forbid(unsafe_code)]

use std::process::ExitCode;

fn main() -> ExitCode {
    let args: Vec<_> = std::env::args_os().skip(1).collect();
    let [path] = args.as_slice() else {
        eprintln!("usage: write_schema <OUTPUT.json>");
        return ExitCode::from(2);
    };
    match std::fs::write(path, iqsim_scenario::scenario_schema_json()) {
        Ok(()) => {
            println!("wrote {}", path.to_string_lossy());
            ExitCode::SUCCESS
        }
        Err(error) => {
            eprintln!("cannot write {}: {error}", path.to_string_lossy());
            ExitCode::FAILURE
        }
    }
}
