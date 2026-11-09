//! Command-line entry point for the simulator.
//! Per the plan section "Interfaces" (CLI), this binary will expose scenario
//! linting, capacity dry runs, execution, and verified replay.
//! `sim lint`, `sim run --offline` and `sim replay` are implemented; `sim dry-run`
//! is not yet.
//!
//! Dependencies: iqsim-runtime, iqsim-control, iqsim-hil, iqsim-scenario, iqsim-contracts.

#![forbid(unsafe_code)]

use std::env;
use std::ffi::OsString;
use std::fs;
use std::process::ExitCode;

const USAGE: &str = "Usage: sim <COMMAND>

Commands:
  lint <FILE>                            Validate a scenario file
  dry-run                                Check scenario capacity (not yet implemented)
  run --offline <FILE> --out <DIR>       Render a scenario as fast as possible into DIR
  replay <MANIFEST>                      Verify a run's inputs, re-render it and compare

Options:
  -h, --help       Print help
  -V, --version    Print version
  --threads <N>    Render threads for run and replay (output does not depend on it)

Exit status: 0 success, 1 lint, run or replay failure, 2 usage error.";

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

/// Splits `--threads N` out of `args`.
fn take_threads(args: &mut Vec<OsString>) -> Result<iqsim_runtime::RunOptions, String> {
    let mut options = iqsim_runtime::RunOptions::default();
    if let Some(i) = args.iter().position(|a| a == "--threads") {
        let value = args
            .get(i + 1)
            .and_then(|v| v.to_str())
            .and_then(|v| v.parse::<usize>().ok())
            .filter(|&n| (1..=256).contains(&n))
            .ok_or("--threads needs a number from 1 to 256")?;
        options.threads = value;
        args.drain(i..i + 2);
    }
    Ok(options)
}

/// `sim run --offline <FILE> --out <DIR>`.
fn run(mut args: Vec<OsString>) -> ExitCode {
    let options = match take_threads(&mut args) {
        Ok(o) => o,
        Err(e) => return usage_error(&format!("run: {e}")),
    };
    let mut offline = false;
    let mut out = None;
    let mut file = None;
    let mut i = 0;
    while i < args.len() {
        let a = &args[i];
        if a == "--offline" {
            offline = true;
        } else if a == "--out" {
            let Some(dir) = args.get(i + 1) else {
                return usage_error("run: --out needs a directory");
            };
            out = Some(dir.clone());
            i += 1;
        } else if a.to_string_lossy().starts_with('-') || file.is_some() {
            return usage_error(&format!("run: unexpected argument {}", a.to_string_lossy()));
        } else {
            file = Some(a.clone());
        }
        i += 1;
    }
    if !offline {
        return usage_error("run: only --offline runs are implemented so far");
    }
    let (Some(file), Some(out)) = (file, out) else {
        return usage_error("run: needs a scenario <FILE> and --out <DIR>");
    };
    match iqsim_runtime::run_offline(file.as_ref(), out.as_ref(), options) {
        Ok(m) => {
            for w in &m.warnings {
                eprintln!("warning: {w}");
            }
            for o in &m.outputs {
                println!(
                    "{}: {} samples, sha256 {}",
                    o.path, o.sample_count, o.sha256
                );
            }
            println!(
                "{}: manifest {}",
                std::path::Path::new(&out)
                    .join("run/manifest.json")
                    .display(),
                m.deterministic_sha256()
            );
            ExitCode::SUCCESS
        }
        Err(e) => {
            eprintln!("sim run: error: {e}");
            ExitCode::FAILURE
        }
    }
}

/// `sim replay <MANIFEST>`.
fn replay(mut args: Vec<OsString>) -> ExitCode {
    let options = match take_threads(&mut args) {
        Ok(o) => o,
        Err(e) => return usage_error(&format!("replay: {e}")),
    };
    let [manifest] = args.as_slice() else {
        return usage_error("replay: takes exactly one <MANIFEST>");
    };
    match iqsim_runtime::replay(manifest.as_ref(), options) {
        Ok(r) => {
            println!(
                "{}: replay ok, {} outputs match (manifest {})",
                manifest.to_string_lossy(),
                r.outputs_matched,
                r.deterministic_sha256
            );
            ExitCode::SUCCESS
        }
        Err(e) => {
            eprintln!("sim replay: error: {e}");
            ExitCode::FAILURE
        }
    }
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
        [command, rest @ ..] if command == "run" => return run(rest.to_vec()),
        [command, rest @ ..] if command == "replay" => return replay(rest.to_vec()),
        _ => {
            return usage_error(&format!(
                "unsupported command or arguments: {}",
                args[0].to_string_lossy()
            ));
        }
    }
    ExitCode::SUCCESS
}
