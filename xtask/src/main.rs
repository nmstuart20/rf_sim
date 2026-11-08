//! Standard-library-only workspace maintenance commands.
//! Implements offline CI, architectural dependency checks, golden blessing,
//! and the future GPU-test entry point from the plan's "Validation strategy".

#![forbid(unsafe_code)]

use std::env;
use std::error::Error;
use std::io::{self, Write};
use std::path::Path;
use std::process::{Command, ExitCode};

const USAGE: &str = "Usage: cargo xtask <COMMAND>

Commands:
  ci                  Check formatting, lints, tests, docs, and dependencies
  deps                Check architectural dependency boundaries
  bless [TEST ARGS]   Regenerate golden hashes and show changed goldens
  gpu-test            Placeholder for CUDA tests (prompt 20)
  help                Print help";

type TaskResult = Result<(), Box<dyn Error>>;

fn workspace_root() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("xtask must be a direct child of the workspace root")
}

fn cargo() -> Command {
    let mut command = Command::new(env::var_os("CARGO").unwrap_or_else(|| "cargo".into()));
    command
        .current_dir(workspace_root())
        .env("CARGO_NET_OFFLINE", "true");
    command
}

fn header(label: &str) -> io::Result<()> {
    println!("\n== {label} ==");
    io::stdout().flush()
}

fn run(label: &str, command: &mut Command) -> TaskResult {
    header(label)?;
    let status = command.status()?;
    if !status.success() {
        return Err(io::Error::other(format!("{label} failed: {status}")).into());
    }
    Ok(())
}

fn ci() -> TaskResult {
    run(
        "cargo fmt --all -- --check",
        cargo().args(["fmt", "--all", "--", "--check"]),
    )?;
    run(
        "cargo clippy --workspace --all-targets --offline --locked -- -D warnings",
        cargo().args([
            "clippy",
            "--workspace",
            "--all-targets",
            "--offline",
            "--locked",
            "--",
            "-D",
            "warnings",
        ]),
    )?;
    run(
        "cargo test --workspace --offline --locked",
        cargo().args(["test", "--workspace", "--offline", "--locked"]),
    )?;
    run(
        "cargo doc --workspace --no-deps --offline --locked (RUSTDOCFLAGS=-D warnings)",
        cargo()
            .args(["doc", "--workspace", "--no-deps", "--offline", "--locked"])
            .env("RUSTDOCFLAGS", "-D warnings"),
    )?;
    header("deps")?;
    deps()
}

fn tree(args: &[&str]) -> Result<String, Box<dyn Error>> {
    let mut command = cargo();
    command.args([
        "tree",
        "--offline",
        "--locked",
        "-e",
        "normal",
        "--prefix",
        "none",
    ]);
    command.args(args);
    let output = command.output()?;
    if !output.status.success() {
        return Err(io::Error::other(format!(
            "cargo tree {} failed: {}\n{}",
            args.join(" "),
            output.status,
            String::from_utf8_lossy(&output.stderr)
        ))
        .into());
    }
    Ok(String::from_utf8(output.stdout)?)
}

// With --prefix none, each line starts with the exact package name; paths,
// versions, and Cargo's repeated-dependency markers follow whitespace.
fn package_names(tree: &str) -> impl Iterator<Item = &str> {
    tree.lines()
        .filter_map(|line| line.split_whitespace().next())
}

fn reject_dependencies(owner: &str, tree: &str, forbidden: &[&str]) -> TaskResult {
    for dependency in package_names(tree) {
        if forbidden.contains(&dependency) {
            return Err(io::Error::other(format!(
                "{owner} must not depend on {dependency} (including transitively)"
            ))
            .into());
        }
    }
    Ok(())
}

fn check_hardware_dependents(hardware: &str, tree: &str) -> TaskResult {
    const ALLOWED: &[&str] = &["iqsim-render-cuda", "iqsim-hil", "iqsim-cli", "xtask"];
    for dependent in package_names(tree) {
        if !ALLOWED.contains(&dependent) {
            return Err(io::Error::other(format!(
                "{dependent} must not depend on {hardware} (including transitively)"
            ))
            .into());
        }
    }
    Ok(())
}

fn deps() -> TaskResult {
    for (owner, forbidden) in [
        ("iqsim-control", ["iqsim-runtime", "iqsim-hil"]),
        ("iqsim-runtime", ["iqsim-control", "iqsim-hil"]),
    ] {
        let output = tree(&["-p", owner])?;
        reject_dependencies(owner, &output, &forbidden)?;
    }
    for hardware in ["iqsim-render-cuda", "iqsim-hil"] {
        let output = tree(&["--workspace", "--invert", hardware])?;
        check_hardware_dependents(hardware, &output)?;
    }
    println!("Dependency boundaries passed (normal dependencies, including transitive edges).");
    Ok(())
}

fn bless(args: &[std::ffi::OsString]) -> TaskResult {
    run(
        "cargo test --workspace --offline --locked (IQSIM_BLESS=1)",
        cargo()
            .args(["test", "--workspace", "--offline", "--locked"])
            .args(args)
            .env("IQSIM_BLESS", "1"),
    )?;
    run(
        "git status --short tests/golden",
        Command::new("git").current_dir(workspace_root()).args([
            "status",
            "--short",
            "tests/golden",
        ]),
    )
}

fn main() -> ExitCode {
    let args: Vec<_> = env::args_os().skip(1).collect();
    let Some(command) = args.first() else {
        println!("{USAGE}");
        return ExitCode::SUCCESS;
    };
    let result = match command.to_str() {
        Some("ci") if args.len() == 1 => ci(),
        Some("deps") if args.len() == 1 => deps(),
        Some("bless") => bless(&args[1..]),
        Some("gpu-test") if args.len() == 1 => {
            println!("CUDA tests arrive in prompt 20.");
            Ok(())
        }
        Some("help" | "--help" | "-h") if args.len() == 1 => {
            println!("{USAGE}");
            Ok(())
        }
        _ => {
            eprintln!(
                "Unknown command or arguments: {}\n\n{USAGE}",
                command.to_string_lossy()
            );
            return ExitCode::from(2);
        }
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("xtask: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{check_hardware_dependents, package_names, reject_dependencies};

    #[test]
    fn tree_scan_uses_exact_names_and_handles_paths_with_spaces() {
        let tree = "iqsim-control v0.1.0 (/a path/crates/iqsim-control)\n\
                    iqsim-core v0.1.0 (/a path/crates/iqsim-core) (*)\n\
                    iqsim-runtime-extra v0.1.0\n";
        assert_eq!(
            package_names(tree).collect::<Vec<_>>(),
            ["iqsim-control", "iqsim-core", "iqsim-runtime-extra"]
        );
        reject_dependencies("iqsim-control", tree, &["iqsim-runtime"]).unwrap();
        let error = reject_dependencies(
            "iqsim-control",
            &format!("{tree}iqsim-runtime v0.1.0\n"),
            &["iqsim-runtime"],
        )
        .unwrap_err();
        assert!(
            error
                .to_string()
                .contains("iqsim-control must not depend on iqsim-runtime")
        );
    }

    #[test]
    fn hardware_reverse_tree_rejects_every_unapproved_dependent() {
        let allowed = "iqsim-hil v0.1.0\niqsim-cli v0.1.0\nxtask v0.1.0\n";
        check_hardware_dependents("iqsim-hil", allowed).unwrap();
        for hardware in ["iqsim-hil", "iqsim-render-cuda"] {
            for dependent in ["iqsim-runtime", "iqsim-control", "iqsim-testkit"] {
                let error = check_hardware_dependents(
                    hardware,
                    &format!("{hardware} v0.1.0\n{dependent} v0.1.0\n"),
                )
                .unwrap_err();
                assert!(
                    error
                        .to_string()
                        .contains(&format!("{dependent} must not depend on {hardware}"))
                );
            }
        }
    }
}
