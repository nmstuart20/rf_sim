//! Records the build identity for run manifests: the compiler release, target, profile,
//! and a fingerprint of the workspace sources. It only reads files and runs the
//! configured `rustc -vV`; it never touches the network.
//!
//! The source fingerprint is SHA-256 over, in sorted path order, each file's
//! workspace-relative path, a NUL, its length (u64 LE) and its bytes, for the root
//! `Cargo.toml`, `rust-toolchain.toml`, every `crates/*/Cargo.toml`, `crates/*/build.rs`
//! and every file under `crates/*/src`. Tests, docs and examples are excluded: they do
//! not change the binary. Any source edit therefore changes the build identity, so
//! replay refuses runs made by different code (contracts §6).

use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

use sha2::{Digest, Sha256};

fn walk(dir: &Path, out: &mut Vec<PathBuf>) {
    let Ok(entries) = fs::read_dir(dir) else {
        return;
    };
    for entry in entries.flatten() {
        let path = entry.path();
        if path.is_dir() {
            walk(&path, out);
        } else {
            out.push(path);
        }
    }
}

fn source_fingerprint(root: &Path) -> String {
    let mut files = vec![root.join("Cargo.toml"), root.join("rust-toolchain.toml")];
    let mut crates: Vec<PathBuf> = fs::read_dir(root.join("crates"))
        .map(|d| {
            d.flatten()
                .map(|e| e.path())
                .filter(|p| p.is_dir())
                .collect()
        })
        .unwrap_or_default();
    crates.sort();
    for c in &crates {
        println!("cargo:rerun-if-changed={}", c.join("src").display());
        files.push(c.join("Cargo.toml"));
        files.push(c.join("build.rs"));
        walk(&c.join("src"), &mut files);
    }
    let mut rel: Vec<(String, PathBuf)> = files
        .into_iter()
        .filter(|p| p.is_file())
        .map(|p| {
            let r = p
                .strip_prefix(root)
                .unwrap_or(&p)
                .to_string_lossy()
                .replace('\\', "/");
            (r, p)
        })
        .collect();
    rel.sort();
    rel.dedup();
    let mut h = Sha256::new();
    for (r, p) in &rel {
        println!("cargo:rerun-if-changed={}", p.display());
        let bytes = fs::read(p).unwrap_or_default();
        h.update(r.as_bytes());
        h.update([0]);
        h.update((bytes.len() as u64).to_le_bytes());
        h.update(&bytes);
    }
    let digest = h.finalize();
    digest.iter().map(|b| format!("{b:02x}")).collect()
}

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-env-changed=RUSTC");
    let rustc = env::var("RUSTC").unwrap_or_else(|_| "rustc".to_owned());
    let release = Command::new(&rustc)
        .arg("-vV")
        .output()
        .ok()
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .and_then(|s| {
            s.lines()
                .find_map(|l| l.strip_prefix("release: ").map(str::to_owned))
        })
        .unwrap_or_else(|| "unknown".to_owned());
    println!("cargo:rustc-env=IQSIM_RUSTC_RELEASE={release}");
    println!(
        "cargo:rustc-env=IQSIM_TARGET={}",
        env::var("TARGET").unwrap_or_default()
    );
    println!(
        "cargo:rustc-env=IQSIM_PROFILE={}",
        env::var("PROFILE").unwrap_or_default()
    );
    let manifest_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap_or_default());
    let root = manifest_dir.join("../..");
    println!(
        "cargo:rustc-env=IQSIM_SOURCE_SHA256={}",
        source_fingerprint(&root)
    );
}
