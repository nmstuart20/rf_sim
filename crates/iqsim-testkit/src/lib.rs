//! Golden-hash support for deterministic validation.
//! Per the plan section "Validation strategy", this crate compares canonical bytes
//! against checked-in SHA-256 digests and explicitly blesses intentional changes.
//! Other iqsim crates may use this support only as a dev-dependency.
//!
//! Dependencies: no iqsim crates.

#![forbid(unsafe_code)]

use std::env;
use std::error::Error;
use std::ffi::OsStr;
use std::fmt;
use std::fs;
use std::io;
use std::path::{Path, PathBuf};

use sha2::{Digest, Sha256};

const BLESS_HINT: &str =
    "if this change is intentional, run `cargo xtask bless` and commit the updated golden file";

/// Select whether a golden digest is verified or deliberately updated.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Mode {
    /// Compare the bytes' digest against an existing, well-formed golden file.
    Compare,
    /// Create or replace the golden file when its content differs.
    Bless,
}

/// Successful result of checking or blessing a golden digest.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Outcome {
    /// The digest matched the existing golden file in comparison mode.
    Matched,
    /// Blessing created or changed the golden file.
    Written,
    /// Blessing found identical content and performed no write.
    Unchanged,
}

/// A validation, comparison, or filesystem failure while checking a golden.
#[derive(Debug)]
pub enum GoldenError {
    /// The name cannot safely identify a portable, relative golden path.
    InvalidName {
        /// The rejected name.
        name: String,
        /// The violated naming rule.
        reason: &'static str,
    },
    /// Comparison requires a golden file that does not exist.
    Missing {
        /// The validated golden name.
        name: String,
        /// The missing file's path.
        path: PathBuf,
    },
    /// The file is not exactly 64 lowercase hexadecimal digits and one newline.
    Malformed {
        /// The validated golden name.
        name: String,
        /// The malformed file's path.
        path: PathBuf,
    },
    /// The bytes' SHA-256 digest differs from the recorded digest.
    Mismatch {
        /// The validated golden name.
        name: String,
        /// The golden file's path.
        path: PathBuf,
        /// The digest recorded in the golden file.
        expected: String,
        /// The digest computed from the supplied bytes.
        actual: String,
        /// The number of supplied bytes.
        byte_len: usize,
    },
    /// A filesystem operation failed for a reason other than a missing golden.
    Io {
        /// The path involved in the failed operation.
        path: PathBuf,
        /// The underlying filesystem error.
        source: io::Error,
    },
}

impl fmt::Display for GoldenError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidName { name, reason } => {
                write!(f, "invalid golden name {name:?}: {reason}")
            }
            Self::Missing { name, path } => write!(
                f,
                "golden file does not exist\nname: {name}\nfile: {}\nrun `cargo xtask bless` to create it",
                path.display()
            ),
            Self::Malformed { name, path } => write!(
                f,
                "malformed golden file\nname: {name}\nfile: {}\nexpected exactly 64 lowercase hexadecimal digits followed by a single newline",
                path.display()
            ),
            Self::Mismatch {
                name,
                path,
                expected,
                actual,
                byte_len,
            } => write!(
                f,
                "golden digest mismatch\nname: {name}\nfile: {}\nexpected digest: {expected}\nactual digest: {actual}\nbyte length: {byte_len}\n{BLESS_HINT}",
                path.display()
            ),
            Self::Io { path, source } => {
                write!(f, "golden file I/O error at {}: {source}", path.display())
            }
        }
    }
}

impl Error for GoldenError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Io { source, .. } => Some(source),
            _ => None,
        }
    }
}

fn validate_name(name: &str) -> Result<(), GoldenError> {
    let invalid = |reason| GoldenError::InvalidName {
        name: name.to_owned(),
        reason,
    };
    if name.is_empty() {
        return Err(invalid("name must not be empty"));
    }
    for segment in name.split('/') {
        if segment.is_empty() {
            return Err(invalid(
                "name must not have leading/trailing '/' or empty path segments",
            ));
        }
        if segment == "." || segment == ".." {
            return Err(invalid("'.' and '..' path segments are forbidden"));
        }
        if !segment
            .bytes()
            .all(|byte| byte.is_ascii_lowercase() || byte.is_ascii_digit() || b"_-".contains(&byte))
        {
            return Err(invalid(
                "each '/'-separated segment may contain only [a-z0-9_-]",
            ));
        }
    }
    Ok(())
}

fn io_error(path: &Path, source: io::Error) -> GoldenError {
    GoldenError::Io {
        path: path.to_owned(),
        source,
    }
}

/// Check or bless `<dir>/<name>.sha256` using the SHA-256 digest of `bytes`.
///
/// `dir` is the golden directory itself. Names must consist of nonempty
/// slash-separated segments containing only lowercase ASCII letters, digits,
/// underscores, and hyphens. Absolute paths and dot segments are rejected.
/// Blessing creates parent directories and writes only when content differs.
///
/// # Errors
///
/// Returns [`GoldenError`] for an invalid name, missing or malformed comparison
/// file, digest mismatch, or filesystem failure. A valid file contains exactly
/// 64 lowercase hexadecimal digits followed by one newline.
pub fn check_golden_in(
    dir: &Path,
    name: &str,
    bytes: &[u8],
    mode: Mode,
) -> Result<Outcome, GoldenError> {
    validate_name(name)?;
    let path = dir.join(format!("{name}.sha256"));
    let actual = hex::encode(Sha256::digest(bytes));
    let content = format!("{actual}\n");
    let existing = match fs::read(&path) {
        Ok(existing) => Some(existing),
        Err(error) if error.kind() == io::ErrorKind::NotFound => None,
        Err(error) => return Err(io_error(&path, error)),
    };

    if mode == Mode::Bless {
        if existing.as_deref() == Some(content.as_bytes()) {
            return Ok(Outcome::Unchanged);
        }
        let parent = path.parent().expect("golden path always has a parent");
        fs::create_dir_all(parent).map_err(|error| io_error(parent, error))?;
        fs::write(&path, content).map_err(|error| io_error(&path, error))?;
        return Ok(Outcome::Written);
    }

    let existing = existing.ok_or_else(|| GoldenError::Missing {
        name: name.to_owned(),
        path: path.clone(),
    })?;
    if existing.len() != 65
        || existing[64] != b'\n'
        || !existing[..64]
            .iter()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(byte))
    {
        return Err(GoldenError::Malformed {
            name: name.to_owned(),
            path,
        });
    }
    let expected = String::from_utf8(existing[..64].to_vec())
        .expect("validated hexadecimal digits are valid UTF-8");
    if expected != actual {
        return Err(GoldenError::Mismatch {
            name: name.to_owned(),
            path,
            expected,
            actual,
            byte_len: bytes.len(),
        });
    }
    Ok(Outcome::Matched)
}

/// Assert that `bytes` match `<workspace root>/tests/golden/<name>.sha256`.
///
/// Set `IQSIM_BLESS=1`, normally through `cargo xtask bless`, to create or update
/// the golden instead. The workspace root is derived from this crate's
/// compile-time manifest directory, so the caller's working directory is irrelevant.
///
/// # Panics
///
/// Panics if the name is invalid, the file cannot be read or written, or the
/// recorded digest is missing, malformed, or different in comparison mode.
pub fn assert_golden(name: &str, bytes: &[u8]) {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(Path::parent)
        .expect("iqsim-testkit must be under <workspace root>/crates");
    let mode = if env::var_os("IQSIM_BLESS").as_deref() == Some(OsStr::new("1")) {
        Mode::Bless
    } else {
        Mode::Compare
    };
    check_golden_in(&root.join("tests/golden"), name, bytes, mode)
        .unwrap_or_else(|error| panic!("{error}"));
    if mode == Mode::Bless {
        eprintln!("blessed {name}: {}", hex::encode(Sha256::digest(bytes)));
    }
}

#[cfg(test)]
mod tests {
    use std::sync::atomic::{AtomicU64, Ordering};

    use super::*;

    const ABC_SHA256: &str = "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad";
    static NEXT_TEMP_DIR: AtomicU64 = AtomicU64::new(0);

    struct TempDir(PathBuf);

    impl TempDir {
        fn new() -> Self {
            loop {
                let sequence = NEXT_TEMP_DIR.fetch_add(1, Ordering::Relaxed);
                let path = env::temp_dir()
                    .join(format!("iqsim-testkit-{}-{sequence}", std::process::id()));
                match fs::create_dir(&path) {
                    Ok(()) => return Self(path),
                    Err(error) if error.kind() == io::ErrorKind::AlreadyExists => continue,
                    Err(error) => panic!("create test directory {}: {error}", path.display()),
                }
            }
        }
    }

    impl Drop for TempDir {
        fn drop(&mut self) {
            fs::remove_dir_all(&self.0).expect("remove test directory");
        }
    }

    #[test]
    fn sha256_known_answer() {
        assert_eq!(hex::encode(Sha256::digest(b"abc")), ABC_SHA256);
    }

    #[test]
    fn comparison_matches_a_known_digest() {
        let dir = TempDir::new();
        fs::write(dir.0.join("abc.sha256"), format!("{ABC_SHA256}\n")).unwrap();
        assert_eq!(
            check_golden_in(&dir.0, "abc", b"abc", Mode::Compare).unwrap(),
            Outcome::Matched
        );
    }

    #[test]
    fn mismatch_reports_digests_name_path_length_and_bless_hint() {
        let dir = TempDir::new();
        let path = dir.0.join("abc.sha256");
        fs::write(&path, format!("{ABC_SHA256}\n")).unwrap();
        let error = check_golden_in(&dir.0, "abc", b"changed", Mode::Compare).unwrap_err();
        let message = error.to_string();
        assert!(matches!(error, GoldenError::Mismatch { byte_len: 7, .. }));
        for expected in [
            "name: abc",
            path.to_str().unwrap(),
            ABC_SHA256,
            &hex::encode(Sha256::digest(b"changed")),
            "byte length: 7",
            BLESS_HINT,
        ] {
            assert!(
                message.contains(expected),
                "missing {expected:?} in {message}"
            );
        }
        assert_eq!(fs::read_to_string(path).unwrap(), format!("{ABC_SHA256}\n"));
    }

    #[test]
    fn missing_file_explains_how_to_bless() {
        let dir = TempDir::new();
        let error = check_golden_in(&dir.0, "missing", b"abc", Mode::Compare).unwrap_err();
        assert!(matches!(error, GoldenError::Missing { .. }));
        let message = error.to_string();
        assert!(message.contains("golden file does not exist"));
        assert!(message.contains("run `cargo xtask bless`"));
        assert!(!dir.0.join("missing.sha256").exists());
    }

    #[test]
    fn malformed_files_are_rejected() {
        let dir = TempDir::new();
        let cases = [
            Vec::new(),
            b"not a digest\n".to_vec(),
            ABC_SHA256.as_bytes().to_vec(),
            format!("{ABC_SHA256}\n\n").into_bytes(),
            format!("{ABC_SHA256}\r\n").into_bytes(),
            format!("{}\n", ABC_SHA256.to_uppercase()).into_bytes(),
            format!("g{}\n", &ABC_SHA256[1..]).into_bytes(),
            [vec![0xff; 64], vec![b'\n']].concat(),
            format!(" {ABC_SHA256}\n").into_bytes(),
        ];
        for content in cases {
            fs::write(dir.0.join("bad.sha256"), content).unwrap();
            let error = check_golden_in(&dir.0, "bad", b"abc", Mode::Compare).unwrap_err();
            assert!(matches!(error, GoldenError::Malformed { .. }));
            assert!(error.to_string().contains("malformed golden file"));
        }
    }

    #[test]
    fn bless_creates_nested_file_then_comparison_matches() {
        let dir = TempDir::new();
        assert_eq!(
            check_golden_in(&dir.0, "nested/abc", b"abc", Mode::Bless).unwrap(),
            Outcome::Written
        );
        assert_eq!(
            fs::read(dir.0.join("nested/abc.sha256")).unwrap(),
            format!("{ABC_SHA256}\n").as_bytes()
        );
        assert_eq!(
            check_golden_in(&dir.0, "nested/abc", b"abc", Mode::Compare).unwrap(),
            Outcome::Matched
        );
    }

    #[test]
    fn bless_unchanged_performs_no_write() {
        let dir = TempDir::new();
        let path = dir.0.join("abc.sha256");
        check_golden_in(&dir.0, "abc", b"abc", Mode::Bless).unwrap();
        let modified = fs::metadata(&path).unwrap().modified().unwrap();
        assert_eq!(
            check_golden_in(&dir.0, "abc", b"abc", Mode::Bless).unwrap(),
            Outcome::Unchanged
        );
        assert_eq!(fs::metadata(path).unwrap().modified().unwrap(), modified);
    }

    #[test]
    fn bless_replaces_malformed_and_different_digests() {
        let dir = TempDir::new();
        let path = dir.0.join("abc.sha256");
        for previous in ["malformed\n".to_owned(), format!("{}\n", "0".repeat(64))] {
            fs::write(&path, previous).unwrap();
            assert_eq!(
                check_golden_in(&dir.0, "abc", b"abc", Mode::Bless).unwrap(),
                Outcome::Written
            );
            assert_eq!(
                fs::read_to_string(&path).unwrap(),
                format!("{ABC_SHA256}\n")
            );
        }
    }

    #[test]
    fn names_are_validated_before_io_in_both_modes() {
        let dir = TempDir::new();
        for name in [
            "", "/abc", "abc/", "abc//def", ".", "..", "a/./b", "a/../b", "ABC", "a.b", "a b",
            "a\\b", "é", "a\nb", "a\0b", "C:/abc",
        ] {
            for mode in [Mode::Compare, Mode::Bless] {
                let error = check_golden_in(&dir.0, name, b"abc", mode).unwrap_err();
                assert!(matches!(error, GoldenError::InvalidName { .. }), "{name:?}");
                assert!(error.to_string().contains("invalid golden name"));
            }
        }
        assert_eq!(fs::read_dir(&dir.0).unwrap().count(), 0);
        for name in ["abc", "a0_-", "nested/self_check-2", "0/_/-"] {
            check_golden_in(&dir.0, name, b"abc", Mode::Bless).unwrap();
            check_golden_in(&dir.0, name, b"abc", Mode::Compare).unwrap();
        }
    }

    #[test]
    fn filesystem_failures_retain_their_source() {
        let dir = TempDir::new();
        let blocked = dir.0.join("file");
        fs::write(&blocked, b"not a directory").unwrap();
        let error = check_golden_in(&blocked, "abc", b"abc", Mode::Bless).unwrap_err();
        assert!(matches!(error, GoldenError::Io { .. }));
        assert!(error.source().is_some());
        assert!(error.to_string().contains("golden file I/O error"));
    }
}
