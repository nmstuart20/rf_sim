use std::fmt;
use std::io;
use std::path::{Path, PathBuf};

/// Invalid recording input or a checked I/O operation that failed.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum IoError {
    /// Input was rejected before any bytes from that operation were written.
    InvalidInput {
        /// Associated file, if the operation targets a path.
        path: Option<PathBuf>,
        /// Explanation of the violated input rule.
        message: String,
    },
    /// An I/O operation or on-disk consistency check failed.
    Io {
        /// Associated file, if the operation targets a path.
        path: Option<PathBuf>,
        /// Original I/O error kind, or `InvalidData` for inconsistent content.
        kind: io::ErrorKind,
        /// Operation context and the underlying error message.
        message: String,
    },
}

impl IoError {
    /// Associated path, when the failure concerns a file.
    pub fn path(&self) -> Option<&Path> {
        match self {
            Self::InvalidInput { path, .. } | Self::Io { path, .. } => path.as_deref(),
        }
    }

    /// Standard I/O error classification (`InvalidInput` for validation errors).
    pub fn kind(&self) -> io::ErrorKind {
        match self {
            Self::InvalidInput { .. } => io::ErrorKind::InvalidInput,
            Self::Io { kind, .. } => *kind,
        }
    }

    pub(crate) fn invalid(message: impl Into<String>) -> Self {
        Self::InvalidInput {
            path: None,
            message: message.into(),
        }
    }

    pub(crate) fn io(error: io::Error) -> Self {
        Self::Io {
            path: None,
            kind: error.kind(),
            message: error.to_string(),
        }
    }

    pub(crate) fn at_path(mut self, file: &Path) -> Self {
        match &mut self {
            Self::InvalidInput { path, .. } | Self::Io { path, .. } => {
                *path = Some(file.to_owned());
            }
        }
        self
    }
}

impl fmt::Display for IoError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Some(path) = self.path() {
            write!(f, "{}: ", path.display())?;
        }
        match self {
            Self::InvalidInput { message, .. } | Self::Io { message, .. } => {
                write!(f, "{}: {message}", self.kind())
            }
        }
    }
}

impl std::error::Error for IoError {}
