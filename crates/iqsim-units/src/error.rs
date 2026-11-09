//! Parse errors for quantities.

use std::fmt;

use crate::mgrs::MgrsError;

/// Describes one kind of quantity for error messages and JSON Schema text.
#[derive(Debug, PartialEq, Eq)]
pub struct QuantityKind {
    /// Human name, e.g. "frequency".
    pub name: &'static str,
    /// Accepted unit spellings, case-sensitive.
    pub units: &'static [&'static str],
    /// A valid example value.
    pub example: &'static str,
}

impl QuantityKind {
    /// The accepted units as a comma-separated list.
    pub fn unit_list(&self) -> String {
        self.units.join(", ")
    }
}

/// Why a quantity string was rejected. Every message names the quantity, and most
/// suggest a valid spelling.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum UnitError {
    /// The value was empty or whitespace.
    Empty {
        /// The expected quantity.
        kind: &'static QuantityKind,
    },
    /// A bare number appeared where a unit is required.
    BareNumber {
        /// The expected quantity.
        kind: &'static QuantityKind,
        /// The number as written (or as TOML gave it).
        number: String,
    },
    /// The text does not start with a number.
    InvalidNumber {
        /// The expected quantity.
        kind: &'static QuantityKind,
        /// The rejected text.
        text: String,
    },
    /// A number with no unit after it.
    MissingUnit {
        /// The expected quantity.
        kind: &'static QuantityKind,
        /// The rejected text.
        text: String,
    },
    /// A unit that this quantity does not accept.
    UnknownUnit {
        /// The expected quantity.
        kind: &'static QuantityKind,
        /// The unit as written.
        unit: String,
    },
    /// The value cannot be represented exactly at the quantity's resolution.
    Inexact {
        /// The expected quantity.
        kind: &'static QuantityKind,
        /// The rejected text.
        text: String,
        /// The finest representable step, e.g. "1 µHz".
        resolution: &'static str,
    },
    /// The value is outside the representable or physical range.
    OutOfRange {
        /// The expected quantity.
        kind: &'static QuantityKind,
        /// The rejected text.
        text: String,
        /// What the range is.
        limit: String,
    },
    /// Any other malformed value (DMS layout, altitude reference, ...).
    Invalid {
        /// The expected quantity.
        kind: &'static QuantityKind,
        /// The rejected text.
        text: String,
        /// What is wrong and how to fix it.
        reason: String,
    },
    /// A malformed MGRS reference.
    Mgrs(MgrsError),
}

impl fmt::Display for UnitError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty { kind } => write!(
                f,
                "empty {}; write a value such as \"{}\"",
                kind.name, kind.example
            ),
            Self::BareNumber { kind, number } => write!(
                f,
                "{} needs a unit: got the bare number {number}; write it as a string such as \"{}\" (units: {})",
                kind.name,
                kind.example,
                kind.unit_list()
            ),
            Self::InvalidNumber { kind, text } => write!(
                f,
                "invalid {} \"{text}\": expected a number followed by a unit, such as \"{}\"",
                kind.name, kind.example
            ),
            Self::MissingUnit { kind, text } => write!(
                f,
                "{} \"{text}\" has no unit; add one of: {} (for example \"{}\")",
                kind.name,
                kind.unit_list(),
                kind.example
            ),
            Self::UnknownUnit { kind, unit } => write!(
                f,
                "unknown {} unit \"{unit}\"; expected one of: {} (units are case-sensitive)",
                kind.name,
                kind.unit_list()
            ),
            Self::Inexact {
                kind,
                text,
                resolution,
            } => write!(
                f,
                "{} \"{text}\" is finer than the {resolution} resolution",
                kind.name
            ),
            Self::OutOfRange { kind, text, limit } => {
                write!(f, "{} \"{text}\" is out of range: {limit}", kind.name)
            }
            Self::Invalid { kind, text, reason } => {
                write!(f, "invalid {} \"{text}\": {reason}", kind.name)
            }
            Self::Mgrs(error) => write!(f, "{error}"),
        }
    }
}

impl std::error::Error for UnitError {}

impl From<MgrsError> for UnitError {
    fn from(error: MgrsError) -> Self {
        Self::Mgrs(error)
    }
}
