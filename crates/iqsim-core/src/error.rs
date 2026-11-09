//! Error type for time, rate, index and clock operations.

use std::fmt;

/// Failure of a time, rate, sample-index, clock or block operation.
///
/// Every integer operation on time, index or rate is checked. Overflow is reported
/// through this type and never wraps or saturates (contracts §1.1).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum TimeError {
    /// A [`TruthTime`](crate::TruthTime) or [`TruthDuration`](crate::TruthDuration)
    /// result does not fit in i64 picoseconds.
    TimeOverflow,
    /// A [`SampleIndex`](crate::SampleIndex) result does not fit in u64.
    IndexOverflow,
    /// A sample rate numerator or denominator was zero.
    ZeroRate,
    /// A reduced sample rate exceeds the accepted bounds (contracts §10, decision 1).
    RateOutOfBounds {
        /// Reduced numerator.
        num: u64,
        /// Reduced denominator.
        den: u64,
    },
    /// The requested truth time precedes the channel epoch, so no sample index maps to it.
    BeforeEpoch,
    /// A block specification had zero length.
    ZeroBlockLength,
    /// A block range had `end < start`.
    InvalidRange,
}

impl fmt::Display for TimeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TimeOverflow => f.write_str("truth time overflow (i64 picoseconds)"),
            Self::IndexOverflow => f.write_str("sample index overflow (u64)"),
            Self::ZeroRate => f.write_str("sample rate numerator and denominator must be nonzero"),
            Self::RateOutOfBounds { num, den } => write!(
                f,
                "sample rate {num}/{den} S/s is outside the accepted bounds \
                 (denominator <= 2^32, numerator <= 2^48, rate <= 1 GS/s)"
            ),
            Self::BeforeEpoch => f.write_str("truth time precedes the channel epoch"),
            Self::ZeroBlockLength => f.write_str("block length must be nonzero"),
            Self::InvalidRange => f.write_str("block range end precedes its start"),
        }
    }
}

impl std::error::Error for TimeError {}

/// Failure of a fixed-point, polynomial or float → integer conversion.
///
/// Wrapping is part of the definition only for carrier phase and polynomial registers
/// (contracts §1.1). Every other range violation is reported here.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum NumericError {
    /// A floating-point input was NaN or ±Inf (contracts §1.4).
    NonFinite,
    /// A value does not fit its fixed-point or integer destination after rounding.
    OutOfRange,
    /// A block polynomial length was zero or above [`MAX_BLOCK_LEN`](crate::MAX_BLOCK_LEN).
    BlockLength {
        /// The rejected length.
        len: u64,
    },
    /// A local sample index was not inside the polynomial's block.
    IndexOutsideBlock {
        /// The requested local index.
        k: u64,
        /// The block length.
        len: u64,
    },
}

impl fmt::Display for NumericError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonFinite => f.write_str("floating-point value is NaN or infinite"),
            Self::OutOfRange => f.write_str("value is outside the destination range"),
            Self::BlockLength { len } => write!(
                f,
                "block polynomial length {len} is outside 1..={}",
                crate::MAX_BLOCK_LEN
            ),
            Self::IndexOutsideBlock { k, len } => {
                write!(f, "local index {k} is outside the block [0, {len})")
            }
        }
    }
}

impl std::error::Error for NumericError {}
