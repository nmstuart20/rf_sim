//! Deterministic numeric and time foundations.
//! Per the plan section "Time, numerics and determinism", this crate owns
//! integer time, fixed-point delay and phase, deterministic sin/cos, and Philox RNG.
//! It also owns the ClockMap trait and shared commit-boundary time types.
//!
//! Implemented so far (prompt 1):
//!
//! - [`TruthTime`] and [`TruthDuration`]: i64 picoseconds, checked arithmetic.
//! - [`SampleRate`]: exact reduced rational with bounded numerator and denominator.
//! - [`SampleIndex`] and [`ChannelTimebase`]: exact index ↔ time conversions.
//! - [`ClockMap`] with [`IdealClockMap`], and the separate [`TimestampClock`].
//! - [`Block`], [`BlockSpec`] and the [`Blocks`] iterator.
//!
//! See `docs/contracts.md` §1.2, §3 and §10 for the contracts these types implement.
//!
//! Dependencies: no iqsim crates.

#![forbid(unsafe_code)]

mod block;
mod clock;
mod error;
mod rate;
mod time;
mod timebase;

pub use block::{Block, BlockSpec, Blocks};
pub use clock::{
    AcqInstant, ClockMap, ClockMapId, IdealClockMap, IdealTimestampClock, ReportedTime,
    TimestampClock,
};
pub use error::TimeError;
pub use rate::{MAX_RATE_DEN, MAX_RATE_HZ, MAX_RATE_NUM, SampleRate};
pub use time::{
    PS_PER_MICROSECOND, PS_PER_MILLISECOND, PS_PER_NANOSECOND, PS_PER_SECOND, TruthDuration,
    TruthTime,
};
pub use timebase::{ChannelTimebase, ExactTime, SampleIndex};
