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
//! Prompt 2 added the numeric core every render backend shares, written so the CUDA
//! backend can port it line by line:
//!
//! - [`wide`]: 128-bit wrapping arithmetic from 64-bit primitives.
//! - [`Delay`] (Q32.32 samples) and [`Phase`] (u64 cycle fraction), with the
//!   ties-to-even float → integer rules ([`f64_to_i64_ties_even`]).
//! - [`BlockPoly`] ([`DelayPoly`], [`PhasePoly`]): Newton-form cubic block polynomials
//!   with exact jump and three-add stepping, up to [`MAX_BLOCK_LEN`] samples.
//! - [`sincos`] over [`SIN_TABLE`].
//! - [`OffsetFrequency::phase_at`]: carrier-offset phase from exact time.
//!
//! See `docs/contracts.md` §1, §3 and §10 for the contracts these types implement.
//!
//! Dependencies: no iqsim crates.

#![forbid(unsafe_code)]

mod block;
mod clock;
mod error;
mod fixed;
mod offset;
mod poly;
mod rate;
mod time;
mod timebase;
mod trig;
pub mod wide;

pub use block::{Block, BlockSpec, Blocks};
pub use clock::{
    AcqInstant, ClockMap, ClockMapId, IdealClockMap, IdealTimestampClock, ReportedTime,
    TimestampClock,
};
pub use error::{NumericError, TimeError};
pub use fixed::{
    DELAY_FRAC_BITS, Delay, POLY_EXTRA_BITS, Phase, f32_to_i64_ties_even, f64_to_i64_ties_even,
    require_finite_f32, require_finite_f64,
};
pub use offset::OffsetFrequency;
pub use poly::{
    BlockPoly, DELAY_RANGE_MARGIN_LSB, DelayPoly, MAX_BLOCK_LEN, PhasePoly, PolyCursor, PolyOutput,
    binomials,
};
pub use rate::{MAX_RATE_DEN, MAX_RATE_HZ, MAX_RATE_NUM, SampleRate};
pub use time::{
    PS_PER_MICROSECOND, PS_PER_MILLISECOND, PS_PER_NANOSECOND, PS_PER_SECOND, TruthDuration,
    TruthTime,
};
pub use timebase::{ChannelTimebase, ExactTime, SampleIndex};
pub use trig::{SIN_TABLE, SIN_TABLE_BITS, sincos};
