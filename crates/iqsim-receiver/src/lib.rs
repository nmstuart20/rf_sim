//! Virtual receiver processing and coherent channel state.
//! Per the plan section "Receiver model", this crate owns gain, shaped noise,
//! LO effects, IQ imbalance, decimation, AGC, and ADC conversion after complete sums.
//! Ordered per-channel state is applied once and never enters antenna-port output.
//!
//! Dependencies: iqsim-core, iqsim-units, iqsim-contracts.
//!
//! # Modules
//!
//! - [`chain`]: the virtual-receiver chain ([`ChainConfig`], [`ReceiverChannelState`]).
//!   It runs once per channel on the complete reduced path sum of the virtual-receiver
//!   target, in block order, and produces the consumer's samples.
//! - [`port`]: the antenna-port output conversion ([`PortConverter`]). The antenna-port
//!   target bypasses every virtual-receiver effect; this is its only conversion.
//!
//! The full stage definitions, number formats and state-ownership rules are in
//! `docs/contracts.md` §5.4 and §5.5.

#![forbid(unsafe_code)]

pub mod chain;
pub mod port;

pub use chain::{
    AdcConfig, AgcConfig, ChainConfig, ChainError, ChannelOutput, DecimationConfig, FrontEnd,
    IqImpairment, LoConfig, NoiseConfig, OutputFormat, PhaseNoise, ReceiverChannelState, Samples,
};
pub use port::{PortConverter, PortError};

/// Boltzmann constant k, J/K (exact, SI 2019).
pub const BOLTZMANN_J_PER_K: f64 = 1.380_649e-23;

/// Reference noise temperature T₀ for noise figure, K (IEEE standard 290 K).
pub const REFERENCE_TEMPERATURE_K: f64 = 290.0;
