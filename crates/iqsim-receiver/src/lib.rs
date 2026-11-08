//! Virtual receiver processing and coherent channel state.
//! Per the plan section "Receiver model", this crate owns gain, shaped noise,
//! LO effects, IQ imbalance, decimation, AGC, and ADC conversion after complete sums.
//! Ordered per-channel state is applied once and never enters antenna-port output.
//!
//! Dependencies: iqsim-core, iqsim-units, iqsim-contracts.

#![forbid(unsafe_code)]
