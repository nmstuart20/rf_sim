//! Waveform generation and immutable signal tables.
//! Per the plan section "Emitters, behaviors and waveforms", this crate owns
//! analytic pulse envelopes, continuous native-rate signals, and IQ-file sources.
//! Tap-specific filtering and fractional-delay tables support deterministic rendering.
//!
//! Dependencies: iqsim-core, iqsim-units.

#![forbid(unsafe_code)]
