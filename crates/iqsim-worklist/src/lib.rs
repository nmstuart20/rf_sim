//! Flat block work lists shared by all render backends.
//! Per the plan section "Architecture", this crate owns pointer-free path,
//! emission, and immutable-table references for each block and output tap.
//! Delay, phase, and retarded waveform coordinates preserve exact render timing.
//!
//! Dependencies: iqsim-core.

#![forbid(unsafe_code)]
