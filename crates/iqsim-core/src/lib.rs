//! Deterministic numeric and time foundations.
//! Per the plan section "Time, numerics and determinism", this crate owns
//! integer time, fixed-point delay and phase, deterministic sin/cos, and Philox RNG.
//! It also owns the ClockMap trait and shared commit-boundary time types.
//!
//! Dependencies: no iqsim crates.

#![forbid(unsafe_code)]
