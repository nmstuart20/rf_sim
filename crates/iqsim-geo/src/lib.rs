//! Reference frames, trajectories, and propagation geometry.
//! Per the plan sections "World model" and "Antennas and propagation", this crate
//! owns WGS84/ECEF/ECI/NED transforms, C2 motion, attitude, geoid altitude,
//! and light-time geometry that retains fractional acquisition time.
//!
//! Dependencies: iqsim-core, iqsim-units.

#![forbid(unsafe_code)]
