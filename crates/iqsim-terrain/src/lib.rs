//! Terrain assets and terrain-dependent propagation paths.
//! Per the plan section "Antennas and propagation", this crate owns DTED/GeoTIFF
//! elevation loading and tile caching, terrain visibility, diffraction,
//! and ground-reflection geometry with deterministic path selection.
//!
//! Dependencies: iqsim-core, iqsim-units, iqsim-geo.

#![forbid(unsafe_code)]
