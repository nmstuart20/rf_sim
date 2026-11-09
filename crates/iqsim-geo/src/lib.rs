//! Reference frames, trajectories, and propagation geometry.
//! Per the plan sections "World model" and "Antennas and propagation", this crate
//! owns WGS84/ECEF/ECI/NED transforms, C2 motion, attitude, geoid altitude,
//! and light-time geometry that retains fractional acquisition time.
//!
//! Dependencies: iqsim-core, iqsim-units.

#![forbid(unsafe_code)]

/// Speed of light in vacuum, m/s (exact, SI).
pub const SPEED_OF_LIGHT_M_S: f64 = 299_792_458.0;

/// A WGS84 Earth-centred, Earth-fixed position, metres.
///
/// Prompt 5 needs only fixed positions and straight-line distance. Frames, light time
/// and motion arrive in phase 3.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Ecef {
    /// x, metres.
    pub x: f64,
    /// y, metres.
    pub y: f64,
    /// z, metres.
    pub z: f64,
}

impl Ecef {
    /// The position (x, y, z) metres.
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Straight-line distance to `other`, metres (IEEE sqrt, fixed operation order).
    pub fn distance(self, other: Ecef) -> f64 {
        let dx = other.x - self.x;
        let dy = other.y - self.y;
        let dz = other.z - self.z;
        let s = dx * dx + dy * dy;
        let s = s + dz * dz;
        s.sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn distance_is_euclidean() {
        let a = Ecef::new(1.0, 2.0, 3.0);
        assert_eq!(a.distance(Ecef::new(4.0, 6.0, 3.0)), 5.0);
        assert_eq!(a.distance(a), 0.0);
    }
}
