//! The Mode C altitude chain (plan, "World model": Mode C altitude; contracts §12.6).
//!
//! ```text
//! ellipsoid height h ──(− N(φ, λ), geoid)──▶ MSL height ──(geopotential, ISA, Δp)──▶
//! pressure altitude ──(nearest 100 ft, Gillham)──▶ Mode C code
//! ```
//!
//! Motion is geometric (ellipsoid heights); a transponder reports pressure altitude.
//! [`AltitudeModel`] holds the geoid and the local atmosphere and runs the chain in
//! both directions; the inverse resolves scenario positions given as flight levels.

use crate::GeoError;
use crate::atmosphere::{Atmosphere, FT_M};
use crate::geodetic::Geodetic;
use crate::geoid::Geoid;
use crate::gillham::GillhamCode;

/// The geoid and atmosphere a run uses.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct AltitudeModel {
    /// Geoid model.
    pub geoid: Geoid,
    /// Local atmosphere.
    pub atmosphere: Atmosphere,
}

impl AltitudeModel {
    /// A model from its parts.
    pub fn new(geoid: Geoid, atmosphere: Atmosphere) -> Self {
        Self { geoid, atmosphere }
    }

    /// Height above MSL of a geodetic position: `h − N(φ, λ)`, metres.
    pub fn msl_height_m(&self, p: Geodetic) -> Result<f64, GeoError> {
        Ok(p.h_m - self.geoid.height_m(p.lat_rad, p.lon_rad)?)
    }

    /// Ellipsoid height of a point `msl_m` metres above MSL at (φ, λ) radians.
    pub fn ellipsoid_height_from_msl_m(
        &self,
        lat_rad: f64,
        lon_rad: f64,
        msl_m: f64,
    ) -> Result<f64, GeoError> {
        Ok(msl_m + self.geoid.height_m(lat_rad, lon_rad)?)
    }

    /// ISA pressure altitude of a geodetic position, metres (geopotential).
    pub fn pressure_altitude_m(&self, p: Geodetic) -> Result<f64, GeoError> {
        self.atmosphere.pressure_altitude_m(self.msl_height_m(p)?)
    }

    /// Pressure altitude of a geodetic position, feet.
    pub fn pressure_altitude_ft(&self, p: Geodetic) -> Result<f64, GeoError> {
        Ok(self.pressure_altitude_m(p)? / FT_M)
    }

    /// The Mode C (Gillham) code a transponder at `p` reports.
    pub fn mode_c_code(&self, p: Geodetic) -> Result<GillhamCode, GeoError> {
        GillhamCode::encode_pressure_altitude_ft(self.pressure_altitude_ft(p)?)
            .map_err(GeoError::Gillham)
    }

    /// Ellipsoid height at (φ, λ) radians of the point whose pressure altitude is
    /// `pressure_alt_m` (geopotential metres), e.g. a flight level.
    pub fn ellipsoid_height_from_pressure_altitude_m(
        &self,
        lat_rad: f64,
        lon_rad: f64,
        pressure_alt_m: f64,
    ) -> Result<f64, GeoError> {
        let msl = self
            .atmosphere
            .msl_height_for_pressure_altitude_m(pressure_alt_m)?;
        self.ellipsoid_height_from_msl_m(lat_rad, lon_rad, msl)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn flight_level_round_trips_to_its_code() {
        let m = AltitudeModel::new(Geoid::Zero, Atmosphere::new(-1_500.0).unwrap());
        let (lat, lon) = (0.7, -1.3);
        for fl in [0, 5, 50, 180, 350, 450, 600, 1267] {
            let hp = f64::from(fl) * 100.0 * FT_M;
            let h = m
                .ellipsoid_height_from_pressure_altitude_m(lat, lon, hp)
                .unwrap();
            let code = m.mode_c_code(Geodetic::new(lat, lon, h)).unwrap();
            assert_eq!(code.altitude_ft(), fl * 100, "FL{fl}");
        }
        // The standard-day FL350 is geometrically higher than 35 000 ft (geopotential).
        let s = AltitudeModel::default();
        let h = s
            .ellipsoid_height_from_pressure_altitude_m(0.0, 0.0, 10_668.0)
            .unwrap();
        assert!((h - 10_685.93).abs() < 0.01, "{h}");
    }
}
