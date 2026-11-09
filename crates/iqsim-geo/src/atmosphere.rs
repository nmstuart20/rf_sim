//! ICAO Standard Atmosphere and pressure altitude (contracts §12.6).
//!
//! **ISA** in the U.S. Standard Atmosphere 1976 form (the same layers as ICAO
//! Doc 7488 below 51 km): `g₀ = 9.806 65 m/s²`, `R = R*/M₀ = 8314.32/28.9644
//! = 287.053 07 J/(kg·K)` (Doc 7488 rounds R to 287.052 87, a 7·10⁻⁷ relative
//! difference, about 8 mm of pressure altitude at FL400), sea level
//! `T₀ = 288.15 K`, `p₀ = 101 325 Pa`, and these layers in **geopotential** height H:
//!
//! | base H (m) | lapse L (K/m) | base T (K) |
//! | ---: | ---: | ---: |
//! | −5 000 | −0.0065 | 320.65 |
//! | 0 | −0.0065 | 288.15 |
//! | 11 000 | 0 | 216.65 |
//! | 20 000 | +0.001 | 216.65 |
//! | 32 000 | +0.0028 | 228.65 |
//! | 47 000 | 0 | 270.65 |
//! | 51 000 (top) | | |
//!
//! (The first row is the troposphere continued below sea level, as Doc 7488's tables
//! do.) Base pressures are computed from p₀ layer by layer with
//! `p = p_b·(T_b/(T_b + L·(H − H_b)))^(g₀/(R·L))`, or `p = p_b·exp(−g₀·(H − H_b)/(R·T_b))`
//! when `L = 0`, and inverted in closed form. Mode C needs −1 200 ft to 126 700 ft
//! (−366 m to 38.6 km), well inside the range.
//!
//! **Geometric vs geopotential.** `H = r₀·h/(r₀ + h)` and `h = r₀·H/(r₀ − H)`, with
//! `r₀ = 6 356 766 m`. ISA altitudes (and so pressure altitudes) are geopotential; a
//! platform's height above MSL is geometric. At FL350 the difference is 59 m.
//!
//! **Local pressure offset.** The local atmosphere is taken as the ISA temperature
//! profile with sea-level pressure `p₀ + Δp` (Δp = QNH − 1013.25 hPa). By hydrostatics
//! that scales the pressure at every height by the same factor, so
//! `p_local(H) = (p₀ + Δp)/p₀ · p_ISA(H)` holds in every layer. Pressure altitude is
//! `H_p = p_ISA⁻¹(p_local(H))`.
//!
//! All transcendentals go through `libm` (`pow`, `exp`, `log`).

use crate::GeoError;

/// Standard gravity, m/s².
pub const G0: f64 = 9.806_65;
/// Specific gas constant of dry air, `R*/M₀`, J/(kg·K) (US 1976).
pub const R_AIR: f64 = 8_314.32 / 28.964_4;
/// ISA sea-level pressure, Pa.
pub const P0_PA: f64 = 101_325.0;
/// ISA sea-level temperature, K.
pub const T0_K: f64 = 288.15;
/// Earth radius for geopotential height, m.
pub const R0_GEOPOTENTIAL_M: f64 = 6_356_766.0;
/// Lowest ISA geopotential height this model covers, m.
pub const ISA_MIN_H_M: f64 = -5_000.0;
/// Highest ISA geopotential height this model covers, m.
pub const ISA_MAX_H_M: f64 = 51_000.0;
/// Largest accepted |Δp|, Pa (QNH 813–1213 hPa).
pub const MAX_PRESSURE_OFFSET_PA: f64 = 20_000.0;
/// One foot, metres.
pub const FT_M: f64 = 0.3048;

/// `(base H, lapse)`; the base temperatures and pressures follow from these.
const LAYERS: [(f64, f64); 6] = [
    (-5_000.0, -0.0065),
    (0.0, -0.0065),
    (11_000.0, 0.0),
    (20_000.0, 0.001),
    (32_000.0, 0.0028),
    (47_000.0, 0.0),
];

/// One layer with its base temperature and pressure.
#[derive(Clone, Copy, Debug)]
struct Layer {
    h_base: f64,
    lapse: f64,
    t_base: f64,
    p_base: f64,
    h_top: f64,
}

/// `p_b·ratio` for height `h` inside the layer.
fn layer_pressure(l: &Layer, h: f64) -> f64 {
    let dh = h - l.h_base;
    if l.lapse == 0.0 {
        l.p_base * libm::exp(-(G0 * dh) / (R_AIR * l.t_base))
    } else {
        let t = l.t_base + l.lapse * dh;
        l.p_base * libm::pow(l.t_base / t, G0 / (R_AIR * l.lapse))
    }
}

fn layer_height(l: &Layer, p: f64) -> f64 {
    let ratio = p / l.p_base;
    if l.lapse == 0.0 {
        l.h_base - (R_AIR * l.t_base / G0) * libm::log(ratio)
    } else {
        let k = libm::pow(ratio, -(R_AIR * l.lapse) / G0);
        l.h_base + (l.t_base / l.lapse) * (k - 1.0)
    }
}

/// The layers with base temperatures and pressures, from sea level outwards (the
/// layer below sea level is computed downward from p₀).
fn layers() -> [Layer; 6] {
    let mut out = [Layer {
        h_base: 0.0,
        lapse: 0.0,
        t_base: 0.0,
        p_base: 0.0,
        h_top: 0.0,
    }; 6];
    let below = LAYERS[0];
    let t_below = T0_K + below.1 * (below.0 - 0.0);
    out[1] = Layer {
        h_base: 0.0,
        lapse: LAYERS[1].1,
        t_base: T0_K,
        p_base: P0_PA,
        h_top: LAYERS[2].0,
    };
    for i in 2..6 {
        let prev = out[i - 1];
        let t_top = prev.t_base + prev.lapse * (prev.h_top - prev.h_base);
        out[i] = Layer {
            h_base: LAYERS[i].0,
            lapse: LAYERS[i].1,
            t_base: t_top,
            p_base: layer_pressure(&prev, prev.h_top),
            h_top: if i + 1 < 6 {
                LAYERS[i + 1].0
            } else {
                ISA_MAX_H_M
            },
        };
    }
    let p_below = layer_pressure(&out[1], below.0);
    out[0] = Layer {
        h_base: below.0,
        lapse: below.1,
        t_base: t_below,
        p_base: p_below,
        h_top: 0.0,
    };
    out
}

/// ISA pressure at geopotential height `h_m`, Pa.
pub fn isa_pressure_pa(h_m: f64) -> Result<f64, GeoError> {
    if !h_m.is_finite() {
        return Err(GeoError::NonFinite("ISA height"));
    }
    if !(ISA_MIN_H_M..=ISA_MAX_H_M).contains(&h_m) {
        return Err(GeoError::AtmosphereRange {
            value: h_m,
            what: "geopotential height, m",
        });
    }
    let ls = layers();
    let l = ls.iter().rev().find(|l| h_m >= l.h_base).unwrap_or(&ls[0]);
    Ok(layer_pressure(l, h_m))
}

/// ISA temperature at geopotential height `h_m`, K.
pub fn isa_temperature_k(h_m: f64) -> Result<f64, GeoError> {
    isa_pressure_pa(h_m)?;
    let ls = layers();
    let l = ls.iter().rev().find(|l| h_m >= l.h_base).unwrap_or(&ls[0]);
    Ok(l.t_base + l.lapse * (h_m - l.h_base))
}

/// The ISA geopotential height at which the pressure is `p_pa` (the pressure
/// altitude of `p_pa`), metres.
pub fn isa_pressure_altitude_m(p_pa: f64) -> Result<f64, GeoError> {
    if !p_pa.is_finite() {
        return Err(GeoError::NonFinite("pressure"));
    }
    let ls = layers();
    let p_top = layer_pressure(&ls[5], ISA_MAX_H_M);
    if !(p_pa >= p_top && p_pa <= ls[0].p_base) {
        return Err(GeoError::AtmosphereRange {
            value: p_pa,
            what: "pressure, Pa",
        });
    }
    let l = ls.iter().rev().find(|l| p_pa <= l.p_base).unwrap_or(&ls[0]);
    Ok(layer_height(l, p_pa))
}

/// Geopotential height of geometric height `h_m`: `r₀·h/(r₀ + h)`.
pub fn geopotential_from_geometric_m(h_m: f64) -> f64 {
    R0_GEOPOTENTIAL_M * h_m / (R0_GEOPOTENTIAL_M + h_m)
}

/// Geometric height of geopotential height `h_m`: `r₀·H/(r₀ − H)`.
pub fn geometric_from_geopotential_m(h_m: f64) -> f64 {
    R0_GEOPOTENTIAL_M * h_m / (R0_GEOPOTENTIAL_M - h_m)
}

/// The local atmosphere: ISA with a sea-level pressure offset (module docs).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Atmosphere {
    pressure_offset_pa: f64,
}

impl Atmosphere {
    /// The ISA standard day (Δp = 0).
    pub const STANDARD: Self = Self {
        pressure_offset_pa: 0.0,
    };

    /// ISA with sea-level pressure `p₀ + pressure_offset_pa`; |Δp| ≤
    /// [`MAX_PRESSURE_OFFSET_PA`].
    pub fn new(pressure_offset_pa: f64) -> Result<Self, GeoError> {
        if !pressure_offset_pa.is_finite() {
            return Err(GeoError::NonFinite("pressure offset"));
        }
        if pressure_offset_pa.abs() > MAX_PRESSURE_OFFSET_PA {
            return Err(GeoError::AtmosphereRange {
                value: pressure_offset_pa,
                what: "pressure offset, Pa",
            });
        }
        Ok(Self { pressure_offset_pa })
    }

    /// Δp, Pa.
    pub const fn pressure_offset_pa(&self) -> f64 {
        self.pressure_offset_pa
    }

    fn factor(&self) -> f64 {
        (P0_PA + self.pressure_offset_pa) / P0_PA
    }

    /// Static pressure at geometric height `h_msl_m` above MSL, Pa.
    pub fn pressure_at_msl_height_pa(&self, h_msl_m: f64) -> Result<f64, GeoError> {
        let p = isa_pressure_pa(geopotential_from_geometric_m(h_msl_m))?;
        Ok(self.factor() * p)
    }

    /// ISA pressure altitude (geopotential metres) of a platform at geometric
    /// height `h_msl_m` above MSL.
    pub fn pressure_altitude_m(&self, h_msl_m: f64) -> Result<f64, GeoError> {
        isa_pressure_altitude_m(self.pressure_at_msl_height_pa(h_msl_m)?)
    }

    /// The geometric height above MSL whose pressure altitude is `pressure_alt_m`
    /// (the inverse of [`Self::pressure_altitude_m`]).
    pub fn msl_height_for_pressure_altitude_m(&self, pressure_alt_m: f64) -> Result<f64, GeoError> {
        let p = isa_pressure_pa(pressure_alt_m)?;
        let h = isa_pressure_altitude_m(p / self.factor())?;
        Ok(geometric_from_geopotential_m(h))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn layer_bases_match_published_values() {
        // U.S. Standard Atmosphere 1976, Table 4 (= ICAO Doc 7488 below 51 km):
        // base pressures at geopotential 11, 20, 32, 47 and 51 km, 7 significant digits.
        for (h, p_pub, t_pub) in [
            (0.0, 101_325.0, 288.15),
            (11_000.0, 22_632.06, 216.65),
            (20_000.0, 5_474.889, 216.65),
            (32_000.0, 868.018_7, 228.65),
            (47_000.0, 110.906_3, 270.65),
            (51_000.0, 66.938_87, 270.65),
        ] {
            let p = isa_pressure_pa(h).unwrap();
            assert!(((p - p_pub) / p_pub).abs() < 2e-7, "{h}: {p} vs {p_pub}");
            let t = isa_temperature_k(h).unwrap();
            assert!((t - t_pub).abs() < 1e-9, "{h}: {t}");
        }
        assert!((isa_temperature_k(-5_000.0).unwrap() - 320.65).abs() < 1e-9);
    }

    #[test]
    fn round_trips_and_continuity() {
        let mut h = -4_999.0;
        while h <= 50_999.0 {
            let p = isa_pressure_pa(h).unwrap();
            let back = isa_pressure_altitude_m(p).unwrap();
            assert!((back - h).abs() < 1e-7, "{h}: {back}");
            h += 37.3;
        }
        // Continuity at the layer boundaries.
        for b in [0.0_f64, 11_000.0, 20_000.0, 32_000.0, 47_000.0] {
            let lo = isa_pressure_pa(b - 1e-6).unwrap();
            let hi = isa_pressure_pa(b).unwrap();
            assert!(((lo - hi) / hi).abs() < 1e-9, "{b}");
        }
        assert!(isa_pressure_pa(51_000.1).is_err());
        assert!(isa_pressure_pa(-5_000.1).is_err());
        assert!(isa_pressure_altitude_m(1e6).is_err());
        assert!(isa_pressure_pa(f64::NAN).is_err());
    }

    #[test]
    fn geopotential_conversion() {
        assert_eq!(geopotential_from_geometric_m(0.0), 0.0);
        let h = 10_668.0;
        let gp = geopotential_from_geometric_m(h);
        assert!((h - gp - 17.88).abs() < 0.01, "{gp}");
        assert!((geometric_from_geopotential_m(gp) - h).abs() < 1e-9);
    }

    #[test]
    fn pressure_offset_shifts_pressure_altitude() {
        // Troposphere closed form, independent of the layer code:
        // H_p = (T₀/L)·(1 − (p/p₀)^(−R·L/g₀)) with L = −0.0065 written as T₀/|L|.
        let atm = Atmosphere::new(1_000.0).unwrap(); // QNH 1023.25 hPa
        let p = atm.pressure_at_msl_height_pa(0.0).unwrap();
        assert!((p - 102_325.0).abs() < 1e-9);
        let hp = atm.pressure_altitude_m(0.0).unwrap();
        let l = 0.0065;
        let want = (T0_K / l) * (1.0 - libm::pow(p / P0_PA, R_AIR * l / G0));
        assert!((hp - want).abs() < 1e-9, "{hp} {want}");
        // About −83 m (−273 ft): 27.3 ft per hPa near sea level.
        assert!((hp / FT_M + 273.0).abs() < 2.0, "{}", hp / FT_M);
        // The inverse returns the MSL height.
        for h in [-300.0, 0.0, 3_000.0, 10_668.0, 38_000.0] {
            let hp = atm.pressure_altitude_m(h).unwrap();
            let back = atm.msl_height_for_pressure_altitude_m(hp).unwrap();
            assert!((back - h).abs() < 1e-6, "{h}: {back}");
        }
        // Standard day: pressure altitude is the geopotential height.
        let s = Atmosphere::STANDARD;
        assert!(
            (s.pressure_altitude_m(10_668.0).unwrap() - geopotential_from_geometric_m(10_668.0))
                .abs()
                < 1e-7
        );
        assert!(Atmosphere::new(20_001.0).is_err());
        assert!(Atmosphere::new(f64::INFINITY).is_err());
    }
}
