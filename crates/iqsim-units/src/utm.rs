//! Inverse Universal Transverse Mercator on WGS84.
//!
//! Karney's (2011) Krüger series to sixth order in the third flattening n
//! ("Transverse Mercator with an accuracy of a few nanometers", J. Geodesy 85),
//! followed by Newton iteration from the conformal latitude. The series truncation
//! error is far below a millimetre within the UTM zones. Every transcendental is
//! from `libm` so results do not depend on the platform libm (AGENTS.md).

/// WGS84 semi-major axis in metres.
pub const WGS84_A: f64 = 6_378_137.0;
/// WGS84 flattening.
pub const WGS84_F: f64 = 1.0 / 298.257_223_563;
/// UTM central scale factor.
pub const UTM_K0: f64 = 0.9996;
/// UTM false easting in metres.
pub const FALSE_EASTING: f64 = 500_000.0;
/// UTM false northing in metres, southern hemisphere.
pub const FALSE_NORTHING_SOUTH: f64 = 10_000_000.0;

/// Hemisphere of a UTM coordinate (selects the false northing).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Hemisphere {
    /// Northern: false northing 0.
    North,
    /// Southern: false northing 10 000 km.
    South,
}

/// Geodetic latitude and longitude of a UTM coordinate, in degrees. The longitude
/// is normalized to [-180, 180). `zone` must be 1..=60.
pub fn utm_to_geodetic(
    zone: u8,
    hemisphere: Hemisphere,
    easting: f64,
    northing: f64,
) -> (f64, f64) {
    debug_assert!((1..=60).contains(&zone));
    let n = WGS84_F / (2.0 - WGS84_F);
    let n2 = n * n;
    let n3 = n2 * n;
    let n4 = n3 * n;
    let n5 = n4 * n;
    let n6 = n5 * n;
    // Rectifying radius A = a/(1+n) · (1 + n²/4 + n⁴/64 + n⁶/256).
    let a_rect = WGS84_A / (1.0 + n) * (1.0 + n2 / 4.0 + n4 / 64.0 + n6 / 256.0);
    // Karney (2011) eq. 36: β coefficients, Horner form written as explicit terms.
    let beta = [
        n / 2.0 - 2.0 * n2 / 3.0 + 37.0 * n3 / 96.0 - n4 / 360.0 - 81.0 * n5 / 512.0
            + 96_199.0 * n6 / 604_800.0,
        n2 / 48.0 + n3 / 15.0 - 437.0 * n4 / 1440.0 + 46.0 * n5 / 105.0
            - 1_118_711.0 * n6 / 3_870_720.0,
        17.0 * n3 / 480.0 - 37.0 * n4 / 840.0 - 209.0 * n5 / 4480.0 + 5569.0 * n6 / 90_720.0,
        4397.0 * n4 / 161_280.0 - 11.0 * n5 / 504.0 - 830_251.0 * n6 / 7_257_600.0,
        4583.0 * n5 / 161_280.0 - 108_847.0 * n6 / 3_991_680.0,
        20_648_693.0 * n6 / 638_668_800.0,
    ];

    let y = match hemisphere {
        Hemisphere::North => northing,
        Hemisphere::South => northing - FALSE_NORTHING_SOUTH,
    };
    let xi = y / (UTM_K0 * a_rect);
    let eta = (easting - FALSE_EASTING) / (UTM_K0 * a_rect);

    let mut xi_p = xi;
    let mut eta_p = eta;
    for (j, b) in beta.iter().enumerate() {
        let k = 2.0 * (j as f64 + 1.0);
        xi_p -= b * libm::sin(k * xi) * libm::cosh(k * eta);
        eta_p -= b * libm::cos(k * xi) * libm::sinh(k * eta);
    }

    let sinh_eta = libm::sinh(eta_p);
    let sin_xi = libm::sin(xi_p);
    let cos_xi = libm::cos(xi_p);
    // Conformal latitude τ' = tan χ.
    let tau_p = sin_xi / (sinh_eta * sinh_eta + cos_xi * cos_xi).sqrt();
    let tau = conformal_to_geodetic_tan(tau_p);

    let lat = libm::atan(tau).to_degrees();
    let central = f64::from(zone) * 6.0 - 183.0;
    let mut lon = libm::atan2(sinh_eta, cos_xi).to_degrees() + central;
    if lon >= 180.0 {
        lon -= 360.0;
    } else if lon < -180.0 {
        lon += 360.0;
    }
    (lat, lon)
}

/// Solves tan φ from the conformal tan χ = τ' by Newton's method (Karney 2011,
/// eqs. 19–21). Converges to full precision in two or three steps.
fn conformal_to_geodetic_tan(tau_p: f64) -> f64 {
    let e2 = WGS84_F * (2.0 - WGS84_F);
    let e = e2.sqrt();
    let one_minus_e2 = 1.0 - e2;
    let mut tau = tau_p;
    for _ in 0..8 {
        let s1 = (1.0 + tau * tau).sqrt();
        let sigma = libm::sinh(e * libm::atanh(e * tau / s1));
        let tau_i = tau * (1.0 + sigma * sigma).sqrt() - sigma * s1;
        let d_tau = (tau_p - tau_i) / (1.0 + tau_i * tau_i).sqrt()
            * (1.0 + one_minus_e2 * tau * tau)
            / (one_minus_e2 * s1);
        tau += d_tau;
        if d_tau.abs() <= 1e-14 * tau.abs().max(1.0) {
            break;
        }
    }
    tau
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn central_meridian_and_equator() {
        let (lat, lon) = utm_to_geodetic(31, Hemisphere::North, FALSE_EASTING, 0.0);
        assert_eq!(lat, 0.0);
        assert_eq!(lon, 3.0);
        let (lat, lon) =
            utm_to_geodetic(31, Hemisphere::South, FALSE_EASTING, FALSE_NORTHING_SOUTH);
        assert_eq!(lat, 0.0);
        assert_eq!(lon, 3.0);
        // Zone 1's central meridian is -177 and zone 60's is 177.
        assert_eq!(
            utm_to_geodetic(1, Hemisphere::North, FALSE_EASTING, 0.0).1,
            -177.0
        );
        assert_eq!(
            utm_to_geodetic(60, Hemisphere::North, FALSE_EASTING, 0.0).1,
            177.0
        );
    }
}
