//! MGRS and UTM against independent published references.
//!
//! - `data/mgrs_geoconvert.txt`: 60 references (hand-picked edge cases plus random
//!   points; encoded by the NGA GEOTRANS-based Python `mgrs` package) decoded by
//!   GeographicLib's GeoConvert 2.5.2 web service, which returns square centres.
//! - `data/utm_pyproj.txt`: 400 random UTM coordinates inverted by PROJ 9.8.1.
//!
//! `data/gen_oracle.py` regenerates both files.

#![forbid(unsafe_code)]

use iqsim_units::utm::{Hemisphere, utm_to_geodetic};
use iqsim_units::{Mgrs, Position};

/// 10⁻¹⁰ degrees is at most 11 µm on the ground. The measured worst cases are
/// 5·10⁻¹² (GeoConvert prints 11 decimals) and 5·10⁻¹³ (PROJ, 12 decimals).
const TOL_DEG: f64 = 1e-10;

fn rows(text: &str) -> impl Iterator<Item = Vec<&str>> {
    text.lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .map(|l| l.split_whitespace().collect())
}

fn zone_hemisphere(s: &str) -> (u8, Hemisphere) {
    let (zone, hemi) = s.split_at(s.len() - 1);
    let hemisphere = match hemi {
        "n" => Hemisphere::North,
        "s" => Hemisphere::South,
        _ => panic!("bad hemisphere {s}"),
    };
    (zone.parse().unwrap(), hemisphere)
}

/// Largest of the latitude and longitude differences, in degrees.
fn error_deg(lat: f64, lon: f64, ref_lat: f64, ref_lon: f64) -> f64 {
    let mut dlon = (lon - ref_lon).abs();
    if dlon > 180.0 {
        dlon = 360.0 - dlon;
    }
    (lat - ref_lat).abs().max(dlon)
}

#[test]
fn mgrs_matches_geoconvert() {
    let mut count = 0;
    let mut worst = 0.0f64;
    for row in rows(include_str!("data/mgrs_geoconvert.txt")) {
        let [reference, lat, lon, zh, e, n] = row[..] else {
            panic!("bad row {row:?}");
        };
        let m = Mgrs::parse(reference).unwrap_or_else(|e| panic!("{e}"));
        let (zone, hemisphere) = zone_hemisphere(zh);
        assert_eq!(
            (m.zone(), m.hemisphere()),
            (zone, hemisphere),
            "{reference}"
        );

        // The square's centre in UTM is exactly GeoConvert's easting/northing.
        let half = f64::from(m.square_size_m()) / 2.0;
        let (sw_e, sw_n) = m.utm_sw_corner();
        assert_eq!(
            f64::from(sw_e) + half,
            e.parse::<f64>().unwrap(),
            "{reference}"
        );
        assert_eq!(
            f64::from(sw_n) + half,
            n.parse::<f64>().unwrap(),
            "{reference}"
        );

        let (clat, clon) = m.center();
        let (rlat, rlon) = (lat.parse().unwrap(), lon.parse().unwrap());
        let error = error_deg(clat, clon, rlat, rlon);
        worst = worst.max(error);
        assert!(
            error <= TOL_DEG,
            "{reference}: got ({clat:.12}, {clon:.12}), GeoConvert ({rlat}, {rlon})"
        );
        count += 1;
    }
    assert_eq!(count, 60);
    println!("worst GeoConvert difference: {worst:e} deg");
}

#[test]
fn utm_inverse_matches_proj() {
    let mut count = 0;
    let mut worst = 0.0f64;
    for row in rows(include_str!("data/utm_pyproj.txt")) {
        let [zh, e, n, lat, lon] = row[..] else {
            panic!("bad row {row:?}");
        };
        let (zone, hemisphere) = zone_hemisphere(zh);
        let (glat, glon) =
            utm_to_geodetic(zone, hemisphere, e.parse().unwrap(), n.parse().unwrap());
        let (rlat, rlon) = (lat.parse().unwrap(), lon.parse().unwrap());
        let error = error_deg(glat, glon, rlat, rlon);
        worst = worst.max(error);
        assert!(
            error <= TOL_DEG,
            "{zh} {e} {n}: got ({glat:.12}, {glon:.12}), PROJ ({rlat}, {rlon})"
        );
        count += 1;
    }
    assert_eq!(count, 400);
    println!("worst PROJ difference: {worst:e} deg");
}

#[test]
fn position_uses_the_south_west_corner() {
    let position: Position =
        toml::from_str::<toml::Table>("p = { mgrs = \"18SUJ2337106519\", alt = \"10 m MSL\" }")
            .unwrap()["p"]
            .clone()
            .try_into()
            .unwrap();
    let m = Mgrs::parse("18SUJ2337106519").unwrap();
    assert_eq!(
        (
            position.latitude().degrees(),
            position.longitude().degrees()
        ),
        m.sw_corner()
    );
    assert_eq!(position.mgrs(), Some(&m));
    // The south-west corner is half a metre south-west of the centre.
    let (lat, lon) = m.sw_corner();
    let (clat, clon) = m.center();
    assert!(lat < clat && lon < clon);
}
