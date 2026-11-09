//! GeographicLib PGM geoid grids: parsing, bilinear interpolation, wrap-around and
//! malformed files. Fixtures are synthetic grids written here.

#![forbid(unsafe_code)]

use iqsim_geo::geoid::{GEOID_FORMAT, GEOID_INTERPOLATION};
use iqsim_geo::{GeoError, Geoid, GeoidGrid};

const D: f64 = core::f64::consts::PI / 180.0;

/// A PGM with `width` columns (Δ = 360/width) and samples `f(row, col)`.
fn pgm(width: usize, offset: &str, scale: &str, f: impl Fn(usize, usize) -> u16) -> Vec<u8> {
    let height = width / 2 + 1;
    let mut out = format!(
        "P5\n# Geoid file in PGM format for the GeographicLib::Geoid class\n# Description test grid\n# Offset {offset}\n# Scale {scale}\n# Origin 90N 0E\n{width} {height}\n65535\n"
    )
    .into_bytes();
    for i in 0..height {
        for j in 0..width {
            out.extend_from_slice(&f(i, j).to_be_bytes());
        }
    }
    out
}

fn sample(i: usize, j: usize) -> u16 {
    (1000 + 37 * i + 101 * j + (i * j) % 7) as u16
}

#[test]
fn nodes_cells_and_wrap() {
    // 8 × 5 grid: Δ = 45°, rows at 90, 45, 0, −45, −90.
    let bytes = pgm(8, "-108", "0.003", sample);
    let g = GeoidGrid::from_pgm_bytes(&bytes).unwrap();
    let m = g.metadata();
    assert_eq!(
        (m.width, m.height, m.offset_m, m.scale_m),
        (8, 5, -108.0, 0.003)
    );
    assert_eq!(
        (m.format, m.interpolation),
        (GEOID_FORMAT, GEOID_INTERPOLATION)
    );
    assert_eq!(m.description, "test grid");
    let n = |i: usize, j: usize| -108.0 + 0.003 * f64::from(sample(i, j));

    // Every node is exact.
    for i in 0..5 {
        for j in 0..8 {
            let lat = 90.0 - 45.0 * i as f64;
            let lon = 45.0 * j as f64;
            let got = g.height_m(lat * D, lon * D).unwrap();
            assert!((got - n(i, j)).abs() < 1e-12, "({i},{j}) {got} {}", n(i, j));
        }
    }

    // Inside a cell: the bilinear weights, written independently.
    let bil = |lat: f64, lon: f64, i: usize, j: usize, j1: usize| {
        let ty = (90.0 - 45.0 * i as f64 - lat) / 45.0;
        let tx = (lon - 45.0 * j as f64) / 45.0;
        let tx = if tx < 0.0 { tx + 8.0 } else { tx };
        (1.0 - ty) * ((1.0 - tx) * n(i, j) + tx * n(i, j1))
            + ty * ((1.0 - tx) * n(i + 1, j) + tx * n(i + 1, j1))
    };
    for (lat, lon, i, j, j1) in [
        (30.0, 10.0, 1, 0, 1),
        (-12.5, 100.0, 2, 2, 3),
        (-80.0, 200.0, 3, 4, 5),
        (89.0, 300.0, 0, 6, 7),
        // Across 360°/0°: between column 7 (315°) and column 0.
        (10.0, 350.0, 1, 7, 0),
    ] {
        let got = g.height_m(lat * D, lon * D).unwrap();
        let want = bil(lat, lon, i, j, j1);
        assert!((got - want).abs() < 1e-9, "({lat},{lon}) {got} {want}");
    }
    // Longitude wraps: −10° is 350°, and 720° + 30° is 30°.
    let a = g.height_m(10.0 * D, -10.0 * D).unwrap();
    let b = g.height_m(10.0 * D, 350.0 * D).unwrap();
    assert!((a - b).abs() < 1e-9);
    let a = g.height_m(-20.0 * D, 750.0 * D).unwrap();
    let b = g.height_m(-20.0 * D, 30.0 * D).unwrap();
    assert!((a - b).abs() < 1e-9);
    // The surface is continuous across the 0° seam.
    let e = 1e-9;
    let w = g.height_m(20.0 * D, -e * D).unwrap();
    let o = g.height_m(20.0 * D, e * D).unwrap();
    assert!((w - o).abs() < 1e-6);
    // The south pole uses the last row.
    let s = g.height_m(-90.0 * D, 45.0 * D).unwrap();
    assert!((s - n(4, 1)).abs() < 1e-12);
    // Non-finite input.
    assert!(matches!(
        g.height_m(f64::NAN, 0.0),
        Err(GeoError::NonFinite(_))
    ));
}

#[test]
fn zero_geoid() {
    assert_eq!(Geoid::Zero.height_m(0.5, -2.0).unwrap(), 0.0);
    assert!(Geoid::Zero.height_m(f64::INFINITY, 0.0).is_err());
    let g = Geoid::Grid(GeoidGrid::from_pgm_bytes(&pgm(4, "0", "1", |_, _| 7)).unwrap());
    assert_eq!(g.height_m(0.3, 0.3).unwrap(), 7.0);
}

/// Review regressions: netpbm allows a comment right after maxval (its newline is
/// then the raster delimiter) and comments terminated by CR alone.
#[test]
fn netpbm_comment_forms() {
    // 6 × 4 grid: 24 samples, sample k = 3k m.
    let raster: Vec<u8> = (0..24u16).flat_map(|v| (v * 3).to_be_bytes()).collect();
    for header in [
        "P5\n# Offset 0\n# Scale 1\n6 4\n65535# tail\n",
        "P5\r# Offset 0\r# Scale 1\r6 4\r65535\r",
        "P5 #c\n6 # Offset 0\n# Scale 1\n4 65535\n",
    ] {
        let bytes = [header.as_bytes(), &raster[..]].concat();
        let g = GeoidGrid::from_pgm_bytes(&bytes).unwrap_or_else(|e| panic!("{header:?}: {e}"));
        assert_eq!((g.metadata().offset_m, g.metadata().scale_m), (0.0, 1.0));
        // Node (row 1, column 2) = sample 8 → 24 m at (30° N, 120° E).
        assert!((g.height_m(30.0 * D, 120.0 * D).unwrap() - 24.0).abs() < 1e-12);
    }
}

#[test]
fn malformed_files_are_rejected() {
    let good = pgm(8, "-108", "0.003", sample);
    assert!(GeoidGrid::from_pgm_bytes(&good).is_ok());
    let text = |b: &[u8]| String::from_utf8_lossy(b).into_owned();
    let cases: Vec<(&str, Vec<u8>)> = vec![
        ("magic", text(&good).replacen("P5", "P2", 1).into_bytes()),
        (
            "no offset",
            text(&good).replacen("# Offset -108\n", "", 1).into_bytes(),
        ),
        (
            "no scale",
            text(&good).replacen("# Scale 0.003\n", "", 1).into_bytes(),
        ),
        (
            "bad scale",
            text(&good)
                .replacen("Scale 0.003", "Scale -1", 1)
                .into_bytes(),
        ),
        (
            "bad offset",
            text(&good)
                .replacen("Offset -108", "Offset abc", 1)
                .into_bytes(),
        ),
        (
            "maxval",
            text(&good).replacen("65535", "255", 1).into_bytes(),
        ),
        ("shape", text(&good).replacen("8 5", "8 6", 1).into_bytes()),
        (
            "odd width",
            text(&good).replacen("8 5", "7 5", 1).into_bytes(),
        ),
        ("truncated", good[..good.len() - 1].to_vec()),
        ("extra", [good.clone(), vec![0]].concat()),
        ("empty", Vec::new()),
    ];
    for (name, bytes) in cases {
        let r = GeoidGrid::from_pgm_bytes(&bytes);
        assert!(matches!(r, Err(GeoError::GeoidFormat(_))), "{name}: {r:?}");
    }
    let r = GeoidGrid::from_pgm_file(std::path::Path::new("/nonexistent/egm96-5.pgm"));
    assert!(matches!(r, Err(GeoError::Io(_))));
}
