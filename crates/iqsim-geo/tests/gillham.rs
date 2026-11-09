//! Gillham encoding checked against the AirSport Avionics reference table.

#![forbid(unsafe_code)]

use iqsim_geo::gillham::{GILLHAM_MAX_FT, GILLHAM_MIN_FT, GillhamBits, GillhamCode, GillhamError};

#[test]
fn every_airsport_reference_row_matches() {
    let table = include_str!("data/airsport_modecascii.txt");
    let mut row_count = 0;
    let mut last_altitude = None;
    for line in table.lines() {
        let columns: Vec<_> = line.split_whitespace().collect();
        let Some(altitude) = columns
            .first()
            .and_then(|column| column.parse::<i32>().ok())
        else {
            // Title, column headings, and blank lines precede the data.
            continue;
        };
        assert_eq!(columns.len(), 14, "malformed table row: {line}");
        assert_eq!(altitude, GILLHAM_MIN_FT + row_count * 100);
        let bit = |index: usize| match columns[index] {
            "0" => false,
            "1" => true,
            other => panic!("invalid table bit {other:?} at {altitude} ft"),
        };
        let bits = GillhamBits {
            a1: bit(1),
            a2: bit(2),
            a4: bit(3),
            b1: bit(4),
            b2: bit(5),
            b4: bit(6),
            c1: bit(7),
            c2: bit(8),
            c4: bit(9),
            d1: bit(10),
            d2: bit(11),
            d4: bit(12),
        };
        let code = GillhamCode::from_altitude_ft(altitude).unwrap();
        assert_eq!(code.bits(), bits, "bits at {altitude} ft");
        assert_eq!(code.squawk_string(), columns[13], "squawk at {altitude} ft");
        assert_eq!(
            code.squawk_octal(),
            u16::from_str_radix(columns[13], 8).unwrap(),
            "octal squawk at {altitude} ft"
        );
        assert_eq!(
            GillhamCode::from_bits(bits).unwrap().altitude_ft(),
            altitude
        );
        // Construct wire order independently from the reference column order.
        let mut wire = 0;
        for index in [10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9] {
            wire = (wire << 1) | u16::from(bit(index));
        }
        assert_eq!(code.wire_order(), wire, "wire order at {altitude} ft");
        assert_eq!(
            GillhamCode::from_wire_order(wire).unwrap().altitude_ft(),
            altitude
        );
        row_count += 1;
        last_altitude = Some(altitude);
    }
    assert_eq!(row_count, 1280);
    assert_eq!(last_altitude, Some(GILLHAM_MAX_FT));
}

#[test]
fn exhaustive_wire_values_round_trip() {
    let mut valid_count = 0;
    for wire in 0..4096 {
        if let Ok(code) = GillhamCode::from_wire_order(wire) {
            let altitude = code.altitude_ft();
            assert!((GILLHAM_MIN_FT..=GILLHAM_MAX_FT).contains(&altitude));
            assert_eq!(altitude % 100, 0);
            assert_eq!(
                GillhamCode::from_altitude_ft(altitude)
                    .unwrap()
                    .wire_order(),
                wire
            );
            assert_eq!(GillhamCode::from_bits(code.bits()).unwrap(), code);
            valid_count += 1;
        }
    }
    assert_eq!(valid_count, 1280);
}

#[test]
fn adjacent_altitudes_change_exactly_one_bit() {
    let mut previous = GillhamCode::from_altitude_ft(GILLHAM_MIN_FT)
        .unwrap()
        .wire_order();
    for altitude in ((GILLHAM_MIN_FT + 100)..=GILLHAM_MAX_FT).step_by(100) {
        let next = GillhamCode::from_altitude_ft(altitude)
            .unwrap()
            .wire_order();
        assert_eq!(
            (previous ^ next).count_ones(),
            1,
            "transition to {altitude} ft"
        );
        previous = next;
    }
}

#[test]
fn pressure_altitude_rounds_and_validates() {
    for (input, expected) in [
        (949.999, 900),
        (950.0, 1000),
        (1049.999, 1000),
        (-1249.9, -1200),
        (-1250.0, -1200),
        (126_749.9, 126_700),
        (-50.0, 0),
        (-0.0, 0),
    ] {
        assert_eq!(
            GillhamCode::encode_pressure_altitude_ft(input)
                .unwrap()
                .altitude_ft(),
            expected,
            "{input} ft"
        );
    }
    for input in [-1250.1, 126_750.0, f64::MIN, f64::MAX] {
        assert!(matches!(
            GillhamCode::encode_pressure_altitude_ft(input),
            Err(GillhamError::OutOfRange { .. })
        ));
    }
    for input in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            GillhamCode::encode_pressure_altitude_ft(input),
            Err(GillhamError::NonFinite)
        );
    }
}

#[test]
fn every_rounding_transition_has_ties_toward_positive_infinity() {
    for lower in ((GILLHAM_MIN_FT - 100)..=GILLHAM_MAX_FT).step_by(100) {
        let transition = f64::from(lower) + 50.0;
        for (input, expected) in [
            (transition.next_down(), lower),
            (transition, lower + 100),
            (transition.next_up(), lower + 100),
        ] {
            let result = GillhamCode::encode_pressure_altitude_ft(input);
            if (GILLHAM_MIN_FT..=GILLHAM_MAX_FT).contains(&expected) {
                assert_eq!(result.unwrap().altitude_ft(), expected, "{input:?} ft");
            } else {
                assert!(
                    matches!(result, Err(GillhamError::OutOfRange { .. })),
                    "{input:?} ft"
                );
            }
        }
    }
}

#[test]
fn illegal_c_patterns_and_d1_are_rejected() {
    for c_bits in [0b000, 0b101, 0b111] {
        let bits = GillhamBits {
            c1: c_bits & 4 != 0,
            c2: c_bits & 2 != 0,
            c4: c_bits & 1 != 0,
            ..GillhamBits::default()
        };
        assert_eq!(
            GillhamCode::from_bits(bits),
            Err(GillhamError::IllegalC { c_bits })
        );
        assert_eq!(
            GillhamCode::from_wire_order(c_bits),
            Err(GillhamError::IllegalC { c_bits })
        );
    }
    let bits = GillhamBits {
        d1: true,
        c4: true,
        ..GillhamBits::default()
    };
    assert_eq!(GillhamCode::from_bits(bits), Err(GillhamError::D1Set));
    assert_eq!(
        GillhamCode::from_wire_order(0x801),
        Err(GillhamError::D1Set)
    );
}

#[test]
fn wire_bits_above_bit_11_are_rejected() {
    // Each high-bit combination is invalid even with an otherwise legal code.
    for high in 1..16 {
        let bits = (high << 12) | 1;
        assert_eq!(
            GillhamCode::from_wire_order(bits),
            Err(GillhamError::WireBitsOutOfRange { bits })
        );
    }
}

#[test]
fn exact_altitude_inputs_are_checked_before_arithmetic() {
    for altitude_ft in [
        i32::MIN,
        GILLHAM_MIN_FT - 100,
        GILLHAM_MAX_FT + 100,
        i32::MAX,
    ] {
        assert_eq!(
            GillhamCode::from_altitude_ft(altitude_ft),
            Err(GillhamError::OutOfRange {
                altitude_ft: f64::from(altitude_ft)
            })
        );
    }
    for altitude_ft in [GILLHAM_MIN_FT + 1, -1, 1, GILLHAM_MAX_FT - 1] {
        assert_eq!(
            GillhamCode::from_altitude_ft(altitude_ft),
            Err(GillhamError::NotMultipleOf100 { altitude_ft })
        );
    }
}
