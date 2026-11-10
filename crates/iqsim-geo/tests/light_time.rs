//! Independent light-time references, fractional time arithmetic and lookback checks.

use std::sync::Arc;

use iqsim_core::{AcqInstant, TruthDuration, TruthTime};
use iqsim_geo::light_time::{
    LIGHT_TIME_MAX_ITERATIONS, LightTime, instant_offset, seconds_to_q64_ps, solve_emit,
    solve_receive,
};
use iqsim_geo::trajectory::{AttitudeMode, FlatGround, Heights, Limits, Primitive, Start};
use iqsim_geo::{
    AltitudeModel, EARTH_ROTATION_RATE_RAD_S as OMEGA, EarthOrientation, Ecef, Eci, Flight,
    GeoError, Geodetic, Lookback, MAX_LOOKBACK, Quat, SPEED_OF_LIGHT_M_S as C, Span, State, Static,
    Trajectory, Vec3, WGS84_A_M as R,
};

fn at(ps: i64) -> AcqInstant {
    AcqInstant::from_truth_time(TruthTime::from_ps(ps))
}

fn ground() -> Static {
    Static::from_ecef(Ecef::new(R, 0.0, 0.0), Quat::IDENTITY)
}

fn equator(longitude: f64) -> Static {
    Static::from_ecef(
        Ecef::new(R * libm::cos(longitude), R * libm::sin(longitude), 0.0),
        Quat::IDENTITY,
    )
}

// RX angle minus TX angle is (lambda_rx + theta(t)) minus
// (lambda_tx + theta(t) - omega*tau) = lambda_rx - lambda_tx + omega*tau.
// The plus sign is essential; reversing the baseline reverses the Sagnac term.
fn equatorial_reference(delta_longitude: f64) -> f64 {
    let residual =
        |tau: f64| 2.0 * R * libm::sin((delta_longitude + OMEGA * tau).abs() / 2.0) - C * tau;
    let (mut lo, mut hi) = (0.0, 2.0 * R / C);
    for _ in 0..100 {
        let mid = lo + (hi - lo) / 2.0;
        if hi - lo <= 1e-19 || mid == lo || mid == hi {
            break;
        }
        if residual(mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    lo + (hi - lo) / 2.0
}

#[test]
fn earth_fixed_equatorial_links_and_sagnac() {
    let mut worst_tau: f64 = 0.0;
    for baseline_m in [30_000.0, 3_000_000.0] {
        let angle = baseline_m / R;
        let no_rotation = 2.0 * R * libm::sin(angle / 2.0) / C;
        for sign in [-1.0, 1.0] {
            let rx = equator(sign * angle);
            let tx = equator(0.0);
            let reference = equatorial_reference(sign * angle);
            assert!((reference - no_rotation) * sign > 0.0);
            for ps in [0, 123_456_789_012_345, 8_640_000_000_000_000_000] {
                let solution = solve_receive(&rx, &tx, at(ps)).unwrap();
                let error = (solution.tau_s - reference).abs();
                worst_tau = worst_tau.max(error);
                assert!(
                    error <= 1e-15,
                    "baseline={baseline_m}, sign={sign}, t={ps}: {error:e}"
                );
                assert!(solution.iterations <= 5);
                check_residual(&solution);
            }
        }
    }
    eprintln!("equatorial reference: worst delay error = {worst_tau:e} s");
}

/// Constant ECI velocity, represented as a fully consistent ECEF trajectory.
#[derive(Debug)]
struct InertialLine {
    p0: Vec3,
    u: Vec3,
}

impl InertialLine {
    fn leo() -> Self {
        Self {
            p0: Vec3::new(7_000_000.0, 0.0, 0.0),
            u: Vec3::new(0.0, 7_500.0, 0.0),
        }
    }

    fn position(&self, t: AcqInstant) -> Vec3 {
        let local_ps =
            t.ps().as_ps() as f64 + t.frac() as f64 * (1.0 / 18_446_744_073_709_551_616.0);
        self.p0 + self.u * (local_ps / 1e12)
    }
}

impl Trajectory for InertialLine {
    fn state(&self, t: AcqInstant) -> Result<State, GeoError> {
        let orientation = EarthOrientation::at(t);
        let w = Vec3::new(0.0, 0.0, OMEGA);
        let (position, velocity) =
            orientation.state_to_ecef(Eci::from_vec(self.position(t)), self.u);
        // R^T a_eci = a_ecef + 2 omega x v_ecef + omega x (omega x r_ecef).
        let acceleration = w.cross(velocity) * -2.0 - w.cross(w.cross(position.vec()));
        Ok(State {
            position,
            velocity,
            acceleration,
            q_ecef_body: Quat::IDENTITY,
            angular_rate_body: Vec3::ZERO,
        })
    }

    fn span(&self) -> Span {
        Span {
            start: None,
            end: None,
        }
    }

    fn knots(&self, _from: TruthTime, _to: TruthTime) -> Vec<TruthTime> {
        Vec::new()
    }
}

// |d + u*tau|^2 = c^2*tau^2 -> A*tau^2 - 2*b*tau - d.d = 0.
// Rationalize the positive root when b < 0 to avoid cancellation.
fn quadratic_reference(rx: &Static, tx: &InertialLine, t: AcqInstant) -> (f64, f64) {
    let rx_state = rx.state(t).unwrap();
    let (r, v) = EarthOrientation::at(t).state_to_eci(rx_state.position, rx_state.velocity);
    let d = r.vec() - tx.position(t);
    let a = C * C - tx.u.dot(tx.u);
    let b = d.dot(tx.u);
    let d2 = d.dot(d);
    let root = libm::sqrt(b * b + a * d2);
    let mut tau = if b < 0.0 {
        d2 / (root - b)
    } else {
        (root + b) / a
    };
    // One Newton polish of the quadratic, independent of the solver's norm iteration.
    let residual = (a * tau - 2.0 * b) * tau - d2;
    tau -= residual / (2.0 * (a * tau - b));
    let d_rate = v - tx.u;
    let tau_rate = (d_rate.dot(tx.u) * tau + d.dot(d_rate)) / (a * tau - b);
    (tau, tau_rate)
}

fn check_residual(solution: &LightTime) -> f64 {
    let range = solution
        .rx
        .position_eci
        .vec()
        .distance(solution.tx.position_eci.vec());
    assert_eq!(range.to_bits(), solution.range_m.to_bits());
    assert_eq!((range / C).to_bits(), solution.tau_s.to_bits());
    let residual = (range - C * solution.tau_s).abs();
    assert!(residual <= 1e-6, "range residual {residual:e} m");
    residual
}

#[test]
fn uniform_eci_motion_delay_and_doppler() {
    // A polar receiver keeps its Earth-sized coordinate fixed under rotation.
    // This lets a 1 us difference resolve the derivative without differencing
    // roundoff from rotating Earth-sized x/y coordinates. Other tests exercise
    // equatorial rotation, and the TX still moves and rotates between frames.
    let rx = Static::from_ecef(
        Ecef::new(0.0, 0.0, R * (1.0 - iqsim_geo::WGS84_F)),
        Quat::IDENTITY,
    );
    let tx = InertialLine {
        p0: Vec3::new(0.0, 0.0, 7_000_000.0),
        u: Vec3::new(7_500.0, 0.0, 0.0),
    };
    let mut worst_tau: f64 = 0.0;
    let mut worst_derivative: f64 = 0.0;
    let mut worst_fd: f64 = 0.0;
    let mut worst_residual: f64 = 0.0;
    let mut iterations = 0;
    for ps in [
        -13_123_456_789_012,
        0,
        500_000_000_000,
        5_000_000_000_000,
        13_123_456_789_012,
    ] {
        let t = AcqInstant::new(TruthTime::from_ps(ps), 1 << 63);
        let solved = solve_receive(&rx, &tx, t).unwrap();
        let (reference, derivative) = quadratic_reference(&rx, &tx, t);
        worst_tau = worst_tau.max((solved.tau_s - reference).abs());
        worst_derivative = worst_derivative.max((solved.tau_rate - derivative).abs());
        let h = 1e-6;
        let before = solve_receive(&rx, &tx, instant_offset(t, -h).unwrap()).unwrap();
        let after = solve_receive(&rx, &tx, instant_offset(t, h).unwrap()).unwrap();
        let finite_difference = (after.tau_s - before.tau_s) / (2.0 * h);
        worst_fd = worst_fd.max((finite_difference - solved.tau_rate).abs());
        worst_residual = worst_residual.max(check_residual(&solved));
        iterations = iterations.max(solved.iterations);
        let evaluated_tx = tx.state(instant_offset(t, -solved.tau_s).unwrap()).unwrap();
        assert!(evaluated_tx.position.distance(solved.tx.state.position) <= 1e-8);
    }
    eprintln!(
        "LEO: delay {worst_tau:e} s, derivative {worst_derivative:e}, 1 us finite difference {worst_fd:e}, residual {worst_residual:e} m, iterations {iterations}"
    );
    assert!(worst_tau <= 1e-15);
    assert!(worst_derivative <= 1e-12);
    assert!(worst_fd <= 1e-12);
    assert!(iterations <= 5);
}

#[test]
fn emission_reception_inverse_and_geo_convergence() {
    let rx = ground();
    let leo = InertialLine::leo();
    let geo = Static::from_ecef(Ecef::new(0.0, 42_164_000.0, 0.0), Quat::IDENTITY);
    let mut worst_inverse_ps: f64 = 0.0;
    let mut worst_residual: f64 = 0.0;
    let mut max_iterations = 0;
    for tx in [&leo as &dyn Trajectory, &geo] {
        for ps in [0, 3_123_456_789_012, 60_000_000_000_000] {
            let emitted = AcqInstant::new(TruthTime::from_ps(ps), 0x1234_5678_9abc_def0);
            let forward = solve_emit(&rx, tx, emitted).unwrap();
            let reverse = solve_receive(&rx, tx, forward.rx.instant).unwrap();
            let delta = reverse.tx.instant.to_q64_64() - emitted.to_q64_64();
            let error_ps = (delta as f64 / 18_446_744_073_709_551_616.0).abs();
            worst_inverse_ps = worst_inverse_ps.max(error_ps);
            assert!(error_ps <= 1e-3, "inverse error {error_ps:e} ps");
            assert!((forward.tau_rate - reverse.tau_rate).abs() <= 1e-14);
            for solution in [&forward, &reverse] {
                worst_residual = worst_residual.max(check_residual(solution));
                max_iterations = max_iterations.max(solution.iterations);
            }
            let again = solve_emit(&rx, tx, reverse.tx.instant).unwrap();
            let delta = again.rx.instant.to_q64_64() - forward.rx.instant.to_q64_64();
            assert!((delta as f64 / 18_446_744_073_709_551_616.0).abs() <= 1e-3);
        }
    }
    eprintln!(
        "LEO/GEO: inverse {worst_inverse_ps:e} ps, residual {worst_residual:e} m, iterations {max_iterations}"
    );
    assert!(max_iterations <= 5);
}

#[test]
fn fractional_offsets_round_trip_and_validate() {
    let mut worst_units = 0;
    for ps in [
        i64::MIN + 1_000_000_000_001,
        -123,
        0,
        5,
        8_640_000_000_000_000_000,
    ] {
        for frac in [1, 1 << 63, u64::MAX] {
            let t = AcqInstant::new(TruthTime::from_ps(ps), frac);
            for s in [-1.0, -1e-4, -1e-12, -1e-24, 0.0, 1e-24, 1e-12, 1e-4, 1.0] {
                let shifted = instant_offset(t, s).unwrap();
                assert_eq!(
                    shifted.to_q64_64() - t.to_q64_64(),
                    seconds_to_q64_ps(s).unwrap()
                );
                let back = instant_offset(shifted, -s).unwrap();
                let error = (back.to_q64_64() - t.to_q64_64()).abs();
                worst_units = worst_units.max(error);
                assert!(error <= 2);
            }
        }
    }
    // Exact tiny offsets straddle ties in integer Q64.64 units. Dividing by 1e12
    // and multiplying back yields these exact binary half-integers.
    let scale = 18_446_744_073_709_551_616.0;
    for (units, rounded) in [(0.5, 0), (1.5, 2), (2.5, 2), (-1.5, -2)] {
        assert_eq!(seconds_to_q64_ps((units / scale) / 1e12).unwrap(), rounded);
    }
    assert_eq!(
        seconds_to_q64_ps(1e6).unwrap(),
        1_000_000_000_000_000_000_i128 << 64
    );
    assert!(seconds_to_q64_ps(-1e6).is_ok());
    for s in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(seconds_to_q64_ps(s), Err(GeoError::NonFinite(_))));
        assert!(matches!(
            instant_offset(at(0), s),
            Err(GeoError::NonFinite(_))
        ));
    }
    for s in [-1_000_001.0, 1_000_001.0] {
        assert!(matches!(seconds_to_q64_ps(s), Err(GeoError::LightTime(_))));
    }
    assert!(instant_offset(AcqInstant::from_q64_64(i128::MAX), 1e-12).is_err());
    assert!(instant_offset(AcqInstant::from_q64_64(i128::MIN), -1e-12).is_err());
    assert_eq!(
        instant_offset(AcqInstant::from_q64_64(i128::MAX), 0.0)
            .unwrap()
            .to_q64_64(),
        i128::MAX
    );
    eprintln!("offset round trip: worst error {worst_units} Q64.64 units");
}

fn flight(start: TruthTime) -> Flight {
    Flight::maneuvers(
        Start {
            time: start,
            position: Geodetic::from_degrees(51.0, -1.0, 3000.0),
            track_rad: Some(0.3),
            speed_m_s: 150.0,
        },
        Limits::default(),
        AttitudeMode::COORDINATED,
        &[Primitive::Straight {
            duration: TruthDuration::from_ps(2_000_000_000_000),
        }],
        Heights {
            model: &AltitudeModel::default(),
            ground: &FlatGround::ELLIPSOID,
        },
    )
    .unwrap()
}

fn state_bits(state: State) -> Vec<u64> {
    let q = state.q_ecef_body;
    state
        .position
        .vec()
        .to_array()
        .into_iter()
        .chain(state.velocity.to_array())
        .chain(state.acceleration.to_array())
        .chain([q.w, q.x, q.y, q.z])
        .chain(state.angular_rate_body.to_array())
        .map(f64::to_bits)
        .collect()
}

#[test]
fn exact_ecef_static_and_unbounded_lookback_pass_through() {
    let p = Ecef::new(-0.0, 6_378_137.123_456_789, 0.25);
    let q = Quat::from_euler_zyx(0.2, -0.3, 0.4);
    let inner = Arc::new(Static::from_ecef(p, q));
    let wrapper = Lookback::new(inner.clone(), MAX_LOOKBACK).unwrap();
    assert_eq!(wrapper.span(), inner.span());
    assert_eq!(
        wrapper.knots(TruthTime::MIN, TruthTime::MAX),
        inner.knots(TruthTime::MIN, TruthTime::MAX)
    );
    for ps in [i64::MIN, -1, 0, 1, i64::MAX] {
        let state = wrapper.state(at(ps)).unwrap();
        assert_eq!(state_bits(state), state_bits(inner.state(at(ps)).unwrap()));
        assert_eq!(state.position.x.to_bits(), p.x.to_bits());
        assert_eq!(state.position.y.to_bits(), p.y.to_bits());
        assert_eq!(state.position.z.to_bits(), p.z.to_bits());
        assert_eq!(state.q_ecef_body, q);
        assert_eq!(state.velocity, Vec3::ZERO);
        assert_eq!(state.acceleration, Vec3::ZERO);
        assert_eq!(state.angular_rate_body, Vec3::ZERO);
    }
}

#[test]
fn flight_lookback_join_span_knots_and_errors() {
    let start = TruthTime::from_ps(1_000_000_000_000);
    let inner = Arc::new(flight(start));
    let wrapper = Lookback::new(inner.clone(), MAX_LOOKBACK).unwrap();
    assert_eq!(
        wrapper.span(),
        Span {
            start: Some(TruthTime::EPOCH),
            end: inner.span().end
        }
    );
    for ps in [
        start.as_ps(),
        start.as_ps() + 1,
        start.as_ps() + 123_456_789,
        inner.span().end.unwrap().as_ps(),
    ] {
        let t = AcqInstant::new(TruthTime::from_ps(ps), 0);
        assert_eq!(
            state_bits(wrapper.state(t).unwrap()),
            state_bits(inner.state(t).unwrap())
        );
    }
    let initial = inner.state_at(start).unwrap();
    let before = wrapper.state(at(start.as_ps() - 1)).unwrap();
    let position_error = initial.position.distance(before.position);
    let velocity_error = initial.velocity.distance(before.velocity);
    let acceleration_error = initial.acceleration.distance(before.acceleration);
    assert!(position_error <= 1e-9);
    assert!(velocity_error <= 1e-9);
    assert!(acceleration_error <= 1e-9);
    for ps in [0, 500_000_000_000, start.as_ps() - 1] {
        let extrapolated = wrapper.state(at(ps)).unwrap();
        let dt = (ps - start.as_ps()) as f64 / 1e12;
        let expected_p = initial.position.vec()
            + (initial.velocity * dt + (initial.acceleration * (0.5 * dt)) * dt);
        assert_eq!(
            state_bits(extrapolated),
            state_bits(State {
                position: Ecef::from_vec(expected_p),
                velocity: initial.velocity + initial.acceleration * dt,
                ..initial
            })
        );
    }
    let knots = wrapper.knots(TruthTime::MIN, TruthTime::MAX);
    assert!(knots.contains(&start));
    assert!(knots.windows(2).all(|pair| pair[0] < pair[1]));
    assert_eq!(wrapper.knots(start, start), vec![start]);
    assert!(
        wrapper
            .knots(TruthTime::EPOCH, TruthTime::from_ps(start.as_ps() - 1))
            .is_empty()
    );
    assert!(wrapper.knots(TruthTime::MAX, TruthTime::MIN).is_empty());
    for ps in [-1, inner.span().end.unwrap().as_ps() + 1] {
        assert_eq!(
            wrapper.state(at(ps)).unwrap_err(),
            GeoError::OutsideSpan {
                t_ps: ps,
                start_ps: 0,
                end_ps: inner.span().end.unwrap().as_ps(),
            }
        );
    }
    for duration in [
        TruthDuration::ZERO,
        TruthDuration::from_ps(-1),
        TruthDuration::from_ps(MAX_LOOKBACK.as_ps() + 1),
    ] {
        assert!(Lookback::new(inner.clone(), duration).is_err());
    }
    let near_min = Arc::new(flight(TruthTime::from_ps(i64::MIN + 1)));
    assert!(Lookback::new(near_min, MAX_LOOKBACK).is_err());
    eprintln!(
        "lookback C2: position {position_error:e} m, velocity {velocity_error:e} m/s, acceleration {acceleration_error:e} m/s^2"
    );
}

#[test]
fn solver_failures_and_trajectory_errors_propagate() {
    let rx = ground();
    assert!(matches!(
        solve_receive(&rx, &rx, at(0)),
        Err(GeoError::LightTime(_))
    ));
    assert!(matches!(
        solve_emit(&rx, &rx, at(0)),
        Err(GeoError::LightTime(_))
    ));
    let nonfinite = Static::from_ecef(Ecef::new(f64::NAN, 0.0, 0.0), Quat::IDENTITY);
    assert!(matches!(
        solve_receive(&rx, &nonfinite, at(0)),
        Err(GeoError::NonFinite(_))
    ));
    let fast = InertialLine {
        p0: Vec3::new(7e6, 0.0, 0.0),
        u: Vec3::new(0.5 * C, 0.0, 0.0),
    };
    assert_eq!(
        solve_receive(&rx, &fast, at(0)).unwrap_err(),
        GeoError::NoConvergence {
            what: "light time",
            iterations: LIGHT_TIME_MAX_ITERATIONS
        }
    );
    assert_eq!(
        solve_emit(&fast, &rx, at(0)).unwrap_err(),
        GeoError::NoConvergence {
            what: "light time",
            iterations: LIGHT_TIME_MAX_ITERATIONS
        }
    );
    let tx = Arc::new(flight(TruthTime::EPOCH));
    let outside = solve_receive(&rx, tx.as_ref(), at(0)).unwrap_err();
    let GeoError::OutsideSpan { t_ps, .. } = outside else {
        panic!("expected OutsideSpan: {outside:?}")
    };
    assert_eq!(outside, tx.state(at(t_ps)).unwrap_err());
    let wrapped = Lookback::new(tx, MAX_LOOKBACK).unwrap();
    assert!(solve_receive(&rx, &wrapped, at(0)).is_ok());
}

#[test]
fn the_first_guess_stays_inside_the_moving_endpoints_span() {
    // Review regression: a receiver flight starting at 0 and a transmitter knot at
    // −1 µs whose signal arrives after 0 (`solve_emit`), and a transmitter that ends
    // just after the emission but before the reception (`solve_receive`).
    let rx = flight(TruthTime::EPOCH);
    let tx = equator(0.0);
    let emitted = at(-1_000_000);
    let lt = solve_emit(&rx, &tx, emitted).unwrap();
    assert!(lt.rx.instant > AcqInstant::from_truth_time(TruthTime::EPOCH));
    let back = solve_receive(&rx, &tx, lt.rx.instant).unwrap();
    let d = back.tx.instant.to_q64_64() - emitted.to_q64_64();
    assert!((d as f64 / 18_446_744_073_709_551_616.0).abs() < 1e-3);
    // A transmitter flight ending before the receive instant but after its emission.
    let tx_flight = flight(TruthTime::EPOCH);
    let end = tx_flight.span().end.unwrap();
    let ground = equator(0.0);
    let t_rx = AcqInstant::from_truth_time(TruthTime::from_ps(end.as_ps() + 1_000_000));
    let lt = solve_receive(&ground, &tx_flight, t_rx).unwrap();
    assert!(lt.tx.instant <= AcqInstant::from_truth_time(end));
}
