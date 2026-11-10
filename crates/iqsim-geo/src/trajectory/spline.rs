//! Natural cubic smoothing splines for recorded positions (Reinsch, 1967).
//!
//! Independently on each axis, minimise
//! `Σ_i (y_i - f(t_i))² + λ ∫ f''(t)² dt` over C2 cubic splines with the
//! supplied knots. Times are converted to local seconds on the control plane;
//! callers should supply ECEF metres relative to the first recorded position.
//! The natural end conditions are `γ_0 = γ_{n-1} = 0`, where γ denotes f''.
//!
//! For interval lengths `h_i`, R has diagonal `(h_j + h_{j+1})/3` and
//! off-diagonal `h_{j+1}/6`. Column j of the n × (n-2) matrix Q has entries
//! `1/h_j, -1/h_j - 1/h_{j+1}, 1/h_{j+1}` in rows j, j+1, j+2. Reinsch's
//! equations are `(R + λ QᵀQ) γ = Qᵀy`, followed by `f = y - λ Qγ`.
//! This symmetric positive definite pentadiagonal system uses O(n) storage and
//! a fixed-order banded LDLᵀ factorisation, shared by all three axes.
//!
//! The smoothing time constant τ gives `λ = τ⁴ / Δt̄` (seconds cubed), where
//! `Δt̄ = (t_{n-1} - t_0)/(n-1)`. For uniform spacing Δt, the continuum
//! sinusoidal gain is `1/(1 + λ Δt ω⁴) = 1/(1 + (ωτ)⁴)`; components with
//! `ω > 1/τ` are attenuated. Discrete sampling and natural boundaries perturb
//! this frequency response. τ = 0 gives the interpolating natural spline;
//! two points always give a straight line.
//!
//! An affine trend is removed before solving and restored afterwards, preserving
//! lines without cancellation in Qᵀy. For λ > 1 the system is scaled by 1/λ
//! and the unknown is λγ, avoiding large matrix entries. Neither step changes
//! the minimisation problem. All products and sums remain separate (no FMA).

use std::fmt;

use iqsim_core::TruthTime;

use crate::vec::Vec3;

/// Why a spline could not be fitted.
#[derive(Clone, Debug, PartialEq)]
pub enum SplineError {
    /// Fewer than two points were supplied.
    TooFewPoints(usize),
    /// Timestamp and position counts differ.
    LengthMismatch {
        /// Number of timestamps.
        times: usize,
        /// Number of positions.
        values: usize,
    },
    /// A timestamp is less than or equal to its predecessor.
    NotIncreasing {
        /// Index of the offending timestamp (at least one).
        index: usize,
    },
    /// A value, smoothing constant, or computed floating-point result is NaN/∞.
    NonFinite(&'static str),
    /// The smoothing time constant is negative.
    NegativeSmoothing,
    /// The banded solve encountered a non-positive pivot.
    Singular,
}

impl fmt::Display for SplineError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TooFewPoints(n) => write!(f, "a spline needs at least two points; got {n}"),
            Self::LengthMismatch { times, values } => {
                write!(f, "spline has {times} timestamps but {values} values")
            }
            Self::NotIncreasing { index } => {
                write!(
                    f,
                    "spline timestamp {index} does not follow its predecessor"
                )
            }
            Self::NonFinite(what) => write!(f, "spline {what} is not finite"),
            Self::NegativeSmoothing => write!(f, "spline smoothing time constant is negative"),
            Self::Singular => write!(f, "spline solve encountered a non-positive pivot"),
        }
    }
}

impl std::error::Error for SplineError {}

/// A fitted piecewise cubic in three axes.
#[derive(Clone, Debug, PartialEq)]
pub struct CubicSpline3 {
    /// Knot times (the input timestamps), strictly increasing, n ≥ 2.
    pub knots: Vec<TruthTime>,
    /// One entry per interval `[knots[i], knots[i+1]]`: for each axis (x, y, z),
    /// coefficients [c0, c1, c2, c3] of `f(s) = c0 + c1 s + c2 s² + c3 s³`,
    /// in local seconds `s = (t - knots[i]) / 1e12 ps`.
    pub coeffs: Vec<[[f64; 4]; 3]>,
}

impl CubicSpline3 {
    /// Fit with smoothing time constant `tau_s` in seconds (≥ 0); 0 interpolates.
    ///
    /// Rejects invalid inputs, non-finite intermediate results and non-positive
    /// factorisation pivots. Extremely ill-conditioned spacing or smoothing can
    /// exceed f64 precision even for a mathematically positive definite system.
    pub fn fit(times: &[TruthTime], values: &[Vec3], tau_s: f64) -> Result<Self, SplineError> {
        if times.len() != values.len() {
            return Err(SplineError::LengthMismatch {
                times: times.len(),
                values: values.len(),
            });
        }
        let n = times.len();
        if n < 2 {
            return Err(SplineError::TooFewPoints(n));
        }
        finite(tau_s, "smoothing time constant")?;
        if tau_s < 0.0 {
            return Err(SplineError::NegativeSmoothing);
        }
        for value in values {
            for component in value.to_array() {
                finite(component, "value")?;
            }
        }
        let mut h = Vec::with_capacity(n - 1);
        for (i, pair) in times.windows(2).enumerate() {
            if pair[1] <= pair[0] {
                return Err(SplineError::NotIncreasing { index: i + 1 });
            }
            h.push(seconds_between(pair[1], pair[0]));
        }
        let lambda = if n == 2 {
            0.0
        } else {
            Self::lambda(times, tau_s)
        };
        finite(lambda, "lambda")?;
        let r_scale = if lambda > 1.0 { 1.0 / lambda } else { 1.0 };
        let q_scale = if lambda > 1.0 { 1.0 } else { lambda };
        let inv_h: Vec<_> = h.iter().map(|v| 1.0 / v).collect();
        let q: Vec<_> = inv_h
            .windows(2)
            .map(|v| [v[0], -v[0] - v[1], v[1]])
            .collect();
        let factor = Factor::new(&h, &q, r_scale, q_scale)?;
        let mut coeffs = vec![[[0.0; 4]; 3]; n - 1];
        let span = seconds_between(times[n - 1], times[0]);

        for axis in 0..3 {
            let origin = values[0].to_array()[axis];
            let slope = (values[n - 1].to_array()[axis] - origin) / span;
            finite(slope, "trend")?;
            let trend: Vec<_> = times
                .iter()
                .map(|&t| origin + slope * seconds_between(t, times[0]))
                .collect();
            let residual: Vec<_> = values
                .iter()
                .zip(&trend)
                .map(|(v, t)| v.to_array()[axis] - t)
                .collect();
            let mut interior: Vec<_> = residual
                .windows(3)
                .enumerate()
                .map(|(j, v)| {
                    // Qᵀy as a difference of slopes, so constant data cancel exactly.
                    (v[2] - v[1]) * inv_h[j + 1] - (v[1] - v[0]) * inv_h[j]
                })
                .collect();
            factor.solve(&mut interior);
            // Restore the affine trend algebraically: f = y - λQγ. Starting
            // from y also keeps the interpolating knot values exact at λ = 0.
            let mut fitted: Vec<_> = values.iter().map(|v| v.to_array()[axis]).collect();
            let mut gamma = vec![0.0; n];
            for (j, (&u, column)) in interior.iter().zip(&q).enumerate() {
                finite(u, "solution")?;
                gamma[j + 1] = r_scale * u;
                for (k, &entry) in column.iter().enumerate() {
                    fitted[j + k] -= q_scale * entry * u;
                }
            }
            for (i, interval) in coeffs.iter_mut().enumerate() {
                let c0 = fitted[i];
                let c1 = (fitted[i + 1] - fitted[i]) / h[i]
                    - h[i] * (2.0 * gamma[i] + gamma[i + 1]) / 6.0;
                let c2 = gamma[i] / 2.0;
                let c3 = (gamma[i + 1] - gamma[i]) / (6.0 * h[i]);
                interval[axis] = [c0, c1, c2, c3];
                for &coefficient in &interval[axis] {
                    finite(coefficient, "coefficient")?;
                }
            }
        }
        Ok(Self {
            knots: times.to_vec(),
            coeffs,
        })
    }

    /// Smoothing parameter λ (s³): `tau_s⁴ / Δt̄`, with
    /// `Δt̄ = (times[n-1] - times[0])/(n-1)` in seconds.
    ///
    /// Returns NaN for fewer than two or non-increasing timestamps, negative τ,
    /// or non-finite τ. Overflow of the computed parameter returns infinity;
    /// [`Self::fit`] checks representability before solving.
    pub fn lambda(times: &[TruthTime], tau_s: f64) -> f64 {
        if times.len() < 2
            || !tau_s.is_finite()
            || tau_s < 0.0
            || times.windows(2).any(|p| p[1] <= p[0])
        {
            return f64::NAN;
        }
        let mean = seconds_between(times[times.len() - 1], times[0]) / (times.len() - 1) as f64;
        let square = tau_s * tau_s;
        square * (square / mean)
    }

    /// Evaluate value, first, second and third derivatives (per axis as Vec3)
    /// at local seconds `s` in interval `i`, using Horner with separate operations.
    ///
    /// Derivative units are metres/s, metres/s² and metres/s³ for metre inputs.
    /// `s` is not clamped, so the polynomial also permits extrapolation.
    ///
    /// # Panics
    /// Panics if `i` is outside [`Self::coeffs`].
    pub fn eval(&self, i: usize, s: f64) -> [Vec3; 4] {
        let mut result = [[0.0; 3]; 4];
        for (axis, &[c0, c1, c2, c3]) in self.coeffs[i].iter().enumerate() {
            result[0][axis] = ((c3 * s + c2) * s + c1) * s + c0;
            result[1][axis] = (3.0 * c3 * s + 2.0 * c2) * s + c1;
            result[2][axis] = 6.0 * c3 * s + 2.0 * c2;
            result[3][axis] = 6.0 * c3;
        }
        result.map(|v| Vec3::new(v[0], v[1], v[2]))
    }
}

fn seconds_between(later: TruthTime, earlier: TruthTime) -> f64 {
    // The difference of two i64 timestamps fits in i128, even across the epoch.
    (i128::from(later.as_ps()) - i128::from(earlier.as_ps())) as f64 / 1e12
}

fn finite(value: f64, what: &'static str) -> Result<(), SplineError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(SplineError::NonFinite(what))
    }
}

// Unit lower triangular L has two subdiagonals; d is the diagonal of D.
struct Factor {
    d: Vec<f64>,
    l1: Vec<f64>,
    l2: Vec<f64>,
}

impl Factor {
    fn new(h: &[f64], q: &[[f64; 3]], r: f64, p: f64) -> Result<Self, SplineError> {
        let m = q.len();
        let mut factor = Self {
            d: vec![0.0; m],
            l1: vec![0.0; m],
            l2: vec![0.0; m],
        };
        for (i, column) in q.iter().enumerate() {
            let norm = (column[0] * column[0] + column[1] * column[1]) + column[2] * column[2];
            let mut pivot = r * ((h[i] + h[i + 1]) / 3.0) + p * norm;
            if i >= 2 {
                let a2 = p * q[i - 2][2] * column[0];
                factor.l2[i] = a2 / factor.d[i - 2];
                pivot -= factor.l2[i] * factor.l2[i] * factor.d[i - 2];
            }
            if i >= 1 {
                let dot = q[i - 1][1] * column[0] + q[i - 1][2] * column[1];
                let mut a1 = r * (h[i] / 6.0) + p * dot;
                if i >= 2 {
                    a1 -= factor.l2[i] * factor.d[i - 2] * factor.l1[i - 1];
                }
                factor.l1[i] = a1 / factor.d[i - 1];
                pivot -= factor.l1[i] * factor.l1[i] * factor.d[i - 1];
            }
            finite(pivot, "pivot")?;
            if pivot <= 0.0 {
                return Err(SplineError::Singular);
            }
            factor.d[i] = pivot;
        }
        Ok(factor)
    }

    fn solve(&self, rhs: &mut [f64]) {
        for i in 0..rhs.len() {
            if i >= 2 {
                rhs[i] -= self.l2[i] * rhs[i - 2];
            }
            if i >= 1 {
                rhs[i] -= self.l1[i] * rhs[i - 1];
            }
        }
        for (v, &d) in rhs.iter_mut().zip(&self.d) {
            *v /= d;
        }
        for i in (0..rhs.len()).rev() {
            if i + 2 < rhs.len() {
                rhs[i] -= self.l2[i + 2] * rhs[i + 2];
            }
            if i + 1 < rhs.len() {
                rhs[i] -= self.l1[i + 1] * rhs[i + 1];
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    fn close(a: f64, b: f64, tolerance: f64) {
        assert!(
            (a - b).abs() <= tolerance * (1.0 + a.abs().max(b.abs())),
            "{a} != {b} (tol {tolerance})"
        );
    }

    fn close_vec(a: Vec3, b: Vec3, tolerance: f64) {
        for (a, b) in a.to_array().into_iter().zip(b.to_array()) {
            close(a, b, tolerance);
        }
    }

    fn track() -> impl Strategy<Value = (Vec<TruthTime>, Vec<Vec3>)> {
        prop::collection::vec(
            (
                100_000_000_000_i64..1_000_000_000_000,
                prop::array::uniform3(-100.0_f64..100.0),
            ),
            2..25,
        )
        .prop_map(|rows| {
            let mut ps = -5_000_000_000_000;
            let mut times = Vec::new();
            let mut values = Vec::new();
            for (step, [x, y, z]) in rows {
                ps += step;
                times.push(TruthTime::from_ps(ps));
                values.push(Vec3::new(x, y, z));
            }
            (times, values)
        })
    }

    fn natural(spline: &CubicSpline3) {
        let last = spline.coeffs.len() - 1;
        close_vec(spline.eval(0, 0.0)[2], Vec3::ZERO, 1e-9);
        close_vec(
            spline.eval(
                last,
                seconds_between(spline.knots[last + 1], spline.knots[last]),
            )[2],
            Vec3::ZERO,
            1e-9,
        );
    }

    fn config() -> ProptestConfig {
        let mut config = ProptestConfig::default();
        if std::env::var_os("PROPTEST_CASES").is_none() {
            config.cases = 2048;
        }
        config
    }

    proptest! {
        #![proptest_config(config())]

        #[test]
        fn spline_interpolates((times, values) in track()) {
            let spline = CubicSpline3::fit(&times, &values, 0.0).unwrap();
            assert_eq!(spline.knots, times);
            assert_eq!(spline.coeffs.len(), times.len() - 1);
            for (i, pair) in times.windows(2).enumerate() {
                close_vec(spline.eval(i, 0.0)[0], values[i], 1e-9);
                close_vec(spline.eval(i, seconds_between(pair[1], pair[0]))[0], values[i + 1], 1e-9);
            }
            natural(&spline);
        }

        #[test]
        fn spline_is_c2((times, values) in track(), tau in 0.0_f64..5.0) {
            let spline = CubicSpline3::fit(&times, &values, tau).unwrap();
            for i in 1..times.len() - 1 {
                let left = spline.eval(i - 1, seconds_between(times[i], times[i - 1]));
                let right = spline.eval(i, 0.0);
                for derivative in 0..3 { close_vec(left[derivative], right[derivative], 1e-9); }
            }
            natural(&spline);
        }

        #[test]
        fn spline_reproduces_lines((times, _) in track(), a in prop::array::uniform3(-10.0_f64..10.0), b in prop::array::uniform3(-10.0_f64..10.0), tau in 0.0_f64..100_000.0) {
            let a = Vec3::new(a[0], a[1], a[2]);
            let b = Vec3::new(b[0], b[1], b[2]);
            let values: Vec<_> = times.iter().map(|&t| a + b * seconds_between(t, times[0])).collect();
            let spline = CubicSpline3::fit(&times, &values, tau).unwrap();
            for (i, pair) in times.windows(2).enumerate() {
                let h = seconds_between(pair[1], pair[0]);
                let s = h * 0.37;
                let result = spline.eval(i, s);
                close_vec(result[0], a + b * (seconds_between(pair[0], times[0]) + s), 1e-9);
                close_vec(result[1], b, 1e-9);
                close_vec(result[2], Vec3::ZERO, 1e-9);
                close_vec(result[3], Vec3::ZERO, 1e-9);
            }
        }
    }

    #[test]
    fn spline_large_tau_approaches_least_squares_line() {
        let seconds = [0, 1, 3, 4, 7, 8, 10, 12];
        let times: Vec<_> = seconds
            .iter()
            .map(|t| TruthTime::from_ps(t * 1_000_000_000_000))
            .collect();
        let values: Vec<_> = [3.0, -2.0, 4.0, 1.0, 8.0, -1.0, 2.0, 6.0]
            .iter()
            .map(|&v| Vec3::new(v, -2.0 * v, 0.5 * v))
            .collect();
        let mean_t = seconds.iter().map(|&t| t as f64).sum::<f64>() / times.len() as f64;
        let mean_y =
            values.iter().copied().fold(Vec3::ZERO, |a, b| a + b) * (1.0 / times.len() as f64);
        let mut numerator = Vec3::ZERO;
        let mut denominator = 0.0;
        for (&t, &y) in seconds.iter().zip(&values) {
            let d = t as f64 - mean_t;
            numerator = numerator + (y - mean_y) * d;
            denominator += d * d;
        }
        let slope = numerator * (1.0 / denominator);
        let spline = CubicSpline3::fit(&times, &values, 100.0).unwrap();
        for (i, pair) in times.windows(2).enumerate() {
            let s = seconds_between(pair[1], pair[0]) * 0.43;
            close_vec(
                spline.eval(i, s)[0],
                mean_y + slope * (seconds[i] as f64 + s - mean_t),
                1e-4,
            );
        }
    }

    #[test]
    fn spline_sinusoid_gain() {
        let times: Vec<_> = (0..2001)
            .map(|i| TruthTime::from_ps(i * 100_000_000_000))
            .collect();
        let omega = 1.0;
        let values: Vec<_> = (0..2001)
            .map(|i| {
                let t = i as f64 * 0.1;
                Vec3::new(
                    libm::sin(omega * t),
                    libm::cos(omega * t),
                    -2.0 * libm::sin(omega * t),
                )
            })
            .collect();
        for (omega_tau, expected) in [(1.0, 0.5), (0.5, 1.0 / (1.0 + 1.0 / 16.0))] {
            let spline = CubicSpline3::fit(&times, &values, omega_tau / omega).unwrap();
            let mut numerator = Vec3::ZERO;
            let mut denominator = Vec3::ZERO;
            for (i, &input) in values.iter().enumerate().take(1334).skip(667) {
                let output = spline.eval(i, 0.0)[0];
                numerator = numerator
                    + Vec3::new(output.x * input.x, output.y * input.y, output.z * input.z);
                denominator = denominator
                    + Vec3::new(input.x * input.x, input.y * input.y, input.z * input.z);
            }
            let gains = Vec3::new(
                numerator.x / denominator.x,
                numerator.y / denominator.y,
                numerator.z / denominator.z,
            );
            println!(
                "omega*tau={omega_tau}: measured gains {:?}, expected {expected:.12}",
                gains.to_array()
            );
            for gain in gains.to_array() {
                assert!((gain / expected - 1.0).abs() < 0.02);
            }
        }
    }

    #[test]
    fn spline_three_point_closed_form() {
        let times = [0, 1_000_000_000_000, 2_000_000_000_000].map(TruthTime::from_ps);
        let spline = CubicSpline3::fit(&times, &[Vec3::ZERO, Vec3::X, Vec3::ZERO], 1.0).unwrap();
        // R = 2/3, QᵀQ = 6, Qᵀy = -2; gamma = -0.3.
        close(spline.eval(0, 0.0)[0].x, 0.3, 1e-14);
        close(spline.eval(1, 0.0)[0].x, 0.4, 1e-14);
        close(spline.eval(1, 1.0)[0].x, 0.3, 1e-14);
        close(spline.eval(1, 0.0)[2].x, -0.3, 1e-14);
        natural(&spline);
    }

    #[test]
    fn spline_two_points_and_extreme_timestamps() {
        for times in [[i64::MIN, i64::MAX], [i64::MAX - 1, i64::MAX], [0, 1]] {
            let times = times.map(TruthTime::from_ps);
            let h = seconds_between(times[1], times[0]);
            let spline =
                CubicSpline3::fit(&times, &[Vec3::ZERO, Vec3::new(1.0, -2.0, 3.0)], f64::MAX)
                    .unwrap();
            close_vec(spline.eval(0, h * 0.5)[0], Vec3::new(0.5, -1.0, 1.5), 1e-14);
            assert_eq!(spline.eval(0, h)[2], Vec3::ZERO);
            assert_eq!(spline.eval(0, h)[3], Vec3::ZERO);
        }
    }

    #[test]
    fn spline_errors_and_lambda() {
        let times = [0, 1_000_000_000_000, 3_000_000_000_000].map(TruthTime::from_ps);
        let values = [Vec3::ZERO; 3];
        assert_eq!(
            CubicSpline3::fit(&[], &[], 0.0),
            Err(SplineError::TooFewPoints(0))
        );
        assert_eq!(
            CubicSpline3::fit(&times[..1], &values[..1], 0.0),
            Err(SplineError::TooFewPoints(1))
        );
        assert_eq!(
            CubicSpline3::fit(&times, &values[..2], 0.0),
            Err(SplineError::LengthMismatch {
                times: 3,
                values: 2
            })
        );
        for bad in [[0, 0, 1], [0, 2, 1]] {
            let index = if bad[1] == 0 { 1 } else { 2 };
            assert_eq!(
                CubicSpline3::fit(&bad.map(TruthTime::from_ps), &values, 0.0),
                Err(SplineError::NotIncreasing { index })
            );
        }
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(
                CubicSpline3::fit(&times, &values, bad),
                Err(SplineError::NonFinite(_))
            ));
            for axis in 0..3 {
                let mut v = [0.0; 3];
                v[axis] = bad;
                assert_eq!(
                    CubicSpline3::fit(
                        &times,
                        &[Vec3::ZERO, Vec3::new(v[0], v[1], v[2]), Vec3::ZERO],
                        0.0
                    ),
                    Err(SplineError::NonFinite("value"))
                );
            }
        }
        assert_eq!(
            CubicSpline3::fit(&times, &values, -1.0),
            Err(SplineError::NegativeSmoothing)
        );
        assert_eq!(
            CubicSpline3::fit(&times, &values, f64::MAX),
            Err(SplineError::NonFinite("lambda"))
        );
        close(CubicSpline3::lambda(&times, 2.0), 16.0 / 1.5, 1e-15);
        assert!(CubicSpline3::lambda(&[], 0.0).is_nan());
        assert!(CubicSpline3::lambda(&times, -1.0).is_nan());
        assert!(CubicSpline3::lambda(&[times[1], times[0]], 1.0).is_nan());
        assert!(matches!(
            Factor::new(&[0.0, 0.0], &[[0.0; 3]], 1.0, 0.0),
            Err(SplineError::Singular)
        ));
        let error: &dyn std::error::Error = &SplineError::Singular;
        assert!(!error.to_string().is_empty());
    }

    #[test]
    fn spline_determinism() {
        let times = [
            0,
            100_000_000_003,
            350_000_000_001,
            1_000_000_000_000,
            1_100_000_000_007,
        ]
        .map(TruthTime::from_ps);
        let values = [
            Vec3::ZERO,
            Vec3::new(1.3, -0.7, 0.3),
            Vec3::new(1.0, 2.0, -3.0),
            Vec3::new(-0.2, 3.1, 0.8),
            Vec3::X,
        ];
        for tau in [0.0, 0.4, 20.0] {
            let a = CubicSpline3::fit(&times, &values, tau).unwrap();
            let b = CubicSpline3::fit(&times, &values, tau).unwrap();
            assert_eq!(a.knots, b.knots);
            for (a, b) in a
                .coeffs
                .iter()
                .flatten()
                .flatten()
                .zip(b.coeffs.iter().flatten().flatten())
            {
                assert_eq!(a.to_bits(), b.to_bits());
            }
        }
    }
}
