//! Band-limited pulse edge tables (contracts §5.3).
//!
//! A pulse edge is the unit-height ramp e(t) of the chosen [`EdgeShape`] (0 for t < 0,
//! 1 after the rise time). The table for an edge, a tap response H and a carrier
//! offset f_off from the channel centre is
//!
//! ```text
//! G(τ) = ∫ e(τ − x)·h(x)·exp(−j2π·f_off·x) dx        (τ in seconds)
//! ```
//!
//! the edge filtered by H shifted to the emitter's offset. The renderer multiplies it by
//! the path's carrier rotation, so a pulse arrives band-limited exactly as if its
//! RF signal had passed through H: `(p·e^{j2πf_off t}) ⊛ h = e^{j2πf_off t}·(p ⊛ h_s)`.
//! G is 0 before the truncated impulse response starts and constant (the plateau
//! ∫h·e^{−j2πf_off x} ≈ H(f_off)) after the edge and the response have both ended.
//!
//! # Construction (deterministic: f64, `libm`, fixed order)
//!
//! Points lie on a grid of 1/R target samples, R = 2^log2_res. Every grid shares the
//! origin `u_lo = ⌊(τ_g − T_h)·F_s·64⌋/64` (the base grid of 1/64 sample), so
//! a finer grid contains every coarser point.
//!
//! 1. Per grid cell [τ_j, τ_{j+1}], 8-point Gauss–Legendre on enough equal panels that
//!    none spans more than half a cycle of |f_off| + 2f_c gives
//!    A_j = ∫g and B_j = ∫(x − τ_g)·g, g(x) = h(x)·e^{−j2πf_off x}. Prefix sums give
//!    S(τ_i) = ΣA and the first moment ΣB.
//! 2. For a linear ramp of duration T_r, with τ_i − T_r = τ_{i−m} − φ·Δ:
//!    `G_i = S(τ_i − T_r) + (1/T_r)·∫_{τ_i−T_r}^{τ_i} (τ_i − x)·g(x) dx`.
//!    The integral is a difference of prefix sums plus one Gauss–Legendre partial cell
//!    of width φ·Δ. A step edge (T_r = 0) is G = S.
//! 3. The plateau is the analytic H(f_off) rounded to f32. Every table of the same
//!    response and offset stores that same value as its last point, so a pulse's rise
//!    and fall cancel exactly after it ends. The quadrature's own end value must
//!    match it to [`PLATEAU_ERROR_BOUND`] (10⁻⁶), which also checks that the quadrature
//!    has converged.
//!
//! # Resolution and error bound
//!
//! The renderer interpolates linearly. Starting at R = 64, the builder also computes the
//! table at 2R and measures the linear-interpolation error at every cell midpoint. It
//! doubles R until that error is ≤ [`INTERP_ERROR_BOUND`] (5·10⁻⁶ of the unit edge
//! height), up to R = 2¹², and fails otherwise. With the truncation tail (< 10⁻⁷) and
//! f32 storage (< 6·10⁻⁸ per component), the table is within 10⁻⁵ of the exact
//! band-limited edge. Scaled by the path amplitude, that is the contracts' 10⁻⁵·V_ref
//! absolute bound. The midpoint is where linear interpolation error peaks for a
//! function whose curvature is nearly constant over a cell, so this is a measured
//! check, not a proof between points.
//!
//! # Keys and cache preparation
//!
//! A table depends only on [`EdgeTableKey`]: the response identity (which covers tap
//! kind, response, output rate and centre), the edge shape and the exact offset in µHz.
//! Tables are immutable. A runtime frequency or response change produces a new key, and
//! its table is built while the work list for the change's effective block is
//! constructed, before that block can render. Old tables stay valid for blocks already
//! built.

use std::fmt;

use iqsim_core::{TapResponse, TapResponseId};

/// Base grid exponent: 64 points per target sample.
pub const BASE_LOG2_RES: u32 = 6;
/// Largest grid exponent: 4096 points per target sample.
pub const MAX_LOG2_RES: u32 = 12;
/// Required linear-interpolation error, in units of the unit edge height.
pub const INTERP_ERROR_BOUND: f64 = 5e-6;
/// Largest accepted difference between the quadrature's end value and the analytic
/// plateau H(f_off).
pub const PLATEAU_ERROR_BOUND: f64 = 1e-6;
/// Longest accepted rise time, ps (1 ms).
pub const MAX_RISE_PS: i64 = 1_000_000_000;

/// The shape of a 0 → 1 pulse edge.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum EdgeShape {
    /// A linear ramp lasting `rise_ps` picoseconds. 0 is an ideal step.
    Linear {
        /// Rise (or fall) duration, ps.
        rise_ps: i64,
    },
}

/// Everything an edge table depends on.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct EdgeTableKey {
    /// Response identity.
    pub response: TapResponseId,
    /// Edge shape.
    pub shape: EdgeShape,
    /// Emitter carrier offset from the channel centre, µHz.
    pub offset_uhz: i64,
}

/// A built table, ready to be copied into a work-list arena.
#[derive(Clone, Debug, PartialEq)]
pub struct BuiltEdgeTable {
    /// Coordinate of point 0, Q32.32 target samples.
    pub u_lo: i64,
    /// Points per target sample = 2^log2_res.
    pub log2_res: u32,
    /// Points (re, im); the last equals the shared plateau.
    pub points: Vec<[f32; 2]>,
    /// Measured worst linear-interpolation error at cell midpoints.
    pub interp_error: f64,
}

/// A table that cannot be built.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum EdgeError {
    /// The rise time is negative or above [`MAX_RISE_PS`].
    Rise,
    /// The interpolation bound was not met at the finest resolution.
    Resolution {
        /// Error achieved at R = 2^[`MAX_LOG2_RES`].
        achieved: f64,
    },
    /// The quadrature's end value does not match the analytic plateau H(f_off) to
    /// [`PLATEAU_ERROR_BOUND`] at the finest resolution.
    Quadrature {
        /// Difference achieved at the finest resolution.
        achieved: f64,
    },
    /// The table extends beyond the coordinate range.
    Range,
}

impl fmt::Display for EdgeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Rise => f.write_str("edge rise time must be 0 to 1 ms"),
            Self::Resolution { achieved } => write!(
                f,
                "edge table interpolation error {achieved:.3e} exceeds {INTERP_ERROR_BOUND:.0e} \
                 at the finest resolution"
            ),
            Self::Quadrature { achieved } => write!(
                f,
                "edge table quadrature misses the analytic plateau by {achieved:.3e} \
                 (bound {PLATEAU_ERROR_BOUND:.0e})"
            ),
            Self::Range => f.write_str("edge table exceeds the coordinate range"),
        }
    }
}

impl std::error::Error for EdgeError {}

/// 8-point Gauss–Legendre nodes and weights on [−1, 1].
const GL_NODES: [f64; 4] = [
    0.183_434_642_495_649_8,
    0.525_532_409_916_329,
    0.796_666_477_413_626_7,
    0.960_289_856_497_536_3,
];
const GL_WEIGHTS: [f64; 4] = [
    0.362_683_783_378_362,
    0.313_706_645_877_887_3,
    0.222_381_034_453_374_5,
    0.101_228_536_290_376_3,
];

#[derive(Clone, Copy, Default)]
struct C64 {
    re: f64,
    im: f64,
}

/// The integrand g(x) = h(x)·e^{−j2πf_off x}.
struct Integrand<'a> {
    resp: &'a TapResponse,
    offset_hz: f64,
    center: f64,
}

impl Integrand<'_> {
    fn g(&self, x: f64) -> C64 {
        let h = self.resp.impulse(x);
        if h == 0.0 {
            return C64::default();
        }
        let (s, c) = libm::sincos(2.0 * core::f64::consts::PI * self.offset_hz * x);
        C64 {
            re: h * c,
            im: -(h * s),
        }
    }

    /// Gauss–Legendre panels per cell of width `dt`: enough that no panel spans more
    /// than half a cycle of the integrand's highest frequency, |f_off| + 2·f_c, so
    /// a large carrier offset is never undersampled.
    fn panels(&self, dt: f64) -> u32 {
        let fmax = self.offset_hz.abs() + 2.0 * self.resp.cutoff_hz();
        let p = libm::ceil(dt * fmax * 2.0);
        if p < 1.0 { 1 } else { p as u32 }
    }

    /// (∫g, ∫(x − center)·g) over [a, b]: `panels` equal panels, each by 8-point
    /// Gauss–Legendre, summed in a fixed order.
    fn cell(&self, a: f64, b: f64, panels: u32) -> (C64, C64) {
        let width = (b - a) / f64::from(panels);
        let half = 0.5 * width;
        let (mut s0, mut s1) = (C64::default(), C64::default());
        for p in 0..panels {
            let mid = a + (f64::from(p) + 0.5) * width;
            for i in 0..4 {
                for sign in [-1.0, 1.0] {
                    let x = mid + sign * GL_NODES[i] * half;
                    let v = self.g(x);
                    let w = GL_WEIGHTS[i];
                    let xt = x - self.center;
                    s0.re += w * v.re;
                    s0.im += w * v.im;
                    s1.re += w * xt * v.re;
                    s1.im += w * xt * v.im;
                }
            }
        }
        (
            C64 {
                re: s0.re * half,
                im: s0.im * half,
            },
            C64 {
                re: s1.re * half,
                im: s1.im * half,
            },
        )
    }
}

/// The f64 table at one resolution.
struct Grid {
    values: Vec<C64>,
    plateau: C64,
}

/// Grid geometry shared by every resolution: base-grid start and end (1/64 sample).
fn base_span(resp: &TapResponse, rise_s: f64, fs_hz: f64) -> Result<(i64, i64), EdgeError> {
    let (x_min, x_max) = resp.support_s();
    let base = f64::from(1u32 << BASE_LOG2_RES);
    let lo = libm::floor(x_min * fs_hz * base);
    let hi = libm::ceil((x_max + rise_s) * fs_hz * base);
    // |coordinate| < 2³⁰ samples, on the 1/64 base grid.
    let limit = (1u64 << (30 + BASE_LOG2_RES)) as f64;
    if !(lo.abs() < limit && hi.abs() < limit) {
        return Err(EdgeError::Range);
    }
    Ok((lo as i64, hi as i64))
}

fn build_grid(f: &Integrand<'_>, rise_s: f64, fs_hz: f64, span: (i64, i64), extra: u32) -> Grid {
    let scale = 1i64 << extra;
    let first = span.0 * scale;
    let cells = ((span.1 - span.0) * scale) as usize;
    let dt = 1.0 / (fs_hz * f64::from(1u32 << (BASE_LOG2_RES + extra)));
    let tau = |i: usize| (first + i as i64) as f64 * dt;
    let panels = f.panels(dt);

    let mut ca = Vec::with_capacity(cells + 1);
    let mut cb = Vec::with_capacity(cells + 1);
    let (mut sa, mut sb) = (C64::default(), C64::default());
    ca.push(sa);
    cb.push(sb);
    for j in 0..cells {
        let (a, b) = f.cell(tau(j), tau(j + 1), panels);
        sa.re += a.re;
        sa.im += a.im;
        sb.re += b.re;
        sb.im += b.im;
        ca.push(sa);
        cb.push(sb);
    }
    let plateau = ca[cells];
    if rise_s == 0.0 {
        return Grid {
            values: ca,
            plateau,
        };
    }
    let rho = rise_s / dt;
    let m = libm::floor(rho);
    let phi = rho - m;
    let m = m as usize;
    let mut values = Vec::with_capacity(cells + 1);
    for i in 0..=cells {
        let ti = tau(i) - f.center;
        let g = if i < m {
            // The ramp window starts before the grid, where g = 0.
            C64 {
                re: (ti * ca[i].re - cb[i].re) / rise_s,
                im: (ti * ca[i].im - cb[i].im) / rise_s,
            }
        } else {
            let j = i - m;
            let (pa, pb) = if phi == 0.0 {
                (C64::default(), C64::default())
            } else {
                f.cell(tau(j) - phi * dt, tau(j), panels)
            };
            let s_lag = C64 {
                re: ca[j].re - pa.re,
                im: ca[j].im - pa.im,
            };
            let wa = C64 {
                re: ca[i].re - ca[j].re + pa.re,
                im: ca[i].im - ca[j].im + pa.im,
            };
            let wb = C64 {
                re: cb[i].re - cb[j].re + pb.re,
                im: cb[i].im - cb[j].im + pb.im,
            };
            C64 {
                re: s_lag.re + (ti * wa.re - wb.re) / rise_s,
                im: s_lag.im + (ti * wa.im - wb.im) / rise_s,
            }
        };
        values.push(g);
    }
    Grid { values, plateau }
}

/// Builds the table for `shape` through `resp` at carrier offset `offset_uhz`.
pub fn build_edge_table(
    resp: &TapResponse,
    shape: EdgeShape,
    offset_uhz: i64,
) -> Result<BuiltEdgeTable, EdgeError> {
    let EdgeShape::Linear { rise_ps } = shape;
    if !(0..=MAX_RISE_PS).contains(&rise_ps) {
        return Err(EdgeError::Rise);
    }
    let rise_s = rise_ps as f64 * 1e-12;
    let rate = resp.rate();
    let fs_hz = rate.num() as f64 / rate.den() as f64;
    let f = Integrand {
        resp,
        offset_hz: offset_uhz as f64 * 1e-6,
        center: resp.group_delay_s(),
    };
    let span = base_span(resp, rise_s, fs_hz)?;
    let mut coarse = build_grid(&f, rise_s, fs_hz, span, 0);
    // The shared plateau: the analytic H(f_off), rounded to f32. The quadrature's own
    // end value must agree with it (truncation < 10⁻⁷ plus quadrature error).
    let (pr, pi) = resp.response(f.offset_hz);
    let plateau = [pr as f32, pi as f32];
    let mut extra = 0;
    loop {
        let fine = build_grid(&f, rise_s, fs_hz, span, extra + 1);
        let mut err = 0.0_f64;
        for i in 0..coarse.values.len() - 1 {
            let (a, b, m) = (
                coarse.values[i],
                coarse.values[i + 1],
                fine.values[2 * i + 1],
            );
            let dr = m.re - 0.5 * (a.re + b.re);
            let di = m.im - 0.5 * (a.im + b.im);
            err = err.max((dr * dr + di * di).sqrt());
        }
        let (dr, di) = (coarse.plateau.re - pr, coarse.plateau.im - pi);
        let plateau_err = (dr * dr + di * di).sqrt();
        if err <= INTERP_ERROR_BOUND && plateau_err <= PLATEAU_ERROR_BOUND {
            let mut points: Vec<[f32; 2]> = coarse
                .values
                .iter()
                .map(|v| [v.re as f32, v.im as f32])
                .collect();
            if let Some(last) = points.last_mut() {
                *last = plateau;
            }
            let log2_res = BASE_LOG2_RES + extra;
            return Ok(BuiltEdgeTable {
                u_lo: span.0 << (32 - BASE_LOG2_RES),
                log2_res,
                points,
                interp_error: err,
            });
        }
        if BASE_LOG2_RES + extra >= MAX_LOG2_RES {
            return Err(if plateau_err > PLATEAU_ERROR_BOUND {
                EdgeError::Quadrature {
                    achieved: plateau_err,
                }
            } else {
                EdgeError::Resolution { achieved: err }
            });
        }
        extra += 1;
        coarse = fine;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use iqsim_core::{ResponseSpec, SampleRate, TapKind};

    fn resp(gd_ps: i64) -> TapResponse {
        TapResponse::new(
            TapKind::VirtualReceiver,
            ResponseSpec {
                passband_edge_uhz: 6_000_000_000_000,
                stopband_edge_uhz: 7_500_000_000_000,
                stopband_atten_mdb: 80_000,
                group_delay_ps: gd_ps,
            },
            SampleRate::from_hz(16_000_000).unwrap(),
            1_090_000_000_000_000,
        )
        .unwrap()
    }

    #[test]
    fn plateau_is_the_response_at_the_offset() {
        let r = resp(1_000_000);
        for off in [0_i64, 3_000_000_000_000, -6_500_000_000_000] {
            let t = build_edge_table(&r, EdgeShape::Linear { rise_ps: 100_000 }, off).unwrap();
            let p = t.points.last().unwrap();
            let (hr, hi) = r.response(off as f64 * 1e-6);
            let (dr, di) = (f64::from(p[0]) - hr, f64::from(p[1]) - hi);
            let err = (dr * dr + di * di).sqrt();
            assert!(err < 3e-7, "offset {off}: plateau error {err}");
            assert!(t.interp_error <= INTERP_ERROR_BOUND);
            assert_eq!(t.points[0], [0.0, 0.0]);
        }
    }

    #[test]
    fn rise_and_fall_tables_share_the_plateau() {
        let r = resp(0);
        let a = build_edge_table(&r, EdgeShape::Linear { rise_ps: 50_000 }, 1_000_000).unwrap();
        let b = build_edge_table(&r, EdgeShape::Linear { rise_ps: 0 }, 1_000_000).unwrap();
        let c = build_edge_table(&r, EdgeShape::Linear { rise_ps: 200_000 }, 1_000_000).unwrap();
        assert_eq!(a.points.last(), b.points.last());
        assert_eq!(a.points.last(), c.points.last());
        assert_eq!(a.u_lo, b.u_lo);
    }

    #[test]
    fn step_is_symmetric_about_the_group_delay() {
        // h is real and even about τ_g, so a step at f_off = 0 crosses ½ at τ_g and
        // G(τ_g + x) + G(τ_g − x) = 1 (up to the plateau's truncation).
        let r = resp(1_000_000);
        let t = build_edge_table(&r, EdgeShape::Linear { rise_ps: 0 }, 0).unwrap();
        let res = 1i64 << t.log2_res;
        let center = 16 * res - (t.u_lo >> (32 - t.log2_res));
        let c = center as usize;
        assert!((f64::from(t.points[c][0]) - 0.5).abs() < 1e-6);
        for d in 1..(c.min(t.points.len() - 1 - c)) {
            let s = f64::from(t.points[c - d][0]) + f64::from(t.points[c + d][0]);
            assert!((s - 1.0).abs() < 2e-6, "d = {d}: {s}");
            assert!(t.points[c + d][1].abs() < 1e-6);
        }
    }

    #[test]
    fn far_offset_plateau_is_not_aliased_by_the_quadrature() {
        // Regression (prompt 5 review): 1 MS/s, 300/450 kHz, 80 dB, 512 MHz offset,
        // 1 µs rise. One GL rule per 15.6 ns cell undersampled the 512 MHz carrier and
        // gave a plateau of 0.355 instead of ≈ 0.
        let r = TapResponse::new(
            TapKind::VirtualReceiver,
            ResponseSpec {
                passband_edge_uhz: 300_000_000_000,
                stopband_edge_uhz: 450_000_000_000,
                stopband_atten_mdb: 80_000,
                group_delay_ps: 0,
            },
            SampleRate::from_hz(1_000_000).unwrap(),
            0,
        )
        .unwrap();
        let t = build_edge_table(
            &r,
            EdgeShape::Linear { rise_ps: 1_000_000 },
            512_000_000_000_000,
        )
        .unwrap();
        let peak = t
            .points
            .iter()
            .map(|p| (f64::from(p[0]) * f64::from(p[0]) + f64::from(p[1]) * f64::from(p[1])).sqrt())
            .fold(0.0, f64::max);
        assert!(peak < 1e-5, "peak {peak}");
        assert_eq!(t.points.last(), Some(&[0.0, 0.0]));
    }

    #[test]
    fn bad_rise_is_rejected() {
        let r = resp(0);
        assert_eq!(
            build_edge_table(&r, EdgeShape::Linear { rise_ps: -1 }, 0),
            Err(EdgeError::Rise)
        );
    }
}
