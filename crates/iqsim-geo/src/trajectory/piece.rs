//! Piecewise quintic kinematics in ECEF (contracts §13.2).
//!
//! A [`Path`] is a sequence of [`Piece`]s with knots at whole picoseconds. Each piece
//! holds, per ECEF axis, the coefficients of `p(τ) = c₀ + c₁τ + … + c₅τ⁵` in local
//! seconds τ from its start. Velocity, acceleration and jerk are its derivatives, so
//! they are always consistent with position. [`Piece::hermite`] builds the unique
//! quintic matching position, velocity and acceleration at both ends, so consecutive
//! pieces built from shared knot values are C2 (to f64 rounding).

use iqsim_core::{AcqInstant, TruthTime};

use crate::GeoError;
use crate::trajectory::seconds_since;
use crate::vec::Vec3;

/// Pieces shorter than this (seconds) take their displacement from the end
/// velocities and accelerations ([`Piece::hermite`]).
pub const SHORT_PIECE_S: f64 = 0.1;

/// Position, velocity, acceleration and jerk at one instant (ECEF, SI units).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Kinematics {
    /// Position, metres.
    pub p: Vec3,
    /// Velocity, m/s.
    pub v: Vec3,
    /// Acceleration, m/s².
    pub a: Vec3,
    /// Jerk, m/s³.
    pub j: Vec3,
}

/// One polynomial piece.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Piece {
    /// First instant of the piece.
    pub start: TruthTime,
    /// Last instant of the piece (the next piece's start, or the path end).
    pub end: TruthTime,
    /// Per axis (x, y, z), coefficients c₀…c₅ in local seconds.
    pub coeffs: [[f64; 6]; 3],
}

impl Piece {
    /// The quintic from `start` to `end` with position, velocity and acceleration
    /// `(p0, v0, a0)` at the start and `(p1, v1, a1)` at the end. `end` must be after
    /// `start`.
    ///
    /// For pieces shorter than [`SHORT_PIECE_S`] the displacement is taken from the
    /// end velocities and accelerations, `Δp = T(v₀ + v₁)/2 + T²(a₀ − a₁)/12` (exact
    /// for quartic motion), instead of `p₁ − p₀`. The rounding of Earth-sized
    /// positions (~10⁻⁹ m) would otherwise be divided by T² and T³. The piece then
    /// ends within that rounding (plus a T⁵ truncation term, < 10⁻¹¹ m) of `p₁`, and
    /// velocity and acceleration still match exactly.
    pub fn hermite(
        start: TruthTime,
        end: TruthTime,
        k0: [Vec3; 3],
        k1: [Vec3; 3],
    ) -> Result<Self, GeoError> {
        let span = i128::from(end.as_ps()) - i128::from(start.as_ps());
        if span <= 0 {
            return Err(GeoError::Trajectory(
                "a trajectory piece must have a positive length".to_owned(),
            ));
        }
        let t = span as f64 / 1.0e12;
        let (t2, t3) = (t * t, t * t * t);
        let (t4, t5) = (t3 * t, t3 * t2);
        let short = t < SHORT_PIECE_S;
        let axis = |f: fn(Vec3) -> f64| {
            let (p0, v0, a0) = (f(k0[0]), f(k0[1]), f(k0[2]));
            let (p1, v1, a1) = (f(k1[0]), f(k1[1]), f(k1[2]));
            let dp = if short {
                t * (0.5 * (v0 + v1)) + t2 * ((a0 - a1) / 12.0)
            } else {
                p1 - p0
            };
            let c3 = (20.0 * dp - (8.0 * v1 + 12.0 * v0) * t - (3.0 * a0 - a1) * t2) / (2.0 * t3);
            let c4 = (-30.0 * dp + (14.0 * v1 + 16.0 * v0) * t + (3.0 * a0 - 2.0 * a1) * t2)
                / (2.0 * t4);
            let c5 = (12.0 * dp - 6.0 * (v1 + v0) * t + (a1 - a0) * t2) / (2.0 * t5);
            [p0, v0, 0.5 * a0, c3, c4, c5]
        };
        Ok(Self {
            start,
            end,
            coeffs: [axis(|v| v.x), axis(|v| v.y), axis(|v| v.z)],
        })
    }

    /// A piece from per-axis cubic coefficients (track splines).
    pub fn cubic(start: TruthTime, end: TruthTime, origin: Vec3, c: [[f64; 4]; 3]) -> Self {
        let o = origin.to_array();
        let axis = |i: usize| [o[i] + c[i][0], c[i][1], c[i][2], c[i][3], 0.0, 0.0];
        Self {
            start,
            end,
            coeffs: [axis(0), axis(1), axis(2)],
        }
    }

    /// Kinematics at local time `tau` seconds (Horner, separate multiplies and adds).
    pub fn eval(&self, tau: f64) -> Kinematics {
        let mut out = [[0.0; 3]; 4];
        for (i, c) in self.coeffs.iter().enumerate() {
            let p = c[0] + tau * (c[1] + tau * (c[2] + tau * (c[3] + tau * (c[4] + tau * c[5]))));
            let v = c[1]
                + tau * (2.0 * c[2] + tau * (3.0 * c[3] + tau * (4.0 * c[4] + tau * (5.0 * c[5]))));
            let a = 2.0 * c[2] + tau * (6.0 * c[3] + tau * (12.0 * c[4] + tau * (20.0 * c[5])));
            let j = 6.0 * c[3] + tau * (24.0 * c[4] + tau * (60.0 * c[5]));
            out[0][i] = p;
            out[1][i] = v;
            out[2][i] = a;
            out[3][i] = j;
        }
        let v = |a: [f64; 3]| Vec3::new(a[0], a[1], a[2]);
        Kinematics {
            p: v(out[0]),
            v: v(out[1]),
            a: v(out[2]),
            j: v(out[3]),
        }
    }
}

/// Consecutive pieces; each piece's `end` is the next one's `start`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Path {
    pieces: Vec<Piece>,
}

impl Path {
    /// The pieces, in time order.
    pub fn pieces(&self) -> &[Piece] {
        &self.pieces
    }

    /// First instant, if any piece exists.
    pub fn start(&self) -> Option<TruthTime> {
        self.pieces.first().map(|p| p.start)
    }

    /// Last instant, if any piece exists.
    pub fn end(&self) -> Option<TruthTime> {
        self.pieces.last().map(|p| p.end)
    }

    /// Appends a piece that starts where the path ends.
    pub(crate) fn push(&mut self, piece: Piece) {
        debug_assert!(self.end().is_none_or(|e| e == piece.start));
        debug_assert!(piece.end > piece.start);
        self.pieces.push(piece);
    }

    /// Kinematics at `t`. A knot belongs to the piece that starts there, except the
    /// path end, which belongs to the last piece.
    pub fn eval(&self, t: AcqInstant) -> Result<Kinematics, GeoError> {
        let (Some(start), Some(end)) = (self.start(), self.end()) else {
            return Err(GeoError::Trajectory("the trajectory is empty".to_owned()));
        };
        if t < AcqInstant::from_truth_time(start) || t > AcqInstant::from_truth_time(end) {
            return Err(GeoError::OutsideSpan {
                t_ps: t.ps().as_ps(),
                start_ps: start.as_ps(),
                end_ps: end.as_ps(),
            });
        }
        let i = self
            .pieces
            .partition_point(|p| AcqInstant::from_truth_time(p.start) <= t)
            .saturating_sub(1);
        let piece = &self.pieces[i];
        Ok(piece.eval(seconds_since(t, piece.start)))
    }

    /// Drops everything after `t` (inclusive of pieces starting at `t`) and shortens
    /// the piece containing `t` to end there, without changing its coefficients, so
    /// every instant before `t` evaluates to the same bits as before.
    pub(crate) fn truncate_at(&mut self, t: TruthTime) {
        self.pieces.retain(|p| p.start < t);
        if let Some(last) = self.pieces.last_mut()
            && last.end > t
        {
            last.end = t;
        }
    }

    /// Piece starts and the path end within `[from, to]`.
    pub fn knots(&self, from: TruthTime, to: TruthTime) -> Vec<TruthTime> {
        self.pieces
            .iter()
            .map(|p| p.start)
            .chain(self.end())
            .filter(|&k| from <= k && k <= to)
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hermite_matches_its_end_conditions() {
        let k0 = [
            Vec3::new(6.4e6, 1.0, -2.0),
            Vec3::new(250.0, -3.0, 1.0),
            Vec3::new(0.1, 2.0, -0.5),
        ];
        let k1 = [
            Vec3::new(6.4e6 + 251.0, -1.5, 0.0),
            Vec3::new(252.0, -1.0, 0.0),
            Vec3::new(0.0, 1.5, 0.3),
        ];
        let p = Piece::hermite(
            TruthTime::EPOCH,
            TruthTime::from_ps(1_000_000_000_000),
            k0,
            k1,
        )
        .unwrap();
        let a = p.eval(0.0);
        assert_eq!((a.p, a.v, a.a), (k0[0], k0[1], k0[2]));
        let b = p.eval(1.0);
        assert!(b.p.distance(k1[0]) < 1e-8);
        assert!(b.v.distance(k1[1]) < 1e-11);
        assert!(b.a.distance(k1[2]) < 1e-11);
        // Jerk is the derivative of acceleration.
        let h = 1e-6;
        let fd = (p.eval(0.5 + h).a - p.eval(0.5 - h).a) * (0.5 / h);
        assert!(fd.distance(p.eval(0.5).j) < 1e-6);
    }

    #[test]
    fn short_pieces_stay_well_conditioned() {
        // Straight flight at 100 m/s through an Earth-sized position: a 1 ms piece
        // whose endpoint positions carry their f64 rounding.
        let p0 = Vec3::new(4_023_305.760_974_29, -70_161.256_903, 4_935_946.151_055);
        let v = Vec3::new(-60.0, 1.0, 79.0);
        let a = Vec3::new(-0.001, 0.0, -0.0012);
        let t = 1.0e-3;
        let p1 = p0 + v * t + a * (0.5 * t * t);
        let k0 = [p0, v, a];
        let k1 = [p1, v + a * t, a];
        let p =
            Piece::hermite(TruthTime::EPOCH, TruthTime::from_ps(1_000_000_000), k0, k1).unwrap();
        for tau in [0.0, 2.5e-4, 5e-4, 1e-3] {
            assert!(p.eval(tau).a.distance(a) < 1e-9, "{:?}", p.eval(tau).a);
        }
        assert!(p.eval(t).p.distance(p1) < 1e-8);
        assert!(Piece::hermite(TruthTime::EPOCH, TruthTime::EPOCH, k0, k1).is_err());
        assert!(Piece::hermite(TruthTime::MAX, TruthTime::MIN, k0, k1).is_err());
    }

    #[test]
    fn truncation_keeps_the_past() {
        let k = [
            Vec3::new(1.0, 2.0, 3.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::ZERO,
        ];
        let k2 = [
            Vec3::new(2.0, 2.0, 3.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::ZERO,
        ];
        let mut path = Path::default();
        let s = 1_000_000_000_000;
        path.push(Piece::hermite(TruthTime::EPOCH, TruthTime::from_ps(s), k, k2).unwrap());
        path.push(Piece::hermite(TruthTime::from_ps(s), TruthTime::from_ps(2 * s), k2, k).unwrap());
        let t = AcqInstant::from_truth_time(TruthTime::from_ps(s / 3));
        let before = path.eval(t).unwrap();
        path.truncate_at(TruthTime::from_ps(s / 2));
        assert_eq!(path.eval(t).unwrap(), before);
        assert_eq!(path.end(), Some(TruthTime::from_ps(s / 2)));
        assert!(
            path.eval(AcqInstant::from_truth_time(TruthTime::from_ps(s)))
                .is_err()
        );
        assert_eq!(
            path.knots(TruthTime::MIN, TruthTime::MAX),
            vec![TruthTime::EPOCH, TruthTime::from_ps(s / 2)]
        );
    }
}
