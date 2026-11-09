//! Three-vectors and 3×3 matrices in f64.
//!
//! Every operation is written out with a fixed order of separate multiplies and adds
//! (no FMA contraction, AGENTS.md), so results are the same bits on every machine.
//! Sums of three products are evaluated left to right: `(a₀b₀ + a₁b₁) + a₂b₂`.

use core::ops::{Add, Mul, Neg, Sub};

/// A 3-vector of f64. The frame is the caller's: the type carries no frame tag.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Vec3 {
    /// First component.
    pub x: f64,
    /// Second component.
    pub y: f64,
    /// Third component.
    pub z: f64,
}

impl Vec3 {
    /// The zero vector.
    pub const ZERO: Self = Self::new(0.0, 0.0, 0.0);
    /// Unit x.
    pub const X: Self = Self::new(1.0, 0.0, 0.0);
    /// Unit y.
    pub const Y: Self = Self::new(0.0, 1.0, 0.0);
    /// Unit z.
    pub const Z: Self = Self::new(0.0, 0.0, 1.0);

    /// `(x, y, z)`.
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// The components as an array.
    pub const fn to_array(self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }

    /// Dot product `(x·x' + y·y') + z·z'`.
    pub fn dot(self, o: Self) -> f64 {
        let s = self.x * o.x + self.y * o.y;
        s + self.z * o.z
    }

    /// Cross product `self × o`.
    pub fn cross(self, o: Self) -> Self {
        Self::new(
            self.y * o.z - self.z * o.y,
            self.z * o.x - self.x * o.z,
            self.x * o.y - self.y * o.x,
        )
    }

    /// Euclidean norm, `sqrt(self·self)` (IEEE sqrt).
    pub fn norm(self) -> f64 {
        self.dot(self).sqrt()
    }

    /// `self / |self|`. A zero vector stays zero.
    pub fn normalized(self) -> Self {
        let n = self.norm();
        if n == 0.0 { self } else { self * (1.0 / n) }
    }

    /// Distance to `o`, `|o − self|`.
    pub fn distance(self, o: Self) -> f64 {
        (o - self).norm()
    }
}

impl Add for Vec3 {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self::new(self.x + o.x, self.y + o.y, self.z + o.z)
    }
}

impl Sub for Vec3 {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Self::new(self.x - o.x, self.y - o.y, self.z - o.z)
    }
}

impl Neg for Vec3 {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(-self.x, -self.y, -self.z)
    }
}

impl Mul<f64> for Vec3 {
    type Output = Self;
    fn mul(self, s: f64) -> Self {
        Self::new(self.x * s, self.y * s, self.z * s)
    }
}

/// A 3×3 matrix, row-major: `m[row][col]`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Mat3(pub [[f64; 3]; 3]);

impl Mat3 {
    /// The identity.
    pub const IDENTITY: Self = Self([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);

    /// The matrix whose rows are `r0`, `r1`, `r2`.
    pub const fn from_rows(r0: Vec3, r1: Vec3, r2: Vec3) -> Self {
        Self([r0.to_array(), r1.to_array(), r2.to_array()])
    }

    /// Row `i` as a vector.
    pub const fn row(&self, i: usize) -> Vec3 {
        Vec3::new(self.0[i][0], self.0[i][1], self.0[i][2])
    }

    /// `M·v`.
    pub fn mul_vec(&self, v: Vec3) -> Vec3 {
        Vec3::new(self.row(0).dot(v), self.row(1).dot(v), self.row(2).dot(v))
    }

    /// `Mᵀ`.
    pub fn transpose(&self) -> Self {
        let m = &self.0;
        Self([
            [m[0][0], m[1][0], m[2][0]],
            [m[0][1], m[1][1], m[2][1]],
            [m[0][2], m[1][2], m[2][2]],
        ])
    }

    /// `self · o`.
    pub fn mul(&self, o: &Self) -> Self {
        let t = o.transpose();
        let mut r = [[0.0; 3]; 3];
        for (i, row) in r.iter_mut().enumerate() {
            for (j, v) in row.iter_mut().enumerate() {
                *v = self.row(i).dot(t.row(j));
            }
        }
        Self(r)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn products() {
        let a = Vec3::new(1.0, 2.0, 3.0);
        let b = Vec3::new(-4.0, 5.0, 0.5);
        assert_eq!(a.dot(b), 7.5);
        assert_eq!(Vec3::X.cross(Vec3::Y), Vec3::Z);
        assert_eq!(Vec3::Y.cross(Vec3::Z), Vec3::X);
        let c = a.cross(b);
        assert_eq!(c.dot(a), 0.0);
        assert_eq!(c.dot(b), 0.0);
        assert_eq!(Vec3::new(3.0, 4.0, 12.0).norm(), 13.0);
        assert_eq!(Vec3::ZERO.normalized(), Vec3::ZERO);
    }

    #[test]
    fn matrices() {
        let m = Mat3([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 10.0]]);
        assert_eq!(m.mul(&Mat3::IDENTITY), m);
        assert_eq!(Mat3::IDENTITY.mul(&m), m);
        assert_eq!(m.transpose().transpose(), m);
        assert_eq!(
            m.mul_vec(Vec3::new(1.0, 0.0, -1.0)),
            Vec3::new(-2.0, -2.0, -3.0)
        );
        let p = m.mul(&m.transpose());
        assert_eq!(p.0[0][1], p.0[1][0]);
        assert_eq!(p.0[0][0], 14.0);
    }
}
