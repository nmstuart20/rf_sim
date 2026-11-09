//! Golden hashes over fixed sets of numeric-core evaluations (prompt 2).
//!
//! Each golden covers 10⁶ evaluations from a fixed SplitMix64 input stream, serialized
//! little-endian. A change to any output bit fails the test; bless only on purpose
//! (`cargo xtask bless`) and explain it in docs/progress.md. The CUDA backend
//! (prompt 20) must reproduce these hashes.

use iqsim_core::wide::U128;
use iqsim_core::{AcqInstant, DelayPoly, OffsetFrequency, Phase, PhasePoly, SIN_TABLE, sincos};
use iqsim_testkit::assert_golden;

const N: usize = 1_000_000;

/// SplitMix64: a fixed test-input generator (not simulation randomness).
struct SplitMix(u64);

impl SplitMix {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    fn u128(&mut self) -> U128 {
        U128::new(self.next(), self.next())
    }
}

#[test]
fn sincos_table_golden() {
    let bytes: Vec<u8> = SIN_TABLE
        .iter()
        .flat_map(|v| v.to_bits().to_le_bytes())
        .collect();
    assert_golden("core/sincos_table", &bytes);
}

#[test]
fn sincos_golden() {
    let mut rng = SplitMix(0x0002_51C0);
    let mut bytes = Vec::with_capacity(N * 8);
    for _ in 0..N {
        let (s, c) = sincos(Phase::from_bits(rng.next()));
        bytes.extend(s.to_bits().to_le_bytes());
        bytes.extend(c.to_bits().to_le_bytes());
    }
    assert_golden("core/sincos_1e6", &bytes);
}

/// 100 polynomials (delay and phase alternating) × 10 000 samples. Even polynomials
/// are evaluated by jumping to each k, odd ones by stepping a cursor.
#[test]
fn block_poly_golden() {
    let mut rng = SplitMix(0x0002_B10C);
    let len = 10_000u64;
    let mut bytes = Vec::with_capacity(N * 8);
    for i in 0..100 {
        let regs = [rng.u128(), rng.u128(), rng.u128(), rng.u128()];
        if i % 4 < 2 {
            let p = DelayPoly::from_registers(regs, len).unwrap();
            let mut cur = p.cursor(0).unwrap();
            for k in 0..len {
                let v = if i % 2 == 0 {
                    p.eval(k).unwrap()
                } else {
                    cur.value()
                };
                bytes.extend(v.to_bits().to_le_bytes());
                cur.step();
            }
        } else {
            let p = PhasePoly::from_registers(regs, len).unwrap();
            let mut cur = p.cursor(0).unwrap();
            for k in 0..len {
                let v = if i % 2 == 0 {
                    p.eval(k).unwrap()
                } else {
                    cur.value()
                };
                bytes.extend(v.to_bits().to_le_bytes());
                cur.step();
            }
        }
    }
    assert_eq!(bytes.len(), N * 8);
    assert_golden("core/block_poly_1e6", &bytes);
}

#[test]
fn offset_phase_golden() {
    let mut rng = SplitMix(0x0002_0FF5);
    let mut bytes = Vec::with_capacity(N * 8);
    for _ in 0..N {
        let f = OffsetFrequency::from_bits(rng.next() as i64);
        let t = AcqInstant::from_q64_64(rng.u128().to_i128());
        bytes.extend(f.phase_at(t).to_bits().to_le_bytes());
    }
    assert_golden("core/offset_phase_1e6", &bytes);
}
