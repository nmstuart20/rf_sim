//! Statistical tests of the Gaussian noise generator (prompt 3).
//!
//! Every test uses fixed seeds, so its outcome is deterministic. The thresholds are
//! fixed before looking at the data: 6 standard errors for moments and correlations,
//! and the 1 % asymptotic Anderson-Darling critical value for a fully specified
//! distribution (3.857).
//!
//! The reference CDFs use `libm` (erfc, expm1, log) in f64. libm is a dev-dependency
//! here, never a simulation input.

use iqsim_core::{NoiseStream, SampleIndex, StreamId};

/// Samples for the moment and distribution tests.
const N_DIST: u64 = 10_000_000;
/// Samples per pair for the cross-correlation tests.
const N_CORR: u64 = 1_000_000;
/// Anderson-Darling 1 % critical value (asymptotic, fully specified distribution).
const AD_CRIT_1PCT: f64 = 3.857;
/// Bound, in standard errors, for moments and correlations.
const K_SIGMA: f64 = 6.0;

const SEED: u64 = 0x5EED_0003;
const ENTITY: u32 = 17;

fn samples(stream: NoiseStream, n: u64) -> Vec<(f32, f32)> {
    (0..n).map(|i| stream.sample(SampleIndex::new(i))).collect()
}

/// Neumaier-compensated sum.
#[derive(Default)]
struct Sum {
    s: f64,
    c: f64,
}

impl Sum {
    fn add(&mut self, x: f64) {
        let t = self.s + x;
        if self.s.abs() >= x.abs() {
            self.c += (self.s - t) + x;
        } else {
            self.c += (x - t) + self.s;
        }
        self.s = t;
    }

    fn get(&self) -> f64 {
        self.s + self.c
    }
}

/// Anderson-Darling A² of `xs` against a CDF given as (ln F(x), ln(1 − F(x))).
///
/// A² = −n − (1/n) Σᵢ (2i − 1)·[ln F(x₍ᵢ₎) + ln(1 − F(x₍ₙ₊₁₋ᵢ₎))], reindexed so each
/// sorted point contributes once, and summed as Σ −(1 + …) so no large cancellation
/// remains.
fn anderson_darling(mut xs: Vec<f64>, ln_cdf: impl Fn(f64) -> (f64, f64)) -> f64 {
    xs.sort_unstable_by(f64::total_cmp);
    let n = xs.len() as f64;
    let mut sum = Sum::default();
    for (i0, &x) in xs.iter().enumerate() {
        let i = i0 as f64 + 1.0;
        let (ln_f, ln_sf) = ln_cdf(x);
        let w_lo = (2.0 * i - 1.0) / n;
        let w_hi = (2.0 * (n - i) + 1.0) / n;
        sum.add(-(1.0 + w_lo * ln_f + w_hi * ln_sf));
    }
    sum.get()
}

/// Standard normal: ln Φ(x) and ln Φ(−x), accurate in both tails.
fn ln_normal_cdf(x: f64) -> (f64, f64) {
    let r = std::f64::consts::FRAC_1_SQRT_2;
    (
        libm::log(0.5 * libm::erfc(-x * r)),
        libm::log(0.5 * libm::erfc(x * r)),
    )
}

/// Exp(1): ln(1 − e⁻ʸ) and −y.
fn ln_exp_cdf(y: f64) -> (f64, f64) {
    (libm::log(-libm::expm1(-y)), -y)
}

/// Moments of 10⁷ samples: means, component variances (½ each), unit total power,
/// I/Q correlation, skewness and excess kurtosis; then Anderson-Darling on I and Q
/// (scaled to unit variance) against N(0, 1) and on |z|² against Exp(1).
#[test]
fn moments_and_anderson_darling_1e7() {
    let s = NoiseStream::new(SEED, ENTITY, StreamId::receiver_thermal_noise(0));
    let z = samples(s, N_DIST);
    let n = N_DIST as f64;

    let mut acc = [(); 6].map(|_| Sum::default());
    let mut acc_q = [(); 4].map(|_| Sum::default());
    for &(i, q) in &z {
        let (i, q) = (i as f64, q as f64);
        acc[0].add(i);
        acc[1].add(i * i);
        acc[2].add(i * i * i);
        acc[3].add(i * i * i * i);
        acc[4].add(i * q);
        acc[5].add(i * i + q * q);
        acc_q[0].add(q);
        acc_q[1].add(q * q);
        acc_q[2].add(q * q * q);
        acc_q[3].add(q * q * q * q);
    }
    let m = |a: &Sum| a.get() / n;
    let power = m(&acc[5]);
    let iq = m(&acc[4]);
    eprintln!("E|z|² = {power:.6}, E[IQ] = {iq:+.2e}");
    // Std of the mean of |z|² (Exp(1)) is 1/√n; of the mean of I·Q it is ½/√n.
    assert!((power - 1.0).abs() < K_SIGMA / n.sqrt(), "E|z|² = {power}");
    assert!(iq.abs() < K_SIGMA * 0.5 / n.sqrt(), "E[IQ] = {iq}");

    for (name, a) in [("I", &acc[..4]), ("Q", &acc_q[..])] {
        let mean = m(&a[0]);
        let var = m(&a[1]) - mean * mean;
        // Central moments from raw moments (mean ≈ 0, so this is well conditioned).
        let m3 = m(&a[2]) - 3.0 * mean * m(&a[1]) + 2.0 * mean * mean * mean;
        let m4 = m(&a[3]) - 4.0 * mean * m(&a[2]) + 6.0 * mean * mean * m(&a[1])
            - 3.0 * mean * mean * mean * mean;
        let skew = m3 / (var * var.sqrt());
        let kurt = m4 / (var * var) - 3.0;
        eprintln!("{name}: mean {mean:+.2e}, var {var:.6}, skew {skew:+.2e}, ex. kurt {kurt:+.2e}");
        // Standard errors for a normal with variance ½: mean √(½/n), variance
        // √(2·¼/n), skewness √(6/n), excess kurtosis √(24/n).
        assert!(
            mean.abs() < K_SIGMA * (0.5 / n).sqrt(),
            "{name} mean {mean}"
        );
        assert!(
            (var - 0.5).abs() < K_SIGMA * (0.5 / n).sqrt(),
            "{name} var {var}"
        );
        assert!(
            skew.abs() < K_SIGMA * (6.0 / n).sqrt(),
            "{name} skew {skew}"
        );
        assert!(
            kurt.abs() < K_SIGMA * (24.0 / n).sqrt(),
            "{name} kurt {kurt}"
        );
    }

    let sqrt2 = std::f64::consts::SQRT_2;
    let a_i = anderson_darling(
        z.iter().map(|&(i, _)| i as f64 * sqrt2).collect(),
        ln_normal_cdf,
    );
    let a_q = anderson_darling(
        z.iter().map(|&(_, q)| q as f64 * sqrt2).collect(),
        ln_normal_cdf,
    );
    let a_p = anderson_darling(
        z.iter()
            .map(|&(i, q)| {
                let (i, q) = (i as f64, q as f64);
                i * i + q * q
            })
            .collect(),
        ln_exp_cdf,
    );
    eprintln!("Anderson-Darling A²: I {a_i:.3}, Q {a_q:.3}, |z|² {a_p:.3}");
    assert!(a_i < AD_CRIT_1PCT, "A²(I) = {a_i}");
    assert!(a_q < AD_CRIT_1PCT, "A²(Q) = {a_q}");
    assert!(a_p < AD_CRIT_1PCT, "A²(|z|²) = {a_p}");
}

/// Complex cross-correlation at lags −2..=2 and power correlation at lag 0.
fn assert_uncorrelated(label: &str, a: &[(f32, f32)], b: &[(f32, f32)]) {
    let n = a.len().min(b.len()) - 2;
    let nf = n as f64;
    // The lagged sums run over k in 2..n: n − 2 products.
    let nl = (n - 2) as f64;
    let mut worst = 0.0f64;
    for lag in -2i64..=2 {
        let (mut re, mut im) = (Sum::default(), Sum::default());
        for k in 2..n {
            let (ai, aq) = a[k];
            let (bi, bq) = b[(k as i64 + lag) as usize];
            let (ai, aq, bi, bq) = (ai as f64, aq as f64, bi as f64, bq as f64);
            // a · conj(b)
            re.add(ai * bi + aq * bq);
            im.add(aq * bi - ai * bq);
        }
        let (re, im) = (re.get() / nl, im.get() / nl);
        let mag = (re * re + im * im).sqrt();
        // Re and Im of the product each have variance ½, so |ρ| is Rayleigh with
        // σ = √(½/n).
        assert!(
            mag < K_SIGMA * (0.5 / nl).sqrt(),
            "{label}: |ρ({lag})| = {mag}"
        );
        worst = worst.max(mag);
    }
    // Pearson correlation of |a|² and |b|² (each Exp(1): mean 1, variance 1).
    let mut c = Sum::default();
    for k in 0..n {
        let (ai, aq, bi, bq) = (a[k].0 as f64, a[k].1 as f64, b[k].0 as f64, b[k].1 as f64);
        let pa = ai * ai + aq * aq;
        let pb = bi * bi + bq * bq;
        c.add((pa - 1.0) * (pb - 1.0));
    }
    let rho_p = c.get() / nf;
    assert!(
        rho_p.abs() < K_SIGMA / nf.sqrt(),
        "{label}: power ρ = {rho_p}"
    );
    eprintln!("{label}: max |ρ| over lags {worst:.2e}, power ρ {rho_p:+.2e}");
}

/// Neighbouring streams share no structure: adjacent channels (stream ids differ by
/// one), adjacent entities, adjacent seeds, and nearby samples of the same stream.
#[test]
fn adjacent_streams_are_uncorrelated() {
    let ch0 = StreamId::receiver_thermal_noise(0);
    let ch1 = StreamId::receiver_thermal_noise(1);
    let base = samples(NoiseStream::new(SEED, ENTITY, ch0), N_CORR);
    let pairs = [
        ("channel+1", NoiseStream::new(SEED, ENTITY, ch1)),
        ("entity+1", NoiseStream::new(SEED, ENTITY + 1, ch0)),
        ("seed+1", NoiseStream::new(SEED + 1, ENTITY, ch0)),
    ];
    for (label, s) in pairs {
        assert_uncorrelated(label, &base, &samples(s, N_CORR));
    }
    // The same stream against itself shifted by 3 samples: the lag window −2..=2
    // covers stream lags 1..=5 (power correlation at lag 3).
    assert_uncorrelated("same stream, lags 1..=5", &base, &base[3..]);
}
