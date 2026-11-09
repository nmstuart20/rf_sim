//! Pre-sampling tap responses: H_rx(f) for the virtual receiver and B_out(f) for an
//! antenna-port output target (contracts §5.1).
//!
//! Every render target declares the response that band-limits its signal before it is
//! sampled, so nothing outside the band can alias into the output. Responses are
//! control-plane objects. They are evaluated in f64 with `libm`, never on the render
//! path. Render kernels only see what a response produces: folded amplitude and phase
//! for narrowband emitters, and edge tables for pulses.
//!
//! # Response family (settled, prompt 5)
//!
//! A linear-phase lowpass at complex baseband (frequency relative to the channel centre)
//! with a Gaussian transition:
//!
//! ```text
//! |H(f)| = ½·[erf((f + f_c)/(√2·σ)) − erf((f − f_c)/(√2·σ))]      f_c = (f_p + f_s)/2
//! H(f)   = |H(f)|·exp(−j·2π·f·τ_g)
//! h(t)   = 2f_c·sinc(2f_c·(t − τ_g))·exp(−2π²σ²·(t − τ_g)²)     sinc x = sin(πx)/(πx)
//! ```
//!
//! |H| is the ideal brick wall of half-width f_c convolved with a Gaussian of standard
//! deviation σ, so h is a sinc times a Gaussian. With z = (f_s − f_c)/(√2σ), σ is
//! chosen so that the passband droop at f_p, `1 − |H(f_p)| = ½·erfc(z) +
//! ½·erfc((f_p + f_c)/(√2σ))`, equals 10^(−A/20). Then:
//!
//! - **Stopband:** |H(f_s)| = ½·erfc(z) − ½·erfc((f_s + f_c)/(√2σ)) ≤ ½·erfc(z) ≤
//!   10^(−A/20), and |H| decreases monotonically beyond f_s. The attenuation is at
//!   least A dB, not a flat floor.
//! - **Passband:** 1 − 10^(−A/20) ≤ |H(f)| ≤ 1 for |f| ≤ f_p, with no ripple. Both
//!   Gaussian tails count toward the droop; the far one matters only for very narrow
//!   passbands (f_p + f_c not much larger than f_s − f_c).
//! - **Transition:** f_p < |f| < f_s, with |H(f_c)| = ½·erf(2f_c/(√2σ)), which is ½
//!   to within the far tail (−6.02 dB for any realistic band).
//! - **Group delay:** τ_g for every frequency (linear phase).
//! - **Impulse truncation:** h is used only on |t − τ_g| ≤ T_h with
//!   T_h = √(c/(2π²σ²)), c = [`TRUNCATION_C`]. Because |2f_c·sinc(2f_c·t)| ≤ 1/(π|t|),
//!   the discarded tail satisfies ∫|h| ≤ e^(−c)/(π·c) < 10⁻⁷ ([`TRUNCATION_TAIL_BOUND`]).
//!
//! A response is accepted only if `0 < f_p < f_s ≤ F_s/2` (anti-alias: everything that
//! could fold into the output is at least A dB down before sampling), `20 ≤ A ≤ 200` dB
//! and `0 ≤ τ_g ≤ 1 ms`. All parameters are integers, so the identity hash is exact.

use std::fmt;

use sha2::{Digest, Sha256};

use crate::SampleRate;

/// Truncation constant c: T_h = √(c/(2π²σ²)). e^(−13)/(13π) ≈ 5.5·10⁻⁸.
pub const TRUNCATION_C: f64 = 13.0;

/// Proven upper bound on ∫|h| over the discarded impulse-response tails.
pub const TRUNCATION_TAIL_BOUND: f64 = 1e-7;

/// Smallest accepted stopband attenuation, milli-dB (20 dB).
pub const MIN_STOPBAND_ATTEN_MDB: u32 = 20_000;
/// Largest accepted stopband attenuation, milli-dB (200 dB).
pub const MAX_STOPBAND_ATTEN_MDB: u32 = 200_000;
/// Largest accepted group delay, ps (1 ms).
pub const MAX_GROUP_DELAY_PS: i64 = 1_000_000_000;

/// Which of a channel's two taps a render target is (contracts §5).
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
#[repr(u32)]
pub enum TapKind {
    /// The antenna-port scene: B_out(f), output-device timebase, no receiver effects.
    AntennaPort = 1,
    /// The virtual receiver: H_rx(f), receiver ClockMap, the full receiver chain.
    VirtualReceiver = 2,
}

impl TapKind {
    /// Stable code used in identities and work lists.
    pub const fn code(self) -> u32 {
        self as u32
    }
}

/// Integer definition of a tap response. Frequencies are one-sided and relative to the
/// channel centre.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct ResponseSpec {
    /// Passband edge f_p, µHz.
    pub passband_edge_uhz: i64,
    /// Stopband edge f_s, µHz.
    pub stopband_edge_uhz: i64,
    /// Stopband attenuation A at and beyond f_s, milli-dB.
    pub stopband_atten_mdb: u32,
    /// Group delay τ_g, ps.
    pub group_delay_ps: i64,
}

/// A rejected response definition.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ResponseError {
    /// The edges do not satisfy `0 < f_p < f_s`.
    Edges,
    /// `f_s` exceeds half the output sample rate, so out-of-band content could alias.
    AliasUnsafe,
    /// The stopband attenuation is outside 20–200 dB.
    Attenuation,
    /// The group delay is negative or above 1 ms.
    GroupDelay,
}

impl fmt::Display for ResponseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::Edges => "tap response needs 0 < passband edge < stopband edge",
            Self::AliasUnsafe => {
                "tap response stopband edge exceeds half the output sample rate \
                 (out-of-band content would alias)"
            }
            Self::Attenuation => "tap response stopband attenuation must be 20 to 200 dB",
            Self::GroupDelay => "tap response group delay must be 0 to 1 ms",
        })
    }
}

impl std::error::Error for ResponseError {}

/// Content identity of a tap response: SHA-256 of the tap kind, the response
/// definition, the output rate and the channel centre (contracts §5).
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct TapResponseId(pub [u8; 32]);

impl fmt::Display for TapResponseId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for b in self.0 {
            write!(f, "{b:02x}")?;
        }
        Ok(())
    }
}

/// A validated tap response for one render target, with its derived constants.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TapResponse {
    kind: TapKind,
    spec: ResponseSpec,
    rate: SampleRate,
    center_uhz: i64,
    /// f_c in Hz.
    cutoff_hz: f64,
    /// σ in Hz.
    sigma_hz: f64,
    /// τ_g in seconds.
    group_delay_s: f64,
    /// T_h in seconds.
    half_width_s: f64,
    id: TapResponseId,
}

const ID_TAG: &[u8] = b"iqsim.tap-response.gauss-lowpass.v1\0";

/// z such that the passband droop ½·erfc(z) + ½·erfc(z·ratio) equals y, by 200 fixed
/// bisection steps on libm erfc. `ratio = (f_p + f_c)/(f_s − f_c)`, so z·ratio is the
/// far Gaussian tail's argument at f_p. The droop decreases in z.
fn droop_inverse(y: f64, ratio: f64) -> f64 {
    let (mut lo, mut hi) = (0.0_f64, 40.0_f64);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        let droop = 0.5 * libm::erfc(mid) + 0.5 * libm::erfc(mid * ratio);
        if droop > y {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

impl TapResponse {
    /// Validates `spec` for a target of `kind` at output `rate` and centre `center_uhz`.
    pub fn new(
        kind: TapKind,
        spec: ResponseSpec,
        rate: SampleRate,
        center_uhz: i64,
    ) -> Result<Self, ResponseError> {
        let fp = spec.passband_edge_uhz;
        let fs = spec.stopband_edge_uhz;
        if fp <= 0 || fs <= fp {
            return Err(ResponseError::Edges);
        }
        // f_s ≤ F_s/2  ⇔  2·f_s[µHz]·den ≤ 10⁶·num, exactly in i128.
        if 2 * i128::from(fs) * i128::from(rate.den()) > 1_000_000 * i128::from(rate.num()) {
            return Err(ResponseError::AliasUnsafe);
        }
        if !(MIN_STOPBAND_ATTEN_MDB..=MAX_STOPBAND_ATTEN_MDB).contains(&spec.stopband_atten_mdb) {
            return Err(ResponseError::Attenuation);
        }
        if !(0..=MAX_GROUP_DELAY_PS).contains(&spec.group_delay_ps) {
            return Err(ResponseError::GroupDelay);
        }
        let fp_hz = fp as f64 * 1e-6;
        let fs_hz = fs as f64 * 1e-6;
        let cutoff_hz = 0.5 * (fp_hz + fs_hz);
        let level = libm::exp10(-(f64::from(spec.stopband_atten_mdb) * 1e-3) / 20.0);
        let z = droop_inverse(level, (fp_hz + cutoff_hz) / (fs_hz - cutoff_hz));
        let sigma_hz = (fs_hz - cutoff_hz) / (core::f64::consts::SQRT_2 * z);
        let a = 2.0 * core::f64::consts::PI * core::f64::consts::PI * sigma_hz * sigma_hz;
        let half_width_s = (TRUNCATION_C / a).sqrt();

        let mut h = Sha256::new();
        h.update(ID_TAG);
        h.update(kind.code().to_le_bytes());
        h.update(fp.to_le_bytes());
        h.update(fs.to_le_bytes());
        h.update(spec.stopband_atten_mdb.to_le_bytes());
        h.update(spec.group_delay_ps.to_le_bytes());
        h.update(rate.num().to_le_bytes());
        h.update(rate.den().to_le_bytes());
        h.update(center_uhz.to_le_bytes());
        Ok(Self {
            kind,
            spec,
            rate,
            center_uhz,
            cutoff_hz,
            sigma_hz,
            group_delay_s: spec.group_delay_ps as f64 * 1e-12,
            half_width_s,
            id: TapResponseId(h.finalize().into()),
        })
    }

    /// The tap this response belongs to.
    pub const fn kind(&self) -> TapKind {
        self.kind
    }

    /// The integer definition.
    pub const fn spec(&self) -> ResponseSpec {
        self.spec
    }

    /// Output sample rate of the target.
    pub const fn rate(&self) -> SampleRate {
        self.rate
    }

    /// Channel centre frequency, µHz.
    pub const fn center_uhz(&self) -> i64 {
        self.center_uhz
    }

    /// Content identity.
    pub const fn id(&self) -> TapResponseId {
        self.id
    }

    /// Half-amplitude frequency f_c, Hz.
    pub const fn cutoff_hz(&self) -> f64 {
        self.cutoff_hz
    }

    /// Gaussian transition standard deviation σ, Hz.
    pub const fn sigma_hz(&self) -> f64 {
        self.sigma_hz
    }

    /// Group delay τ_g, seconds.
    pub const fn group_delay_s(&self) -> f64 {
        self.group_delay_s
    }

    /// Impulse-response truncation half-width T_h, seconds.
    pub const fn half_width_s(&self) -> f64 {
        self.half_width_s
    }

    /// Lag interval `[τ_g − T_h, τ_g + T_h]` (seconds) outside which h is taken as 0.
    pub fn support_s(&self) -> (f64, f64) {
        (
            self.group_delay_s - self.half_width_s,
            self.group_delay_s + self.half_width_s,
        )
    }

    /// |H(f)| at offset `f_hz` from the channel centre.
    ///
    /// Inside f_c it sums two positive erf terms. Outside f_c it subtracts two erfc
    /// terms whose second is much smaller, so the far stopband keeps full relative
    /// precision instead of cancelling to zero.
    pub fn magnitude(&self, f_hz: f64) -> f64 {
        let f = f_hz.abs();
        let s = core::f64::consts::SQRT_2 * self.sigma_hz;
        if f <= self.cutoff_hz {
            0.5 * (libm::erf((f + self.cutoff_hz) / s) + libm::erf((self.cutoff_hz - f) / s))
        } else {
            0.5 * (libm::erfc((f - self.cutoff_hz) / s) - libm::erfc((f + self.cutoff_hz) / s))
        }
    }

    /// Phase of H(f) in cycles: −f·τ_g.
    pub fn phase_cycles(&self, f_hz: f64) -> f64 {
        -(f_hz * self.group_delay_s)
    }

    /// The complex response H(f) as (re, im).
    pub fn response(&self, f_hz: f64) -> (f64, f64) {
        let m = self.magnitude(f_hz);
        let angle = 2.0 * core::f64::consts::PI * self.phase_cycles(f_hz);
        let (s, c) = libm::sincos(angle);
        (m * c, m * s)
    }

    /// The truncated impulse response h(t), 1/s. Zero outside [`Self::support_s`].
    pub fn impulse(&self, t_s: f64) -> f64 {
        let x = t_s - self.group_delay_s;
        if x.abs() > self.half_width_s {
            return 0.0;
        }
        let two_fc = 2.0 * self.cutoff_hz;
        let arg = core::f64::consts::PI * two_fc * x;
        let sinc = if arg == 0.0 {
            1.0
        } else {
            libm::sin(arg) / arg
        };
        let a = 2.0 * core::f64::consts::PI * core::f64::consts::PI * self.sigma_hz * self.sigma_hz;
        two_fc * sinc * libm::exp(-(a * x * x))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spec(fp_mhz: f64, fs_mhz: f64, atten_db: u32, gd_ps: i64) -> ResponseSpec {
        ResponseSpec {
            passband_edge_uhz: (fp_mhz * 1e12) as i64,
            stopband_edge_uhz: (fs_mhz * 1e12) as i64,
            stopband_atten_mdb: atten_db * 1000,
            group_delay_ps: gd_ps,
        }
    }

    fn rx() -> TapResponse {
        let rate = SampleRate::from_hz(16_000_000).unwrap();
        TapResponse::new(
            TapKind::VirtualReceiver,
            spec(6.0, 7.5, 80, 1_000_000),
            rate,
            1_090_000_000_000_000,
        )
        .unwrap()
    }

    #[test]
    fn truncation_constant_meets_its_bound() {
        let c = TRUNCATION_C;
        assert!(libm::exp(-c) / (core::f64::consts::PI * c) < TRUNCATION_TAIL_BOUND);
    }

    #[test]
    fn band_edges_meet_the_declared_attenuation() {
        let r = rx();
        let level = 1e-4;
        let at_fs = r.magnitude(7.5e6);
        assert!(at_fs <= level * (1.0 + 1e-9), "{at_fs}");
        assert!(at_fs >= level * (1.0 - 1e-3), "{at_fs}");
        assert!((1.0 - r.magnitude(6.0e6)) <= level * (1.0 + 1e-9));
        assert!((r.magnitude(6.75e6) - 0.5).abs() < 1e-12);
        assert_eq!(r.magnitude(-3e6), r.magnitude(3e6));
        // Monotone beyond f_s and never above 1.
        let mut prev = at_fs;
        for i in 1..200 {
            let m = r.magnitude(7.5e6 + f64::from(i) * 50e3);
            assert!(m <= prev, "not monotone at step {i}");
            prev = m;
        }
        for i in 0..=100 {
            assert!(r.magnitude(f64::from(i) * 60e3) <= 1.0);
        }
        // The far stopband keeps relative precision (no cancellation to zero) ...
        assert!(r.magnitude(8.0e6) > 0.0 && r.magnitude(8.0e6) < 1e-6);
        // ... and the phase is linear with the group delay.
        assert_eq!(r.phase_cycles(2e6), -2e6 * 1e-6);
    }

    #[test]
    fn narrow_passband_still_meets_the_droop_and_stopband() {
        // Regression (prompt 5 review): fp = 1 Hz, fs = 1 kHz, A = 20 dB. Solving σ
        // from the near tail alone gave |H(fp)| = 0.8009, below the promised 0.9.
        let rate = SampleRate::from_hz(16_000).unwrap();
        let spec = ResponseSpec {
            passband_edge_uhz: 1_000_000,
            stopband_edge_uhz: 1_000_000_000,
            stopband_atten_mdb: 20_000,
            group_delay_ps: 0,
        };
        let r = TapResponse::new(TapKind::AntennaPort, spec, rate, 0).unwrap();
        assert!(
            r.magnitude(1.0) >= 0.9 * (1.0 - 1e-12),
            "{}",
            r.magnitude(1.0)
        );
        assert!(
            r.magnitude(1000.0) <= 0.1 * (1.0 + 1e-12),
            "{}",
            r.magnitude(1000.0)
        );
        for f in [0.0, 0.5, 1.0] {
            assert!(r.magnitude(f) >= 0.9 * (1.0 - 1e-12));
        }
    }

    #[test]
    fn impulse_response_transforms_to_the_magnitude() {
        // Numerical Fourier transform of the truncated h (Simpson) against H(f).
        let r = rx();
        let (lo, hi) = r.support_s();
        let n = 40_000;
        let dt = (hi - lo) / f64::from(n);
        for f in [0.0, 2.5e6, 6.0e6, 6.75e6, 7.5e6] {
            let (mut re, mut im) = (0.0, 0.0);
            for i in 0..=n {
                let t = lo + f64::from(i) * dt;
                let w = if i == 0 || i == n {
                    1.0
                } else if i % 2 == 1 {
                    4.0
                } else {
                    2.0
                };
                let (s, c) = libm::sincos(-2.0 * core::f64::consts::PI * f * t);
                re += w * r.impulse(t) * c;
                im += w * r.impulse(t) * s;
            }
            re *= dt / 3.0;
            im *= dt / 3.0;
            let (hr, hi_) = r.response(f);
            let err = ((re - hr) * (re - hr) + (im - hi_) * (im - hi_)).sqrt();
            assert!(err < 2e-7, "f = {f}: |FT(h) − H| = {err}");
        }
    }

    #[test]
    fn validation_rejects_bad_definitions() {
        let rate = SampleRate::from_hz(16_000_000).unwrap();
        let k = TapKind::AntennaPort;
        let bad = |s| TapResponse::new(k, s, rate, 0).unwrap_err();
        assert_eq!(bad(spec(6.0, 6.0, 80, 0)), ResponseError::Edges);
        assert_eq!(bad(spec(0.0, 6.0, 80, 0)), ResponseError::Edges);
        assert_eq!(bad(spec(6.0, 8.000001, 80, 0)), ResponseError::AliasUnsafe);
        assert!(TapResponse::new(k, spec(6.0, 8.0, 80, 0), rate, 0).is_ok());
        assert_eq!(bad(spec(6.0, 7.0, 19, 0)), ResponseError::Attenuation);
        assert_eq!(bad(spec(6.0, 7.0, 201, 0)), ResponseError::Attenuation);
        assert_eq!(bad(spec(6.0, 7.0, 80, -1)), ResponseError::GroupDelay);
        assert_eq!(
            bad(spec(6.0, 7.0, 80, MAX_GROUP_DELAY_PS + 1)),
            ResponseError::GroupDelay
        );
    }

    #[test]
    fn identity_covers_kind_definition_rate_and_centre() {
        let rate = SampleRate::from_hz(16_000_000).unwrap();
        let s = spec(6.0, 7.5, 80, 0);
        let base = TapResponse::new(TapKind::VirtualReceiver, s, rate, 5)
            .unwrap()
            .id();
        let other_kind = TapResponse::new(TapKind::AntennaPort, s, rate, 5)
            .unwrap()
            .id();
        let other_centre = TapResponse::new(TapKind::VirtualReceiver, s, rate, 6)
            .unwrap()
            .id();
        let rate2 = SampleRate::from_hz(20_000_000).unwrap();
        let other_rate = TapResponse::new(TapKind::VirtualReceiver, s, rate2, 5)
            .unwrap()
            .id();
        let mut s2 = s;
        s2.group_delay_ps = 1;
        let other_spec = TapResponse::new(TapKind::VirtualReceiver, s2, rate, 5)
            .unwrap()
            .id();
        let ids = [base, other_kind, other_centre, other_rate, other_spec];
        for i in 0..ids.len() {
            for j in i + 1..ids.len() {
                assert_ne!(ids[i], ids[j]);
            }
        }
        assert_eq!(
            base,
            TapResponse::new(TapKind::VirtualReceiver, s, rate, 5)
                .unwrap()
                .id()
        );
    }
}
