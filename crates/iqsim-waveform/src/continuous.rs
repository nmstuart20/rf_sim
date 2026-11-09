//! Continuous sources and the stage-1 resampling filter (contracts §5.6).
//!
//! A continuous source is a sequence x\[n\] of complex samples at a native rate F_n
//! (an IQ file now, generated waveforms later), normalized so that |x| = 1 is full
//! scale. Sample n is emitted at `start + n/F_n` on the emitter's carrier.
//!
//! # Signal model
//!
//! The continuous emitted envelope is `x(t) = Σ_n x[n]·φ(t − n)` (t in native
//! samples), where φ is the **source reconstruction response**, a zero-group-delay
//! member of the tap-response family (§5.1, [`TapKind::SourceReconstruction`]) with
//! stopband edge ≤ F_n/2, so every image of the native spectrum is at least its
//! attenuation down.
//!
//! For a render target with response H centred Δf = f_e − f_c below the emitter
//! carrier, the band-limited envelope a path delivers is `z(w) = Σ_n x[n]·K(w − n)`
//! with the **combined kernel**
//!
//! ```text
//! K(t) = ∫ Φ(ν·F_n)·H(ν·F_n + Δf)·e^{j2πνt} dν        (ν in cycles per native sample)
//! ```
//!
//! i.e. φ convolved with h shifted to the emitter's offset, as for pulse edge tables.
//! The renderer multiplies z by the path's carrier rotation, which applies the
//! frequency translation exactly at the target rate, so a low native rate never has to
//! translate its own samples (no intermediate aliasing). K is computed from this
//! frequency integral over the finite band where both responses are non-negligible
//! (8-point Gauss–Legendre, a convergence check by doubling), then truncated to the
//! sum of both impulse supports, which costs at most ~10⁻⁷·(‖φ‖₁ + ‖h‖₁).
//!
//! # Stage 1: filtering at the native rate
//!
//! [`ResampleKernel::filter_window`] evaluates z on a grid of U = 2^log2_up points per
//! native sample, `z[m] = Σ_i x[⌊m/U⌋ − i]·K((m mod U)/U + i)`, in f64 with a fixed
//! order (increasing i), and rounds each point once to f32. Each point is a pure
//! function of m and the source, so windows of different block splits agree bit for bit
//! where they overlap; the filter "history" is recomputed from the immutable source
//! instead of being carried.
//!
//! U is the smallest power of two (≤ 32) for which the 6-tap Lagrange Farrow kernel
//! that reads the window (stage 2, in the renderer) errs by at most
//! [`FARROW_ERROR_BOUND`] at the highest frequency z can contain,
//! `ν_max = min(f_s,φ, |Δf| + f_s,H)/(U·F_n)` ([`lagrange6_error`]).
//!
//! # Boundaries ([`EndPolicy`])
//!
//! Before sample 0 the source is zero for every policy, so the onset is band-limited by
//! K. After the last sample it is zero (`Stop`, `ZeroPad`) or periodic,
//! `x[n mod N]` (`Loop`): sample N − 1 is followed directly by sample 0, with no gap or
//! crossfade, and K band-limits the seam. `Stop` is strict: the caller rejects a block
//! whose retarded coordinate passes the last sample ([`SourceTimeline::exhausted`]); the
//! filter halo past the end still reads zeros.

use std::fmt;

use iqsim_core::{ResponseError, ResponseSpec, SampleRate, TapKind, TapResponse, TapResponseId};

use crate::edge::{GL_NODES, GL_WEIGHTS};

/// Largest stage-1 oversampling exponent (U = 32), as in the work list.
pub const MAX_LOG2_UP: u32 = 5;

/// Bound on the stage-2 (Farrow) interpolation error relative to a spectral component's
/// amplitude.
pub const FARROW_ERROR_BOUND: f64 = 2e-6;

/// Bound on the change of any kernel tap between quadrature resolutions, relative to
/// the largest tap.
pub const KERNEL_QUADRATURE_BOUND: f64 = 1e-9;

/// Largest number of native taps per polyphase branch.
pub const MAX_KERNEL_TAPS: i64 = 1 << 16;

/// What happens after a source's last sample.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum EndPolicy {
    /// The source must cover the run: rendering past the last sample is an error.
    Stop,
    /// Silence after the last sample.
    ZeroPad,
    /// The samples repeat with period N.
    Loop,
}

/// A failure while preparing or filtering a continuous source.
#[derive(Clone, Debug, PartialEq)]
pub enum ContinuousError {
    /// A response is invalid.
    Response(ResponseError),
    /// The kernel quadrature did not converge.
    Quadrature,
    /// The kernel is longer than [`MAX_KERNEL_TAPS`].
    KernelLength,
    /// The source is empty or a window request is out of range.
    Range,
    /// Reading the source failed.
    Source(String),
}

impl fmt::Display for ContinuousError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Response(e) => e.fmt(f),
            Self::Quadrature => f.write_str("resampling kernel quadrature did not converge"),
            Self::KernelLength => f.write_str("resampling kernel is too long"),
            Self::Range => f.write_str("source window request is out of range"),
            Self::Source(e) => write!(f, "source read failed: {e}"),
        }
    }
}

impl std::error::Error for ContinuousError {}

impl From<ResponseError> for ContinuousError {
    fn from(e: ResponseError) -> Self {
        Self::Response(e)
    }
}

/// A sampled source: N normalized complex samples, readable by index.
pub trait ContinuousSource {
    /// Number of samples N (> 0).
    fn len(&self) -> u64;
    /// Whether the source is empty (never, for a valid source).
    fn is_empty(&self) -> bool {
        self.len() == 0
    }
    /// Fills `out` with samples `[first, first + out.len())`, all inside `[0, N)`.
    fn fill(&mut self, first: u64, out: &mut [[f64; 2]]) -> Result<(), String>;
}

/// An in-memory source (tests and generated waveforms).
#[derive(Clone, Debug, PartialEq)]
pub struct MemorySource(pub Vec<[f64; 2]>);

impl ContinuousSource for MemorySource {
    fn len(&self) -> u64 {
        self.0.len() as u64
    }
    fn fill(&mut self, first: u64, out: &mut [[f64; 2]]) -> Result<(), String> {
        let a = usize::try_from(first).map_err(|_| "range".to_owned())?;
        let src = self.0.get(a..a + out.len()).ok_or("range")?;
        out.copy_from_slice(src);
        Ok(())
    }
}

/// A source's length and end policy: the extension of x\[n\] to every integer n.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SourceTimeline {
    /// N.
    pub len: u64,
    /// After the last sample.
    pub policy: EndPolicy,
}

impl SourceTimeline {
    /// Whether a retarded coordinate up to `w_max` native samples (Q32.32 bits, relative
    /// to sample 0) passes the last sample under [`EndPolicy::Stop`].
    pub fn exhausted(&self, w_max_bits: i128) -> bool {
        self.policy == EndPolicy::Stop && w_max_bits > (i128::from(self.len) - 1) << 32
    }

    /// Whether every sample in `[first, last]` is a padding zero.
    pub fn all_zero(&self, first: i64, last: i64) -> bool {
        last < 0 || (self.policy != EndPolicy::Loop && i128::from(first) >= i128::from(self.len))
    }

    /// Fills `out` with the extended samples `[first, first + out.len())`.
    pub fn read(
        &self,
        src: &mut dyn ContinuousSource,
        first: i64,
        out: &mut [[f64; 2]],
    ) -> Result<(), ContinuousError> {
        if self.len == 0 || src.len() != self.len {
            return Err(ContinuousError::Range);
        }
        let n = self.len as i128;
        let mut i = 0usize;
        while i < out.len() {
            let idx = i128::from(first) + i as i128;
            let rest = out.len() - i;
            // Length of the run that maps contiguously, and where it starts in x.
            let (run, start) = if idx < 0 {
                ((-idx).min(rest as i128) as usize, None)
            } else if idx < n || self.policy == EndPolicy::Loop {
                let pos = idx.rem_euclid(n);
                (((n - pos).min(rest as i128)) as usize, Some(pos as u64))
            } else {
                (rest, None)
            };
            let dst = &mut out[i..i + run];
            match start {
                Some(s) => src.fill(s, dst).map_err(ContinuousError::Source)?,
                None => dst.fill([0.0, 0.0]),
            }
            i += run;
        }
        Ok(())
    }
}

/// The defining responses of a kernel: everything K depends on.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ResampleKey {
    /// The source reconstruction response (covers its rate).
    pub source: TapResponseId,
    /// The target's tap response (covers tap kind, rate and centre).
    pub target: TapResponseId,
    /// Emitter carrier minus target centre, µHz.
    pub offset_uhz: i64,
}

/// The source reconstruction response for native rate `rate`: passband edge
/// `bandwidth/2`, stopband edge `stopband/2` (µHz, full widths), attenuation in milli-dB.
pub fn source_response(
    rate: SampleRate,
    spec: ResponseSpec,
) -> Result<TapResponse, ContinuousError> {
    Ok(TapResponse::new(
        TapKind::SourceReconstruction,
        spec,
        rate,
        0,
    )?)
}

/// Worst error over μ ∈ [0, 1] (513 points) of the 6-point Lagrange interpolator with
/// nodes −2..3 on a unit complex tone of `nu` cycles per grid sample.
pub fn lagrange6_error(nu: f64) -> f64 {
    let mut worst = 0.0_f64;
    for k in 0..=512 {
        let mu = f64::from(k) / 512.0;
        let (mut re, mut im) = (0.0, 0.0);
        for i in -2..=3_i32 {
            let mut l = 1.0;
            for j in -2..=3_i32 {
                if j != i {
                    l *= (mu - f64::from(j)) / f64::from(i - j);
                }
            }
            let (s, c) = libm::sincos(2.0 * core::f64::consts::PI * nu * (f64::from(i) - mu));
            re += l * c;
            im += l * s;
        }
        let dr = re - 1.0;
        worst = worst.max((dr * dr + im * im).sqrt());
    }
    worst
}

/// The combined kernel K for one (source, target) pair, sampled at its U polyphase
/// branches, with the constants the builder needs.
#[derive(Clone, Debug, PartialEq)]
pub struct ResampleKernel {
    key: ResampleKey,
    log2_up: u32,
    /// Smallest native tap index i of any branch.
    i_lo: i64,
    /// Taps per branch.
    taps_per_branch: usize,
    /// `taps[j·taps_per_branch + (i − i_lo)] = K(j/U + i)`.
    taps: Vec<[f64; 2]>,
    l1_max: f64,
    farrow_error: f64,
    nu_max: f64,
}

/// One point of the band integrand: (ν, weight·Φ·H re, im).
struct Node {
    nu: f64,
    re: f64,
    im: f64,
}

/// Gauss–Legendre nodes over [lo, hi] (cycles per native sample) with `panels` panels,
/// each carrying weight·Φ(ν F_n)·H(ν F_n + Δf).
fn band_nodes(
    phi: &TapResponse,
    h: &TapResponse,
    rate_hz: f64,
    offset_hz: f64,
    lo: f64,
    hi: f64,
    panels: u32,
) -> Vec<Node> {
    let width = (hi - lo) / f64::from(panels);
    let half = 0.5 * width;
    let mut out = Vec::with_capacity(panels as usize * 8);
    for p in 0..panels {
        let mid = lo + (f64::from(p) + 0.5) * width;
        for i in 0..4 {
            for sign in [-1.0, 1.0] {
                let nu = mid + sign * GL_NODES[i] * half;
                let f = nu * rate_hz;
                let a = phi.magnitude(f);
                let (hr, hi_) = h.response(f + offset_hz);
                let w = GL_WEIGHTS[i] * half;
                out.push(Node {
                    nu,
                    re: w * a * hr,
                    im: w * a * hi_,
                });
            }
        }
    }
    out
}

/// K(t) from band nodes, summed in node order.
fn kernel_at(nodes: &[Node], t: f64) -> [f64; 2] {
    let (mut re, mut im) = (0.0, 0.0);
    for n in nodes {
        let (s, c) = libm::sincos(2.0 * core::f64::consts::PI * n.nu * t);
        re += n.re * c - n.im * s;
        im += n.re * s + n.im * c;
    }
    [re, im]
}

/// Frequency (Hz) beyond which a response is below erfc(8)/2 ≈ 10⁻²⁹.
fn band_edge(r: &TapResponse) -> f64 {
    r.cutoff_hz() + 8.0 * core::f64::consts::SQRT_2 * r.sigma_hz()
}

impl ResampleKernel {
    /// Builds K for source response `phi` (native rate) into target response `h` at
    /// carrier offset `offset_uhz` (emitter carrier − target centre).
    pub fn build(
        phi: &TapResponse,
        h: &TapResponse,
        offset_uhz: i64,
    ) -> Result<Self, ContinuousError> {
        let rate = phi.rate();
        let rate_hz = rate.num() as f64 / rate.den() as f64;
        let offset_hz = offset_uhz as f64 * 1e-6;
        let key = ResampleKey {
            source: phi.id(),
            target: h.id(),
            offset_uhz,
        };
        // Oversampling: the highest |frequency| z can carry above the responses'
        // stopbands, as cycles per native sample.
        let f_stop_phi = phi.spec().stopband_edge_uhz as f64 * 1e-6;
        let f_stop_h = h.spec().stopband_edge_uhz as f64 * 1e-6;
        let f_max = f_stop_phi.min(offset_hz.abs() + f_stop_h);
        let mut log2_up = 0;
        let (nu_max, farrow_error) = loop {
            let nu = f_max / rate_hz / f64::from(1u32 << log2_up);
            let e = lagrange6_error(nu);
            if e <= FARROW_ERROR_BOUND || log2_up == MAX_LOG2_UP {
                break (nu, e);
            }
            log2_up += 1;
        };
        let up = 1i64 << log2_up;

        // Time support of K, native samples: supp φ + supp h.
        let (h_lo, h_hi) = h.support_s();
        let t_phi = phi.half_width_s();
        let t_lo = (h_lo - t_phi) * rate_hz;
        let t_hi = (h_hi + t_phi) * rate_hz;
        let i_lo = libm::floor(t_lo) as i64 - 1;
        let i_hi = libm::ceil(t_hi) as i64;
        let taps_per_branch = i_hi - i_lo + 1;
        if taps_per_branch > MAX_KERNEL_TAPS {
            return Err(ContinuousError::KernelLength);
        }

        // Integration band: where both responses are non-negligible.
        let lo = (-band_edge(phi)).max(-offset_hz - band_edge(h)) / rate_hz;
        let hi = band_edge(phi).min(-offset_hz + band_edge(h)) / rate_hz;
        let mut taps = vec![[0.0, 0.0]; (up * taps_per_branch) as usize];
        if hi > lo {
            let t_abs = t_lo.abs().max(t_hi.abs());
            // Half a cycle of e^{j2πνt} per panel, and at least 8 panels per σ.
            let sigma = phi.sigma_hz().min(h.sigma_hz()) / rate_hz;
            let base = libm::ceil((hi - lo) * t_abs * 2.0).max(libm::ceil(8.0 * (hi - lo) / sigma));
            let mut panels = (base as u32).max(16);
            let eval = |panels: u32| -> Vec<[f64; 2]> {
                let nodes = band_nodes(phi, h, rate_hz, offset_hz, lo, hi, panels);
                let mut v = Vec::with_capacity(taps.len());
                for j in 0..up {
                    for i in i_lo..=i_hi {
                        let t = j as f64 / up as f64 + i as f64;
                        v.push(if t < t_lo || t > t_hi {
                            [0.0, 0.0]
                        } else {
                            kernel_at(&nodes, t)
                        });
                    }
                }
                v
            };
            let mut coarse = eval(panels);
            let mut converged = false;
            for _ in 0..6 {
                let fine = eval(2 * panels);
                let peak = fine
                    .iter()
                    .map(|v| (v[0] * v[0] + v[1] * v[1]).sqrt())
                    .fold(0.0, f64::max);
                let diff = coarse
                    .iter()
                    .zip(&fine)
                    .map(|(a, b)| {
                        let (dr, di) = (a[0] - b[0], a[1] - b[1]);
                        (dr * dr + di * di).sqrt()
                    })
                    .fold(0.0, f64::max);
                coarse = fine;
                panels *= 2;
                if diff <= KERNEL_QUADRATURE_BOUND * peak.max(1e-300) {
                    converged = true;
                    break;
                }
            }
            if !converged {
                return Err(ContinuousError::Quadrature);
            }
            taps = coarse;
        }
        let l1_max = (0..up as usize)
            .map(|j| {
                taps[j * taps_per_branch as usize..(j + 1) * taps_per_branch as usize]
                    .iter()
                    .map(|v| (v[0] * v[0] + v[1] * v[1]).sqrt())
                    .sum::<f64>()
            })
            .fold(0.0, f64::max);
        Ok(Self {
            key,
            log2_up,
            i_lo,
            taps_per_branch: taps_per_branch as usize,
            taps,
            l1_max,
            farrow_error,
            nu_max,
        })
    }

    /// The kernel's identity.
    pub fn key(&self) -> ResampleKey {
        self.key
    }

    /// log2 of the stage-1 oversampling U.
    pub fn log2_up(&self) -> u32 {
        self.log2_up
    }

    /// max over branches of Σ|K|: |z| ≤ l1_max·max|x|.
    pub fn l1_max(&self) -> f64 {
        self.l1_max
    }

    /// The Farrow error bound at the selected U ([`lagrange6_error`] at ν_max).
    pub fn farrow_error(&self) -> f64 {
        self.farrow_error
    }

    /// ν_max, cycles per z sample.
    pub fn nu_max(&self) -> f64 {
        self.nu_max
    }

    /// Native samples of lookahead and history a z point needs: z\[m\] reads
    /// x\[⌊m/U⌋ − i\] for i ∈ [`Self::tap_range`].
    pub fn tap_range(&self) -> (i64, i64) {
        (self.i_lo, self.i_lo + self.taps_per_branch as i64 - 1)
    }

    /// K(j/U + i) for branch j and native tap i (0 outside the support).
    pub fn tap(&self, j: u32, i: i64) -> [f64; 2] {
        let k = i - self.i_lo;
        if (j as usize) >= (1usize << self.log2_up) || k < 0 || k >= self.taps_per_branch as i64 {
            return [0.0, 0.0];
        }
        self.taps[j as usize * self.taps_per_branch + k as usize]
    }

    /// Native sample range `[first, last]` that z points `[m_lo, m_lo + len)` read.
    pub fn source_span(&self, m_lo: i64, len: usize) -> (i64, i64) {
        let (i_lo, i_hi) = self.tap_range();
        let m_hi = m_lo + len as i64 - 1;
        let n_lo = m_lo.div_euclid(1 << self.log2_up);
        let n_hi = m_hi.div_euclid(1 << self.log2_up);
        (n_lo - i_hi, n_hi - i_lo)
    }

    /// Stage 1: z points `[m_lo, m_lo + len)` (absolute grid index, U per native
    /// sample), from the extended source, each summed in f64 over increasing native tap i
    /// and rounded once to f32. `scale` multiplies every point (applied in f64 before
    /// rounding).
    pub fn filter_window(
        &self,
        timeline: &SourceTimeline,
        src: &mut dyn ContinuousSource,
        m_lo: i64,
        len: usize,
        scale: f64,
    ) -> Result<Vec<[f32; 2]>, ContinuousError> {
        if len == 0 {
            return Err(ContinuousError::Range);
        }
        let (first, last) = self.source_span(m_lo, len);
        let mut x = vec![
            [0.0, 0.0];
            usize::try_from(last - first + 1).map_err(|_| ContinuousError::Range)?
        ];
        if !timeline.all_zero(first, last) {
            timeline.read(src, first, &mut x)?;
        }
        let up = 1i64 << self.log2_up;
        let (i_lo, _) = self.tap_range();
        let mut out = Vec::with_capacity(len);
        for d in 0..len as i64 {
            let m = m_lo + d;
            let n = m.div_euclid(up);
            let j = m.rem_euclid(up) as usize;
            let row = &self.taps[j * self.taps_per_branch..(j + 1) * self.taps_per_branch];
            let (mut re, mut im) = (0.0_f64, 0.0_f64);
            for (k, t) in row.iter().enumerate() {
                // x index n − i with i = i_lo + k, as an offset into the buffer.
                let xi = (n - (i_lo + k as i64) - first) as usize;
                let v = x[xi];
                re += v[0] * t[0] - v[1] * t[1];
                im += v[0] * t[1] + v[1] * t[0];
            }
            let (re, im) = (re * scale, im * scale);
            if !(re.is_finite() && im.is_finite())
                || re.abs() > f64::from(f32::MAX)
                || im.abs() > f64::from(f32::MAX)
            {
                return Err(ContinuousError::Range);
            }
            out.push([re as f32, im as f32]);
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rate(hz: u64) -> SampleRate {
        SampleRate::from_hz(hz).unwrap()
    }

    fn spec(bw_hz: i64, stop_hz: i64) -> ResponseSpec {
        ResponseSpec {
            passband_edge_uhz: bw_hz * 500_000,
            stopband_edge_uhz: stop_hz * 500_000,
            stopband_atten_mdb: 80_000,
            group_delay_ps: 0,
        }
    }

    fn phi_2m() -> TapResponse {
        source_response(rate(2_000_000), spec(1_600_000, 2_000_000)).unwrap()
    }

    fn rx_16m(gd_ps: i64) -> TapResponse {
        TapResponse::new(
            TapKind::VirtualReceiver,
            ResponseSpec {
                group_delay_ps: gd_ps,
                ..spec(12_000_000, 15_000_000)
            },
            rate(16_000_000),
            1_090_000_000_000_000,
        )
        .unwrap()
    }

    #[test]
    fn lagrange_error_matches_the_design_table() {
        // Values computed independently in the planning notes.
        assert!((lagrange6_error(0.5 / 16.0) - 2.78e-7).abs() < 0.02e-7);
        assert!((lagrange6_error(0.5 / 8.0) - 1.75e-5).abs() < 0.02e-5);
        assert!(lagrange6_error(0.0) < 1e-15);
    }

    #[test]
    fn kernel_spectrum_is_phi_times_shifted_h() {
        // A tone through stage 1 has gain Φ(f)·H(f + Δf): DTFT of each branch at ν.
        let phi = phi_2m();
        let h = rx_16m(1_000_000);
        let offset = 5_000_000_000_000; // +5 MHz
        let k = ResampleKernel::build(&phi, &h, offset).unwrap();
        assert!(k.log2_up() <= 4 && k.farrow_error() <= FARROW_ERROR_BOUND);
        let (i_lo, i_hi) = k.tap_range();
        for f in [-700_000.0, 0.0, 300_000.0, 790_000.0] {
            let nu = f / 2e6;
            for j in [0u32, 1] {
                let frac = f64::from(j) / f64::from(1u32 << k.log2_up());
                // z(n + frac) for x[n] = e^{j2πνn}: Σ_i e^{j2πν(n−i)}K(frac+i)
                // = e^{j2πν(n+frac)}·G with G = Σ_i K(frac+i)e^{−j2πν(i+frac)}.
                let (mut re, mut im) = (0.0, 0.0);
                for i in i_lo..=i_hi {
                    let t = k.tap(j, i);
                    let (s, c) =
                        libm::sincos(-2.0 * core::f64::consts::PI * nu * (i as f64 + frac));
                    re += t[0] * c - t[1] * s;
                    im += t[0] * s + t[1] * c;
                }
                let (hr, hi) = h.response(f + 5e6);
                let a = phi.magnitude(f);
                let (er, ei) = (re - a * hr, im - a * hi);
                assert!((er * er + ei * ei).sqrt() < 2e-7, "f = {f}, j = {j}");
            }
        }
    }

    #[test]
    fn out_of_band_kernel_is_zero() {
        let phi = phi_2m();
        let h = rx_16m(0);
        // 60 MHz away: no overlap at all.
        let k = ResampleKernel::build(&phi, &h, -60_000_000_000_000).unwrap();
        assert_eq!(k.l1_max(), 0.0);
    }

    #[test]
    fn timeline_padding_and_loop() {
        let mut src = MemorySource((0..5).map(|i| [f64::from(i), 0.0]).collect());
        let mut out = vec![[9.0, 9.0]; 14];
        let lp = SourceTimeline {
            len: 5,
            policy: EndPolicy::Loop,
        };
        lp.read(&mut src, -2, &mut out).unwrap();
        let re: Vec<f64> = out.iter().map(|v| v[0]).collect();
        assert_eq!(re, [0., 0., 0., 1., 2., 3., 4., 0., 1., 2., 3., 4., 0., 1.]);
        let zp = SourceTimeline {
            len: 5,
            policy: EndPolicy::ZeroPad,
        };
        zp.read(&mut src, 3, &mut out).unwrap();
        assert_eq!(out[0][0], 3.0);
        assert!(out[2..].iter().all(|v| *v == [0.0, 0.0]));
        assert!(zp.all_zero(5, 100) && zp.all_zero(-9, -1) && !lp.all_zero(5, 6));
        let st = SourceTimeline {
            len: 5,
            policy: EndPolicy::Stop,
        };
        assert!(!st.exhausted(4 << 32) && st.exhausted((4 << 32) + 1));
        assert!(!zp.exhausted(1 << 40));
    }

    #[test]
    fn windows_agree_bit_for_bit_where_they_overlap() {
        let phi = phi_2m();
        let h = rx_16m(250_000);
        let k = ResampleKernel::build(&phi, &h, 1_234_000_000).unwrap();
        let x: Vec<[f64; 2]> = (0..400)
            .map(|n| {
                let t = f64::from(n);
                [libm::cos(0.3 * t + 0.001 * t * t), libm::sin(0.2 * t)]
            })
            .collect();
        let mut src = MemorySource(x);
        let tl = SourceTimeline {
            len: 400,
            policy: EndPolicy::Loop,
        };
        let a = k.filter_window(&tl, &mut src, -40, 3000, 1.0).unwrap();
        let b = k.filter_window(&tl, &mut src, 1111, 777, 1.0).unwrap();
        assert_eq!(&a[1151..1151 + 777], &b[..]);
        // Looping: z repeats with period N·U once both points are past the onset.
        let p = 400 * (1usize << k.log2_up());
        let c = k
            .filter_window(&tl, &mut src, 1111 + p as i64, 777, 1.0)
            .unwrap();
        assert_eq!(b, c);
    }
}
