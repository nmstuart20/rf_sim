# iqsim contracts

Shared contracts every crate and backend must honor. They come from `project_plan.md` ("Time, numerics and determinism", "Receiver model", "Compute backends and scaling", "Interfaces", "Validation strategy").

**Status:** numeric contract version **0 (proposed)**. Values marked *proposed* are starting points. Phase 1 (prompts 1–4) fixes them before any render code lands, and each change is recorded here and in `docs/progress.md`. The version becomes 1 when §10 is fully settled. After that, any change that can alter output bits increments it, and the version is written into every run manifest.

Contents:

1. Numeric representation and operation rules
2. Amplitude normalization, accumulation and overflow policy
3. Clocks: ClockMap, timestamp clock, LO model
4. Exclusive commit boundary
5. Tap points: responses and state ownership
6. Run manifest and replay
7. Transport choice, framing and delivery guarantees
8. Resource budgeting
9. Crate placement and dependency graph
10. Proposed tolerances and decisions phase 1 must settle
11. Scenario quantities and files

---

## 1. Numeric representation and operation rules

### 1.1 Integer quantities

| Quantity | Type | Unit / scale | Range | Overflow |
| --- | --- | --- | --- | --- |
| Truth time `TruthTime` | i64 | picoseconds from scenario epoch | ±2⁶³ ps ≈ ±106.7 days | checked; error |
| Sub-ps fraction (`AcqInstant.frac`) | u64 | 2⁻⁶⁴ ps, floor of the exact instant (§3.1) | [0, 1) ps | n/a |
| Exact index-time remainder (`ExactTime.rem`) | u64 over the reduced rate numerator | fraction of 1 ps | [0, 1) ps | n/a |
| Sample index `SampleIndex` | u64 | samples, per receiver channel, from channel epoch | 2⁶⁴ | checked; error |
| Sample rate `SampleRate` | (u64 num, u64 den), reduced, nonzero | samples per second = num/den | den ≤ 2³², num ≤ 2⁴⁸, num/den ≤ 10⁹ | rejected at construction |
| Path delay (output) `Delay` | i64 Q32.32 | output samples | ±2³¹ samples (2.1 s at 1 GS/s) | range-checked at fit time |
| Carrier phase (output) `Phase` | u64 | 2⁻⁶⁴ cycle | one cycle | **wraps by design** |
| Carrier-offset frequency `OffsetFrequency` | i64 | 2⁻³² Hz | ±2³¹ Hz ≈ ±2.1 GHz | checked at construction; error |
| Delay/phase polynomial coefficients | 128-bit, Newton forward differences D0..D3 | output scale + 64 extra fractional bits | n/a | **wraps by design** (§1.3) |
| Accumulator | i64 per I and Q | per-target LSB voltage (§2) | ±(2⁶³−1) | **forbidden**; error |

Wrapping is allowed only for carrier phase and the Newton-form polynomial registers, where arithmetic mod 2¹²⁸ (or mod 2⁶⁴ for phase) is the definition. Every other integer operation on time, index, delay output or accumulators is checked, and overflow returns an error. It never wraps or saturates silently. Release builds keep `overflow-checks = true`.

The split is visible in the code (settled, prompt 2):

- Intentional wrapping lives only in `iqsim_core::wide` (`U128`, every method named `wrapping_*`), in `BlockPoly`/`PolyCursor`, in `Phase::wrapping_add`/`wrapping_sub`, and in `OffsetFrequency::phase_at`'s final reduction mod 2⁶⁴.
- `Delay` has no arithmetic operators. Converting to it is range-checked (`NumericError::OutOfRange`).
- Accumulators use plain i64 with `overflow-checks = true` on CPU, and the headroom proof of §2 on every backend. A wrapping op on an accumulator is a contract violation.

### 1.2 Time and index conversions

- Ideal map: `t(n) = channel_epoch + n · den / num` seconds, evaluated in ps as `n · den · 10¹² / num` with u128 intermediates. The quotient is the integer ps, and the remainder `(n · den · 10¹²) mod num` over `num` is the exact sub-ps fraction.
- Time → index: floor (`ChannelTimebase::index_at_or_before`: largest n with exact t(n) ≤ t; a time before the epoch is `BeforeEpoch`). A ceiling companion (`index_at_or_after`) returns the smallest n with t(n) ≥ t, and 0 for times before the epoch.
- Index → time (`ChannelTimebase::index_to_time`): returns `ExactTime { floor, rem, rem_den = num }`, exact for every index. It is an integer ps time when `rem = 0`. Converting to an `AcqInstant` floors `rem / num` onto the 2⁻⁶⁴ ps grid. The fraction is never dropped on the geometry path.
- **Rate bounds (settled, prompt 1):** after reduction, `den ≤ 2³²`, `num ≤ 2⁴⁸` and `num/den ≤ 10⁹` (1 GS/s). `SampleRate::from_ratio` rejects anything else, and zero numerators or denominators.
- **No-overflow proof (prompt 1).** Let `M = den · 10¹² < 2⁷²`.
  - Index → time: `n · den < 2⁹⁶` always. The product with 10¹² is checked. For runs ≤ 100 days at ≤ 1 GS/s, `n ≤ 8.64·10¹⁵ < 2⁵³`, so it is `< 2¹²⁵` and never overflows. If it does overflow (far larger n), the product is ≥ 2¹²⁸ and the quotient ≥ 2¹²⁸ / 2⁴⁸ = 2⁸⁰, which is not a valid i64 ps time anyway. So the check never rejects a representable result, and the error is `TimeOverflow` either way. The epoch is added in i128 and range-checked into i64.
  - Time → index: `d = t − epoch` is computed in i128, `|d| < 2⁶⁴`. `d · num < 2¹¹²`. The quotient is ≤ d / 1000 < 2⁵⁴, so it always fits in u64.
  - Remainder → fraction: `rem < num ≤ 2⁴⁸`, so `rem · 2⁶⁴ < 2¹¹²`.
  - A property test checks that `index_to_time` fails exactly when the true time does not fit in i64, for arbitrary u64 indices, epochs and accepted rates.

### 1.3 Block polynomials (delay and phase)

- Order ≤ 3 in local sample index k ∈ [0, block_len). Stored as Newton forward differences with 64 fractional bits beyond the output format: delay coefficients are Q32.96, phase coefficients are 0.128 cycle fractions (`DelayPoly`, `PhasePoly` = `BlockPoly<Delay>`, `BlockPoly<Phase>`).
- Jump: `p(k) = D0 + C(k,1)·D1 + C(k,2)·D2 + C(k,3)·D3` mod 2¹²⁸, with `U128 × u64` wrapping multiplies. Step: the cursor holds `(p, Δ1, Δ2, D3)` with `Δ1(k) = D1 + C(k,1)·D2 + C(k,2)·D3` and `Δ2(k) = D2 + C(k,1)·D3`, and a step is `p += Δ1; Δ1 += Δ2; Δ2 += D3`, in that order. Jump and step agree bit for bit at every k.
- **Binomials (settled, prompt 2):** `C(k,1) = k`, `C(k,2) = k(k−1)/2`, and `C(k,3) = C(k,2)·((k−2)/3)` if 3 | k−2, otherwise `(C(k,2)/3)·(k−2)`. Every division is exact and no intermediate exceeds C(k,3), so all of it is u64 arithmetic.
- Output extraction keeps the top 64 bits with an arithmetic shift (floor toward −∞). It does not round.
- **Maximum block length (settled, prompt 2):** `MAX_BLOCK_LEN = 4 801 281`. The largest local index is 4 801 280, the largest k with C(k,3) < 2⁶⁴. That is 4.8 ms at 1 GS/s and 19.2 ms at 250 MS/s. `BlockPoly` rejects longer blocks. `BlockSpec` itself stays unbounded, and the render configuration must choose `block_len ≤ MAX_BLOCK_LEN`.
- **Coefficient quantization bound (settled, prompt 2):** rounding each Dj to the nearest register value errs by ≤ ½ register unit (2⁻⁶⁵ output LSB). At index k the register error is ≤ ½·(1 + C(k,1) + C(k,2) + C(k,3)) units. That is < 2⁶³ + 2⁴³ units, i.e. ½ + 2⁻²¹ output LSB, for every k < `MAX_BLOCK_LEN`. Floor extraction adds < 1 LSB, so the output is within 1.5 + 2⁻²¹ LSB (under 2 LSB) of the exact polynomial:
  - delay: < 2·2⁻³² ≈ 4.7·10⁻¹⁰ samples (budget 10⁻³)
  - phase: < 2·2⁻⁶⁴ cycle ≈ 6.9·10⁻¹⁹ rad (budget 10⁻²)
  - For a 1 ms block at 500 MS/s, C(k,3) < 2⁵⁵, so the register error is < 2⁻⁹ LSB.

  Quantization therefore uses none of the §10 delay or phase budgets in practice. All of the budget goes to the geometry fit (phase 3). A realistic 61.44 MS/s delay polynomial measures 1.0 LSB worst case, including the floor.
- `BlockPoly::from_f64_differences` rounds each Dj ties-to-even: Q32.96 with a range check for delay, and mod one cycle for phase (exact, because C(k,j) is an integer).
- The 128-bit arithmetic is written with explicit 64-bit primitives (`mul_hi_u64` from four 32 × 32 products, wrapping mul, `add_with_carry`, `sub_with_borrow`) in `iqsim_core::wide`, so the CUDA C port mirrors it line by line (`__umul64hi`, `*`, carry compares). It is tested against native u128/i128.
- **Delay range margin (settled, prompt 2):** the geometry engine checks that the delay polynomial's true range over the block stays inside [−2³¹ + m, 2³¹ − m) samples, with m = `DELAY_RANGE_MARGIN_LSB` = 2 Q32.32 LSB. The margin matters because quantization (< 2 LSB) can push a true delay that is just inside the range across −2³¹, and the wrapped output then reads as nearly `i64::MAX`. Leaving the margin is a fit error.

### 1.4 Floating point

Floating point is allowed only where its rounding is fixed by construction.

- **Render path:** f32 only, for amplitude polynomials, sin/cos results, waveform samples and per-path complex products. Never f64 (T1000 FP64 is 1:32). Never time, delay or carrier phase.
- **Control plane:** f64 for geometry, fitting, link budgets and attitude. It produces integer or f32 coefficients for the work list. Transcendentals come from the `libm` crate, never std (clippy enforces this). Results must be identical on x86-64 (RHEL 9 target) and aarch64 (development). Only IEEE basic operations plus libm are allowed. No x87, no `target-cpu=native`.
- **Operation order:** the source order is the order. No reassociation, no fast-math, no FMA contraction. `mul_add` is banned. CUDA kernels build with `--fmad=false`.
- **Rounding mode:** IEEE 754 round-to-nearest-even for all arithmetic.
- **Division and sqrt:** correctly rounded IEEE on both backends. CUDA uses `-prec-div=true -prec-sqrt=true`.
- **Subnormals (settled, prompt 2):** preserved (IEEE gradual underflow) on both backends. Nothing ever flushes. CPU code never sets FTZ/DAZ (Rust's default floating-point environment), and CUDA uses `-ftz=false`. The cost is negligible because render-path values stay far from the subnormal range: `sincos` cannot produce one (smallest nonzero intermediate ≈ 5·10⁻²⁹), and amplitudes in volts at −200 dBm are ≈ 2·10⁻¹¹. Preserving is also the only choice that needs no mode bit to match between backends. If a later kernel can produce subnormals on a hot path, it must measure the cost and propose a change here.
- **Float → integer (settled, prompt 2):** ties-to-even, never a bare Rust `as` (which truncates and saturates). The rule is `round_ties_even()`, then a range check, then `as`, which is exact at that point. On CUDA, the equivalents are `__float2ll_rn`/`__double2ll_rn` plus the same range check. Out-of-range is an error. Implementations: `f64_to_i64_ties_even`, `f32_to_i64_ties_even` (the f32 → f64 widening is exact), `Delay::from_samples_f64`, `Phase::from_cycles_f64` and `BlockPoly::from_f64_differences`. Scaling by 2ⁿ (0 ≤ n ≤ 128) before rounding is an exact multiply, because scaling up never discards significand bits, even for subnormal inputs or products. An overflow to ±Inf is rejected. Wrapping (cycle) quantities remove whole cycles first with `x − trunc(x)`, which is exact.
- **Integer → float:** round to nearest even (Rust `as`; CUDA `__ll2float_rn`/`__ull2float_rn` or a plain C cast). Render code converts only values with ≤ 24 significant bits to f32, so these conversions are exact. Example: the 24-bit sincos remainder.
- **Non-finite values:** NaN and ±Inf are rejected at every input validation and at every float → integer boundary (`require_finite_f64`/`require_finite_f32`, `NumericError::NonFinite`). A non-finite value never reaches an accumulator.
- **Compiler settings (settled, prompt 2):**
  - CPU: the toolchain is pinned (`rust-toolchain.toml`), with no nightly features. Rust never contracts a·b + c into an FMA, even when the target has FMA, and stable Rust has no fast-math flag. `mul_add` is banned by clippy. `-C target-cpu=native` is banned, and the release target is plain `x86_64-unknown-linux-gnu`.
  - CUDA: `nvcc --fmad=false -ftz=false -prec-div=true -prec-sqrt=true`, never `--use_fast_math`, and no `__f*` fast intrinsics (`__sinf`, `__fdividef`, …). Fatbins are built offline for sm_75 and sm_86 with these flags, and their hashes go in the run manifest.

### 1.5 sin/cos, ln, Gaussian

- `sincos(phase: Phase) -> (f32, f32)` (settled, prompt 2). The full method is in `iqsim-core/src/trig.rs`.
  - The top 10 phase bits index `SIN_TABLE` (1024 f32 entries; cos uses the entry a quarter cycle ahead). The next 24 bits form δ = r·(2π_f32·2⁻³⁴), which is exact in f32. The low 30 bits are dropped, an angle error < 3.7·10⁻¹⁰ rad.
  - Correction, with every operation separate and in source order: `sd = δ − δ·(δ²·(1/6))`, `h = δ²·0.5`, `sin = s + (c·sd − s·h)`, `cos = c − (s·sd + c·h)`.
  - Table construction runs at compile time and uses only integers, no float library. sin((π/2)·j/256) for j = 0..256 comes from an integer Taylor series in Q2.62 starting at ⌊(π/2)·2⁶²⌋ = `0x6487ED5110B4611A`. Each value is rounded to f32 ties-to-even by bit manipulation, and quarter-wave symmetry fills the table, with zero always +0.0. A test checks that every entry is correctly rounded and that none is within 2⁸ Q2.62 units of a rounding tie.
  - Measured max absolute error over all 2²⁴ values of the top 24 bits, with random low bits: 5.98·10⁻⁸ (sin) and 5.97·10⁻⁸ (cos), about 1 ulp. The target is ≤ 3·10⁻⁷.
  - Golden hashes: the table bytes (`core/sincos_table`) and 10⁶ evaluations (`core/sincos_1e6`).
- `ln(x: f32) -> Result<f32, NumericError>` (settled, prompt 3). It uses no std or libm call, so the GPU port matches. The full method is in `iqsim-core/src/log.rs`.
  - Range reduction on the bit pattern, with integer operations only. A subnormal is first scaled by 2²³ (exact). `x = 2ᵉ·m`, and m is halved when its bit pattern is above `0x3FB504F3`, so m ∈ [0.7071068, 1.4142135].
  - `f = m − 1` (exact), `s = f/(2 + f)`, `z = s²`, `hfsq = 0.5·f²`, `R = z·(2/3 + z·(2/5 + z·(2/7 + z·2/9)))` (Taylor coefficients rounded to f32, Horner form, every operation separate).
  - `ln x = e·LN2_HI − ((hfsq − (s·(hfsq + R) + e·LN2_LO)) − f)`. `LN2_HI` is ln 2 rounded to f32 with the low 12 mantissa bits cleared, so `e·LN2_HI` is exact for every exponent. `LN2_LO = (ln 2 − LN2_HI)` rounded to f32.
  - Error, measured against an f64 reference over every f32 in 13 binades plus a sample of subnormals: **≤ 0.86 ulp** (asserted ≤ 1 ulp, and ≤ 1 ulp in a property test over random positive finite f32). On [2⁻⁶⁵, 1], the Box-Muller range, the absolute error is ≤ 1.932·10⁻⁶, asserted ≤ 2·10⁻⁶ (half an ulp of |ln 2⁻⁶⁵| ≈ 45). `ln(1) = +0`. Zero, negative and non-finite inputs are errors.
  - Golden: 10⁶ evaluations (`core/ln_1e6`).
- **Philox (settled, prompt 3).** Philox4x32-10 exactly as in Random123 (`iqsim-core/src/philox.rs`): multipliers `0xD2511F53`, `0xCD9E8D57`, Weyl key bumps `0x9E3779B9`, `0xBB67AE85`, ten rounds and nine bumps. It passes the six philox4x32 known-answer vectors (7 and 10 rounds) of Random123's `tests/kat_vectors`, which are checked in verbatim with their provenance. The 32 × 32 → 64-bit products are plain u64 multiplies (CUDA: `__umulhi` and `*`).
- **Stream key (settled, prompt 3).** `mix(z)` is the SplitMix64 finalizer (`z ^= z>>30; z *= 0xBF58476D1CE4E5B9; z ^= z>>27; z *= 0x94D049BB133111EB; z ^= z>>31`, mod 2⁶⁴), a bijection on u64. With `KEY_DOMAIN = 0x697173696D2D7031` (ASCII `iqsim-p1`):
  - `k = mix(mix(seed ⊕ KEY_DOMAIN) ⊕ (entity << 32 | stream))`, and the key is `(k mod 2³², k >> 32)`.
  - For a fixed seed the map (entity, stream) → key is injective, so no two streams of one scenario share a key. A property test checks this through an explicit inverse of `mix`.
- **Counter (settled, prompt 3).** `(c0, c1, c2, c3) = (index_lo, index_hi, sub_lo, sub_hi)`: the u64 sample index at the rate where the draw happens, then a u64 sub-index. Each model documents which sub-indices it uses.
- **Stream-id registry (settled, prompt 3).** `StreamId = kind << 16 | instance`, where `StreamKind` is a u16 code. `StreamId` can only be built from a registered kind, so two models never share a stream. Codes are part of this contract: new kinds are appended, and no code is ever renumbered or reused, because that would change output bits. Code 0 is reserved. Registered so far:

  | Code | Kind | Instance | Entity | Sub-indices |
  | --- | --- | --- | --- | --- |
  | 1 | `ReceiverThermalNoise` | channel index | receiver | 0 (`NOISE_SUB_INDEX`) |

- **Gaussian noise (settled, prompt 3).** `noise(seed, entity, stream, n) -> (f32, f32)` (or `NoiseStream::sample`) is circularly symmetric complex Gaussian with **unit complex power**: `E[I²] = E[Q²] = ½`, `E[I·Q] = 0`, `E|z|² = 1`. Scale it by the RMS voltage of the wanted noise power (§2). Sample n uses one block at counter (n, 0):
  1. `x = w1·2³² + w0`. The radius uniform `u` is `x·2⁻⁶⁴` truncated to 24 significant bits, built as f32 from the bit pattern (leading-zero count, shift, exponent field). `x = 0` gives `u = 2⁻⁶⁵`. So `u ∈ [2⁻⁶⁵, 1 − 2⁻²⁴]`.
  2. `r = sqrt(−ln u)` (IEEE sqrt). `r ≤ sqrt(65·ln 2) ≈ 6.71`, about 9.5 σ per component.
  3. `(s, c) = sincos(Phase(w2 << 32))`: the angle is w2 as a u32 phase.
  4. Output `(r·c, r·s)`. Word w3 is unused.

  Quantization of u and θ and the ln and sincos errors bias the variance by under 10⁻⁶ relative. Tests over fixed seeds:
  - 10⁷ samples: means, variances, I/Q correlation, skewness and excess kurtosis all within 6 standard errors.
  - Anderson-Darling A² against the 1 % critical value 3.857: I 1.153, Q 0.282 (scaled to N(0, 1)) and |z|² 0.390 (against Exp(1)).
  - Complex cross-correlation at lags −2..2, and power correlation, for adjacent channels, entities and seeds, and for the same stream at lags 1..5. All are within 6 standard errors over 10⁶ samples.
  - Golden: 10⁶ samples (`core/noise_1e6`).

### 1.6 Carrier-offset phase

Offset phase `φ(t) = frac(Δf · t)` comes from exact time with modular arithmetic. Nothing ever computes a large floating-point `f · t` and discards its integer cycles.

**Settled (prompt 2).**

- Δf is `OffsetFrequency`, an i64 in units of 2⁻³² Hz: range ±2³¹ Hz, resolution 2.3·10⁻¹⁰ Hz. It is built from whole Hz (exact), from integer µHz (`from_uhz`: `uhz·2³²/10⁶` rounded ties-to-even in integers, the required path for scenario frequencies, prompt 5), or from f64 Hz (ties-to-even, range-checked). Going µHz → f64 Hz → `OffsetFrequency` adds a rounding that can reach 0.04 rad over 100 days at 7 MHz, so the f64 path is only for values that are already f64.
- t is the Q64.64 ps `AcqInstant` T. The result is `φ = ⌊F·T / (2⁴⁴·5¹²)⌋ mod 2⁶⁴`, in units of 2⁻⁶⁴ cycle, floored toward −∞. It is exact for every i64 frequency and every instant, with error in [0, 1) unit.
- `OffsetFrequency::phase_at` uses 64-bit primitives only, so the CUDA port can mirror it:
  1. a 192-bit product |F|·|T|
  2. a shift right by 44
  3. schoolbook division by 5¹² over 32-bit digits (each partial dividend is < 2⁶⁰)
  4. for a negative product, `−⌈A/D⌉` mod 2⁶⁴
- Tests:
  - an independent bit-serial big-integer reference
  - additivity in t and in Δf to within one floor unit, with t up to ±104 days
  - for every whole-Hz offset, φ(100 days + δ) equals φ(δ) exactly
  - a golden over 10⁶ evaluations (`core/offset_phase_1e6`)

## 2. Amplitude normalization, accumulation and overflow policy

**Physical normalization (settled, prompt 2).** A sample value is the complex-envelope voltage at a reference impedance R₀ = 50 Ω, using the RMS convention `P = |v|² / R₀`. A CW tone of power P has `|v| = √(P·R₀)`. Path amplitude polynomials carry path loss, antenna gains, polarization mismatch and the selected tap response. They are in volts at the render target's reference plane.

**Per-target integer scale.** Each render target (one tap of one channel) declares:

- reference plane: antenna port of the receiver, or the virtual ADC input
- full-scale voltage `V_fs`
- LSB voltage `q = V_fs / 2^S`, with S = 40 (settled, prompt 2). One LSB is −240 dBFS, so the f32 contribution's 24-bit precision, not the integer grid, limits small signals.

Contribution component value = `round_ties_even(component_volts / q)` as i64 (`f32_to_i64_ties_even` after the f32 scaling, §1.4). The antenna-port scale is independent of any virtual ADC full scale.

**Scale factor (settled, prompt 5).** `TargetScale::new(V_fs)` computes `lsb_per_volt = f32(2⁴⁰ / V_fs)` (one f64 division, one rounding to f32) and requires it to be a normal f32. **The f32 factor defines the scale**: q is exactly `1 / lsb_per_volt`, and every backend multiplies by that same f32. Each target (antenna port, virtual receiver) declares its own V_fs. Implementation: `iqsim_core::amplitude`.

**Supported peak overlap (settled, prompt 5).** `FULL_SCALE_PATHS = 2⁶² / 2⁴⁰ = 2²²` simultaneous full-scale contributions per sample. A path's bound is `max|a_p|·max|s_p|·lsb_per_volt`:
- `max|a_p|`: Σ|c_j| of the f32 amplitude cubic in x ∈ [0, 1], widened by 2⁻²⁰ relative for the f32 Horner evaluation.
- `max|s_p|`: taken as at least 1, so the proof also bounds the kernel's intermediate `m = a·lsb_per_volt` (a tiny pulse amplitude on a strong path could otherwise pass while m overflows f32). 1 for CW. Pulse amplitudes are limited to |A| ≤ `MAX_PULSE_AMPLITUDE` = 10⁶, which keeps the envelope sums finite. For a pulse train, the largest sum of `A_i·(peak_rise + peak_fall)` over the pulses whose extents overlap one point. A pulse's extent runs from the earlier of its two shifted table starts to the **later** of its two table ends (`EdgeTableArena::pulse_extent`), because a long rise can outlast the fall's table. The peak of a table also bounds every interpolated value, because linear interpolation is a convex combination.

**Culling (settled, prompt 5).** A path is dropped from a block's work list when its peak contribution is below the target's cull threshold, in LSB. The default `LOSSLESS_CULL_LSB = ¼` is **lossless**: every computed component is then below ¼·(1 + 2⁻²⁰) < ½ LSB, which rounds to 0, so output bits are identical with and without culling. A test checks this. A larger threshold may be configured, and it is then an approximation of at most that many LSB per culled path. Narrowband paths far outside the band (|H(f_off)| ≈ 0) are culled this way.

**Headroom proof.** Before rendering, and again whenever the path set or amplitude bounds change, the control plane computes, per target and block:

```
Σ_paths  ceil(max|a_p| · max|s_e| / q)  ≤  L = 2⁶²
```

`max|a_p|` comes from the block's amplitude polynomial bound. `max|s_e|` is the waveform's declared peak (pulses: edge-table peak; continuous: declared crest factor). **Guard (settled, prompt 2):** the limit is L = 2⁶², one full guard bit below the i64 range. That slack absorbs the f32 rounding of each contribution (relative ≤ 2⁻²³ per operation, a few operations per path) and any looseness in the amplitude bound, with a factor of 2 to spare. With S = 40 and contributions bounded by full scale, 2²² full-scale paths fit per sample. A configuration that fails the proof is rejected (`sim lint` and run start). Within the bound, integer addition is associative, so partitioning and thread order cannot change the sum. The CPU backend also uses checked adds, so a proof bug fails loudly instead of wrapping.

**Reduction and finalization.** Emitter partitions return raw i64 I/Q sums. The host, or one designated device, reduces them exactly. Receiver effects and output conversion run once on the complete sum (§5).

**Output conversion.** sc16 or cf32 per sink. i64 → sc16 is `round_ties_even` of the scaled value with clip detection. Clipping is a reported receiver effect on the virtual ADC. On the antenna port, clipping is an error.

**Amplitude tolerance and null floor (settled, prompt 2, with the user).** Each link has a reference voltage `V_ref`, its nominal main-beam amplitude. For a block, `V_ref` is the amplitude the link's direct path would have at the block-start range, at the target's reference plane, with:

- both antennas at their peak gain
- matched polarization
- no blockage or diffraction loss
- the tap response's passband peak

It is computed once per block on the control plane and logged with the fit residuals. Tolerances:

- Where `|a| ≥ floor`, the amplitude error must be ≤ 0.01 dB relative (|Δa| ≤ (10^(0.01/20) − 1)·|a| ≈ 1.152·10⁻³·|a|). The floor is `V_ref · 10⁻²`, i.e. −40 dB.
- Below the floor, the absolute complex-gain error must be ≤ `10⁻⁵ · V_ref`.
- The two bounds are continuous: at the floor the relative bound is 1.152·10⁻⁵·V_ref, so the tolerance never loosens as |a| crosses the floor downward. The originally proposed 10⁻³·V_ref was rejected because it is 100× looser just below the floor than just above it.
- Phase is tested only above the floor.

## 3. Clocks: ClockMap, timestamp clock, LO model

Three distinct models. No code may use one as a stand-in for another.

### 3.1 ClockMap (sampling oscillator)

A ClockMap maps a sample index to the actual truth-time acquisition instant. Geometry, propagation delay and carrier phase are evaluated at these instants.

```text
trait ClockMap {
    /// Actual acquisition instant of sample n. Strictly increasing in n.
    fn instant(&self, n: SampleIndex) -> Result<AcqInstant>;
    /// Smallest n with instant(n) >= t (and the companion "largest n with instant(n) <= t").
    fn index_at_or_after(&self, t: AcqInstant) -> Result<SampleIndex>;
    /// Stable identity (content hash of the model + parameters) for cache keys and manifests.
    fn id(&self) -> ClockMapId;
}
AcqInstant = { ps: TruthTime, frac: u64 }   // ps + frac·2⁻⁶⁴ ps, i.e. signed Q64.64 ps (settled, prompt 1)
```

- **Ideal map:** the exact rational map of §1.2, floored onto the 2⁻⁶⁴ ps grid (`IdealClockMap`). This is the v1 implementation, and the first render code uses the trait, not a shortcut. Its identity is SHA-256 of `"iqsim.clockmap.ideal.v1\0" ‖ epoch_ps (i64 LE) ‖ num (u64 LE) ‖ den (u64 LE)`.
- **AcqInstant representation (settled, prompt 1):** `ps + frac · 2⁻⁶⁴` picoseconds, floor of the exact instant. A fixed binary grid works the same for the ideal and perturbed maps, and it is a plain 128-bit integer for the GPU port. The quantization error is < 2⁻⁶⁴ ps ≈ 5.4·10⁻²⁰ ps. At 12 GHz that is < 5·10⁻²¹ rad of carrier phase, so it uses none of the §10 phase or delay budget. Rounding keeps the map strictly increasing, because samples are ≥ 1000 ps apart.
- **Inverse:** `index_at_or_after(t)` returns the smallest n with `instant(n) ≥ t`, and `index_at_or_before(t)` the largest n with `instant(n) ≤ t` (`None` before sample 0). Both compare against the floored instants, so `index_at_or_after(instant(n)) = index_at_or_before(instant(n)) = n` exactly. The inverse returns only the index. A caller that needs the sub-sample offset computes `t − instant(n)` exactly in Q64.64. When no representable sample is at or after t, `index_at_or_after` is `TimeOverflow`.
- **Perturbed maps (prompt 15):** frequency offset, drift, aging, optional GPSDO discipline, and phase noise driven by Philox. They are deterministic and strictly monotone, so the inverse is well defined (integer search over a monotone function). Cumulative random processes are integrated on the control plane in block order and frozen into per-block coefficients. Evaluating a block never depends on partitioning or on which device renders it.
- **Coherent channels:** the M channels of one receiver share one oscillator. Each channel's index → instant goes through the shared oscillator time at shared physical times, so channels with different rates stay coherent.
- **Antenna-port tap:** uses the declared output-device timebase, which is the ideal map at the device's nominal rate. It bypasses every virtual clock impairment.

### 3.2 Timestamp clock

Maps an acquisition instant to the time reported to the consumer (stream headers, VITA-49 timestamps, SigMF `core:datetime`). It has its own offset and drift. It relabels only and never moves sampling instants. Truth exports record the truth instant, the reported timestamp and the sample index together.

### 3.3 LO model

Controls mixing phase and frequency error, i.e. the `e^{j2π(f_c − f_LO,r)t}` term and phase noise, evaluated with §1.6 modular arithmetic at actual acquisition instants. Coupling to the sampling oscillator (common reference) is explicit in the receiver profile. It is never implied.

## 4. Exclusive commit boundary

- **Horizons.** Render progress R (per channel, in sample indices) and the commit boundary C (truth time) are distinct. Telemetry reports R, C and command latency separately.
- **Meaning of C.** C is exclusive. Every truth-time dependency before C is immutable. C covers everything already consumed by behavior (including lookahead), geometry, work-list construction, receiver processing and filter lookahead.
- **Reservation.** A stage reserves its dependency interval `[start, end)` before advancing. Advancing C rounds fractional endpoints upward to the event-time grid (*proposed* grid: 1 ps).
- **Serialization.** Command admission and stage advancement are serialized. In v1 this is one commit actor, with no rollback.
- **Admission.** A command with requested time `t_req` gets effective time `t_eff = max(t_req, C)`. It is appended to the durable event log, fsync'd, and only then acknowledged. The ack carries `t_eff` and the ingress sequence. All of this completes before any stage advances across `t_eff`. Every command pending at a boundary is ordered before the advance.
- **Ordering.** Behavior events are ordered by `(truth time, entity id, sequence)`. *Proposed*: an operator command at time t precedes entity events at t, and commands are ordered among themselves by ingress sequence.
- **Latency.** Behavior lookahead counts toward the same latency budget as pipeline depth. *Proposed* defaults: block 1 ms, depth 50–100 ms.
- **In-flight signals.** A signal already emitted keeps its emission-time parameters. A command changes future emissions only.
- **Log records.** Start, stop, accepted, rejected (with reason), duplicate (returns the original outcome) and pending-at-stop each have an explicit record. Every record holds: requested time, effective time, ingress sequence, client/session id, request id and outcome. Stop may drain committed output and reports the final output boundary per channel.
- **Replay.** Replay applies logged effective times and ingress order directly. It never recomputes them from replay speed or wall clock.

## 5. Tap points: responses and state ownership

Every receiver channel has two taps.

| | Antenna-port tap | Virtual-receiver tap |
| --- | --- | --- |
| Purpose | Drive real hardware (HIL TX) | Feed software consumers |
| Pre-sampling response | B_out(f), output band selection for the DAC/output channel | H_rx(f), analog front end + anti-alias |
| Timebase | output-device timebase (ideal) | receiver ClockMap (§3.1) |
| Receiver effects | none | chain 1–7: gain/NF, kTBF noise, LO offset + phase noise, IQ imbalance + DC, decimation, AGC, ADC quantization/clip |
| Scale | its own `V_fs` (calibrated output) | virtual ADC full scale |

- **Response identity.** `TapResponseId` = content hash of (tap kind, response definition, output rate, center frequency). Work lists, edge tables and narrowband response samples are keyed by it. Changing any virtual-receiver setting cannot change an antenna-port key. This is the tap-isolation check.
- **Separate scenes.** Antenna-port and virtual-receiver work lists are separate whenever their response, clock map, rate, center or scale differ. Sharing is an optimization, allowed only when all of those are identical.
- **Response application.** Narrowband emitters: the response at the emitter offset folds into the amplitude and phase polynomials, and paths below the in-band threshold are culled. Analytic pulses: complex band-limited edge tables (edge ⊛ response, shifted to the offset), built deterministically before their effective boundary. Wideband emitters: filtered inside the resampler before the fractional delay.
- **Stateless kernels.** Render kernels see only the work list and immutable tables. They own no state across blocks.
- **Stateful processes.** Decimation filter history, AGC, phase noise, clock processes and resampler source history live in one owner per channel (`ReceiverChannelState`, or the waveform source on the control plane). They are handed off in block order: block b+1 starts from the state after block b. Parallel or multi-GPU rendering never resets or duplicates that state. Finalization of a channel happens on exactly one owner.
- **Noise.** Thermal noise is added once, after the complete path sum. Its Philox stream is keyed by (seed, receiver entity, `StreamId::receiver_thermal_noise(channel)`), with counter = sample index at the rate where noise is injected (§1.5).
- **Truth.** Truth and consumer-knowledge exports leave from the control plane. They never enter IQ samples, IQ headers or context packets.

### 5.1 Tap response family (settled, prompt 5)

`iqsim_core::TapResponse`. A linear-phase lowpass at complex baseband, f relative to the channel centre, with a Gaussian transition:

```text
|H(f)| = ½·[erf((f + f_c)/(√2σ)) − erf((f − f_c)/(√2σ))]      f_c = (f_p + f_s)/2
H(f)   = |H(f)|·exp(−j2π·f·τ_g)
h(t)   = 2f_c·sinc(2f_c(t − τ_g))·exp(−2π²σ²(t − τ_g)²)
```

- **Definition** (integers, so the identity hash is exact): passband edge f_p (µHz), stopband edge f_s (µHz), stopband attenuation A (milli-dB) and group delay τ_g (ps).
- **σ:** solved so that the passband droop at f_p, ½·erfc(z) + ½·erfc((f_p + f_c)/(√2σ)) with z = (f_s − f_c)/(√2σ), equals 10^(−A/20), by 200 fixed bisection steps on libm `erfc`. Both Gaussian tails count. The far one matters only for very narrow passbands (review fix: f_p = 1 Hz, f_s = 1 kHz, A = 20 dB had given |H(f_p)| = 0.80).
- **Stopband:** |H| ≤ 10^(−A/20) at and beyond f_s, decreasing monotonically. A is a minimum attenuation, not a flat floor.
- **Passband:** droop ≤ 10^(−A/20) at f_p, with no ripple.
- **Transition:** f_p < |f| < f_s, with |H(f_c)| = ½·erf(2f_c/(√2σ)), which is ½ to within the far tail.
- **Group delay:** τ_g at every frequency.
- **Validation:** `0 < f_p < f_s ≤ F_s/2` (anti-alias, compared exactly), 20 ≤ A ≤ 200 dB and 0 ≤ τ_g ≤ 1 ms.
- **Evaluation:** inside f_c, |H| uses the sum of two erf terms. Outside f_c it uses a difference of erfc terms, so the far stopband does not cancel to zero.
- **Truncation:** h is used on |t − τ_g| ≤ T_h = √(13/(2π²σ²)). Because |2f_c·sinc(2f_c t)| ≤ 1/(π|t|), the discarded tail has ∫|h| ≤ e⁻¹³/(13π) < 10⁻⁷. A unit test checks the numerical Fourier transform of the truncated h against H(f) to 2·10⁻⁷.
- **Why this family:** a raised-cosine transition's impulse response decays only as 1/t³. For a 1 MHz transition and 10⁻⁵ accuracy it needs ≈ 100 µs of support, while the Gaussian transition needs ≈ 4 µs.
- **Identity:** `TapResponseId` = SHA-256 of `"iqsim.tap-response.gauss-lowpass.v1\0" ‖ tap kind (u32) ‖ f_p ‖ f_s (i64) ‖ A (u32) ‖ τ_g (i64) ‖ rate num ‖ den (u64) ‖ centre µHz (i64)`, all little-endian. The same definition on the two taps has different identities.
- **Narrowband fold:** a CW path's amplitude is multiplied by |H(f_off)|, and −f_off·τ_g cycles are added to its phase D0.

### 5.2 Work list and CPU kernel (settled, prompt 5)

`iqsim_worklist::BlockWorkList` is one render target and one block. Element types are `#[repr(C)]` with explicit padding, contain no pointers, and use `u32` indices only. Sizes are pinned by a test:

| Type | Bytes | Fields |
| --- | --- | --- |
| `PolyRegisters` | 64 | D0..D3 as `[lo, hi]` u64 pairs |
| `PathSegment` | 240 | delay, coord and phase registers; block u64; `AmplitudePoly` (4 × f32 + `inv_len`); target, coord units and coord source u32; `EmissionRef` (kind, index) |
| `Pulse` | 32 | start and fall_start (Q32.32 target samples, relative to the emission time); amplitude f32; rise and fall table u32; pad |
| `PulseTrain` | 8 | first pulse, count (sorted by start) |
| `CwEmission` | 8 | constant complex envelope |
| `EdgeTableHeader` | 40 | u_lo (Q32.32); first, len, log2_res; plateau (2 × f32); peak f32; pads |

- **Waveform coordinate** `w(k) = F_s·(t_acq(k) − τ(k) − t_ref)`, Q32.32 target samples. `t_acq` comes from the ClockMap, so clock error moves pulse timing. `t_ref` is the origin of the source segment (for a pulse train, its emission time).
- **Rebasing:** |w| and |pulse offsets| < 2³⁰ samples (`MAX_COORD_SAMPLES`). A long-running source gets a later origin, so `w − offset` never overflows i64.
- **Coordinate shape (v1):** `validate()` accepts only linear, non-decreasing coordinates (D2 = D3 = 0, D1 ≥ 0), so the endpoint range check bounds every interior value. Endpoint checks alone would miss an interior wrap of a curved polynomial. Phase 3 must add an extremum and monotonicity check before it accepts cubic coordinates.
- **Kernel:** each worker takes a contiguous chunk, jumps the coordinate and phase polynomials to the chunk start, then steps. Per sample, in this order with separate f32 operations:
  1. `a = c0 + x(c1 + x(c2 + x·c3))`, x = f32(k)·inv_len, from scratch
  2. the envelope: CW constant, or `Σ A_i·(G_rise(w − start_i) − G_fall(w − fall_i))`
  3. `sincos(phase)`
  4. the complex product
  5. `m = a·lsb_per_volt`
  6. `I, Q = (pr·m, pi·m)`, ties-even to i64, checked add
- **Determinism:** every contribution depends only on k and w, never on the chunk, so output is independent of thread count, chunk size and path order. Errors are reported for the lowest failing chunk. Output is identical across 1, 2, 8 threads and chunk lengths 1 and 777 (golden `render/static_scene`).
- **Table lookup:** `d = sat_sub(u, u_lo)`. Below 0 the value is 0. At or past the last point it is the plateau. Otherwise linear interpolation, with the weight taken from the top ≤ 24 fraction bits, which is exact in f32.

### 5.3 Pulse edge tables (settled, prompt 5)

`iqsim_waveform::edge`. For a linear ramp edge e of duration T_r (0 = step), response h and carrier offset f_off:

```text
G(τ) = ∫ e(τ − x)·h(x)·e^{−j2π f_off x} dx
```

The renderer multiplies G by the path's carrier rotation, which equals filtering the RF pulse by H.

- **Construction:** f64 and libm in a fixed order.
  - 8-point Gauss–Legendre on equal panels per grid cell gives ∫g and ∫(x − τ_g)g, followed by prefix sums. Each cell has enough panels that none spans more than half a cycle of |f_off| + 2f_c, so a large carrier offset is never undersampled.
  - `G_i = S(τ_i − T_r) + (1/T_r)·∫_{τ_i−T_r}^{τ_i}(τ_i − x)g dx`, from prefix-sum differences plus one partial Gauss–Legendre cell.
- **Grid:** 1/R target samples with R = 2^log2_res. Every resolution shares the origin `u_lo = ⌊(τ_g − T_h)·F_s·64⌋/64`.
- **Plateau:** the analytic H(f_off), rounded to f32. It is stored as the last point of every table with the same (response, offset), so a finished pulse cancels exactly. The quadrature's own end value must match it to 10⁻⁶ (`PLATEAU_ERROR_BOUND`), which also checks that the quadrature converged; otherwise R is doubled and, at the limit, the build fails (`EdgeError::Quadrature`). (Review fix: taking the plateau from an undersampled base grid gave 0.355 instead of ≈ 0 for a 512 MHz offset at 1 MS/s.)
- **Resolution:** start at R = 64. Build at 2R and measure the linear-interpolation error at every cell midpoint. Double R until the error is ≤ 5·10⁻⁶ of the unit edge height, up to R = 2¹², and fail otherwise.
- **Error budget:** truncation < 10⁻⁷, interpolation ≤ 5·10⁻⁶ (measured, not proven between midpoints) and f32 storage < 6·10⁻⁸. The total is below 10⁻⁵·V_ref, which is §2's below-floor bound.
- **Measured:** the 16 MS/s reference responses need R = 128 for both in-band and −60 MHz offsets (16 691 points). A pulse matches an untabulated Simpson reference to 3.4·10⁻⁶·V in band and 4.3·10⁻⁶·V at −60 MHz.
- **Key:** `EdgeTableKey = (TapResponseId, EdgeShape, offset µHz)`. The output rate and tap are part of the response id.
- **Cache preparation:** tables are immutable. A runtime frequency or response change creates a new key, whose table is built while the work list for the change's effective block is constructed, before that block can render. Blocks already built keep their tables. Eviction is the runtime's job (phase 5).

## 6. Run manifest and replay

A versioned document written alongside every run. It is never part of consumer-visible output.

**Fields (initial).**

- `manifest_schema_version`, `numeric_contract_version`.
- **build:** package version, git commit and dirty flag, rustc version, target triple, profile, enabled features, `Cargo.lock` SHA-256, CUDA toolkit version and fatbin SHA-256s when the GPU is used.
- **scenario:** the resolved configuration in canonical serialization with its SHA-256, plus the source file paths and their SHA-256s. Canonical form: *proposed* RFC 8785 JSON canonicalization, decided when the manifest is implemented (prompt 9 or 10).
- `seed`, run mode (offline or real time), block length and pipeline depth.
- **backend and partitioning:** recorded for diagnostics. By contract they cannot change sample or truth hashes.
- **event log:** path, SHA-256, record count, first and last ingress sequence, final C.
- **assets:** one entry per external input: `{kind, logical name, resolved path, byte length, SHA-256, format/version}`. Kinds: library item, terrain tile, geoid grid, antenna pattern, IQ input, calibration table, live recording, live timing/gap record. Live recordings and their timing/gap records become assets when finalized.
- **outputs:** per stream: stream id, tap, sample format, first and last sample index, sample count, SHA-256 of the sample bytes, finalization status. Also the canonical truth-record SHA-256.
- **nondeterministic section** (excluded from every deterministic hash): wall-clock start and end, host, telemetry summary, performance counters, container-file metadata.

**Replay rules.** `sim replay <manifest>` verifies the build, numeric contract version, resolved scenario, event log and every asset hash before rendering anything. Missing or changed content is rejected. Different builds are not replay. They are compared against stored golden outputs with stated tolerances. Hash algorithm everywhere: SHA-256.

## 7. Transport choice, framing and delivery guarantees

Handlers and framing are independent of transport adapters (§9). The adapters are thin, so an alternative transport can be added later without touching handlers.

**Decision for v1 and alternatives** (from the design's transport table):

| Use | v1 choice | Guarantee | Alternatives kept open |
| --- | --- | --- | --- |
| Commands and status | ZMQ ROUTER (sim) / DEALER (clients) | acknowledged, deduplicated, logged | HTTP/WebSocket (phase 12, same handlers); gRPC if a consumer requires it |
| Operator telemetry | ZMQ PUB | best effort, sequence-numbered, snapshot resync | — |
| Strict software IQ | framed ZMQ ROUTER/DEALER over TCP or IPC | end-to-end acknowledged or the run fails | plain framed TCP, shared memory, DDS: revisit for a concrete integration or measured bottleneck |
| Receiver-interop IQ | VITA-49.2 over UDP | paced production; the consumer detects discontinuities | — |
| Preview | raw ZMQ PUB, sc16/fc32 (GNU Radio ZMQ SUB compatible) | best effort, no header, losses uncounted | — |
| Recording | SigMF files | checked writes, finalized count and hash | — |

This is a design choice to validate with the real consumer and target host, not a throughput claim.

**Control envelope.** JSON: `{schema_version, type, client_id, session_id, request_id, payload}`. The JSON Schema is generated from the Rust types and published.

- The ROUTER uses `ROUTER_MANDATORY`, bounded queues (HWM) and send deadlines. Delivery failure is reported.
- Ack = after durable log append. It returns `t_eff`, the ingress sequence and the outcome.
- Retries are deduplicated by `(client_id/session_id, request_id)`. A retry returns the original outcome and never schedules the action again.
- Clients can query command and run status.
- Telemetry carries sequence numbers. A status/snapshot request resynchronizes a client.

**Strict framed IQ.**

- Session per configured receiver-channel stream, with a dedicated queue per consumer. Fan-out is explicit.
- Handshake: consumer READY.
- Binary message header (versioned): stream/session id, block sequence, first sample index, sample count, datatype and byte order, receiver-reported first-sample timestamp (§3.2), context version. Context carries RF reference frequency, sample rate, bandwidth and gain. No truth fields and no manifest.
- Consumers send cumulative acknowledgements of contiguous blocks. The window is bounded in outstanding bytes and blocks. End-of-stream is an acknowledged final boundary.
- Queueing a message is not delivery.
- Real time: a disconnect, a missing ack, or backpressure past the configured deadline fails the sink and the run. Reconnect never silently resumes a stream as successful.
- Offline: the sender may wait up to its configured operation timeout.
- A reference receiver and a GNU Radio adapter are supplied.

**VITA-49 over UDP.** The interoperability profile is pinned when implemented (phase 5): VITA 49.2 signal data packets with stream id, 4-bit packet count, integer and fractional timestamps (units pinned with the profile, *proposed* UTC seconds + picoseconds), and payload sample count. Context packets carry RF reference frequency, sample rate, gain and bandwidth at a pinned cadence and on every change. MTU is pinned (*proposed* 1500 bytes, jumbo optional). UDP send success proves nothing about delivery. The supplied receiver checks continuity from packet count, timestamp and payload count together, handling count wrap. It reports loss, reordering and duplicates, and fails strict consumption on any gap.

**SigMF.** Sequential writes with every error checked. The sample count and SHA-256 are finalized at close. Truth annotations go to a separate sidecar file, never into the consumer's `.sigmf-meta`. In real time, recording shares the sink deadline and resource budget. Offline recording is unpaced.

**Failure classes** (kept distinct in status, logs and exit codes):

1. Producer underrun: hard error.
2. Strict-sink deadline or delivery failure: hard error.
3. Preview loss: allowed and never blocks required sinks.
4. UDP loss detected by the consumer: reported by the consumer.

## 8. Resource budgeting

- Every sink and device declares resource identities (NIC, disk, USB controller, PCIe link/root) with usable capacities.
- `sim lint` and `sim dry-run` sum concurrent traffic along every shared resource: payload, packet/framing/ack overhead, fan-out copies, and i64 emitter-partition transfers (16 B per complex sample per partition).
- Each resource gets a configured headroom allowance.
- Rough payload limits: sc16 about 30 MS/s per 1 GbE and 300 MS/s per 10 GbE. fc32 halves both.
- Two 20 MS/s sc16 sinks on one 1 GbE link need 1.28 Gb/s before overhead, so they fail the real-time check even though each fits alone.
- Unknown capacity is reported as *unverified*, never as passing.

## 9. Crate placement and dependency graph

Shared contracts sit low in the graph, so no runtime ↔ control cycle can form:

```text
iqsim-core ─┬─ iqsim-units ─┬─ iqsim-contracts ─┬─ iqsim-control      (adapters → handler traits)
            │               │                   ├─ iqsim-io ── iqsim-hil
            │               │                   ├─ iqsim-receiver
            │               │                   ├─ iqsim-scenario ── iqsim-behavior
            │               │                   └─ iqsim-runtime     (implements handler traits)
            │               ├─ iqsim-geo ── iqsim-antenna, iqsim-terrain
            │               └─ iqsim-waveform
            └─ iqsim-worklist ── iqsim-render ── iqsim-render-cuda
iqsim-cli = runtime + control + hil + scenario + contracts (wiring only)
```

- **iqsim-core:** TruthTime, SampleIndex, SampleRate, `ClockMap` trait and ideal map, commit-boundary time types, fixed point, sin/cos, Philox, the amplitude scale and headroom proof (§2), tap responses and their identities (§5.1), `NUMERIC_CONTRACT_VERSION`. `libm` is a normal dependency since prompt 5, for the control-plane response math only.
- **iqsim-worklist:** the `#[repr(C)]` work-list and edge-table arena types (§5.2). **iqsim-waveform:** the edge-table builder (§5.3). **iqsim-render:** the `Backend` trait and the CPU kernel.
- **iqsim-contracts:** pure data and pure encode/decode, with no sockets, files or simulation logic. It holds:
  - the control envelope, commands, acks and status
  - the handler traits (`CommandHandler`, `StatusProvider`)
  - telemetry records
  - the IQ frame header and context, and the VITA-49 profile constants
  - the run-manifest schema
- **iqsim-control:** ZMQ adapters that decode envelopes and call handler traits.
- **iqsim-io:** sinks that use contracts framing.
- **iqsim-runtime:** implements the handlers and owns the commit actor, pipeline and replay. It never depends on control or HIL.
- **iqsim-cli:** wires runtime to adapters.
- `cargo xtask deps` enforces the forbidden edges.

## 10. Proposed tolerances and decisions phase 1 must settle

**Proposed tolerances**

| Item | Proposed | Settled by |
| --- | --- | --- |
| Carrier phase error (above null floor) | ≤ 0.01 rad (**settled**). Time representation contributes < 10⁻²⁰ rad, polynomial quantization < 6.9·10⁻¹⁹ rad (§1.3), so the full budget goes to fitting. | prompts 1–2 |
| Delay error | ≤ 0.001 samples (**settled**). Time representation contributes < 2⁻⁶⁴ ps, polynomial quantization < 4.7·10⁻¹⁰ samples (§1.3). | prompts 1–2 |
| Acquisition-instant representation error | < 2⁻⁶⁴ ps, floor (**settled**) | prompt 1 |
| Amplitude error (above null floor) | ≤ 0.01 dB (**settled**) | prompt 2 |
| Null floor | V_ref · 10⁻² (−40 dB) (**settled**; V_ref defined in §2) | prompt 2 |
| Absolute complex-gain error below floor | ≤ 10⁻⁵ · V_ref (**settled**; continuous with the relative bound at the floor) | prompt 2 |
| sincos max abs error | ≤ 3·10⁻⁷ (**settled**; measured 5.98·10⁻⁸) | prompt 2 |
| Power check | ±0.1 dB vs closed-form link budget | phase 2 |
| Noise floor check | ±0.2 dB vs kTBF·NF | phase 2 |
| Diffraction check | ±0.1 dB vs independent P.526 | phase 8 |
| Out-of-band in-band level | ±0.1 dB vs selected H_rx/B_out | phase 2 |
| Determinism | identical hashes (CPU, GPU models, split/unsplit) | always |

**Decisions phase 1 must settle** (record the outcome here, then bump to contract version 1):

1. **Settled (prompt 1).** `SampleRate` bounds: den ≤ 2³², num ≤ 2⁴⁸, rate ≤ 1 GS/s, enforced at construction after reduction. The u128 no-overflow proof is in §1.2.
2. **Settled (prompt 1).** `AcqInstant` is a fixed 2⁻⁶⁴ ps fraction (floor). `ChannelTimebase::index_to_time` also returns the exact rational remainder. The inverse map returns an index only, and the caller gets the offset by exact Q64.64 subtraction (§3.1).
3. **Settled (prompt 2).** `MAX_BLOCK_LEN = 4 801 281`. Quantization error < 2 output LSB: < 4.7·10⁻¹⁰ samples of delay and < 6.9·10⁻¹⁹ rad of phase (§1.3). The true delay keeps a 2 LSB margin from the Q32.32 range ends.
4. **Settled (prompt 2).** sin/cos uses a 1024-entry f32 table, a 24-bit remainder, sin δ to the cubic term and cos δ to the quadratic, and an integer Taylor table builder (§1.5).
5. **Settled (prompt 2).** Subnormals are preserved on both backends (CPU default, CUDA `-ftz=false`) (§1.4).
6. **Settled (prompt 2).** `OffsetFrequency` is an i64 in 2⁻³² Hz, with exact floor-mod phase from Q64.64 time (§1.6).
7. **Settled (prompt 2).** R₀ = 50 Ω RMS, S = 40, headroom limit L = 2⁶², null floor V_ref·10⁻² with V_ref defined per link and block, and below-floor absolute bound 10⁻⁵·V_ref (§2).
8. **Settled (prompt 3).** Philox4x32-10 with the SplitMix64-finalizer key mix (injective within a seed), counter (sample index, sub-index), the `StreamKind` registry, the in-house f32 `ln` and the unit-power Box-Muller noise (§1.5).
9. Event-time grid for rounding C, and command vs entity-event ordering at equal times (before phase 4–5).
10. Error-budget sampling and refinement policy, including derivative bounds vs adaptive refinement and maximum depth (before phase 3).

Later phases pin these with the implementation: the manifest canonical encoding, the VITA-49 profile (timestamp units, MTU, context cadence) and the strict-sink window defaults.

## 11. Scenario quantities and files

Scenario files are TOML, read by `iqsim-scenario` into types from `iqsim-units`. `schema/scenario.schema.json` is generated from those types (`cargo xtask schema`), and a test fails when it is stale.

### 11.1 Quantity text

- A quantity is a string: a number, optional whitespace, then a **case-sensitive** unit (`"1090 MHz"`, `"1090MHz"`; `"1090 mhz"` is an error because `mHz` and `MHz` differ).
- The number is `[+-]digits[.digits][e[+-]digits]` with at most 38 significant digits. It is parsed as an **exact decimal**, never through an intermediate f64.
- A bare TOML number where a unit is required is rejected, and the message suggests a spelling. Latitude and longitude are the exception: a bare number is decimal degrees.

### 11.2 Representations

| Type | Units | Representation | Rule |
| --- | --- | --- | --- |
| `Frequency` | Hz, kHz, MHz, GHz | i64 µHz (±9.2 THz) | exact; finer than 1 µHz is an error |
| `SampleRate` | S/s, kS/s, MS/s, GS/s | `iqsim_core::SampleRate` | exact reduced rational; core bounds (§10 decision 1) |
| `Duration` | s, ms, us/µs, ns, ps | `TruthDuration` (i64 ps) | exact; finer than 1 ps is an error |
| `DataRate` | b/s, kb/s, Mb/s, Gb/s, B/s, kB/s, MB/s, GB/s | u64 bit/s | exact; SI prefixes (×1000); B = 8 b |
| `Power` | W, mW, dBm, dBW | f64 W, ≥ 0 | dB via `libm::exp10` |
| `Gain` | dB, dBi | f64 dB | both stored as dB |
| `Distance` | m, km, ft, nmi | f64 m | ft = 0.3048 m, nmi = 1852 m |
| `Speed` | m/s, km/h, kt | f64 m/s | kt = 1852/3600 m/s |
| `Angle` | deg, °, rad, DMS | f64 deg | |
| `Altitude` | m, km, ft, nmi + reference, or `FLnnn` | f64 m + `AltitudeRef` | reference required |

- `Display` output parses back to the same value for every value parsed from text. A `SampleRate` prints as an exact decimal (`2.5 S/s`). Only a rate built in code with a non-decimal denominator prints as `num/den S/s`. A flight level prints as `FLnnn`.
- Float types apply decimal unit factors exactly (km, mW, ft, nmi) and round once to f64 with the correctly rounded std parser. km/h and kt then divide once by 3600. The f64 bits are therefore the same on every platform.
- **Altitude references:** `MSL`, `AGL`, `ellipsoid` (alias `HAE`). `FLnnn` is a pressure altitude in the ISA standard atmosphere (1013.25 hPa datum) and stores `AltitudeRef::PressureStd` with nnn × 30.48 m. It is not a geometric height. Converting it to geometry needs the ISA model and any local pressure offset (plan, "Mode C altitude"), so it accepts no other reference.
- **DMS:** `D°M'S"` (also `′ ″`, `º`) or ASCII `DdMmSs`. Minutes and seconds are below 60, only the last component may have a fraction, and only the degrees may carry a sign. Latitude and longitude also take an **upper-case** hemisphere letter as a prefix or suffix (N/S, E/W), which cannot be combined with a sign. Lower case is excluded because `s` is the seconds marker.

### 11.3 Positions and MGRS

- A position is a table `{ lat, lon, alt }` or `{ mgrs, alt }` on WGS84. Giving both forms, or only one of lat/lon, is an error.
- An MGRS reference names a grid square. **A scenario position is the square's south-west corner** (NGA convention). `Mgrs::center()` is also available, and it is GeographicLib's representative point.
- The decoder uses WGS84 "AA" lettering. Column sets follow the zone mod 3, and row letters start at F in even zones. The 2000 km northing cycle is chosen by requiring the **100 km block** (not the digit square, so validity does not depend on precision, as in GeographicLib) to meet its latitude band, with 0.01° slack. A southern-hemisphere block must lie south of the equator. Zones 32X, 34X and 36X are rejected. Whitespace is ignored and letters may be lower case.
- The UTM inverse uses Karney's (2011) Krüger series to n⁶ and Newton iteration, all through `libm`. It agrees with GeographicLib GeoConvert to 5·10⁻¹² ° and with PROJ to 5·10⁻¹³ °, which is the print precision of the reference files.
- Polar UPS references (bands A, B, Y, Z) are rejected as not supported.

### 11.4 Scenario file rules

- Every table uses `deny_unknown_fields`. The top level requires `schema_version = 1`, `metadata.name`, `epoch`, `seed` and `run_mode` (`offline` or `realtime`).
- **Epoch:** a TOML offset date-time (or the same text as a string) in UTC (`Z` or `+00:00`), with fractional seconds to 1 ns. It is stored as POSIX seconds plus nanoseconds. A leap second (`:60`) is rejected.
- **Seed:** a TOML integer from 0 to 2⁶³−1.
- **Names:** entity names (platforms, emitters, receivers, antenna ports, sinks, resources) share one namespace and match `[A-Za-z0-9][A-Za-z0-9_-]*`. Mount names are unique within their platform, and channel names within their receiver. Library references (antennas, behaviors, waveforms) may also contain `.` and `/`. They are not resolved yet.
- **Lint checks:** references name an entity of the right kind, and a mount reference names a mount on the referenced platform. Emitter and centre frequencies are positive. Each channel's H_rx bandwidth and each antenna port's B_out bandwidth is positive and at most its complex sample rate, compared exactly as µHz·den ≤ 10⁶·num. A receiver has at least one channel.
- **Diagnostics:** `file:line:col: error: message` with 1-based lines and columns counted in characters, then the source line with a caret under the span. Duplicates add a `note:` at the first definition. `sim lint` exits 0 when clean, 1 on errors or an unreadable file, and 2 on usage errors.
- Clock and LO models (`{ model = "ideal" }`), trajectories (`{ kind = "static" }`), mounts, sinks and resources are placeholders that later prompts extend. Resource capacity is optional, and a missing capacity is reported as unverified (§8).
