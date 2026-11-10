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
12. Frames, geodesy and altitude
13. Trajectories and attitude
14. Light time and per-block polynomial fitting

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

- **Entity ids (prompt 6).** A receiver's entity id is the first four bytes, little-endian, of SHA-256(`"iqsim.entity.v1\0"` ‖ name). It does not depend on file order, so reordering a scenario keeps every noise stream. Two receivers whose ids collide are rejected at resolution. Other entity kinds get ids the same way when they first draw random numbers.
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

**Continuous sources (prompt 6a).** The envelope peak of a continuous path is `FARROW_LEBESGUE·max|z|` over its block's stage-1 window, plus a 10⁻⁵ allowance for the f32 Farrow evaluation (§5.6). It is not a declared crest factor.

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

**Output conversion (settled, prompt 6).** sc16 or cf32 per sink.

- **Antenna port** (`iqsim_receiver::PortConverter`), integers only: sc16 `code = round_ties_even(acc / 2²⁵)`, so code 32768 is 2⁴⁰ LSB, the **effective full scale** `2⁴⁰ / lsb_per_volt` (`TargetScale::full_scale_effective_volts`, within 2⁻²⁴ relative of the declared V_fs, and reported in the truth record). A code outside [−32768, 32767] is an error (`PortError::Clip`): clipping would corrupt the hardware drive. cf32 is `f32(acc)·2⁻⁴⁰` (volts / effective full scale), whose only rounding is the i64 → f32 conversion.
- **Virtual receiver:** the ADC stage of the chain (§5.4). Clipping there is a reported receiver effect, counted per stream in the manifest (`adc_clipped_components`).

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
- **Identity:** `TapResponseId` = SHA-256 of `"iqsim.tap-response.gauss-lowpass.v1\0" ‖ tap kind (u32) ‖ f_p ‖ f_s (i64) ‖ A (u32) ‖ τ_g (i64) ‖ rate num ‖ den (u64) ‖ centre µHz (i64)`, all little-endian. Tap kind codes: 1 antenna port, 2 virtual receiver, 3 source reconstruction (φ of a sampled source, §5.6; not a tap). The same definition on two kinds has different identities.
- **Narrowband fold:** a CW path's amplitude is multiplied by |H(f_off)|, and −f_off·τ_g cycles are added to its phase D0.

### 5.2 Work list and CPU kernel (settled, prompt 5)

`iqsim_worklist::BlockWorkList` is one render target and one block. Element types are `#[repr(C)]` with explicit padding, contain no pointers, and use `u32` indices only. Sizes are pinned by a test:

| Type | Bytes | Fields |
| --- | --- | --- |
| `PolyRegisters` | 64 | D0..D3 as `[lo, hi]` u64 pairs |
| `PathSegment` | 248 | delay, coord and phase registers; block u64; `AmplitudePoly` (4 × f32 + `inv_len`); target, coord units and coord source u32; `EmissionRef` (kind, index); `first`, `len` u32 (since prompt 9) |
| `Pulse` | 32 | start and fall_start (Q32.32 target samples, relative to the emission time); amplitude f32; rise and fall table u32; pad |
| `PulseTrain` | 8 | first pulse, count (sorted by start) |
| `CwEmission` | 8 | constant complex envelope |
| `EdgeTableHeader` | 40 | u_lo (Q32.32); first, len, log2_res; plateau (2 × f32); peak f32; pads |

- **Waveform coordinate** `w(k) = F_s·(t_acq(k) − τ(k) − t_ref)`, Q32.32 target samples. `t_acq` comes from the ClockMap, so clock error moves pulse timing. `t_ref` is the origin of the source segment (for a pulse train, its emission time).
- **Rebasing:** |w| and |pulse offsets| < 2³⁰ samples (`MAX_COORD_SAMPLES`). A long-running source gets a later origin, so `w − offset` never overflows i64.
- **Segments (prompt 9):** a `PathSegment` covers the block samples `[first, first + len)` (`len ≥ 1`, `first + len ≤` block length), and its delay, coordinate, phase and amplitude polynomials are in the segment-local index j = k − first over `len` samples (amplitude x = f32(j)·f32(1/len)). A path is one segment per knot-free, budget-meeting interval (§14); segments of one path never overlap. The work-list hash is `iqsim.worklist.v2` and covers `first` and `len`.
- **Coordinate shape (prompt 9):** cubic coordinates are accepted when they are non-decreasing and wrap-free over the segment (`coordinate_is_monotone`): the exact value `D0 + j·D1 + C(j,2)·D2 + C(j,3)·D3` fits i128 at both ends, and the step `Δ1(j) = D1 + j·D2 + C(j,2)·D3` is ≥ 0 at j = 0, at j = len − 2 and, when D3 > 0, at the two integers around the vertex `j* = ½ − D2/D3`, all in checked i128 (overflow rejects). A wrap-free, monotone value between in-range endpoints never leaves the range, so the endpoint checks bound the interior. (v1 accepted only linear coordinates; endpoint checks alone would miss an interior wrap of a curved polynomial.)
- **Kernel:** each worker takes a contiguous chunk and, per segment, the part of the chunk inside `[first, first + len)`; it jumps the coordinate and phase polynomials to where the chunk enters the segment (j = max(k0, first) − first), then steps. Per sample, in this order with separate f32 operations:
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

### 5.4 Virtual-receiver chain (settled, prompt 6)

`iqsim_receiver::chain`. It runs once per channel on the complete reduced i64 sum of the channel's virtual-receiver target. H_rx is already applied to every path before sampling, so the chain never filters the signal by H_rx again. Each stage is switched on by an `Option` in `ChainConfig`. Every stage is f32 with each operation separate and in source order; constants come from f64 + libm on the control plane, rounded once to f32.

| # | Stage | Definition |
| --- | --- | --- |
| 0 | volts | `x = f32(acc) / lsb_per_volt` at the receiver input |
| 1 | front end | `x ·= g`, `g = 10^(G/20)`. The noise figure F feeds stage 2. Off: G = 0 dB, F = 0 dB |
| 2 | thermal noise | `x += Σ_m t_m·w[n + m]`, w = unit-power Philox noise of `ReceiverThermalNoise(channel)` at the channel sample index. `t_m = σ_w·T·h(τ_g + (m − c)T)`, m = 0..L, `L = 2⌈T_h/T⌉ + 1`, c = (L − 1)/2, `σ_w = √(N₀·F_s·R₀)`, `N₀ = k·T₀·F·g²` (k = 1.380649·10⁻²³ J/K, T₀ = 290 K). Output power is N₀·B_n with B_n = ∫\|H_rx\|² df |
| 3 | LO offset | `x ·= e^{−j2π·Δf·t_acq(n)}`: `sincos(−phase_at(Δf, ClockMap.instant(n)))`, exact modular phase (§1.6). Phase noise: stub (`PhaseNoise::None`) |
| 4 | IQ, DC | `I' = I + dc_i`, `Q' = a·(Q·cos φ − I·sin φ) + dc_q`, a = 10^(ε/20), \|φ\| < 45° |
| 5 | decimation | FIR from the §5.1 family at the input rate (τ_g = 0, f_s ≤ output Nyquist), sampled as in stage 2 with unit scale. Output j = Σ_m t_m·x[jM − m] (causal, zero history at start), so it adds (L − 1)/2 input samples of delay. Output rate F_s/M |
| 6 | AGC | per output sample: `y = g·x`, `p += α(\|y\|² − p)`, `g += β·g·(T² − p)/T²`, clamp. Starts at g = 1, p = T². Off by default |
| 7 | ADC | `code = clamp(round_ties_even(x·2^(b−1)/V_adc))`, 2 ≤ b ≤ 16, left-justified (`code << (16 − b)`) in sc16; clips counted. cf32 output is code/32768, or with the ADC off the unquantized `x / V_adc` |

- **Full scales.** V_adc is the per-component clip voltage: a CW tone of power V_adc²/R₀ just reaches it (scenario `adc.full_scale`, a power). The virtual target's accumulator full scale is V_adc/g, so accumulator and ADC share headroom.
- **State ownership and handoff.** One `ReceiverChannelState` per channel owns the decimation history (L − 1 inputs) and phase, and the AGC gain and power. `process(block, first, acc, clock)` requires block b+1 after block b and contiguous sample indices, and rejects anything else (`OutOfOrder`, `Gap`) without changing state; it never resets. Noise carries no state: it is a correlation over counter-indexed white noise, so each block draws its own L − 1 sample lookahead from the counter. The LO phase comes from each sample's acquisition instant. The output is therefore identical for any block split (tested with splits down to 1 sample and with every stage on).
- **Determinism note for the CUDA port.** The chain runs on the host after the global reduction in v1. A device port must keep the f32 operation order above, including the noise and FIR sums in increasing m.

### 5.5 Antenna-port output (settled, prompt 6)

The antenna-port target renders through B_out with its own V_fs (scenario `full_scale`, a power: code 32768 ↔ a CW tone of that power) and the ideal output-device timebase (epoch 0, nominal rate). Its accumulator goes straight to `PortConverter` (§2): no gain, noise, LO, IQ, decimation, AGC or ADC, and no state. A test changes H_rx, gain, NF, noise, seed, clock, LO, ADC, IQ, decimation and AGC one at a time: the antenna-port hash never changes, and the virtual-receiver hash always does.

### 5.6 Continuous sources and the two-stage resampler (settled, prompt 6a)

`iqsim_waveform::continuous` (stage 1), `iqsim_render::farrow_value` (stage 2), `iqsim_runtime::scene` (coordinates and windows; `static_scene` until prompt 9), `iqsim_io::iq_input` (SigMF IQ files).

**Model.** A source is N complex samples x[n] at native rate F_n, normalized so |x| = 1 is full scale (ci16: i/32768; cf32 as stored). Sample n is emitted at `start + n/F_n` on the emitter carrier f_e. The continuous envelope is `x(t) = Σ x[n]·φ(t − n)` (t in native samples), with φ the **source reconstruction response**: the §5.1 family with tap kind `SourceReconstruction` (code 3), zero group delay, default passband 0.8·F_n and stopband F_n (full widths), 80 dB. Its stopband edge ≤ F_n/2 puts every image of the native spectrum at least A dB down. For a target with response H and Δf = f_e − f_c, a path delivers `a·z(w)·e^{jφ(k)}` with

```text
z(w) = Σ_n x[n]·K(w − n),   K(t) = ∫ Φ(ν·F_n)·H(ν·F_n + Δf)·e^{j2πνt} dν
w(k) = F_n·(t_acq(k) − τ − start)          (retarded native coordinate, ClockMap instants)
```

K is φ convolved with h shifted to the offset, exactly as an edge table is the edge convolved with it (§5.3), and it includes H's group delay. The carrier rotation `f_off·t_acq − f_e·τ` is the path's ordinary phase polynomial, so **frequency translation happens at the target rate**, never on native samples: a low native rate can sit anywhere in a wider channel without aliasing.

**Kernel construction** (`ResampleKernel::build`, f64 + libm, fixed order). The frequency integral runs over the band where both responses exceed erfc(8)/2, with 8-point Gauss–Legendre on equal panels (half a cycle of e^{j2πνt} at the largest |t| and at least 8 per σ). Panels double until no tap moves by more than 10⁻⁹ of the largest (`KERNEL_QUADRATURE_BOUND`), up to six doublings, else `Quadrature`. K is truncated to supp φ + supp h, which costs at most ~10⁻⁷·(‖φ‖₁ + ‖h‖₁). Key: `ResampleKey = (φ id, H id, offset µHz)`; ids cover rates, tap kind and centre.

**Stage 1** (`filter_window`, control plane, per source × target × block). z is evaluated on a grid of U = 2^log2_up points per native sample: `z[m] = Σ_i x[⌊m/U⌋ − i]·K((m mod U)/U + i)`, summed in f64 over increasing i (complex products as `xr·kr − xi·ki`, `xr·ki + xi·kr`), each point multiplied by the source's level factor (1, or 1/RMS) and rounded once to f32, so z is near unit scale whatever the recording's raw level. U is the smallest power of two ≤ 32 with `lagrange6_error(ν_max) ≤ 2·10⁻⁶` (`FARROW_ERROR_BOUND`), `ν_max = min(f_s,φ, |Δf| + f_s,H)/(U·F_n)`; `lagrange6_error` is the worst |interpolated − 1| of a unit tone over 513 values of μ. U = 16 suffices for any content up to native Nyquist (2.8·10⁻⁷).
- **History and halos.** Each z[m] is a pure function of m and the source, so windows of different block splits are bit-identical where they overlap. No filter state is carried: each window reads its own halo (`source_span`) from the immutable source.
- **Boundary padding** (`SourceTimeline`): zeros before sample 0 for every policy (the onset is band-limited by K). After sample N − 1: zeros (`stop`, `zero_pad`) or `x[n mod N]` (`loop`, the seam is sample N − 1 followed directly by sample 0, no gap or crossfade, band-limited by K). `stop` is strict: a block whose retarded coordinate passes N − 1 (`w_max > (N − 1)·2³²`) over the samples the run keeps fails with `SourceExhausted` (the builder's `with_sample_limit` excludes the discarded tail of a final partial block); its filter halo past the end still reads zeros.
- **Group delay.** The target's τ_g only (φ is zero-phase); stage 2 adds none.
- **Lookahead.** A z point reads `i ∈ tap_range()` native samples around it (supp φ + supp h), plus 3 cells of Farrow lookahead. For a recorded file this is known in advance. A live or generated source must count it as filter lookahead for the commit boundary (§4).

**Stage 2** (render kernel, f32). Work-list emission kind `EMISSION_CONTINUOUS` (3) points to a `ContinuousEmission { window }`; windows live in the block's `SourceWindowArena` (`SourceWindowHeader { first, len, log2_up, peak }`, 24 bytes, ≤ 2²⁶ points, log2_up ≤ 5). The path coordinate uses `COORD_NATIVE_SAMPLES`: Q32.32 native samples relative to the window's z[0]. Per sample, every operation separate, per component:
1. `u = w << log2_up`, `m = u >> 32`, `μ = f32((u mod 2³²) >> 8)·2⁻²⁴` (exact).
2. For j = 1..5: `c_j = FARROW[j][0]·z[m−2] + FARROW[j][1]·z[m−1] + … + FARROW[j][5]·z[m+3]`, accumulated left to right from the first product.
3. `y = c5`; then `y = c_j + μ·y` for j = 4, 3, 2, 1; then `out = z[m] + μ·y`.

`FARROW` is the power-basis matrix of the 6-point Lagrange interpolator on nodes −2..3, each exact rational rounded to f32 (`iqsim_render::FARROW`). The CUDA port (prompt 20) must use this table, this order, `--fmad=false`, and the same μ extraction. Stage 1 stays on the host in f64; a device port of it must keep the f64 sum order above.
- `validate()` requires native units for continuous paths (and target units for the others), that `D0 + (len − 1)·D1` of every linear coordinate does not leave i128 (no wrap mod 2¹²⁸, so the endpoints bound the interior), and, from the endpoints, `2 ≤ ⌊w·U⌋` and `⌊w·U⌋ + 4 ≤ len`, so every tap is inside the window.

**Coordinates and rebasing.** At k = 0 the builder computes `F_n·(t_acq(0) − start)` exactly from the Q64.64 instant and the rational rate in i128 (`q64_interval_to_q32`, within one Q32.32 unit), subtracts `F_n·τ` rounded once to Q32.32, and sets that as D0 of a slope-only polynomial. A group of links that share a source, kernel, start, timeline (length and end policy) and level factor gets one window covering all their coordinates plus the Farrow taps and one cell of margin. Each coordinate is then rebased to that window's z[0]. A window above 2²⁶ points is rejected (`WindowTooLong`) before anything is allocated or read. Coordinates stay below a block's worth of samples, so neither 100-day runs nor many loops can overflow.

**Level, headroom and culling.** Path amplitude = Friis volts; the normalization factor (1 for `full_scale`; 1/√(mean |x|²) for `rms`) is applied in stage 1. The envelope bound is `FARROW_LEBESGUE·max|z|·(1 + 10⁻⁵)` from the window itself (max|z| computed in f64 and rounded up to f32), with `FARROW_LEBESGUE` = 1.3911 ≥ max_μ Σ|ℓ_i(μ)| and 10⁻⁵ for the f32 evaluation; it enters the §2 proof like any envelope peak. Culling is lossless as in §2. A path is dropped before stage 1 when even `Λ·max_n Σ|K|·factor·max|x|` cannot reach ¼ LSB, or when its window's whole source span is padding zeros.

**Error budget** (relative to a component's amplitude): K truncation ~10⁻⁷, quadrature 10⁻⁹, f32 storage 6·10⁻⁸, Farrow ≤ 2·10⁻⁶ at ν_max. **Measured** (`crates/iqsim-runtime/tests/iq_file.rs`, both taps): tone level within 4·10⁻⁷ dB of Friis·|Φ(f_t)|·|H(Δf + f_t)|; frequency within 3·10⁻¹⁰ Hz; native-rate images ≤ 1.6·10⁻⁸; chirp against an independent direct sum with a time-domain Simpson kernel ≤ 2.2·10⁻⁷ of peak; envelope delay ≤ 7.8·10⁻⁷ samples; block splits of 16 000, 4000 and 1777 samples bit-identical.

**IQ files** (`iqsim_io::iq_input`). `.sigmf-meta`: `core:datatype` exactly `ci16_le` or `cf32_le`; `core:sample_rate` finite and > 0, made exact as the shortest round-trip decimal parsed by `iqsim_units` (a scenario `sample_rate` must equal it); `core:version` present; `core:num_channels` 1 if present; `core:sha512` checked if present; captures are objects with non-decreasing `core:sample_start` < N. `.sigmf-data`: one streaming pass in 1 MiB pages (non-empty, a whole number of samples, otherwise "truncated"; cf32 finite) computes SHA-256, per-page SHA-256, peak |x| and Σ|x|² (f64, file order). Rendering reads through an LRU cache of at most 64 pages, and each page is re-hashed when loaded, so a file changed mid-run is an error, never different samples.

## 6. Run manifest and replay

A versioned document written alongside every run. It is never part of consumer-visible output.

**Fields.** The initial list is below; what prompt 6 implements is under *Implemented*.

- `manifest_schema_version`, `numeric_contract_version`.
- **build:** package version, git commit and dirty flag, rustc version, target triple, profile, enabled features, `Cargo.lock` SHA-256, CUDA toolkit version and fatbin SHA-256s when the GPU is used.
- **scenario:** the resolved configuration in canonical serialization with its SHA-256, plus the source file paths and their SHA-256s. Canonical form: *proposed* RFC 8785 JSON canonicalization, decided when the manifest is implemented (prompt 9 or 10).
- `seed`, run mode (offline or real time), block length and pipeline depth.
- **backend and partitioning:** recorded for diagnostics. By contract they cannot change sample or truth hashes.
- **event log:** path, SHA-256, record count, first and last ingress sequence, final C.
- **assets:** one entry per external input: `{kind, logical name, resolved path, byte length, SHA-256, format/version}`. Kinds: library item, terrain tile, geoid grid, antenna pattern, IQ input, calibration table, live recording, live timing/gap record. Live recordings and their timing/gap records become assets when finalized.
- **outputs:** per stream: stream id, tap, sample format, first and last sample index, sample count, SHA-256 of the sample bytes, finalization status. Also the canonical truth-record SHA-256.
- **nondeterministic section** (excluded from every deterministic hash): wall-clock start and end, host, telemetry summary, performance counters, container-file metadata.

**Implemented (settled, prompt 6).** `iqsim_contracts::RunManifest`, schema version 1, unknown fields rejected:

- `numeric_contract_version` (`iqsim_core::NUMERIC_CONTRACT_VERSION`, still 0).
- `build`: package version, rustc release (from `rustc -vV` in `iqsim-runtime/build.rs`), target triple, profile, features (none yet), the SHA-256 of the workspace `Cargo.lock` compiled into the binary, and a **source fingerprint** (SHA-256 over the root `Cargo.toml`, `rust-toolchain.toml`, every crate's `Cargo.toml`, `build.rs` and `src/` files, by sorted relative path, with lengths). Any source edit, committed or not, is a different build. No git commit is recorded: the fingerprint covers what it would, without git at build time. Replay requires every field to match.
- `scenario`: absolute source path and file SHA-256, the resolved scenario (every default filled in, `iqsim_runtime::ResolvedScenario`) and its canonical SHA-256.
- `seed`, `run_mode`, `block_ps`, `duration_ps`, `assets` (`AssetRegistry`: kind, logical name, path, length, SHA-256, format, metadata SHA-256). Since prompt 6a, each IQ-file emitter registers `<emitter>.data` (kind `iq_input`, format `sigmf/<datatype>`, metadata = the canonical JSON of the datatype, exact rate, N, start, end policy, normalization and factor, reconstruction response and meta SHA-256) and `<emitter>.meta` (format `sigmf-meta`). Large files are hashed in a streaming pass (`register_digest`) and verified the same way (`verify_digests`), `event_log` (always absent until prompt 16), `outputs` (stream, sink, tap, format, output rate, centre, first index, count, SHA-256 of the sample bytes, path, finalized, ADC clip count) and `truth` (path and SHA-256 of the canonical truth record). The truth record (schema 1) has per stream the render and output rates and counts, the output → render index factor, the decimation delay in render samples, the first sample's acquisition instant (ClockMap) and reported label (timestamp clock), and the effective full scale; and per (stream, emitter) link the range, delay (s and render samples), carrier offset, response gain and power at the target plane.
- **Since prompt 7 (manifest schema 2, resolver version 2):** `warnings` (deterministic resolution warnings, e.g. the zero-geoid fallback; replay requires them to match), replay also requires the asset list to equal the one the re-resolved scenario produces (no omitted, extra or edited records), the resolved scenario's `world` (geoid model and file identity, pressure offset), and the `geoid` asset when a grid is configured (§12.5).
- **Since prompt 9 (resolver version 3, truth schema 2; manifest schema unchanged):** the resolved scenario has `platforms` (name, trajectory kind, span, knot count and `trajectory_sha256`, §14.5), and emitters and streams name their `platform` (`position_ecef_m` is the platform's position at truth time 0, or at its trajectory start if later). The built trajectories travel with the resolved scenario but are not serialized; replay re-resolves them from the verified scenario and requires the same fingerprints. The truth record's links are evaluated by light time at the stream's first sample (`range_m` = c·τ, `range_rate_m_s`, `doppler_hz`, `emission_time`, delay, response gain at the Doppler-shifted offset, power) and carry `geometry_fit` (blocks, segments, worst residuals and every split).
- `backend` (name, threads) and `nondeterministic` (wall-clock start/end, host, render seconds) are excluded from `deterministic_sha256`.
- **Canonical JSON (decided, prompt 6):** "iqsim canonical JSON v1": no whitespace, object keys sorted by UTF-8 bytes, serde_json string escaping, decimal integers and serde_json's shortest round-trip floats (`float_roundtrip` is enabled so parsing gives the same bits back). It is RFC 8785 except for float formatting outside 10⁻⁶..10²¹; exact quantities are integers, so they never reach that difference.
- **Run directory:** `consumer/<sink>/<stream>.sigmf-{data,meta}` for consumers; `run/manifest.json` and `run/truth.json` are privileged and never under `consumer/`.

**Replay rules.** `sim replay <manifest>` verifies the build, numeric contract version, resolved scenario, event log and every asset hash before rendering anything. Missing or changed content is rejected. Prompt 6 checks, in this order and before rendering: manifest schema, numeric contract version, build identity, no event log, scenario file hash, the recorded resolved scenario against its own hash, every asset, and the scenario re-resolved by this build against the recorded hash (plus seed, block and duration). It then re-renders into hashing sinks and requires every output record and the truth hash to match. Different builds are not replay. They are compared against stored golden outputs with stated tolerances. Hash algorithm everywhere: SHA-256.

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

**SigMF.** `iqsim_io::sigmf` (prompt 6). Data goes to `<stem>.sigmf-data` (never overwritten), and `<stem>.sigmf-meta` is written only by a successful `finish`, after flush, `sync_all` and a check that the file length equals count × sample size. The meta holds only `core:datatype`, `core:sample_rate`, `core:version` (1.2.0), `core:num_channels`, `core:recorder`, `core:description` and per capture `core:sample_start`, `core:frequency` and `core:datetime`; the API has no extension fields. `core:datetime` is the timestamp clock's label of the first sample's acquisition instant (§3.2; ideal until prompt 15), as UTC with 12 fractional digits, floored to the picosecond. Sequential writes with every error checked. The sample count and SHA-256 are finalized at close. Truth annotations go to a separate sidecar file, never into the consumer's `.sigmf-meta`. In real time, recording shares the sink deadline and resource budget. Offline recording is unpaced.

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
- **Continuous sources (prompt 6a):** the source trait, timelines and the stage-1 kernel in **iqsim-waveform** (`continuous`); window types in **iqsim-worklist**; the Farrow kernel in **iqsim-render**; the SigMF IQ reader in **iqsim-io** (`iq_input`); coordinates, grouping and the reader adapter in **iqsim-runtime**. No dependency edges were added: the runtime adapts the io reader to the waveform trait.
- **iqsim-geo (prompt 7):** vectors, geodetic ↔ ECEF, pseudo-inertial ECI, quaternions, NED/body/antenna frames, the PGM geoid reader, ISA, Gillham and the altitude chain (§12). It reads a geoid file but hashes nothing; the runtime hashes and registers it. No dependency edges were added (proptest is a dev-dependency).
- **Trajectories (prompt 8):** the `Trajectory` trait, maneuvers, routes, tracks (the smoothing spline), derived attitude and appends in **iqsim-geo** (`trajectory`, §13); the scenario → trajectory builder in **iqsim-runtime** (`trajectory`). No dependency edges were added.
- **Light time and fitting (prompt 9):** the light-time solver and the `Lookback` wrapper in **iqsim-geo** (`light_time`, `trajectory::Lookback`); the exact path model, the fitter, refinement and the fit log in **iqsim-runtime** (`geometry`), the work-list builder in **iqsim-runtime** (`scene`, which replaces `static_scene`); segment ranges in **iqsim-worklist**; the segment kernel in **iqsim-render**. No dependency edges were added.
- **iqsim-receiver:** the virtual-receiver chain and the antenna-port conversion (§5.4, §5.5). **iqsim-io:** SigMF and UTC labels (sha2, serde_json). **iqsim-runtime:** scenario resolution, the offline pipeline, truth record, manifest and replay (sha2, serde, serde_json, hex, libm; a `build.rs` that only runs `rustc -vV`).
- **iqsim-contracts:** pure data and pure encode/decode, with no sockets, files or simulation logic (serde, serde_json, sha2, hex since prompt 6). It holds:
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
| Power check | ±0.1 dB vs closed-form link budget (**tested**, prompt 6: worst 0.0004 dB over six tones and both taps) | phase 2 |
| Noise floor check | ±0.2 dB vs kTBF·NF (**tested**, prompt 6: total −0.037 dB, Welch passband −0.036 dB) | phase 2 |
| Diffraction check | ±0.1 dB vs independent P.526 | phase 8 |
| Out-of-band in-band level | ±0.1 dB vs selected H_rx/B_out (**tested**, prompt 6) | phase 2 |
| Doppler and phase | moving CW emitter's phase within the carrier-phase bound of the analytic solution (**tested**, prompt 9: 1.8·10⁻⁷ rad worst over six blocks through closest approach; frequency within 2.2·10⁻⁵ Hz) | phase 3 |
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
10. **Settled (prompt 9).** Error-budget sampling and refinement: cubic interpolation on Chebyshev–Lobatto nodes, 15 check points per interval against the exact model, acceptance at a quarter of every budget, conservative halving up to `MAX_REFINE_DEPTH` = 12 with an explicit failure beyond it, and splits first at trajectory knots mapped to receive time (§14.3). Point checks are diagnostics, not a proof between points; the factor 4 is the margin for f⁽⁴⁾ varying across an interval.

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
| `Speed` | m/s, km/h, kt, ft/min | f64 m/s | kt = 1852/3600 m/s; ft/min = 0.3048/60 m/s |
| `AngularRate` | deg/s, °/s, rad/s | f64 deg/s | |
| `Acceleration` | m/s², m/s^2, g | f64 m/s² | g = 9.806 65 m/s² |
| `Jerk` | m/s³, m/s^3 | f64 m/s³ | |
| `Angle` | deg, °, rad, DMS | f64 deg | |
| `Pressure` | Pa, hPa, mbar, inHg | f64 Pa, signed | inHg = 3386.389 Pa (conventional); mbar = hPa |
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
- **Lint checks:** references name an entity of the right kind, and a mount reference names a mount on the referenced platform. Emitter and centre frequencies are positive. Each channel's H_rx and each antenna port's B_out has 0 < bandwidth < stopband ≤ complex sample rate, compared exactly as µHz·den against 10⁶·num; the stopband defaults to the sample rate, so a bandwidth equal to the rate is an error (prompt 6; prompt 4 allowed it). Attenuation is 20–200 dB and group delay 0–1 ms. A `fixed_offset` clock or LO has an `offset` and an `ideal` one has none. Render duration and block are positive. A receiver has at least one channel.
- **Render and receiver fields (prompt 6).** `[render] duration, block` (default 1 ms; a whole number of samples at every rate, at most `MAX_BLOCK_LEN`). Per channel and antenna port: `response = { stopband, attenuation, group_delay }` (passband edge = bandwidth/2, stopband edge = stopband/2, floored to µHz; defaults rate, 80 dB, 0). Per channel: `front_end = { gain, noise_figure }`, `thermal_noise` (default true), `clock = { model = "fixed_offset", offset }` (an ideal map with a later epoch), `lo = { model = "fixed_offset", offset }`, `iq = { gain_imbalance, phase_skew, dc_i, dc_q }` (DC as a fraction of full scale), `decimation = { factor, bandwidth, stopband, attenuation }` (defaults 0.8 and 1.0 × output rate, 80 dB), `agc = { target, attack, step, min_gain, max_gain }` (defaults −12 dB FS, 0.01, 0.01, ±40 dB) and `adc = { full_scale, bits, enabled }` (defaults 0 dBm, 16, true). Per antenna port: `full_scale` (default 0 dBm). Per sink: `format = "ci16_le" | "cf32_le"`.
- **IQ-file emitters (prompt 6a).** `waveform = "iq_file"` with `iq_file = { path, start, end, normalization, sample_rate, response }`. `path` names the `.sigmf-meta` or `.sigmf-data` file, or their common stem, relative to the scenario file. `start` is the emission time of sample 0 (default 0, not negative). `end` is `stop` (strict: the recording must cover the run), `zero_pad` or `loop`, and it is required. `normalization` is `full_scale` (`power` is the power of a |x| = 1 sample) or `rms` (`power` is the recording's mean power), and it is also required. `sample_rate` is optional and must equal the recording's. `response = { bandwidth, stopband, attenuation }` is φ (defaults 0.8 × rate, 1 × rate, 80 dB). The emitter `frequency` is the RF frequency of the recording's 0 Hz. Lint requires an `iq_file` table exactly when `waveform = "iq_file"`, a non-empty path, a non-negative start, and a valid response (checked against `sample_rate` when given). The recording itself is checked at resolution.
- **Earth model (prompt 7).** `[world] geoid = { path }` (a GeographicLib PGM grid, relative to the scenario file) and `pressure_offset` (a `Pressure`, QNH − 1013.25 hPa, default 0). Lint requires a non-empty path and |offset| ≤ 200 hPa; the grid itself is read at resolution (§12.5, §12.6).
- **Trajectories (prompt 8).** `[platforms.trajectory] kind = "static" | "maneuvers" | "route" | "track" | "ground"` (default static; §13). Fields:
  - `start` (initial truth time, default 0 s), `heading` (initial track) and `speed`, for maneuvers, route and ground;
  - `attitude = { heading, pitch, roll }`, static only;
  - `angle_of_attack` and `wind = { speed, from }`, for the airborne kinds (maneuvers, route, track);
  - `limits = { bank, roll_rate, accel, jerk, vertical_accel, vertical_jerk, climb_rate }`, for moving kinds;
  - `[[…maneuvers]] { kind = straight | turn | climb | speed, duration, rate, bank, direction, heading, by, vertical_rate, angle, altitude, speed }`;
  - `[[…waypoints]] { lat, lon, alt, speed }`;
  - `points = [{ time, lat, lon, alt }]` and `smoothing` (τ), for tracks.

  Lint rejects:
  - a field its kind does not use;
  - a missing required field: heading and speed for maneuvers; speed for route and ground; at least one maneuver or waypoint; a ground trajectory needs exactly one of maneuvers and waypoints;
  - a `position` on a track (it starts at its first point), or a missing `position` on any other kind (`position` is optional in the schema since prompt 8);
  - a maneuver without exactly one rate (rate/bank or vertical_rate/angle) and one end (duration/heading/by or altitude/duration), or a turn by duration or angle without `direction`;
  - a non-positive limit, a bank outside (0°, 80°), or a speed below 1 m/s;
  - for ground trajectories, a non-AGL altitude or a climb;
  - fewer than two track points, or times that do not strictly increase.

  Geometric feasibility (a waypoint inside the turn circle, an unreachable altitude) is checked when the trajectory is built at resolution.
- **What `sim run` supports (prompt 6, updated in prompts 7, 8 and 9).** Static and moving platforms of every trajectory kind (since prompt 9): every trajectory is built at resolution, extended straight and level to the run end if its plan ends earlier, and rendered through light time (§14); MSL heights add the geoid height (0 without a grid), flight levels resolve through the inverse pressure-altitude chain, and AGL heights are still taken as ellipsoidal (no terrain yet); mounts and antennas are ignored (isotropic at the platform origin); emitters must use `behavior = "always_on"` and the built-in `waveform = "cw"` or `"iq_file"` (prompt 6a); only SigMF sinks are written, and a stream's sinks must agree on the format. Anything else is a resolution error, never a silent default.
- **Diagnostics:** `file:line:col: error: message` with 1-based lines and columns counted in characters, then the source line with a caret under the span. Duplicates add a `note:` at the first definition. `sim lint` exits 0 when clean, 1 on errors or an unreadable file, and 2 on usage errors.
- Clock and LO models (`{ model = "ideal" }`), mounts, sinks and resources are placeholders that later prompts extend. Resource capacity is optional, and a missing capacity is reported as unverified (§8).

## 12. Frames, geodesy and altitude (settled, prompt 7)

All in `iqsim-geo`, f64, every transcendental through `libm`, every product and sum a separate operation in a fixed order (no FMA), so the bits are the same on every machine. These are control-plane computations; the render path still sees only fixed-point polynomials.

### 12.1 Frames and Earth rotation

- **WGS84** a = 6 378 137 m, f = 1/298.257 223 563. **ECEF** is WGS84 Earth-fixed.
- **ECI is pseudo-inertial:** it coincides with ECEF at the scenario epoch (truth time 0), and ECEF spins about the shared z axis at ω_E = 7 292 115·10⁻¹¹ rad/s (WGS84). No precession, nutation, polar motion or UT1 (plan: below this fidelity). `r_eci = R_z(θ)·r_ecef`, `v_eci = R_z(θ)·(v_ecef + ω × r_ecef)`, θ = ω_E·t.
- **θ is reduced in integers** (`earth_rotation_angle`): for `t = ps + frac·2⁻⁶⁴` ps, `S = ps·7 292 115·2⁴¹ + ⌊frac·7 292 115/2²³⌋` (i128, units 10⁻²³·2⁻⁴¹ rad, |S| < 2¹²⁷ over all i64 ps), `r = S mod round(2π·10²³·2⁴¹)`, `F = ⌊r·2⁶⁴/T⌋` by exact long division (a 64-bit fraction of a turn), then `θ = F·2⁻⁶⁴·2π` with a double-double 2π and Dekker's two-product (no FMA). θ ∈ [0, 2π] is within 0.5 ulp + 3.5·10⁻¹⁹ rad of the exact angle at any time in the truth-time range (≤ 4.5·10⁻¹⁶ rad: 1.9·10⁻⁸ m at GEO radius, 1.7·10⁻⁷ m at lunar distance); a naive f64 ω·t has a 1.1·10⁻¹³ rad ulp at 100 days. Checked within 1 ulp (or 3.5·10⁻¹⁹ rad) of mpmath at 212 instants, including both ends of the i64 range. After rounding, θ can equal 2π just before a whole turn.

### 12.2 Geodetic ↔ ECEF

- **Forward:** the closed form `N = a/√(1 − e² sin²φ)`, `((N+h) cos φ cos λ, (N+h) cos φ sin λ, (N(1−e²)+h) sin φ)`. Its bits are unchanged since prompt 5.
- **Inverse:** Bowring's iteration on the parametric latitude, `β₀ = atan2(z, (1−f)p)`, `φ = atan2(z + e′²b sin³β, p − e²a cos³β)`, `β = atan2((1−f) sin φ, cos φ)`, stopping when |Δφ| ≤ 10⁻¹⁵ rad, capped at 8 steps (not converging is an error). Then `h = p cos φ + z sin φ − a√(1 − e² sin²φ)`. Domain 1000 km ≤ |r| ≤ 10¹² m (the evolute, where latitude is not unique, is within 43 km of the centre; the upper limit keeps `hypot` and h finite). On the polar axis φ = ±90°, λ = 0, h = |z| − b. λ = atan2(y, x) in (−180°, 180°]. Measured: at most 3 steps from −11 km to 4·10⁸ m.
- **Tolerances** (separate per dimension, tested against PROJ and a 60-digit mpmath reference, and by 2048-case round trips):
  - ECEF position and height: **10⁻⁶ m**, with a scale-aware allowance max(10⁻⁶ m, 4 ulp(|r|)) that only exceeds 1 µm beyond 2³⁰ m ≈ 1.07·10⁹ m. Measured worst: forward 7.9·10⁻⁸ m, inverse 1.4·10⁻⁷ m.
  - Latitude and longitude as arc lengths, |Δφ|·(N + h) and |Δλ|·p, within the same metre bound (≈ 1.6·10⁻¹³ rad at the surface, scaling with distance). Longitude is unconstrained on the polar axis.
  - PROJ's own geocentric inverse is only accurate to ~10⁻⁵ m at the surface and degrades with height, so it is a loose cross-check (2·10⁻⁵ m below |r| = 6400 km).
- **Why 1 µm, not 1 nm.** The highest carrier in the plan is Ka band (uplinks to ~31 GHz, λ ≈ 9.7 mm). The §10 phase budget of 0.01 rad is 15 µm of path there; 1 µm costs ≤ 6.5·10⁻⁴ rad (6.5 % of the budget), and 2·10⁻⁵ rad at 1.09 GHz. The delay budget (0.001 samples, ≥ 1 ps = 0.3 mm at the 1 GS/s maximum rate) is 300× looser. Besides, a position is computed once from its definition and used consistently, so the inverse's error only enters where a round trip feeds geometry (NED-defined trajectories, prompt 8), and there it is a fixed offset, not noise. A universal 1 nm bound is finer than one f64 ulp beyond 2²³ m ≈ 8400 km.

### 12.3 NED, body and antenna frames

- **NED** at geodetic (φ, λ, h): origin at that point; axes north `(−sφcλ, −sφsλ, cφ)`, east `(−sλ, cλ, 0)`, down `(−cφcλ, −cφsλ, −sφ)` in ECEF (`C_ecef_ned` has them as columns). Down is the inward **geodetic** normal.
- **Body:** x forward, y right, z down. Attitude is `q_ned_body` relative to the NED frame at the platform's own position.
- **Antenna:** mount orientation `q_body_ant`; the lever arm is the phase-centre offset in body coordinates (`PlatformPose::antenna_position`). The light-time solver starts using it in prompt 10.

### 12.4 Quaternions

- **Hamilton** (ij = k), stored `(w, x, y, z)`. `q_a_b` maps frame-b coordinates to frame a by `v_a = q ⊗ v_b ⊗ q*`; this is simultaneously the active rotation that turns a's axes onto b's and the passive change of coordinates b → a. Composition `q_a_c = q_a_b ⊗ q_b_c` (the `*` operator; right factor applies first). `to_mat3` gives `C_a_b` with `v_a = C_a_b·v_b`. `normalized` picks w ≥ 0.
- **Euler:** aerospace 3-2-1, `q_ned_body = q_z(yaw) ⊗ q_y(pitch) ⊗ q_x(roll)`; positive yaw turns the nose east of north, positive pitch raises it, positive roll lowers the right wing. Pinned by tests (heading 090 → nose east; pitch +10° → nose (cos10°, 0, −sin10°); roll +30° → right wing (0, cos30°, sin30°)). `to_euler_zyx` takes pitch as atan2(−C₂₀, hypot(C₀₀, C₁₀)); when cos(pitch) < 10⁻⁸ (gimbal lock) roll is 0 and yaw = atan2(−C₀₁, C₁₁). The rebuilt rotation matches to ~2·10⁻⁸ near the lock and ~10⁻¹⁵ elsewhere.

### 12.5 Geoid

- **Format:** GeographicLib PGM (`egm96-5.pgm`, `egm2008-1.pgm`, …): `P5`, `# Offset` and `# Scale` comments (required; netpbm comment rules: `#` to the next CR or LF, anywhere in the header including directly after maxval), width × (width/2 + 1) big-endian u16 samples, maxval 65535, rows from 90° N to 90° S inclusive, columns from 0° E at 360°/width spacing. N = Offset + Scale·v. Anything else (shape, maxval, missing or extra bytes) is an error.
- **Interpolation:** bilinear, in GeographicLib's order: `a = (1−fx)v₀₀ + fx·v₀₁`, `b = (1−fx)v₁₀ + fx·v₁₁`, `c = (1−fy)a + fy·b`, longitude wrapped to [0°, 360°) with the last column joined to the first, and the south-pole row folded into the cell above it. The model error is the grid's own (EGM96 5′ bilinear: 0.14 m max).
- **Configuration:** `[world] geoid = { path }`. **Without a grid, N = 0 everywhere** (MSL = ellipsoid): the resolved scenario records `world.geoid = { model = "zero" }` explicitly, and the run manifest carries the warning `no geoid grid configured: …`; `sim run` prints it.
- **Manifest:** a grid registers asset `geoid` (kind `geoid_grid`, format `geographiclib-pgm`), with the file SHA-256 and length, and a metadata hash over width, height, Offset, Scale and interpolation. The resolved scenario also carries the file hash, so a changed grid fails replay at asset verification, before any rendering.

### 12.6 Pressure altitude and Mode C

- **ISA**, U.S. Standard Atmosphere 1976 form (the ICAO Doc 7488 layers below 51 km): g₀ = 9.806 65, R = R*/M₀ = 8314.32/28.9644 (Doc 7488 rounds R to 287.052 87, a 7·10⁻⁷ relative difference, ≈ 8 mm of pressure altitude at FL400), T₀ = 288.15 K, p₀ = 101 325 Pa; geopotential layer bases −5, 0, 11, 20, 32, 47 km (top 51 km) with lapse rates −6.5, −6.5, 0, +1, +2.8, 0 K/km. Base pressures are computed from p₀ and match the published 22 632.06, 5474.889, 868.0187, 110.9063, 66.938 87 Pa within 2·10⁻⁷ relative. Closed-form inverse per layer.
- **Geopotential:** H = r₀h/(r₀ + h), r₀ = 6 356 766 m. ISA (and pressure) altitudes are geopotential; platform heights are geometric.
- **Local pressure offset** Δp = QNH − 1013.25 hPa (`[world] pressure_offset`, |Δp| ≤ 200 hPa, default 0): the local atmosphere is the ISA temperature profile with sea-level pressure p₀ + Δp, which scales pressure at every height by (p₀ + Δp)/p₀.
- **Chain** (`AltitudeModel`): ellipsoid h → MSL (h − N) → geopotential → local pressure → ISA pressure altitude → Gillham. The inverse resolves `FLnnn` scenario positions (pressure altitude nnn × 100 ft) to ellipsoid heights.
- **Gillham:** 500 ft index n = ⌊(alt + 1200)/500⌋ as a reflected Gray code on D2 D4 A1 A2 A4 B1 B2 B4; C1 C2 C4 = 001, 011, 010, 110, 100 for the 100 ft step m = 0…4 within the band, reversed when n is odd; D1 always 0; C patterns 000, 101, 111 illegal. Range −1200 ft to 126 700 ft, the full published table (ICAO's nominal range starts at −1000 ft; the two lower codes exist and are accepted). Encoding rounds to the nearest 100 ft with transitions at ±50 ft and ties toward +∞. All 1280 rows of the AirSport Avionics published table match.

## 13. Trajectories and attitude (settled, prompt 8)

All in `iqsim-geo::trajectory`, f64 through `libm`, separate multiplies and adds in a fixed order, as in §12. A trajectory is control-plane truth: prompt 9 fits render polynomials to it.

### 13.1 Trait and state

- `Trajectory::state(t: AcqInstant) -> State` is pure and deterministic, for any instant in `span()` (fractional picoseconds allowed); outside it is `GeoError::OutsideSpan`. `knots(from, to)` lists every instant where jerk (and so angular rate) may be discontinuous; prompt 9 splits fits there.
- `State`: ECEF position (m), Earth-relative velocity and acceleration in ECEF axes (m/s, m/s²), `q_ecef_body` (§12.4 convention; `q_ned_body()` and `pose()` derive the NED form), and the angular velocity of the body **relative to ECEF**, in **body axes** (rad/s). Add Earth rotation for inertial rates.
- Inside a piece, local time is the exact Q64.64 difference `t − start`, rounded once to f64 ps and divided by 10¹², so no absolute time meets a floating-point product.

### 13.2 Representation

- Every moving trajectory (`Flight`) is a `Path` of `Piece`s with **knots at whole picoseconds**. A piece is a quintic per ECEF axis in local seconds (tracks use its cubic subset). Velocity, acceleration and jerk are its derivatives, so they are always consistent with position.
- Maneuver pieces are quintic Hermite interpolants of (p, v, a) at both ends, and consecutive pieces share their knot values, so motion is **C2** (to f64 rounding: over the 356 knots of the test flight the largest jumps are 0 m, 2·10⁻¹² m/s and 6·10⁻¹² m/s²). Jerk is discontinuous at knots.
- Nodes are placed at:
  - segment ends;
  - every channel breakpoint, rounded to the nearest ps;
  - extra nodes so spacing is at most **1 s**, and at most **0.1 rad of turn** (spacing ≤ 0.1/|ω|max, for fast turns at low speed).
- Nodes closer than **1 ms** merge, which leaves that breakpoint just inside a piece.
- **Segments** last at least 1 ms: shorter transitions are held to 1 ms, and an explicit duration under 1 ms is an error. A segment with nothing to do (a change to the current speed or heading) adds no piece.
- **Short pieces** (< 0.1 s) take their displacement from the end velocities and accelerations, `Δp = T(v₀ + v₁)/2 + T²(a₀ − a₁)/12` (exact for quartic motion), instead of `p₁ − p₀`. Otherwise the ~10⁻⁹ m rounding of Earth-sized positions would be divided by T² and T³. They end within that rounding (plus < 10⁻¹¹ m truncation) of the next knot's position, with velocity and acceleration matched exactly. Pieces of 0.1 s or more use `p₁ − p₀`, which keeps acceleration rounding noise below 6·10⁻⁶ m/s².
- Piece lengths are computed in i128, and a non-positive length is an error.
- The Hermite error against the ODE solution is O(h⁶) (≈ 2·10⁻⁹ m for a 3°/s turn at 1 s spacing). It is not a model error: the pieces *are* the truth.

### 13.3 Horizontal motion (n-vector ODE)

- State `(n, d)`:
  - `n` is the geodetic up normal; with height h, `point(n, h) = (N n_x, N n_y, N(1 − e²) n_z) + h n`.
  - `d` is the horizontal unit track direction.
- Neither has a pole singularity.
- ODE:
  - `ṅ = K(n, h)·V_h d`, with `K x = x/(N + h) + c (w·x) w`, `w = ẑ − n_z n` and `c = a e²/(W³(M + h)(N + h))`. This equals the meridian and prime-vertical rates and is smooth at the poles.
  - `ḋ = −ω (n × d) − (d·ṅ) n`.
  - ω > 0 turns right (clockwise seen from above).
  - ω = 0 is parallel transport: **"straight" means zero lateral acceleration**, a geodesic of the height surface, not a rhumb line. The right-hand lateral acceleration is exactly `V_h ω + ḣ (ṅ·right)`; the second term is ~10⁻⁹ of the first.
- Integration: classic RK4 with steps ≤ min(0.125 s, 0.02 rad/|ω|max), renormalizing n and re-orthogonalizing d after every step.
- `V_h = √(V² − ḣ²)` (V is the total speed). It must stay ≥ **1 m/s**.

### 13.4 Channels and primitives

- **Channels.** Turn rate ω is continuous and piecewise linear, ramping at `g₀·roll_rate/V` (V at the segment start). Speed V and vertical rate ḣ are S-curves: their derivative is continuous, piecewise linear, bounded by the acceleration limit, and changes at the jerk limit. h is the exact integral of ḣ (piecewise cubic).
- **S-curve** from (x₀, ẋ₀) to x₁ with ẋ = 0 at the end:
  - The direction comes from the stopping point `x₀ + ẋ₀|ẋ₀|/(2J)`.
  - The peak is `p² = J·D + b₀²/2`, clamped to A with a hold of `(D − (2A² − b₀²)/(2J))/A`.
  - An initial |ẋ₀| above A raises the limit to |ẋ₀|.
- **Move by Δ** (turn by angle, climb to altitude): rise to a peak, hold, fall to 0, with total integral Δ. Either the peak is ±peak_max with a closed-form hold, or it is a root of `integral(peak) = Δ`.
  - The integral is not monotone in the peak when the start has non-zero rate and acceleration. So the range is scanned in 64 steps from ±peak_max towards ∓peak_max, and the first sign change is refined by 200 bisection steps (the root with the largest peak towards Δ).
  - No root is an error: the start state overshoots Δ whatever the peak.
- **Primitives:**
  - `straight {duration}`
  - `turn {rate | bank; duration (total) | heading | by; direction}`
  - `climb {vertical_rate | angle; altitude | duration}`
  - `speed {speed}`

  Each primitive is one segment that commands all three channels. Channels it does not name return to rest: ω → 0, ḣ → 0, and V̇ → 0 at the speed reached. A segment lasts until its channels finish, or its duration if that is later. Turns and climbs by duration fail if the duration is shorter than their transitions.
  - A bank-specified turn uses `ω = g₀ tan φ / V` at the start speed. Banks of 80° or more are rejected.
  - A climb by angle uses `ḣ = V sin γ` at the start speed.
  - Altitude targets convert to ellipsoid heights through the `AltitudeModel` (and the `Ground` for AGL) **at the point where the climb starts**. The platform then holds constant ellipsoid height, so a flight level is not followed across geoid or pressure changes.
  - **Turn to heading** shoots on the turn angle, `Δ ← Δ + wrap(target − final track)`, to 10⁻¹² rad, because the track also drifts over the curved Earth (`χ̇ = ω + λ̇ sin φ`). The default direction is the shortest way. **Turn by** is relative to a straight path, so the final track differs from start + Δ by that drift.
- **Defaults** (`Limits`; `[platforms.trajectory.limits]`):

  | Limit | Default |
  | --- | --- |
  | route bank | 25° |
  | roll rate | 10°/s |
  | accel | 1 m/s² |
  | jerk | 0.5 m/s³ |
  | vertical accel | 0.5 m/s² |
  | vertical jerk | 0.25 m/s³ |
  | route climb rate | 1500 ft/min = 7.62 m/s |

### 13.5 Derived attitude

- One function of (p, v, a, j) for every moving source, so attitude always agrees with the kinematics:
  - n comes from the geodetic inverse, with `ḣ = v·n`, `v_h = v − ḣn`, `ṅ = K v_h`, `ḧ = a·n + v·ṅ`;
  - with a constant wind w (NED, rotated with the local frame), the air-relative horizontal velocity is `u = v_h − w`;
  - the heading frame H is x along u, y right, z down;
  - **pitch** `θ = atan2(ḣ, |u|) + α`;
  - **bank** `φ = atan2(a·y_H, g₀)` with g₀ = 9.806 65 m/s², i.e. `tan φ = V ω / g₀` in a coordinated turn;
  - `q_ecef_body = q_ecef_H ⊗ q_y(θ) ⊗ q_x(φ)`, which is heading (track plus crab), pitch, roll in the aerospace 3-2-1 sequence;
  - the angular rate is `Ω_H + θ̇ y_H + φ̇ x_body`, with `Ω_H = ½ Σ hᵢ × ḣᵢ`, all in closed form (φ̇ uses the jerk).
- Approximations:
  - α is applied as a pitch offset before bank;
  - bank uses the Earth-fixed lateral acceleration, so Coriolis (≈ 0.04 m/s² at 250 m/s) and the curvature's vertical acceleration (V²/R ≈ 0.01 m/s²) are excluded, as the design's formula says.
- Errors:
  - the heading is undefined below 0.1 m/s air-relative horizontal speed (an error at evaluation);
  - wind is rejected within 0.01° of a pole.
- Modes:
  - `Derived {angle_of_attack, wind}`: the airborne kinds.
  - `Ground`: no bank and no α; pitch follows the along-track slope.
  - `Fixed(q_ned_body)`: static platforms only.
- Measured:
  - bank within 0.01° of `atan(Vω/g₀)` at five latitudes up to 89.5° and three speeds;
  - the angular rate matches quaternion differences to 2·10⁻⁶ rad/s.

### 13.6 Routes

- Waypoints `{lat, lon, alt, speed?}` are flown over. Each leg is one segment from the previous waypoint (or the start):
  - turn by Δ at the route bank;
  - climb or descend at the route climb rate to the waypoint altitude (converted at the waypoint), starting at the leg start and level by the waypoint;
  - change to the waypoint speed (kept if not given);
  - straight for the rest.
- Shooting on (Δ, leg duration T):
  - From the end state, the along-track miss `e_a` and the cross-track miss `e_c` update `T += e_a/V_h` and `Δ += e_c/L`, with L the distance flown after the turn.
  - The misses are the chord projected on the local horizontal when the normals are within 10⁻³ rad. Beyond that they are the surface arc `angle(n, n_w)·(N + h)` along the initial direction to the waypoint, so a far waypoint is never mistaken for an arrival. An antipodal waypoint (no defined direction) is an error.
  - It stops at **10⁻⁴ m** in both, and at 10⁻³ m in 3-D distance. It fails after 40 iterations, or when the waypoint is too close for the transitions.
- The initial heading defaults to the bearing of the first waypoint.
- Tested passing within 1 cm horizontally and 1 cm vertically, with the speed reached, in three cases: a four-waypoint route in England, a dateline crossing, and a route over the North Pole.

### 13.7 Track replay

- **Reinsch** natural cubic smoothing spline per ECEF axis (relative to the first point):
  - It minimizes `Σ |y_i − f(t_i)|² + λ ∫ |f″|² dt`.
  - The solve is a banded LDLᵀ in a fixed order.
  - Knots are at the timestamps.
  - The span is [t₀, tₙ].
- **Smoothing parameter.** `smoothing = τ` (default 0 s, which interpolates exactly), with `λ = τ⁴/Δt̄`, where Δt̄ is the mean spacing. For uniform sampling the spline acts on motion at angular frequency ω as the zero-phase gain `1/(1 + (ωτ)⁴)`. Measured: 0.500 000 03 at ωτ = 1 and 0.941 176 48 at ωτ = 0.5.
- Points convert with the `AltitudeModel` at each point.
- Replay timing:
  - at τ = 0, positions at the timestamps match within 1 µm;
  - at τ = 2 s, the closest approach to each point is at its timestamp within 10 ms (mean offset < 0.2 ms), with local shifts only where the recorded motion's jerk steps.

### 13.8 Static, ground and appends

- **`Static`:** fixed position and attitude relative to NED (`{heading, pitch, roll}`, default 0), zero velocity, acceleration and rate, unbounded span.
- **Ground-clamped (placeholder):** maneuvers without climbs, or routes whose heights are AGL over a `Ground` surface. Until terrain (prompt 22) the surface is `FlatGround`, which the runtime sets to the ellipsoid, consistent with AGL positions (§11.4). A terrain model must be C2 along the path for the clamp to stay C2.
- **Append** (`Flight::append(t, primitives)`, `append_route(t, waypoints)`, used by runtime maneuvers in prompt 16):
  - It replaces everything planned after t.
  - The piece containing t keeps its coefficients and only ends at t, so **every instant before t evaluates to the same bits**.
  - The new channel and ODE state is the exact inverse of the node formulas applied to the (p, v, a) at t (`ω = (a·right − ḣ ṅ·right)/V_h`, …).
  - The first new piece starts from that exact (p, v, a), so the join is C2. The geodetic round trip (~10⁻⁷ m) is absorbed inside that first piece.
  - Appends work after track prefixes too.
  - On error the flight is unchanged.
- **`extend_to(t)`** appends straight and level until at least t.

## 14. Light time and per-block polynomial fitting (settled, prompt 9)

`iqsim_geo::light_time` (the solver), `iqsim_geo::trajectory::Lookback`, `iqsim_runtime::geometry` (exact model, fit, checks, refinement, log) and `iqsim_runtime::scene` (work lists). Control plane, f64 through `libm`, fixed operation order; the render path still sees only fixed-point polynomials.

### 14.1 Light time

- Every path solves `|r_rx(t) − r_tx(t − τ)| = c·τ` in the pseudo-inertial ECI frame (§12.1), each endpoint rotated at its own instant (`EarthOrientation::at`), velocities including ω × r. Earth rotation during the light time (the Sagnac term) is therefore included: ≈ 155 ps on a 30 km east–west link at the equator.
- `solve_receive(rx, tx, t)`: fixed-point iteration from τ₀ = |r_rx(t) − r_tx(t)|/c, `τₙ₊₁ = |r_rx(t) − r_tx(t − τₙ)|/c`, with `t − τ` formed as Q64.64 ps (`instant_offset`, τ rounded once to f64 ps then scaled exactly by 2⁶⁴; |τ| ≤ 10⁶ s), so no absolute time meets a float product. **Stopping rule:** stop when `|τₙ₊₁ − τₙ| ≤ max(10⁻¹⁶ s, 4·ulp(τₙ₊₁))`, at most 10 updates (`NoConvergence` beyond); the transmitter is then evaluated once more at `t − τₙ₊₁`, and τ is reported as that distance / c, so τ, the range and the endpoint states agree. The map contracts with factor β ≈ |v_tx|/c (3·10⁻⁵ for LEO), so the remaining error is ≈ β/(1−β) times the last step (< 10⁻²⁰ s); LEO and GEO converge in 3 updates (tests assert ≤ 5). Coincident endpoints are an error.
- `tau_rate` = dτ/dt = `û·(v_rx − v_tx)/(c − û·v_tx)` (ECI, û from transmitter to receiver). The received carrier's Doppler shift is `−f_e·τ̇`.
- `solve_emit(rx, tx, t_e)` solves `|r_rx(t_e + τ) − r_tx(t_e)| = c·τ` the same way; it maps a transmitter knot (emission time) to its receive time.
- The first guess evaluates the moving endpoint at the fixed instant **clamped into its span** (review fix: a transmitter knot before the receiver's start, or a reception after the transmitter's end, made the first guess fail although the solution was valid).
- **Phase centres:** platform origins until prompt 10 (antenna phase centres and lever arms).
- **Tested** (`crates/iqsim-geo/tests/light_time.rs`): two Earth-fixed equator points against an independent bisection of `cτ = 2R·sin(|Δλ − ω_E·τ|/2)`, both directions (Sagnac of opposite sign) at 30 km and 3000 km and at 100 days; uniform ECI motion (7.5 km/s) against the closed-form quadratic root and its derivative: worst 5.6·10⁻¹⁸ s (bound 10⁻¹⁵ s), τ̇ within 2.8·10⁻¹³ of a finite difference; `solve_emit` inverts `solve_receive` within 10⁻³ ps for LEO and GEO (measured 0).

### 14.2 Look-back before a trajectory start

Light time evaluates the transmitter before t = 0, which is before a `Flight`'s span start. **Decision (with the user): extrapolate.** Every transmitter trajectory is wrapped in `Lookback` (1 s, `MAX_LOOKBACK`): for t in [start − 1 s, start) the state is the second-order Taylor extrapolation of the start state, `p₀ + (v₀·dt + (a₀·(½dt))·dt)`, `v₀ + a₀·dt`, `a₀`, with attitude and angular rate held, so the join is C2; the start itself is reported as a knot (a jerk step), so no fit crosses it. A trajectory without a start (`Static`) passes through unchanged. Before start − 1 s it is `OutsideSpan`. Receivers are evaluated only at their own acquisition instants and are not wrapped.

### 14.3 Exact model, fit and refinement

For sample k of a block at acquisition instant t (`ClockMap.instant`, the Q64.64 fraction kept), with τ and τ̇ from §14.1:

| Quantity | Exact model | Representation |
| --- | --- | --- |
| delay | τ·F_s target samples | f64 |
| phase | `phase_at(f_off, t′) − carrier_phase(f_e, τ(t′))` | u64, 2⁻⁶⁴ cycle, modular |
| coordinate | `q64_interval_to_q32(t − t_ref, F) − round(F·τ·2³²)` | i128, Q32.32 in coordinate units |
| amplitude | `V_ref·g`, `V_ref = √(P·R₀)·λ/(4π·c·τ)` | f64 volts |

t′ = t, except for CW through a response with group delay, where t′ = t − τ_g (below).

- `carrier_phase(f_e, τ) = ⌊f_µHz·d/10¹⁸⌋ mod 2⁶⁴` with d = τ as Q64.64 ps: an exact 64 × 128-bit product in three limbs and long division by 10¹⁸, so `f_e·τ` keeps its fraction for any link (prompt 5's 2⁴⁰-cycle limit is gone; τ's own f64 precision is what remains).
- **CW fold:** `g = |H(f_inst)|` at the instantaneous offset `f_inst = f_off − f_e·τ̇`, and the response's group delay is applied **exactly**: the whole geometry (phase, τ, amplitude) is evaluated at `t − τ_g`, because a linear-phase response delays an in-band tone's trajectory by τ_g. (Review fix: the first-order fold `−f_inst·τ_g` omitted ½τ_g²·φ̈, measured 0.029 rad for a 31 GHz pass at 7.5 km/s and 550 km with τ_g = 1 ms; the exact form is within 1.0·10⁻⁶ rad there.)
- **Pulses and continuous sources are filtered at the Doppler-shifted offset** (review fix: they were filtered at the nominal offset, which for a 31 GHz emitter closing at 7.5 km/s is 775 kHz off). Each segment's edge tables or resampling kernel are built at `offset_uhz = f_off + round(D/grid)·grid`, D = the mean of f_inst over the nodes minus f_off, on a grid of half the tolerance `doppler_tolerance_hz` = min(¼·A_budget·√(2π)·σ, ¼·P_budget/(2π·τ_g)): max|H′| = 1/(√(2π)σ) bounds the magnitude change, 2π·δ·τ_g the group-delay phase. With no Doppler the key is exactly f_off. Tables and kernels are cached by key and built when a segment first needs them, before its block renders; a pulse's activity is tested on its nominal-offset tables with one sample of margin. The Doppler compression of the source's own spectrum (ν·τ̇, ≤ 25 Hz for 1 MHz at LEO speeds) is not modelled in the key. `g = 1` for them.
- The coordinate origin `t_ref` is the pulse train's emission time (F = F_s) or the continuous source's sample-0 time (F = F_n, rebased to the window per §5.6).
- **Nodes:** an interval of n samples is interpolated through the samples at `u = 0, ¼, ¾, 1` (`j = round(u·(n − 1))`, duplicates dropped, so the degree is ≤ min(3, n − 1)): the degree-3 Chebyshev–Lobatto points, which include the ends, so adjacent segments and blocks meet at exactly evaluated values. Each quantity's variation from the first node is interpolated (Newton divided differences, expanded to the power basis in j, then forward differences `D1 = a1 + a2 + a3`, `D2 = 2a2 + 6a3`, `D3 = 6a3`). Delay D0 = τ₀·F_s; phase D0 = the exact phase at j = 0 (high word); the coordinate's exact Q32.32 start is kept as an i128 beside a D0 = 0 variation polynomial and placed into D0 after rebasing (|w| < 2³⁰ samples, §5.2). Amplitude: an f32 cubic in the kernel's `x = f32(j)·f32(1/n)` through the same nodes (constant when all node values are equal).
- **Checks:** at 15 points, `u = 0, 1/16, ½ − √(5/32), 3/16, ¼, 5/16, 7/16, ½, 9/16, 11/16, ¾, 13/16, ½ + √(5/32), 15/16, 1`: the ends, the inner nodes, the three extrema of the nodal polynomial `u(u − ¼)(u − ¾)(u − 1)` (where a smooth function's interpolation error peaks; |ω| = 1/64 at ½ and 9/1024 at the outer two) and the odd sixteenths. The fixed-point polynomials are evaluated exactly as the kernel does (jump, floor extraction, f32 Horner), against the exact model:
  - delay: |Δ| in target samples;
  - phase: the wrapped difference of the two u64 phases, radians;
  - coordinate: |Δ| in target samples (native units × F_s/F);
  - amplitude: |Δa|/a above the null floor; below it |Δa| over 10⁻⁵·V_ref/budget, so acceptance there is the §2 absolute bound (review fix: a deep-stopband value that underflows f32 had counted as a 100 % error);
  - response (`Residuals::response`): for pulses and continuous sources |f_inst − key|/tolerance; for CW the change of f_inst between adjacent check points over σ/8, so |H(f_inst)| is resolved by the checks (review fix: a sweep through a sharp edge could vanish between checks). Accepted at ≤ 1.
- **Acceptance:** every residual ≤ budget / `FIT_MARGIN` (4): delay 2.5·10⁻⁴ samples, phase 2.5·10⁻³ rad, coordinate 2.5·10⁻⁴ target samples, amplitude 2.88·10⁻⁴ (¼ of 0.01 dB). Phase is checked at every point (null floors arrive with antenna gains in prompt 10).
- **Why this is conservative.** On an interval without known discontinuities the interpolation error is ≈ f⁽⁴⁾(ξ)/24·ω(u)·T⁴, and the checks include ω's extrema, so the error between checks exceeds the largest checked residual only through the variation of f⁽⁴⁾ across the interval. The factor 4 covers f⁽⁴⁾ changing by up to 4× within one interval. Point checks remain diagnostics, not a proof between points (plan, "Error budget").
- **Refinement:** an interval that misses acceptance is halved and both halves are fitted and checked again, depth first, left half first, at most `MAX_REFINE_DEPTH` = 12 halvings (4096 segments per knot-free interval); missing it there is `FitError::Budget` with the interval and residuals, and the block (and the run) fails. Halving an interval cuts the cubic interpolation error by ≈ 16.
- **Knot splits first:** before fitting, the block is split at every receiver knot (receive time) and every transmitter knot mapped to receive time by `solve_emit`, at the first sample whose instant is at or after it (`index_at_or_after`); a split at a block edge needs none. The transmitter knots searched are those between the emission instants of the block's first and last samples. Jerk is discontinuous at knots, so no fitted interval crosses one. Future effective parameter changes and blockage transitions (prompt 10) will enter the same list.
- **Delay range:** the fitted delay cubic must stay inside ±(2³¹ − 1) samples over the whole interval, its interior extrema included (review fix: node checks alone missed an interior crossing that wraps in Q32.32), else `FitError::DelayRange`.
- **Coordinates far from their origin:** the interval `t − t_ref` is taken as whole ps (i128) plus the fraction (review fix: two instants ±100 days apart overflowed one Q64.64 i128).
- **Measured** (`crates/iqsim-runtime/tests/geometry.rs`, against independent light-time references that solve the closed-form ECI motion directly):
  - moving CW (uniform ECI motion at 7.5 km/s passing 50 km over the receiver, 2.2 GHz, 2 MS/s, 10 ms blocks; Doppler −12.5 to +26.4 kHz): rendered phase within **1.8·10⁻⁷ rad** of `f_off·t − f_e·τ(t)` at every sample of six blocks through closest approach; frequency over 1 ms windows within 2.2·10⁻⁵ Hz;
  - LEO pass (circular orbit at 550 km passing overhead, 10 GHz, 1 MS/s, 4 s blocks, a CW and a pulse emitter): one cubic per block misses the budget; 3090 refinement splits over three blocks, most of them on the pulse path, whose 14 kHz Doppler sweep per block must stay within the 8.7 Hz tolerance of a 400/500 kHz response; at ≥ 98 points per segment the worst errors are delay 5.2·10⁻⁹ samples, phase 3.3·10⁻⁴ rad, coordinate 1.0·10⁻⁹ samples, amplitude 7.3·10⁻⁸;
  - review regressions: a CW sweep at ≈ 1 kHz/s through a 100/200 Hz response edge refines until every sample's amplitude is within 0.23 of its budget; a pulse from a 31 GHz emitter closing at 7.5 km/s, nominally 775.5 kHz off, arrives at the channel centre at its Friis level within 2·10⁻⁷ dB; a moving receiver whose plan ends inside the final block, with a 100 µs clock offset, renders;
  - an aircraft whose turn knots fall inside blocks: each transmitter-knot split is the first sample whose emission time is at or after the knot; output is identical across 1, 3 and 8 threads and chunk lengths 777 and 1.

### 14.4 Work lists, culling, headroom and the fit log

- Each fitted interval becomes one `PathSegment` (§5.2). Culling is per segment (lossless at the default threshold, as in §2: every component of a culled segment rounds to 0). The headroom proof adds, per path, the **largest** bound over its kept segments, because segments of one path never overlap in time.
- Pulse activity uses the coordinate over the whole block (the first segment's start to the last segment's end; the coordinate is monotone). A continuous path's window covers all of its segments' coordinates, and each segment's start is rebased to the window origin.
- Pulses get one `PulseTrain` per segment, with the tables at that segment's offset, and continuous segments are grouped into windows by their segment's kernel.
- **`stop` sources:** exhaustion is checked on the **exact** retarded coordinate of the last kept sample (review fix: the fitted one can fall short by the fit residual).
- **Fit log:** `SceneBuilder::build_block_logged` returns every path's segments (with residuals) and splits (`rx_knot`, `tx_knot`, `refine` with the halved interval's residuals). The offline runner aggregates them per (stream, emitter) into the truth record's `geometry_fit`: blocks, segments, worst residuals, and every split with its block, sample and reason. The log is deterministic and covered by the truth hash.

### 14.5 Resolution

- Every platform's trajectory is built at resolution and, if its plan ends earlier, extended straight and level (`trajectory::build_until`, `Flight::extend_to`) to the last instant the run can evaluate: duration + one block (the final block is fitted whole) + the largest receiver clock offset + 1 ms (review fix: extending to the nominal duration left the final block's tail outside the span).
- The resolved scenario's `platforms` record kind, span, knot count and `trajectory_sha256` = SHA-256(`"iqsim.trajectory-fingerprint.v1\0"` ‖ per instant: ps (i64 LE) ‖ position, velocity, acceleration, q_ecef_body (w, x, y, z), angular rate as f64 bits LE) over max(span start, −1 s), every knot before that last instant, and the last instant (or the span end).
