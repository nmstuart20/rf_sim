# Claude Code prompts

Thirty prompts, run in order (0–6, 6a, then 7–28), take the design from an empty repository to phase 12. Prompts 0–25 and 6a contain concrete implementation steps; 26–28 plan first and get refined once earlier phases land.

## How to run these

- Keep `project_plan.md` and `agent_prompts.md` at the repository root. `project_plan.md` is the canonical design; every prompt cites its sections by heading.
- Run one prompt per fresh Claude Code session (or `/clear` between prompts), in order. Start each in plan mode and approve the plan before any code is written.
- Prompt 0 writes `CLAUDE.md` with the project invariants and a Definition of Done. Later prompts point to it instead of repeating it.
- Prompt 0 also writes `docs/contracts.md` with the numeric, clock, commit-boundary, replay, tap and transport contracts. Every session reads it; phase 1 fixes the proposed tolerances before rendering begins.
- Each prompt ends by appending to `docs/progress.md`: what was built, deviations from the design, open questions. The next session reads it first.
- After each prompt: review the diff, run `cargo xtask ci`, commit. Golden hashes change only through `cargo xtask bless`, and only on purpose.
- `iqsim` is a working name for the workspace and `sim` for the binary; rename them before prompt 0 if you like.
- Prompts 20–21 need an NVIDIA GPU and 24–25 need an SDR. The rest run on any Linux machine.

## Setup

**Prompt 0: Repository bootstrap**

```text
Read project_plan.md end to end. You are starting a new Rust workspace, "iqsim", for the real-time scenario IQ simulator it describes. Do not implement simulation logic yet.

1. Create a Cargo workspace (stable toolchain pinned in rust-toolchain.toml) with these library crates, each with a crate-level doc stating its responsibility from the design:
   iqsim-core (time, fixed-point, phase, sin/cos, RNG), iqsim-units, iqsim-scenario, iqsim-geo, iqsim-antenna, iqsim-terrain, iqsim-behavior, iqsim-waveform, iqsim-worklist, iqsim-render (backend trait + CPU backend), iqsim-render-cuda, iqsim-receiver, iqsim-io (sinks and exports), iqsim-control, iqsim-runtime, iqsim-hil.
   Add a binary crate iqsim-cli that builds `sim`, and an xtask crate.
2. xtask commands: `ci` (fmt check, clippy with -D warnings, all tests), `bless` (regenerate golden hashes), `gpu-test` (placeholder for now).
3. Golden-test helper: assert_golden(name, bytes) compares a SHA-256 against tests/golden/<name>.sha256 and fails with a clear message; `cargo xtask bless` rewrites them.
4. Write CLAUDE.md containing:
   - Project summary; read project_plan.md, docs/contracts.md and docs/progress.md at the start of every session after bootstrap.
   - Invariants:
     * No floating point for time, delay or carrier phase on the render path.
     * Integer (i64) accumulation with documented physical units, scale and checked headroom. Emitter partitions return raw sums; receiver effects and final conversion occur once after global reduction.
     * No std f64 transcendentals (sin, cos, exp, ln, atan2, powf) in control-plane code; use the libm crate so results match across machines and libc versions.
     * No fast-math, no FMA contraction, no nightly features.
     * #![forbid(unsafe_code)] everywhere except iqsim-render-cuda and iqsim-hil.
     * Deterministic iteration: no HashMap/HashSet iteration order in outputs; use BTreeMap or sorted Vecs.
     * Every random draw comes from the Philox streams in iqsim-core.
     * Actual sampling instants, reported timestamps and LO errors are distinct models. Use the ClockMap interface from the first render implementation.
     * Antenna-port rendering uses B_out; virtual-receiver rendering uses H_rx. No virtual receiver effects enter an antenna-port output.
     * Commands are admitted atomically at the shared exclusive commit boundary, covering behavior lookahead and all committed dependencies; replay uses logged effective times.
     * Replay validates a run manifest and every external input hash. Truth and manifests stay separate from consumer metadata.
     * Strict delivery uses acknowledged framed streams; raw PUB is best effort. Budget aggregate traffic on shared resources.
   - Dependency policy: few, well-maintained crates that vendor cleanly for an air-gapped RHEL 9 build; justify each new dependency in your summary; no network access at build time.
   - Definition of Done: `cargo xtask ci` passes; new behavior has unit tests and, where project_plan.md names a validation check, a test for that check; public items documented; docs/progress.md appended with what was built, deviations from the design and why, and open questions.
5. Write docs/contracts.md from the corresponding sections of project_plan.md: numeric operation/rounding rules, amplitude normalization and overflow policy; ideal and perturbed ClockMap interfaces; exclusive commit-boundary semantics; tap responses and state ownership; run-manifest fields; transport choice/alternatives and delivery guarantees. Record proposed tolerances and the decisions phase 1 must settle, including the amplitude null floor. Keep protocols independent of simulation handlers. A dependency graph must keep these shared contracts usable without a runtime/control dependency cycle.
6. Create docs/progress.md with a "Prompt 0" entry.

Stop when the workspace builds, `cargo xtask ci` passes and CLAUDE.md and docs/contracts.md are complete. Summarize the layout and any decisions you made.
```

## Phase 1: Foundations

**Prompt 1: Time and rate types**

```text
Read CLAUDE.md, docs/progress.md, and project_plan.md section "Time, numerics and determinism". Implement in iqsim-core:

- TruthTime: i64 picoseconds with checked arithmetic. Overflow is an error, never a wrap.
- SampleRate: an exact rational (u64 numerator and denominator, reduced). Constructors from integer Hz and from rationals; reject zero.
- SampleIndex: u64 per receiver channel. Exact conversions between SampleIndex and TruthTime given a channel epoch and SampleRate, with documented rounding (floor for time to index; exact for index to time when representable, otherwise floor plus the returned remainder).
- ClockMap interface: sample index -> actual acquisition truth time with its sub-picosecond remainder, and a documented inverse. Implement the ideal rational map now; prompt 15 adds monotone oscillator-error models behind this interface. Keep the timestamp-clock mapping separate. Geometry must retain the remainder rather than floor every acquisition instant.
- Block: a half-open range of sample indices for one channel; BlockSpec (block length in samples) and an iterator over blocks.
- Use u128/i128 intermediates where needed, and document why no intermediate can overflow for rates up to 1 GS/s and runs up to 100 days.
- Fix the time/phase tolerances in docs/contracts.md and validate or explicitly bound rational numerator/denominator sizes so the intermediate-overflow argument covers every accepted SampleRate.

Tests: property tests (proptest) for round trips, monotonicity and block tiling (no gaps, no overlaps), plus edge cases at 61.44 MS/s and 250 MS/s. Follow the Definition of Done in CLAUDE.md.
```

**Prompt 2: Fixed-point delay, carrier phase, polynomials and sin/cos**

```text
Read CLAUDE.md, docs/progress.md, and project_plan.md "Time, numerics and determinism". In iqsim-core, implement the numeric core every render backend shares. The CUDA backend (prompt 20) will port this code line by line, so write it with that in mind.

- Delay output: Q32.32 fixed-point samples in i64.
- Phase output: u64 fraction of a cycle (2^64 = one cycle), wrapping arithmetic.
- BlockPoly: a polynomial of order <= 3 in local sample index k (0..block_len), stored in Newton forward-difference form (D0..D3) as 128-bit fixed point, with 64 extra fractional bits beyond the outputs above.
  * Jump: p(k) = D0 + C(k,1)*D1 + C(k,2)*D2 + C(k,3)*D3 in wrapping integer arithmetic, exact for any k.
  * Step: advance one sample with three 128-bit adds.
  * Document the maximum block length (C(k,3) must fit in u64), and show that coefficient quantization error stays far below 0.001 samples of delay and 0.01 rad of phase over a full block.
- Implement the 128-bit arithmetic with explicit 64-bit primitives (mul_hi_u64, wrapping mul, add with carry) in one small module rather than native i128, so the CUDA C version can mirror it exactly. Test it against an i128 reference.
- sincos(phase) -> (f32, f32): the top phase bits index a table built by a documented deterministic routine; a fixed-order polynomial corrects the remainder. Write every operation explicitly (separate mul and add). Report max absolute error over a sweep of the top 24 bits with random low bits; target <= 3e-7.
- Pin the numeric contract in docs/contracts.md: ties-to-even conversion, finite-value validation, subnormal handling, precise division/sqrt and compiler flags on CPU/CUDA. Separate intentional polynomial/phase wrapping from forbidden accumulator overflow. Settle amplitude tolerances and the null floor before prompt 5. Provide modular carrier-offset phase from exact time without a large floating-point absolute-time product.

Tests: jumping and stepping agree bit for bit over full blocks, including phase wrap and negative delay rates; modular offset phase retains precision at times near 100 days; rounding boundary cases and non-finite rejection; golden hash over a fixed set of 1e6 evaluations. Follow the Definition of Done in CLAUDE.md.
```

**Prompt 3: Philox and Gaussian noise**

```text
Read CLAUDE.md, docs/progress.md, and project_plan.md "Time, numerics and determinism" (Noise). In iqsim-core:

- Implement Philox4x32-10 as specified in the Random123 paper. Take the Philox4x32-10 known-answer vectors from the Random123 distribution's kat_vectors file and include them as tests. Do not invent expected values; if you cannot obtain the file, stop and ask.
- Stream keying: key derived from (scenario seed: u64, entity id: u32, stream id: u32) with a documented, deterministic mixing function; counter = (sample index, sub-index).
- Gaussian: Box-Muller where the angle is a u32 treated as a phase and passed to the sincos from prompt 2, and the radius uses an in-house ln (range reduction plus a fixed polynomial, explicit operations, no std or libm call) and IEEE sqrt, so the GPU port matches bit for bit. Document the error of ln.
- noise(seed, entity, stream, sample_index) -> (f32, f32) with documented variance.

Tests: known-answer vectors; moments plus an Anderson-Darling or chi-square test over 1e7 samples; cross-correlation between adjacent streams; golden hash of 1e6 samples. Follow the Definition of Done in CLAUDE.md.
```

**Prompt 4: Units, scenario schema skeleton and `sim lint`**

```text
Read CLAUDE.md, docs/progress.md, and project_plan.md "Interfaces" (Scenario files, CLI) and "World model". Scenario files are TOML.

1. iqsim-units: typed quantities parsed from strings with units and deserializable with serde:
   Frequency (Hz, kHz, MHz, GHz); Power (W, mW, dBm, dBW); Gain (dB, dBi); Distance (m, km, ft, nmi); Speed (m/s, km/h, kt); Altitude (m, ft, FL, with an explicit reference: ellipsoid, MSL or AGL); Angle (deg, rad, DMS); Duration (s, ms, us, ns); Position (latitude/longitude in decimal degrees, DMS or MGRS, plus altitude).
   Reject unitless numbers where a unit is required, with a helpful error. For MGRS, implement it or pick a vendorable crate, and test against published reference points.
2. iqsim-scenario: serde types for the top level of a scenario: metadata, epoch (UTC), seed, run mode, receivers (channels with center frequency, sample rate, bandwidth and separate clock/LO placeholders), antenna-port output targets, platforms (position, with trajectory and mounts as placeholders), emitters (behavior and waveform references as placeholders), library search paths, and sink/resource identities for shared NICs/disks/USB/PCIe links. Use deny_unknown_fields. Generate JSON Schema with schemars into schema/scenario.schema.json through an xtask command and check it in.
3. `sim lint <file>`: parse and validate; errors point at file:line:column using toml spans; non-zero exit code on error.
4. examples/minimal.toml (one static CW emitter, one receiver) lints clean. Add negative tests for common mistakes: missing unit, unknown field, bad MGRS.

Follow the Definition of Done in CLAUDE.md.
```

## Phase 2: Static render on CPU

**Prompt 5: Work list and CPU render backend**

```text
Read CLAUDE.md, docs/progress.md, and project_plan.md "Architecture", "Time, numerics and determinism", "Receiver model" (Band-limiting before sampling) and "Emitters, behaviors and waveforms" (Rendering paths).

1. iqsim-worklist: plain-data #[repr(C)] types a GPU can consume unchanged, with no pointers (indices into flat arrays only):
   PathSegment (render target/tap, block, delay BlockPoly, retarded waveform-coordinate BlockPoly with declared units/source segment, phase BlockPoly, f32 amplitude polynomial, emission index); PulseTrain (pulse start offsets in Q32.32 target sample units relative to emission time, widths, edge-table references, per-pulse amplitude); CwEmission. A BlockWorkList holds everything one block needs. Waveform coordinates come from ClockMap acquisition instants minus propagation delay, preserving fractional time, rather than ideal sample index minus nominal delay. Define source rebasing so fixed-point coordinates cannot overflow during long runs.
2. Define physical complex-voltage/RMS normalization, reference impedance, output full scale, supported peak overlap and the fixed-point accumulator scale with guard bits. Keep antenna-port units independent of virtual ADC configuration. Document these in iqsim-core, CLAUDE.md and docs/contracts.md. Prove headroom from configured peak path contributions and reject unsafe configurations; i64 accumulator overflow is an error.
3. Responses: add separate H_rx(f) for the virtual receiver and B_out(f) for an antenna-port output target, each with documented passband, transition band, stopband attenuation and group delay. The antenna-port target excludes H_rx and all virtual receiver effects. Build separate tap-specific work lists when responses differ.
4. Band-limiting before sampling, applied per path:
   * Narrowband emitters (CW now, modulated later): fold the selected tap's response at the emitter's offset frequency, gain and phase, into the path's amplitude and phase. Cull paths whose in-band contribution is below a documented threshold.
   * Pulses: complex edge tables (the edge shape convolved with the selected response shifted to the emitter's offset), built deterministically before use (libm, fixed order), keyed by tap, response identity, edge shape and offset. A pulse renders as the difference of two shifted edge responses. Choose the table resolution so interpolation error stays under a documented bound. Define cache invalidation/preparation at future runtime frequency or response changes.
5. iqsim-render: a Backend trait, render(&BlockWorkList, accumulators per channel as &mut [i64]). CPU backend: each worker takes a chunk of output samples, jumps every overlapping path's polynomials to the chunk start, then steps sample by sample: band-limited envelope at the exact delayed time, carrier rotation from the phase polynomial and sincos, amplitude scaling to accumulator units, round to i64 and add. Integer accumulation makes the result independent of chunking and path order.
6. A temporary static-geometry builder using the ideal ClockMap: emitter and receiver at fixed ECEF positions, free-space loss, isotropic antennas, constant delay and carrier-offset phase per block, and retarded waveform-coordinate polynomials.

Tests: a single pulse lands at the expected fractional sample (estimated time of arrival within 0.001 samples after removing the selected filter's group delay); a CW tone lands at the expected frequency and level for each tap response; an emitter 60 MHz outside the target band leaves no aliased image above the declared threshold; output hashes are identical across 1, 2 and 8 threads; invalid accumulator headroom is rejected. Follow the Definition of Done in CLAUDE.md.
```

**Prompt 6: Receiver chain, SigMF output and first validation checks**

```text
Read CLAUDE.md, docs/progress.md, and project_plan.md "Receiver model", "Interfaces" (IQ outputs) and "Validation strategy".

1. iqsim-receiver: the chain from the design, in order, each stage switchable: front-end gain and noise figure; thermal noise (kTBF) from the Philox noise generator, shaped by a digital equivalent of H_rx(f) from prompt 5; LO frequency offset (phase noise as a documented stub); IQ imbalance and DC offset; optional digital decimation filter; AGC off by default; ADC quantization and clipping to sc16. H_rx is already applied to each virtual-receiver path before sampling and is not applied again here. Specify per-channel state ownership, filter history and ordered block handoff. Apply the chain once after the complete path sum. The separately rendered antenna-port target uses B_out and bypasses all virtual receiver effects; give it its own output-level conversion.
2. iqsim-io: a SigMF writer (ci16_le and cf32_le) with core:sample_rate, core:frequency and core:datetime from the ideal receiver timestamp-clock interface; prompt 15 adds impairments. Check write failures, final sample count and sample hash. Truth annotations go to a separate sidecar file, never into the consumer's .sigmf-meta.
3. Versioned run manifest and asset registry: resolved scenario, seed, build/numeric-contract identities, available asset hashes, event-log identity when present and deterministic sample/canonical-truth hashes. External asset types register here as later prompts add them. Keep manifests with privileged run artifacts, separate from consumer files. Hash source bytes and meaningful metadata, and reject missing/changed assets before replay. Runtime measurements are separate from deterministic truth.
4. iqsim-runtime and iqsim-cli: `sim run --offline <scenario> --out <dir>`, wiring scenario, static geometry, both render targets, CPU backend, receiver chain, SigMF writer and run manifest. `sim replay <run-manifest>` initially supports runs without live commands; prompt 16 extends it with the event log.
5. Validation tests from the design's table: delay; power within 0.1 dB of a closed-form link budget; noise floor within 0.2 dB of kTBF plus noise figure over the channel's noise bandwidth (Welch PSD); out-of-band rejection (a 1030 MHz emitter into a 1090 MHz channel leaves no aliased image, and in-band levels match the selected response within 0.1 dB). Change virtual H_rx, gain, noise and clock configuration while holding B_out/output timebase fixed: antenna-port hashes must stay unchanged. Golden hash for examples/minimal.toml, replay equivalence and input-mismatch rejection.

Follow the Definition of Done in CLAUDE.md.
```

**Prompt 6a: IQ-file emitters and CPU resampling**

```text
Read CLAUDE.md, docs/contracts.md, docs/progress.md, and project_plan.md "Emitters, behaviors and waveforms" (IQ-file emitter, Rendering paths), "Receiver model" and "Time, numerics and determinism".

1. Add a continuous-source interface and SigMF ci16_le/cf32_le IQ-file emitter, with native rate/frequency, explicit physical level normalization, start time and stop/zero-pad/loop policy. Define loop-boundary behavior and validate malformed metadata, non-finite values and truncated inputs.
2. Generate/read each source once at its native rate and share immutable windows across paths. On CPU, select and filter its spectrum to each render target before sampling, with the target's H_rx or B_out, then apply the rate conversion/Farrow fractional delay at retarded native-sample coordinates. Retain deterministic filter history/halos across blocks and specify boundary padding and group delay. A common low native rate must still support frequency translation into another channel without aliasing intermediate samples.
3. Extend work-list source/table references and scenario TOML; manifest entries hash the IQ bytes and source metadata. Large files are streamed through bounded caches, with source-relative rebasing for long runs.
4. examples/iq_file.toml and tests: independently generated tone/chirp fixtures verify frequency, power, fractional delay and out-of-band rejection for both taps; split-block processing preserves source/filter continuity; stop/loop policies are deterministic; changed IQ bytes make replay fail. Record the numeric operation order for the CUDA port in prompt 20.

Follow the Definition of Done in CLAUDE.md.
```

## Phase 3: Motion and channel

**Prompt 7: Frames, geodesy and altitude**

```text
Read CLAUDE.md, docs/progress.md, and project_plan.md "World model" (Frames, Mode C altitude). In iqsim-geo, all in f64 using libm:

- WGS84 geodetic <-> ECEF (document the method and convergence), ECEF <-> ECI as a fixed-rate rotation about z from the scenario epoch, local NED at a reference point, body frame from attitude quaternions, antenna frame from mount orientation. One quaternion type with documented conventions (Hamilton; active vs passive), pinned by tests.
- Geoid: load an EGM96 or EGM2008 grid from a file path in config (no download at build time), bilinear interpolation. If no grid is configured, geoid height is 0 and a warning goes into run metadata.
- Register geoid content/metadata hashes in the run manifest, including the explicit zero-geoid fallback when selected.
- Pressure altitude: ellipsoid height -> MSL -> ISA pressure altitude, with an optional local pressure offset. Gillham (Mode C) encoding in 100 ft steps, tested against published Gillham tables.

Tests: published geodetic/ECEF reference points and round trips with separate dimensional tolerances: an initial 1e-6 m absolute ECEF tolerance for terrestrial/LEO cases, scale-aware ULP allowances, and angular/altitude tolerances stated separately. Justify them against the highest carrier's phase budget in docs/contracts.md rather than requiring a universal 1 nm bound. Test poles, dateline and altitude extremes; NED and body rotations against hand-computed cases; changed geoid data rejects replay. Follow the Definition of Done in CLAUDE.md.
```

**Prompt 8: Trajectories and attitude**

```text
Read CLAUDE.md, docs/progress.md, and project_plan.md "World model" (Motion sources, Runtime maneuvers). In iqsim-geo:

- Trajectory trait: state(t) -> position, velocity, acceleration, attitude, angular rate, for any TruthTime; deterministic and pure.
- Maneuver primitives: straight and level, coordinated turn (by rate or bank), climb or descent (by rate or angle), speed change. Jerk-limited transitions keep position C2 continuous.
- Derived attitude as in the design: bank from tan(phi) = V*omega/g, heading from track, pitch from flight-path angle; optional angle of attack and wind crab.
- Route compiler: waypoints with speeds and altitudes -> primitive sequence.
- Track replay: timestamped positions -> C2 smoothing spline, with a documented smoothing parameter.
- Static and ground-clamped placeholders (terrain arrives in prompt 22).
- An append API that adds primitives from the current state at a given time (runtime control uses it in prompt 16).
- Scenario TOML for all of the above.

Tests: C2 continuity by finite differences at every joint; coordinated-turn bank within 0.01 deg of the formula; routes pass through waypoints within tolerance; replay preserves timing. Follow the Definition of Done in CLAUDE.md.
```

**Prompt 9: Light time and per-block polynomial fitting**

```text
Read CLAUDE.md, docs/progress.md, and project_plan.md "Time, numerics and determinism" (Error budget) and "Antennas and propagation" (Light time).

1. For each path, solve |r_rx(t) - r_tx(t - tau)| = c*tau in ECI by fixed-point iteration with a documented stopping rule. Use platform origins for now; prompt 10 switches to antenna phase centers.
2. Evaluate geometry at ClockMap acquisition instants with their retained fractional time. Fit delay, retarded waveform-coordinate and phase polynomials (order <= 3, documented nodes such as Chebyshev), then convert to fixed-point BlockPoly. Include propagation phase -f_c*tau and transmitter/receiver carrier offsets evaluated at actual acquisition time. Use the modular phase routine from prompt 2; avoid large floating-point absolute-time phase products.
3. Check endpoints and interior residuals against exact solutions, with default bounds 0.001 samples and 0.01 rad, plus waveform-coordinate error under the equivalent timing bound. Define derivative bounds or conservative adaptive refinement between checks. Add amplitude/complex-gain checks as prompt 10 supplies them. Split first at trajectory knots and known discontinuities, mapped from emit time to receive time per path, then refine smooth intervals as needed. Log every split with block, path and residual; fail explicitly at the refinement limit.
4. Replace the static-geometry builder from prompt 5 with this geometry engine.

Tests: a moving CW emitter's measured frequency and phase trajectory match the analytic solution within the bounds; a LEO-like fast mover forces splits and still meets the bounds; hashes are stable across thread counts. Follow the Definition of Done in CLAUDE.md.
```

**Prompt 10: Antennas, mounts and airframe blockage**

```text
Read CLAUDE.md, docs/progress.md, and project_plan.md "Antennas and propagation" (Pattern representation through When gains are evaluated) and "World model" (Platforms and mounts). In iqsim-antenna:

- PatternGrid: complex gain on a theta/phi grid with polarization as a Jones vector, bilinear interpolation, documented frame conventions.
- Parametric library: isotropic; dipole; monopole and blade over a finite ground plane (document the approximation); SSR large-vertical-aperture antenna with sum, difference and control patterns (parameters: azimuth beamwidth, elevation shaping, control-pattern fill); parabolic dish from diameter and efficiency; circularly polarized patch; arrays as element pattern x array factor.
- Mounts: lever arm and orientation in the body frame. The light-time solver now uses antenna phase centers.
- Polarization mismatch from TX and RX Jones vectors rotated into a common frame.
- TX gain evaluated at emit time, RX gain at receive time; antenna phase goes into the phase polynomial.
- Airframe blockage: a fuselage cylinder plus wing and tail plates per platform. A path whose ray from the antenna hits a primitive gets shadow loss (document the loss model). Gains, losses and mismatch fold into the amplitude polynomial.
- Extend prompt 9's refinement to the amplitude and near-null complex-gain bounds in docs/contracts.md. Detect blockage/pattern discontinuity crossings, refine their times and create explicit segment boundaries; never smooth across a shadow step or demand phase at a zero. Carry future effective parameter-change boundaries through the same API. Register external pattern/airframe assets with the manifest.
- Scenario TOML for mounts, patterns and airframes.

Tests: each parametric pattern integrates to its expected directivity within tolerance; linear-to-circular mismatch is 3 dB; a belly antenna loses sight of a ground site on the outside of a banked turn; a lever arm shifts delay by the expected amount; a null and a blockage crossing inside a block meet the gain budgets with explicit boundaries. Follow the Definition of Done in CLAUDE.md.
```

## Phase 4: Mode A/C end to end, offline

**Prompt 11: Behavior engine and truth log**

```text
Read CLAUDE.md, docs/progress.md, and project_plan.md "Emitters, behaviors and waveforms" (Behavior engine) and "Interfaces" (Exports). In iqsim-behavior:

- A discrete-event engine: events ordered by (TruthTime, entity id, sequence number). Handlers are pure functions of state and event; they may schedule new events and produce Emission records. No wall-clock access.
- A link query API that behaviors use to ask what a given antenna on another platform would see from an emission: arrival time, received power and per-component gains, using the geometry and antenna code. Link-budget level, not IQ.
- Lookahead: the engine hands each block every emission that can reach any receiver within it, accounting for maximum path delay, emission duration and filter support. It may run ahead of rendering only inside the overall latency budget, reserving its consumed dependency interval at the shared exclusive commit boundary C before advancing. State before C is immutable. Future commands at a boundary are processed before advancing across it; future scheduled events are not themselves processed state. Keep rendering progress distinct from C. Provide an atomic admission API for prompt 16; v1 uses no rollback.
- Truth log records: emissions, causal links (the event that caused each emission), and entity states at a configurable rate.
- Wire emissions into the work-list builder.

Tests: ties in event ordering resolve deterministically; lookahead never misses an emission (adversarial test with long delays); an admitted boundary event cannot be overtaken by behavior advancement; truth hashes are identical across runs. Follow the Definition of Done in CLAUDE.md.
```

**Prompt 12: SSR interrogator and Mode A/C transponder**

```text
Read CLAUDE.md, docs/progress.md, and project_plan.md "Emitters, behaviors and waveforms" (SSR interrogator, Mode A/C transponder, Rendering paths). Implement in iqsim-behavior and iqsim-waveform, and cite the source of every timing constant in a comment:

- SSR interrogator: rotating antenna using the sum and control patterns from prompt 10; configurable scan period and start azimuth; PRF with stagger or jitter patterns drawn from Philox; Mode A/C interlace pattern; P1/P2/P3 pulses 0.8 us wide, P1-P3 spacing 8 us (Mode A) or 21 us (Mode C), P2 2 us after P1 and radiated on the control pattern.
- Mode A/C transponder: reception decided at link-budget level through the link query (P1/P3 power vs sensitivity; P2 vs P1 for sidelobe suppression). Reply 3.0 +/- 0.5 us after P3 with jitter; reply pulses 0.45 us wide on a 1.45 us grid; framing pulses 20.3 us apart; the 12 code pulses in standard positions; SPI. Per-unit carrier offset within 1090 +/- 3 MHz; dead time and suppression windows. Mode A code from the scenario (settable at runtime later); Mode C from the Gillham encoder in prompt 7.
- Pulse shapes with parameterized rise and fall times.
- Scenario TOML for both, with defaults matching typical equipment.

Tests: pulse positions and widths match the standard timing; a transponder in the interrogator's sidelobes does not reply; replies carry the right code bits; dead time is enforced; all randomness is reproducible. Follow the Definition of Done in CLAUDE.md.
```

**Prompt 13: Reference scenario end to end, offline**

```text
Read CLAUDE.md, docs/progress.md, and project_plan.md "Purpose and scope" (Reference scenario), "Interfaces" (Exports) and "Validation strategy".

1. examples/reference_modeac.toml: 3 aircraft with Mode A/C transponders on maneuvering routes; 4 fixed SSR interrogators with different scan periods and PRFs; 1 fixed receiver with two channels at 1030 and 1090 MHz.
2. Exports: the truth log as Parquet (entity states, emissions, causal links, receiver sample indices) and the consumer-knowledge JSON (interrogator positions with optional survey error, nominal PRI and scan schedules, receiver site). Truth never appears in SigMF metadata.
3. An independent Mode A/C decoder under tests/tools that does not import iqsim-waveform or iqsim-behavior: pulse detection on IQ, framing-pair search, code extraction. Decode the reference run with it.
4. An independent baseline bistatic-ellipse estimator under tests/tools, importing neither geometry/light-time code nor behavior implementations. Its inputs are the two receiver IQ streams, ordinary receiver metadata and the consumer-knowledge JSON. Associate received interrogations/replies using observable timing and codes; never use truth identifiers, true stagger schedules or causal links to choose associations. Declare the altitude/ambiguity assumptions and choose a well-conditioned test geometry. Truth is used only after estimation for scoring.
5. Tests: decoded squawks and altitudes match truth; reply-delay statistics match the model; noiseless geolocation and seeded timing/noise cases meet an explicit, independently derived position-error budget and report ambiguity/missed fixes; two runs give identical sample and canonical truth hashes; `sim run --offline` and manifest replay complete the reference scenario. Add impaired-clock geolocation cases in prompt 15.
6. Record reference duration, channel rates/bandwidths, render time against simulated time, estimator error statistics and acceptance thresholds in docs/progress.md. Add all resolved inputs to the manifest.

Follow the Definition of Done in CLAUDE.md.
```

## Phase 5: Real time and control

**Prompt 14: Real-time pipeline and pacer**

```text
Read CLAUDE.md, docs/progress.md, and project_plan.md "Time, numerics and determinism" (Blocks and latency) and "Architecture".

In iqsim-runtime, build the block pipeline: behavior and geometry stages run ahead, then render, global path reduction, tap-specific finalization/receiver chain and sinks. Stages connect through bounded queues sized from the latency budget. Offline mode runs the same pipeline, unpaced. Stateful receiver stages have one ordered owner per channel, with explicit history handoff.

- Real-time mode: a wall-clock pacer releases blocks against the configured output timebase, and each block must be ready before its release time. A late block is a hard error with a diagnostic (which stage, by how much), never a silent drop or zero-fill. Required sinks have bounded delivery deadlines; optional preview queues cannot stall them.
- Coordinate the shared exclusive commit boundary C from prompt 11 across behavior, geometry, work-list preparation and receiver/filter lookahead. A stage reserves its dependencies before consuming them. Atomic command admission at max(requested time, C) inserts the event before any further advancement across that time. All lookahead counts toward the overall latency budget; no rollback in v1.
- Telemetry counters: render margin per block, queue depths, stage timings, commit boundary, admission latency and required-sink delivery margin. Keep these nondeterministic measurements outside canonical truth hashes.
- A pacing abstraction, so a hardware sink (prompt 24) can replace the wall clock as the pacer.
- Config: block length, latency budget, thread counts.

Tests: offline and real-time modes give identical sample hashes for the same scenario; an artificially slowed stage triggers the underrun error; injected test commands cannot change a dependency already committed by a later stage; optional preview stalls do not block required output; the reference scenario runs in real time on this machine, with its margin recorded in docs/progress.md. Follow the Definition of Done in CLAUDE.md.
```

**Prompt 15: Receiver clock and streaming outputs**

```text
Read CLAUDE.md, docs/progress.md, and project_plan.md "Receiver model" (Clock model) and "Interfaces" (IQ outputs).

1. Implement distinct sampling-oscillator, timestamp-clock and LO-error models, with explicit coupling in receiver profiles. Extend ClockMap with offset, frequency error, drift/aging and an optional documented Philox-driven GPSDO process. Require monotonicity and a stable inverse, preserve fractional time, and evaluate geometry, retarded waveform coordinates and carrier phase at actual acquisition instants. Generate stream timestamps/SigMF datetimes from the reported clock. Shared-channel errors are evaluated at common physical times even at different rates. Truth records acquisition time, reported time and sample index separately; antenna-port rendering uses its independent output-device timebase.
2. Strict framed ZMQ IQ sink over TCP or IPC: ROUTER/DEALER sessions with ROUTER_MANDATORY, bounded queues and explicit per-consumer routing/fan-out. Define a versioned header (stream/session id, block sequence, first sample index, count, datatype/endianness, reported first-sample timestamp, context version) plus binary IQ payload and RF context. Include no truth fields. READY handshake, acknowledgements of contiguous received blocks, bounded outstanding-byte/ack windows, connection-loss handling and a final acknowledged end-of-stream boundary. Required-sink deadline expiry fails a real-time run; successful enqueue alone never proves delivery. A failed session cannot silently reconnect/resume as a successful stream. Prefer the zmq crate with vendored libzmq; verify an offline build.
3. VITA-49 UDP sink: implement and document a specific VITA 49.2 profile, signal-data stream IDs/packet counts, integer/fractional reported timestamps, context fields/cadence and MTU. Its receiver detects gaps/reordering/duplicates using packet count, timestamp and payload length together, including count wrap. UDP has no receiver delivery guarantee; document that sender-visible end-to-end confirmation requires the strict framed sink. Validate interoperability with an independent parser or an authoritative packet fixture, as well as local round trips.
4. Raw ZMQ PUB preview: sc16/fc32 items directly readable by GNU Radio ZMQ SUB Source, one socket per channel. Explicitly best effort under startup/overflow, with no claim of gap detection or strict TDOA timing. Keep its queues independent of required sinks. Document that per-subscriber drops are not all observable at the publisher.
5. Reference strict and VITA receivers plus a GNU Radio adapter for framed streams. Publish docs/transport.md with the alternative/selection rationale from the design, wire formats, startup, failure and shutdown semantics. Raw preview uses the stock SUB Source; framed metadata needs the adapter. Preserve this decision behind replaceable sink interfaces.
6. Measure strict framed ZMQ and VITA-49 throughput/latency on this machine at the reference rates, including slow-consumer behavior, memory/copy cost and ACK overhead. Record results and required headroom in docs/progress.md; if the reference workload cannot meet its latency budget, resolve the sink design or configured workload before declaring phase 5 complete. Prompt 21 repeats capacity measurement on the target host.

Tests: two receivers observing analytic CW/pulses show the predicted sampling-error shift in IQ and a separately configured timestamp offset in metadata; coherent channels at differing rates share the correct error process; run impaired-clock geolocation cases from prompt 13. Strict streams hash exactly after handshake and either finish with final acknowledgement or fail on slow/missing consumers, disconnects and lost ACKs; use protocol synchronization instead of arbitrary startup sleeps. UDP tests inject loss/reordering/duplicates and count wrap; preview late-join/overflow is allowed and never stalls required output. Re-run antenna-port isolation checks with clock/LO impairments enabled. Follow the Definition of Done in CLAUDE.md.
```

**Prompt 16: Control plane, event log and replay**

```text
Read CLAUDE.md, docs/progress.md, and project_plan.md "Interfaces" (Control plane) and "World model" (Runtime maneuvers). In iqsim-control:

- Transport-independent command handlers and JSON messages in a versioned envelope {schema_version, msg_type, client_session_id, id, payload}, defined as Rust types. Generate docs/control.schema.json with schemars through xtask and check it in; it is the contract the existing range tool will adapt to and the later HTTP/WebSocket adapter will reuse.
- A ROUTER socket for DEALER command clients, with ROUTER_MANDATORY, bounded queues, send deadlines and acks carrying the effective time and ingress sequence after durable event-log append. Deduplicate retries by client/session and request id and return the original outcome. Add a status/snapshot query. PUB telemetry (entity states, link status, render/commit/delivery margins, underruns) is best effort, includes sequence numbers and supports snapshot-based resynchronization.
- Commands: maneuver (through the trajectory append API), add and remove entity, enable and disable emitter, set parameters (squawk, power, frequency), noise or jamming injection, start, stop.
- Effective time = max(requested truth time, shared exclusive commit boundary C), using the atomic admission API from prompts 11/14. A stage cannot advance across that time before insertion. No rollback; preserve already emitted/in-flight signals. Prepare frequency-dependent response tables and other command resources before accepting their effective boundary. The append-only deterministic event log records request identity, requested/effective times, ingress sequence, command and acceptance/rejection outcome, plus start/stop/final-output boundaries. Define how future pending commands are handled at stop.
- Extend `sim replay <run-manifest>` to verify all inputs/event-log content and replay accepted events at their logged effective times and ordering; never recompute effective times from replay speed. Canonical truth hashes exclude wall-clock telemetry. Register the finalized event-log hash with the manifest and report incomplete runs explicitly.

Tests: a scripted client sends commands during a real-time run, and replay matches sample/canonical-truth hashes. Exercise behavior already ahead of rendering, concurrent admission at a boundary, long-delay in-flight emissions, equal-time commands, a lost ACK followed by retry, runtime frequency/table changes and stop with pending commands. No command executes twice or changes committed history; malformed requests get clear errors. Changed/missing input assets or event-log bytes reject replay before rendering. Follow the Definition of Done in CLAUDE.md.
```

## Phase 6: Authoring and packaging

**Prompt 17: Libraries and templates**

```text
Read CLAUDE.md, docs/progress.md, and project_plan.md "Interfaces" (Scenario files) and "World model" (Platforms and mounts).

- Library system: platform templates (airframe primitives, mounts, default emitters), antennas, waveform presets, receiver profiles and interrogator presets. Each is a TOML file found on the library search paths, referenced by name, with parameter overrides in the scenario. Document and test resolution order and conflict rules.
- Save fully resolved configuration and hash every resolved library asset in the manifest. Replay must not silently pick a changed search-path item; test changed-content rejection and replay from an explicit matching asset bundle.
- Ship a starter library: generic transport and business-jet airframes with top and bottom L-band blades; generic terminal and en-route SSR presets; a two-channel 1030/1090 MHz receiver profile; a wideband receiver profile.
- `sim new --template <name> <dir>` creates a scenario from a template.
- Convert examples/reference_modeac.toml to use the library. Its output hash must not change.

Follow the Definition of Done in CLAUDE.md.
```

**Prompt 18: Full lint and dry-run**

```text
Read CLAUDE.md, docs/progress.md, and project_plan.md "Interfaces" (CLI, IQ outputs).

- `sim lint`: everything from prompt 4, plus library resolution, unattached antennas, emitters on platforms with no mount, overlapping receiver channels, and plausibility warnings (for example a 1090 MHz transponder at 10 W).
- `sim dry-run`: no IQ rendering. Run behavior and geometry over the scenario and report per-link visibility over time, predicted SNR per link, emitters outside every receiver's band, terrain coverage (stub until prompt 22), polynomial split counts, and estimated render cost per second of scenario. Terminal summary plus --json.
- Output/resource budget: map every sink/device to named shared NICs, disks, USB controllers and PCIe links/roots with usable capacity and headroom. Sum concurrent sample traffic, framing/packet/ACK overhead and fan-out copies along each shared resource, plus the two-i64-per-complex-sample partition transfers introduced in prompt 21. Show per-sink rates and aggregate resource totals. Two 20 MS/s sc16 streams sharing 1 GbE need at least 1.28 Gb/s and fail even though each fits alone. Reject known real-time overloads; offline runs warn where pacing can avoid them. Report unknown topology/capacity as unverified, with an actionable measurement/configuration step.
- Every lint and dry-run message is actionable: what is wrong, where (file:line), and how to fix it.

Tests: golden text output for the reference scenario and deliberately broken scenarios, including an over-budget sink, individually valid sinks sharing an overloaded NIC/disk, explicit fan-out copies and unknown capacity. Extend PCIe cases when prompt 21 lands. Follow the Definition of Done in CLAUDE.md.
```

**Prompt 19: Packaging for air-gapped RHEL 9**

```text
Read CLAUDE.md, docs/progress.md, and project_plan.md "Requirements" (Easy installation).

- `cargo xtask package`: a release build, then an RPM (spec file or cargo-generate-rpm) and a tarball containing the binary, starter library, schemas, examples, transport reference clients/GNU Radio adapter and docs. Include a way to bundle a manifest with its permitted input assets for replay on an air-gapped machine.
- Vendored build: a `cargo vendor` flow, a .cargo/config.toml template that points at either vendored sources or an internal registry, and a script that verifies the build with networking disabled.
- The default build needs no CUDA. The CUDA backend sits behind a feature and loads the driver at runtime.
- docs/install.md and docs/quickstart.md: install, first lint, first dry-run, first offline run and verified replay, first strict real-time stream into a GNU Radio flowgraph through the supplied adapter, and optional raw preview using SUB Source. Explain delivery guarantees and enumerate native build/runtime prerequisites for vendored libzmq; later CUDA/UHD features extend packaging verification.
- If podman is available, add an xtask that installs the RPM in a clean RHEL 9-compatible container and runs the quickstart.

Follow the Definition of Done in CLAUDE.md.
```

## Phase 7: CUDA and multi-GPU

**Prompt 20: CUDA backend**

```text
Read CLAUDE.md, docs/progress.md, and project_plan.md "Compute backends and scaling", "Time, numerics and determinism" and "Receiver model" (Band-limiting before sampling). Implement iqsim-render-cuda behind a cargo feature.

- build.rs: when the feature is on, compile the kernels with nvcc into a fatbin for sm_75 and sm_86 with -fmad=false and no fast-math flags, and embed it with include_bytes!. If nvcc is missing, fail with a clear message.
- Host side on cudarc with dynamic loading. At startup, if no driver or device is found, fall back to the CPU backend and log it.
- Extend prompt 19's air-gapped package checks/documentation for CUDA-enabled artifacts: pin compatible toolkit/driver versions, include the fatbins and distinguish build-time nvcc from the runtime driver prerequisite. Verify the same shipped CUDA-enabled binary falls back on a machine without a CUDA driver.
- Kernels port the iqsim-core arithmetic line by line under docs/contracts.md: Newton-form delay, retarded waveform-coordinate and phase polynomials with the same explicit 64-bit primitives (each thread jumps once to its chunk start, then steps with 128-bit adds); the sincos table in constant memory; band-limited edge tables in global memory; the CPU IQ-file filter/resampler from prompt 6a; Philox, ln and Box-Muller. Each render thread owns a chunk of output samples and accumulates path contributions into i64; no atomics and no receiver noise/effects in partial-sum output.
- Expose separate raw-sum and finalization operations. A complete channel may finalize on its owning GPU and convert to sc16 there. An emitter partition always returns raw i64 I/Q; global reduction precedes exactly one noise/receiver-chain pass and one final output conversion. The designated finalizer owns filter/AGC/LO history in channel order; apply its history or deterministic halos rather than resetting per CUDA chunk. Antenna-port finalization excludes virtual effects.
- `cargo xtask gpu-test`: thousands of seeded random work lists rendered on CPU and GPU must hash identically, including raw sums and final sc16/fc32 outputs, IQ-file paths, both taps, rounding/subnormal boundaries, phase wrap, clock errors and GPU chunks differing from CPU chunks. Test finalization across block boundaries with noise shaping/decimation enabled; reject unsafe accumulation bounds.

GPU tests skip with a clear message when no GPU is present. Record CPU vs GPU throughput on this machine in docs/progress.md. Follow the Definition of Done in CLAUDE.md.
```

**Prompt 21: Multi-GPU partitioning and benchmark**

```text
Read CLAUDE.md, docs/progress.md, and project_plan.md "Compute backends and scaling".

- Partitioner: assigns complete receiver/output channels to devices first. Split one channel's emitters only when capacity and shared-link budgets justify it; each partition returns raw i64 I/Q with no noise, receiver effects, clipping or sc16 conversion. Reduce exactly on the host, then finalize once on the CPU or one designated device, preserving per-channel state. A device can be a GPU or the CPU backend.
- Heterogeneous test: GPU plus CPU and two CUDA streams on one GPU match the raw sums and final output of an unsplit run. Enable noise, IQ imbalance, shaping/decimation and clipping (AGC too when implemented) so duplicate effects or early conversion cannot pass unnoticed. Cover both taps and multiple consecutive blocks.
- Memory: size buffers from each device's reported memory, with 4 GB T1000s as the minimum target.
- `sim bench`: synthetic dense and sparse workloads; reports path-samples per second per device, PCIe transfer rates and real-time margin; writes a JSON report.
- Extend prompt 18's resource graph with actual topology/traffic: full-channel sc16 at 100 MS/s is 400 MB/s, whereas each i64 I/Q partition transfers 1.6 GB/s, plus return traffic if GPU finalization is used. Include shared roots, output copies and the active strict/VITA/preview sinks. Reject a known real-time aggregate overload. Benchmark each supported sink at the reference and target rates; record bytes, copy/ACK overhead and observed delivery margin instead of assuming ZMQ or UDP is faster.
- docs/performance.md: how to run the benchmark on the 4x T1000 host and what to record.

Follow the Definition of Done in CLAUDE.md.
```

## Phase 8: Terrain

**Prompt 22: Terrain data**

```text
Read CLAUDE.md, docs/progress.md, and project_plan.md "Antennas and propagation" (Terrain). In iqsim-terrain:

- An in-house DTED reader for Levels 0, 1 and 2: UHL, DSI and ACC headers, data records, checksums, latitude-dependent longitude spacing. Validate against the DTED specification.
- A GeoTIFF elevation reader (a pure-Rust TIFF crate is fine; parse the GeoKeys you need in-house).
- A tile index over a directory, a memory-mapped tile cache with a size limit, bilinear elevation queries, and terrain-profile extraction between two points with documented sampling.
- Run metadata records the terrain source, level and resolution. Register all selected tile bytes and relevant CRS/vertical-datum metadata in the run manifest; hash before first use and require matching assets on replay. Record the explicit smooth-Earth fallback when no terrain is loaded.
- Ground-clamped trajectories from prompt 8 now use terrain, and `sim dry-run` reports terrain coverage.

Tests use small synthetic DTED and GeoTIFF files that the test suite generates. Follow the Definition of Done in CLAUDE.md.
```

**Prompt 23: Propagation paths**

```text
Read CLAUDE.md, docs/progress.md, and project_plan.md "Antennas and propagation" (Path set per link).

- Write the reference first: a Python script under tests/reference, written from the ITU-R P.526 text before any Rust code, that generates expected knife-edge diffraction losses for a set of test profiles.
- Line of sight over the terrain profile with Earth curvature and a refraction k-factor (default 4/3).
- Diffraction: P.526 multiple knife-edge (state which method you implement) as an amplitude loss on the direct path.
- Ground reflection: find the specular point on the terrain; compute its own delay and phase, Fresnel reflection coefficients by polarization with ground-constant presets from the ITU-R P.527 categories, and a roughness reduction factor. It becomes its own path in the work list.
- Optional Longley-Rice (ITM) excess loss for beyond-line-of-sight links, using the NTIA C++ implementation behind a feature, or a documented hook if you defer it.
- Path updates run at the geometry rate and flow through the block polynomials. Path appearance/disappearance, reflection topology changes and visibility/blockage crossings create explicit segment boundaries through prompts 9/10's refinement API, with amplitude/complex-gain tolerances applied on smooth intervals. Register external ground-constant tables in the manifest.

Tests: flat-earth two-ray interference against the closed form; diffraction within 0.1 dB of the Python reference; line-of-sight masking behind a ridge. Follow the Definition of Done in CLAUDE.md.
```

## Phase 9: Hardware in the loop

**Prompt 24: SDR TX through uhd-rs**

```text
Read CLAUDE.md, docs/progress.md, and project_plan.md "Hardware in the loop" (including Device limits) and "Receiver model". In iqsim-hil:

- An SDR TX sink on uhd-rs (libuhd FFI) for the X410 and B205mini, one or more channels: render the antenna-port target using its B_out, output-device timebase and level normalization. Exclude virtual H_rx, gain, clock/LO impairments, noise and ADC effects. Use timed start with UHD time specs and shared 10 MHz/PPS references when configured. Verify actual negotiated device rates against the scenario before starting; do not silently accept rate coercion.
- The device becomes the pacer through UHD backpressure, using the pacing abstraction from prompt 14. Underflows reported by UHD are hard errors.
- Device capabilities: a built-in table per supported radio (channels, tuning range, bandwidth and sample-rate limits, host link), with values cited from the vendor specifications listed in project_plan.md. `sim lint` rejects over-limit mappings and known aggregate real-time link overloads using prompt 18's shared-resource budget (for example four X410 channels at 400 MHz without sufficient host links). Unknown capacity is unverified. Warn when several radios share no 10 MHz and PPS reference.
- Calibration: a table per device serial, channel and frequency that maps digital full scale to dBm at the port, plus per-path losses (matrix ports, cables, external converters). `sim calibrate` steps frequency and level and records values entered from an analyzer or power meter.
- Record calibration bytes, mapping, device-rate/timebase configuration and relevant external converter/loss assets in the manifest. Extend the air-gapped package documentation/checks with pinned UHD runtime and device image prerequisites; CUDA packaging checks from prompt 20 must also remain current.
- Scenario TOML: a hil.tx section mapping simulated receivers or channels to devices and ports.

Write a manual test plan in docs/hil.md: a B205mini transmitting into a second SDR through attenuators; verify power with an analyzer; decode the reference scenario's replies from the captured RF with the decoder from prompt 13. Tests: lint accepts valid mappings and rejects over-limit/aggregate-overload cases; rate coercion fails clearly; antenna-port IQ remains identical when virtual receiver H_rx/gain/noise/clock settings change. Follow the Definition of Done in CLAUDE.md.
```

**Prompt 25: SDR RX injection**

```text
Read CLAUDE.md, docs/progress.md, and project_plan.md "Hardware in the loop" (In: SDR RX -> simulation).

- An SDR RX source on uhd-rs that captures IQ from a device and exposes it as a scenario emitter placed on a platform. Reuse the continuous-source/filter/resampler path from prompt 6a and its CUDA port from prompt 20: the selected H_rx or B_out acts before sampling, with retarded source coordinates and fractional delay. Map source device time to the simulator's truth timeline explicitly, while keeping output-reported timestamps separate.
- Record the exact IQ samples accepted into simulation as SigMF, with source rate/frequency/level metadata, device timestamps, timeline mapping and gap/overrun/stop decisions. Add finalized bytes and timing records to the manifest. Replay uses the same source interface and recorded decisions; it verifies all hashes before rendering and never substitutes recaptured hardware input.
- The output runs a fixed latency behind real time; capture jitter beyond the latency budget is a hard error.

Tests: replaying a recorded live-input run matches the original sample/canonical-truth hashes; changed recording bytes or timestamp/gap records reject replay; with no hardware, a file-backed fake device exercises timing alignment, capture jitter and the same error decisions. Follow the Definition of Done in CLAUDE.md.
```

## Phases 10–12: Waveforms, satellite and editor

These prompts open with a planning step and span several sessions; refine them once phases 1–9 have landed.

**Prompt 26: Waveform expansion**

```text
Read CLAUDE.md, docs/progress.md, and project_plan.md "Emitters, behaviors and waveforms" (Library roadmap). Work one family per session, in this order, each with an independent decoder or demodulator test. Before coding, propose a sub-plan for this session's family and wait for approval.

1. Extend the continuous-waveform path from prompts 6a/20 for generated waveforms: select/filter the source spectrum with the selected tap's H_rx or B_out in source-frequency coordinates, then perform output-frame frequency translation and retarded-coordinate rate conversion/Farrow delay without an aliased intermediate. Generate each emitter once at its native rate and share source windows across paths. Preserve bit-exact CPU/CUDA behavior, filter history, clock mapping and the two tap contracts; extend the existing implementation rather than adding a second resampler.
2. Mode S and ADS-B: 1090ES extended squitters with CRC, and Mode S replies; DF17 position and velocity encoded from truth.
3. DME interrogator and transponder pulse pairs.
4. Generic pulsed radar: PRI patterns, pulse-compression codes, scan.
5. AM and FM; generic PSK and QAM with root-raised-cosine shaping.
```

**Prompt 27: Satellite**

```text
Read CLAUDE.md, docs/progress.md, and project_plan.md "Antennas and propagation" (Light time, Relay nodes) and "Emitters, behaviors and waveforms" (Library roadmap). Plan first and wait for approval; split into sessions as needed.

1. SGP4 trajectories from TLE and OMM (in-house or a vendorable crate), converted from TEME into the simulator's frames.
2. Relay nodes: the link graph becomes a per-block DAG. Bent-pipe transponder with receive filter, frequency translation, amplifier nonlinearity (Saleh or similar), and uplink noise and interference carried through.
3. DVB-S2-structured carriers: PL framing, pilots, symbol rate, rolloff, random payload. Tests detect the PL headers with an independent detector.
4. Dense-rendering performance on the T1000 host, compared with the capacity estimate in project_plan.md.
```

**Prompt 28: Map-based scenario editor**

```text
Read CLAUDE.md, docs/progress.md, project_plan.md and docs/control.schema.json. Plan first and wait for approval.

Build a web scenario editor served by `sim serve`. It works fully offline: a DTED shaded-relief basemap (no tile server); place and edit platforms; draw routes; attach library templates; save TOML that passes `sim lint`; run dry-run and show its report; start a run and show live telemetry. Add an HTTP/JSON and WebSocket adapter over the transport-independent handlers from prompt 16; browser clients reuse command IDs, acknowledgements, commit-boundary admission and status/snapshot semantics. They do not connect to native ZMQ sockets or create a second command implementation. Keep operator truth telemetry separate from consumer-IQ endpoints and label raw spectrum previews as best effort.

Acceptance: the reference scenario can be built and run without editing files by hand.
```
