# Claude Code prompts

Twenty-nine prompts, run in order, take the design from an empty repository to phase 12. Prompts 0–25 are ready to run; 26–28 plan first and get refined once earlier phases land.

## How to run these

- Export the design doc's main tab as Markdown and commit it as `docs/design.md`. Every prompt cites its sections by heading.
- Run one prompt per fresh Claude Code session (or `/clear` between prompts), in order. Start each in plan mode and approve the plan before any code is written.
- Prompt 0 writes `CLAUDE.md` with the project invariants and a Definition of Done. Later prompts point to it instead of repeating it.
- Each prompt ends by appending to `docs/progress.md`: what was built, deviations from the design, open questions. The next session reads it first.
- After each prompt: review the diff, run `cargo xtask ci`, commit. Golden hashes change only through `cargo xtask bless`, and only on purpose.
- `iqsim` is a working name for the workspace and `sim` for the binary; rename them before prompt 0 if you like.
- Prompts 20–21 need an NVIDIA GPU and 24–25 need an SDR. The rest run on any Linux machine.

## Setup

**Prompt 0: Repository bootstrap**

```text
Read docs/design.md end to end. You are starting a new Rust workspace, "iqsim", for the real-time scenario IQ simulator it describes. Do not implement simulation logic yet.

1. Create a Cargo workspace (stable toolchain pinned in rust-toolchain.toml) with these library crates, each with a crate-level doc stating its responsibility from the design:
   iqsim-core (time, fixed-point, phase, sin/cos, RNG), iqsim-units, iqsim-scenario, iqsim-geo, iqsim-antenna, iqsim-terrain, iqsim-behavior, iqsim-waveform, iqsim-worklist, iqsim-render (backend trait + CPU backend), iqsim-render-cuda, iqsim-receiver, iqsim-io (sinks and exports), iqsim-control, iqsim-runtime, iqsim-hil.
   Add a binary crate iqsim-cli that builds `sim`, and an xtask crate.
2. xtask commands: `ci` (fmt check, clippy with -D warnings, all tests), `bless` (regenerate golden hashes), `gpu-test` (placeholder for now).
3. Golden-test helper: assert_golden(name, bytes) compares a SHA-256 against tests/golden/<name>.sha256 and fails with a clear message; `cargo xtask bless` rewrites them.
4. Write CLAUDE.md containing:
   - Project summary; read docs/design.md and docs/progress.md at the start of every session.
   - Invariants:
     * No floating point for time, delay or carrier phase on the render path.
     * Integer (i64) accumulation into receiver buffers.
     * No std f64 transcendentals (sin, cos, exp, ln, atan2, powf) in control-plane code; use the libm crate so results match across machines and libc versions.
     * No fast-math, no FMA contraction, no nightly features.
     * #![forbid(unsafe_code)] everywhere except iqsim-render-cuda and iqsim-hil.
     * Deterministic iteration: no HashMap/HashSet iteration order in outputs; use BTreeMap or sorted Vecs.
     * Every random draw comes from the Philox streams in iqsim-core.
   - Dependency policy: few, well-maintained crates that vendor cleanly for an air-gapped RHEL 9 build; justify each new dependency in your summary; no network access at build time.
   - Definition of Done: `cargo xtask ci` passes; new behavior has unit tests and, where docs/design.md names a validation check, a test for that check; public items documented; docs/progress.md appended with what was built, deviations from the design and why, and open questions.
5. Create docs/progress.md with a "Prompt 0" entry.

Stop when the workspace builds, `cargo xtask ci` passes and CLAUDE.md is complete. Summarize the layout and any decisions you made.
```

## Phase 1: Foundations

**Prompt 1: Time and rate types**

```text
Read CLAUDE.md, docs/progress.md, and docs/design.md section "Time, numerics and determinism". Implement in iqsim-core:

- TruthTime: i64 picoseconds with checked arithmetic. Overflow is an error, never a wrap.
- SampleRate: an exact rational (u64 numerator and denominator, reduced). Constructors from integer Hz and from rationals; reject zero.
- SampleIndex: u64 per receiver channel. Exact conversions between SampleIndex and TruthTime given a channel epoch and SampleRate, with documented rounding (floor for time to index; exact for index to time when representable, otherwise floor plus the returned remainder).
- Block: a half-open range of sample indices for one channel; BlockSpec (block length in samples) and an iterator over blocks.
- Use u128/i128 intermediates where needed, and document why no intermediate can overflow for rates up to 1 GS/s and runs up to 100 days.

Tests: property tests (proptest) for round trips, monotonicity and block tiling (no gaps, no overlaps), plus edge cases at 61.44 MS/s and 250 MS/s. Follow the Definition of Done in CLAUDE.md.
```

**Prompt 2: Fixed-point delay, carrier phase, polynomials and sin/cos**

```text
Read CLAUDE.md, docs/progress.md, and docs/design.md "Time, numerics and determinism". In iqsim-core, implement the numeric core every render backend shares. The CUDA backend (prompt 20) will port this code line by line, so write it with that in mind.

- Delay output: Q32.32 fixed-point samples in i64.
- Phase output: u64 fraction of a cycle (2^64 = one cycle), wrapping arithmetic.
- BlockPoly: a polynomial of order <= 3 in local sample index k (0..block_len), stored in Newton forward-difference form (D0..D3) as 128-bit fixed point, with 64 extra fractional bits beyond the outputs above.
  * Jump: p(k) = D0 + C(k,1)*D1 + C(k,2)*D2 + C(k,3)*D3 in wrapping integer arithmetic, exact for any k.
  * Step: advance one sample with three 128-bit adds.
  * Document the maximum block length (C(k,3) must fit in u64), and show that coefficient quantization error stays far below 0.001 samples of delay and 0.01 rad of phase over a full block.
- Implement the 128-bit arithmetic with explicit 64-bit primitives (mul_hi_u64, wrapping mul, add with carry) in one small module rather than native i128, so the CUDA C version can mirror it exactly. Test it against an i128 reference.
- sincos(phase) -> (f32, f32): the top phase bits index a table built by a documented deterministic routine; a fixed-order polynomial corrects the remainder. Write every operation explicitly (separate mul and add). Report max absolute error over a sweep of the top 24 bits with random low bits; target <= 3e-7.

Tests: jumping and stepping agree bit for bit over full blocks, including phase wrap and negative delay rates; golden hash over a fixed set of 1e6 evaluations. Follow the Definition of Done in CLAUDE.md.
```

**Prompt 3: Philox and Gaussian noise**

```text
Read CLAUDE.md, docs/progress.md, and docs/design.md "Time, numerics and determinism" (Noise). In iqsim-core:

- Implement Philox4x32-10 as specified in the Random123 paper. Take the Philox4x32-10 known-answer vectors from the Random123 distribution's kat_vectors file and include them as tests. Do not invent expected values; if you cannot obtain the file, stop and ask.
- Stream keying: key derived from (scenario seed: u64, entity id: u32, stream id: u32) with a documented, deterministic mixing function; counter = (sample index, sub-index).
- Gaussian: Box-Muller where the angle is a u32 treated as a phase and passed to the sincos from prompt 2, and the radius uses an in-house ln (range reduction plus a fixed polynomial, explicit operations, no std or libm call) and IEEE sqrt, so the GPU port matches bit for bit. Document the error of ln.
- noise(seed, entity, stream, sample_index) -> (f32, f32) with documented variance.

Tests: known-answer vectors; moments plus an Anderson-Darling or chi-square test over 1e7 samples; cross-correlation between adjacent streams; golden hash of 1e6 samples. Follow the Definition of Done in CLAUDE.md.
```

**Prompt 4: Units, scenario schema skeleton and `sim lint`**

```text
Read CLAUDE.md, docs/progress.md, and docs/design.md "Interfaces" (Scenario files, CLI) and "World model". Scenario files are TOML.

1. iqsim-units: typed quantities parsed from strings with units and deserializable with serde:
   Frequency (Hz, kHz, MHz, GHz); Power (W, mW, dBm, dBW); Gain (dB, dBi); Distance (m, km, ft, nmi); Speed (m/s, km/h, kt); Altitude (m, ft, FL, with an explicit reference: ellipsoid, MSL or AGL); Angle (deg, rad, DMS); Duration (s, ms, us, ns); Position (latitude/longitude in decimal degrees, DMS or MGRS, plus altitude).
   Reject unitless numbers where a unit is required, with a helpful error. For MGRS, implement it or pick a vendorable crate, and test against published reference points.
2. iqsim-scenario: serde types for the top level of a scenario: metadata, epoch (UTC), seed, run mode, receivers (channels with center frequency, sample rate, bandwidth), platforms (position, with trajectory and mounts as placeholders), emitters (behavior and waveform references as placeholders), library search paths. Use deny_unknown_fields. Generate JSON Schema with schemars into schema/scenario.schema.json through an xtask command and check it in.
3. `sim lint <file>`: parse and validate; errors point at file:line:column using toml spans; non-zero exit code on error.
4. examples/minimal.toml (one static CW emitter, one receiver) lints clean. Add negative tests for common mistakes: missing unit, unknown field, bad MGRS.

Follow the Definition of Done in CLAUDE.md.
```

## Phase 2: Static render on CPU

**Prompt 5: Work list and CPU render backend**

```text
Read CLAUDE.md, docs/progress.md, and docs/design.md "Architecture", "Time, numerics and determinism", "Receiver model" (Band-limiting before sampling) and "Emitters, behaviors and waveforms" (Rendering paths).

1. iqsim-worklist: plain-data #[repr(C)] types a GPU can consume unchanged, with no pointers (indices into flat arrays only):
   PathSegment (receiver channel, block, delay BlockPoly, phase BlockPoly, f32 amplitude polynomial, emission index); PulseTrain (pulse start offsets in Q32.32 samples relative to emission time, widths, edge-table references, per-pulse amplitude); CwEmission. A BlockWorkList holds everything one block needs.
2. Choose the fixed-point accumulator scale (fractional guard bits below the ADC LSB), document it in iqsim-core and in CLAUDE.md, and justify it against the dynamic range in docs/design.md.
3. Channel response: add the channel's analog front-end and anti-alias response H(f) to the receiver channel schema as a documented parameterized model (passband, transition band, stopband attenuation, group delay).
4. Band-limiting before sampling, applied per path:
   * Narrowband emitters (CW now, modulated later): fold H at the emitter's offset frequency, gain and phase, into the path's amplitude and phase. Cull paths whose in-band contribution is below a documented threshold.
   * Pulses: complex edge tables (the edge shape convolved with H shifted to the emitter's offset), built deterministically at scenario start (libm, fixed order) and shared by every path with the same edge shape, channel and offset. A pulse renders as the difference of two shifted edge responses. Choose the table resolution so interpolation error stays under a documented bound.
5. iqsim-render: a Backend trait, render(&BlockWorkList, accumulators per channel as &mut [i64]). CPU backend: each worker takes a chunk of output samples, jumps every overlapping path's polynomials to the chunk start, then steps sample by sample: band-limited envelope at the exact delayed time, carrier rotation from the phase polynomial and sincos, amplitude scaling to accumulator units, round to i64 and add. Integer accumulation makes the result independent of chunking and path order.
6. A temporary static-geometry builder: emitter and receiver at fixed ECEF positions, free-space loss, isotropic antennas, constant delay and phase per block.

Tests: a single pulse lands at the expected fractional sample (estimated time of arrival within 0.001 samples after removing the filter's group delay); a CW tone lands at the expected frequency and at the level H predicts; an emitter 60 MHz outside the channel leaves no aliased image above the noise floor; output hashes are identical across 1, 2 and 8 threads. Follow the Definition of Done in CLAUDE.md.
```

**Prompt 6: Receiver chain, SigMF output and first validation checks**

```text
Read CLAUDE.md, docs/progress.md, and docs/design.md "Receiver model", "Interfaces" (IQ outputs) and "Validation strategy".

1. iqsim-receiver: the chain from the design, in order, each stage switchable: front-end gain and noise figure; thermal noise (kTBF) from the Philox noise generator, shaped by a digital equivalent of the channel response H(f) from prompt 5; LO frequency offset (phase noise as a documented stub); IQ imbalance and DC offset; optional digital decimation filter; AGC off by default; ADC quantization and clipping to sc16. The analog anti-alias response is not applied here, because prompt 5 applies it to each path before sampling. Also expose the antenna-port tap, taken before any receiver effects.
2. iqsim-io: a SigMF writer (ci16_le and cf32_le) with core:sample_rate, core:frequency and core:datetime (truth time for now; the clock model arrives in prompt 15). Truth annotations go to a separate sidecar file, never into the consumer's .sigmf-meta.
3. iqsim-runtime and iqsim-cli: `sim run --offline <scenario> --out <dir>`, wiring scenario, static geometry, CPU backend, receiver chain and SigMF writer.
4. Validation tests from the design's table: delay; power within 0.1 dB of a closed-form link budget; noise floor within 0.2 dB of kTBF plus noise figure over the channel's noise bandwidth (Welch PSD); out-of-band rejection (a 1030 MHz emitter into a 1090 MHz channel leaves no aliased image, and in-band levels match H within 0.1 dB). Golden hash for examples/minimal.toml.

Follow the Definition of Done in CLAUDE.md.
```

## Phase 3: Motion and channel

**Prompt 7: Frames, geodesy and altitude**

```text
Read CLAUDE.md, docs/progress.md, and docs/design.md "World model" (Frames, Mode C altitude). In iqsim-geo, all in f64 using libm:

- WGS84 geodetic <-> ECEF (document the method and convergence), ECEF <-> ECI as a fixed-rate rotation about z from the scenario epoch, local NED at a reference point, body frame from attitude quaternions, antenna frame from mount orientation. One quaternion type with documented conventions (Hamilton; active vs passive), pinned by tests.
- Geoid: load an EGM96 or EGM2008 grid from a file path in config (no download at build time), bilinear interpolation. If no grid is configured, geoid height is 0 and a warning goes into run metadata.
- Pressure altitude: ellipsoid height -> MSL -> ISA pressure altitude, with an optional local pressure offset. Gillham (Mode C) encoding in 100 ft steps, tested against published Gillham tables.

Tests: round trips within 1e-9 m; published geodetic/ECEF reference points; NED and body rotations against hand-computed cases. Follow the Definition of Done in CLAUDE.md.
```

**Prompt 8: Trajectories and attitude**

```text
Read CLAUDE.md, docs/progress.md, and docs/design.md "World model" (Motion sources, Runtime maneuvers). In iqsim-geo:

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
Read CLAUDE.md, docs/progress.md, and docs/design.md "Time, numerics and determinism" (Error budget) and "Antennas and propagation" (Light time).

1. For each path, solve |r_rx(t) - r_tx(t - tau)| = c*tau in ECI by fixed-point iteration with a documented stopping rule. Use platform origins for now; prompt 10 switches to antenna phase centers.
2. Fit per-block delay and phase polynomials (order <= 3, documented nodes such as Chebyshev) to exact solutions, then convert to fixed-point BlockPoly. Phase = -f_c * tau in cycles, plus any carrier offset.
3. Check residuals at interior points against the exact solution. If delay or phase error exceeds the configured bounds (defaults 0.001 samples and 0.01 rad), split the block and refit. Log every split with block, path and residual.
4. Replace the static-geometry builder from prompt 5 with this geometry engine.

Tests: a moving CW emitter's measured frequency and phase trajectory match the analytic solution within the bounds; a LEO-like fast mover forces splits and still meets the bounds; hashes are stable across thread counts. Follow the Definition of Done in CLAUDE.md.
```

**Prompt 10: Antennas, mounts and airframe blockage**

```text
Read CLAUDE.md, docs/progress.md, and docs/design.md "Antennas and propagation" (Pattern representation through When gains are evaluated) and "World model" (Platforms and mounts). In iqsim-antenna:

- PatternGrid: complex gain on a theta/phi grid with polarization as a Jones vector, bilinear interpolation, documented frame conventions.
- Parametric library: isotropic; dipole; monopole and blade over a finite ground plane (document the approximation); SSR large-vertical-aperture antenna with sum, difference and control patterns (parameters: azimuth beamwidth, elevation shaping, control-pattern fill); parabolic dish from diameter and efficiency; circularly polarized patch; arrays as element pattern x array factor.
- Mounts: lever arm and orientation in the body frame. The light-time solver now uses antenna phase centers.
- Polarization mismatch from TX and RX Jones vectors rotated into a common frame.
- TX gain evaluated at emit time, RX gain at receive time; antenna phase goes into the phase polynomial.
- Airframe blockage: a fuselage cylinder plus wing and tail plates per platform. A path whose ray from the antenna hits a primitive gets shadow loss (document the loss model). Gains, losses and mismatch fold into the amplitude polynomial.
- Scenario TOML for mounts, patterns and airframes.

Tests: each parametric pattern integrates to its expected directivity within tolerance; linear-to-circular mismatch is 3 dB; a belly antenna loses sight of a ground site on the outside of a banked turn; a lever arm shifts delay by the expected amount. Follow the Definition of Done in CLAUDE.md.
```

## Phase 4: Mode A/C end to end, offline

**Prompt 11: Behavior engine and truth log**

```text
Read CLAUDE.md, docs/progress.md, and docs/design.md "Emitters, behaviors and waveforms" (Behavior engine) and "Interfaces" (Exports). In iqsim-behavior:

- A discrete-event engine: events ordered by (TruthTime, entity id, sequence number). Handlers are pure functions of state and event; they may schedule new events and produce Emission records. No wall-clock access.
- A link query API that behaviors use to ask what a given antenna on another platform would see from an emission: arrival time, received power and per-component gains, using the geometry and antenna code. Link-budget level, not IQ.
- Lookahead: the engine runs ahead of the render horizon by the configured margin and hands each block every emission that can reach any receiver within it, accounting for maximum path delay and emission duration.
- Truth log records: emissions, causal links (the event that caused each emission), and entity states at a configurable rate.
- Wire emissions into the work-list builder.

Tests: ties in event ordering resolve deterministically; lookahead never misses an emission (adversarial test with long delays); truth logs are identical across runs. Follow the Definition of Done in CLAUDE.md.
```

**Prompt 12: SSR interrogator and Mode A/C transponder**

```text
Read CLAUDE.md, docs/progress.md, and docs/design.md "Emitters, behaviors and waveforms" (SSR interrogator, Mode A/C transponder, Rendering paths). Implement in iqsim-behavior and iqsim-waveform, and cite the source of every timing constant in a comment:

- SSR interrogator: rotating antenna using the sum and control patterns from prompt 10; configurable scan period and start azimuth; PRF with stagger or jitter patterns drawn from Philox; Mode A/C interlace pattern; P1/P2/P3 pulses 0.8 us wide, P1-P3 spacing 8 us (Mode A) or 21 us (Mode C), P2 2 us after P1 and radiated on the control pattern.
- Mode A/C transponder: reception decided at link-budget level through the link query (P1/P3 power vs sensitivity; P2 vs P1 for sidelobe suppression). Reply 3.0 +/- 0.5 us after P3 with jitter; reply pulses 0.45 us wide on a 1.45 us grid; framing pulses 20.3 us apart; the 12 code pulses in standard positions; SPI. Per-unit carrier offset within 1090 +/- 3 MHz; dead time and suppression windows. Mode A code from the scenario (settable at runtime later); Mode C from the Gillham encoder in prompt 7.
- Pulse shapes with parameterized rise and fall times.
- Scenario TOML for both, with defaults matching typical equipment.

Tests: pulse positions and widths match the standard timing; a transponder in the interrogator's sidelobes does not reply; replies carry the right code bits; dead time is enforced; all randomness is reproducible. Follow the Definition of Done in CLAUDE.md.
```

**Prompt 13: Reference scenario end to end, offline**

```text
Read CLAUDE.md, docs/progress.md, and docs/design.md "Purpose and scope" (Reference scenario), "Interfaces" (Exports) and "Validation strategy".

1. examples/reference_modeac.toml: 3 aircraft with Mode A/C transponders on maneuvering routes; 4 fixed SSR interrogators with different scan periods and PRFs; 1 fixed receiver with two channels at 1030 and 1090 MHz.
2. Exports: the truth log as Parquet (entity states, emissions, causal links, receiver sample indices) and the consumer-knowledge JSON (interrogator positions with optional survey error, nominal PRI and scan schedules, receiver site). Truth never appears in SigMF metadata.
3. An independent Mode A/C decoder under tests/tools that does not import iqsim-waveform or iqsim-behavior: pulse detection on IQ, framing-pair search, code extraction. Decode the reference run with it.
4. Tests: decoded squawks and altitudes match truth; reply-delay statistics match the model; two runs give identical output and truth hashes; `sim run --offline` completes the reference scenario.
5. Record render time against simulated time in docs/progress.md.

Follow the Definition of Done in CLAUDE.md.
```

## Phase 5: Real time and control

**Prompt 14: Real-time pipeline and pacer**

```text
Read CLAUDE.md, docs/progress.md, and docs/design.md "Time, numerics and determinism" (Blocks and latency) and "Architecture".

In iqsim-runtime, build the block pipeline: behavior and geometry stages run ahead, then render, receiver chain and sinks. Stages connect through bounded queues sized from the latency budget. Offline mode runs the same pipeline, unpaced.

- Real-time mode: a wall-clock pacer releases blocks to sinks at the sample rate, and each block must be ready before its release time. A late block is a hard error with a diagnostic (which stage, by how much), never a silent drop or zero-fill.
- Telemetry counters: render margin per block, queue depths, stage timings.
- A pacing abstraction, so a hardware sink (prompt 24) can replace the wall clock as the pacer.
- Config: block length, latency budget, thread counts.

Tests: offline and real-time modes give identical sample hashes for the same scenario; an artificially slowed stage triggers the underrun error; the reference scenario runs in real time on this machine, with its margin recorded in docs/progress.md. Follow the Definition of Done in CLAUDE.md.
```

**Prompt 15: Receiver clock and streaming outputs**

```text
Read CLAUDE.md, docs/progress.md, and docs/design.md "Receiver model" (Clock model) and "Interfaces" (IQ outputs).

1. Receiver clock model: offset, drift and aging, or a GPSDO-disciplined model with a documented, Philox-driven error process. Output timestamps and SigMF datetimes now come from this clock, never from truth time. The truth log records both.
2. ZMQ PUB IQ sink: raw sc16 or fc32 items readable by a GNU Radio ZMQ SUB Source (no tags), one socket per channel, configurable high-water mark. Prefer the zmq crate with libzmq built from source (vendored) so no system package is needed; confirm it builds offline.
3. VITA-49 sink over UDP, implemented in-house: signal data packets with integer and fractional timestamps from the receiver clock, and context packets (RF reference frequency, sample rate, bandwidth, gain, timestamp) on change and at a configurable interval. Document which subset of VITA 49.2 is implemented.
4. A small test consumer for each sink.

Tests: VITA-49 packets round-trip through an in-house parser; timestamps follow the clock model, not truth; the ZMQ consumer receives the exact sample stream (hash). Follow the Definition of Done in CLAUDE.md.
```

**Prompt 16: Control plane, event log and replay**

```text
Read CLAUDE.md, docs/progress.md, and docs/design.md "Interfaces" (Control plane) and "World model" (Runtime maneuvers). In iqsim-control:

- JSON messages in a versioned envelope {schema_version, msg_type, id, payload}, defined as Rust types. Generate docs/control.schema.json with schemars through xtask and check it in; it is the contract the existing range tool will adapt to.
- A ROUTER socket for commands, with acks carrying the effective time; a PUB socket for telemetry (entity states, link status, render margin, underruns).
- Commands: maneuver (through the trajectory append API), add and remove entity, enable and disable emitter, set parameters (squawk, power, frequency), noise or jamming injection, start, stop.
- Effective time = max(requested time, render horizon). An append-only event log with a deterministic encoding records each request, its effective time and the command.
- `sim replay <scenario> <event-log>` re-renders a run.

Tests: a scripted client sends commands during a real-time run, and the replay is bit-identical (sample and truth hashes); malformed and out-of-order commands get clear errors. Follow the Definition of Done in CLAUDE.md.
```

## Phase 6: Authoring and packaging

**Prompt 17: Libraries and templates**

```text
Read CLAUDE.md, docs/progress.md, and docs/design.md "Interfaces" (Scenario files) and "World model" (Platforms and mounts).

- Library system: platform templates (airframe primitives, mounts, default emitters), antennas, waveform presets, receiver profiles and interrogator presets. Each is a TOML file found on the library search paths, referenced by name, with parameter overrides in the scenario. Document and test resolution order and conflict rules.
- Ship a starter library: generic transport and business-jet airframes with top and bottom L-band blades; generic terminal and en-route SSR presets; a two-channel 1030/1090 MHz receiver profile; a wideband receiver profile.
- `sim new --template <name> <dir>` creates a scenario from a template.
- Convert examples/reference_modeac.toml to use the library. Its output hash must not change.

Follow the Definition of Done in CLAUDE.md.
```

**Prompt 18: Full lint and dry-run**

```text
Read CLAUDE.md, docs/progress.md, and docs/design.md "Interfaces" (CLI, IQ outputs).

- `sim lint`: everything from prompt 4, plus library resolution, unattached antennas, emitters on platforms with no mount, overlapping receiver channels, and plausibility warnings (for example a 1090 MHz transponder at 10 W).
- `sim dry-run`: no IQ rendering. Run behavior and geometry over the scenario and report per-link visibility over time, predicted SNR per link, emitters outside every receiver's band, terrain coverage (stub until prompt 22), polynomial split counts, and estimated render cost per second of scenario. Terminal summary plus --json.
- Output budget: for every sink, compute the required rate (sample rate x bytes per sample x channels, plus packet overhead for VITA-49 and ZMQ) and compare it with the link capacity configured for that sink (for example 1, 10, 25 or 100 GbE, or a disk write rate). Real-time runs that exceed it fail lint; offline runs only warn.
- Every lint and dry-run message is actionable: what is wrong, where (file:line), and how to fix it.

Tests: golden text output for the reference scenario and for a deliberately broken one, including an over-budget sink. Follow the Definition of Done in CLAUDE.md.
```

**Prompt 19: Packaging for air-gapped RHEL 9**

```text
Read CLAUDE.md, docs/progress.md, and docs/design.md "Requirements" (Easy installation).

- `cargo xtask package`: a release build, then an RPM (spec file or cargo-generate-rpm) and a tarball containing the binary, starter library, schemas, examples and docs.
- Vendored build: a `cargo vendor` flow, a .cargo/config.toml template that points at either vendored sources or an internal registry, and a script that verifies the build with networking disabled.
- The default build needs no CUDA. The CUDA backend sits behind a feature and loads the driver at runtime.
- docs/install.md and docs/quickstart.md: install, first lint, first dry-run, first offline run, first real-time stream into a GNU Radio flowgraph.
- If podman is available, add an xtask that installs the RPM in a clean RHEL 9-compatible container and runs the quickstart.

Follow the Definition of Done in CLAUDE.md.
```

## Phase 7: CUDA and multi-GPU

**Prompt 20: CUDA backend**

```text
Read CLAUDE.md, docs/progress.md, and docs/design.md "Compute backends and scaling", "Time, numerics and determinism" and "Receiver model" (Band-limiting before sampling). Implement iqsim-render-cuda behind a cargo feature.

- build.rs: when the feature is on, compile the kernels with nvcc into a fatbin for sm_75 and sm_86 with -fmad=false and no fast-math flags, and embed it with include_bytes!. If nvcc is missing, fail with a clear message.
- Host side on cudarc with dynamic loading. At startup, if no driver or device is found, fall back to the CPU backend and log it.
- Kernels port the iqsim-core arithmetic line by line: Newton-form polynomials with the same explicit 64-bit primitives (each thread jumps once to its chunk start, then steps with 128-bit adds); the sincos table in constant memory; band-limited edge tables in global memory; Philox, ln and Box-Muller; noise shaping. Each thread owns a chunk of output samples and accumulates every overlapping path into i64; no atomics. Elementwise receiver-chain stages run on the GPU, and output is converted to sc16 on the device.
- `cargo xtask gpu-test`: thousands of seeded random work lists rendered on CPU and GPU must hash identically, including edge cases at block boundaries, phase wrap, and GPU chunk sizes that differ from the CPU's.

GPU tests skip with a clear message when no GPU is present. Record CPU vs GPU throughput on this machine in docs/progress.md. Follow the Definition of Done in CLAUDE.md.
```

**Prompt 21: Multi-GPU partitioning and benchmark**

```text
Read CLAUDE.md, docs/progress.md, and docs/design.md "Compute backends and scaling".

- Partitioner: assigns receiver channels to devices first. It splits one channel's emitters across devices only when a single device cannot keep up; partial sums come back as i64 and are reduced on the host. A device can be any GPU or the CPU backend.
- Heterogeneous test: the same run split across GPU plus CPU, and across two CUDA streams on one GPU, hashes identically to an unsplit run.
- Memory: size buffers from each device's reported memory, with 4 GB T1000s as the minimum target.
- `sim bench`: synthetic dense and sparse workloads; reports path-samples per second per device, PCIe transfer rates and real-time margin; writes a JSON report.
- docs/performance.md: how to run the benchmark on the 4x T1000 host and what to record.

Follow the Definition of Done in CLAUDE.md.
```

## Phase 8: Terrain

**Prompt 22: Terrain data**

```text
Read CLAUDE.md, docs/progress.md, and docs/design.md "Antennas and propagation" (Terrain). In iqsim-terrain:

- An in-house DTED reader for Levels 0, 1 and 2: UHL, DSI and ACC headers, data records, checksums, latitude-dependent longitude spacing. Validate against the DTED specification.
- A GeoTIFF elevation reader (a pure-Rust TIFF crate is fine; parse the GeoKeys you need in-house).
- A tile index over a directory, a memory-mapped tile cache with a size limit, bilinear elevation queries, and terrain-profile extraction between two points with documented sampling.
- Run metadata records the terrain source, level and resolution. With no terrain, use a smooth Earth and warn.
- Ground-clamped trajectories from prompt 8 now use terrain, and `sim dry-run` reports terrain coverage.

Tests use small synthetic DTED and GeoTIFF files that the test suite generates. Follow the Definition of Done in CLAUDE.md.
```

**Prompt 23: Propagation paths**

```text
Read CLAUDE.md, docs/progress.md, and docs/design.md "Antennas and propagation" (Path set per link).

- Write the reference first: a Python script under tests/reference, written from the ITU-R P.526 text before any Rust code, that generates expected knife-edge diffraction losses for a set of test profiles.
- Line of sight over the terrain profile with Earth curvature and a refraction k-factor (default 4/3).
- Diffraction: P.526 multiple knife-edge (state which method you implement) as an amplitude loss on the direct path.
- Ground reflection: find the specular point on the terrain; compute its own delay and phase, Fresnel reflection coefficients by polarization with ground-constant presets from the ITU-R P.527 categories, and a roughness reduction factor. It becomes its own path in the work list.
- Optional Longley-Rice (ITM) excess loss for beyond-line-of-sight links, using the NTIA C++ implementation behind a feature, or a documented hook if you defer it.
- Path updates run at the geometry rate and flow through the block polynomials.

Tests: flat-earth two-ray interference against the closed form; diffraction within 0.1 dB of the Python reference; line-of-sight masking behind a ridge. Follow the Definition of Done in CLAUDE.md.
```

## Phase 9: Hardware in the loop

**Prompt 24: SDR TX through uhd-rs**

```text
Read CLAUDE.md, docs/progress.md, and docs/design.md "Hardware in the loop" (including Device limits) and "Receiver model". In iqsim-hil:

- An SDR TX sink on uhd-rs (libuhd FFI) for the X410 and B205mini, one or more channels: continuous streaming of the antenna-port tap (no simulated receiver noise or ADC effects); timed start with UHD time specs; shared 10 MHz and PPS references when configured.
- The device becomes the pacer through UHD backpressure, using the pacing abstraction from prompt 14. Underflows reported by UHD are hard errors.
- Device capabilities: a built-in table per supported radio (channels, tuning range, bandwidth and sample-rate limits, host link), with values cited from the vendor specifications listed in docs/design.md. `sim lint` rejects HIL mappings whose bandwidth, sample rate, center frequency or channel count exceed the mapped radio. It warns when the host link cannot carry the total rate (for example four X410 channels at 400 MHz without dual 100 GbE) and when several radios share no 10 MHz and PPS reference.
- Calibration: a table per device serial, channel and frequency that maps digital full scale to dBm at the port, plus per-path losses (matrix ports, cables, external converters). `sim calibrate` steps frequency and level and records values entered from an analyzer or power meter.
- Scenario TOML: a hil.tx section mapping simulated receivers or channels to devices and ports.

Write a manual test plan in docs/hil.md: a B205mini transmitting into a second SDR through attenuators; verify power with an analyzer; decode the reference scenario's replies from the captured RF with the decoder from prompt 13. Tests: lint accepts valid mappings and rejects each kind of over-limit mapping. Follow the Definition of Done in CLAUDE.md.
```

**Prompt 25: SDR RX injection**

```text
Read CLAUDE.md, docs/progress.md, and docs/design.md "Hardware in the loop" (In: SDR RX -> simulation).

- An SDR RX source on uhd-rs that captures IQ from a device and exposes it as a scenario emitter placed on a platform. It passes through the simulated channel like any other emitter: band-limited to each channel's H(f), then the Farrow fractional delay (implement a minimal CPU version here if prompt 26 has not landed).
- All live input is recorded as SigMF with device timestamps, so `sim replay` reproduces the run bit for bit.
- The output runs a fixed latency behind real time; capture jitter beyond the latency budget is a hard error.

Tests: replaying a recorded live-input run matches the original output hash; with no hardware, a file-backed fake device exercises the same path. Follow the Definition of Done in CLAUDE.md.
```

## Phases 10–12: Waveforms, satellite and editor

These prompts open with a planning step and span several sessions; refine them once phases 1–9 have landed.

**Prompt 26: Waveform expansion**

```text
Read CLAUDE.md, docs/progress.md, and docs/design.md "Emitters, behaviors and waveforms" (Library roadmap). Work one family per session, in this order, each with an independent decoder or demodulator test. Before coding, propose a sub-plan for this session's family and wait for approval.

1. Continuous-waveform path, per path: shift the emitter's baseband to the channel frame, filter it to the channel band with H(f) (polyphase), then apply the Farrow fractional delay. Generate each emitter once at its native rate and share it across paths. Bit-exact between CPU and CUDA (extend the minimal version from prompt 25 if it exists).
2. Mode S and ADS-B: 1090ES extended squitters with CRC, and Mode S replies; DF17 position and velocity encoded from truth.
3. DME interrogator and transponder pulse pairs.
4. Generic pulsed radar: PRI patterns, pulse-compression codes, scan.
5. AM and FM; generic PSK and QAM with root-raised-cosine shaping.
```

**Prompt 27: Satellite**

```text
Read CLAUDE.md, docs/progress.md, and docs/design.md "Antennas and propagation" (Light time, Relay nodes) and "Emitters, behaviors and waveforms" (Library roadmap). Plan first and wait for approval; split into sessions as needed.

1. SGP4 trajectories from TLE and OMM (in-house or a vendorable crate), converted from TEME into the simulator's frames.
2. Relay nodes: the link graph becomes a per-block DAG. Bent-pipe transponder with receive filter, frequency translation, amplifier nonlinearity (Saleh or similar), and uplink noise and interference carried through.
3. DVB-S2-structured carriers: PL framing, pilots, symbol rate, rolloff, random payload. Tests detect the PL headers with an independent detector.
4. Dense-rendering performance on the T1000 host, compared with the capacity estimate in docs/design.md.
```

**Prompt 28: Map-based scenario editor**

```text
Read CLAUDE.md, docs/progress.md, docs/design.md and docs/control.schema.json. Plan first and wait for approval.

Build a web scenario editor served by `sim serve`. It works fully offline: a DTED shaded-relief basemap (no tile server); place and edit platforms; draw routes; attach library templates; save TOML that passes `sim lint`; run dry-run and show its report; start a run and show live telemetry over the control plane.

Acceptance: the reference scenario can be built and run without editing files by hand.
```
