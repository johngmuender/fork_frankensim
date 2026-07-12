# GUM Physics Core — Gap Analysis v4: the Scale Frontier

Four-surveyor measured survey (full structured output:
`scale_survey_v4.json`); successor to GAP_ANALYSIS_v3 after the ROADMAP_v3
build-out. Everything below is measured on this box (4-core Cascade
Lake-class Xeon, AVX-512 present, 33 MB L3, 15 GB RAM) or verified by
running code — not guessed.

## The measured facts that set the design

1. **Everything is single-threaded and the descent loop is 99% of all
   wall-clock.** All gate binaries show user ≈ real; the three ANF
   descents are 147 of gum_statics_gates' 149 s run-1. The only thing
   worth parallelizing is the sweep inside `engine.rs::eval()`.
2. **The kernels are compute-bound by ~20×, and scalar-issue-bound at
   that**: ~700 FLOP/point against 0.3–0.5 GB/s DRAM demand (11.7 GB/s
   available single-thread). Achieved 1.9–3.1 GFLOP/s vs 11.3 GFLOP/s
   measured SSE2 peak. Cache behavior benign to N=192 (5-plane stencil
   window ≪ L3).
3. **Compiler flags are a dead end — measured**: `-C target-cpu=native`
   makes the real kernels 10–25% *slower* (the per-point `[f64;4]`
   local-array style defeats the autovectorizer) while speeding up
   synthetic FMA loops 1.7×. The ~15× single-box headroom (≈4× SIMD ×
   3.96× threads) lives in *restructured k-contiguous sweeps*, not flags.
4. **Multiprocess scaling is already perfect** (4 copies: 3.96×
   throughput) — task-farm workloads need no new architecture.
5. **Two determinism hazards block naive threading**, both fixable by
   construction: (a) fixed-order reduction sums → per-tile partials in a
   thread-count-independent tile decomposition combined by a fixed tree
   (bit-identical serial↔threaded BY CONSTRUCTION; costs one deliberate,
   documented golden bump since the tiled order ≠ the current flat
   order); (b) the 13-point adjoint SCATTER is a data race → invert to
   GATHER, which simultaneously enables SIMD. (c) ANF's arrest cascade
   amplifies 1 ulp into divergent trajectories → cross-BACKEND goldens
   must be endpoint tolerance-band/metric goldens (the v2 assessment's
   sanctioned fallback); bit-identity is reserved for same-backend and
   CPU-serial↔CPU-threaded classes.
6. **The GPU rung's f64 question is answered concretely and TESTABLY
   HERE**: wgpu/naga support f64 WGSL via `SHADER_F64`/
   `Capabilities::FLOAT64`, and the surveyor verified it working
   headlessly on **lavapipe** (software Vulkan) in this container. At
   AI ≈ 5 FLOP/B the sweep kernel sits near the memory roofline of an
   A100-class part: ~3000× this box's serial throughput. Device-resident
   field, ~12 scalars read back per iteration, 4 per-component buffers
   to stay under binding limits to N ≥ 256.
7. **fs-exec remains unusable here** (asupersync sibling absent) — its
   `pairwise_fold`+`Compensated` (~200 LOC) can be vendored with the
   documented-swap-back precedent; std::thread::scope suffices (no new
   deps).

## What the scale buys (physics demand, ranked by value/FLOP)

- **G*-vs-𝔠₀ precision question** (2-D spectral, single-core, sub-hour):
  does the halo-saturated objective G* equal 𝔠₀ = 128√42/105π exactly?
  Currently agree to ~1.2×10⁻⁴ with G* still drifting down. Potentially
  a new exact identity for the F-R5 record. Zero architecture risk — do
  it first.
- **The I.2 two-knot bond-equation loop — the flagship** (never executed
  anywhere; grep-verified): two knots on 144×96×96 at the frozen
  h = 0.09375, 3 orientations × separations, gradient-flow relax, fit
  ℬ_eff, close the bond equation x₀ = 1.90 ± 0.05 from 𝔟 = 42 vs the
  corpus's measured 1.92 ± 0.08 (stretch: the ⟨r6⟩ locked bond
  2.42 ± 0.12). ~10¹⁴ FLOP ≈ 2.4 h as a 4-process task farm on the
  EXISTING serial engine — no new kernel needed; per-run bit-identity
  intact; the fitted observable is exactly the tolerance-band golden
  class the parallel rungs need anyway.
- **N=192 F-R5 halo endpoint** (8× cells vs the campaign's N=96): the
  threaded kernel's proof workload; hardens the DISCHARGE_PACKAGE.
- **Genuinely GPU-class** (unchanged): blue-fog condensate/disclination
  boxes, tower-depth studies at scale.

## Verdict

The scale gap is real but narrower and cheaper than the v2 assessment
priced it: one dep-free kernel crate (tiled deterministic sweeps,
scatter→gather, serial/threaded/wgpu backends), one deliberate golden
bump, tolerance-band goldens only at backend boundaries — and the
science payloads (two-knot, G*-identity, N=192 endpoint) are all
reachable on this box. ROADMAP_v4_SCALE.md executes it.
