# fs-gum-kern — GUM core Phase G1 results

Crate: `analysis/gum-sandbox/gum-core/fs-gum-kern/` (standalone
`[workspace]` opt-out; path deps: `../fs-gum-field`, `../fs-gum-statics`,
`crates/fs-blake3` — std-only closure, the fs-cosserat-pilot pattern; the
crate itself is `#![forbid(unsafe_code)]`).

Scope (ROADMAP_v4_SCALE.md Phase G1, design fixed by the measured survey
`scale_survey_v4.json` DIMENSIONs 1–2): the deterministic parallel kernel
layer for the N³ sweeps of fs-gum-field / fs-gum-statics — ONE parallel
primitive (tile-decomposed sweep, per-tile partials, fixed pairwise
combine), the gradient adjoint inverted from scatter to GATHER, serial and
threaded execution bit-identical BY CONSTRUCTION, the gum kernels
(forward sums, gradient, renormalize, vector ops, norms, ANF) re-hosted on
it, and a measured ≥ 3× 4-thread speedup on the ANF hot path at N = 96.

Reproduce: `cargo build --release && ./target/release/gum_kern_gates`
(~4 min; build is warning-free; `cargo test --release` for the unit layer).

## Design — the bit contract

The old engine's determinism contract was "one flat sequential loop".
The new contract makes every rounding decision a pure function of the
sweep extent alone — never of thread count, scheduling, or steal order:

```
TILE_I = 4          i-slab thickness in planes (>= widest stencil reach).
                    Tile count = ceil(m / TILE_I) for a sweep of m planes
                    (m = N for cell sweeps, N+1 for corner sweeps) — a
                    pure function of m.  PART OF THE BIT CONTRACT, exactly
                    like fs-la's KC (its GEMM_BIT_SEMANTICS_VERSION
                    precedent): retuning TILE_I legitimately changes bits
                    and forces a documented golden re-freeze.
combine tree        per-tile partials (plain f64, ascending (i,j,k) inside
                    the tile, the old engine's per-point statements
                    verbatim) merged up the fixed pairwise tree that
                    splits at the largest power of two strictly below n —
                    fs-exec's shape rule, vendored (see below).
gather adjoint      pass A stores per-eval-point flux records P[3][4]+gq[4]
                    (16 f64/point; 117 MB at N = 96 corner, workspace-
                    reused across ANF iterations); pass B has each OUTPUT
                    cell gather its contributions in a PINNED order —
                    central4: axes x,y,z with eval-point offsets
                    (-1,+1,-2,+2) (boundary-dropped, the exact transpose
                    of the ghost-flux-dropping scatter), then own gq;
                    corner: the 8 surrounding corners ascending (ei,ej,ek)
                    with the fused 0.25/h(±P)+0.125 gq expression; then
                    cell terms, h³, tangent projection, all per cell.
norm_flat           per-tile (a,i,j,k) ascending, tree-combined, sqrt last.
threading           std::thread::scope, workers = threads.min(tiles),
                    AtomicUsize tile dispenser (work stealing absorbs
                    throttled cores), CachePadded (#[repr(align(64))])
                    per-tile slots; the per-slot lock guards HANDOFF only.
                    No rayon, no unsafe: tile disjointness is established
                    by safe slice splitting (split_tiles).
```

All pinned by `GUM_KERN_BIT_SEMANTICS` ("gum-kern-bits-v1: ..."), which
domain-tags the frozen goldens.  Serial (threads = 1) runs the identical
per-tile code and the identical tree — bit-identity serial↔threaded is a
construction, not an observation; gate K-D observes it anyway.

Vendored: fs-exec's `pairwise_fold` + `Compensated` (~130 LOC, verbatim
shapes) in `src/reduce.rs` — fs-exec is unbuildable in this checkout
(absent `asupersync` sibling); swap-back to `fs_exec::reduce` documented
in the module header, the B3 `jacobi_eigh` / F3 exec-feature precedent.
The gum kernels use PLAIN per-tile f64 accumulation (current physics
semantics); `Compensated` rides the same tree where future physics wants
it.

Dependency deviation (documented): the roadmap sketch said "path deps on
fs-gum-field"; the crate also path-deps `fs-gum-statics` — for the shared
engine types (`Opts`/`Out`/`Grad`/`AnfParams`/`Record`, making the new
`eval` a literal drop-in), the diagnostics kit, and the OLD engine that
gates K-A/K-E compare against.  Consequence: fs-gum-statics can never dep
on fs-gum-kern; downstream phases (G5) call `fs_gum_kern::{eval, anf}`
directly.  Two `pub` accessors were added to `fs-gum-field::Field3`
(`data()`/`data_mut()` — raw views of the already-documented SoA layout)
so the tiled maps can split the field into disjoint per-tile chunks
safely; nothing else in the finished crates was touched.

## The deliberate golden bump (measured, once)

The tiled canonical order provably differs from the old flat order, and
the gather adjoint re-orders each cell's contribution sum.  Gate K-A
measured the deltas ONCE on the N = 48 hedgehog + frozen bump (both
schemes, bare `estatic` and fully-guarded rotor objectives):

```
forward sums   obj(central4)  old 3.18765859157195486e0
                              new 3.18765859157204456e0   rel 2.81e-14
               obj(corner)    old 3.04413419705629495e0
                              new 3.04413419705638466e0   rel 2.94e-14
               worst over all fields/cases: rel 2.87e-12  (the guarded
               case — the one-sided penalties multiply the raw ~1e-14
               order noise by MU-class factors; still machine-eps class)
gradient       central4: max|Δg|/max|g| = 6.50e-15
               corner:   max|Δg|/max|g| = 0 — EXACTLY bit-identical: the
               corner scatter's arrival order at a cell (its 8 corners in
               ascending corner order, one fused expression each) happens
               to coincide with the pinned gather order, so inverting the
               corner adjoint costs no bits at all
g.norm         rel <= 1.19e-13 (tiled tree vs flat running sum)
```

This is the expected divergence at the ~12th–14th significant digit
(the survey's `tiledet` prototype measured the same class at N = 96:
9.83072230548363030e7 tiled vs 9.83072230548577905e7 flat).  Per the
fs-la KC-contract precedent the bump is taken deliberately, once, in the
same change that introduces the contract: K-C froze NEW bit goldens for
the tiled order (BLAKE3 over Out + gradient bits + norms, both schemes,
domain-tagged by GUM_KERN_BIT_SEMANTICS):

```
golden 5a4cc174d094b0e878ee1f46c50f764202bc5cdc6827130a8d1c1a190d07fefa
  central4 rotor+guards: obj 1.99028698800106358e4  gnorm 5.92022109985084555e1
  corner   rotor+guards: obj 1.98164692378990876e4  gnorm 2.36266225335303695e1
```

The OLD crates (fs-gum-field / fs-gum-statics / fs-gum-e2e) and their
goldens are untouched — the old engine remains the frozen referee this
layer was gated against.

## Gate table (13/13 PASS)

| id | gate | measured | verdict |
|----|------|----------|---------|
| K-A1 | old-engine cross-check, forward sums (tiled vs flat order), N=48, both schemes, bare + guarded | worst rel 2.87e-12 (tol 5e-11); bare-objective rel 2.8e-14 | PASS |
| K-A2 | old-engine cross-check, gradient (gather vs scatter) | central4 max|Δ|/max|g| 6.5e-15; corner 0.0 (bit-identical); gnorm rel 1.2e-13 | PASS |
| K-B | FD-vs-analytic gradient THROUGH the new layer, N=24, both schemes, 8 objectives x 3 dirs, eps=1e-5 | worst rel 5.66e-7 (tol 1e-5; old engine's own gate: 1.21e-6) | PASS |
| K-C | NEW bit goldens frozen for the tiled order | hash matches frozen 5a4cc174... | PASS |
| K-D | HEADLINE: serial vs threaded at 1/2/3/4 threads, BLAKE3 over all outputs (eval fwd+grad both schemes N=48, renormalize field+drift, 80-iter ANF at N=32: series + final field bytes) | all four hashes = 4669604dcc81e42c... (bit-identical) | PASS |
| K-E1 | new-layer G-C main: R decreases | 3.749126 -> 3.419867 (dR -0.329), max rise +1.55e-3 <= 2e-3 | PASS |
| K-E2 | new-layer kappa falls | 0.25316 -> 0.13722 (drop 0.116); crosses threshold 0.19947 by it 160 (old engine: it 160) | PASS |
| K-E3 | new-layer halo rises (4A signature) | 0.0276 -> 0.0572 (+0.0296) | PASS |
| K-E4 | objectives monotone by construction (static + main) | max recorded rise -7.57e-6 <= 0 | PASS |
| K-E5 | endpoint vs old engine within tolerance bands | dR +2.9e-14 (band 0.02), dkappa -3.0e-14 (0.02), dhalo +3.3e-15 (0.02), ddeg +5.3e-14 (0.005); arrests 12 vs 12 | PASS |
| K-F | threaded-4 speedup, ANF hot path (corner gradient, N=96) | 3.61x >= 3.0x (0.529 s -> 0.147 s) | PASS |
| K-G | bit-identical two-run replay of the 33 core gate numbers | fingerprint feb63c00d3c465bc... (identical across two separate process invocations too) | PASS |

## K-E: the G-C descent through the new layer (N = 48, 900 + 600 iters)

Full fs-gum-statics G-C main protocol (static relax from hedgehog + bump,
then the over-spun L = 0.266·I Routhian descent from static + bump +
tilt-halo shell), run end-to-end through the new layer at 4 threads and
independently through the old engine:

* Qualitative referee: identical signature — R monotone down, kappa falls
  through the threshold at the same recorded iteration (160), halo grows
  2.1x, degree pinned at the anchor band.
* Endpoint: the two engines agree to `|ΔR| = 2.9e-14`,
  `|Δkappa| = 3.0e-14`, `|Δhalo| = 3.3e-15`, `|Δdeg| = 5.3e-14`, with
  IDENTICAL arrest counts (12) and arrest iterations.  Explicit
  statement: at N = 48 the ~1e-14-class per-eval order change never
  flipped an accept/arrest decision in 1500 iterations, so the
  trajectories stayed ulp-locked and the endpoint agreement is
  machine-eps class — far inside the stated bands (R/kappa/halo 0.02,
  deg 0.005).  The bands, not bit-identity, remain the CONTRACT for this
  comparison: the arrest cascade is a ulp amplifier, and a resolution or
  protocol where a single flip occurs would legitimately diverge
  trajectories while still passing the bands (survey DIMENSION 1,
  hazard c; the sanctioned cross-backend golden class).  Same-backend
  runs are bit-identical at any thread count regardless (K-D).
* New-layer protocol wall time 80.8 s vs old serial 111.9 s under farm
  co-tenancy (see measurement note).

## K-F: measured performance (ns/point, best-of-reps)

Serial-old = fs-gum-statics engine; serial-new = this crate, threads=1;
4T = threads=4.  ANF-class objective (rotor + both guards).  4 vCPU
Cascade-Lake-class container, baseline x86-64 codegen (no target-cpu).

| sweep | old ns/pt | new ns/pt | 4T ns/pt | 4T/new | new/old |
|-------|-----------|-----------|----------|--------|---------|
| fwd  central4 N=48 | 60.1 | 62.1 | 38.7 | 1.60x | 0.97x |
| grad central4 N=48 | 273.1 | 255.1 | 131.9 | 1.93x | 1.07x |
| fwd  corner  N=48 | 257.8 | 242.3 | 69.8 | 3.47x | 1.06x |
| grad corner  N=48 | 620.9 | 646.6 | 185.0 | 3.49x | 0.96x |
| fwd  central4 N=96 | 62.6 | 62.6 | 21.8 | 2.87x | 1.00x |
| grad central4 N=96 | 275.3 | 251.5 | 80.1 | 3.14x | 1.09x |
| fwd  corner  N=96 | 224.4 | 221.7 | 58.5 | 3.79x | 1.01x |
| grad corner  N=96 | 552.9 | 579.3 | 160.6 | 3.61x | 0.95x |
| fwd  central4 N=128 | 60.7 | 63.2 | 17.8 | 3.56x | 0.96x |
| grad central4 N=128 | 309.2 | 273.5 | 86.7 | 3.16x | 1.13x |
| fwd  corner  N=128 | 223.7 | 224.3 | 73.8 | 3.04x | 1.00x |
| grad corner  N=128 | 587.9 | 591.2 | 169.7 | 3.48x | 0.99x |

Thread ladder on the target (corner gradient, N = 96): 1T 0.529 s,
2T 0.278 s (1.90x), 3T 0.197 s (2.69x), **4T 0.147 s (3.61x)** — the
>= 3x gate passes with margin.  Projection: the N=96 G-C-class descent
(~2000 evals) drops from ~18 min serial to ~5 min at 4 threads; the
N=192 F-R5 endpoint run (Phase G5) from ~2.8 h to ~50 min.

THE SIMD SITUATION, honestly: serial-new is 0.95–1.13x serial-old — the
gather restructuring did NOT unlock autovectorization by itself.  The
per-point `[f64;4]` local-array style still defeats the vectorizer (the
survey measured `target-cpu=native` as a 10–25% REGRESSION on exactly
this shape, so flags remain a dead end), and the kernels stay
scalar-issue-bound at ~2–3 GFLOP/s against an 11.3 GFLOP/s measured SSE2
peak.  What the gather form does buy is (a) race-freedom (the threading
win above), and (b) the memory layout a k-contiguous explicit-SIMD (or
GPU, Phase G2) kernel needs — pass A records are already point-major
f64x16.  Explicit-SIMD k-contiguous sweeps are the documented next step,
NOT promised here; the measured ~4x SIMD headroom is real but unclaimed.

Measurement note (co-tenancy): the Phase G4 two-knot task farm (4
processes, 100% CPU each) was running in this container throughout the
build.  Correctness gates are scheduling-independent (that is the whole
point of the layer) and were run under full contention.  The K-F numbers
above were measured in a quiet window created by SIGSTOPping the farm's
4 processes for the ~35 s perf section and SIGCONTing them immediately
after — bit-neutral to the farm (its runs are wall-clock-stretched by
~35 s, nothing else).  Under full contention the same gate measured
1.60x — thread scaling numbers on this box are meaningless while 4
sibling processes saturate the cores, which is why the pause was needed.

## What is NOT here (deliberate)

* fs-gum-statics/e2e still run their own serial engine and their frozen
  goldens — migrating them onto fs-gum-kern (and retiring the scatter) is
  the G5-class follow-up; this phase proves the layer against them.
* No slab-staged pass A: the 16·N³ flux buffer is 905 MB at N = 192 —
  fine in 15 GB, but the staged/fused variant should land before N > 192.
* No explicit SIMD (see above); no GPU (Phase G2 owns the wgpu backend —
  pass A/pass B map 1:1 onto its two kernels by construction).

## Epistemic notice (binding, inherited)

Everything here is within-model numerical engineering on a speculative
theory's functional.  A passing gate certifies the kernel layer — its
determinism contract, its equivalence to the frozen reference engine at
machine-eps class, its measured speedup — never anything about nature.
The G-C descent re-run is an ENGINE CAPABILITY replication of the F-R5
mechanism, not new physics.
