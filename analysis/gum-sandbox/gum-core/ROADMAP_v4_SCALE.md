# GUM Physics Core — Roadmap v4: Scale (CPU → GPU)

Successor to ROADMAP_v3 (executed in full). Closes the scale gap of
GAP_ANALYSIS_v4_SCALE.md. Design rules carried forward: machine-checked
gates against frozen referees; deterministic-by-construction wherever
possible; deliberate, documented golden bumps when a canonical order
changes; tolerance-band/metric goldens only at backend boundaries (GPU,
and ANF endpoints across backends); standalone crates under gum-core/.

## Phase G1 — `fs-gum-kern`: the deterministic parallel kernel layer ▶ EXECUTING
The one parallel primitive: a domain-decomposed sweep (stencil-read →
per-point physics → per-point flux records) with (a) per-tile partial
reductions in a thread-count-INDEPENDENT tile decomposition (i-slabs,
TILE_I = 4, part of the bit contract) combined by a fixed pairwise tree,
and (b) the adjoint inverted from scatter to GATHER. Backends: serial
(the new canonical order — one deliberate golden bump, documented and
re-frozen) and threaded (std::thread::scope, tile dispenser; bit-identical
to serial BY CONSTRUCTION — gate asserts hash equality at 1/2/3/4
threads). Vendors fs-exec's pairwise_fold/Compensated (~200 LOC,
swap-back documented — the B3/F3 precedent). Gates: equivalence to the
old engine at machine-eps-class once, then re-frozen bit goldens;
measured speedup table (target ≥ 3× at N=96 on 4 threads); fs-gum-statics
descent gates re-run through the new layer.

## Phase G2 — wgpu f64 backend on lavapipe ▶ EXECUTING
The GPU rung, TESTED HERE via software Vulkan: WGSL f64 compute kernels
(SHADER_F64 / naga FLOAT64) for the forward sums + flux + gather sweep,
device-resident field, per-component storage buffers, fixed-workgroup
reduction with a deterministic combine on readback. Goldens: the
tolerance-band/metric class at kernel level (sector sums vs CPU within
stated bands) and at ANF endpoints (deg/E_static/κ/halo bands) — the
sanctioned fallback, stated per gate. Lavapipe throughput is irrelevant
(it is a correctness rig); the deliverable is a working, gated backend
that a real GPU inherits unchanged.

## Phase G3 — the G*-vs-𝔠₀ precision push ▶ EXECUTING
Single-core, 2-D spectral, sub-hour: drive the halo-saturated objective
G* = min(E_static − I/16π) to ≤ 1e-7-class precision (basis/grid
ladders, Richardson) and adjudicate whether G* = 𝔠₀ = 128√42/105π
exactly (currently 1.2×10⁻⁴ apart, drifting down). Either outcome is a
result: a new exact identity for the F-R5 record, or a pinned gap.

## Phase G4 — the I.2 two-knot bond-equation flagship ▶ EXECUTING
The corpus protocol never executed anywhere: two B=1 knots at
separations, product-ansatz seed, guarded relaxation on 144×96×96
(h = 0.09375 frozen), 3 relative orientations; extract the interaction
energy vs separation, fit the screened-dipole ℬ_eff form, close the
bond equation and compare x₀ = 1.90 ± 0.05 (𝔟 = 42) vs the corpus's
measured 1.92 ± 0.08. Runs as a 4-process task farm on the EXISTING
serial engine (measured 3.96× multiprocess scaling) — per-run
bit-identity intact; the fitted x₀ is a tolerance-band observable by
nature. Reduced protocol if budget forces it (fewer separations,
documented), full protocol as the follow-on.

## Phase G5 — the N=192 F-R5 endpoint (after G1)
The threaded kernel's proof workload: the halo-descent referee at 8×
the campaign's cells, refining the F-R5 endpoint and hardening the
DISCHARGE_PACKAGE. Scheduled once G1's speedup lands.

## Out of envelope (unchanged)
Blue-fog condensate/disclination boxes and tower-depth studies at scale
(real-GPU class — the G2 backend is their on-ramp); gravitation;
nucleation.

## Execution record
- Phase G-A (measured scaling survey, 4 surveyors): ✅ this document +
  GAP_ANALYSIS_v4_SCALE.md + scale_survey_v4.json.
- G1–G4: launched as parallel agents (G4 farms the existing serial
  engine and does not wait on G1).
- G5: queued behind G1.
