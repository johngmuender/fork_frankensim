# S4 — Population-scale certification screener (P2 extension) (F-T8-S4)

**Goal.** Extend the P2 exact-field certification (432 stratified
configurations) toward the engineering open "population-scale
certification": screen a ~23x-refined population drawn from the SAME
P2/O1 parametrization with a validated cheap float64 pipeline, locate the
margin distribution and its minimum, and re-certify the K = 24
lowest-margin configurations through the unmodified P2 200-bit ladder.

**Campaign:** Tier 8 Phase S (ROADMAP_v10_ALTERNATIVES.md, S4).
**Date:** 2026-07-18.  **Code:** `s4_screen.py` (phases bench / validate /
screen / ladder / gates / fig, all checkpointed; screen chunk files in
`s4_chunks/`, ladder per-config checkpoints in `s4_ladder/`).
**Raw output:** `s4_results.json` + intermediates (`_bench.json`,
`_validate.json`, `_validate_attempt1.json`, `_worstk.json`).
**Figure:** `s4_fig.png`.  **Log:** `s4_run.log`.
**Prior art (mandatory reads, honored):** tier7-program/P2/RESULTS.md +
p2_certify.py — imported VERBATIM (field spectrum, velocity law +
regularization, margin definition, classification semantics, DOP853
reference integrator, and the 200-bit gmpy2/MPFR ladder `mp_run_path`).

*Within-model; nothing here bears on nature.*

## 0. Verdict in one line

PLACEHOLDER

## 1. Population (parametrization reuse, pre-registered)

The P2/O1 parametrization per family n_y in {+1, +0.75} is t0 in the
kill bins [6.4, 7.2] (O1: dt = 0.01) x z0 in [0.0025, 0.9975] (O1:
dz = 0.005); the -n_y mirrors are numerically identical under z -> 1-z
(re-verified in P2 at 3.5e-8) and are not re-run, exactly as in P2.
S4 refines to N = 10,000 NEW configurations:

- t0: NT = 100 midpoints 6.4 + 0.8 (i + 1/2)/100 (dt = 0.008 — finer
  than O1's 0.01 and everywhere OFF the O1 t-grid);
- z0: NZ = 50 uniform 0.0025 + j (0.995/49) (dz ~ 0.0203; includes the
  extreme wall rows 0.0025 / 0.9975 exactly);
- families {+1, +0.75}: N = 2 x 100 x 50 = 10,000.

No scale-down was needed: the bench-projected cost fit the budget at
full N (the 10,000-target refinement ~23x the 432-path P2 subset).

## 2. The screener and its calibration (bench)

Float64 fixed-step classical RK4 on the exact 2048-mode sum — the SAME
integration scheme as the P2 mp ladder, in float64, batched (per-mode
time phases advanced by half-step recurrences; spatial factors by
cumulative-product power ladders + BLAS; evaluator agrees with
p2_certify.fields_exact at ~9e-15).  Margins refined by local h/20
re-integration on p2_certify.rhs_factory + 3-point parabola (the mp
ladder's own refinement scheme); s_pc by h/20 bracket + linear
interpolation; classification, monitors (relative-rho floor, clip
headroom), and box-exit semantics mirrored from p2_certify.

**Bench-frozen h policy** (calibrated against p2.integrate_path, DOP853
rtol 1e-12, 10 reference configs spanning the z rows at t0 = 6.404;
frozen in `_bench.json` BEFORE validate/screen ran):

| z-row band (wall distance dw) | engine | measured margin error |
|---|---|---|
| dw >= 0.05 (interior) | RK4, h = 1e-3 (R_INT = 2) | <= 1.85e-7 |
| 0.02 <= dw < 0.05 | RK4, h = 5e-4 (R_NEAR = 4) | <= 2.4e-8 |
| dw < 0.02 (wall band) | p2.integrate_path (DOP853) at rtol 1e-11 | <= 1.0e-7 |

**Calibration history (printed, coordinator-owned policy choices):**
(i) a pure fixed-step screener FAILS on the extreme wall rows — the
cot(pi z) drift gives |v_x| ~ 400 at z0 = 0.0025 against a field
wavelength ~0.06, and convergence degrades to ~first order (rel err
1.5e-5 at 10x steps, 8.7e-6 at 20x); the wall band therefore runs the
P2 DOP853 pipeline itself (hybrid screener).  (ii) The wall band's
rtol was first relaxed to 1e-10 for cost; the G1 attempt-1 run then
measured 1.392e-6 on one long wall path (ny = 0.75, t0 = 6.855,
z0 = 0.9975) — an honest first-attempt FAIL, archived in
`_validate_attempt1.json` — after which the wall rtol was tightened to
1e-11 (measured 1.0e-7 on that same config) and the seeded validation
re-run in full.  The production screen had NOT started at either
calibration step; the G1 gate bar (1e-6) was never moved.

## 3. Gates

PLACEHOLDER_GATES

## 4. Key numbers

PLACEHOLDER_KEY

## 5. Honest scope (S4-G4)

SCREENED POPULATION COVERAGE, NOT A FORMALLY-VERIFIED KERNEL.  The
N = 10,000 screened margins are float64 values validated to ~1e-7
against the P2 pipeline; the worst-K margins are precision-certified by
200-bit cross-validation; but NO margin here is a machine-checked
inequality.  Formal certification would still require validated
interval/Taylor-model enclosures of the backward flow with
outward-rounded mode sums over the full population — and S3's covering
pilot measured exactly that technology's present cost (wrapping-effect
dominated).  The formal-kernel open of the Phase-R closeout REMAINS
OPEN; the analytic T4-W5 proof remains the terminal open.  Within-model;
nothing here bears on nature.

## 6. Caveats

PLACEHOLDER_CAVEATS
