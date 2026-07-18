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

**4/4 gates PASS.**  The P2 stability margin survives at population scale:
all 10,000 screened configurations classify pre-crossed (9,991) or
box-exit (9, the P2 semantics' benign class), the population minimum
margin is 0.62156 at (n_y = 0.75, t0 = 6.404, z0 = 0.6929) — consistent
with P2's own reference minimum 0.62068 — and the 24 lowest-margin
configurations all re-certify **STABLE** through the unmodified 200-bit
ladder with 9–10 certified digits and zero UNSTABLE findings.

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

| gate | requirement | measured | verdict |
|---|---|---|---|
| S4-G1 | float64 screener reproduces P2 margins on 24 seeded (seed 20260718) class-(a) P2 configs to ≤ 1e-6 relative | max 2.28e-7, p95 1.53e-7, median 3.54e-9; all 24 classifications preserved; attempt-1 FAIL at 1.392e-6 archived, wall rtol tightened 1e-10 → 1e-11 BEFORE production, bar never moved | **PASS** |
| S4-G2 | population screen complete at final N with margin distribution + minimum located; nohup + chunk checkpoints | N = 10,000 (full target, no scale-down), 100 resume-safe chunks, ~96 min wall; classes 9,991 A_PRECROSSED / 9 BOX_EXIT; margin percentiles (0/5/50/95/100) = 0.622 / 1.324 / 4.850 / 17.52 / 28.82; minimum 0.621562 at (0.75, 6.404, 0.6929); clips never bind; min relative ρ 1.95e-18 above floor | **PASS** |
| S4-G3 | K = 24 lowest-margin configs through the unmodified P2 200-bit ladder; STABLE/UNSTABLE per config with certified digits; any UNSTABLE reported | 24/24 STABLE, 0 UNSTABLE; screener-vs-mp relative median 8.9e-10, max 4.3e-9; certified digits 9–10; step-halving shifts ≤ 3.5e-11 (3 lowest); DOP853 rtol-1e-12 cross-checks ≤ 3.8e-11 (3 lowest) | **PASS** |
| S4-G4 | honest scope statement: screened coverage ≠ formally-verified kernel | §5 verbatim; formal-kernel open REMAINS OPEN | **PASS** |

## 4. Key numbers

| quantity | value |
|---|---|
| population N (families × NT × NZ) | 10,000 = 2 × 100 × 50 (all t0 off the O1 grid) |
| class counts | 9,991 A_PRECROSSED / 9 BOX_EXIT / 0 other |
| margin percentiles 0/1/5/25/50/75/95/100 | 0.6216 / 0.8550 / 1.3242 / 3.1013 / 4.8496 / 8.6533 / 17.5214 / 28.8196 |
| population minimum (screener) | 0.621561690406 at (n_y = 0.75, t0 = 6.404, z0 = 0.69291) |
| population minimum (200-bit mp) | 0.621561690483489914… (certified digits 9) |
| P2 reference minimum (432 configs) | 0.6206790 — population min sits 0.14% above it |
| worst-K certification | 24/24 STABLE, 0 UNSTABLE |
| screener vs mp relative (worst-K) | median 8.9e-10, max 4.3e-9 |
| G1 validation (24 seeded P2 configs) | max 2.28e-7, median 3.54e-9; attempt-1 1.392e-6 FAIL archived |
| bench-frozen h policy | interior RK4 h = 1e-3 (≤ 1.85e-7); near-wall h = 5e-4 (≤ 2.4e-8); wall band DOP853 rtol 1e-11 (≤ 1.0e-7) |
| batched evaluator vs fields_exact | ~9e-15 |
| wall time | screen ~96 min (4 workers, 100 chunks) + ladder ~40 min (30 items) + gates/fig ~3 s |

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

1. **The screener is a hybrid by measured necessity** (printed spec
   amendment, coordinator-owned): a pure fixed-step float64 RK4 cannot
   hold 1e-6 on the extreme wall rows (cot(πz) drift |v_x| ~ 400 vs field
   wavelength ~0.06 degrades convergence to ~first order), so the wall
   band (dw < 0.02, ~8% of configs) runs the verbatim P2 DOP853 pipeline
   at rtol 1e-11 instead.  This is a method finding of the workstream,
   not a silent deviation.
2. **The 9 BOX_EXIT configs** inherit P2's classification semantics
   verbatim (trajectory leaves the integration box before the horizon —
   handled as its own class, not a stability failure); all 9 sit at the
   extreme wall rows.
3. The G1 attempt-1 calibration FAIL (1.392e-6 > 1e-6 at wall rtol 1e-10)
   is preserved in `_validate_attempt1.json`; the fix (rtol 1e-11) was
   made before the production screen started and the gate bar was never
   moved.
4. Certified digits for worst-K margins are limited by the
   screener-vs-mp agreement (9–10 digits), not by the 200-bit ladder
   itself; step-halving and DOP853 cross-checks on the 3 lowest configs
   agree to ≤ 3.8e-11.
5. Population coverage is still a *sample* of the continuous
   (t0, z0) rectangle — 23× the P2 grid, minimum located interior to the
   z-range and at the t0 lower edge (the kill-bin lower edge, where P2
   also found its minimum band).  Nothing here is a machine-checked
   inequality (§5).
