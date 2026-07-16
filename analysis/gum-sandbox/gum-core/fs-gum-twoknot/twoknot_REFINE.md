# fs-gum-twoknot — GUM core Phase G5b refine (longer-cap attractive channel)

Companion to `twoknot_RESULTS.md` (the G4 record, kept intact).  Scope:
tighten G4's two weakest numbers — the attractive-channel well depth
(2.1 sigma at cap 400) and the direct x0 band (±0.15) — by 3x-longer-cap
reruns of the ATTRACTIVE channel only in the well window, plus the far
anchor and the single-knot reference.  No physics code was changed; the
G4 binary was reused as-is (gates re-verified ALL PASS, 20/20, after the
box restart).

Reproduce: `sh farm_refine.sh runs1200 1200 && python3 analyze_refine.py`
(~43 min wall on 3-4 free cores; each run bit-deterministic for fixed args).

## Protocol as run

```
runs        attract m in {16, 18, 20, 23}  (x = d/R* in {1.6837, 1.8942,
            2.1046, 2.4203} — brackets x0 ~ 1.9-2.0) + far anchor m = 38
            (x = 3.9988) + single-knot reference; grid, seeding, guards,
            frozen dial ALL identical to G4 (twoknot_run.rs untouched).
cap         1200 UNIFORM (3x G4's 400).  The roadmap plan was 1600; it was
            cut to the sanctioned fallback after a container restart killed
            the first farm (~55%-done runs lost, no checkpointing — reruns
            are deterministic so nothing but wall time was lost) and G5a's
            4-thread n192_fr5 co-tenanted the 4-core box.
farm        farm_refine.sh: 2 batches of 3 processes (3-wide because of the
            G5a co-tenant); banking order — single + far anchor + the well
            point (m = 18) FIRST so the key observable survives another
            restart.  All 6 runs end at iter_cap (expected; the descent is
            guard-arrested, not converged — that is what the differential
            extrapolation below is for).
wall        42m34s farm (94m14s user); ~0.99 s/iter batch 1 (n192_fr5
            exited early, so nearly uncontended), ~1.13 s/iter batch 2.
            + ~40 min sunk in the restart-killed first attempt.
            Analysis ~5 s.  Total G5b wall ~1.5 h.
```

## THE CAP SYSTEMATIC (cap 400 -> 1200, far-anchored E_int, corner)

The headline measured fact: **G4's cap-400 E_int values in the bond zone
were dominated by under-convergence, not by the ±2.5e-4 differential
noise G4 quoted.**  The shift IS the cap systematic:

| m | x = d/R* | E_int (cap 400) | E_int (cap 1200) | shift 1200-400 |
|---|------|-------------|-------------|-----------|
| 16 | 1.6837 | +1.8422e-01 | +1.9626e-01 | +1.20e-02 (+6.5%) |
| 18 | 1.8942 | -5.2665e-04 | **-6.0011e-03** | **-5.47e-03** |
| 20 | 2.1046 | +2.2217e-04 | -3.8885e-03 | -4.11e-03 |
| 23 | 2.4203 | -2.2501e-04 | +4.8589e-04 | +7.11e-04 |

rms 6.9e-3, max |shift| 1.2e-2.  At the well point the shift is ~22x G4's
sigma_dc: the 400-cap runs had NOT yet expressed the pair's mutual
relaxation (the well deepens 11x when the descent is given 3x the
iterations; d_eff at m=18 pulls in further, 3.3552 -> 3.3419).  G4's
LOCATION of the well was right; its DEPTH was a large underestimate.
Guard penalties at the 1200 endpoints are uniform across separations
(epen 8.4-9.1e-4, efpen 1.0-1.3e-3, degrees 1.9664-1.9669 — the same
class as G4 and the certified single-knot runs), so the deepening is
descent, not guard drift.

Anchor cross-check E(far) - 2 E(single): +3.129e-2 at cap 1200 vs
+2.657e-2 at cap 400 — the B=2-vs-B=1 protocol offset GROWS with cap,
re-confirming G4's decision to define E_int against the far anchor
(which shares the pair protocol) rather than against 2 E_single.

## Convergence discipline (step-2 instrument)

Absolute E_static tails are still descending at cap 1200 (residual slope
~ -2e-5 to -3e-5 per iteration, every run iter_cap) — the absolute
energies remain protocol quantities.  The E_int instrument is the
DIFFERENTIAL series E_m(it) - E_far(it) on the common instrumentation
grid (every 10 iters): the shared descent cancels iteration-by-iteration,
leaving physics + arrest-desynchronization noise.  Cap extrapolation =
Aitken delta^2 on three block medians of the last-half tail (medians
because arrest-desync spikes are outliers — measured at cap 400, where
m=23's endpoint sat at -2.25e-4 vs a -6e-6 tail plateau):

| m | x | method | tail ratio r | E_int extrapolated | sigma |
|---|------|--------|------|-------------|--------|
| 16 | 1.6837 | median (no contraction, r=1.16) | — | +1.9330e-01 | 6.2e-03 |
| 18 | 1.8942 | **Aitken, r = 0.40** | 0.40 | **-6.5303e-03** | **1.30e-03** |
| 20 | 2.1046 | median (marginal, r=0.97) | — | -2.8130e-03 | 3.3e-03 |
| 23 | 2.4203 | Aitken, r = 0.03 | 0.03 | +2.1959e-04 | 1.8e-04 |

The m=23 extrapolation (+2.2e-4 where the seed estimator bounds physics
at <= 5e-6) directly measures the residual differential noise floor at
this cap: ~2-5e-4.  The m=18 differential tail is a clean geometric
contraction — the well depth extrapolates with a defensible error bar.

Refreshed noise model (single channel only at 1200, so G4's cross-channel
scatter is not re-measurable): sigma_ref = 2.30e-3 adopted from
max(|E_int(m23)|, min(median sigma_extrap, quiet-zone cap shift)) =
max(4.9e-4, min(2.3e-3, 4.1e-3)) — deliberately the CONSERVATIVE reading;
same exponential-decay form as G4, floored at 10% relative.

## The well (the 2.1-sigma number, revisited)

| estimator | E_int(x = 1.8942) | sigma | depth/sigma |
|---|-------------|--------|------|
| G4, cap 400 | -5.27e-04 | 2.5e-04 | 2.1 |
| cap 1200 endpoint, conservative noise model | -6.00e-03 | 2.30e-03 | **2.6** |
| cap-extrapolated (differential Aitken) | **-6.53e-03** | **1.30e-03** | **5.0** |

**The attractive well is real at 5.0 sigma on the cap-extrapolated
estimator (the instrument built for exactly this question); the raw
endpoint with the deliberately conservative single-channel noise model
gives 2.6 sigma.**  Depth ~1.0e-3 of the two-knot energy (G4: 1.7e-4).
The m=20 point also went negative (-3.9e-3): the well is now resolved as
a two-point feature, not a single-point excursion.

## x0 (bond location) — all estimators, corner

| estimator | x0 |
|---|---|
| direct parabolic minimum, cap 1200 (primary) | **1.997**, band [1.992, 2.002] |
| direct parabolic minimum, cap-extrapolated | 1.996, band [1.989, 2.002] |
| pair-law + core fit argmin, cap 1200 | 1.864 |
| **quoted (ensemble centre ± estimator spread)** | **x0 = 1.95 ± 0.07** (band [1.86, 2.00]) |
| G4 quoted | 2.0 ± 0.15 |
| corpus measured | 1.92 ± 0.08 |
| corpus bond equation at b = 42 | 1.90 ± 0.05 |

**x0 = 1.95 ± 0.07 vs corpus 1.92 ± 0.08: consistent at 0.3 combined
sigma; the band halved from G4's ±0.15.**  The direct estimator is now
stable under cap (1.997 at 1200 vs 1.996 extrapolated vs 2.00 at G4) and
its formal band is honest (the sigma model grew to match the measured
noise).  The residual spread is still protocol choice (direct vs fitted),
same class as G4.

b = 42 closed loop, reading B (running estimator crossing +42 from the
repulsive side): **x0 = 1.81 (endpoint) / 1.81 (extrapolated)** vs G4's
1.88 and the corpus's 1.90 ± 0.05.  This reading DEGRADED slightly: the
deeper-relaxed amplitude is larger, so the crossing moves down.  It sits
1.8 sigma below the corpus value — reported as-is.  (Reading A remains
uninformative for the reasons in the G4 record.)

Pair-law + core fit on the four 1200 points (dof = 1): C_d = 4.26e8, nu
at its lower bound 1.20 mu, chi2/dof = 1.87, b_eff = 7.5e8.  The relaxed
amplitude moved TOWARD the seed-ansatz value (G4: relaxed 7.9e6 vs seed
1.7e8) — consistent with the G4 diagnosis that the relaxed cap-400
amplitude was suppressed by under-convergence.  The ~7-order mismatch vs
the corpus's 42 is unchanged (G4 caveat 5 stands).

## What improved and what didn't

* IMPROVED: well depth significance 2.1 -> 5.0 sigma (extrapolated
  estimator; 2.6 sigma on the most conservative reading) and the well is
  now a two-point feature; x0 band ±0.15 -> ±0.07, central value 1.95,
  0.3 sigma from the corpus's 1.92 ± 0.08; the cap systematic of every
  G4 bond-zone number is now MEASURED instead of assumed.
* NOT improved: the b_eff amplitude normalization (still ~7 orders above
  42; deeper relaxation made it larger, not smaller); the b=42 running
  crossing moved from 1.88 to 1.81 (1.8 sigma below corpus); absolute
  energies still iter_cap-limited (slope ~ -2e-5/iter at 1200).

## Honest caveats (binding)

1. Every run still ends at iter_cap: 1200 is 3x deeper, not converged.
   The m=20 differential tail is only marginally contracting (r = 0.97);
   its extrapolated value carries the largest well-window error bar.
2. Single channel at cap 1200: the noise model rests on the m=23
   noise-dominated point + the extrapolation error bars, not on G4's
   three-channel scatter.  The conservative sigma_ref (2.3e-3) was chosen
   over the aggressive one (4.9e-4); with the aggressive choice the
   endpoint well would be >10 sigma — not quoted, deliberately.
3. The Aitken extrapolation assumes a single dominant exponential mode in
   the differential tail; the r = 0.40 contraction at m=18 supports this
   but 81 tail samples cannot exclude slower modes below the error bar.
4. Cap 1200 instead of the planned 1600 (restart + co-tenancy, coordinator
   sanctioned).  A cap-1600-vs-1200 shift was therefore NOT measured; the
   1200-vs-400 shift and the extrapolation bracket it instead.
5. All G4 caveats (guards in the B=2 sector, product-ansatz bias,
   central4-on-relaxed unusable, epistemic notice) carry over unchanged.
   Every number is a within-model measured fact at eps = 0.05.

## Deviations from the G5b plan (with reasons)

1. Cap 1200, not 1600: container restart killed the first (1600) farm at
   ~55% of batch 1; G5a's 4-thread n192_fr5 then owned the box, capping
   this farm at 3 concurrent processes.  1200 kept the 2-batch farm
   inside the wall-time budget; the cap-systematic analysis (400 vs
   1200) is intact.
2. 3-wide farm batches (not 4), banking order anchor-pair + well first —
   restart resilience.
3. The E_int cap extrapolation uses the DIFFERENTIAL series, not per-run
   absolute tail fits: absolute tails are still non-contracting at 1200
   (measured, table above) while the differential series is the quantity
   that actually enters E_int; this is the G4 sigma-model lesson
   (absolute descent cancels between runs) applied to extrapolation.

## Files

```
fs-gum-twoknot/
  farm_refine.sh                 G5b farm (3-wide, banking order, cap 1200)
  analyze_refine.py              G5b analysis (diff-series extrapolation)
  runs1200/                      6 per-run JSONs + logs (G4's runs/ intact)
  runs1200_farm_console.txt      farm wall-time record
  twoknot_refine_results.json    full numbers (per-run, tails, systematic,
                                 noise model, well, bond, fit)
  twoknot_refine_fig.png         E_int 400 vs 1200 vs extrap + diff series
  analysis_refine_console.txt    analysis stdout
  twoknot_REFINE.md              this file
```
