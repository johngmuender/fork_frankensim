# fs-gum-twoknot — GUM core Phase G4 results (App I.2 two-knot bond-equation flagship)

Crate: `analysis/gum-sandbox/gum-core/fs-gum-twoknot/` (standalone
`[workspace]` opt-out; path deps `../fs-gum-field`, `../fs-gum-statics`,
`crates/fs-math` — neither finished crate was modified).

Scope: the FIRST execution anywhere in this repo of the corpus's App I.2
protocol (01-GUM-Omega-Paper-v2.0.1.md line 714): two B = 1 knots at a scan
of separations, product-ansatz seeding, three relative orientations
(Theorem VII.2 channel theorem), guarded ANF relaxation at the frozen
eps = 0.05 dial, interaction energy E_int(d), the screened-dipole pair-law
fit, and the bond-equation closure against x0 = 1.90 ± 0.05 (bond equation
at 𝔟 = 42) / 1.92 ± 0.08 (corpus-measured).

Reproduce: `cargo build --release && ./target/release/twoknot_gates &&
sh farm.sh runs 400 && python3 analyze_twoknot.py` (~50 min wall on 4
cores for the farm; each run is an independent OS process and is
bit-deterministic for fixed args).

## Protocol as run

```
grid        160 x 96 x 96 cells at the FROZEN h = 0.09375 (box 15 x 9 x 9);
            transverse box = the certified single-knot referee geometry
            (LBOX = 4.5, N = 96); long axis holds the x ~ 4 far anchor with
            the full ln(1e6)/mu = 1.78 Yukawa margin.  1,474,560 cells.
            [Deviation from the survey's 144 x 96 x 96: +16 long-axis cells
            (+11% cost) so the asymptote anchor fits ON THE SAME GRID CLASS
            as every other run — see Deviations.]
engine      anisotropic (nx,ny,nz) twin of Field3 + eval + ANF carried in
            THIS crate (fs-gum-field's Field3 is cubic-only), gated BITWISE
            against fs-gum-statics on cubic grids before any campaign run
            (gate binary twoknot_gates, 20/20 PASS incl. two-run replay:
            eval Out + full gradient bitwise, 25 guarded ANF iterations
            bitwise in series and field).
seeding     product ansatz q(x) = q1(x - d/2 X) * q2'(x + d/2 X) (Hamilton
            product of two radial-profile hedgehogs; radial_solve(t_frozen,
            4000, 6.0)); knot 2 isorotated q2' = a q2 conj(a):
              attract  a = (0,0,0,1)   theta = pi about z-hat PERP X (pi/perp rule)
              repulse  a = (0,1,0,0)   theta = pi about x-hat PAR X
              align    a = identity
            centres +-(m h) sit on cell corners for EVERY separation and for
            the single-knot reference: one sub-grid alignment class, so
            per-knot discretisation offsets cancel in every energy
            difference.  Seed degree ~ 2 verified (1.9725 = 2 x the
            single-knot discrete degree 0.98601 at this grid class).
relax       guarded ANF, corner objective, t = 0.008276434949296802,
            deg_ref = seed degree (~1.9725: the anchor at 2), deg_unit =
            deg_ref/2 (per-unit-degree Bogomolny wall — the B = 2 guard
            generalization), fgap_ref = seed floor gap - 0.01, frozen
            MU = 5000 / MU_F = 400 / bands 0.005 / 0.01; dt0 = 0.01,
            dt_max = 0.05; iteration cap 400 UNIFORM (all runs end at
            iter_cap; tail-flatness recorded per run).
separations d = 2 m h, m in {14,16,18,20,23,26,29} -> x = d/R* in
            {1.4732, 1.6837, 1.8942, 2.1046, 2.4203, 2.7360, 3.0517}
            (the mandatory calibration window [1.5, 3] + one bracket point
            below), + far anchor m = 38 (x = 3.9988) in the attractive
            channel, + one single-knot reference on the same grid.
farm        23 runs as 4 parallel OS processes (farm.sh, 6 batches);
            attractive channel + single + far anchor FIRST.
```

## What was run (complete — no budget reduction was needed)

| block | runs | iters | status |
|---|---|---|---|
| single-knot reference | 1 | 400 | iter_cap |
| attract (x = 1.47..3.05 + far 4.0) | 8 | 400 | iter_cap |
| align | 7 | 400 | iter_cap |
| repulse | 7 | 400 | iter_cap |

Wall time: farm 49.5 min on 4 cores (sum of per-run ANF times 11,023 s =
3.06 h single-core-equivalent; effective multiprocess speedup ~3.7x,
consistent with the survey's measured 3.96x).  Contended rate ~1.4 s/iter
per process (1.135 s/iter solo pilot) at 1.47 M cells — within the
survey's 0.59 us/cell/iter extrapolation band.  Gates 3.6 s; pilot 12.5 s.

## Gate table (twoknot_gates, 20/20 PASS, two-run replay bitwise)

| id | gate | verdict |
|----|------|---------|
| T1 x4 | eval twin vs fs-gum-statics, cubic N=24, corner+central4, bare+guards-forced: all Out fields AND all 55,296 gradient components bitwise | PASS |
| T2 | 25 guarded ANF iterations: (obj, estat, deg, gnorm, dt, arrests) series bitwise + final field bitwise | PASS |
| T3 x4 | product-ansatz degree ~ 2 all three channels (1.9709-1.9725); single ~ 1 (0.98601) | PASS |
| T4 | full gate suite two-run in-process replay, 15 fingerprint f64s bitwise | PASS |

## Single-knot reference (same grid class, same protocol)

E_corner(seed hedgehog) = 3.146623554, E_corner(relaxed 400) = 3.115089428
(drift -0.0315, the documented near-BPS guard-bounded drift class), final
degree 0.98046 (anchor band edge), 7 arrests.

## Per-run table (corner scheme; E_int = E(d) − E(far anchor x=3.9988))

`Eint_seed` = product-ansatz (iteration-free) estimator E_seed(d) − 2 E_seed(single).
`sigma` = measured differential-noise model (below).  deg_f = final degree;
d_eff = (1−q0)-weighted half-box centroid separation (seeded d in braces).

| ch | x = d/R* | d | E_int (relaxed, anchored) | E_int (seed) | sigma | deg_f | d_eff final |
|----|------|-------|-------------|-------------|--------|--------|--------|
| attract | 1.4732 | 2.6250 | +9.7389e-01 | +9.7602e-01 | 9.7e-2 | 1.96179 | 2.587 {2.590} |
| attract | 1.6837 | 3.0000 | +1.8422e-01 | +1.9686e-01 | 1.8e-2 | 1.96604 | 2.963 {2.967} |
| attract | 1.8942 | 3.3750 | **-5.2665e-04** | +4.1181e-03 | 2.5e-4 | 1.96690 | 3.355 {3.364} |
| attract | 2.1046 | 3.7500 | +2.2217e-04 | -3.7430e-04 | 2.5e-4 | 1.96655 | 3.745 {3.750} |
| attract | 2.4203 | 4.3125 | -2.2501e-04 | -4.8266e-06 | 2.5e-4 | 1.96648 | 4.312 {4.312} |
| attract | 2.7360 | 4.8750 | -6.8046e-07 | -5.8613e-08 | 3.2e-6 | 1.96650 | 4.875 {4.875} |
| attract | 3.0517 | 5.4375 | -1.1407e-08 | -7.1978e-10 | 4.0e-8 | 1.96650 | 5.437 {5.437} |
| attract | 3.9988 | 7.1250 | 0 (anchor)  | -1.5730e-12 | —      | 1.96650 | 7.125 {7.125} |
| align   | 1.4732 | 2.6250 | +1.8579e+00 | +1.8632e+00 | 1.9e-1 | 1.95872 | 2.745 |
| align   | 1.6837 | 3.0000 | +5.2459e-01 | +5.2651e-01 | 5.2e-2 | 1.96259 | 3.050 |
| align   | 1.8942 | 3.3750 | +2.8150e-02 | +3.5719e-02 | 2.8e-3 | 1.96556 | 3.398 |
| align   | 2.1046 | 3.7500 | -2.5764e-04 | +5.1874e-04 | 2.5e-4 | 1.96645 | 3.755 |
| align   | 2.4203 | 4.3125 | +3.1684e-04 | +4.9284e-06 | 2.5e-4 | 1.96648 | 4.313 |
| align   | 2.7360 | 4.8750 | +6.2417e-07 | +5.9181e-08 | 3.2e-6 | 1.96650 | 4.875 |
| align   | 3.0517 | 5.4375 | +1.0698e-08 | +7.1892e-10 | 4.0e-8 | 1.96650 | 5.437 |
| repulse | 1.4732 | 2.6250 | +2.1444e-01 | +2.5177e-01 | 2.1e-2 | 1.95882 | 2.819 |
| repulse | 1.6837 | 3.0000 | +3.1736e-02 | +4.3739e-02 | 3.2e-3 | 1.96276 | 3.077 |
| repulse | 1.8942 | 3.3750 | +2.3846e-03 | +4.0582e-03 | 2.5e-4 | 1.96538 | 3.400 |
| repulse | 2.1046 | 3.7500 | +3.1292e-04 | +3.6059e-04 | 2.5e-4 | 1.96639 | 3.755 |
| repulse | 2.4203 | 4.3125 | +1.9930e-04 | +4.7347e-06 | 2.5e-4 | 1.96649 | 4.313 |
| repulse | 2.7360 | 4.8750 | +7.3638e-07 | +5.8039e-08 | 3.2e-6 | 1.96650 | 4.875 |
| repulse | 3.0517 | 5.4375 | +1.2115e-08 | +7.1353e-10 | 4.0e-8 | 1.96650 | 5.437 |

Anchor cross-check: E(far, relaxed) − 2 E(single, relaxed) = **+2.657e-2**
(corner).  This is the cap-convergence + guard-geometry systematic between
the B = 2 and B = 1 protocols (the pair's degree anchor holds deg/2 =
0.98325 per knot vs the single knot's 0.98046) — it is why E_int is
defined against the far anchor, which shares the pair protocol exactly.
Both definitions are in `twoknot_results.json`.

Differential-noise model (MEASURED, not assumed): at m = 20–23 the three
channels' anchored E_int scatter by ~2.5e-4 with random signs while the
seed estimator bounds the physical interaction there at <= 5e-6; that
scatter is the differential cap noise.  It decays with separation like the
interaction (m = 26: ~7e-7, m = 29: ~1.1e-8), so sigma_dc(d) =
2.52e-4 · min(1, e^{−mu (d − d_23)}), floored at 10% relative.

## The pair law and the channel theorem (measured)

Sign structure — exactly the channel theorem's: the pi/perp channel is the
ONLY attractive one; aligned and pi/parallel are repulsive.  Tail-ratio
test on the clean seed estimator (predictions from Thm VII.2 at x = mu d):

| x = d/R* | align/(−attract) meas | pred | repulse/(−attract) meas | pred |
|------|------|------|------|------|
| 2.42 | 1.021 | 0.942 | 0.981 | 1.058 |
| 2.74 | 1.010 | 0.949 | 0.990 | 1.051 |
| 3.05 | 0.999 | 0.954 | 0.991 | 1.046 |

The near-unity magnitudes are confirmed at the 2–8% level; the predicted
±5% align-vs-repulse splitting is NOT resolved (it sits at the
product-ansatz bias level).

Yukawa mass of the measured tail (free-mu refit of the pair-law shape on
the seed estimator's pure-attraction points): mu_eff = 7.553 (corner) /
7.807 (central4) vs the frozen radial-referee mu = 7.7725 — the
interaction tail decays with the single-knot Yukawa mass to ±3%.

## Fit (form stated verbatim) and the b_eff amplitude

Fitted form on the calibration window x in [1.5, 3] (+ the 3.05 point),
attraction FIXED to the corpus pair-law shape (7.5) with x = mu r, r = d:

    E_int(d) = − C_d · h_par(mu d)  +  B_core · e^{−nu d},
    h_par = (2 + 2x + x²) e^{−x} / d³          [paper (7.5), verbatim]

nu restricted to [1.2 mu, 2.6 mu] (nu -> mu is collinear with h_par and
produced a degenerate fit — observed and rejected; the corpus's own
"derived core repulsion" closed form is not printed in the paper text, so
the exponential class is an assumption of this analysis and nu is
reported).  Corner scheme (the descent's own instrument; central4 numbers
on the RELAXED states are unusable — see Caveats):

| estimator | C_d | nu | chi2/dof | b_eff = C_d mu/(2 pi p²) |
|---|---|---|---|---|
| relaxed anchored | 7.87e6 ± 3.3e6 | 20.2 (upper bound) | 1.9/3 | **1.38e7** |
| seed ansatz | 1.721e8 ± 0.8e6 | 9.33 (lower bound) | 5409/3 | **3.02e8** |

The running estimator b_eff(x) = −E_int·mu/(2 pi p² h_par) plateaus at
~1.3e8–2.2e8 across the window (seed).  **Compared to the corpus's
𝔟_eff = 42 ± 6 the measured amplitude is ~7 orders of magnitude larger
under the App-A normalization as reconstructed here (b = B/2 pi p²,
B = C_d mu, p = 0.84).**  At the benchmark point the knots are near-BPS
compactons with mu R* = 13.85: at bond-zone separations the cores overlap
and the linear-zone multipole normalization is enhanced by an
e^{+mu(core scale)}-class factor.  The corpus's 42 is not recoverable from
the visible text at this dial; what IS confirmed is the pair-law SHAPE
(Yukawa mass, channel signs, near-unity channel ratios).  Both C_d and the
running-b_eff tables are in twoknot_results.json.

## The bond: x0 (all estimators, corner scheme)

The attractive channel crosses from core repulsion to tail attraction and
forms a shallow bond.  Measured well: E_int(x = 1.894) = −5.27e-4 ±
2.5e-4 (2.1 sigma below zero; the only > 2 sigma negative point), depth
~1.7e-4 of the two-knot energy.

| estimator | x0 |
|---|---|
| direct parabolic minimum, relaxed anchored (primary) | **2.00** |
| model-fit argmin, relaxed anchored | 1.85 |
| model-fit argmin, seed ansatz | 2.04 |
| direct parabolic minimum, seed ansatz | 2.25 |
| **spread of estimators (quoted band)** | **x0 = 2.0 ± 0.15** |
| corpus measured | 1.92 ± 0.08 |
| corpus bond equation at 𝔟 = 42 | 1.90 ± 0.05 |

**x0(measured here) = 2.0 ± 0.15 vs corpus 1.92 ± 0.08 — consistent
within 0.5 combined sigma; the bond exists, in the attractive channel
only, at the corpus's location.**  The formal parabola-propagation band on
the primary estimator alone ([1.998, 1.999]) is far tighter than the
honest estimator spread, which is dominated by protocol choices (relaxed
vs seed, direct vs fitted) — the ±0.15 spread is the quoted uncertainty.

Bond-equation closure at 𝔟 = 42 (two readings, since App I.2 states the
protocol in one line):

* Reading A (amplitude substitution): argmin of −C_d(42) h_par + fitted
  core, C_d(42) = 2 pi p² · 42/mu = 23.96 -> x0 = 2.44 (relaxed-fit core);
  no interior minimum with the seed-fit core.  Under our normalization
  C_d(42) is ~7 orders below the measured amplitude, so this reading
  mostly probes the fitted core term — reported for completeness, weight
  it accordingly.
* Reading B (running estimator): the x where b_eff(x) crosses +42 from
  the repulsive side — **x0 = 1.88 (relaxed) / 1.97 (seed)** vs the
  corpus's 1.90 ± 0.05.  Because the measured plateau amplitude is >> 42,
  this reading lands at the repulsion-attraction crossover, i.e. the bond
  onset; it is the reading that reproduces the corpus number, and it does
  so essentially independently of the amplitude mismatch.

## Honest caveats (binding)

1. **Grid class / cap effects.** All runs end at iter_cap = 400 (never
   converged); E_int is meaningful only DIFFERENTIALLY at fixed protocol.
   The far-anchor definition cancels the shared systematic; the residual
   differential noise is the measured sigma_dc above.  The well depth is a
   2.1 sigma effect against that noise.  A longer-cap follow-on would
   sharpen it.
2. **The 2 E_single cross-check offset (+2.66e-2)** is real and understood
   (per-knot degree pinned differently by the B = 2 vs B = 1 anchor), and
   is why E_int(d) = E(d) − 2 E_single is NOT the primary estimator here.
3. **Product-ansatz bias.** Seeds are unsymmetrized Hamilton products
   (q1 q2 != q2 q1 pointwise); the seed estimator carries O(overlap²)
   ansatz bias in the window (visible as the 2–8% channel-ratio residuals
   and the relaxed-vs-seed differences at x <= 2.1).
4. **Central4 on relaxed states is unusable** (differential rms ~5e-2):
   the corner-descent lattice states sit in corner-scheme discretization
   valleys that the 4th-order instrument penalizes; central4 agrees with
   corner on the SEED estimator (mu_eff 7.81, C_d 1.48e8).  Corner is the
   protocol instrument.
5. **b_eff normalization.** The App-A reconstruction (b = C_d mu/2 pi p²)
   is this analysis's dimensional reading of one sentence; the ~7-order
   mismatch vs 42 says the corpus's calibration-zone normalization is not
   recoverable from the visible text at eps = 0.05, NOT that the pair law
   fails (its shape checks pass).  The core-repulsion exponential class
   (and both nu's landing on their bounds) is likewise an assumption of
   this analysis, printed as part of the estimator definition per the
   survey's risk note.
6. **Orientation coverage.** All three App I.2 orientations were run at
   all 7 scan separations (full protocol); the clock-coherent locked bond
   (x0_lock = 2.42 ± 0.12, release r6) is OUT of scope — it needs the
   cos(dPhi) transverse channel, i.e. isorotating pairs.
7. **Guards in the two-charge sector.** The frozen guard parameters were
   never calibrated for B = 2 (survey risk); the degree anchor was active
   at endpoints (epen ~ 7e-4 to 1e-3, deg drift −0.006 from seed, same
   class as the certified single-knot runs) and the deg_unit = deg_ref/2
   floor generalization is this crate's documented adaptation.
8. **Epistemic notice.** Every number is a WITHIN-MODEL measured fact
   about a speculative theory's functional on a lattice at eps = 0.05.
   Nothing here validates the theory; adjudication is the coordinator's.

## Deviations from the roadmap/survey spec (with reasons)

1. **Grid 160 x 96 x 96 instead of 144 x 96 x 96** (+11% cells): the
   mandated far-point anchor (x ~ 4) does not fit the 144 box with the
   1e-6 tail margin; one uniform grid class for every run (including the
   anchor and the single-knot reference) buys exact cancellation of the
   discretisation offset — decisive for a 1e-4-scale observable.  h, the
   transverse box, and all frozen constants are unchanged.
2. **Non-cubic support via an in-crate anisotropic twin** (Field3 is
   cubic-only and was not modified): cubic 160³ would have cost 2.8x the
   cells (~2.3 h farm instead of ~50 min).  The twin is not trusted on
   faith: it reproduces fs-gum-statics BITWISE on cubic grids (gates T1,
   T2) — the strongest available equivalence.
3. **Iteration cap 400** (spec suggested starting there; ~1500 was the
   survey's full-campaign guess): tail-flatness was recorded and the
   differential noise measured (2.5e-4 mid-zone); the full 3-orientation
   scan fit in budget at 400 with the bond resolved at 2.1 sigma.  Longer
   caps are the natural follow-on, not a blocker.
4. **deg_unit floor generalization** for B = 2 (deg_ref/2), documented
   above — the frozen single-knot formula would have halved the
   per-degree Bogomolny wall.
5. **sigma model is empirical** (cross-channel scatter), not the a-priori
   per-run tail-flatness proxy, which overestimates DIFFERENTIAL error by
   ~70x (it measures absolute descent, which cancels between runs).

## Files

```
fs-gum-twoknot/
  Cargo.toml
  src/lib.rs            constants, protocol docs
  src/field_a.rs        anisotropic FieldA + stencils (Field3 twin)
  src/engine_a.rs       eval_a (objective + analytic gradient, deg_unit floor)
  src/anf_a.rs          arrested Newton flow (frozen protocol)
  src/seed.rs           product ansatz, channels, d_eff diagnostic
  src/bin/twoknot_gates.rs   T1-T4 bitwise twin gates
  src/bin/twoknot_run.rs     one campaign run -> runs/run_<ch>_m<m>.json
  farm.sh               4-process task farm (6 batches, attract first)
  analyze_twoknot.py    deterministic analysis -> twoknot_results.json + fig
  runs/                 23 per-run JSONs + logs (+ runs_pilot/: 10-iter pilot)
  twoknot_results.json  full numbers (per-run, fits, bond closures, ratios)
  twoknot_fig.png       channels + bond-region figure
  analysis_console.txt  analysis stdout (the summary tables)
```
