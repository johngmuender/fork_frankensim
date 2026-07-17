# L5′ — A3 arrival times, smoothed truncation + slower boost (F-T6-L5-EXEC-2)

**Goal:** L5′ is the pre-registered follow-up to L5 (F-T6-L5-EXEC), designed
in advance — no scanning — to convert L5's two spec-artifact PARTIALs:
(a) G2's single non-convergent trajectory born on the discontinuous ±3σ
truncation edge → fixed by a C²-smooth compact taper (S2 smoothstep, W = 1
for |x−x0| ≤ 2.5σ, 0 at 3.5σ); (b) G3's 2.8 % axial tail fraction (under the
5 % bar, compressed by convective transport at k0 = 4) → fixed by a slower
packet, k0 = 2, with windows T_near = 16 (near-field analysis) and T = 24
(far field; extended to T = 40 in the G4 coverage retry). Same physics,
same engines, same gate battery, full ensemble, no exclusions.

Within-model; nothing here bears on nature.

## Gates

| id | spec | measured | verdict |
|----|------|----------|---------|
| G1 | norm conserved to 1e-6; wall amplitude < 1e-10; energy drift ≤ 1e-6 | norm dev 2.2e-16; wall max 2.8e-16 (diag) / 1.1e-15 (traj run); energy drift 3.7e-16 (rel); z-mode leakage 0.0; 1-D energy now grid-converged: E(Nx=2048/4096/8192) = 2.2503989 ± 5e-9 (L5's log-divergent truncation energy is gone) | **PASS** |
| G2 | transverse near-field hard cutoff: zero arrivals in (τ_max, T], gap ≥ 2·IQR, τ_max stable ≤ 2 % under double-Nx AND halve-dt separately, full ensemble incl. every edge trajectory; sample-max trajectory convergent over ≥ 3 resolutions | τ_max = 5.13095; zero arrivals in (5.131, 16] and in (5.131, 40] on the retry horizon; gap 10.869 ≥ 2·IQR = 2.432; stability 8.3e-8 (Nx×2), 5.6e-10 (dt/2); all 2000 trajectories converge (max per-trajectory |Δτ| 3.7e-5 under Nx×2, 1.7e-7 under dt/2, zero arrival-status changes); sample max = idx 1754, relative spread 8.6e-7 across the 6-run ladder | **PASS** |
| G3 | axial arrivals at d_near strictly beyond transverse τ_max: fraction > 5 % | 219/1976 = **11.08 %** (window T = 16; 11.58 % on T = 24), axial support to 15.79 (window) / 23.39 (T = 24) vs τ_max = 5.13; zero transverse arrivals in that region | **PASS** |
| G4 | far-field two-sample KS (transverse vs axial) at d_far = 25: p > 0.05 | T = 24 window: stat 0.0486, p = 0.025 — but coverage was asymmetric (non-arrivals 30 vs 290), so one coverage retry (T = 40, box [−15, 225], same dx): balanced coverage (0 vs 88), stat 0.0925, **p = 9.9e-8** (dt-stable to 9 digits). KS(d) decays monotonically 0.126 (d = 1) → 0.057 (d = 35, p = 3.7e-3): the null is being approached, not violated; d = 25 is not asymptotic at k0 = 2 | **FAIL-AT-d25** (physics of configuration, honest) |
| G5 | (report-only) mean arrival ratio transverse/axial at d_near, bootstrap error | 0.8569 ± 0.0119 (means 2.9589 vs 3.4530; transverse earlier at ≈ 12σ from 1; L5: 0.9767 ± 0.0063) | REPORT-ONLY |

**L5's two PARTIALs both convert to PASS. The pre-registered configuration
additionally exposes a real finite-distance effect at G4 (below).**

## L5 vs L5′ comparison

| quantity | L5 (k0 = 4, hard cut) | L5′ (k0 = 2, smooth taper) |
|---|---|---|
| τ_max (transverse cutoff) | 2.28503 | 5.13095 |
| empty gap (T_window − τ_max) | 5.715 (≥ 2·IQR = 0.764) | 10.869 (≥ 2·IQR = 2.432); holds to horizon 40 (gap 34.87) |
| τ_max stability (dt/2, Nx×2) | 7.1e-7, 1.6e-4 | 5.6e-10, 8.3e-8 |
| edge trajectory | idx 1249 non-convergent over 7 resolutions (arrivals 1.67…7.34/never) — G2 PARTIAL | sample max idx 1754 converged, rel spread 8.6e-7; **failure mode gone** |
| axial fraction beyond τ_max | 2.80 % (< 5 %) — G3 PARTIAL | **11.08 %** (> 5 %) — PASS |
| far-field KS at d = 25 | stat 0.036, p = 0.307 — PASS | stat 0.0925, p = 9.9e-8 (balanced coverage) — FAIL-AT-d25 |
| KS(d) trend | not measured | 0.126 → 0.057 monotone over d = 1…35, p rising to 3.7e-3 |

**Bench-design tradeoff (recorded for the Layer-1 envelope):** slowing the
boost sharpens the near-field spin contrast (G3: 2.8 % → 11.1 %; G5 ratio
0.98 → 0.86) but pushes the asymptotic far-field regime outward — at k0 = 2
the detector at d = 25 sits effectively nearer in packet-width units, so the
same d that satisfied the far-field null at k0 = 4 no longer does. Near/far
detector placement must be co-designed with the boost.

## Key numbers

- Transverse near-field cutoff: **τ_max = 5.130950** (reference: exact-in-z
  engine, dx = 0.0293 ≡ Nx = 2048 on the L5 box, dt = 2.5e-4).
  Stability: 5.1309500 at Nx×2 (8.3e-8 rel), 5.1309505 at dt/2 (5.6e-10 rel);
  5.1309462 on the pre-registered box (8.3e-7); 5.130417 on the independent
  2-D engine at matched coarse resolution (1.0e-4 rel).
- Full-ensemble convergence (all 2000 transverse trajectories, no
  exclusions): per-trajectory |Δτ| median/p99/max = 2.5e-6 / 2.2e-5 / 3.7e-5
  under Nx×2 and 1.8e-9 / 2.7e-8 / 1.7e-7 under dt/2; arrival-status changes: 0.
- Axial arrivals beyond τ_max: 219/1976 = 11.08 % (window 16); axial last
  arrival 15.787 (window) / 23.389 (T = 24, dt-stable to 1e-8); zero
  transverse arrivals beyond τ_max out to t = 40.
- IQR of transverse near-field arrivals: 1.2159; empty gap 16 − τ_max = 10.869.
- Non-arrivals (of 2000): near-field window 16 — transverse 0, axial 24;
  by T = 24 — 0 / 13; by T = 40 — 0 / 8. Far field d = 25 by T = 24 —
  30 / 290 (the coverage asymmetry that triggered the retry); by T = 40 —
  0 / 88. Frozen box escapees (reference run): 1 transverse, 0 axial.
- Far-field KS at d = 25, T = 40, balanced coverage: stat 0.09247,
  p = 9.9e-8 (identical to 9 digits at dt = 5e-4). KS(d), T = 40:

  | d | 1 | 10 | 15 | 20 | 25 | 30 | 35 |
  |---|---|----|----|----|----|----|----|
  | KS stat | 0.1261 | 0.1161 | 0.1179 | 0.1049 | 0.0925 | 0.0788 | 0.0571 |
  | p | 2.7e-14 | 3.9e-12 | 1.9e-12 | 6.7e-10 | 9.9e-8 | 1.1e-5 | 3.7e-3 |
  | non-arr T/A | 0/8 | 0/24 | 0/42 | 0/60 | 0/88 | 0/118 | 0/160 |

- Cross-engine validation (2-D grid engine vs exact-in-z engine, matched
  Nx, dt = 2e-3, T = 24): median |Δτ| = 1.9e-4, p99 = 2.7e-3, KS p ≈ 1.0;
  far-field arrival counts and KS statistic bit-identical across engines
  (same 1970/1710 samples, same rank configuration).
- G5 ratio T/A at d_near: 0.8569 ± 0.0119 (bootstrap, 2000 resamples).
- 2-D total energy 7.18520 = E_x (2.25040) + π²/2; z-mode leakage exactly 0.

## Method summary

Identical to L5 except the two pre-registered changes. 2-D (x, z)
waveguide, hard walls at z = 0, 1 (odd extension, sine basis), exact
spectral propagation (G1 at machine precision by construction, and
measured). Initial state: Gaussian × S2-taper in x — amplitude
exp(−(x−x0)²/(2σ_x²))·W(x)·e^{ik0x}, σ_x = 1, x0 = −5, k0 = 2, compact
support |x−x0| ≤ 3.5σ retained — × ground mode sin(πz). Pauli guidance as
in L5 (transverse σ_y: j_spin = (−½∂_zρ, +½∂_xρ); axial: j_spin = 0).
N = 2000 from |φ(0)|² (analytic inverse-CDF, grid-independent, seed
20260717), identical ensemble for both spin cases, all resolutions, and
both engines. Primary engine: factorized exact-in-z (validated in L5;
revalidated here). RK4 trajectories on stored velocity fields at
t, t+dt/2, t+dt; arrival = first crossing of x = d.

Reference: dx = 60/2048 (pre-registered resolution) on an extended box
[−15, 105] (Nx = 4096), dt = 2.5e-4, T = 24; refinements Nx = 8192 and
dt = 1.25e-4 run separately; ladder of 6 resolutions total for the
edge-trajectory clause. The extended box is a containment fix: at k0 = 2
the packet front reaches x ≈ 93 by T = 24 and would wrap the periodic
FFT box [−15, 45]; a control run on the pre-registered box gives
τ_max identical to 8.3e-7 and the same G3 fraction (11.26 %), so no gate
conclusion depends on the box. G4 coverage retry: T = 40 on box
[−15, 225] (Nx = 8192, same dx), seven detector planes d = 1…35.

## Realization

Both L5 spec artifacts are confirmed as artifacts: with the C² taper the
truncation-edge diffraction ringing is gone — every one of the 2000
trajectories, including the latest arrival (idx 1754, born mid-packet at
x0 = −5.22, z0 = 0.222, not on the support edge), converges across the
resolution ladder to ≤ 3.7e-5 in |Δτ| — and with k0 = 2 the axial late
tail widens to 11.08 %, comfortably over the 5 % bar. The transverse hard
cutoff itself is strengthened: the arrival density ends at τ_max = 5.131
and stays empty for the following 34.9 time units of evolution, while the
axial density continues smoothly past it — the Das–Dürr mechanism
(the σ_y spin current advances the slow half of the guide and removes the
would-be-late arrivals) realized cleanly at slow boost.

G4 is the honest cost of the same pre-registered choice: with balanced
coverage the far-field distributions at d = 25 still differ (KS 0.0925,
p = 9.9e-8), and the monotone KS(d) decay (0.126 → 0.057 across d = 1…35,
p climbing seven orders of magnitude to 3.7e-3) shows the far-field null
being *approached, not violated* — at k0 = 2 the detector at d = 25 is
effectively nearer in packet widths, i.e. not yet asymptotic. This is a
physics-of-configuration FAIL at the pre-registered d, not a numerical
artifact (dt-stable to 9 digits; both engines bit-consistent) and not a
contradiction of L5's k0 = 4 PASS at the same d.

**Connection.** Same S3/A3 observable of GUM Sec. IX.A (Layer 1) as L5;
the hard cutoff is the non-POVM candidate statistic, now with its two
harness defects discharged. Filed as F-T6-L5-EXEC-2.

## Caveats

- "This reproduces the phenomenon; it neither proves nor tests
  POVM-exclusion (the corpus's T4-W5 theorem gap) - it builds the testbed
  any such argument must clear."
- The pre-registered amplitude formula exp(−(x−x0)²/(2σ_x²)) was followed
  literally; it gives position-density std σ_x/√2 ≈ 0.707, narrower than
  L5's code convention exp(−(x−x0)²/(4σ_x²)) (std σ_x). A sensitivity run
  with the L5 convention (+ taper) gives τ_max = 4.457, axial-beyond
  11.66 %, far-field KS p (T = 24 window) = 0.004: G2/G3 conclusions are
  robust to the convention; the G4 finite-distance failure is too.
- G4's verdict is bound to the pre-registered d_far = 25: KS(d) indicates
  p > 0.05 would be reached at larger d (extrapolating the measured trend),
  but no such run was scored — the gate is graded FAIL at its own spec.
- τ_max is the endpoint of the sampled ensemble (N = 2000), not a proved
  support edge of the continuum density.
- The x-domain is periodic (FFT); reference runs use an extended box for
  containment (see method). On the pre-registered box the wrapped fast tail
  re-enters and artificially pushes the 13 slowest axial trajectories
  across d_near (2000 vs 1987 arrivals) — gate values are unaffected
  (control deltas: τ_max 8.3e-7, G3 +0.2 pp, far KS p 0.079 vs 0.025 on the
  T = 24 window), but the extended box is the meaningful configuration.
- Far-field KS compares arrived sub-populations (2000/1912 at T = 40; the
  88 axial stragglers are the slow tail, recorded honestly); both ensembles
  share the identical initial points.
- Both σ_y eigenstates give mirror trajectories under z → 1 − z with
  identical statistics (z-symmetric initial state); M_y = +ρ/2 w.l.o.g.
