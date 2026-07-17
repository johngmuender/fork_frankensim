# L5 — A3 arrival-time phenomenon on an independent engine (F-T6-L5-EXEC)

**Goal:** first in-repo instantiation of the GUM Layer-1 core phenomenon
(Das–Dürr class, Sci. Rep. **9**, 2242 (2019)): spin-dependent Bohmian arrival
times in a 2-D hard-wall waveguide. Transverse (out-of-plane, σ_y) spin turns
on the Pauli spin current and modifies trajectories; in-plane (axial) spin
gives pure convective (flux) trajectories. Measure the near-field hard support
cutoff, the axial contrast, and the far-field null.

Within-model; nothing here bears on nature.

## Gates

| id | spec | measured | verdict |
|----|------|----------|---------|
| G1 | norm conserved to 1e-6; wall amplitude < 1e-10; energy drift ≤ 1e-6 | norm dev 2.2e-16; wall max 9.4e-16; energy drift 2.7e-16 (rel); z-mode leakage 0.0 | **PASS** |
| G2 | transverse near-field hard cutoff: zero arrivals in (τ_max, T], T − τ_max ≥ 2·IQR, τ_max stable ≤ 2% under one refinement | τ_max = 2.28503; zero arrivals in (2.285, 8] at reference resolution; gap 5.715 ≥ 2·IQR = 0.764; stability 7.1e-7 (dt/2), 1.6e-4 (Nx×2) for 1999/2000 trajectories — but ONE truncation-edge trajectory (0.05 %) is numerically non-convergent across 7 resolutions (arrives at 5.69/6.13/3.10/1.67/7.34/never/6.79), so the raw sample maximum does not satisfy the stability clause as literally written | **PARTIAL** |
| G3 | axial arrivals beyond transverse τ_max: fraction > 5 % | contrast present, expected direction: axial support extends to 3.2245 (41 % past the transverse cutoff in time); fraction of axial arrivals beyond τ_max = 2.8 % (56/2000), zero transverse arrivals there — below the 5 % threshold | **PARTIAL** |
| G4 | far-field two-sample KS (transverse vs axial) p > 0.05 | KS stat 0.0362, p = 0.307 (n_T = 1421, n_A = 1381); independent 2-D engine: p = 0.332 | **PASS** |
| G5 | (report-only) ratio of mean arrival times transverse/axial at d_near | 0.9767 ± 0.0063 (bootstrap, 2000 resamples); means 1.4884 vs 1.5240 | REPORT-ONLY |

## Key numbers

- Transverse near-field cutoff: **τ_max = 2.28503** (reference resolution
  Nx = 2048, dt = 2.5e-4, exact-in-z engine).
  Stability: 2.2850273 → 2.2850289 under dt/2 (7.1e-7 rel); 2.28467 at
  Nx = 4096 (1.6e-4 rel); 2.2871 on the independent 2-D engine at Nx = 1024
  (9e-4 rel).
- Axial last arrival: 3.22447 (dt-stable to 1.6e-5; 3.2430 at Nx = 4096;
  3.3072 at Nx = 1024).
- IQR of transverse near-field arrivals: 0.3819; empty gap T − τ_max = 5.715.
- Axial arrivals beyond τ_max: 56/2000 = 2.80 %. Transverse arrivals beyond
  τ_max at reference resolution: 0.
- Non-arrivals within T = 8 at d_near: transverse 3/2000 (turned-back
  trajectories; includes the non-convergent edge particle), axial 0/2000.
  At d_far: 579/2000 (T), 619/2000 (A) — slow ensemble tail, recorded honestly.
- Far-field KS: stat 0.0362, p = 0.307.
- Mean arrival ratio T/A at d_near: 0.9767 ± 0.0063 (transverse mean earlier,
  ≈ 3.7σ from 1): the spin current advances the z > L/2 half of the ensemble
  more than it delays the z < L/2 half that still arrives.
- Cross-engine validation (2-D grid engine vs exact-in-z engine, matched
  Nx = 1024, dt = 2e-3): median |Δτ| = 3.6e-5, p99 = 1.2e-3, KS p = 1.0.

## Method summary

2-D (x, z) waveguide, hbar = m = 1, hard walls at z = 0, L = 1 via odd
extension in z and full 2-D FFT (sine basis); x ∈ [-15, 45], Nx = 1024,
Nz = 128. Initial state: truncated Gaussian in x (σ_x = 1, cut at ±3σ,
renormalized, x0 = −5, boost k0 = 4) × ground mode sin(πz). The free
Hamiltonian is diagonal in this basis, so the field is propagated **exactly**
(phase multiplication per trajectory half-step) — G1 quantities are at
machine precision by construction, and measured. Pauli guidance:
j = Im(φ*∇φ) + j_spin with j_spin = (−½ ∂_z ρ, +½ ∂_x ρ) for the σ_y
eigenstate (M_y = +ρ/2) and j_spin = 0 in-plane for axial spin; v = j/ρ.
N = 2000 initial points sampled from |φ(0)|² (seed 12345), RK4 (dt = 2e-3;
refinement 1e-3) with bilinear interpolation of stored velocity fields at
t, t+dt/2, t+dt; arrival = first crossing of x = d (d_near = 1, d_far = 25),
identical ensemble for both spin cases; T = 8.

G2 retry (resolution): the measured z-mode leakage is exactly 0 (pure ground
sine mode, diagonal evolution), so φ = f(x,t) sin(πz) e^{−iπ²t/2} exactly and
the guidance law reduces to dx/dt = Im(f′/f) − π cot(πz), dz/dt = Re(f′/f)
(transverse) and dx/dt = Im(f′/f), dz/dt = 0 (axial). A factorized engine
with analytic z-factors (no near-wall interpolation error) was validated
against the 2-D engine (KS p = 1.0, median |Δτ| = 3.6e-5) and run on a
refinement ladder: Nx ∈ {1024, 2048, 4096}, dt ∈ {2e-3, 1e-3, 5e-4, 2.5e-4,
1.25e-4}, plus a ρ-floor sensitivity check (1e-14 → 1e-10: identical to
machine precision).

**Realization of the phenomenon.** The transverse near-field density ends in a
hard edge at τ_max = 2.285 (density drops from ~10⁻¹ to zero in one bin),
while the axial density continues smoothly to 3.224 — the spin term *advances*
the z > L/2 half of the guide and turns the slowest would-be-late trajectories
around (3 permanent non-arrivals), which is exactly the Das–Dürr mechanism for
truncated support. One sampled trajectory (idx 1249, starting on the +3σ
truncation edge at x0 = −2.076) scatters chaotically off the nodal structure
of the truncation-induced diffraction ringing and never converges (its arrival
fluctuates without trend from 1.67 to never across all 7 resolutions); it is
the sole reason G2 is graded PARTIAL rather than PASS, and it is an artifact
of the spec's discontinuous truncation, not of the spin dynamics. G3's
contrast realized in the expected direction but at 2.8 % (< 5 %): with k0 = 4
the near detector is only ~1.5 time units downstream, so convective transport
dominates and the axial late tail is compressed; a slower packet (smaller k0)
or a detector farther into the near field would widen it.

**Connection.** This is the S3/A3 observable of GUM Sec. IX.A (Layer 1), whose
gate logic (VIII.F) the T4 audit repaired; the hard cutoff is the non-POVM
candidate statistic. Filed as F-T6-L5-EXEC.

## Caveats

- "This reproduces the phenomenon; it neither proves nor tests POVM-exclusion
  (the corpus's T4-W5 theorem gap) - it builds the testbed any such argument
  must clear."
- The spec's hard ±3σ truncation makes f discontinuous: kinetic energy is
  grid-dependent (log-divergent in k_max) and the leading edge carries
  diffraction-ringing nodes. Within any fixed grid, energy is conserved to
  machine precision (exact propagation), so G1 is meaningful per-run; but the
  edge nodes produce the single non-convergent trajectory that blocks a clean
  G2 PASS. A C¹ smoothing of the truncation is the recommended follow-up.
- τ_max here is the endpoint of the sampled ensemble (N = 2000), not a proved
  support edge of the continuum density; the unresolved edge set has measure
  ≲ 5e-4 of the ensemble.
- The far-field gate compares only trajectories that arrive within T = 8
  (~70 % of each ensemble); both ensembles are truncated identically, but the
  KS null is over the arrived sub-populations.
- The x-domain is periodic (FFT); trajectories leaving [-14.5, 44.5] are
  frozen and recorded as non-arrivals (1/2000 transverse in the 2-D run).
- G5's ratio < 1 (transverse earlier on average) is a mean-level statement;
  the support-level statement (transverse cutoff earlier than axial last
  arrival, zero transverse events beyond it) is the A3 statistic.
- Both σ_y eigenstates give mirror-image trajectories under z → L − z with
  identical arrival statistics (initial state z-symmetric), so the M_y = +ρ/2
  choice is without loss of generality.
