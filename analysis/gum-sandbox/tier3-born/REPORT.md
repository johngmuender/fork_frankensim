# Tier 3 — Born-Rule Relaxation Pilot: REPORT & ADJUDICATION

**Status: COMPLETE.** Valentini–Westman-class de Broglie–Bohm relaxation in a
2-D box; the numerical claim under test is the corpus's ⟨V-SIM-4⟩ summary
(GUM-Ω §III.D): *"near-exponential H̄ decay on system dynamical times for
M ≳ 4; power-law tails and residual pockets for small-M near-integrable
flows."* Run: `born.py` (seeded: phases 42, sampling 12345), N = 20,000
trajectories, RK4 dt = 2×10⁻³, t_final = 4π (the exact ψ-revival period),
runtime 829 s. Deliverables: `hbar_decay.png`, `snapshots_M16.png`,
`hdata.csv`, `results.json`, `run.log`.

## Setup
ψ(t) = M⁻¹ᐟ² Σ e^{iθ_mn} φ_mn e^{−iE_mn t} on [0,π]² (exact evolution — no
PDE error); trajectories integrate v = Im(∇ψ/ψ) analytically from the mode
sum; initial ensemble ρ₀ = |φ₁₁|² (non-equilibrium); coarse-grained
H̄(t) = Σ P̄ ln(P̄/Q̄) at 32×32 (and 16×16 for the grain-dependence caveat).
Noise floors (occupied-cells/2N): 0.0256 (32²), 0.0064 (16²).

## Results (32×32 coarse-graining)

| Run | H̄(0) | H̄(4π) | fitted τ | fit r² | verdict |
|-----|-------|--------|----------|--------|---------|
| M = 4 | 0.443 | 0.377 | 83.0 | **0.43** | weak, non-exponential; large residual |
| M = 9 | 0.796 | 0.217 | 10.2 | 0.91 | near-exponential |
| M = 16 | 0.557 | **0.046** | 4.85 | **0.99** | clean exponential, → noise floor |
| control (ρ₀ = \|ψ(0)\|²) | 0.024 | 0.025 | — (max 0.030) | — | **flat at the noise floor** |

16×16 graining: same ordering (τ = 53.4 / 7.8 / 4.0; final H̄ 0.303 / 0.148
/ 0.018), confirming grain-robustness of the qualitative picture.

## Adjudication

**Replicated:**
1. **Relaxation to Born equilibrium is real and monotone-in-M**: τ falls
   83 → 10 → 4.9 as M grows 4 → 9 → 16, with M = 16 reaching the finite-N
   noise floor within one revival period — the corpus's qualitative engine
   ("faster/cleaner with more modes") confirmed.
2. **Equivariance (the control)**: an ensemble started AT ρ = |ψ|² stays at
   the noise floor for the entire run (max 0.030 vs floor 0.026) — the
   H-theorem's fixed point verified, and with it the padlock premise that
   equilibrium ensembles stay equilibrium.
3. **The particle distribution does not revive** at t = 4π even though ψ
   does (exactly) — relaxation is a property of the trajectories, not the
   wave: visible in `hbar_decay.png` at the marked revival time.

**Mild tension, adjudicated:**
The corpus says "near-exponential for **M ≳ 4**"; this pilot's M = 4 case is
*not* near-exponential (r² = 0.43, residual 0.38). Ruling: **attributable to
the mode set, not a corpus defect** — the first-2×2 set has energies
E ∈ {1, 2.5, 2.5, 4}: one degenerate pair and *commensurate* gaps (1.5, 1.5),
i.e. very few distinct beat frequencies — precisely the corpus's own
"small-M near-integrable" caveat case, which it predicts shows "power-law
tails and residual pockets." Our M = 4 result IS that predicted behavior;
the boundary between "pocketed" and "near-exponential" simply falls between
this particular M = 4 set and M = 9. A replication at M = 4 with a
non-degenerate mode choice (e.g. {(1,2),(2,3),(3,1),(1,4)}) would test
whether "M ≳ 4" holds for generic sets — recorded as follow-up, not a kill.

**Caveats (standing):** single phase seed; fit window t ≤ 2π; coarse-grain
dependence quantified above (τ shifts ~30% between grains); RK4 near nodal
lines guarded by substepping but not error-controlled; N = 2×10⁴ noise floor
limits the deepest observable H̄.

## Epistemic notice
This replicates the *dynamical-relaxation* ingredient GUM imports from the
pilot-wave literature (Valentini's H-theorem + mixing rates) — a real
property of de Broglie–Bohm dynamics, and the corpus's [IM+CAL] grading of
it is fair. It is not evidence for GUM's substrate ontology; it verifies
that the borrowed engine behaves as the corpus says it does.
