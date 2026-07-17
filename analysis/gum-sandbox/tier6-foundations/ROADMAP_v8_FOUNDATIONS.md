# ROADMAP v8 — FOUNDATIONS (Tier 6)

Successor to ROADMAP_v7_HARDENING (theory-audit/). Source analysis:
`FOUNDATIONS_ANALYSIS.md` (findings F-T6-1…6). Method unchanged:
replication-first, within-model, independent code; every workstream has
named gates and a RESULTS.md; nothing bears on nature.

Units throughout: ħ = m = 1 unless stated.

## L1 — Contextuality exhibit: the flipped Stern–Gerlach [F-T6-2a]
**Goal.** First in-repo execution of the dBB/GUM guidance-sector
demonstration that "measurements don't measure": same Ψ(0), same X(0),
reversed field gradient ⇒ same spatial motion, flipped outcome label.
**Model.** 1-D two-component spinor; impulsive SG coupling H_int = λzσ_z
(kick ∓λτ_imp), then free split-step evolution; guidance
dZ/dt = Im(Ψ†∂_zΨ)/Ψ†Ψ.
**Gates.** G1: nodal-line theorem — for z-symmetric Ψ(0), zero median
crossings across all trajectories. G2: outcome flip — under λ → −λ with
identical (Ψ(0), Z(0)) ensemble, spatial exit side unchanged for every
trajectory, spin label flipped for every trajectory. G3: Born check —
upper-exit fraction = ½ within 3σ binomial. G4: determinism control —
trajectory reproducibility under halved dt at 10⁻⁶.
**Deliverables.** `L1/l1_sg.py`, `L1/l1_results.json`, `L1/RESULTS.md`,
trajectory figure.

## L2 — Momentum time-of-flight: dBB exact + Nelson leg [F-T6-2b/6]
**Goal.** Execute Bricmont App. 1 exactly, then the GUM-specific
extension: the same ToF statistic under *Nelson* dynamics (GUM's actual
Sec. III kinematics) at finite diffusion.
**Model.** Ψ(x,0) = π^{-1/4}e^{-x²/2} (S ≡ 0; dBB particles at rest).
dBB leg: verify X(t) = X(0)√(1+t²) (analytic vs integrated guidance);
p := X(T)/T at large T distributed as |Ψ̂(p)|² = π^{-1/2}e^{-p²}.
Nelson leg: dX = (v + u)dt + √(2ν)dW, u = ν∂_x ln ρ_ψ, ν = ½;
equivariance (ensemble stays |ψ_t|²) and the same ToF distribution.
**Gates.** G1: dBB trajectory law to 10⁻⁶ (RMS, T = 100). G2: dBB ToF
KS-test vs |Ψ̂|², p > 0.1, N ≥ 20,000. G3: Nelson equivariance — ensemble
vs |ψ_t|² sup-norm of binned densities within 4/√N-class floor at three
times. G4: Nelson ToF KS vs |Ψ̂|², p > 0.1. G5: ν-dial — ToF distribution
invariant across ν ∈ {0.1, 0.5, 1.0} (pairwise KS p > 0.05).
**Deliverables.** `L2/l2_tof.py`, `L2/l2_results.json`, `L2/RESULTS.md`.

## L3 — Lock-2 conservation lemma: the leakage-exponent ladder [F-T6-5]
**Goal.** Toy-grade discharge of the T4-W3 debt: the (c/c_L)⁵ leakage
exponent presupposes vanishing monopole/dipole source channels; prove +
measure that conservation supplies exactly that.
**Model.** Scalar wave equation (1/c_L²)∂_t²p − ∇²p = s(x,t) in 3-D,
compact oscillating sources, far-field radiated power P(c_L) measured;
sources: (a) monopole-breathing (∫s dV oscillates — non-conserved),
(b) dipole (∫s dV = 0, ∫x s dV oscillates — momentum-violating),
(c) quadrupole (both moments static — conserved-source class).
Analytics: multipole radiation P ∝ (c/c_L)^{2ℓ+1} per channel at fixed
source frequency; lemma text: s built from conserved density/current ⇒
ℓ = 0, 1 channels cancel identically.
**Method.** Spherical-harmonic reduction to radial 1-D per ℓ (exact
multipole solve), plus a direct 3-D FFT check at one point.
**Gates.** G1: measured exponent for (a) = 1 ± 0.1 over ≥ 1.5 decades of
c_L. G2: (b) = 3 ± 0.15. G3: (c) = 5 ± 0.25. G4: lemma statement with
proof sketch printed in RESULTS.md; the GUM mapping (s from knot
mass/momentum densities; Lock 2 fallback (c/c_L)³ ⇒ restored (c/c_L)⁵)
stated with its residual assumption named (source built from locally
conserved densities — still an assumption about the constraint sector,
now isolated).
**Deliverables.** `L3/l3_leakage.py`, `L3/l3_results.json`,
`L3/RESULTS.md`, exponent-ladder figure.

## L5 — A3 phenomenon on an independent engine [Sec. 4 of analysis]
**Goal.** First in-repo instantiation of the Layer-1 core phenomenon:
spin-dependent Bohmian arrival times in a waveguide; transverse spin ⇒
hard support cutoff τ_max; axial ⇒ flux-form; far-field convergence.
**Model.** 2-D Pauli guidance: ψ(x,z,t) spinor, waveguide ground state
(hard-wall or truncated) in z × right-moving Gaussian in x; current
j = Im(ψ†∇ψ) + ½∇×(ψ†σψ) (spin term is the discriminator); split-step
FFT evolution; trajectories to a detector plane at x = d (near field
d ~ few widths; far field d ≫).
**Gates.** G1: norm + energy conservation ≤ 10⁻⁶. G2: transverse-spin
near-field arrival density has hard support — zero arrivals beyond
measured τ_max, with τ_max stable under grid/dt refinement to 2%.
G3: axial-spin arrival density strictly positive over the same window
(no cutoff), matching quantum-flux form within binned tolerance.
G4: far-field null — transverse vs axial arrival distributions converge
(KS p > 0.05) at d_far. G5: the X-Vb1-style internal consistency ratio
reported with error.
**Deliverables.** `L5/l5_arrival.py`, `L5/l5_results.json`,
`L5/RESULTS.md`, arrival-histogram figure.
**Scope caveat (printed in advance).** This reproduces the *phenomenon*;
it neither proves nor tests POVM-exclusion (the T4-W5 theorem gap) — it
builds the testbed any such argument must clear.

## L4/L6 — Analysis-only (folded into FOUNDATIONS_ANALYSIS.md)
L4: GUM↔GRWm discriminator ledger (F-T6-3) — done, §2.3.
L6: PLC selection principle applied to GUM's effective-m/VI.D P1
premise — exploration note, §2.4; full transcription deferred (needs
the induced-gravity sector, out of tier envelope).

## Execution
Workstreams L1, L2, L3, L5 are independent → one parallel workflow
fan-out (independent agents, independent code), then a synthesis
adjudication (`TIER6_ADJUDICATION.md`) by the coordinator. Gate
verdicts are the agents'; adjudication and any F-T6 promotions are the
coordinator's. All artifacts commit to branch `claude/analyze-gum-po`.
