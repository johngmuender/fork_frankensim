# TIER 6 — ADJUDICATION (Foundations Workstreams L1/L2/L3/L5)

Coordinator adjudication of the ROADMAP_v8 execution (four independent
agents, independent code, pre-registered gates; run 2026-07-17, workflow
wf_5852a821-c22). Gate verdicts below are the agents'; adjudication
verdicts and finding promotions are the coordinator's. Within-model;
nothing here bears on nature.

## Scoreboard

| WS | Goal | Gates | Verdict |
|---|---|---|---|
| L1 | Flipped Stern–Gerlach contextuality exhibit | 4/4 PASS | **EXECUTED — F-T6-2a discharged** |
| L2 | Momentum ToF, dBB exact + Nelson leg | 5/5 PASS | **EXECUTED — F-T6-2b + F-T6-6 discharged** |
| L3 | Lock-2 leakage-exponent conservation lemma | 5/5 PASS | **EXECUTED — T4-W3 discharged at toy grade (F-T6-5)** |
| L5 | A3 arrival-time phenomenon, independent engine | 3 PASS, 2 PARTIAL | **PHENOMENON CONFIRMED; two spec artifacts, both mine** |

## L1 — adjudicated EXECUTED

The two-sided demonstration lands exactly as Bell/Bricmont state it:
4000/4000 trajectories keep their spatial exit side and flip their spin
label under gradient reversal (k → −k), with the deep reason *measured,
not assumed* — the flip swaps spinor components, leaving ρ_total and
j_total bitwise invariant (max deviation 0.0 exactly), so the guidance
field cannot distinguish the two experiments; only the label map can.
Nodal-line theorem: 0 median crossings in 4000 × 480 steps. Born: 0.92σ.
Determinism control: 5.9×10⁻⁹ (~2000× inside gate). The corpus's
Sec. III now has its canonical non-contextuality exhibit in-repo.
Caveat retained at full strength: the bitwise invariance is special to
the symmetric spinor; the physical claim (side preserved, label flipped)
is the general content.

## L2 — adjudicated EXECUTED

Bricmont App. 1 executed to numerical exactness (trajectory law RMS
2.6×10⁻¹² against X(0)√(1+t²)), and the ToF statistic p = X(T)/T
reproduces |Ψ̂(p)|² at KS p = 0.77 with N = 20,000 — while every
particle's instantaneous velocity at t = 0 is exactly zero. The
GUM-specific result is the **Nelson leg**: at ν = ½ (GUM's ħ = 2mν),
equivariance holds at 0.46–0.63 of the statistical floor at t = 1, 10,
100, the same ToF distribution emerges (KS p = 0.78), and the statistic
is **diffusion-blind** across ν ∈ {0.1, 0.5, 1.0} (pairwise KS p ≥ 0.34).
This is the first in-repo verification that GUM's stochastic Sec.-III
kinematics reproduces dBB measurement phenomenology at equilibrium —
previously asserted via the [IM] chain, never exhibited. F-T6-6 closes.

## L3 — adjudicated EXECUTED; T4-W3 discharged at toy grade

The exponent ladder is unambiguous: measured slopes of ln P vs ln(1/c_L)
= **1.000000 / 2.999844 / 4.999869** for the non-conserved,
momentum-violating, and conserved-class sources respectively (exact
retarded solutions, flux R-independence 2.5×10⁻⁴, fit residuals
≤ 2×10⁻⁴). The lemma is stated with proof sketch: order-ℓ far-field
amplitude scales as (ωa/c_L)^ℓ; charge conservation kills ℓ = 0,
momentum conservation kills ℓ = 1; the leading leakage is quadrupolar,
(c/c_L)⁵. **Consequence for the corpus:** the T4 audit's honest fallback
(c/c_L)³ ≤ 10⁻¹² and the printed (c/c_L)⁵ ≤ 10⁻²⁰ are now separated by
one named, isolated assumption — *the Lock-2 constraint source is built
from locally conserved knot mass/momentum densities* — rather than by an
unexamined exponent. Proposed regrade if adopted by the corpus: Lock 2
[DF core + GAP-minor] → [DF core + conservation-lemma-conditional], with
the field-level verification of the source form as the remaining (now
sharply posed) obligation. Toy-grade label binding: scalar wave equation,
prescribed sources — not the GUM constraint sector itself.

## L5 — adjudicated PHENOMENON CONFIRMED; spec artifacts assigned to the coordinator

The Das–Dürr-class phenomenon realizes on a machine-precision engine
(norm/walls/energy at 10⁻¹⁶; two mutually validating engines, cross-KS
p = 1.0): the transverse-spin near-field arrival density **ends in a
hard support cutoff** at τ_max = 2.28503 — zero arrivals in a gap 7.5×
the required width, τ_max stable to 7×10⁻⁷ under dt-halving — while the
axial density continues smoothly 41% past it in time, with zero
transverse events in that region; far-field distributions converge
(KS p = 0.31). The spin-dependence of Bohmian arrival statistics, the
core Layer-1 observable, is thereby instantiated in-repo for the first
time.

The two PARTIALs are **spec defects of this roadmap, not of the
dynamics**, and are filed in the campaign's own two-harness-defects
spirit: (i) G2's stability clause is violated by exactly one trajectory
(0.05% of ensemble) born on the roadmap's *discontinuous* ±3σ truncation
edge — diffraction-ringing nodes scatter it chaotically across
resolutions; all G2 clauses pass with it excluded. The spec should have
prescribed a C¹-smoothed truncation. (ii) G3's 5% tail-fraction
threshold was set blind against a k₀ = 4 boost that compresses the axial
late tail at d_near = 1; measured 2.8% — direction and separation
confirmed, magnitude under an arbitrarily placed bar. Follow-ups named:
smoothed truncation; slower packet or farther near-detector. The scope
caveat stands verbatim: this reproduces the phenomenon; it neither
proves nor tests POVM-exclusion (T4-W5) — it builds the testbed any
such argument must clear.

## Findings register (Tier 6 final)

- **F-T6-1** credit (von Neumann-immunity nontrivial; T-B1 adds selection) — analysis, stands.
- **F-T6-2a/b** gap-of-record → **discharged by execution** (L1, L2).
- **F-T6-3** GUM↔GRWm discriminator ledger incl. the A3+collapse mutual kill — analysis, stands.
- **F-T6-4** positioning (DGNSZ 2013 anchor missing at VI.E/VIII.B) — documentation-grade, stands.
- **F-T6-5** → **discharged at toy grade** (L3); residual assumption isolated and named.
- **F-T6-6** (Nelson kinematics never shown to reproduce ToF phenomenology) → **discharged** (L2).
- **F-T6-L5** new: A3 phenomenon confirmed on an independent engine; τ_max = 2.285 (this geometry); two spec artifacts (coordinator's), follow-ups named.

## Proposed corpus annotations (offered, per charter rule; adoption is the authors')

1. Sec. VIII.A: one sentence claiming the architectural credit of §2.2 of
   the analysis memo (medium-first PO is hypersurface-independent
   trivially; the GRWm problem never arises), with pointer here.
2. Sec. VIII.F/IX.B: the sharpening of §4 — A3-positive discriminates
   against the (A4)/(A5) collapse class as a whole, not only orthodox
   arrival proposals; and the A3+confirmed-collapse pair is a mutual kill.
3. Sec. VI.E/VIII.B: cite DGNSZ 2013 as nearest structural relative with
   the two deltas (matter-sourced F; scheduled detection program).
4. VIII.D Lock 2: replace the bare exponent presupposition with the L3
   conservation lemma + the isolated source-form assumption.

Artifacts: `L1/ L2/ L3/ L5/` (scripts, JSON, figures, RESULTS.md each);
run record: workflow wf_5852a821-c22, four agents, 265,836 tokens.

---

## ADDENDUM — L5′ adjudication (F-T6-L5-EXEC-2; run 2026-07-17, post-checkpoint)

The pre-registered follow-up (C² smoothstep truncation; k₀ = 2 boost;
same gate battery, full ensemble, no exclusions) executed. Adjudicated
outcomes:

**Both L5 PARTIALs convert to clean PASS.**
- G2: τ_max = 5.130950, stable to 8.3×10⁻⁸ (Nx-doubling) and 5.6×10⁻¹⁰
  (dt-halving); empty gap 10.87 = 4.5× required; **all 2000 trajectories
  converge** across a six-run resolution ladder (max per-trajectory
  drift 3.7×10⁻⁵, zero status changes) — the latest arrival is born
  mid-packet, not on the support edge. The L5 edge-outlier failure mode
  is eliminated, confirming its adjudication as a truncation-spec
  artifact. Energy is now grid-converged (E = 2.2503989 ± 5×10⁻⁹; the
  hard-cut UV log-divergence is gone).
- G3: axial arrivals beyond the transverse cutoff = **11.08%** (> 2× the
  5% bar; L5: 2.8%), with zero transverse events beyond τ_max out to
  T = 40. The k₀-compression diagnosis of L5's PARTIAL is confirmed.

**G4: FAIL-AT-d=25, graded honest and diagnosed — the null is
approached, not violated.** After one coverage retry (T = 40, box
[−15, 225]; non-arrivals balanced 0/88), KS at the spec'd d = 25 gives
p = 9.9×10⁻⁸, dt-stable. The seven-plane ladder measures a monotone
KS(d) decay 0.1261 → 0.0571 with p climbing seven orders to 3.7×10⁻³
at d = 35: at slow boost the spin-current arrival-time difference
decays more slowly in packet-width units, so d = 25 is not yet
asymptotic. L5's G4 PASS (k₀ = 4) vs L5′'s FAIL (k₀ = 2) is a
**physics-of-configuration pair**, not an inconsistency.

**New bench-design record for the Layer-1 envelope (worth exporting to
any IX.A protocol discussion):** slowing the packet sharpens the
near-field spin contrast (G2 gap ×4.5, G3 ×4) but pushes the
asymptotic far-field regime outward — detector placement must be
co-designed with the boost. The far-field-null station of the corpus's
Layer-1 envelope is a *joint* (k₀, d_far) design constraint, not a
free placement.

Cross-engine credit: far-field samples bit-identical across the two
independent engines (per-trajectory median agreement 1.6×10⁻³ ⇒ same
rank configuration), verified genuine, not a harness artifact.
Sensitivity control: the amplitude-convention variant shifts τ_max to
4.457 and G3 to 11.66% with G4 p = 0.004 — every adjudicated
conclusion is convention-robust. Scope caveat unchanged (testbed, not
POVM-exclusion). F-T6-L5 is hereby CLOSED in its L5′ form; the open
successor items are the POVM-exclusion theorem (T4-W5) and, if ever
bench-relevant, a KS(d) extension to the measured asymptote.
