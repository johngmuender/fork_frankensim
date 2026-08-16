# U1 — Junction–antijunction (J–J̄) spectrum from printed parameters [F-T10-U1]

**Phase:** Tier 10 (Phase U), ROADMAP_v12_JUNCTION.md workstream U1
(pre-registered; gates frozen before computation).
**Date:** 2026-08-16.  **Code:** `u1_dumbbell.py` (deterministic, no RNG,
self-checking: Airy validation 7.8e-10, T1-limit reproduction exact to
machine, dumbbell grid convergence 2.5e-7; ~6 s wall).
**Outputs:** `u1_results.json`, `u1_fig.png`, `u1_run.log`.
**Sources:** T1 frozen machinery reused verbatim (variable-mass radial
Schrödinger, Zhu–Kroemer central / BenDaniel–Duke systematic, frozen
cutoff f\* = 1.9724 read from `tier9-glueball/T1/t1_results.json` — never
recalibrated); σ = 0.19 GeV² [IM] (corpus Q-5 inversion); junction-mass
anchors: M_J/√σ = 0.1355(36) 2+1D SU(3) lattice (JHEP 12 (2025) 019,
arXiv:2508.00608 [IM]); f_q central 0.17 GeV (corpus inversion window
0.14–0.20); state class per Csörgő–Gyulassy–Kharzeev J. Phys. G 30, L17
(2004) and Frenklakh–Kharzeev–Rossi–Veneziano JHEP 07 (2024) 262 [IM];
`JUNCTION_PAPER_DIGEST.md`, `JUNCTION_CONTEXT_ANALYSIS.md`,
`u_context_agents.json`.

**Epistemic frame (binding).  Within-model; nothing here bears on
nature.**  STAR/lattice/Regge values are [IM] imported anchors.  Under
seal q-θ every mass number below that confronts a hadron mass is
reportable ONLY at V.F grade: *dimensionally secure / structurally
plausible / quantitatively unclaimed*.  **No fitting to STAR observables
anywhere** — the three STAR values appear in no computation here (the
code performs a mechanical token self-scan on its own source and on the
results JSON; result: clean).  The junction interpretation of the STAR
data is itself contested (CGC-saturation rival, neutron-skin and
strangeness confounds — context §4); nothing below leans on it.
**Numerology pre-emption (U1-G4):** the corpus's dimensionless closure
benchmark 𝔠 = 2.37 ± 0.09 shares digits with GeV-denominated numbers by
man-made-unit coincidence only; no digit identity is cited anywhere in
U1 as structure, and may never be.  Junction energy is **[CJ-new]**: the
corpus prints no junction energy/mass/inertia; m_J is a scanned explicit
parameter, never fitted.

---

## 1. Question and answer

**Question (pre-registered).**  The minimal knot-free junction state —
two order-3 junctions joined by three tubes (the "baryonium glueball"
class) — at the corpus tension σ = 0.19 GeV² [IM]: where does its ground
state land dimensionlessly, and how does the answer depend on the
(unprinted) junction energy?

**Answer.**  Under the spec-frozen extension of the validated T1
machinery, the J–J̄ ground state lands at

- m/√σ = **2.59** (envelope 2.45–2.88) at m_J = 0,
- m/√σ = **2.70** (2.55–2.86) at m_J/√σ = 0.1355 (2+1D lattice [IM]),
- m/√σ = **2.94** (2.67–3.12) at m_J/√σ = 0.39 (f_q/√σ),

i.e. **BELOW the loop sector's 0⁺⁺ head (3.583) at every pre-registered
junction mass** — the junction sector, so computed, opens a lighter floor
under the Tier-9 closed-loop spectrum, and its excitations interleave
with the loop ladder rather than duplicating it (sector-relation verdict,
§5).  The dependence on the junction energy is weak and monotone: the
full m_J scan moves the ground by only +0.35 in m/√σ (+14 %), because
2m̂_J enters statics and inertia with partially cancelling effects.

**The honest counterweight (computed, printed):** the transfer of T1's
inertia *pattern* to the dumbbell is a genuine ambiguity.  Under the
physically derived stretch inertia (σL/4 + m_J/2 instead of the
spec-frozen 3σL + 2m_J) the ground rises to m/√σ ≈ 5.47–5.77 — above the
loop 0⁺⁺ head.  The **downward-extension half of the verdict is owned by
the frozen inertia convention; the interleaving half survives both
conventions.**  Both statements are within-model, V.F, quantitatively
unclaimed.

## 2. Method and provenance (frozen; deviations printed)

Dumbbell configuration coordinate: junction separation L; dimensionless
x = √σ·L, ε = E/√σ, m̂ = m_J/√σ.

```
V_N(x) = 3x + 2m̂ + (πN − π/4)·(1 − exp[−f*x])/x,      μ(x) = 3x + 2m̂
```

- **3x** — three tubes at leading order parallel along the J–J̄ axis
  (energy exactly 3σL); Y-law end geometry idealized (A1 below).
- **πN/x** — transverse phonons of an open tube with both ends fixed on
  junctions (Dirichlet), ω_n = nπ/L; N = Σn over 3 tubes × 2 transverse
  polarizations (T1 census conventions transferred).
- **−π/4·(1/x)** — ζ-regularized zero point: 6 polarization channels ×
  (π/2L)ζ(−1) = −π/(24L) each = three times the single open-string D = 4
  Lüscher term −π/(12L).
- **(1 − e^{−f\*x})** — T1's short-distance smoothing transferred as
  junction-end smoothing at the FROZEN f\* = 1.9724361517 (read from
  T1's results file; not recalibrated).
- **μ(x) = 3x + 2m̂** — configuration inertia = full configuration mass,
  imitating T1's μ(ρ) = 2πσρ per the frozen spec (A2 below).
- Radial Schrödinger equation in x, Zhu–Kroemer ordering central,
  BenDaniel–Duke as the ordering systematic; Dirichlet at 0 and x_max;
  scheme envelope = T1's frozen set {ZK × f ∈ [f\*/1.5, f\*, 1.5f\*]} ∪
  {BDD × f ∈ [8, 32]}.
- m̂ ∈ {0, 0.1355, 0.39}: pre-registered, never fitted.  **Zero fitted
  parameters in this workstream** (f\* was frozen by T1's validation
  before Tier 10 existed; no Tier-10 anchor entered any calibration).

**Coordinator-owned amendments (printed by the code, verbatim in
`u1_results.json`):**

- **[U1-A1] Y-law geometry idealization.**  At leading order in the
  adiabatic separation coordinate the three tubes run parallel along the
  J–J̄ axis (energy exactly 3σL); the Y-law 120° opening (a banked corpus
  structural consistency, lattice-favored [IM]) lives in the transverse
  junction structure collapsed to a point here.  Transverse junction
  shape, tube–tube interaction, and junction-core excitation are beyond
  this leading order.
- **[U1-A2] Inertia convention.**  μ = 3σL + 2m_J (T1-pattern, full
  configuration mass) is the spec-frozen central choice.  The physically
  derived symmetric-stretch inertia (ends at ±L/2, linear velocity
  profile) is σL/4 + m_J/2; for the ring breathing mode the two
  prescriptions coincide, for the dumbbell they do not.  The variant is
  computed and printed as a sensitivity (ground ε: 2.59 → 5.77 at m̂ = 0)
  **outside** the frozen envelope — the single largest systematic in this
  workstream, priced, not hidden.
- **[U1-A3] Phonon census transfer.**  Open-tube Dirichlet modes,
  3 × 2 polarization channels, ζ-regularized zero point −π/(4L), T1
  smoothing factor with frozen f\* as junction-end smoothing; scheme
  envelope transferred unchanged.  No new parameter; nothing
  recalibrated.
- **[U1-A4] Junction-mass status.**  m_J enters statics (+2m_J) and
  inertia; 0.1355 is a **2+1D** SU(3) lattice value — the 3+1D
  counterpart is unmeasured and dimensional transfer is an untested
  assumption, printed; 0.39 = f_q/√σ is a corpus scale *analogy*, not a
  junction mass.  Junction energy remains [CJ-new].

## 3. Validation (U1-G1) — before any junction number entered

| check | requirement | measured | verdict |
|---|---|---|---|
| Airy: two-body linear potential (μ = 1/2, k = 1 ⇒ E_n = \|a_n\|), same solver code path, 6 levels, h²-Richardson | ≤ 1e-6 rel. | max 7.8e-10 | **PASS** |
| ZK/BDD agreement at constant mass (orderings must coincide) | — | 0.0 | PASS |
| T1 limit: generic solver fed T1's loop V, μ at frozen f\*, T1 grid — 0⁺⁺ | ≤ 1e-3 rel. of 3.5826743623 | 0.0 (machine-identical) | **PASS** |
| T1 limit — 0⁻⁺ (M = 4) | ≤ 1e-3 rel. of 8.2749381892 | 0.0 | PASS |
| dumbbell grid convergence (npts, x_max doubling) | < 5e-6 | 2.5e-7 | PASS |

## 4. Spectrum at the three pre-registered junction masses (U1-G2)

σ = 0.19 GeV² [IM]; √σ = 0.43589 GeV.  ε = m/√σ; bands = frozen scheme
envelope (ordering × cutoff), model systematics only.  GeV values V.F.
Lowest five configuration levels per m̂ (adiabatic model assignment;
ground = scalar-class head; the N = 1 sextet decomposes under S₃ tube
permutation × polarization J_z = ±1 × end swap into singlet + doublet
structure per polarization — exact degeneracy is an adiabatic artifact,
same caveat class as T1's M = 2 level):

**m̂ = m_J/√σ = 0:**

| level | ε | envelope | m [GeV, V.F] |
|---|---|---|---|
| ground (N=0) | **2.588** | 2.453–2.884 | 1.128 |
| separation excitation (n_r=1) | 4.728 | 4.665–4.925 | 2.061 |
| first phonon level (N=1, ×6 channels) | 5.615 | 5.243–6.234 | 2.447 |
| phonon + separation (N=1, n_r=1) | 7.037 | 6.818–7.491 | 3.067 |
| second phonon level (N=2) | 8.223 | 7.594–8.680 | 3.584 |

**m̂ = 0.1355 (2+1D lattice anchor [IM]; dimensional-transfer caveat):**

| level | ε | envelope | m [GeV, V.F] |
|---|---|---|---|
| ground | **2.695** | 2.546–2.857 | 1.175 |
| separation excitation | 4.807 | 4.735–4.898 | 2.096 |
| first phonon level | 5.816 | 5.421–6.436 | 2.535 |
| phonon + separation | 7.195 | 6.959–7.653 | 3.136 |
| second phonon level | 8.473 | 7.820–8.927 | 3.693 |

**m̂ = 0.39 (= f_q central 0.17 GeV / √σ; scale analogy only):**

| level | ε | envelope | m [GeV, V.F] |
|---|---|---|---|
| ground | **2.943** | 2.674–3.122 | 1.283 |
| separation excitation | 5.001 | 4.874–5.104 | 2.180 |
| first phonon level | 6.216 | 5.782–6.852 | 2.710 |
| phonon + separation | 7.520 | 7.253–8.007 | 3.278 |
| second phonon level | 8.950 | 8.257–9.401 | 3.901 |

**Inertia-convention sensitivity (printed, outside the envelope):**
stretch-inertia ground ε = 5.77 / 5.60 / 5.47 at m̂ = 0 / 0.1355 / 0.39 —
a factor ~2.2 on the ground state.  The junction-mass dependence under
the stretch variant is *inverted* (heavier junction → lower ε: the
kinetic zero-point reduction outweighs +2m̂), another marker that the
dumbbell's absolute placement is convention-owned.

## 5. Confrontation (U1-G3; V.F grade locked)

| quantity | m̂ = 0 | m̂ = 0.1355 | m̂ = 0.39 | anchors / T1 |
|---|---|---|---|---|
| ground ε (envelope) | 2.59 (2.45–2.88) | 2.70 (2.55–2.86) | 2.94 (2.67–3.12) | loop 0⁺⁺ 3.58; lattice 0⁺⁺ 3.5 [IM] |
| ground below loop 0⁺⁺? | yes | yes | yes | — |
| ground vs lattice 0⁺⁺ 3.5 | −26 % | −23 % | −16 % | not claimed as any state |
| levels touching 5.41–5.9 (X(2370) [IM] ∪ lattice 0⁻⁺ [IM]) | N=1 phonon (5.61, band 5.24–6.23) | N=1 phonon (5.82, band 5.42–6.44) | N=1 phonon band edge (5.78) | X 5.41–5.45; lat 0⁻⁺ 5.5–5.9 |
| levels interleaving T1 loop gaps | 4 | 3 | 3 | T1 ladder 3.58 … 8.27 |
| stretch-inertia ground (sensitivity) | 5.77 | 5.60 | 5.47 | above loop 0⁺⁺; convention-owned |

**Sector-relation verdict (computed, printed by the code):**

1. **EXTENDS DOWNWARD** — at every pre-registered junction mass the
   J–J̄ ground state lies below the loop sector's 0⁺⁺ head: the junction
   sector opens a new, lighter floor under the closed-loop spectrum.
2. **INTERLEAVES** — 10 dumbbell level placements across the m̂ set fall
   inside gaps of the T1 loop ladder; the two sectors would overlay as
   one interleaved spectrum, not two separated bands, and no dumbbell
   level duplicates a loop level to within 5 %.
3. **CONVENTION-FRAGILE** — under the physically derived stretch inertia
   the ground rises to ε ≈ 5.5–5.8, above the loop 0⁺⁺ head; the
   downward-extension half of the verdict is owned by the spec-frozen
   T1-pattern inertia, not by the dumbbell physics itself.  The
   interleaving half survives both conventions.  **Either answer was
   pre-registered as a finding; the finding is that the sector relation
   is convention-split, and the split is now priced.**

No dumbbell state is identified with X(2370), any lattice glueball, or
any STAR-discussed exotic; the N = 1 phonon level's overlap with the
5.41–5.9 window is reported as geometry of the table, not as a match
(its J^PC content is an adiabatic sextet, not an established 0⁻⁺).

**Tetraquark-kill note (n = 3 only).**  Every network here is built from
order-3 junctions — the corpus's Q-1 inventory class with the tetraquark
kill armed and junction order owned by WS1-H3 (an ansatz, [CAL]-at-best;
edition-status seam F-K0-1 carried).  No 4-fold junction and no
baryonium-*tetraquark* (knot-terminated) state is constructed or priced;
the J–J̄ dumbbell is the knot-FREE minimal network.  Nothing here bears
on the tetraquark kill's standing.

## 6. Gates (pre-registered; may not move)

| gate | requirement | measured | verdict |
|---|---|---|---|
| U1-G1 | Airy ≤ 1e-6; T1-limit reproduction | Airy 7.8e-10; T1 0⁺⁺ and 0⁻⁺ reproduced to 0.0 (≤ 1e-3 gate) | **PASS** |
| U1-G2 | spectrum at all three m_J values, lowest 3+ states, zero fitted parameters, uncertainty = scheme envelope | 5 levels × 3 m̂ values; envelope = frozen T1 scheme set; f\* frozen, m̂ pre-registered, nothing fitted; inertia sensitivity additionally printed | **PASS** |
| U1-G3 | dimensionless confrontation table + sector-relation verdict, all V.F | table §5; verdict: extends-downward + interleaves, with the computed convention-fragility printed as part of the finding | **PASS (finding = convention-split sector relation)** |
| U1-G4 | register: seal q-θ; numerology pre-emption; no grade motion; [CJ-new] status carried | header + §2 A4; mechanical STAR-token self-scan clean; no stake, no forbidden sentences, no grade motion | **PASS** |

## 7. Key numbers

| quantity | value |
|---|---|
| σ [IM] / √σ | 0.19 GeV² / 0.43589 GeV |
| f\* (frozen by T1) | 1.9724361517 |
| zero point (3 tubes × 2 pol, ζ-reg) | −π/4 per 1/L (= 3 × open-string Lüscher, D = 4) |
| ground ε at m̂ = 0 / 0.1355 / 0.39 | 2.588 / 2.695 / 2.943 |
| ground envelopes | 2.45–2.88 / 2.55–2.86 / 2.67–3.12 |
| ground m [GeV, V.F] | 1.128 / 1.175 / 1.283 |
| T1 loop references | IP 0⁺⁺ 3.583, 0⁻⁺ 8.275; NG N=1 4.800 |
| m_J scan lever arm on ground | +0.355 in ε over full m̂ range (+13.7 %) |
| stretch-inertia sensitivity (ground) | 5.77 / 5.60 / 5.47 (convention-owned; outside envelope) |
| interleaving placements (Σ over m̂) | 10; duplications within 5 %: 0 |
| validation | Airy 7.8e-10; T1-limit 0.0; grid 2.5e-7 |

## 8. Honest caveats (U1-G4 register)

1. **The inertia convention is the dominant systematic and it flips the
   headline.**  The spec-frozen μ = 3σL + 2m_J (T1 pattern) puts the
   ground at ε ≈ 2.6–2.9; the physically derived stretch inertia puts it
   at ε ≈ 5.5–5.8.  T1's pattern was physically exact for the ring
   breathing mode; for the dumbbell it is an analogy.  Both values are
   printed; neither is claimed.  This is exactly why nothing here rises
   above V.F.
2. **The junction-mass anchor is dimensionally transplanted.**  0.1355
   is a 2+1D SU(3) result; no 3+1D lattice value exists.  The corpus
   itself prints NO junction energy ([CJ-new]) — this entire parameter
   direction is imported scaffolding, and the weak m̂ dependence (+14 %
   over the scan) is the useful output: the answer is *not* hostage to
   the unprinted junction energy within the scanned window under the
   frozen convention (it IS hostage to the inertia convention).
3. **Adiabatic degeneracies are artifacts.**  The N = 1 level is a
   6-channel sextet at one energy; real dynamics would split it, and its
   J^PC content is model-assigned, not derived.  Its overlap with the
   X(2370)/lattice-0⁻⁺ window is table geometry, not a candidate claim.
4. **Y-law geometry enters only as an idealization** (parallel tubes at
   leading order); transverse junction structure, tube–tube interaction,
   junction-core excitations, and any knot-terminated (quark-analog)
   admixture are all outside the model space.
5. **The state class itself is NOT-CONSTRUCTED in the corpus.**  Course
   20.1's "exactly" inventory sentence excludes knot-free
   junction-network composites as printed; this workstream computes an
   extension the corpus never wrote down — offer-class input to U3's
   annotations, not corpus content.
6. **No absolute-mass claim above V.F.**  All GeV values are
   √σ-rescalings with σ [IM] (inversion-provenanced, κ_q never valued).
   The ground-state GeV range 1.13–1.28 is not identified with any
   hadron, glueball candidate, or STAR-discussed state.
7. **Forbidden-sentence compliance:** this workstream does not say and
   does not imply "GUM predicted the junction result" or "STAR
   confirms/refutes GUM"; no stake is issued; no grade moves; STAR's
   1.84/0.64/1.04 appear only in this compliance sentence and the
   code's self-scan disclaimer, never as computation inputs; rivals to
   the junction interpretation (CGC saturation, neutron skin,
   strangeness) are carried at equal weight in the Tier-10 context and
   nothing here adjudicates among them.
