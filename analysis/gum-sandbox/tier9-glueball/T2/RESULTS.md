# T2 — Worldsheet-mode census on the GUM tube: the chirality hook (F-T9-T2)

**Campaign:** Tier 9 (Phase T), ROADMAP_v11_GLUEBALL.md workstream T2 (frozen spec;
gates pre-registered). **Date:** 2026-08-15. **Code:** `t2_census.py` (python3 +
sympy, deterministic, runtime 2.6 s, 16/16 internal checks PASS). **Raw output:**
`t2_results.json`. **Figure:** none — the workstream is symbolic/structural; no
computed spectrum exists to plot (declared, not skipped silently). **Sources:**
corpus3/01-GUM-Omega-Paper-v4.3-ext.md (Sec. II.A:114, II.C:132–137, II.E:147,
II.F:151, II.H:159–167, VII.E:563, VII.I:599, VII.J:615–617);
theory-audit/h21_RESULTS.md §2, h25_RESULTS.md §3; tier7-program/M3/RESULTS.md §3,
N2/RESULTS.md §3; GLUEBALL_CONTEXT_ANALYSIS.md; t_context_agents.json
(selection-rules + color-sector archaeology agents).

**Epistemic frame (binding).** *Within-model; nothing here bears on nature.* All
lattice/BESIII statements are [IM] anchors. **Seal q-θ:** this workstream produces
**no mass numbers at all** — every internal-mode scale is symbol-only (M_gap is
never valued in print; context analysis §1) — so the V.F protocol is satisfied
trivially and no quantitative claim exists to grade. **Numerology pre-emption
(mandatory):** the corpus's dimensionless closure number 𝔠 = 2.37 ± 0.09 and the
imported m_X ≈ 2.37 GeV share digits by man-made-units coincidence; that
coincidence is **not structure, is not used anywhere below, and may not be cited
as structure** in any reading of this memo. Grade language never rises; all
corpus-side changes are offers.

## 0. Verdict in one line

**The GUM T2 mismatch tube, as printed, is a PLAIN BOSONIC STRING: the only
worldsheet degrees of freedom derivable from printed structure are the two
transverse translation Goldstones (validated massless, checks B1/C1–C3) — no
pseudoscalar worldsheet mode emerges from print, so the tube inherits the known
3+1D Isgur–Paton/Nambu–Goto J^PC problems (the 0⁻ sector that lattice torelon
spectroscopy resolves only with a massive worldsheet axion [IM: ADLT]). The
chirality hook is real but undischarged: W_χ is the corpus's UNIQUE parity-odd
energy class (check A4), any core-localized axial rotor/twist mode is FORCED to
be a worldsheet pseudoscalar by the printed tensor characters (check A5), and the
axion-type vertex ϑ·ε_ab∂X_a∂X_b is symmetry-allowed with W_χ as its only
printed source (check D2) — but the mode's existence, localization, and gap are
all unprinted. The "tube-core axion" is filed as a [CJ-new]/[DW, testbed-grade]
proposal with mass ∝ M_gap, symbol-only.**

## 1. Gates

| id | requirement (pre-registered) | measured | verdict |
|---|---|---|---|
| T2-G1 | tube background reproduced from campaign machinery with citations; mode-enumeration table complete with printed-evidence column (h25 branch cross-check) | h21 straight-static-tube background (radius a, uniform core; h21_RESULTS.md §2, reused verbatim by M3 §3/N2 §3) specialized to a single tube (amendment A-T2-1, §2); 8-row census (§3), every row carrying file:line evidence or a [CJ-new] flag; h25 cross-check performed: the audited branch content (B1, B2±, B3, B4, knot band-edges) contains **no locking-stratum mediator branch at all** | PASS |
| T2-G2 | quadratic worldsheet action derived symbolically for ≥ the 2 transverse modes (validation: MUST come out massless — Goldstone check) and every candidate internal mode; [CJ-new] flags printed at every unprinted step | transverse modes massless by an exact operator identity for the full printed energy class (check B1) AND by explicit halo reduction to Nambu–Goldstone form with tension σ_hal = πν²f_q²ln(R/a) (checks C1–C3); internal twist-mode worldsheet action extracted mechanically (checks E1–E5) with its uniform-profile limit reproducing the h25-audited B4 dispersion Jω² = (α+β)k² + 4μ_c + m_V² **exactly**; internal masses symbol-only; ansatz steps flagged [CJ-new] inline (§4) | PASS |
| T2-G3 | census verdict with parity classification; if NO pseudoscalar from printed structure, that is the finding; if one exists, structure printed at [CJ-new]/[DW, testbed-grade] | verdict delivered (§5): **no pseudoscalar mode emerges from printed structure — plain bosonic string as printed** (the finding); the forced-parity tube-core-axion candidate is printed at [CJ-new]/[DW, testbed-grade] proposal class, mass ∝ M_gap symbolic | PASS |
| T2-G4 | within-model; SILENT-corpus verdicts carried from archaeology; no grade motion | register held: no numbers, no grade motion, no stake, forbidden sentences absent; the archaeology's SILENT verdicts (tube-core chirality; cross-stratum chirality transfer; compact-loop quantization) carried verbatim (§3, §6) | PASS |

## 2. Background: what is reused, and the one printed amendment

- **Geometry (mandate honored):** the campaign's straight-static-tube
  parametrization — a straight tube along ẑ of radius a with uniform core —
  is h21's (`h21_RESULTS.md` §2: "static network of straight tubes, radius a,
  uniform … inside, 0 outside"; `h21_mc.py`), reused verbatim by M3
  (`M3/RESULTS.md` §3, gate condition) and N2 (`N2/RESULTS.md` §3). T2 needs the
  single-tube member of that family, not the Poisson web. **Amendment A-T2-1
  (coordinator-owned, printed):** specialization to one straight tube, plus an
  orientation *halo* θ₀ = ν·(polar angle), ν = 1/3, outside the core — the halo
  is forced by Q-1's fractional frame winding 2π/n, n = 3 (corpus3/01:615), but
  its profile normalization is a closest-feasible-variant choice the corpus does
  not print.
- **Printed energy terms used:** W₂ (01:133), W_χ = χ₁e_kkΓ_ll + χ₂e_(ij)Γ_(ij) +
  χ₃e_[ij]Γ_[ij] (01:134), the chiral gate (e tensor, Γ pseudotensor; 01:114),
  the locking/potential sector normalized to the h25-audited quadratic operators
  (B4 row: Jω² = (α+β)k² + 4μ_c + m_V², "exactly KG"; h25_RESULTS.md §3), and the
  Q-5 orientation field with stiffness f_q (01:617).
- **The h25 cross-check (gate T2-G1), stated sharply:** the corpus's audited
  quadratic branch spectrum enumerates B1, B2±, B3, B4 and knot band-edges —
  **and no locking-stratum mediator branch**. The stratum whose gap is supposed
  to make the T2 tube (Q-2: T ~ M_gap²) has no printed bulk field in the
  enumerated content (GLUEBALL_CONTEXT_ANALYSIS.md §1, gap 1). Every internal-mode
  candidate below must therefore be built either from enumerated substrate
  branches (localization on the tube unprinted) or from a never-enumerated
  mediator field (the field itself unprinted). This single printed fact already
  bounds what the census can find.

## 3. T2-G1 — the mode-enumeration table

Worldsheet symmetry frame: unbroken group of the straight tube = O(2) rotations
about the axis × worldsheet Poincaré; P_ws = reflection through a plane
containing the axis (the transverse-parity used in torelon spectroscopy [IM]).
Tensor characters are the printed ones (01:114), mechanically verified (checks
A1/A2/A5).

| # | candidate mode | field content | printed evidence (file:line) or flag | O(2) × P_ws | worldsheet mass | status |
|---|---|---|---|---|---|---|
| 1 | **X₁, X₂ transverse translations** (2) | tube collective coordinates | tube existence + tension: 01:615 (Q-1/Q-2), 01:617 (Q-5); background: h21 §2 | O(2) vector (X₁ even, X₂ odd) | **0 exactly** (checks B1, C1–C3) | **PRINTED-DERIVABLE** |
| 2 | core axial rotor ϑ (rotation of core orientation texture about the axis) | locking-stratum orientation field, stiffness f_q | field printed 01:617 (Q-5); **existence as a tube-localized mode NOT printed → [CJ-new]**; gap stratum: 01:615 (Q-2, T ~ M_gap², M_gap symbol-only) | O(2) scalar, **P_ws-ODD (pseudoscalar — forced, check A5)** | 𝒪(M_gap) or f_q-scaled, **symbol-only** | [CJ-new] existence; parity forced *if* it exists |
| 3 | core-bound axial twist (B4-class descendant) | substrate micro-rotation φ_z | branch printed 01:147; audited operator h25 §3; **binding to the T2 core NOT printed → [CJ-new]** | O(2) scalar, **P_ws-ODD (pseudoscalar)** | branch gap (4μ_c+m_V²)/J — substrate stratum; the hadronic-scale (M_gap) reading is itself [CJ-new] | [CJ-new] localization; reduction validated (E3–E5) |
| 4 | core-bound relative-rotation doublet (B3-class, m = ±1 channels) | substrate ψ_⊥ (B3) | branch printed 01:147; h25 §3; binding [CJ-new] | O(2) vector (parity-even pair) | gapped, branch scale, symbolic | [CJ-new]; NON-axionic |
| 5 | breathing / radius mode | core profile deformation | **core profile nowhere printed → [CJ-new] entirely** | O(2) scalar, P_ws-EVEN | 𝒪(M_gap) symbolic | [CJ-new] |
| 6 | phason (T3) admixture | soft-sector phason, f ∈ [4, 60] MeV | T3 printed 01:615 (Q-2), 01:167 (II.H); **no printed dynamical cross-stratum coupling to T2** — the only printed cross-stratum relation is the stiffness law f_q = √ε_q·𝔪_Sk (V15.7); stratum discipline = the F-Q1 lesson | (would-be scalar) | — | **NOT-SUPPORTED-BY-PRINT** (any admixture [CJ-new]) |
| 7 | longitudinal u₃ | worldsheet reparametrization | standard bookkeeping | — | — | GAUGE (not physical) |
| 8 | ℤ₃ label k ∈ {1, 2} | discrete winding class 2πk/3 | 01:615 (Q-1) | — | — | discrete superselection label, not a mode |

**Census bottom line:** rows 2–5 — every internal candidate — carry a [CJ-new]
flag on the load-bearing step (existence/localization); row 6 is affirmatively
unsupported by print. Only row 1 is printed-derivable.

## 4. T2-G2 — the symbolic quadratic reduction (what was derived, and how)

All derivations are executable in `t2_census.py`; each numbered check is an
exact sympy identity (no numerics, no tuning).

1. **Parity characters of the printed energy (check A, mechanical).** With u a
   polar vector and φ axial (the printed reading that makes e a tensor and Γ a
   pseudotensor, 01:114 — both verified as polynomial identities, A1/A2): W₂ is
   parity-even (A3); **W_χ is parity-ODD (A4)**; W₄ ~ L⁴ even; b_P ~ ε(L)³ odd
   but enters the printed action squared (½Λ²b_P²) hence even; 𝒱(σ_P) even.
   **W_χ is therefore the corpus's unique parity-odd energy class** — the entire
   chirality hook must run through χ₁/χ₂/χ₃.
2. **Goldstone check (check B — the gate's required validation).** For the
   general two-component energy class containing the printed W₂ + W_χ structure
   at quadratic reduction (quadratic gradients, symmetric and χ-type
   antisymmetric cross-gradients, field–gradient couplings of the ε_ijkφ_k-in-e
   type, generic polynomial locking potential), the linearized fluctuation
   operator applied to the translation direction (∂_xψ⁰) equals ∂_x(EL[ψ⁰])
   identically — so on **any** static background solving the field equations,
   translations are exact zero modes. **The two transverse worldsheet modes are
   massless for any tube profile — printed or unprinted.** The validation the
   spec demands PASSES at operator level, ansatz-independently.
3. **Explicit X-sector reduction (check C).** On the Q-1-forced orientation halo
   (θ₀ = νφ, ν = 1/3, stiffness f_q), the collective-coordinate reduction gives
   the Nambu–Goldstone worldsheet action
   S_X = ∫dt dz [½μ_ℓ Ẋ_a² − ½σ_hal X′_a²], with **zero mass term**,
   σ_hal = πν²f_q² ln(R_IR/a_core), μ_ℓ = (J_m/f_q²)σ_hal, cone c_q² = f_q²/J_m.
   The tension reproduces the **printed form** σ = πf_q² ln κ_q (01:617) with the
   ν² = 1/9 winding bookkeeping absorbable into the unvalued ln κ_q — consistent
   with (and a small sharpening of) the context analysis's finding that π ln κ_q
   is load-bearing and never valued.
4. **Internal-mode reduction (check E, on the flagged ansatz).** On the
   [CJ-new] rigid-core ansatz (u_a = X_a(z,t)g(r), φ₃ = ϑ(z,t)h(r); winding axis
   aligned with the tube axis — a further [CJ-new] alignment assumption), the
   mechanical reduction of W₂-wryness + locking + W_χ gives:
   - **all X–ϑ cross terms, including every W_χ bilinear, vanish on angular
     integration** (E1) — the mechanical confirmation of the O(2) selection rule
     (check D1: no invariant ϑ–X bilinear exists);
   - **all background-wryness tadpoles vanish** (E2) — the straight tube is a
     consistent background at this order;
   - the internal-mode worldsheet action is
     S_ϑ = ∫dt dz ½[J𝓘₂ ϑ̇² − (α+β)𝓘₂ ϑ′² − (𝓚_⊥ + 𝓥_lock)ϑ²], with
     𝓘₂ = 2π∫h²r dr, 𝓚_⊥ = 2π∫[(β+γ)/2]h′(r)² r dr, and 𝓥_lock the locking-well
     integral (well shape [CJ-new]);
   - **validation against the audited branch spectrum:** in the uniform-profile
     limit the reduction reproduces the h25 B4 row **exactly**:
     Jω² = (α+β)k² + 4μ_c + m_V² (E3–E5). The reduction machinery is therefore
     anchored to the campaign's audited operators, not invented.
5. **Mass scale, symbol-only.** m_ϑ² = (𝓚_⊥ + 𝓥_lock)/(J𝓘₂). A mode below the
   bulk continuum exists iff the core well supports a bound state — unprinted
   [CJ-new]. On the locking-stratum reading its scale is ∝ M_gap (symbol-only in
   the entire archive); on the f_q reading, f_q-scaled. **No number is claimed at
   any grade** (seal q-θ: nothing here reaches even V.F, because no number is
   produced).

## 5. T2-G3 — parity classification and the census verdict

**The parity question, answered from print.** Under P_ws, any core-localized
mode that is an *axial rotation angle about the tube axis* (row 2's ϑ) or an
*axial micro-rotation component* (row 3's φ_z) flips sign — this is forced by
the printed axial-vector character (mechanical check A5), not by any new
assumption. So GUM's candidate worldsheet-axion **class** is well-defined: it is
the tube-core axial rotor/twist. Furthermore:

- the axion-type vertex ϑ·ε_ab∂X_a∂X_b — the coupling structure the worldsheet
  axion needs [IM: Athenodorou–Dubovsky–Luo–Teper; Axionic String Ansatz] — is
  O(2)×P_ws invariant (check D2), and
- among the printed energy classes **only W_χ can source it** (check A4:
  unique parity-odd class); its coefficient would be χ-class × (core overlap
  integral) [CJ-new].

**But no printed statement supplies the mode itself.** The three independent
blockers, each carried from the archaeology as SILENT verdicts:
(i) no locking-stratum mediator branch exists in the enumerated field content
(h25 §3 cross-check); (ii) no statement localizes any orientation/twist mode on
the T2 core, and the compact-loop/tube quantization machinery does not exist
(the h4/i4 moduli are properly-embedded lines); (iii) no cross-stratum chirality
transfer is printed — and the corpus's printed chirality consequences elsewhere
are uniformly *suppression* theorems (photon achirality II.3 and Ω-2, 01:151;
knot-matter achirality gate χ* = 0, 01:563), though none of them applies to the
T2 core (the core is SILENT, in both directions).

**Verdict (pre-registered categories):** **GUM tube = plain bosonic string as
printed** — two transverse modes, nothing else derivable — so the corpus's
closed-tube sector, taken as printed, **inherits the known 3+1D Isgur–Paton/
Nambu–Goto J^PC problems** in exactly the 0⁻ sector where the [IM] anchors put
X(2370). That is the finding. The constructive residue is filed at proposal
class: **the tube-core axion [CJ-new]/[DW, testbed-grade]** — a core-bound axial
rotor/twist mode, worldsheet-pseudoscalar by forced parity, axion vertex sourced
by the printed χ-class, mass ∝ M_gap (symbol-only). What would discharge it:
(a) exhibit the locking-stratum mediator field and its core-bound spectrum
(new construction), or (b) exhibit a core-localized mode of an enumerated branch
(B4-class) with a computed binding — either is corpus-side work; both are offers.

## 6. T2-G4 — register statement

Within-model throughout; nothing bears on nature. No mass number is produced
anywhere in this workstream (all scales symbolic), so the seal q-θ V.F protocol
has nothing to grade and nothing was tuned — the pipeline was frozen by the
roadmap before execution. The numerology hazard (𝔠 = 2.37 vs m_X ≈ 2.37 GeV) is
pre-empted in the header and cited nowhere as structure. SILENT-corpus verdicts
are carried, not upgraded. No stake is issued, no board clock moves, no grade
rises. The sentences "GUM predicted X(2370)", "X(2370) confirms GUM", "X(2370)
refutes GUM" are unavailable and absent.

## 7. Key numbers / symbols

| quantity | value | provenance |
|---|---|---|
| σ (printed tube tension) | 0.19 GeV² | [IM] inversion anchor, 01:617 (used only as form-comparison; no spectrum computed here) |
| f_q | 0.14–0.20 GeV | inversion output, 01:617 |
| winding ν | 1/3 (2π/3 frame winding) | 01:615 (Q-1, n = 3) |
| halo tension (derived form) | σ_hal = πν²f_q² ln(R_IR/a_core) | check C2; matches printed π f_q² ln κ_q with ν² absorbed into unvalued κ_q |
| transverse-mode mass | 0 (exact) | checks B1, C3 (Goldstone validation) |
| internal twist stiffness (longitudinal) | α+β | check E3 = h25 B4 row |
| internal twist gap (uniform limit) | (4μ_c + m_V²)/J | check E4 = h25 B4 row |
| tube-core-axion mass | ∝ M_gap, **symbol-only** (or f_q-scaled) | M_gap never valued in print (context §1); no number claimed |
| checks | 16/16 PASS | t2_results.json |

## 8. Honest caveats

1. **The rigid-core ansatz (g, h profiles, winding-axis alignment) is [CJ-new].**
   Its two mechanical conclusions that matter are ansatz-robust: the ϑ–X
   quadratic decoupling is the O(2) selection rule (check D1), and transverse
   masslessness is the operator identity (check B1) — neither depends on the
   profile. The extracted 𝓚_⊥/𝓥_lock integrals do depend on it and are quoted
   symbol-only.
2. **The halo normalization (A-T2-1)** reproduces the printed tension *form*;
   the ν² prefactor bookkeeping is absorbed into the unvalued κ_q. If the corpus
   ever values κ_q, this absorption becomes checkable.
3. **Two substrate levels are in play** (micropolar u/φ vs the Q-5 orientation
   field), because the corpus itself distributes the color mediator's roles
   across them without printing the mediator branch (h25 cross-check). The
   census is honest about which level each row lives on; a corpus-side
   unification of the two readings would sharpen rows 2–3 into one.
4. **The parity classification uses P_ws** (reflection through a plane
   containing the axis), the torelon-spectroscopy convention [IM]. Strict 3D
   parity is not a symmetry of the substrate at all (chiral constitutive
   choice, 01:114); any emergent parity label carries O(χ) corrections — which
   is precisely why a pseudoscalar mode here would be *natural* rather than
   fine-tuned, and equally why its absence from print is a genuine structural
   silence, not an oversight this workstream may fill.
5. **Cubic order was classified, not derived.** The axion vertex's existence is
   established as symmetry-allowed with W_χ its only possible printed source;
   its coefficient (χ-class × overlap) would require the unprinted core profile.
   No claim is made that it is nonzero.
6. **What this workstream does not adjudicate:** whether the closed loop of this
   tube reproduces any glueball-analog spectrum (T1's job); the seven-criteria
   scorecard (T3); the two-sector dictionary (T4). Promotion of the tube-core
   -axion proposal is the coordinator's call.

## 9. Deliverables

- `t2_census.py` — symbolic census: parity characters, Goldstone operator
  identity, halo reduction, O(2)×P_ws selection rules, rigid-core reduction
  with h25 B4 cross-check (deterministic, 2.6 s, 16/16 checks).
- `t2_results.json` — machine-readable: checks, 8-row census, symbolic
  worldsheet actions, gates, verdict.
- This memo. **F-T9-T2 status: census executed — plain bosonic string as
  printed (no pseudoscalar from printed structure: the finding); tube-core
  axion filed at [CJ-new]/[DW, testbed-grade] proposal class, mass ∝ M_gap
  symbol-only.**

— end of memo —
