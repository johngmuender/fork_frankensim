# T1 — Closed-tube spectrum from printed corpus parameters [F-T9-T1]

**Phase:** Tier 9 (Phase T), ROADMAP_v11_GLUEBALL.md workstream T1
(pre-registered; gates frozen before computation).
**Date:** 2026-08-15.  **Code:** `t1_spectrum.py` (deterministic, no RNG,
self-checking: SHO solver regression 5.5e-6, IP grid convergence 9.9e-7,
NG identity checks exact/4e-16; ~40 s wall).
**Outputs:** `t1_results.json`, `t1_fig.png`.
**Sources:** σ = 0.19 GeV² [IM-inversion anchor] from corpus Theorem Q-5
(σ = πf_q² ln κ_q, the only printed tension; see
`GLUEBALL_CONTEXT_ANALYSIS.md` §2); Isgur–Paton, PRD 31, 2910 (1985)
(lowest 0⁺⁺ = 1.52 GeV; snippet-pinned via Crede–Meyer arXiv:0812.0600 and
Mathieu–Kochelev–Vento arXiv:0810.4453); Nambu–Goto closed-string levels
per Athenodorou–Bringoltz–Teper arXiv:1007.4720; closed Regge slope α′/2
per Sonnenschein–Weissman arXiv:1507.01604; [IM] anchors per
`GLUEBALL_PAPER_DIGEST.md` / `t_context_agents.json`.

**Epistemic frame (binding).  Within-model; nothing here bears on nature.**
BESIII/lattice values are [IM] imported anchors.  Under seal q-θ, every
mass number below that confronts a hadron mass is reportable ONLY at V.F
grade: *dimensionally secure / structurally plausible / quantitatively
unclaimed*.  The pipeline was frozen by the roadmap before execution; no
parameter was tuned toward 2359/2376 MeV or any lattice value.
**Numerology pre-emption (T1-G4):** the corpus's dimensionless closure
benchmark 𝔠 = 2.37 ± 0.09 and m_X ≈ 2.37 GeV share digits by man-made-unit
coincidence only; that digit identity is cited nowhere in this workstream
as structure, and may never be.

---

## 1. Question and answer

**Question (pre-registered).**  Taking the corpus's T2 tube at face value
(σ = 0.19 GeV² [IM], √σ = 0.43589 GeV), what do the two standard
closed-flux-tube quantizations give for the lowest glueball-analog states,
and — the honest headline — does the corpus-parametrized closed tube put
its **lightest pseudoscalar** in the X(2370)/lattice-0⁻⁺ class
(m/√σ ≈ 5.4–5.9)?

**Answer: NO — in neither route, and the two routes miss in opposite
directions.**  (Either answer was pre-registered as a finding.)

1. **Route A (Isgur–Paton loop, validated to IP's own 1.52 GeV first):**
   lightest pseudoscalar 0⁻⁺ at m/√σ = **8.27** (systematics envelope
   7.56–9.73) — **40–70 % above** the window.  The IP-route states that DO
   fall in/near the window are the wrong-J^PC ones: the model's known
   pathological 1⁻⁺ orbital excitation (5.62) and the 2⁺⁺/2⁻⁺ phonon
   doublet (5.97).
2. **Route B (free Nambu–Goto closed string, exact):** the only
   pseudoscalar candidate (antisymmetric level-1 state) sits at
   m/√σ = **4.80**, degenerate with the 0⁺⁺ and 2⁺⁺ of the same level —
   **11–12 % below** the window, and with the ratio 0⁻⁺/0⁺⁺ = 1 instead of
   the lattice 1.50.  This is precisely the sector where the lattice says
   free NG fails and a **massive worldsheet axion** is required
   (Athenodorou–Dubovsky–Luo–Teper [IM]) — the 0⁻ caveat carried from the
   context analysis, and the exact question T2 interrogates on the GUM
   tube.
3. What the corpus tube DOES land on: the IP-route 0⁺⁺ head at
   m/√σ = 3.58 sits on the lattice scalar anchor m(0⁺⁺)/√σ ≈ 3.5 (within
   2.4 %; band 3.28–3.88).  Note honestly: the *dimensionless* 0⁺⁺ value is
   fixed by the validation calibration to IP's published 1.52 GeV, so this
   is IP's success re-expressed at the corpus tension, not a new GUM
   result; the corpus-specific content is that σ = 0.19 GeV² [IM] maps it
   to 1.56 GeV.
4. **V.F statement of record:** every number here is dimensionally secure
   (masses scale as √σ exactly), structurally plausible (both routes are
   the field's standard closed-tube treatments), and quantitatively
   unclaimed (scheme envelopes span tens of percent; no grade motion).

## 2. Method and provenance (frozen; deviations printed)

### 2.1 Route A — Isgur–Paton adiabatic loop quantization

Closed circular flux loop of radius ρ; transverse phonons of frequency
m/ρ, m ≥ 2 (m = 0 is the collective radial coordinate, m = 1 the spurious
translation), each carrying angular momentum ±m about the loop axis.
Adiabatic effective potential in the review-pinned rendition
(Crede–Meyer 0812.0600; Mathieu et al 0810.4453):

```
E_tot(ρ) = 2πσρ + c₀ + (M − 13/12)/ρ · (1 − exp[−f√σ ρ]),
M = Σ_m m(n_m⁺ + n_m⁻),   −13/12 = ζ-reg Σ_{m≥2} m,
```

with f the IP short-distance "fudge factor" (finite tube thickness).
Radial Schrödinger equation in ρ with the loop's own inertia
μ(ρ) = 2πσρ (the ring's kinetic energy is ½·2πσρ·ρ̇²); in dimensionless
variables x = √σρ every eigenvalue ε = E/√σ depends only on (M, n, L, f,
scheme) — the σ = 0.19 evaluation is exact rescaling, so nothing was
re-solved at the corpus tension (no tuning surface exists).

**Coordinator-owned amendments (printed by the code, verbatim in
`t1_results.json`):**

- **[A1]** IP's own radial-quantization details (operator ordering, c₀,
  f value) are not pinnable from accessible literature (all direct fetches
  egress-blocked; WebSearch snippets only — the context sweep's unverified
  list already flagged exactly these items).  Closest published variant
  implemented: review-pinned E_tot + variable-mass radial Schrödinger
  equation.  Central scheme = **Zhu–Kroemer ordering** (von Roos family),
  Dirichlet BCs, with the single cutoff **f\* = 1.9724** calibrated ONCE to
  IP's published 0⁺⁺ = 1.52 GeV at b = 0.18 GeV², then frozen.  The
  BenDaniel–Duke ordering validates at 0.62 % only in its large-f limit
  and is carried as the ordering systematic.  Calibrating to the model's
  own published benchmark is validation, not tuning: no X(2370)/lattice
  anchor entered the calibration.
- **[A2]** c₀ = 0: no accessible source values it for the glueball
  sector; a nonzero c₀ would be an un-cited second parameter.
- **[A3]** The −13/12 zero point counts one polarization pair per mode m;
  the J^PC census (below) uses the full in-plane/out-of-plane phonon
  content of a loop in 3D, which reproduces IP's signature features (the
  odd-J PC = +− "oddballs"; a 0⁻ in the low spectrum).  IP's internal
  polarization bookkeeping is not pinnable; the energy validation is
  unaffected.

**J^PC census (computed, not assumed).**  Operators derived from loop
geometry — P: a†(m,pol,c) → (−1)^m·(pol = z ? −1 : +1)·a†(m,pol,−c);
C (flux-orientation reversal): c → −c; J_z = Σ c·m; band heads J = |Λ|,
parity doublets (Λ-doubling) counted as distinct states:

| level | phonon content | J^PC (multiplicity) |
|---|---|---|
| M = 0 | none | 0⁺⁺ |
| M = 2 | one m = 2 | 2⁺⁺, 2⁻⁺, 2⁺⁻, 2⁻⁻ (1 each) |
| M = 3 | one m = 3 | 3⁻⁺, **3⁺⁻ (oddball)**, 3⁺⁺, 3⁻⁻ |
| M = 4 | two m = 2 / one m = 4 | **0⁻⁺ (1)**, 0⁺⁻ (1), 0⁺⁺ (2), 4⁺⁺(3), 4⁻⁺(2), 4⁺⁻(2), 4⁻⁻(3) |

The **lightest pseudoscalar is structural**: a 0⁻⁺ needs Λ = 0 with mixed
in-plane/out-of-plane polarization, first available at M = 4 (two m = 2
phonons, antisymmetric in polarization).  Since J ≥ |K| for the symmetric
top, no rotor band on a Λ = 2 head can produce J = 0 — the M = 4
assignment is robust to the rotor freedom.

### 2.2 Route B — free Nambu–Goto closed string (exact)

```
E²(N_L,N_R; q, ℓ) = (σℓ)² + 8πσ[(N_L+N_R)/2 − (D−2)/24] + (2πq/ℓ)²,
level matching N_L − N_R = q;  glueball reading: ℓ → 0, q = 0,
N_L = N_R = N, D = 4  ⇒  E² = 8πσ(N − 1/12).
```

Machine-precision identity checks (T1-G1b): formula reproduction 0.0
absolute error; D = 26 level-1 masslessness exact; Regge slope from level
differences = 1/(4πσ) to 4.0e-16 = α′_open/2 exactly (Sonnenschein–
Weissman closed slope).  Level content: N = 1 is (D−2)² = 4 states =
0⁺⁺ (trace) ⊕ 2⁺⁺ helicity ±2 (its J = 2 completion is the known
non-critical-D consistency gap) ⊕ **0⁻⁺ candidate** (antisymmetric
ε_ij α^i ᾱ^j, the axion-like state; its P/C reading follows the closed-
string glueball literature [IM] and is exactly where the lattice finds
free-NG failure).  N = 0 is the closed-string tachyon (E² = −2πσ/3 < 0),
excluded from the glueball reading as in all standard treatments —
printed, not hidden.

## 3. Validation (T1-G1) — before any GUM number entered

| check | requirement | measured | verdict |
|---|---|---|---|
| IP lowest 0⁺⁺ at b = 0.18 GeV² | ≤ 3 % of 1.52 GeV | 1.520000 GeV (1.2e-9 rel.; f\* = 1.9724, calibrated then frozen) | PASS |
| IP ordering cross-check (BDD, large-f) | — | 1.5294 GeV (0.62 %) | PASS (also within gate) |
| honesty line | — | f = 1 raw gives 1.722 GeV (13.3 % off): calibration was necessary and is printed | — |
| NG formula identities | machine precision | 0.0 abs; slope err 4.0e-16; D=26 massless exact | PASS |
| solver regression (SHO) | — | 5.5e-6 rel. | PASS |
| grid convergence | — | 9.9e-7 rel. | PASS |

## 4. Spectrum at the corpus tension (T1-G2) — zero fitted parameters

σ = 0.19 GeV² [IM]; √σ = 0.43589 GeV.  Bands = envelope over the frozen
scheme set {ZK × f ∈ [f\*/1.5, f\*, 1.5f\*]} ∪ {BDD × f ∈ [8, 32]} — the
models' own systematics (IP cutoff freedom + ordering freedom), nothing
fitted.  GeV values are V.F-graded throughout.

**Route A (Isgur–Paton):**

| state (model-assigned J^PC) | level | m/√σ | envelope | m [GeV, V.F] | band [GeV] |
|---|---|---|---|---|---|
| 0⁺⁺ ground | M=0 | 3.583 | 3.28–3.88 | 1.562 | 1.43–1.69 |
| 1⁻⁺ orbital (**known 3+1D pathology**) | L=1 | 5.616 | 5.43–5.82 | 2.448 | 2.37–2.54 |
| 2⁺⁺/2⁻⁺ (+2⁺⁻/2⁻⁻ flagged) | M=2 | 5.971 | 5.73–7.00 | 2.603 | 2.50–3.05 |
| 0⁺⁺\* radial | M=0,n=1 | 6.665 | 6.50–6.85 | 2.905 | 2.83–2.99 |
| 3⁻⁺/3⁺⁻ oddballs (+3⁺⁺/3⁻⁻) | M=3 | 7.135 | 6.65–8.44 | 3.110 | 2.90–3.68 |
| **0⁻⁺** (+0⁺⁻, 0⁺⁺′, 4's) | M=4 | **8.275** | 7.56–9.73 | 3.607 | 3.30–4.24 |

**Route B (free Nambu–Goto, D = 4):**

| level | content (model-assigned) | m/√σ | m [GeV, V.F] | σ/M² (short-string proxy) |
|---|---|---|---|---|
| N=0 | tachyon, E² = −2πσ/3 | — | excluded | — |
| N=1 | 0⁺⁺, 2⁺⁺, **0⁻⁺ cand.** (degenerate) | 4.800 | 2.092 | 0.043 |
| N=2 | 0⁺⁺, 2⁺⁺, 4⁺⁺ head, tower | 6.941 | 3.025 | 0.021 |
| N=3 | level-3 tower (J ≤ 6) | 8.562 | 3.732 | 0.014 |

Short-string validity: the glueball (contracted-loop) reading sits where
the NG derivative expansion is formally uncontrolled (winding expansion
parameter 8π(N−1/12)/(σℓ²) → ∞ as ℓ → 0); the printed σ/M² proxy uses
ℓ_eff = M/σ.  Lattice finds free NG nonetheless accurate at short lengths
**except the 0⁻ sector** [IM] — the caveat is load-bearing exactly on the
state this workstream asks about.

## 5. Confrontation table (T1-G3; V.F grade locked)

| quantity | IP route | NG route | [IM] anchor | reading |
|---|---|---|---|---|
| m(0⁺⁺)/√σ | 3.58 (3.28–3.88) | 4.80 | ≈ 3.5 (lattice) | IP on-anchor (inherited from validation); NG 37 % high |
| m(0⁻⁺)/√σ | **8.27 (7.56–9.73)** | **4.80** | 5.5–5.9 (lattice); 5.41–5.45 (X(2370)) | **IP 40–70 % HIGH; NG 11–12 % LOW; neither in window** |
| 0⁻⁺/0⁺⁺ | 2.31 (1.95–2.97) | 1.00 | 1.50 (lattice, both standard studies) | IP too large, NG too small — the ratio brackets but misses |
| 2⁺⁺ | 5.97 (5.73–7.00) | 4.80 | ~5.4–5.5 (lattice 2400/440, not pre-registered; noted) | near-window both sides |
| states actually in the 5.4–5.9 window | 1⁻⁺ (5.62, pathological), 2⁺⁺ edge (5.97) | none | — | the window is populated only by wrong-J^PC content |

**Pre-registered question, answered:** the corpus-parametrized closed tube
does **NOT** put its lightest pseudoscalar in the X(2370)/lattice-0⁻⁺
class.  Route A's 0⁻⁺ (8.27) is far above; Route B's candidate (4.80) is
below and level-degenerate with the scalar.  The two standard treatments
bracket the window from opposite sides and *both* fail the 0⁻⁺/0⁺⁺ = 1.50
ratio anchor — which is the known state of the art for axion-less closed
bosonic strings in 3+1D [IM].  **Within-model consequence (structural,
not graded):** a GUM closed-tube 0⁻⁺ in the X class is not delivered by
the printed tension + off-the-shelf quantization; it would require
tube-internal structure (a worldsheet pseudoscalar mode) — precisely what
T2's chirality census interrogates.  If T2 finds none, the corpus tube is
a plain bosonic string and inherits both failure directions above.

## 6. Gates (pre-registered; may not move)

| gate | requirement | measured | verdict |
|---|---|---|---|
| T1-G1 | IP reproduces published 0⁺⁺ ≈ 1.52 GeV at their σ to ≤ 3 % before GUM numbers; NG machine-precision identities | 1.2e-9 rel. (calibrated f\* = 1.97, frozen; BDD cross 0.62 %); NG identities ≤ 4.0e-16 | **PASS** (with printed amendments A1–A3: prescription pinned to review rendition, not the original paper) |
| T1-G2 | lowest 3+ states per route at σ = 0.19 GeV², model-assigned J^PC, bands printed, zero fitted parameters | 6 IP levels + 3 NG levels + tachyon disclosure; J^PC from computed census; envelope bands; parameters cited or amendment-flagged | **PASS** |
| T1-G3 | dimensionless table vs [IM] anchors; pre-registered pseudoscalar question answered honestly; V.F language | table in §5; answer: **NO in both routes** (miss in opposite directions); V.F locked | **PASS (finding = NO)** |
| T1-G4 | seal q-θ compliance; numerology pre-emption; within-model framing; no grade motion | header + §1.4; 𝔠 = 2.37 digit coincidence pre-empted, cited nowhere as structure; no stake, no forbidden sentences | **PASS** |

## 7. Key numbers

| quantity | value |
|---|---|
| σ (corpus, [IM]) / √σ | 0.19 GeV² / 0.43589 GeV |
| IP validation: f\*, m(0⁺⁺) at b = 0.18 | 1.972436, 1.520000 GeV (1.2e-9) |
| IP ε = m/√σ: 0⁺⁺, 1⁻⁺, 2⁺⁺/2⁻⁺, 0⁺⁺\*, 3's, 0⁻⁺ | 3.583, 5.616, 5.971, 6.665, 7.135, 8.275 |
| IP 0⁻⁺ envelope | 7.560–9.730 (m: 3.30–4.24 GeV, V.F) |
| IP 0⁻⁺/0⁺⁺ | 2.310 (1.949–2.967) |
| NG m/√σ: N = 1, 2, 3 | 4.7998, 6.9405, 8.5618 |
| NG N=1 mass at σ = 0.19 (V.F) | 2.092 GeV |
| NG 0⁻⁺/0⁺⁺ (level degeneracy) | 1.000 |
| anchors [IM] | X: 5.41–5.45; lat 0⁻⁺: 5.5–5.9; lat 0⁺⁺: 3.5; lat ratio: 1.50 |
| NG tachyon level | E² = −2πσ/3 (excluded, printed) |
| self-checks | SHO 5.5e-6; grid 9.9e-7; NG identities ≤ 4.0e-16 |

## 8. Honest caveats (T1-G4 register)

1. **The IP prescription is a reconstruction.**  The original PRD 31, 2910
   text was not accessible in this environment; the energy formula is
   review-pinned, but ordering, c₀, f, and the polarization bookkeeping
   are amendment-class (A1–A3).  The scheme envelope (tens of percent on
   excited states, e.g. 0⁻⁺/0⁺⁺ ∈ [1.95, 2.97]) is the honest price, and
   is why nothing here rises above V.F.  The *sign* of the headline
   finding (IP 0⁻⁺ far above the window) survives every scheme in the
   envelope; the finding is robust, its magnitude is not claimed.
2. **The IP-route 0⁺⁺ "hit" on the lattice 3.5 anchor is inherited.**
   Dimensionlessly it is IP's own published success (1.52/√0.18 = 3.58),
   frozen in by validation; evaluating at σ = 0.19 only rescales GeV
   values.  It is not evidence for the corpus.
3. **Rotor band-head corrections omitted** for Λ > 0 heads (symmetric-top
   Λ/(2I_⊥) + Λ²/(2I_z) terms); these shift high-Λ heads upward by a few
   to ~15 % and do not affect the Λ = 0 states carrying the headline
   (0⁺⁺, 0⁻⁺).
4. **Degeneracies are adiabatic artifacts.**  The M = 2 quadruplet
   (2⁺⁺/2⁻⁺/2⁺⁻/2⁻⁻) is exactly degenerate because E_tot depends only on
   M; real splittings (and the exotic 2⁻⁻/0⁺⁻ content, which the lattice
   puts heavy) are beyond the model — a known IP-class J^PC weakness,
   carried as printed.
5. **The 0⁻ caveat cuts both ways and is the bridge to T2.**  Free NG has
   no massive pseudoscalar mode; the lattice-established worldsheet axion
   is extra structure.  GUM's printed chirality (χ-couplings) is the
   corpus's only candidate source of such a mode on the T2 tube; whether
   it delivers one is T2's question, and nothing in T1 presumes the
   answer.
6. **No absolute-mass claim is made anywhere above V.F grade.**  All GeV
   numbers are √σ-rescalings of dimensionless model output with σ [IM];
   the corpus's own σ provenance is inversion, not derivation (κ_q never
   valued — context analysis §1.2).
7. **Forbidden-sentence compliance:** this workstream does not say and
   does not imply "GUM predicted X(2370)" or "X(2370) confirms/refutes
   GUM"; no stake is issued; no grade moves; all corpus-side implications
   are offers routed through T3/T4.
