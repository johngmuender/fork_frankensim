# U2 — Two-carrier transport sign test (F-T10-U2)

**Campaign:** Tier 10 (Phase U), ROADMAP_v12_JUNCTION.md workstream U2 (frozen
spec; gates pre-registered). **Date:** 2026-08-16. **Code:** `u2_transport.py`
(python3 + numpy, deterministic, runtime ~3 s; exact-arithmetic ledger via
`fractions.Fraction`). **Raw output:** `u2_results.json`. **Figure:**
`u2_fig.png` (matplotlib Agg). **Sources:** JUNCTION_CONTEXT_ANALYSIS.md §§1–3;
JUNCTION_PAPER_DIGEST.md §§1–2, 5; u_context_agents.json
(tier10-tube-junction-dynamics-transport-reader: printed-scale inventory,
transport-silence verification, K-1 dictionary evidence); corpus evidence as
cited by that agent — Primer ex. 15.1 chain corpus3/03:807, 811, 827 (snap +
minting, ~1e5 N × 1e-15 m ~ 100-MeV-class); corpus3/01:615, 617 (Q-6′, Q-5,
junction order conditional), 01:779 (kill registry); WS-K census
(corpus/WS-K-Q1-Q2-Charge-Census-Crisis-v0_1.md:11, 14, 17, 34–38 — K-1
dictionary, taping rule, reconnection rows); tier9-glueball/T2/RESULTS.md
(two massless transverse Goldstones — the light-side census).

**Epistemic frame (binding).** *Within-model; nothing here bears on nature.*
This is a **[CJ-new] structural-genericity test, not a prediction**: the corpus
prints statics (tensions, mass floors, junction order, termination/reconnection
legality) and **no transport law of any kind** (context §2 — the dynamics agent
verified zero printed statements on mobility, drag, rapidity, or momentum
partition; "Regge" has zero corpus hits). Every dynamical element below is
therefore a declared construction on printed *scale inputs* only.
**Sign/monotonicity ONLY — magnitudes are UNCLAIMABLE** (no α_B value, no
transport ratio, no rapidity slope is produced, and none could be graded above
V.F if it were). **Seal q-θ:** no hadron-mass-confronting number is produced
anywhere in this workstream; the GeV quantities that appear (W_V, W_J) are
*inputs assembled from printed windows*, not computed masses — the V.F protocol
("dimensionally secure / structurally plausible / quantitatively unclaimed") is
satisfied with the third clause doing all the work. **No fitting:** the three
STAR observables appear in this memo exactly once, as [IM] anchors in §6, and
enter no formula, no scan, no gate (mechanical scan, §1 gate U2-G4). **The
junction interpretation of the STAR data is itself contested** — CGC gluon
saturation reproduces the α_B systematics with valence quarks only
(García-Montero–Schlichting), neutron-skin geometry contaminates the peripheral
isobar trend, and strangeness asymmetry shifts the ratio (Ross–Lin); the rivals
are carried at equal weight (§6) and nothing below leans on the junction
reading of the data. **Numerology pre-emption (mandatory):** the corpus's
dimensionless closure number 𝔠 = 2.37 ± 0.09 and any GeV-valued number sharing
its digits are related by man-made-units coincidence only; that coincidence is
not structure, is not used anywhere below, and may not be cited as structure in
any reading of this memo. Grade language never rises; all corpus-side changes
are offers.

## 0. Verdict in one line

**Over the entire printed parameter window — M_knot ∈ [1.7, 172.5] GeV, f_q ∈
[0.14, 0.20] GeV, m_J/√σ ∈ {0, 0.1355, 0.39}, snap scale 100–140 MeV — the
junction-led channel is strictly the cheaper carrier of baryon number (min
margin W_V/W_J = 1.128 in the most conservative tube-length reading, 2.90 in
the compact reading, rising to ~343 at the top of the knot window), so for ANY
monotone-decreasing transport-cost function the toy produces B-transport ≥
Q-transport — the STAR-*shaped* sign — without the junction ever carrying a
charge: the exact-arithmetic ledger shows B moving across the gap with zero
original knots moving and net signed core count conserved, i.e. K-1
(knot-carried B) is COMPATIBLE with junction-led transport. Carrier = knots;
enforcer = taping rule/junction constraint; transporter = junction migration +
pair minting. Magnitudes unclaimable; the sign is the whole result.**

## 1. Gates (pre-registered; may not move)

| id | requirement (pre-registered) | measured | verdict |
|---|---|---|---|
| U2-G1 | toy pre-registered with both validation limits exact: minting disabled → junction channel closed → asymmetry exactly 0; tube+junction+minting cost → 0 → maximal junction dominance | (a) minting off: P_J == 0.0, B − Q == 0.0, B/Q == 1.0 by exact float equality across all 3 cost families × 3 rapidity gaps × 4 scale choices; (b) W_J = 0: P_J == 1.0 exactly (= the supremum of every normalized decreasing g, verified on a 101-pt probe), and in the sharp-suppression limit the junction share == 1.0 exactly (g(W_V) underflows to 0.0, P_J = 1.0) — maximal dominance | PASS |
| U2-G2 | sign result over the FULL printed windows with robustness map | 4D grid 241 × 61 × 3 × 9 (M_knot × f_q × m_J × E_snap), both tube-length readings: junction channel strictly cheaper at **100 % of grid points in both readings** (no sign flip anywhere in-window); min margin 1.128 (conservative) / 2.902 (compact); robustness map = `u2_fig.png`; ansatz-robustness verified mechanically at 42 768 (point × family × scale) combinations: dominance sign == cost sign everywhere, and B − Q = P_J ≥ 0 everywhere minting is on — the B ≥ Q **sign needs no cost model at all**, only channel-J openness + Q-blind minting | PASS |
| U2-G3 | reconciliation statement as computed: carrier/enforcer/transporter trichotomy — is K-1 compatible with junction-led transport? | exact rational-arithmetic ledger over all 2³ minted-flavor assignments: B moves with **0 original knots moving**; net signed core count 3 → 3 (+3 original, +3 minted knots, −3 minted antiknots); B_total = 1 before and after; B at destination = 1 (forced by the junction constraint/taping rule); net minted charge in the destination region = 0 for every assignment. Verdict: **COMPATIBLE** — carrier (knots, K-1 intact) / enforcer (taping rule — it *forces* the signed destination assembly) / transporter (junction migration + minting) are three distinct roles; §5 | PASS |
| U2-G4 | register: no STAR fitting — mechanical scan of deliverables for the anchor values used as targets; V.F on any mass; forbidden sentences absent | regex scan of `u2_transport.py`, `u2_results.json`, `RESULTS.md` for the three anchor values: **zero hits in code and JSON**; RESULTS.md hits only on the [IM]-anchor-context line in §6 (scan verifies the context mechanically; result stored in JSON `G4_anchor_scan.pass = true`); no mass claimed at any grade (no mass computed); forbidden sentences absent | PASS |

## 2. The toy, as pre-registered — and every [CJ-new] step, printed

**Printed inputs (the only numbers that enter):** σ = 0.19 GeV² [IM import],
√σ = 0.4359 GeV; M_knot window [m_Sk = 1.7 GeV floor, m̃_t-class 172.5 GeV];
f_q ∈ [0.14, 0.20] GeV; m_J ∈ {0, 0.1355, 0.39}·√σ (U1's pre-registered set —
0.1355 is the 2+1D lattice anchor [IM] with the dimensional-transfer caveat;
0.39·√σ = f_q-central, the corpus's only junction-adjacent scale); snap scale
E_snap ∈ [100, 140] MeV — the **printed 100-MeV-class scale from Primer
ex. 15.1**, whose chain has the student compute ~1e5 N × 1e-15 m ~ 100
MeV-class against the 140 MeV pion (corpus3/03:807, 811, 827: "a snapped cord
has two fresh ends, which must terminate on fresh fractional knots… you pay
for two new ones"). Junction order n = 3 tubes — **Primer-asserted (03:811),
paper-conditional on WS1-H3 (01:615, 779), execution record absent**: the
F-K0-1 seam is carried, not resolved, and the whole construction inherits it.

**Channel inertias ([CJ-new] construction c2):**

- Channel V (valence-led): the three original knots traverse the gap;
  transported inertia **W_V = 3·M_knot**. Q rides only on knots (V.D winding;
  K-1), so channel V co-transports the baryon's charge with its B.
- Channel J (junction-led): the junction + three tube stubs migrate and pay
  2× pair-minting at destination; **W_J = m_J + 3σℓ_tube + 2·E_snap**. The
  corpus prints no structure size, so BOTH printed hadronic-length readings
  are scanned: *compact* ℓ_tube = 1/√σ (0.45 fm) and *conservative*
  ℓ_tube = 1/f_q (0.99–1.41 fm — the longest printed hadronic length, which
  maximizes W_J and is the hardest case for the junction channel).

**Suppression ansatz ([CJ-new] c1), printed:** the probability of moving
transported inertia W across a rapidity gap dy is P(W) = g(W) with g ANY
normalized monotone-decreasing cost function, g(0) = 1. Three families are
run — Boltzmann-type g = exp(−W·dy/Λ), power-law g = (1 + W·dy/Λ)^(−p),
Gaussian g = exp(−(W·dy/Λ)²) — with the [CJ-new] scale Λ scanned over
{0.5, 1, 5, 50} GeV and dy over {0.5, 1, 3} **only to demonstrate that the
sign never depends on them**. No value of Λ, dy, or p is claimed or claimable.

**Channel bookkeeping ([CJ-new] c4):** equal prior channel weights, no tuning.
Per unit of transported B, channel V delivers one unit of knot charge; channel
J delivers zero *net* charge (minting is local and pairwise net-zero at the
destination; the ensemble is C-blind). So B_trans ∝ g(W_V) + g(W_J),
Q_trans ∝ g(W_V), and the asymmetry sign is sign(g(W_J)).

## 3. U2-G1 — validation limits, exact in the code

1. **Minting disabled → junction channel closed.** A snapped tube must
   terminate on a fresh fractional core (the taping rule); with minting off,
   channel J cannot deliver, so P_J is set identically 0 and the code verifies
   **B − Q == 0.0 and B/Q == 1.0 by exact float equality** at the window
   corner, for all three g-families × all dy × all Λ (36 combinations). The
   toy is exactly valence-only in this limit.
2. **Tube + junction + minting cost → 0 → maximal junction dominance.** With
   W_J = 0, P_J = g(0) **== 1.0 exactly**, which is the supremum of every
   normalized monotone-decreasing g (verified against a 101-point probe of
   g(w) for each family/scale); in the sharp-suppression limit (Λ = 1e-3 GeV,
   exponential family) g(W_V) underflows to exactly 0.0 and the junction
   share P_J/(P_V + P_J) **== 1.0 exactly** — the junction channel carries
   all transported B.

## 4. U2-G2 — the sign over the full printed windows

**Grid:** M_knot 241 log-spaced × f_q 61 × m_J 3 × E_snap 9, both tube-length
readings (396 927 points per reading; margins computed exactly as
W_V/W_J ratios, no ansatz needed for the sign).

| reading | W_J range (GeV) | W_V range (GeV) | min margin W_V/W_J | max margin | junction-dominant fraction |
|---|---|---|---|---|---|
| conservative (ℓ = 1/f_q) | 3.050 – 4.521 | 5.1 – 517.5 | **1.128** (at M_knot = 1.7, f_q = 0.14, m_J = 0.39√σ, E_snap = 0.14) | 169.7 | **100 %** |
| compact (ℓ = 1/√σ) | 1.508 – 1.758 | 5.1 – 517.5 | **2.902** (same corner) | 343.2 | **100 %** |

**Sign result.** (i) *B ≥ Q:* B_trans − Q_trans = g(W_J) > 0 at every grid
point with minting enabled, for every cost family — this sign requires **no
cost model at all**, only that channel J is open (minting legal) and minting
is Q-blind. (ii) *Junction dominance* (the STAR-cartoon-relevant statement,
g(W_J) > g(W_V)): holds wherever W_J < W_V, which is **everywhere in the
printed window, in both readings** — and for ANY monotone-decreasing g, since
W_J < W_V ⟺ g(W_J) ≥ g(W_J′) ordering is preserved by monotonicity. This was
additionally verified mechanically at 42 768 point × family × scale
combinations (dominance sign == cost sign at every one; `u2_results.json`
`ansatz_robustness`).

**Where the sign WOULD flip (honest fragility edge).** At the hardest corner
(M_knot at the m_Sk floor, m_J and E_snap maximal) the flip boundary sits at
ℓ_flip = 8.16 GeV⁻¹ = 1.61 fm of migrating tube-stub length — only a factor
**×1.14** above the conservative reading's 1/f_q = 1.41 fm (×3.56 above the
compact reading). A junction-led mode that must drag ≳1.6 fm of tube through
the gap at the corner point would lose its cost advantage against the lightest
printed knots. The corner margin is thin and the tube-length reading is
[CJ-new]; the robust statement is the *window-wide* sign plus the strong
growth of the margin with M_knot (linear in M_knot; ~343× at the top of the
window). Figure `u2_fig.png`: panel A — worst-case margin map over
(M_knot, f_q); panel B — margin vs M_knot at the hardest corner, both
readings, with the flip boundary marked.

## 5. U2-G3 — the reconciliation ledger and the trichotomy verdict

Exact rational arithmetic (`fractions.Fraction`), all 2³ minted-flavor
assignments {u, d}³ enumerated:

- **Event:** junction + three tube stubs migrate; each tube snaps at the
  destination, minting one knot–antiknot pair locally; the junction
  constraint (taping rule: a fractional line terminates only on a fractional
  core; thirds sum to integers) forces the three destination tube ends onto
  the three minted **knots** — signed assembly B = +1 at destination is
  *forced by the enforcer*, not chosen.
- **Ledger, every assignment:** original knots moved = **0**; signed core
  count 3 → 3 (+3 original at beam, +3 minted knots −3 minted antiknots at
  destination); B_total = 1 before and after (K-1 exact throughout);
  B(destination) = 1; **net minted charge in the destination region = 0 for
  every flavor assignment** (pairs are local and net-zero — minting is
  Q-blind), while the destination baryon's own charge varies with assignment
  (2 for uuu through −1 for ddd) exactly as pair-symmetric assembly permits.
- **Trichotomy verdict, as computed:** **carrier = knots** (K-1 intact at
  every step; the junction never carries a unit of B or any charge in the
  ledger); **enforcer = the web/taping rule + junction constraint** (it is
  exactly this constraint that forces B = +1 at destination); **transporter =
  the junction+tube migration mode** (B crosses the rapidity interval with
  zero original knots moving). **K-1 (knot-carried B) is COMPATIBLE with
  junction-led transport** — the apparent GUM-vs-STAR tension dissolves into
  three distinct structural roles. Status: [CJ-new] within-model
  construction; the corpus prints the *legality* of every ingredient
  (snap+minting existence, reconnection net-conservation, taping rule) but no
  dynamics — this ledger tests structural consistency, not a prediction.

## 6. Confrontation context ([IM] anchors only — not targets, not fitted)

STAR [IM] anchors, for orientation only:

- [IM] anchor: isobar ⟨B⟩/ΔQ × ΔZ/A = 1.84 ± 0.02 ± 0.09 ± 0.16 (0–10 % central)
- [IM] anchor: α_B(γ+Au) = 1.04 ± 0.22 (junction Regge window [0.42, 1])
- [IM] anchor: α_B(Au+Au) = 0.64 ± 0.05, centrality-independent

None of these numbers enters any formula, scan, or gate above (mechanically
verified, gate U2-G4). The toy's sign statement is *structurally shaped like* the STAR
sign (B transport ≥ Q transport via a light junction-led mode) — and no more
than that: **the toy cannot discriminate the junction interpretation from its
rivals**, because CGC gluon saturation reproduces the same α_B systematics
with valence quarks only (García-Montero–Schlichting), the peripheral isobar
trend is part neutron-skin geometry (TRENTO reproduces the trend without
junctions), and strangeness asymmetry shifts the isobar ratio (Ross–Lin). A
sign-only toy that agrees with a contested sign confirms nothing and is
claimed to confirm nothing.

## 7. Key numbers

| quantity | value | status |
|---|---|---|
| σ (tension) | 0.19 GeV² | [IM] import, printed |
| √σ | 0.4359 GeV | arithmetic |
| M_knot window | 1.7 – 172.5 GeV | printed (m_Sk floor unprovenanced in GeV; m̃_t-class input) |
| f_q window | 0.14 – 0.20 GeV | printed ([CAL]-class inversion) |
| m_J set | {0, 0.0591, 0.1700} GeV = {0, 0.1355, 0.39}·√σ | pre-registered (U1); 0.1355 = lattice [IM], transfer caveat |
| E_snap window | 0.100 – 0.140 GeV | printed 100-MeV-class (Primer ex. 15.1 chain) |
| W_V range | 5.1 – 517.5 GeV | toy input assembly |
| W_J range | 1.508 – 4.521 GeV (both readings pooled) | toy input assembly, [CJ-new] structure |
| min margin W_V/W_J (conservative / compact) | 1.128 / 2.902 | computed, sign-relevant only |
| max margin | 343.2 (compact, top of knot window) | computed, sign-relevant only |
| junction-dominant fraction of 4D grid | 100 % (both readings) | computed |
| flip boundary at hardest corner | ℓ_flip = 1.61 fm (×1.14 above 1/f_q; ×3.56 above 1/√σ) | computed fragility edge |
| ansatz-robustness checks | 42 768 combinations, dominance sign == cost sign at all | computed |
| ledger: original knots moved / B conserved / net minted Q at destination | 0 / exact (1 → 1) / 0, all 8 flavor assignments | exact rational arithmetic |
| B-vs-Q asymmetry magnitude | **UNCLAIMABLE** | corpus prints no transport law |

## 8. Honest caveats

1. **Everything dynamical here is [CJ-new].** The corpus prints no transport
   law, no mobility, no rate for minting or reconnection, and no junction
   energetics; W_J's very form (junction + finite tube stub + 2 snaps) is a
   constructed reading. The result is a *structural-genericity* statement
   about printed scales, not a corpus prediction — and it could not have
   failed to be one: had the sign come out opposite, that too would have been
   a [CJ-new] toy result.
2. **The conservative-reading corner margin is thin (1.128).** A modestly
   longer migrating structure (≳1.6 fm at the corner) flips the corner sign.
   The window-wide claim survives because the margin grows linearly in
   M_knot, but the corner fragility is real and printed in §4.
3. **Junction order n = 3 is Primer-asserted, paper-conditional** (WS1-H3
   execution record absent; F-K0-1 seam). The 2× minting count and the
   three-tube structure both inherit this provenance caveat.
4. **K-1 itself is v1-worksheet-stratum only** (zero operative-edition hits
   for a B dictionary — the unaffixed delta-pack stubs are a filed fold
   debt); the trichotomy verdict cites K-1 at census grade, not paper grade.
5. **The toy has no kinematics.** Rapidity appears only as an abstract gap
   variable inside a monotone cost; there is no collision geometry, no
   energy dependence, no centrality — none of the structure that makes the
   STAR measurements discriminating. Agreement in sign is a consistency
   observation at the cartoon level.
6. **The rivals are live.** CGC saturation, neutron skin, and strangeness
   confounds each reproduce part of the anchor systematics without
   junctions; this workstream neither adjudicates between them nor could.
7. **Ensemble vs event.** "Channel J transports zero net charge" is an
   ensemble-mean statement (C-blind minting); event-by-event the destination
   baryon carries assignment-dependent charge, compensated locally by the
   minted antiknots. The toy addresses net transport only, which is what the
   [IM] anchor observables measure.

## Finding

**F-T10-U2.** Within-model, [CJ-new] construction, sign-only: the corpus's
printed statics (knots ≥ 1.7 GeV, tubes at σ = 0.19 GeV² with two massless
transverse Goldstones, junctions as unweighted constraint points, minting at
the printed 100-MeV-class snap scale) *generically* make the junction-led
channel the cheaper carrier of B across a rapidity gap over the entire printed
parameter window, for any monotone-decreasing transport cost — producing
B-transport ≥ Q-transport with the junction carrying no charge, and doing so
in exact consistency with K-1: B is knot-carried, taping-rule-enforced,
junction-transported. Magnitudes unclaimable; no STAR observable fitted or
confronted beyond [IM]-anchor orientation; the junction interpretation of the
anchors remains contested and this finding does not bear on that contest.
