#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
T3 — Seven-criteria scorecard + the closure obligation.  [F-T9-T3]
ROADMAP_v11_GLUEBALL.md workstream T3 (frozen spec; gates pre-registered).

This workstream is analysis, not heavy compute: the script assembles the
machine-readable scorecard (t3_scorecard.json), validates that every row
verdict is one of the four PRE-REGISTERED categories, that every corpus
cell carries file:line evidence, and computes the two small numbers the
relic-census note uses (Gamma/m; branching sums).  Deterministic, no RNG.

Epistemic register: within-model; nothing here bears on nature.  BESIII /
lattice values are [IM] imported anchors.  Seal q-theta: every mass number
confronting a hadron mass is V.F-graded (dimensionally secure /
structurally plausible / quantitatively unclaimed).  Numerology
pre-emption: the corpus's dimensionless c-frak = 2.37 +/- 0.09 and the
imported m_X ~ 2.37 GeV share digits by man-made-units coincidence; that
coincidence is not structure and is cited nowhere as structure.
"""

import json, math, sys

CATEGORIES = ["SUPPORTED-STRUCTURALLY", "V.F-CLASS (number, unclaimed)",
              "SILENT-NEEDS-NEW-WORK", "TENSION"]

# ------------------------------------------------------------------ [IM] anchors
IM = {
    "m_X_combined_MeV": {"value": "2359 (+13/-14)", "source": "arXiv:2605.26495 via Letter arXiv:2607.20366; GLUEBALL_PAPER_DIGEST.md Sec.1"},
    "Gamma_X_MeV": {"value": "170 (+44/-29)", "source": "same"},
    "JPC": {"value": "0-+ (>9.8 sigma over weakest alternative 2-+; >10.8 over 1++)", "source": "PRL 132, 181901 (2024); digest Sec.3"},
    "BR_production": {"value": "B[J/psi->gamma X] > 1e-3 (est.); LQCD pure-G 2.31(90)e-4 + mixing -> 1e-3-class", "source": "digest Sec.3-4; Gui et al 2019; arXiv:2605.01757"},
    "flavor_singlet": {"value": "K*(892)Kbar null: B < 2.7e-6, R < 0.081 (90% CL); partial width < 2 MeV", "source": "arXiv:2607.20366 (this Letter); digest Sec.1-2"},
    "narrow_partial_widths": {"value": "quasi-two-body modes at B ~ 1-10% => few-MeV partial widths (sqrt-OZI, Robson 1977)", "source": "digest Sec.4 crit.6"},
    "gamma_omega_phi": {"value": "B(gamma omega) < 0.04e-6, B(gamma phi) < 0.11e-6 vs eta(1405) ~ 3.5 both", "source": "digest Sec.3 [29]; PRD 111, 052011 (2025)"},
    "lattice_0mp_window_GeV": {"value": "2.3-3.0 (quenched); continuum 0-+ 2.395(14) (Gui19); MP99/Chen06 2.56-2.59", "source": "t_context_agents.json glueball-theory-literature"},
    "anchor_ratios": {"value": "m_X/sqrt(sigma) = 5.41-5.45; lattice 0-+/0++ ~ 1.50; m(0++)/sqrt(sigma) ~ 3.5; m(0-+)/sqrt(sigma) ~ 5.5-5.9", "source": "GLUEBALL_CONTEXT_ANALYSIS.md Sec.2-3"},
}

# small computed numbers for the census note (T3-G3)
gamma_over_m = 170.0 / 2359.0
IM["Gamma_over_m"] = {"value": round(gamma_over_m, 4),
                      "source": "arithmetic on [IM] 170/2359 (campaign-side)"}

# ------------------------------------------------------------------ T1/T2 computed cells
# T1/T2 run in parallel; their RESULTS.md are read before finalizing.
# Cells below are filled from T1/RESULTS.md and T2/RESULTS.md.
import os
HERE = os.path.dirname(os.path.abspath(__file__))
T1_PATH = os.path.join(HERE, "..", "T1", "RESULTS.md")
T2_PATH = os.path.join(HERE, "..", "T2", "RESULTS.md")
T1_PRESENT = os.path.exists(T1_PATH)
T2_PRESENT = os.path.exists(T2_PATH)

T1_MASS_CELL = ("PENDING-T1" if not T1_PRESENT else
    "T1 (RESULTS.md Secs.1,4,5; t1_results.json; V.F-graded): validation-first PASS "
    "(IP reproduces its own published 0++ = 1.52 GeV at b = 0.18 GeV^2 to 1.2e-9 after "
    "one-time cutoff calibration f* = 1.97, frozen before the corpus tension entered; "
    "NG identities machine-precision). At sigma = 0.19 GeV^2 [IM] the pre-registered "
    "question is answered NO IN BOTH ROUTES, missing in OPPOSITE directions: IP lightest "
    "0-+ at m/sqrt(sigma) = 8.27 (envelope 7.56-9.73; 3.61 GeV-class) is 40-70% ABOVE "
    "the window, NG's only 0-+ candidate at 4.80 (2.09 GeV-class) is 11-12% BELOW and "
    "level-degenerate with 0++/2++ (ratio 0-+/0++ = 1.00 vs lattice 1.50; IP ratio 2.31); "
    "anchors: m_X/sqrt(sigma) = 5.41-5.45, lattice 0-+/sqrt(sigma) = 5.5-5.9. The 5.4-5.9 "
    "window is populated only by wrong-J^PC content (IP's known-pathological 1-+ orbital "
    "at 5.62; the 2++ edge at 5.97). The IP 0++ head (3.58 vs lattice 3.5) is inherited "
    "from validation, not corpus evidence (T1 Sec.8.2). The miss is localized in exactly "
    "the 0- sector where lattice torelon spectroscopy requires the massive worldsheet "
    "axion [IM: ADLT]. All numbers dimensionally secure / structurally plausible / "
    "quantitatively UNCLAIMED (V.F; seal q-theta)")

T1_JPC_CELL = ("PENDING-T1" if not T1_PRESENT else
    "T1: J^PC labels are MODEL-ASSIGNED ([IM] IP phonon census computed from loop "
    "geometry / NG level content), not corpus-derived; free NG cannot split 0-+ from "
    "0++ (N=1 degeneracy) and IP's known 3+1D pathologies surface as computed (low 1-+ "
    "orbital; exact adiabatic degeneracies) — both routes carry the printed 0- caveat "
    "[IM: ADLT]")

T2_JPC_CELL = ("T2 (RESULTS.md Sec.0,5): the GUM tube AS PRINTED is a plain bosonic "
    "string — only the two transverse Goldstones derivable (validated massless); NO "
    "pseudoscalar worldsheet mode emerges from printed structure, so the closed-tube "
    "sector inherits the known 3+1D IP/NG J^PC problems in exactly the 0- sector; "
    "the tube-core axion (forced worldsheet-pseudoscalar IF it exists; W_chi the unique "
    "printed parity-odd source) is filed at [CJ-new]/[DW, testbed-grade], mass ~ M_gap symbol-only")

# ------------------------------------------------------------------ the scorecard
SCORECARD = [
 {
  "criterion": "1. mass (LQCD 0-+ window 2.3-3.0 GeV; X at 2359 MeV [IM])",
  "corpus_as_printed": (
    "NO hadron mass is printed anywhere; the color sector's only quantitative anchors are "
    "sigma = 0.19 GeV^2 [IM import, 'inversion, not derivation'] and f_q = 0.14-0.20 GeV "
    "(corpus3/01-GUM-Omega-Paper-v4.3-ext.md:617, Thm Q-5); T2 tension symbolic T ~ M_gap^2 "
    "with M_gap never valued (01:615, Q-2); F-Q2-prime's 'hadronic spectrum' constraint is "
    "anchor-free — no state list (01:617, 767); m_X/sqrt(sigma) = 5.41-5.45 is campaign "
    "arithmetic — the corpus never forms the ratio (t_context_agents.json mass-anchor reader); "
    "the glueball-analog class itself is named once, v1 dated record only, zero quantitative "
    "content (corpus/WS-D-Q1-Q2-Relic-Census-v0_1.md:11)"),
  "T1_T2_computed": T1_MASS_CELL,
  "verdict": "TENSION" if T1_PRESENT else "SILENT-NEEDS-NEW-WORK",
  "verdict_note": ("TENSION at proposal class, V.F language: the corpus prints no mass, so no "
    "printed claim is hit — what misses the anchor window is the NEW [CJ]-class construction "
    "(printed sigma + frozen axion-less [IM] quantizations), whose lightest pseudoscalar lands "
    "outside 5.4-5.9 in both routes while the scalar channel lands on lattice 3.5; the tension "
    "is localized in the 0- sector and its identified repair is T2's [CJ-new] tube-core axion; "
    "no kill fires (no stake/row exists), nothing is refuted, numbers quantitatively unclaimed "
    "(seal q-theta V.F)" if T1_PRESENT else
    "PENDING-T1 at assembly time — amendment printed"),
 },
 {
  "criterion": "2. spin-parity 0-+ (>9.8 sigma [IM])",
  "corpus_as_printed": (
    "No parity operator, C operation, or J^PC label is printed for ANY composite or tube state "
    "(selection-rules agent, negative sweep; corpus3/01:815, 831); the only composite "
    "quantum-number result is the B=2 deuteron analog J=0-forbidden/J=1-ground (01:831 App F.5; "
    "sector-conditional post F-R17, theory-audit/i4_topology.md:361-393); printed rules deliver "
    "for a knot-free loop exactly 'flavor-singlet boson, J^PC otherwise UNDETERMINED' (ladder at "
    "K=0, corpus3/02:980 + 01:577); the h4/i4 monodromy machinery covers properly-embedded "
    "lines, not compact loops (theory-audit/h4_completion.md:63-99); strict P is not even a "
    "substrate symmetry (chiral constitutive choice, 01:114)"),
  "T1_T2_computed": (T1_JPC_CELL + " || " + T2_JPC_CELL),
  "verdict": "SILENT-NEEDS-NEW-WORK",
  "verdict_note": ("corpus silent on J^PC machinery for loops; T2 sharpens the silence into a "
    "structural finding: as printed the tube CANNOT supply the 0- sector the criterion needs "
    "(plain bosonic string), and the candidate repair (tube-core axion) is [CJ-new] — a "
    "TENSION-flavored silence carried at model grade, not upgraded to TENSION because the "
    "corpus asserts nothing that the criterion contradicts"),
 },
 {
  "criterion": "3. production rate (B[J/psi->gamma X] > 1e-3 [IM]; LQCD 2.31e-4 + mixing)",
  "corpus_as_printed": (
    "Nothing. No J/psi analog, no charmonium, no radiative-production machinery, no gg-fusion "
    "analog; 'gluon' appears twice, both popular-history prose (corpus3/03:80; archaeologist "
    "zero-hit sweep); the h25-audited branch spectrum contains NO locking-stratum mediator "
    "branch at all (theory-audit/h25_RESULTS.md:53-59) — the color mediator is invoked in VII.J "
    "but absent from the enumerated field content; no mixing machinery of any kind "
    "(no eta/eta-prime, no charm sector: archaeologist silences)"),
  "T1_T2_computed": "neither T1 nor T2 computes production (out of both frozen scopes by design)",
  "verdict": "SILENT-NEEDS-NEW-WORK",
  "verdict_note": "the [IM] structural corollary cuts deep here: the 0-+ glueball is only ever "
    "observable as the glueball-dominant member of a G-ccbar mixture (Letter/2605.01757), and "
    "GUM has no mixing machinery and no heavy-quarkonium sector to mix with",
 },
 {
  "criterion": "4. eta_c decay-pattern similarity (same J^PC, gluon-mediated, no dominant mode [IM])",
  "corpus_as_printed": (
    "No decay machinery of any kind exists for line-neutral composites: no hadronic decay "
    "table, no channel selection rules (selection-rules agent exhaustive inventory — the "
    "printed rules are substrate-level: Omega-1 01:151, VIII.D.1/2 01:649, VIII-prime.1 "
    "01:683, none forbids or weights a composite->composite channel); the sole decay-adjacent "
    "sentence for this class: 'stratum composites (glueball-analogs) -> decay to "
    "hadrons/phasons [check]' — instability asserted, zero channel content "
    "(corpus/WS-D-Q1-Q2-Relic-Census-v0_1.md:11, v1 dated record)"),
  "T1_T2_computed": "not computed by T1/T2 (out of scope); no decay-width machinery exists to reuse",
  "verdict": "SILENT-NEEDS-NEW-WORK",
 },
 {
  "criterion": "5. flavor-singlet (K*Kbar null; first flavor-singlet light hadron > 1 GeV [IM])",
  "corpus_as_printed": (
    "The one criterion the printed structure supports: every printed flavor label attaches to "
    "knot cores — base scales m-tilde (corpus3/01-GUM-Omega-Paper-v4.3-ext.md:617), family "
    "class p (corpus3/01-GUM-Omega-Paper-v4.3-ext.md:587-595, VII.H), tube-alignment zeta "
    "(01:617) — a knot-free closed tube has no cores, hence no flavor labels: flavor-singlet "
    "BY CONSTRUCTION ([DF-structural reading]; selection-rules agent deliverable 4); caveats "
    "printed with it: no flavor SYMMETRY exists in print (no SU(3) analog, no isospin "
    "multiplets — 'flavor' occurs once in the paper, an analog-platform figure caption 01:867), "
    "and the discriminating Lipkin-style test itself (0-+ singlet forbidden -> K*Kbar) has NO "
    "in-model counterpart — no G-parity analog is printed; obstacle carried: Course 20.1's "
    "completeness word 'finite-energy states are exactly the line-neutral composites: "
    "knot-antiknot pairs (mesons) and n-knot junction closures (baryons)' "
    "(corpus3/02-The-Substrate-Course-v3-ext.md:1088, verified this session) reads literally "
    "as excluding the knot-free loop"),
  "T1_T2_computed": ("T2 mode census row 8: the Z3 label k in {1,2} is a discrete superselection "
    "label on the tube, not a flavor label (T2/RESULTS.md Sec.3); T1 spectra carry no flavor "
    "structure at all (consistent: nothing to carry)"),
  "verdict": "SUPPORTED-STRUCTURALLY",
  "verdict_note": "supported at the label-absence level ONLY — the corpus can say 'no flavor "
    "content', it cannot reproduce the K*Kbar selection-rule test that made the [IM] "
    "measurement decisive; and the 'exactly' sentence needs a scope reading or amendment "
    "(offer-class, see T4)",
 },
 {
  "criterion": "6. narrow partial widths (few MeV per exclusive mode; sqrt-OZI [IM])",
  "corpus_as_printed": (
    "No width machinery for composites at any grade (WS-N RN-3 certifies 'no binding energy, "
    "level, or radius is computed anywhere in Sec. VII.J' — corpus/WS-N-Hadronic-Pass-P7-v0_1.md:9); "
    "no OZI analog (no quark-line rule exists — no quark lines exist); the only closed-loop "
    "decay dynamics ever printed is the T3 web-loop obituary — radiate phasons, die in tens of "
    "periods (corpus/WS-D-Q1-Q2-Relic-Census-v0_1.md:19, (c1)) — the WRONG stratum: the T2 "
    "closed-loop lifetime is nowhere priced (archaeologist silences)"),
  "T1_T2_computed": "not computed by T1/T2 (widths out of both frozen scopes)",
  "verdict": "SILENT-NEEDS-NEW-WORK",
  "verdict_note": "GUM-specific hook noted at proposal class only: the census line's phason "
    "channel would be an invisible partial width with no QCD counterpart — unmeasured, "
    "unconstructed, [CJ]-class if ever built",
 },
 {
  "criterion": "7. gamma-omega / gamma-phi suppression (both < 0.11e-6 [IM])",
  "corpus_as_printed": (
    "Nothing connects any composite state to a radiative decay channel; no vector-meson "
    "spectrum exists (no omega/phi analogs, no quarkonium content tags — archaeologist "
    "silences: no meson spectrum anywhere); the photon sector's printed theorems are "
    "substrate-level suppression/creation results (Omega-1/Omega-2, "
    "corpus3/01-GUM-Omega-Paper-v4.3-ext.md:151) that say nothing about hadronic radiative "
    "transitions (selection-rules agent exhaustive inventory)"),
  "T1_T2_computed": "not computed by T1/T2 (out of scope)",
  "verdict": "SILENT-NEEDS-NEW-WORK",
 },
]

# ------------------------------------------------------------------ closure obligation (T3-G2)
CLOSURE = {
 "name": "the knotless-closure obligation (proposed corpus-side open problem, F-Q2-prime style; offer-class)",
 "requirement_from_print": {
   "Q6prime_language_lock": ("'no fractionally-wound species can close as an asymptotic state at "
     "the vacuum's hbar. Confinement removes the requirement — only line-neutral composites owe "
     "closure, and they close at the emergent level' + '(Language lock F-A14-2: confined "
     "constituents are internal moduli; only asymptotic states owe closure.)' — "
     "corpus3/01-GUM-Omega-Paper-v4.3-ext.md:615 (Thm Q-6-prime), verified this session"),
   "closure_machinery": ("Sec. IV.I Thm T-B1: the spin-clock closure w*j(1-j) = 1, solvable iff "
     "0 < j < 1 — built for KNOT rotors (rotor number j; j = 1/2 the unique admissible elementary "
     "rotor) — corpus3/01:361-367, verified this session"),
   "the_state": ("a knot-free closed T2 tube loop is line-neutral (terminates nothing) and "
     "asymptotic (finite energy = tension x finite length; nothing confines it) — so under the "
     "printed language lock it OWES closure at the vacuum's hbar"),
 },
 "branches": [
   {"id": "(i) emergent-level closure mechanism for a knotless state",
    "content": ("exhibit the loop's hbar-clock: Q-6-prime's blanket sentence says line-neutral "
      "composites 'close at the emergent level', but every exhibited closure mechanism runs "
      "through knot rotors (T-B1's w*j(1-j)=1; the B=2 composite via constituent moduli, "
      "01:831 App F.5) — a knotless loop has no constituents to aggregate; the candidate clock "
      "is the loop's own collective spectrum (T1's rotational/vibrational modes), but no printed "
      "statement extends the closure condition to loop collective coordinates — new construction, "
      "[CJ-new] at every step")},
   {"id": "(ii) principled excusal",
    "content": ("print a principle under which the loop does not owe closure — e.g. only STABLE "
      "asymptotic states owe closure (the census asserts the class decays: WS-D:11), making a "
      "finite-width resonance exempt; no such principle is printed; the corpus's own genre "
      "precedent for a legal evasal exists (the heliknoton 'evading the mass floor legally, "
      "because its stabilizer is the background helix' — 01:831 App F.6) but no closure analog "
      "of it is printed")},
   {"id": "(iii) cannot-close => contradiction with the relic census",
    "content": ("if the loop provably cannot close at the emergent level, Q-6-prime + T-B1 "
      "forbid the asymptotic state — but the WS-D relic census PRESUPPOSES the class exists "
      "('stratum composites (glueball-analogs) -> decay to hadrons/phasons [check]', WS-D:11 — "
      "the no-WIMP theorem D-0 needs the class in its exhaustive inventory), an internal "
      "contradiction between the operative closure sector and the v1 dated record; within-model "
      "structural cost only — nothing here bears on nature")},
 ],
 "print_adjudication": ("NEW-OPEN. Print contains a blanket ASSERTION ('they close at the "
   "emergent level', 01:615) whose exhibited support covers only knot-bearing composites; for a "
   "knotless loop no mechanism, no excusal, and no impossibility proof is printed. Citing the "
   "blanket sentence as discharge would be grade inflation; citing branch (iii) as fired would "
   "be overclaim. The obligation is open in exactly the F-Q2-prime sense: sharply posed, "
   "unadjudicated, owned by the sector."),
 "proposed_obligation_text": ("#8-class NEW OPEN (proposed, offer-class, strike-able): "
   "**F-Q8 (the knotless-closure obligation):** exhibit the emergent-level closure mechanism "
   "(the hbar-clock) for a knot-free closed T2 tube loop, or print the principle that excuses "
   "it (dual-constrained: Q-6-prime's language lock makes the loop closure-owing as an "
   "asymptotic line-neutral state; the WS-D D-0 relic census presupposes the class exists — a "
   "cannot-close verdict contradicts the census's exhaustive inventory). "
   "[F-T9-T3; tier9-glueball/T3/RESULTS.md]"),
}

# ------------------------------------------------------------------ relic census note (T3-G3)
CENSUS_NOTE = {
 "printed": ("WS-D Theorem D-0 census row (v1 dated record only; zero hits in corpus2/corpus3): "
   "'stratum composites (glueball-analogs) -> decay to hadrons/phasons [check]' — existence "
   "presupposed, instability asserted, no rate/channel/mass content (WS-D-Q1-Q2-Relic-Census-v0_1.md:11)"),
 "IM_measured": ("X(2370): Gamma = 170 (+44/-29) MeV, all OBSERVED modes hadronic (KKbar-pi, "
   "pi-pi-eta, pi-pi-eta-prime, KKbar-eta-prime, a0(980)pi0); Gamma/m = 0.072 — an unstable "
   "hadronically-decaying state [IM]"),
 "consistency": ("QUALITATIVE CONSISTENCY AT CENSUS GRADE: the census asserts exactly two "
   "properties (the class exists; it decays to hadrons/phasons) and the [IM] anchor exhibits "
   "an existing state class decaying hadronically. No overclaim either direction: the census "
   "prints no rate, so the 170 MeV number adjudicates nothing about it; the phason channel is "
   "an invisible width no measurement addresses; the census line is dated-record, not operative "
   "text; and the [IM] identification is itself 'dominant component', never 'pure' "
   "(community caveat). The sentence 'the census is confirmed' is unavailable; so is 'the "
   "census is contradicted'. What can be said: the v1 census's qualitative requirements on the "
   "class are not in tension with the [IM] anchor."),
 "gamma_over_m": round(gamma_over_m, 4),
}

# ------------------------------------------------------------------ stakes constraints (T3-G4)
STAKES = {
 "verbatim_constraints": [
   "no stake, kill row, or watch row anywhere names glueballs, closed flux-tube loops, gluonia, or X(2370) (stakes agent, exhaustive sweep; the only hadron-spectrum kill: 'tetraquark-like ground-state baryons (if WS1-H3 returns 3-fold)', corpus3/01:779)",
   "seal q-theta (corpus/WS-N-Hadronic-Pass-P7-v0_1.md:9): 'hadron-spectrum numerology — any matching of specific hadron or nuclear masses to corpus constants — SEALED; forced arrivals run the V.F protocol'",
   "V.F protocol (corpus3/01:437): 'dimensionally secure / structurally plausible / quantitatively unclaimed'",
   "the only legal verdict categories (stakes agent assessment, flagged coordinator-input): consistent -> offered 'banked structural consistency' / 'compliance, uncounted' (credit zero); inconsistent -> 'Q-5's identification takes the hit at inheritance grade' (T-N5); silent -> ABSENT-CERTIFIED with new constructions at [CJ]/[DW, testbed-grade]",
   "unavailable sentences in every case: 'GUM predicted X(2370)', 'X(2370) confirms GUM', 'X(2370) refutes GUM', any grade rise, any campaign-issued stake ('stakes: zero' discipline, WS-N:37)",
   "the event moves no board clock: X(2370) touches no board row, so the 30-day law does not engage; at most a quarterly-memo line, corpus-side (Watch-Mode constitution, v1.0:33-49)",
   "grades never rise by replication (REVISION_CHARTER_v2.md rule 4); offers, not adoptions (rule 5); 'no external sentence may exceed' App. J (01:853)",
 ],
 "compliance": {
   "stake_issued": False, "board_clock_moved": False,
   "forbidden_sentences_present": False, "grade_risen": False,
   "all_numbers_VF_graded": True,
 },
}

# ------------------------------------------------------------------ gates
GATES = [
 {"id": "T3-G1", "requirement": "scorecard complete, every corpus cell carrying file:line evidence, every verdict one of the four pre-registered categories",
  "measured": None, "verdict": None},
 {"id": "T3-G2", "requirement": "closure obligation stated with the three branches and a print-answerable adjudication; proposed obligation text drafted (offer-class)",
  "measured": "stated (requirement from 01:615 + 01:361-367, both verified by session reads); three branches enumerated; adjudication NEW-OPEN; obligation text drafted in offer language",
  "verdict": "PASS"},
 {"id": "T3-G3", "requirement": "relic-census consistency note at census grade (no overclaim in either direction)",
  "measured": "note delivered: qualitative consistency at census grade; both overclaim directions explicitly barred; Gamma/m = %.3f arithmetic only" % gamma_over_m,
  "verdict": "PASS"},
 {"id": "T3-G4", "requirement": "stakes-machinery constraints honored verbatim: no stake, no board clock, unavailable sentences absent, V.F on all numbers",
  "measured": "constraints quoted verbatim from stakes-agent record; compliance block all-clear; forbidden sentences absent from all T3 deliverables",
  "verdict": "PASS"},
]

def validate():
    ok = True
    n_cited = 0
    for row in SCORECARD:
        if row["verdict"] not in CATEGORIES:
            print("BAD VERDICT:", row["criterion"], row["verdict"]); ok = False
        # file:line evidence check: require a ':NNN' pattern in the corpus cell
        import re
        if re.search(r"\.md:\d+", row["corpus_as_printed"]):
            n_cited += 1
        else:
            print("MISSING file:line in corpus cell:", row["criterion"]); ok = False
    pend = sum(1 for row in SCORECARD if "PENDING-T1" in str(row["T1_T2_computed"]))
    g1 = GATES[0]
    g1["measured"] = ("7/7 rows; %d/7 corpus cells carry file:line; verdicts all in the four "
                      "pre-registered categories; %d T1 cells pending" % (n_cited, pend))
    g1["verdict"] = "PASS" if (ok and pend == 0) else ("PARTIAL" if ok else "FAIL")
    return ok

def main():
    ok = validate()
    out = {
        "workstream": "T3", "finding": "F-T9-T3",
        "register": "within-model; nothing bears on nature; BESIII/lattice = [IM] anchors; seal q-theta V.F on all numbers; numerology hazard (c-frak 2.37 vs m_X 2.37 GeV) pre-empted, cited nowhere as structure",
        "T1_present_at_assembly": T1_PRESENT, "T2_present_at_assembly": T2_PRESENT,
        "IM_anchors": IM,
        "scorecard": SCORECARD,
        "closure_obligation": CLOSURE,
        "relic_census_note": CENSUS_NOTE,
        "stakes_constraints": STAKES,
        "gates": GATES,
        "forbidden_sentence_scan": None,
        "validation_ok": ok,
    }
    # forbidden-sentence self-scan over this very JSON
    blob = json.dumps(out)
    forbidden = ["GUM predicted X(2370)", "confirms GUM", "refutes GUM",
                 "X(2370) confirms", "X(2370) refutes"]
    hits = [s for s in forbidden if s in blob.replace("'", "'")]
    # the list above appears inside the *quotation of the ban itself*; scan outside quotes:
    # operational rule: the strings may appear ONLY inside the stakes 'unavailable sentences' row.
    blob_wo_ban = blob.replace(STAKES["verbatim_constraints"][4], "")
    hits = [s for s in forbidden if s in blob_wo_ban]
    out["forbidden_sentence_scan"] = {"hits_outside_ban_quotation": hits, "clean": not hits}
    if hits:
        print("FORBIDDEN SENTENCES FOUND:", hits); sys.exit(2)
    with open(os.path.join(HERE, "t3_scorecard.json"), "w") as fh:
        json.dump(out, fh, indent=1, ensure_ascii=False)
    print("T1 present:", T1_PRESENT, "| T2 present:", T2_PRESENT)
    for g in GATES:
        print(g["id"], "->", g["verdict"], "|", g["measured"])
    print("wrote t3_scorecard.json; validation_ok =", ok)

if __name__ == "__main__":
    main()
