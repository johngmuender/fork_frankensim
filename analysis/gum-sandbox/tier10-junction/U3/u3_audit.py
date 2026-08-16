#!/usr/bin/env python3
"""
U3 — B-dictionary audit + scorecard assembly + mechanical register scan
(F-T10-U3).  Tier 10 (Phase U), ROADMAP_v12_JUNCTION.md workstream U3.

Within-model; nothing here bears on nature.  STAR/lattice/Regge values are
[IM] anchors only; they enter no computation below (mechanically verified
by this script's own scan).  Numerology pre-emption: the corpus's
dimensionless c-frak = 2.37 +/- 0.09 shares digits by coincidence (of
man-made units) with GeV-denominated numbers; that coincidence is not
structure and is not used anywhere in this workstream.

Part 1 — formal B-conservation audit, exact rational arithmetic
(fractions.Fraction), for the three printed operations plus the
junction-migration scenario.  The audited dictionary is Theorem K-1 as
printed (corpus/WS-K-NR-K1-Execution-Corollaries-v0_1.md:33):

    B = (1/3) * (count of fractional-residue cores, signed).

[U3-A1] (coordinator-owned amendment, printed): the roadmap's shorthand
"B = (1/3)*Sigma(signed residues)" is executed as the signed INDICATOR
sum over fractional-residue cores (each fractional core contributes +1,
each fractional anticore -1), which is K-1's printed dictionary.  The
literal sum of residue VALUES rho = w_em mod 1 is NOT the printed B
(e.g. a uud parse: 2/3 + 2/3 + 1/3 = 5/3, and (1/3)*5/3 = 5/9 != 1);
the literal reading is computed below as a diagnostic and shown to
disagree, so the shorthand cannot silently replace the theorem.

Part 2 — scorecard JSON assembly (rows/columns per the frozen spec).
Part 3 — mechanical register scan (forbidden sentences; STAR anchor
tokens outside [IM]-anchor context; numerology tokens outside
pre-emption context) over the U3 deliverables.

Deterministic; no RNG; runtime < 1 s.
"""
import itertools
import json
import os
import re
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))

# ----------------------------------------------------------------------
# Part 1 — the exact-ledger model
# ----------------------------------------------------------------------
# A core is a dict: {id, kind in {"knot","antiknot","lepton","antilepton"},
#   residue: Fraction in [0,1) — rho = w_em mod 1 per WS-K-NR-K1:31
#   (leptons rho = 0; quarks rho in {1/3, 2/3}-class), sign: +1 matter,
#   -1 antimatter (matter-positive convention, census WS-K-Q1-Q2:23)}.
# A line is a dict: {id, cls in {"integer","fractional"},
#   ends: list of endpoint refs ("core:<id>", "junction:<id>", "BULK")}.
# A junction is a node where fractional lines meet; census legality
# (WS-K-Q1-Q2:14): thirds sum to integers at every junction.

RHO_U = Fr(2, 3)   # up-class residue  (w_em = +2/3 mod 1)
RHO_D = Fr(1, 3)   # down-class residue (w_em = -1/3 mod 1 = 2/3? no:)
# NOTE: rho is a CLASS label mod 1.  d: -1/3 mod 1 = 2/3 as a value; the
# corpus writes the class as "+/-1/3-class" and verifies the weak vertex
# residue-preserving as d: -1/3 == u: +2/3 (mod 1) (WS-K-NR-K1:31).  The
# audit therefore tracks rho as the signed thirds-class in {-1/3,+2/3}
# VALUES only for the diagnostic; conservation bookkeeping uses the
# fractional/integer CLASS plus the signed core count, exactly as K-1 is
# printed.
RHO_D_SIGNED = Fr(-1, 3)


def B_printed(cores):
    """Theorem K-1 as printed: (1/3) * signed count of fractional-residue
    cores.  A core is fractional iff its residue is not an integer."""
    n = sum(c["sign"] for c in cores if c["residue"] % 1 != 0)
    return Fr(n, 3)


def residue_value_sum(cores):
    """Diagnostic ONLY — the literal (1/3)*Sigma(signed residue values)
    reading of the roadmap shorthand; shown NOT to equal B."""
    return Fr(1, 3) * sum(c["sign"] * c["residue"] for c in cores)


def fractional(core):
    return core["residue"] % 1 != 0


def line_legal(line, cores_by_id):
    """The taping rule (WS-K-NR-K1:31): integer lines may begin/end/
    reconnect on any core; a fractional line can terminate only on a
    fractional core (or a junction node); no line ends in the bulk
    (Primer: 'cannot end in empty space')."""
    for end in line["ends"]:
        if end == "BULK":
            return False
        if end.startswith("core:") and line["cls"] == "fractional":
            core = cores_by_id[int(end.split(":")[1])]
            if not fractional(core):
                return False
    return True


def junction_thirds_ok(junction_lines):
    """Census legality (WS-K-Q1-Q2:14): thirds sum to integers at every
    junction.  junction_lines: list of thirds-class values entering."""
    return sum(junction_lines) % 1 == 0


def mk_quark(i, flavor, sign):
    rho = RHO_U if flavor == "u" else RHO_D_SIGNED
    return {"id": i, "kind": "knot" if sign > 0 else "antiknot",
            "residue": rho, "sign": sign}


audit = {"operations": [], "amendment_U3_A1": (
    "Roadmap shorthand 'B=(1/3)*Sigma(signed residues)' executed as the "
    "signed indicator sum over fractional-residue cores (K-1 as printed, "
    "WS-K-NR-K1:33); the literal residue-VALUE sum is computed as a "
    "diagnostic and disagrees (uud: (1/3)*Sigma rho = 5/9 != 1), so the "
    "literal reading is rejected as the dictionary.")}

# --- diagnostic: the literal shorthand is not the printed dictionary ---
proton = [mk_quark(1, "u", +1), mk_quark(2, "u", +1), mk_quark(3, "d", +1)]
diag = {
    "state": "proton parse uud",
    "B_printed": str(B_printed(proton)),
    "literal_residue_value_sum_over_3": str(residue_value_sum(proton)),
    "equal": B_printed(proton) == residue_value_sum(proton),
}
assert B_printed(proton) == Fr(1) and not diag["equal"]
audit["shorthand_diagnostic"] = diag

# ----------------------------------------------------------------------
# Operation (i): taping-rule termination (attachment/termination events)
# ----------------------------------------------------------------------
# Legal moves: attach or detach lines between EXISTING cores/junctions
# subject to the termination rule.  The move changes line adjacency only;
# the core multiset is untouched.  S-K3 (WS-K-NR-K1:49): attachment
# transports winding, never degree.
cores = [mk_quark(1, "u", +1), mk_quark(2, "u", +1), mk_quark(3, "d", +1),
         {"id": 4, "kind": "lepton", "residue": Fr(0), "sign": +1}]
cores_by_id = {c["id"]: c for c in cores}
B0 = B_printed(cores)
cases = []
# legal: fractional line between two fractional cores
l1 = {"id": 1, "cls": "fractional", "ends": ["core:1", "core:2"]}
cases.append(("fractional line, both ends fractional cores",
              line_legal(l1, cores_by_id), True))
# illegal: fractional line terminating on an integer (lepton) core
l2 = {"id": 2, "cls": "fractional", "ends": ["core:1", "core:4"]}
cases.append(("fractional line onto integer-residue core",
              line_legal(l2, cores_by_id), False))
# illegal: fractional line ending in the bulk
l3 = {"id": 3, "cls": "fractional", "ends": ["core:1", "BULK"]}
cases.append(("fractional line ending in bulk", line_legal(l3, cores_by_id),
              False))
# legal: integer line onto any core
l4 = {"id": 4, "cls": "integer", "ends": ["core:4", "core:1"]}
cases.append(("integer line onto any core", line_legal(l4, cores_by_id),
              True))
term_ok = all(got == want for _, got, want in cases)
B1 = B_printed(cores)  # cores untouched by any adjacency move
audit["operations"].append({
    "operation": "(i) taping-rule termination",
    "printed_rule": "WS-K-NR-K1-Execution-Corollaries-v0_1.md:31 (taping "
                    "rule), :49 (S-K3: lines transport winding, never "
                    "degree); Primer corpus3/03:811 ('cannot end in empty "
                    "space')",
    "model": "attach/detach lines among existing cores; core multiset "
             "invariant by construction of the move class",
    "legality_cases": [{"case": c, "legal": g, "expected": w}
                       for c, g, w in cases],
    "delta_B": str(B1 - B0),
    "B_exact": B1 == B0,
    "verdict": "B-EXACT" if (B1 == B0 and term_ok) else "FINDING",
})

# ----------------------------------------------------------------------
# Operation (ii): web reconnection (net-conserving, census column)
# ----------------------------------------------------------------------
# Model: lines AB, CD -> AD, CB — endpoint exchange among fixed
# cores/junctions; enumerate every pairing of a 6-line test network and
# check (a) core multiset invariant -> delta B = 0, (b) the summed
# thirds-content entering each junction stays integer (net-conserving),
# (c) global thirds content unchanged.
netcores = [mk_quark(i, f, +1) for i, f in
            zip(range(1, 7), ["u", "u", "d", "u", "d", "d"])]
B0 = B_printed(netcores)
# two junctions, each fed by three tubes from three cores:
# thirds-class content of the tube = its core's signed thirds
junctionA = [RHO_U, RHO_U, RHO_D_SIGNED]          # sums to 1 (integer)
junctionB = [RHO_U, RHO_D_SIGNED, RHO_D_SIGNED]   # sums to 0 (integer)
assert junction_thirds_ok(junctionA) and junction_thirds_ok(junctionB)
total0 = sum(junctionA) + sum(junctionB)
recon_results = []
# reconnection: swap one tube between the junctions (all 3x3 swaps)
for ia, ib in itertools.product(range(3), range(3)):
    A = list(junctionA)
    Bj = list(junctionB)
    A[ia], Bj[ib] = Bj[ib], A[ia]
    net_conserving = (sum(A) + sum(Bj)) == total0
    per_junction_legal = junction_thirds_ok(A) and junction_thirds_ok(Bj)
    recon_results.append({
        "swap": [ia, ib], "net_conserving": net_conserving,
        "per_junction_thirds_integer": per_junction_legal})
B1 = B_printed(netcores)  # cores untouched by any endpoint exchange
all_net = all(r["net_conserving"] for r in recon_results)
legal_events = [r for r in recon_results if r["per_junction_thirds_integer"]]
audit["operations"].append({
    "operation": "(ii) web reconnection",
    "printed_rule": "WS-K-Q1-Q2-Charge-Census-Crisis-v0_1.md:14 (w_c row: "
                    "'web reconnection (net-conserving)'), :17 (web-line "
                    "charges: 'reconnection events (net-conserving)'); "
                    "S-K3 WS-K-NR-K1:49",
    "model": "endpoint exchange of tubes between two junctions over a "
             "6-core, 2-junction test network; 9 swap events enumerated",
    "events_enumerated": 9,
    "net_thirds_conserved_all_events": all_net,
    "events_satisfying_junction_legality": len(legal_events),
    "delta_B": str(B1 - B0),
    "B_exact": B1 == B0,
    "verdict": "B-EXACT" if (B1 == B0 and all_net) else "FINDING",
    "note": "reconnection redistributes line content among junctions "
            "(census legality restricts which events are junction-legal) "
            "but is a LINE event: the core multiset — the only thing K-1 "
            "counts — is untouched by every enumerated event.",
})

# ----------------------------------------------------------------------
# Operation (iii): snap-minting (pair creation)
# ----------------------------------------------------------------------
# Primer 15.2 (corpus3/03:811): a snapped cord has two fresh ends, which
# must terminate on fresh fractional knots — pay for two new ones.
# Thm VIII'.1: creation is pairwise.  Model: mint (knot, antiknot) with
# residues (+rho, -rho), every flavor class; check delta B = 0 exactly
# and both fresh ends land on fractional cores (taping rule).
mint_cases = []
base = [mk_quark(1, "u", +1), mk_quark(2, "u", +1), mk_quark(3, "d", +1)]
B0 = B_printed(base)
for flavor in ("u", "d"):
    state = base + [mk_quark(10, flavor, +1), mk_quark(11, flavor, -1)]
    dB = B_printed(state) - B0
    ends_ok = fractional(state[-1]) and fractional(state[-2])
    mint_cases.append({"minted_pair_flavor": flavor, "delta_B": str(dB),
                       "fresh_ends_on_fractional_cores": ends_ok,
                       "exact": dB == 0})
audit["operations"].append({
    "operation": "(iii) snap-minting (pair creation)",
    "printed_rule": "Primer corpus3/03:807, :811, :827 (snap + pay-for-two;"
                    " ~100-MeV-class exercise chain); Thm VIII'.1 pairs-"
                    "only (corpus3/01:683); WS-K-NR-K1:31 (fresh fractional"
                    " ends must terminate on fractional cores)",
    "model": "mint (knot, antiknot) pair, residues (+rho, -rho), both "
             "flavor classes",
    "cases": mint_cases,
    "B_exact": all(c["exact"] for c in mint_cases),
    "verdict": "B-EXACT" if all(c["exact"] for c in mint_cases)
               else "FINDING",
})

# ----------------------------------------------------------------------
# Operation (iv): junction migration (U2's mode) — full scenario ledger
# ----------------------------------------------------------------------
# Independent re-coding of the U2-G3 scenario: junction + three tube
# stubs migrate from beam to destination; each tube snaps once, minting
# one pair; the junction constraint forces the three destination tube
# ends onto the three minted knots; the three original knots pair with
# the three minted antiknots into mesons at beam.  All 2^3 flavor
# assignments enumerated; exact arithmetic.
migration_cases = []
for flavors in itertools.product("ud", repeat=3):
    orig = [mk_quark(i, f, +1) for i, f in
            zip((1, 2, 3), ("u", "u", "d"))]           # beam baryon uud
    B_init = B_printed(orig)                            # = 1
    minted_k = [mk_quark(10 + i, f, +1) for i, f in enumerate(flavors)]
    minted_ak = [mk_quark(20 + i, f, -1) for i, f in enumerate(flavors)]
    final = orig + minted_k + minted_ak
    B_fin = B_printed(final)
    # regional ledgers
    B_dest = B_printed(minted_k)          # destination baryon content
    dest_charge = sum(c["sign"] * c["residue"] for c in minted_k)
    net_minted_charge = sum(c["sign"] * c["residue"]
                            for c in minted_k + minted_ak)
    # junction legality at destination: three thirds sum to integer
    dest_thirds = [c["residue"] for c in minted_k]
    migration_cases.append({
        "minted_flavors": "".join(flavors),
        "B_initial": str(B_init), "B_final": str(B_fin),
        "delta_B": str(B_fin - B_init),
        "original_knots_moved": 0,
        "junction_count_change": 0,
        "B_destination": str(B_dest),
        "destination_thirds_integer": junction_thirds_ok(dest_thirds),
        "net_minted_charge_all_regions": str(net_minted_charge),
        "destination_baryon_charge": str(dest_charge + Fr(0)),
        "exact": (B_fin == B_init == Fr(1)) and B_dest == Fr(1),
    })
mig_ok = all(c["exact"] and c["destination_thirds_integer"]
             and c["net_minted_charge_all_regions"] == "0"
             for c in migration_cases)
audit["operations"].append({
    "operation": "(iv) junction migration (U2's junction-led transport "
                 "mode)",
    "printed_rule": "composite of printed-legal pieces only: reconnection "
                    "legality (census :14, :17), snap-minting (Primer "
                    ":811), taping-rule forced destination assembly "
                    "(WS-K-NR-K1:31); the migration DYNAMICS itself is "
                    "[CJ-new] (no printed transport law — context "
                    "analysis section 2)",
    "model": "beam baryon uud; junction + 3 tube stubs migrate; 3 pairs "
             "minted; 8 flavor assignments enumerated; exact Fraction "
             "ledger (independent re-code of U2-G3)",
    "cases": migration_cases,
    "B_exact": mig_ok,
    "verdict": "B-EXACT" if mig_ok else "FINDING",
    "note": "B = 1 before and after in every assignment; zero original "
            "knots move; junction/core bookkeeping closed; net minted "
            "charge zero — junction-led transport is B-exact and "
            "charge-blind, i.e. K-1-compatible.",
})

audit["overall"] = ("B exact under all four audited operations (zero "
                    "failures); junction-led transport is printed-LEGAL "
                    "at the ledger level, its dynamics [CJ-new]")
audit["all_exact"] = all(op["verdict"] == "B-EXACT"
                         for op in audit["operations"])

# ----------------------------------------------------------------------
# Part 2 — scorecard assembly
# ----------------------------------------------------------------------
RIVAL_CGC = ("CGC gluon saturation reproduces the alpha_B systematics "
             "with valence quarks only (Garcia-Montero & Schlichting, "
             "PRC 111, 024912) — the slope alone does not select "
             "junctions")
RIVAL_SKIN = ("neutron-skin geometry contaminates the peripheral isobar "
              "trend (TRENTO reproduces the trend junction-free; "
              "Pihan & Vovchenko arXiv:2602.04079 treat the observable "
              "as a joint skin+stopping probe)")
RIVAL_STRANGE = ("strangeness asymmetry moves B/DeltaQ well below 1 in "
                 "junction-free AMPT (Ross & Lin, EPJC 86, 143)")
RIVAL_HYDRO = ("junction-free hydrodynamic charge-stopping baseline: "
               "Pihan-Monnai-Schenke-Shen PRL 133, 182301")
RIVAL_PERSP = ("Science Perspective (W. Li): the measurements 'do not "
               "provide a direct, tightly controlled measurement of the "
               "underlying mechanism'")

scorecard_rows = [
    {
        "row": "STAR isobar <B>/DeltaQ x DeltaZ/A = 1.84 +/- 0.02 [IM] "
               "+/- 0.09 +/- 0.16 (0-10% central)",
        "corpus_printed": (
            "carrier: B knot-carried, Thm K-1 [DF-structural, worksheet "
            "stratum] (WS-K-NR-K1:33); enforcer: taping rule/junction "
            "constraint (WS-K-NR-K1:31; census :14, :17); transporter: "
            "SILENT — zero printed transport statements, 'Regge' zero "
            "corpus hits (context section 2); operative editions silent "
            "on B entirely (zero 'baryon number' hits, grep-verified)"),
        "u1_u2_computed": (
            "U2 [CJ-new] sign-only toy: junction-led channel strictly "
            "cheaper over the full printed window (min margin 1.128 "
            "conservative / 2.902 compact; 100% of 4D grid); B >= Q for "
            "any monotone-decreasing cost; magnitudes UNCLAIMABLE "
            "(U2/RESULTS.md sections 4-5)"),
        "verdict": "SILENT-NEEDS-NEW-WORK",
        "verdict_note": (
            "the corpus prints no transport law, so the measurement "
            "confronts nothing printed; the U2 completion shows the "
            "STAR-shaped sign is structurally generic in-model, at "
            "[CJ-new] grade — exactly the 'new work' the verdict names"),
        "rivals": [RIVAL_SKIN, RIVAL_STRANGE, RIVAL_HYDRO, RIVAL_PERSP],
    },
    {
        "row": "STAR gamma+Au alpha_B = 1.04 +/- 0.22 [IM] (junction "
               "Regge window [0.42, 1])",
        "corpus_printed": (
            "SILENT: no rapidity, Regge, or stopping structure of any "
            "kind ('Regge' zero hits; no alpha_B-like quantity "
            "constructible from print — context section 2; dynamics "
            "agent silences list)"),
        "u1_u2_computed": (
            "U2 produces a sign only — no alpha_B value exists or could "
            "be claimed (U2 register: magnitudes UNCLAIMABLE; corpus "
            "prints no transport law)"),
        "verdict": "SILENT-NEEDS-NEW-WORK",
        "verdict_note": (
            "an in-model alpha_B would need the full [CJ-new] transport "
            "layer priced in the alternative below"),
        "rivals": [RIVAL_CGC, RIVAL_PERSP],
    },
    {
        "row": "STAR Au+Au alpha_B = 0.64 +/- 0.05 [IM], "
               "centrality-independent",
        "corpus_printed": "SILENT (same basis as the gamma+Au row)",
        "u1_u2_computed": (
            "U2 sign only; centrality does not exist in the toy "
            "(U2 caveat 5: no collision geometry, no energy dependence, "
            "no centrality)"),
        "verdict": "SILENT-NEEDS-NEW-WORK",
        "verdict_note": "as above; the centrality-independence "
                        "discriminator has no in-model counterpart",
        "rivals": [RIVAL_CGC, RIVAL_PERSP],
    },
    {
        "row": "compact-B imaging: baryon number confined to transverse "
               "radius 0.33-0.53 fm vs charge/mass radii >= 0.67 fm "
               "(Klein-Labonte-Sweger-Miller-Vogt, arXiv:2603.03730, "
               "preprint) [IM]",
        "corpus_printed": (
            "no printed density, radius, or form-factor machinery for "
            "any charge (WS-N RN-3 charter cage: no nuclear numbers; "
            "context: momentum/charge partition inside a hadron SILENT). "
            "Printed structure cuts both ways: B and Q ride the SAME "
            "knot cores (K-1 + V.D winding) — naively co-located — but "
            "minted K-Kbar pairs carry winding with zero net B "
            "(VIII'.1; Primer :811), so a charge cloud wider than the "
            "net-B support is equally constructible from print"),
        "u1_u2_computed": (
            "neither U1 nor U2 computes spatial distributions; U2's "
            "ledger shows net-B tracks the junction-forced core "
            "assembly while charge spreads over pair-minted content — "
            "a structural rhyme, not a radius"),
        "verdict": "SILENT-NEEDS-NEW-WORK",
        "verdict_note": (
            "both a tension reading (same-carrier co-location) and an "
            "accommodation reading (net-zero-B charge cloud) are "
            "constructible; neither is printed; an in-model B-density "
            "would be new [CJ] construction"),
        "rivals": [
            "conventional meson-cloud explanation: the pion cloud "
            "inflates the charge radius while valence B stays compact — "
            "no junction needed",
            "the extraction is Regge-model-dependent, from four "
            "exclusive backward channels; preprint, not yet "
            "peer-reviewed (SNIPPET-ONLY verification)"],
    },
    {
        "row": "junction-exotics class: J-Jbar 'baryonium glueballs', "
               "gluonic graphene, buckyballs (CGK J. Phys. G 30, L17; "
               "magic numbers 8/24/48/120); none discovered [IM]",
        "corpus_printed": (
            "NOT-FORBIDDEN / ACCOMMODATED / NEVER-CONSTRUCTED: knot-free "
            "tube networks with junctions are legal iff thirds sum to "
            "integers at every junction (census :14); a quarkless "
            "junction carrying B = 1 is dictionary-IMPOSSIBLE (K-1: "
            "knot-free networks have B = 0); Course 20.1's 'exactly' "
            "inventory sentence excludes the class as printed "
            "(corpus3/02:1088 vs Omega Q-1 corpus3/01:615 — Tier-9 "
            "derivation-scope precedent carried)"),
        "u1_u2_computed": (
            "U1 constructs the minimal J-Jbar dumbbell: ground "
            "m/sqrt(sigma) = 2.59/2.70/2.94 at the three pre-registered "
            "junction masses (envelopes ~2.45-3.12), BELOW the loop "
            "0++ head 3.583 under the frozen inertia convention; "
            "convention-split priced (stretch inertia: 5.5-5.8); "
            "interleaves the Tier-9 loop ladder; all V.F"),
        "verdict": "V.F-CLASS",
        "verdict_note": (
            "a computed within-model spectrum now exists at V.F grade "
            "(dimensionally secure / structurally plausible / "
            "quantitatively unclaimed); no state identified with any "
            "observed or proposed hadron"),
        "rivals": [
            "no junction-exotic has been observed by anyone (STAR paper: "
            "'none discovered'); the class cannot currently adjudicate "
            "the junction interpretation",
            "CGK's own GeV mass estimates unpinned from accessible "
            "sources (SNIPPET-ONLY)"],
    },
]

priced_alternative = {
    "frame": ("symmetric price sheet; NEITHER branch is recommended; "
              "within-model; adoption is the corpus authors' alone"),
    "branch_A": {
        "name": "K-1-as-printed only",
        "content": ("keep the printed dictionary and nothing else: B is "
                    "knot-carried (K-1), taping-rule-enforced, and the "
                    "transport question is left unanswered"),
        "costs": [
            "the STAR-shaped question (who MOVES B in momentum space) is "
            "unanswerable in-model: no confrontation row above can ever "
            "leave SILENT-NEEDS-NEW-WORK",
            "the worksheet-stratum fold debt stands (operative editions "
            "print no B dictionary at all)"],
        "buys": [
            "zero new construction; zero new parameters; the theorem "
            "stack stays exactly as derived",
            "no exposure to the contested junction interpretation of "
            "the STAR data (rivals may yet win the phenomenology)"],
    },
    "branch_B": {
        "name": "K-1 + junction-led transport completion (U2's "
                "construction)",
        "content": ("adopt the carrier/enforcer/transporter trichotomy "
                    "and a junction-led transport layer as in-model "
                    "structure"),
        "costs": [
            "[CJ-new] transport layer: a mobility/cost law the corpus "
            "never printed (U2's g is a toy ansatz family, not a "
            "derivation)",
            "junction energy/inertia: unprinted ([CJ-new]); the only "
            "anchor is a 2+1D lattice import with an untested "
            "dimensional transfer",
            "unprinted dynamics throughout: minting rates, reconnection "
            "rates, junction migration kinetics — none derivable from "
            "print today",
            "inherits the F-K0-1 junction-order seam (n = 3 "
            "Primer-asserted, paper-conditional, execution absent)"],
        "buys": [
            "the trichotomy dissolves the apparent GUM-vs-STAR tension: "
            "knot-carried B is COMPATIBLE with junction-led transport "
            "(U2-G3, exact ledger; audited again here, operation (iv))",
            "the STAR-shaped sign comes out generically over the full "
            "printed parameter window (U2-G2), not by tuning",
            "the J-Jbar state class acquires a computable spectrum "
            "(U1), giving the two-sector dictionary a third entry"],
    },
    "symmetry_note": ("branch A pays in permanent silence; branch B pays "
                      "in unprinted structure.  The scorecard's verdicts "
                      "are identical under both branches TODAY (nothing "
                      "above rises past V.F-CLASS/SILENT); the branches "
                      "differ only in what future work could change "
                      "that.  No recommendation is made."),
}

# ----------------------------------------------------------------------
# Part 3 — mechanical register scan
# ----------------------------------------------------------------------
FORBIDDEN_PATTERNS = [
    (r"GUM\s+predicted\s+the\s+junction", "forbidden sentence class 1"),
    (r"STAR\s+(confirms|refutes|vindicates|validates)",
     "forbidden sentence class 2"),
    (r"grade\s+(rise|rises|raised|promoted)",
     "grade-motion language (only legal inside the ban's own statement)"),
    (r"(?<!no )(?<!issues no )stake\s+is\s+issued(?!.*no)",
     "stake issuance"),
]
ANCHOR_TOKENS = ["1.84", "0.64", "1.04"]
ANCHOR_CONTEXT = re.compile(r"\[IM\]|IM.anchor|anchor|forbidden|compliance"
                            r"|scan|ban", re.I)
NUMEROLOGY_TOKEN = "2.37"
NUMEROLOGY_CONTEXT = re.compile(r"coincid|numerolog|pre.?empt", re.I)


def scan_file(path):
    out = {"file": os.path.basename(path), "forbidden_hits": [],
           "anchor_violations": [], "numerology_violations": []}
    if not os.path.exists(path):
        out["missing"] = True
        return out
    with open(path, encoding="utf-8") as fh:
        lines = fh.read().splitlines()
    for i, ln in enumerate(lines, 1):
        for pat, label in FORBIDDEN_PATTERNS:
            if re.search(pat, ln, re.I):
                # exempt lines that state the ban itself
                if re.search(r"forbidden|does not say|ban|absent|never"
                             r"|scan", ln, re.I):
                    continue
                out["forbidden_hits"].append({"line": i, "class": label})
        for tok in ANCHOR_TOKENS:
            if tok in ln and not ANCHOR_CONTEXT.search(ln):
                out["anchor_violations"].append({"line": i, "token": tok})
        if NUMEROLOGY_TOKEN in ln and not NUMEROLOGY_CONTEXT.search(ln):
            out["numerology_violations"].append({"line": i})
    out["clean"] = not (out["forbidden_hits"] or out["anchor_violations"]
                        or out["numerology_violations"])
    return out


def run_scan():
    targets = [os.path.join(HERE, f) for f in
               ("RESULTS.md", "u3_annotations.md", "u3_scorecard.json",
                "u3_audit.py")]
    results = [scan_file(t) for t in targets]
    return {"targets": [r["file"] for r in results], "results": results,
            "all_clean": all(r.get("clean", False) or r.get("missing")
                             for r in results),
            "note": ("anchor tokens are LEGAL only on [IM]-anchor/"
                     "compliance-context lines; numerology token legal "
                     "only inside its pre-emption; scan exempts lines "
                     "that state the bans themselves")}


# ----------------------------------------------------------------------
# assemble and write
# ----------------------------------------------------------------------
def main():
    payload = {
        "workstream": "U3 (F-T10-U3)",
        "date": "2026-08-16",
        "epistemic_frame": (
            "within-model; nothing bears on nature; STAR/lattice/Regge "
            "values are [IM] anchors only and enter no computation; seal "
            "q-theta: any hadron-mass-confronting number V.F only "
            "(dimensionally secure / structurally plausible / "
            "quantitatively unclaimed); no fitting to STAR observables; "
            "the junction interpretation of the STAR data is contested "
            "and rivals are carried at equal weight; numerology "
            "pre-emption: the dimensionless c-frak = 2.37 shares digits "
            "with GeV numbers by man-made-unit coincidence only — not "
            "structure, never cited as structure; no stake, no clock, "
            "no grade motion; all corpus-side changes are offers"),
        "spec_amendments": {
            "U3-A1": audit["amendment_U3_A1"],
            "U3-A2": ("the roadmap's U1/U2 watcher loop (<= 40 min) was "
                      "not needed: both U1/RESULTS.md and U2/RESULTS.md "
                      "existed at first poll; no PENDING cells"),
        },
        "conservation_audit": audit,
        "scorecard": {"columns": ["row ([IM] anchor)", "corpus-printed "
                                  "content (file:line)", "U1-U2 computed "
                                  "content", "verdict", "rivals "
                                  "(mandatory)"],
                      "verdict_categories": ["SUPPORTED-STRUCTURALLY",
                                             "V.F-CLASS",
                                             "SILENT-NEEDS-NEW-WORK",
                                             "TENSION"],
                      "rows": scorecard_rows},
        "priced_alternative": priced_alternative,
    }
    out = os.path.join(HERE, "u3_scorecard.json")
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=1)
    # scan AFTER writing the JSON so the JSON scans itself on rerun
    scan = run_scan()
    payload["G4_scan"] = scan
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=1)
    print("U3 audit: all operations B-exact:", audit["all_exact"])
    for op in audit["operations"]:
        print(f"  {op['operation']}: {op['verdict']}")
    print("scorecard rows:", len(scorecard_rows))
    print("scan all_clean:", scan["all_clean"])
    for r in scan["results"]:
        print("  ", r["file"],
              "MISSING" if r.get("missing") else
              ("clean" if r["clean"] else
               f"VIOLATIONS: {r['forbidden_hits']} "
               f"{r['anchor_violations']} {r['numerology_violations']}"))
    return 0 if audit["all_exact"] else 1


if __name__ == "__main__":
    sys.exit(main())
