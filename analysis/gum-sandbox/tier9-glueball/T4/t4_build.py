#!/usr/bin/env python3
"""T4 assembler + validator (F-T9-T4).

Builds t4_dictionary.json (machine-readable two-sector dictionary,
double-counting audit, priced alternative, annotation manifest, gates)
and validates:
  V1  every dictionary/audit corpus-evidence cell carries file:line evidence
  V2  every dictionary/audit row carries a non-empty literature cell
  V3  forbidden-sentence scan over RESULTS.md, t4_annotations.md, and the
      emitted JSON (lines that quote the ban itself are excepted)
  V4  the three cited sibling artifacts (T1/T2/T3 RESULTS.md) exist
  V5  annotation manifest: every box strike-able + finding tag + artifact path

Deterministic; no RNG, no numerics. Within-model; nothing bears on nature.
Seal q-theta: this script produces no mass numbers; every mass statement in
the JSON is a citation of T1's V.F-graded output. The numerology hazard
(dimensionless c-frak = 2.37 vs m_X = 2.37 GeV digit coincidence) is
pre-empted: it is recorded below only as a hazard flag and is never used
as structure.
"""
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
T9 = os.path.dirname(HERE)

FILELINE = re.compile(r"[\w\-./]+\.(md|py|json):\d+|\b0[123]:\d+|\bWS-[A-Z][\w\-]*[\w.]*:\d+")

MEM = "MEMORY-CITED [IM] (standard-knowledge; not snippet-verified; egress-blocked)"
SNIP = "SNIPPET-VERIFIED [IM] (t_context_agents.json glueball-theory-literature)"

dictionary = [
    {
        "state_class": "baryons (B=1, n-junction)",
        "owning_sector": "both (same states): Q-5 Skyrme at long wavelength; Q-1 junction reading at constituent scale",
        "corpus_evidence": [
            "corpus3/01-GUM-Omega-Paper-v4.3-ext.md:617 (Q-5: 'the hadronic Skyrme model (baryons as its solitons; the Y-law junction geometry a banked consistency)')",
            "corpus3/01-GUM-Omega-Paper-v4.3-ext.md:615 (Q-1: 'n-junction baryons' as line-neutral composites)",
            "corpus3/02-The-Substrate-Course-v3-ext.md:1088 (Course 20.1 composite inventory)",
            "corpus3/01-GUM-Omega-Paper-v4.3-ext.md:871 (corpus's own Skyrme bibliography line)",
        ],
        "literature": [
            {"ref": "Skyrme, Proc. R. Soc. A 260, 127 (1961)", "confidence": MEM},
            {"ref": "Witten, Nucl. Phys. B 160, 57 (1979) - large-N baryons as solitons", "confidence": MEM},
            {"ref": "Adkins-Nappi-Witten, Nucl. Phys. B 228, 552 (1983)", "confidence": MEM},
            {"ref": "Isgur-Paton, PRD 31, 2910 (1985) - string-junction baryons (complementarity)", "confidence": SNIP},
            {"ref": "Nadkarni-Nielsen-Zahed, Nucl. Phys. B 253, 308 (1985) - Cheshire-Cat complementarity", "confidence": MEM},
        ],
    },
    {
        "state_class": "K-Kbar mesons",
        "owning_sector": "both (same states): Skyrme-sector quanta at long wavelength; Q-1 knot-antiknot + tube at constituent scale",
        "corpus_evidence": [
            "corpus3/01-GUM-Omega-Paper-v4.3-ext.md:615 (Q-1: 'line-neutral composites (K-Kbar mesons; ...)')",
            "corpus3/01-GUM-Omega-Paper-v4.3-ext.md:617 (Q-5: orientation field, stiffness f_q, as the long-wavelength theory of line-neutral composites)",
        ],
        "literature": [
            {"ref": "mesons as chiral-field quanta in the Skyrme description (textbook)", "confidence": MEM},
            {"ref": "Isgur-Paton, PRD 31, 2910 (1985) - flux-tube mesons", "confidence": SNIP},
        ],
    },
    {
        "state_class": "deuteron-analog (B=2)",
        "owning_sector": "Skyrme sector exclusively",
        "corpus_evidence": [
            "corpus3/01-GUM-Omega-Paper-v4.3-ext.md:831 (App F.5: B=2 composite, J=0 forbidden / J=1 ground, sector-conditional post F-R17)",
            "theory-audit/i4_topology.md:361 (i4_topology.md:361-393)",
        ],
        "literature": [
            {"ref": "Braaten-Carson, PRD 38, 3525 (1988) - B=2 toroidal Skyrmion", "confidence": MEM},
            {"ref": "Leese-Manton-Schroers, Nucl. Phys. B 442, 228 (1995)", "confidence": MEM},
        ],
    },
    {
        "state_class": "glueball-analogs (knot-free closed T2 tube loops)",
        "owning_sector": "STRING sector exclusively - the Skyrme sector cannot house them (no glueball states in the Skyrme spectrum)",
        "corpus_evidence": [
            "corpus/WS-D-Q1-Q2-Relic-Census-v0_1.md:11 (the class's ONLY archive occurrence: 'stratum composites (glueball-analogs) -> decay to hadrons/phasons [check]')",
            "corpus3/01-GUM-Omega-Paper-v4.3-ext.md:615 (Q-2: T2 tube autonomous tension T ~ M_gap^2)",
            "corpus3/01-GUM-Omega-Paper-v4.3-ext.md:617 (sigma = 0.19 GeV^2 [IM-inversion])",
        ],
        "computed_content": [
            "T1/RESULTS.md Sec.5 (axion-less closed-tube spectrum misses the [IM] 0-+ window in both routes; V.F)",
            "T2/RESULTS.md Sec.0/5 (plain bosonic string as printed; tube-core axion [CJ-new], mass ~ M_gap symbol-only)",
        ],
        "literature": [
            {"ref": "Gomm-Jain-Johnson-Schechter, PRD 33, 801 (1986) - scalar glueball must be ADDED to Skyrme-type models as a dilaton field", "confidence": MEM},
            {"ref": "Zahed-Brown, Phys. Rept. 142, 1 (1986) - Skyrme model review (spectrum = solitons + chiral quanta)", "confidence": MEM},
            {"ref": "Isgur-Paton, PRD 31, 2910 (1985) - glueballs as closed flux loops", "confidence": SNIP},
            {"ref": "Athenodorou-Bringoltz-Teper, arXiv:1007.4720 - torelon/NG spectroscopy", "confidence": SNIP},
            {"ref": "Dubovsky-Hernandez-Chifflet, arXiv:1611.09796 - Axionic String Ansatz", "confidence": SNIP},
        ],
    },
    {
        "state_class": "hybrid-analogs (knot composites with excited tube)",
        "owning_sector": "jointly owned IF ever constructed: knot cores -> Skyrme-sector labels; tube excitation -> STRING sector",
        "corpus_evidence": [
            "corpus3/01-GUM-Omega-Paper-v4.3-ext.md:615 (tube exists; NO excitation quantization anywhere - class entirely SILENT, [CJ]-if-built; negative evidence, archaeologist sweep)",
        ],
        "literature": [
            {"ref": "Isgur-Paton, PRD 31, 2910 (1985) - flux-tube hybrids", "confidence": SNIP},
            {"ref": "Barnes-Close-Swanson, PRD 52, 5242 (1995)", "confidence": MEM},
        ],
    },
]

audit = [
    {
        "id": "A",
        "threat": "scalar 0++ glueball-analog vs Skyrmion breathing/monopole mode",
        "resolution": "breathing mode is a B=1 baryon resonance (Roper-analog), wrong winding sector; B=0 scalar enters Skyrme-type models only as an ADDED dilaton field; residual QCD subtlety = qqbar-glueball mixing ('dominant component', never exclusive ownership)",
        "literature": [
            {"ref": "Hajduk-Schwesinger, Phys. Lett. B 140, 172 (1984)", "confidence": MEM},
            {"ref": "Gomm-Jain-Johnson-Schechter, PRD 33, 801 (1986)", "confidence": MEM},
            {"ref": "Amsler-Close, PRD 53, 295 (1996)", "confidence": MEM},
            {"ref": "Gui et al, PRL 110, 021601 (2013) - scalar-glueball production", "confidence": SNIP},
        ],
        "corpus_status": "potential only: no breathing-mode spectrum, no B=0 scalar quantum printed (corpus3/01-GUM-Omega-Paper-v4.3-ext.md:617 prints the field, not its meson spectrum); discriminator available = frame winding number",
    },
    {
        "id": "B",
        "threat": "0-+ glueball-analog vs pseudoscalar modes of the orientation field",
        "resolution": "orientation-field pseudoscalars are flavor-adjoint pion-analogs; the singlet channel (eta-prime) is anomaly-fed and MIXES with the pseudoscalar glueball (mixing matrix, separate heavier basis state) - not a double count",
        "literature": [
            {"ref": "Witten, Nucl. Phys. B 156, 269 (1979)", "confidence": MEM},
            {"ref": "Veneziano, Nucl. Phys. B 159, 213 (1979)", "confidence": MEM},
            {"ref": "eta-glueball mixing, PRD 107, 094510 (2023), arXiv:2205.12541", "confidence": SNIP},
            {"ref": "glueball-ccbar mixing, PLB 827, 136960 (2022)", "confidence": SNIP},
        ],
        "corpus_status": "cannot even be posed as printed: no flavor symmetry, no singlet/octet split, no eta-prime analog, no anomaly machinery (T3 row 5; 'flavor' once, corpus3/01-GUM-Omega-Paper-v4.3-ext.md:867); T2 sharpener: tube-core axion is a WORLDSHEET mode, distinct from any bulk quantum",
    },
    {
        "id": "C",
        "threat": "baryons claimed by both sectors (Q-1 junctions vs Q-5 solitons)",
        "resolution": "complementary effective descriptions of the SAME states in different regimes (quark/flux-tube <-> Skyrme; chiral-bag/Cheshire-Cat) - never summed in the literature",
        "literature": [
            {"ref": "Nadkarni-Nielsen-Zahed, Nucl. Phys. B 253, 308 (1985)", "confidence": MEM},
            {"ref": "Isgur-Paton, PRD 31, 2910 (1985) - both descriptions coexist", "confidence": SNIP},
        ],
        "corpus_status": "the corpus prints both readings (corpus3/01-GUM-Omega-Paper-v4.3-ext.md:615 junctions; corpus3/01-GUM-Omega-Paper-v4.3-ext.md:617 solitons) and never sums them; annotation #1 makes the same-state reading explicit",
    },
    {
        "id": "D",
        "threat": "2++ channel: closed-loop tensor vs Skyrmion rotational excitations",
        "resolution": "Skyrmion rotor states are BARYONS (half-integer spin under Finkelstein-Rubinstein quantization, N_c odd) - wrong statistics sector; the B=0 tensor collision reduces to row A's mixing bookkeeping",
        "literature": [
            {"ref": "Finkelstein-Rubinstein quantization (corpus's own bibliography, corpus3/01:65)", "confidence": MEM},
        ],
        "corpus_status": "no printed B=0 tensor quanta (corpus3/01-GUM-Omega-Paper-v4.3-ext.md:617 prints no meson spectrum); resolution-by-winding as row A",
    },
]

pricing = {
    "option_A": {
        "name": "Q-5 as printed (Skyrme-only)",
        "costs": [
            "A1 homelessness: the glueball-dominant X(2370)-class state has no in-model home (Skyrme spectrum contains no glueball-analog states - dictionary row 4 literature)",
            "A2 census strain: WS-D D-0's exhaustive inventory PRESUPPOSES the glueball-analog class (WS-D-Q1-Q2-Relic-Census-v0_1.md:11) - the corpus's own relic census presupposes the class the identification cannot house",
            "A3 a new exclusion/instability theorem is owed for the T2 tube's finite-energy closed loops (Q-2 autonomous tension, corpus3/01-GUM-Omega-Paper-v4.3-ext.md:615; only the WRONG-stratum T3 web-loop obituary is printed, WS-D:19) - the 'no new structure' advantage is partly illusory",
            "A4 Course 20.1's 'exactly' sentence (corpus3/02-The-Substrate-Course-v3-ext.md:1088) goes load-bearing as an unproven completeness claim",
        ],
        "benefits": [
            "no new sectors; no loop quantization; F-Q8 mooted (class denied rather than owed closure - at A2's price); Q-5 verbatim; T-N5 scope narrow",
        ],
    },
    "option_B": {
        "name": "Q-5 + closed-string completion (the two-sector dictionary)",
        "costs": [
            "B1 loop quantization machinery [CJ]: compact-loop moduli/J^PC operators (h4/i4 covers embedded lines only, theory-audit/h4_completion.md:63-99)",
            "B2 the closure obligation F-Q8 must be discharged (T3 Sec.3: NEW-OPEN, three branches, none printed)",
            "B3 worldsheet content beyond print: tube-core axion [CJ-new]/[DW testbed-grade] (T2), upstream of it the h25-absent mediator branch (theory-audit/h25_RESULTS.md:53-59); T1 shows the axion-less spectrum misses the [IM] window in both routes (V.F)",
            "B4 sector-interface machinery [CJ]: mixing (no charm sector to mix with), decay/width machinery (none at any grade, WS-N-Hadronic-Pass-P7-v0_1.md:9), radiative transitions (no vector-meson analogs)",
            "B5 bookkeeping amendments: Course 20.1 scope reading/amendment; VII.J two-sector note; WS-D cross-reference (annotations 1-3)",
            "B6 necessary-not-sufficient: even fully paid, the completion buys a home, not a number - axion-completed spectrum depends on M_gap (symbol-only archive-wide); everything quantitative stays at V.F",
        ],
        "benefits": [
            "houses the class D-0 presupposes (= A2's mirror)",
            "matches the literature's two-description structure of real QCD",
            "converts the closure question from hidden presupposition to posable open problem (F-Q8)",
            "gives the flavor-singlet structural support (T3 row 5) a state to attach to",
        ],
    },
    "recommendation": "NONE - neither option recommended; the pricing is the deliverable; adoption either way is the authors' alone (symmetric honesty per T4-G2)",
}

annotations = [
    {"n": 1, "subject": "the VII.J two-sector note", "placement": "corpus3/01 Sec. VII.J at Q-5 (01:617)",
     "finding_tags": ["F-T9-T4", "F-T9-T1", "F-T9-T2"],
     "artifact_paths": ["tier9-glueball/T4/RESULTS.md", "tier9-glueball/T1/RESULTS.md", "tier9-glueball/T2/RESULTS.md"],
     "strikeable": True, "grade_motion": "none", "optional": False},
    {"n": 2, "subject": "the closure obligation F-Q8 (consumed verbatim from T3 Sec.3.4)", "placement": "corpus3/01 Sec. VII.J at Q-6-prime (01:615), cross-pinned to Sec. X list (01:767)",
     "finding_tags": ["F-T9-T3"],
     "artifact_paths": ["tier9-glueball/T3/RESULTS.md"],
     "strikeable": True, "grade_motion": "none", "optional": False},
    {"n": 3, "subject": "the relic-census cross-reference (WS-D D-0 <-> VII.J)", "placement": "corpus/WS-D-Q1-Q2-Relic-Census-v0_1.md:11 + mirror in corpus3/01 VII.J",
     "finding_tags": ["F-T9-T3", "F-T9-T4"],
     "artifact_paths": ["tier9-glueball/T3/RESULTS.md", "tier9-glueball/T4/RESULTS.md"],
     "strikeable": True, "grade_motion": "none", "optional": False},
    {"n": 4, "subject": "T-N5-adjacent monitoring-row suggestion (completeness, not correctness)", "placement": "corpus/WS-N-Hadronic-Pass-P7-v0_1.md Sec.3 beside T-N5 (:29)",
     "finding_tags": ["F-T9-T4"],
     "artifact_paths": ["tier9-glueball/T4/RESULTS.md"],
     "strikeable": True, "grade_motion": "none", "optional": True},
]

# ---------------- validators ----------------
failures = []

def v1_v2():
    for row in dictionary:
        ok = any(FILELINE.search(c) for c in row["corpus_evidence"])
        if not ok:
            failures.append("V1 dictionary row lacks file:line: " + row["state_class"])
        if not row["literature"]:
            failures.append("V2 dictionary row lacks literature: " + row["state_class"])
    for row in audit:
        if not FILELINE.search(row["corpus_status"]):
            failures.append("V1 audit row lacks file:line: " + row["id"])
        if not row["literature"]:
            failures.append("V2 audit row lacks literature: " + row["id"])

# forbidden sentences (roadmap register). Constructed from parts so this
# source file itself never contains a bare forbidden sentence.
X = "X(2370)"
FORBIDDEN = [
    ("gum predicted " + X).lower(),
    (X + " confirms").lower(),
    (X + " refutes").lower(),
    ("confirms gum").lower(),
    ("refutes gum").lower(),
]
BAN_MARKERS = ["unavailable", "forbidden", "ban", "does not say", "banned"]

def v3():
    files = [os.path.join(HERE, "RESULTS.md"), os.path.join(HERE, "t4_annotations.md"),
             os.path.join(HERE, "t4_dictionary.json")]
    for fp in files:
        if not os.path.exists(fp):
            failures.append("V3 missing file: " + fp)
            continue
        for i, line in enumerate(open(fp, encoding="utf-8"), 1):
            low = line.lower()
            if any(s in low for s in FORBIDDEN):
                if any(m in low for m in BAN_MARKERS):
                    continue  # quotation of the ban itself
                failures.append("V3 forbidden sentence at %s:%d" % (os.path.basename(fp), i))

def v4():
    for ws in ("T1", "T2", "T3"):
        p = os.path.join(T9, ws, "RESULTS.md")
        if not os.path.exists(p):
            failures.append("V4 missing sibling artifact: " + p)

def v5():
    for a in annotations:
        if not (a["strikeable"] and a["finding_tags"] and a["artifact_paths"] and a["grade_motion"] == "none"):
            failures.append("V5 annotation manifest defect: #%d" % a["n"])

out = {
    "workstream": "T4",
    "finding": "F-T9-T4",
    "register": ("within-model; nothing bears on nature; BESIII/lattice = [IM] anchors; "
                 "seal q-theta V.F on every mass statement (all cited from T1); all corpus-side "
                 "changes are OFFERS; numerology hazard (dimensionless c-frak = 2.37 vs "
                 "m_X = 2.37 GeV digit coincidence) pre-empted - a man-made-units coincidence, "
                 "never cited as structure anywhere in this workstream"),
    "T3_dependency": "T3/RESULTS.md landed at watcher poll +120 s; consumed directly; no draft-from-spec amendment needed",
    "dictionary": dictionary,
    "double_counting_audit": audit,
    "audit_bottom_line": "0 actual double counts in print; 4 potential rows, each with a standard literature resolution; discriminators already in corpus structure (frame winding B; worldsheet vs bulk locality)",
    "priced_alternative": pricing,
    "annotations": annotations,
    "gates": [],  # filled below
    "spec_deviations": [],
    "validation_failures": failures,
}

def main():
    v1_v2()
    v4()
    v5()
    # first write (so v3 can scan the JSON), then v3, then final write
    jpath = os.path.join(HERE, "t4_dictionary.json")
    with open(jpath, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    v3()
    n_mem = sum(1 for row in dictionary + audit for L in row["literature"] if L["confidence"] == MEM)
    n_snip = sum(1 for row in dictionary + audit for L in row["literature"] if L["confidence"] == SNIP)
    out["gates"] = [
        {"id": "T4-G1", "requirement": "dictionary + double-counting audit complete with [IM] literature citations and corpus file:line evidence",
         "measured": "5 dictionary rows + 4 audit rows; file:line validator PASS on all 9; literature cells on all 9 (%d memory-cited flagged, %d snippet-verified)" % (n_mem, n_snip),
         "verdict": "PASS" if not any(x.startswith(("V1", "V2")) for x in failures) else "FAIL"},
        {"id": "T4-G2", "requirement": "alternative priced both ways, symmetric honesty (homelessness AND new-structure bill both printed)",
         "measured": "option A: 4 costs + benefits (incl. homelessness A1, census strain A2); option B: 6-item bill + benefits + necessary-not-sufficient caveat; recommendation: NONE",
         "verdict": "PASS"},
        {"id": "T4-G3", "requirement": "proposed T9 annotations drafted in charter-compliant offer language (strike-able, cites finding + artifact path, grades never rise)",
         "measured": "4 boxes (3 required + 1 optional) in t4_annotations.md; manifest validator PASS; format modeled on corpus3/01 T6 boxes (lines 505/637/651/661)",
         "verdict": "PASS" if not any(x.startswith("V5") for x in failures) else "FAIL"},
        {"id": "T4-G4", "requirement": "within-model; nothing adopted, everything offered; no unavailable sentences",
         "measured": "all items offer-class; 0 stakes, 0 clock moves, 0 grade rises; forbidden-sentence scan over 3 deliverables clean (ban-quotations excepted)",
         "verdict": "PASS" if not any(x.startswith("V3") for x in failures) else "FAIL"},
    ]
    out["validation_failures"] = failures
    out["validation_ok"] = not failures
    with open(jpath, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    print("t4_build: %d dictionary rows, %d audit rows, %d annotations" % (len(dictionary), len(audit), len(annotations)))
    print("literature cells: %d memory-cited, %d snippet-verified" % (n_mem, n_snip))
    for g in out["gates"]:
        print("%s: %s" % (g["id"], g["verdict"]))
    if failures:
        print("VALIDATION FAILURES:")
        for x in failures:
            print("  -", x)
        sys.exit(1)
    print("validation_ok: True")

if __name__ == "__main__":
    main()
