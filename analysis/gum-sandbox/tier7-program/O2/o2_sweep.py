#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
O2 — g_lim archaeology (F-T7-O2): exhaustive archive sweep for any printed
Majoron-mode 0nubb rate-level bound (g_lim / T_1/2-equivalent / quantitative
detectability line) — the F-T7-N2 blocking number.

Protocol (pre-registered, ROADMAP_v9 Phase-O addendum O2):
  - Scope: corpus/ corpus2/ corpus3/ theory-audit/ DISCHARGE_PACKAGE/
    gum-core/ tier0-gauntlet/ ... tier7-program/ + top-level *.md of
    analysis/gum-sandbox/.  Supplementary coverage pass (beyond mandate):
    simulator/ substrate-suite/ audit-logs/ + the top-level html.
  - >= 2 independent pattern families:
      family A (mode/observable), family B (bounds/limits/experiments),
      family C (numerics near mode context).  Every pattern logged with its
      raw line count; context filters logged.
  - Every context-passing hit recorded (path, line, verbatim quote) and
    classified:
      class 1 — Majoron-mode RATE-level bound (g_lim or T_1/2-equivalent or a
                quantitative detectability line): the blocking number;
      class 2 — non-rate g-constraint (free-streaming, BBN, SN cooling,
                recoupling, lab perimeter — binds g itself, not the mode rate);
      class 3 — stake/kill language about 0nubb OCCURRENCE (S6 frame; funnels
                as occurrence adjudicators; incl. the mass-mode T_1/2 data
                context, flagged near_miss);
      class 4 — other/incidental (campaign analysis layer self-references,
                unrelated detectability language, code, etc.).
  - Decision rule: class-1 found -> run the F-T7-N2 conditional-kill
    arithmetic against it; only class-2/3 -> catalog + why they cannot
    adjudicate; nothing -> certify ABSENT with the protocol printed.

Deterministic, stdlib-only.  Output: o2_results.json + stdout protocol log.
"""
import json, os, re, sys, unicodedata
from collections import OrderedDict

ROOT = "/home/user/fork_frankensim/analysis/gum-sandbox"
OUT  = os.path.join(ROOT, "tier7-program", "O2", "o2_results.json")

MANDATED_DIRS = ["corpus", "corpus2", "corpus3", "theory-audit",
                 "DISCHARGE_PACKAGE", "gum-core",
                 "tier0-gauntlet", "tier1-spectrum", "tier2-closure",
                 "tier3-born", "tier4-field", "tier5-family",
                 "tier6-foundations", "tier7-program",
                 # promoted to FULL scope after the first pass showed it holds
                 # the v2.0.1 corpus documents (Omega paper, Primer, Course):
                 "substrate-suite"]
SUPPLEMENTARY_DIRS = ["simulator", "audit-logs"]

TEXT_EXT = {".md", ".txt", ".py", ".json", ".rs", ".toml", ".sh", ".qmd",
            ".html", ".yml", ".yaml", ".csv", ".cfg", ".lock", ""}
SKIP_EXT = {".png", ".pdf", ".npz", ".npy", ".pyc", ".gz", ".zip", ".ico"}

SELF_PREFIX = os.path.join("tier7-program", "O2")  # exclude this workstream's own outputs

# ---------------------------------------------------------------- patterns
# family A — mode / observable
FAM_A = OrderedDict([
    ("A1_0nubb",       r"0nubb|0νββ"),
    ("A2_neutrinoless",r"neutrinoless|double[-\s]?beta"),
    ("A3_majoron",     r"majoron"),
    ("A4_halflife",    r"T₁/₂|T_?1/2|half[-\s]?life"),
    ("A5_Pnu4",        r"P-ν4|P-nu4|N-ν4|N-nu4"),
])
# family B — bounds / limits / experiments (independent of family A)
FAM_B = OrderedDict([
    ("B1_glim",        r"\bg_?lim\b"),   # \b: 'glimpse(d)' in the Primer prose is not g_lim
    ("B2_g_ineq",      r"\bg\s*[<≤≲]|\bg²\s*[<≤≲]|\bg\^?2\s*[<≤≲]"),
    ("B3_limit_on_g",  r"limit on g|bound on g"),
    ("B4_experiments", r"KamLAND|EXO-?200|nEXO|GERDA|LEGEND-?1000|NEMO|CUORE"
                       r"|Majorana Demonstrator|¹³⁶Xe|136Xe|Xe-136|⁷⁶Ge|76Ge|Ge-76"),
    ("B4b_LEGEND_cs",  r"LEGEND"),          # case-SENSITIVE (else figure legends)
    ("B5_detectab",    r"detectab|sensitivit"),
    ("B6_funnel",      r"funnel"),
])
CASE_SENSITIVE = {"B4b_LEGEND_cs"}
# family C — numerics near mode context
FAM_C = OrderedDict([
    ("C1_1em9",        r"10⁻⁹|10\^-?9|10\^\{-9\}|1e-0?9|e-09"),
    ("C2_1e26yr",      r"10²⁶|10\^26|10\^\{26\}|3\.8e26"),
])

# hits from these patterns are kept only if a mode token appears within
# +/- CONTEXT_WINDOW lines (else the pattern drowns in unrelated text —
# e.g. 'undetectable medium', matplotlib legends, code inequalities);
# raw counts are still logged for every pattern.
CONTEXT_FILTERED = {"B2_g_ineq", "B5_detectab", "B6_funnel", "C1_1em9"}
CONTEXT_WINDOW = 5
MODE_TOKEN = re.compile(
    r"majoron|0νββ|0nubb|neutrinoless|double[-\s]?beta|ΔL\s?=\s?2|\bS6\b|phason",
    re.IGNORECASE)

# ---------------------------------------------------------------- classifier
def layer_of(rel):
    if rel.startswith(("corpus/", "corpus2/", "corpus3/", "substrate-suite/")):
        return "corpus"   # substrate-suite = the authors' v2.0.1 release documents
    return "campaign-analysis"

RULES_NOTE = {
 1: "Majoron-mode RATE-level bound — the blocking number",
 2: "non-rate g-constraint (binds the coupling itself, not the mode-rate observable)",
 3: "0nubb OCCURRENCE stake/kill frame (S6) — occurrence is not a rate number",
 4: "other/incidental",
}

SIG_NOTE = ("Majoron dictionary + battery stamps (class-2: bind g itself); the same line "
            "carries P-nu4 'Majoron-mode 0nubb at g ~ 1e-9 (far-horizon)' — a QUALITATIVE "
            "signature clause: no detectability threshold, no rate normalization printed")

# Manual adjudications for the load-bearing / ambiguous prints (file, line) -> (class, flag)
OVERRIDES = {
 ("corpus2/01-GUM-Omega-Paper-v3.0-ext.md", 567): (2, SIG_NOTE),
 ("corpus3/01-GUM-Omega-Paper-v4.1-ext.md", 591): (2, SIG_NOTE),
 ("substrate-suite/01-GUM-Omega-Paper-v2.0.1.md", 760): (2,
   "V15.5 (soft sector / Majoron) verification row — battery arithmetic stamps only; "
   "no rate-level bound"),
 ("corpus/WS-nu-Q1-Reading-v0_1.md", 27): (4,
   "P-nu4 located at its corpus grade: '(Majoron-mode 0nubb, far-horizon)' — signature "
   "language only, no number"),
 ("corpus2/deltas/DELTA-NR-K2b-K3a-Two-Scale-Preview-v0_1.md", 11): (3,
   "corpus-side dated delta: the K3a gate returned ADVERSE with R = 1e14-1e18; four-cell "
   "occurrence table re-opened; contains the R statistic, NOT a bound"),
 ("corpus/WS-K4-Ledger-Perimeter-Audit-v0_1.md", 17): (3,
   "NEAR-MISS: contains the archive's ONLY experimental 0nubb half-life datum "
   "T1/2(136Xe) > 3.8e26 yr (KamLAND-Zen) + m_bb < 28-122 meV + LEGEND-1000/nEXO reach "
   "— all MASS-MODE (0nu peak) data context for the S6 occurrence stake row; "
   "not a Majoron-mode bound; no in-archive T1/2(g) conversion exists"),
 ("corpus/WS-nu-Neutrino-Ledger-Prompt-Pack-v0_1.md", 50): (3,
   "NOT a bound — the OPPOSITE: quarantine rule (q-kappa) '0nubb half-life numerology "
   "sealed (no T-half arithmetic against corpus constants)' — the corpus structurally "
   "FORBIDS itself half-life arithmetic; surfaced by the class-1 candidate detector, "
   "adjudicated by hand"),
 ("corpus/WS-K4-Ledger-Perimeter-Audit-v0_1.md", 41): (3,
   "NEAR-MISS: signed-margin table row r3 '0nubb | stake row (S6) | funnels: > 3.8e26 yr | —' "
   "— the funnel datum is the mass-mode half-life; margin column deliberately EMPTY "
   "(stake row, not a compliance row)"),
}

def classify(rel, lineno, line, patterns):
    """Ordered rules; returns (cls, note)."""
    key = (rel, lineno)
    if key in OVERRIDES:
        return OVERRIDES[key]
    low = line.lower()
    layer = layer_of(rel)
    majoronish = bool(re.search(r"majoron|phason|P-ν4|P-nu4", line, re.I))
    nubbish = bool(re.search(r"0νββ|0nubb|neutrinoless|double[-\s]?beta|ΔL\s?=\s?2", line, re.I))

    # ---- class-1 candidate detector (surfaced for manual adjudication) ----
    # a Majoron-context line carrying a half-life / g_lim / experiment-bound
    # numeric would land here; every candidate is re-examined by hand and
    # either promoted to class 1 or dispatched by an OVERRIDE above.
    cand = (majoronish and
            re.search(r"g_?lim|T₁/₂|T_?1/2|half[-\s]?life|×10²⁶|e26|yr\b", line, re.I))
    if cand and layer == "corpus":
        return (1, "CLASS-1 CANDIDATE — manual adjudication required")

    # ---- class 2: the battery + perimeter (binds g, not the mode rate) ----
    if (majoronish or "b-ν1" in low or "b-nu1" in low) and re.search(
        r"free.?stream|B-ν1|B-nu1|BBN|ΔN_eff|dneff|SN.?cool|recoupl|"
        r"g²\s*≲|g\^?2\s*(≲|<|<=)|PNC|perimeter|ε_e ≥|1\.5×10⁻⁶|window|"
        r"∈ \[8×10⁻¹⁰|∈ \[1\.7×10⁻⁹|dictionary", line, re.I):
        return (2, "battery/perimeter/window print: constrains g itself "
                   "(astro/cosmo/lab emission or admission window), not the 0nubb mode rate")

    # ---- class 3: occurrence stake frame ----
    if re.search(r"\bS6\b|stake|occurr|funnel|kill|LEGEND|nEXO|KamLAND|hash", line) and \
       (nubbish or majoronish or re.search(r"funnel|LEGEND|nEXO|KamLAND|\bS6\b", line)):
        return (3, "S6 occurrence stake / funnel-adjudicator / kill language "
                   "— stakes occurrence, never a half-life (WS-nu-Q1 verbatim)")

    # ---- class 4 subtypes ----
    if layer == "campaign-analysis":
        return (4, "campaign analysis layer (audit/tier/roadmap/code) — reports on or "
                   "searches for the bound; not a corpus print")
    if re.search(r"undetectab|detectable by uniform motion|foliation|Unruh|aether", line, re.I):
        return (4, "unrelated detectability language (medium/foliation/Unruh)")
    if majoronish or nubbish:
        return (4, "incidental mode mention without bound/rate/stake content")
    return (4, "incidental pattern hit outside the 0nubb/Majoron question")

# ---------------------------------------------------------------- sweep
def iter_files(dirs, add_toplevel_md=True, add_toplevel_html=False):
    files = []
    for d in dirs:
        p = os.path.join(ROOT, d)
        if not os.path.isdir(p):
            continue
        for base, _, names in os.walk(p):
            for n in names:
                fp = os.path.join(base, n)
                rel = os.path.relpath(fp, ROOT)
                if rel.startswith(SELF_PREFIX):
                    continue
                ext = os.path.splitext(n)[1].lower()
                files.append((rel, ext not in SKIP_EXT))
    if add_toplevel_md:
        for n in sorted(os.listdir(ROOT)):
            if n.endswith(".md") and os.path.isfile(os.path.join(ROOT, n)):
                files.append((n, True))
    if add_toplevel_html:
        for n in sorted(os.listdir(ROOT)):
            if n.endswith(".html"):
                files.append((n, True))
    return files

def read_lines(rel):
    try:
        with open(os.path.join(ROOT, rel), "r", encoding="utf-8", errors="replace") as f:
            return f.read().splitlines()
    except (OSError, UnicodeError):
        return None

def run_sweep(file_list, families):
    """Returns (raw_counts, hits) where hits: (rel, lineno) -> {line, patterns}."""
    raw_counts = OrderedDict((k, 0) for fam in families for k in fam)
    hits = OrderedDict()
    compiled = {}
    for fam in families:
        for k, pat in fam.items():
            flags = 0 if k in CASE_SENSITIVE else re.IGNORECASE
            compiled[k] = re.compile(pat, flags)
    n_scanned = 0
    for rel, is_text in file_list:
        if not is_text:
            continue
        lines = read_lines(rel)
        if lines is None:
            continue
        n_scanned += 1
        for i, line in enumerate(lines, 1):
            matched = [k for k in compiled if compiled[k].search(line)]
            if not matched:
                continue
            for k in matched:
                raw_counts[k] += 1
            kept = []
            for k in matched:
                if k in CONTEXT_FILTERED:
                    lo, hi = max(0, i - 1 - CONTEXT_WINDOW), min(len(lines), i + CONTEXT_WINDOW)
                    if not any(MODE_TOKEN.search(l) for l in lines[lo:hi]):
                        continue
                kept.append(k)
            if kept:
                key = (rel, i)
                if key not in hits:
                    hits[key] = {"line": line, "patterns": []}
                hits[key]["patterns"] = sorted(set(hits[key]["patterns"] + kept))
    return raw_counts, hits, n_scanned

def main():
    families = [FAM_A, FAM_B, FAM_C]
    mandated = iter_files(MANDATED_DIRS, add_toplevel_md=True)
    counts_by_dir = OrderedDict()
    for d in MANDATED_DIRS:
        counts_by_dir[d] = sum(1 for rel, _ in mandated if rel.startswith(d + "/") or rel.startswith(d + os.sep))
    counts_by_dir["<top-level *.md>"] = sum(1 for rel, _ in mandated if os.sep not in rel and "/" not in rel)

    raw_counts, hits, n_scanned = run_sweep(mandated, families)

    # supplementary coverage pass (beyond mandate): family A only
    supp = iter_files(SUPPLEMENTARY_DIRS, add_toplevel_md=False, add_toplevel_html=True)
    supp_counts, supp_hits, supp_scanned = run_sweep(supp, [FAM_A])

    # classify
    catalog = []
    tallies = {1: 0, 2: 0, 3: 0, 4: 0}
    candidates = []
    for (rel, lineno), h in sorted(hits.items()):
        cls, note = classify(rel, lineno, h["line"], h["patterns"])
        if cls == 1:
            candidates.append((rel, lineno))
        tallies[cls] += 1
        q = h["line"].strip()
        entry = OrderedDict([
            ("file", rel), ("line", lineno), ("layer", layer_of(rel)),
            ("patterns", h["patterns"]), ("class", cls), ("note", note),
            ("quote", q),   # verbatim, untruncated (archaeology grade)
        ])
        catalog.append(entry)

    # -------- verdict per pre-registered decision rule --------
    class1 = [e for e in catalog if e["class"] == 1]
    verdict = OrderedDict()
    if class1:
        verdict["branch"] = "class-1 FOUND -> conditional-kill arithmetic required"
        verdict["entries"] = class1
    else:
        verdict["branch"] = ("class-1 ABSENT -> certify ABSENT: the archive prints NO "
                             "Majoron-mode 0nubb rate-level bound (no g_lim, no Majoron-mode "
                             "T_1/2, no quantitative detectability line, no T_1/2(g) rate "
                             "normalization for ANY 0nubb mode). The F-T7-N2 conditional kill "
                             "(fires for any g_lim < 8e-3) remains ARMED but cannot fire "
                             "in-archive.")
        verdict["nearest_misses"] = [
            OrderedDict([
              ("file", "corpus/WS-K4-Ledger-Perimeter-Audit-v0_1.md"), ("lines", [17, 41]),
              ("datum", "T1/2(136Xe) > 3.8e26 yr (KamLAND-Zen); m_bb < 28-122 meV; "
                        "LEGEND-1000/nEXO reach the inverted-hierarchy funnel"),
              ("why_not_class_1", [
                "wrong mode: it is the MASS-mechanism (0nu two-electron peak) limit; the "
                "Majoron mode (0nubb+chi, continuous spectrum) has separate experimental "
                "limits, none printed anywhere in the archive",
                "no conversion printed: the archive contains no T_1/2(g) normalization "
                "(no nuclear matrix element, no phase-space factor) for any mode, so no "
                "half-life can be converted to a g-bound in-archive — and the corpus prints "
                "no predicted T_1/2 at its g for the honest rate R x printed to be compared "
                "against 3.8e26 yr",
                "role in the print: it is data context for the S6 OCCURRENCE stake "
                "('stake row — adjudication owned by S6; nothing moved'; margin column "
                "printed EMPTY), and the corpus states verbatim it stakes 'occurrence, "
                "never a half-life' (WS-nu-Q1-Reading line 31)"]),
            ])
        ]
        verdict["consequence_for_N2"] = (
            "F-T7-N2's verdict (iii) INDETERMINATE-FROM-ARCHIVE is CERTIFIED at archaeology "
            "grade: the blocking number g_lim is confirmed unprinted across the entire "
            "archive; the conditional kill stays filed, adjudicable only by an "
            "archive-external Majoron-mode rate bound.")
    verdict["class_2_3_cannot_adjudicate"] = (
        "class-2 items (B-nu1 free-streaming f >= 1.4 MeV / floor eps_e >= 1.5e-6, BBN "
        "DeltaN_eff = 0.0268, SN-cooling band [edges unprinted], late recoupling "
        "[unprinted], lab perimeter g^2 <~ 1e-16, Cs-PNC, the window g in [8e-10, 1.3e-8]) "
        "bind the COUPLING g via emission/thermalization/admission — they are blind to the "
        "0nubb-insert second moment R, which multiplies the mode-rate observable only; "
        "class-3 items (S6, four-cell table, funnels, 3.8e26 yr) adjudicate OCCURRENCE of "
        "standard-mode 0nubb — occurrence is not a rate number, and the only half-life "
        "printed is for the wrong mode with no in-archive conversion. Neither class can "
        "decide whether R x printed exceeds a Majoron-mode limit.")

    out = OrderedDict([
      ("workstream", "F-T7-O2"),
      ("date", "2026-07-18"),
      ("protocol", OrderedDict([
        ("scope_mandated_dirs", MANDATED_DIRS + ["<top-level *.md>"]),
        ("files_per_dir", counts_by_dir),
        ("files_scanned_text", n_scanned),
        ("files_listed_total", len(mandated)),
        ("pattern_families", OrderedDict([
          ("A_mode_observable", dict(FAM_A)),
          ("B_bounds_limits",   dict(FAM_B)),
          ("C_numerics",        dict(FAM_C))])),
        ("case_sensitive_patterns", sorted(CASE_SENSITIVE)),
        ("context_filtered_patterns", OrderedDict([
          ("patterns", sorted(CONTEXT_FILTERED)),
          ("rule", "kept only if MODE_TOKEN within +/- %d lines" % CONTEXT_WINDOW),
          ("mode_token_regex", MODE_TOKEN.pattern)])),
        ("raw_line_counts_per_pattern", raw_counts),
        ("supplementary_pass", OrderedDict([
          ("dirs", SUPPLEMENTARY_DIRS + ["<top-level *.html>"]),
          ("families", ["A only"]),
          ("files_scanned", supp_scanned),
          ("raw_counts", supp_counts),
          ("hits", [OrderedDict([("file", r), ("line", l),
                                 ("quote", h["line"].strip()[:300])])
                    for (r, l), h in sorted(supp_hits.items())])])),
      ])),
      ("class_definitions", {str(k): v for k, v in RULES_NOTE.items()}),
      ("tallies", OrderedDict([("class_%d" % k, tallies[k]) for k in (1, 2, 3, 4)]
                              + [("total_catalog", len(catalog))])),
      ("class1_candidates_surfaced", candidates),
      ("catalog", catalog),
      ("verdict", verdict),
    ])
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)

    # ----------------- stdout protocol log (G1) -----------------
    print("O2 g_lim archaeology — sweep protocol")
    print("scope files per dir:", json.dumps(counts_by_dir))
    print("text files scanned: %d of %d listed (binaries skipped by extension)"
          % (n_scanned, len(mandated)))
    print("supplementary pass: %d files, family-A hits: %d" % (supp_scanned, len(supp_hits)))
    print("raw line counts per pattern:")
    for k, v in raw_counts.items():
        print("  %-16s %5d" % (k, v))
    print("catalog: %d unique (file,line) hits; tallies: %s"
          % (len(catalog), {("class_%d" % k): tallies[k] for k in (1, 2, 3, 4)}))
    print("class-1 candidates surfaced for manual adjudication:", candidates or "NONE")
    print("VERDICT branch:", verdict["branch"].split("->")[0].strip())
    return 0

if __name__ == "__main__":
    sys.exit(main())
