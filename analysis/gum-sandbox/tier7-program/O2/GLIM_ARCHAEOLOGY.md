# O2 — g_lim archaeology: the F-T7-N2 blocking number (F-T7-O2)

**Goal.** F-T7-N2 (tier7-program/N2) closed with verdict (iii) INDETERMINATE-FROM-ARCHIVE
and a conditional kill: *any Majoron-mode 0νββ rate-level bound g_lim < 8×10⁻³ excludes
the honest second-moment rate at every anchor and every printed g*. N2's claim that the
archive prints NO such bound rested on J3 §0's source census. This workstream certifies
or refutes that claim at archaeology grade: an exhaustive, pattern-logged sweep of the
entire archive, every hit quoted verbatim and classified, the pre-registered decision
rule executed.

**Campaign:** Tier 7 Phase O (ROADMAP_v9_PROGRAM.md addendum, O2).
**Date:** 2026-07-18. **Code:** `o2_sweep.py` (deterministic, stdlib-only, ~2 s).
**Raw output:** `o2_results.json` (machine-readable catalog, 370 entries, verbatim
untruncated quotes). **No figure** (per spec — archival workstream).
**Prior art (honored):** N2/RESULTS.md (the conditional kill and the blocking number),
theory-audit/j3_RESULTS.md §0 (the battery source census), theory-audit/h21_RESULTS.md.

*Within-model; nothing here bears on nature.*

## 0. Verdict in one line

**ABSENT — CERTIFIED. The archive prints no Majoron-mode 0νββ rate-level bound of any
form: no g_lim, no Majoron-mode T₁/₂, no quantitative detectability line, and no
T₁/₂(g) rate normalization for ANY 0νββ mode that could convert one; 370 catalogued
hits across 658 scanned files contain 0 class-(1) items. The nearest miss — the
archive's only experimental 0νββ rate number, T₁/₂(¹³⁶Xe) > 3.8×10²⁶ yr (KamLAND-Zen,
WS-K4 lines 17/41) — is the MASS-mode limit, printed as data context for the S6
occurrence stake, and is triply disqualified (wrong mode; no in-archive conversion in
either direction; occurrence-stake role under the corpus's own q-κ seal: "0νββ
half-life numerology sealed"). The pre-registered decision rule therefore executes its
catalog branch: the class-(2) g-constraints and class-(3) occurrence stakes are
catalogued below with the precise reason each cannot adjudicate the honest-rate
question. F-T7-N2's conditional kill remains ARMED and in-archive UNFIREABLE — it can
only be adjudicated by an archive-external Majoron-mode rate bound.**

## 1. Gates

| id | spec | measured | verdict |
|---|---|---|---|
| G1 | sweep coverage: file counts per directory + full pattern log; ≥ 2 independent pattern families | 15 directories + top-level *.md: 705 files listed, 658 text files scanned (47 binaries skipped by extension, list in JSON); 3 pattern families, 14 patterns, every raw line count logged (§3); supplementary family-A pass over simulator/ + audit-logs/ + top-level html (21 files, 9 hits, all campaign-layer echoes); 2 protocol refinements logged (§3, r1–r2) | PASS |
| G2 | the catalog: every hit quoted + classified; the expected battery items each located and classified | 370 unique (file, line) hits catalogued with verbatim untruncated quotes + class + note in `o2_results.json`; tallies: class 1 = **0**, class 2 = 89, class 3 = 95, class 4 = 186; class-1 candidate detector surfaced 2 lines for manual adjudication, both dispatched (§4.3); every expected item located: B-ν1 f ≥ 1.4 MeV ✓, SN band (edges unprinted) ✓, BBN ✓, recoupling ✓, g-window [8×10⁻¹⁰, 1.3×10⁻⁸] ✓, lab perimeter g² ≲ 10⁻¹⁶ ✓, P-ν4 far-horizon ✓, S6/LEGEND-1000/nEXO stake rows ✓, the 3.8×10²⁶ yr funnel datum ✓ (§4) | PASS |
| G3 | verdict per pre-registered decision rule | class-(1) empty → ABSENT certified with reproducible protocol (§3); class-(2)/(3) catalogued with the observable distinction stated per item class (§5); the F-T7-N2 arithmetic is NOT executed (nothing to execute it against) and the kill's arming condition is restated unchanged | PASS |

## 2. Key numbers

| quantity | value |
|---|---|
| files listed / text-scanned (mandated + promoted scope) | 705 / 658 |
| patterns run (families A/B/C) + supplementary | 14 + 5 |
| unique catalogued hits (context-passing) | 370 (99 distinct files; 216 corpus-layer, 154 campaign-analysis-layer) |
| class 1 — Majoron-mode rate-level bound | **0** |
| class 2 — non-rate g-constraints | 89 |
| class 3 — occurrence stake/kill frame | 95 |
| class 4 — other/incidental | 186 |
| nearest-miss datum (mass-mode, S6 context) | T₁/₂(¹³⁶Xe) > 3.8×10²⁶ yr; m_ββ < 28–122 meV |
| the unfired kill threshold (N2) | g_lim < 8×10⁻³ (min g_equiv 8.2×10⁻³) |
| corpus's own printed coupling window (for scale) | g ∈ [8×10⁻¹⁰, 1.3×10⁻⁸] |

## 3. Search protocol (G1 — reproducible)

**Scope (file counts).** corpus/ 105 · corpus2/ 82 · corpus3/ 12 · theory-audit/ 94 ·
DISCHARGE_PACKAGE/ 5 · gum-core/ 175 · tier0-gauntlet/ 6 · tier1-spectrum/ 6 ·
tier2-closure/ 30 · tier3-born/ 11 · tier4-field/ 13 · tier5-family/ 4 ·
tier6-foundations/ 66 · tier7-program/ 77 (O2's own outputs excluded) ·
substrate-suite/ 8 · top-level *.md 11. Total 705; 658 text files read line-by-line
(binaries skipped by extension: .png/.pdf/.npz/...). Supplementary pass (beyond
mandate): simulator/, audit-logs/, top-level .html — family A only, 21 files.

**Pattern log (regex, case-insensitive unless noted; raw line counts in parentheses).**
Family A (mode/observable): A1 `0nubb|0νββ` (148) · A2 `neutrinoless|double[-\s]?beta`
(6) · A3 `majoron` (137) · A4 `T₁/₂|T_?1/2|half[-\s]?life` (14) · A5
`P-ν4|P-nu4|N-ν4|N-nu4` (59).
Family B (bounds/limits/experiments): B1 `\bg_?lim\b` (20) · B2
`\bg\s*[<≤≲]|\bg²\s*[<≤≲]|\bg\^?2\s*[<≤≲]` (35) · B3 `limit on g|bound on g` (0) · B4
`KamLAND|EXO-?200|nEXO|GERDA|LEGEND-?1000|NEMO|CUORE|Majorana Demonstrator|¹³⁶Xe|136Xe|
Xe-136|⁷⁶Ge|76Ge|Ge-76` (18) · B4b `LEGEND` case-SENSITIVE (18) · B5
`detectab|sensitivit` (125) · B6 `funnel` (71).
Family C (numerics near mode context): C1 `10⁻⁹|10^-9|1e-9|e-09` variants (381) · C2
`10²⁶|10^26|3.8e26` (10).

**Context filter (logged).** B2/B5/B6/C1 hits kept only if a mode token
(`majoron|0νββ|0nubb|neutrinoless|double beta|ΔL = 2|S6|phason`) occurs within ±5
lines — else those patterns drown in unrelated text (the Unruh sector's "undetectable
medium", matplotlib legends, code inequalities). Raw counts above are UNFILTERED; only
the catalog is filtered. All other patterns catalogued unfiltered.

**Protocol refinements (2, logged per campaign rule).** r1: B1 tightened to
word-boundary `\bg_?lim\b` after the first pass matched "glimpsed" in the Primer prose
(10 false positives — verified by hand to be the word "glimpse(d)"). r2:
substrate-suite/ promoted from supplementary to FULL scope after the first pass showed
it holds the authors' v2.0.1 release documents (Ω paper v2.0.1, Primer, Course) —
corpus-class content that the mandate's directory list omitted.

**Candidate-surfacing rule.** Any corpus-layer line matching a Majoron token AND a
rate-bound token (`g_lim|T₁/₂|half-life|×10²⁶|e26|yr`) is auto-flagged "CLASS-1
CANDIDATE — manual adjudication required" and the run fails loud until each is
hand-adjudicated in the override table (visible in `o2_sweep.py`, with reasons). Two
lines surfaced across the runs; both dispatched (§4.3). Final run: zero unadjudicated.

## 4. The catalog (G2 — condensed; complete in `o2_results.json`)

Distinct prints, deduplicated across document versions (v2.0.1 / v3.0-ext / v4.1-ext
carry verbatim-identical rows; every copy is separately catalogued in the JSON).

### 4.1 Class 2 — g-constraints that bind the coupling, not the mode rate (89 hits)

| item | verbatim quote (source) | why it cannot adjudicate |
|---|---|---|
| B-ν1 free-streaming | "free-streaming ✓ (B-ν1: f ≥ 1.4 MeV — passed with margin; inverted, the **cosmological floor B-ν1′: ε_e ≥ 1.5×10⁻⁶**)" (Ω VII.I, corpus2/01:567 = corpus3/01:591; formula B-ν1 = (T·M_Pl·m_ν²)^{1/4} in WS-nu-P4:19) | binds f (hence g) via cosmological free-streaming of the phason — says nothing about the 0νββ-insert rate observable |
| BBN | "BBN ✓; … ΔN_eff = (4/7)(10.75/106.75)^{4/3} = 0.0268 [CJ-thermal]" (same lines; WS-nu-P4:19) | thermalization history of the phason bath; g-side, rate-blind |
| SN-cooling | "SN-cooling band evaded from below ✓" (same lines) | band edges NEVER printed (J3 §3.1 standing); binds g via stellar emission |
| late recoupling | "late recoupling harmless ✓" (same lines) | no formula or number printed anywhere; binds g |
| the g-window | "The Majoron dictionary [DW]. g = m_ν/f ∈ [8×10⁻¹⁰, 1.3×10⁻⁸] across f ∈ [4, 60] MeV (II.H)" (Ω VII.I); operative variant "g ∈ [1.7×10⁻⁹, 1.3×10⁻⁸]" (Course 19.3, corpus2/02:1050) | an ADMISSION window for the coupling, not an experimental bound on the mode rate |
| lab perimeter | "static PV/spin forces … at g² ≲ 10⁻¹⁶ — negligible"; Cs-PNC "(λZ/a₀)² = 10⁻⁴–3×10⁻³ per unit admixture x" (WS-K2:37, b/c faces) | fifth-force/PNC pricing of the ACTUAL coupling in the lab; not a decay-rate observable |
| P-ν4 (rides on the dictionary line) | "P-ν4: Majoron-mode 0νββ at g ∼ 10⁻⁹ (far-horizon)." (Ω VII.I, all three versions) | the mode's ONLY signature line: qualitative "far-horizon" stamp — **no detectability threshold, no rate normalization, no T₁/₂(g)**; classified with its host line, flagged in the note |

All 89 class-2 hits are instances/echoes of these seven items (the J3/i2/N2
analysis-layer restatements included, layer-tagged in the JSON).

### 4.2 Class 3 — occurrence stakes and kill language (95 hits; occurrence is not a rate number)

| item | verbatim quote (source) |
|---|---|
| S6 stake row | "S6 \| 0νββ occurs \| LEGEND-1000/nEXO era \| full-funnel exclusion (kills VII.I) \| conditional" (Ω registry, all three versions; Primer form "definitive full-funnel exclusion … LIVE, conditional on Ch. 14") |
| the funnel datum — **the near-miss** | "Data ⬛ for the row's frame: T₁/₂(¹³⁶Xe) > 3.8×10²⁶ yr (KamLAND-Zen), m_ββ < 28–122 meV; LEGEND-1000/nEXO reach the inverted-hierarchy funnel. Row status: **stake row — adjudication owned by S6; … nothing moved.**" (WS-K4:17); margin table "r3 0νββ \| stake row (S6) \| funnels: > 3.8×10²⁶ yr \| — \|" (WS-K4:41, margin column deliberately EMPTY) |
| the occurrence-only clause | "**any m_ββ value** (q-κ's justification: the corpus stakes *occurrence*, never a half-life)" (WS-nu-Q1:31) |
| the q-κ seal | "**(q-κ)** 0νββ half-life numerology sealed (no T½ arithmetic against corpus constants)" (WS-nu prompt pack:50) — the corpus structurally FORBIDS itself the very arithmetic a class-1 bound would require |
| four-cell table | branch-(b) cell: "line-contact; rate ∼10⁻¹⁸ below funnels → corpus-level contradiction / GUM *predicted* the null…" (WS-K5:40-43) — rate-CLASS language, no absolute number |
| adjudicating-data row | "0νββ funnels (LEGEND/nEXO era — S6's own adjudicators)" (WS-K-Ledger prompt pack:38) |
| gate return (corpus2 deltas) | "the 0νββ insertion rate is quadratic in the amplitude … R = **10¹⁴–10¹⁸** … The 'negligible' clause fails by 13–18 orders" (DELTA-NR-K2b-K3a:11; REVISION_MAP:1166) — carries the R statistic, NOT a bound |

### 4.3 The two auto-surfaced class-1 candidates, adjudicated by hand

1. **WS-K4:17 / :41 (the 3.8×10²⁶ yr datum)** — NOT class 1, on three independent
   grounds: **(a) wrong mode** — it is KamLAND-Zen's MASS-mechanism limit (0ν
   two-electron peak, quoted with its m_ββ conversion 28–122 meV); the Majoron mode
   (continuum spectrum) has separate experimental limits, none printed anywhere in the
   archive; **(b) no in-archive conversion, either direction** — the archive prints no
   T₁/₂(g) normalization (no nuclear matrix element, no phase-space factor) for any
   mode, so this half-life cannot be converted to a g-bound, and equally the corpus
   prints no predicted T₁/₂ at its g, so N2's honest rate R × printed cannot be
   converted to a half-life to compare against 3.8×10²⁶ yr; **(c) role in the print** —
   it is data context for the S6 OCCURRENCE stake ("stake row — adjudication owned by
   S6; nothing moved"; margin column printed empty), under the corpus's own q-κ seal.
2. **WS-nu prompt pack:50** — surfaces on `half-life` but is the OPPOSITE of a bound:
   the q-κ quarantine seal (quoted in §4.2).

### 4.4 Class 4 (186 hits)

Campaign-analysis-layer self-references (N2/J3/i2/roadmaps/DISCHARGE/tier files
discussing the absence of the bound — including every `g_lim` hit in the archive:
all 20 B1 hits live in N2, ROADMAP_v9, TIER7_ADJUDICATION and the tier6 checkpoint —
zero corpus-layer),
unrelated detectability language (the Unruh sector's "the medium is undetectable by
uniform motion", foliation detectability), signature-location lines without numbers,
and incidental numeric coincidences (C-family). Each entry carries its note in the JSON.

## 5. Decision-rule execution (G3)

**Class-(1): EMPTY → the ABSENT branch executes.** The F-T7-N2 conditional-kill
arithmetic is not executable in-archive: there is nothing to run it against. The kill
stays filed exactly as N2 printed it — *for ANY archive-external g_lim < 8×10⁻³ the
honest second moment converts P-ν4 from far-horizon signature to already-excluded
prediction at every anchor and every printed g* (min g_equiv = 8.2×10⁻³) — with its
arming condition now CERTIFIED rather than asserted: the certification is this memo's
§3 protocol, reproducible by `python3 o2_sweep.py`.

**Why class-(2) cannot adjudicate (the observable distinction, stated precisely).**
Every class-2 item binds the coupling g through emission, thermalization, free-streaming
or laboratory forces — channels whose amplitudes involve the FIRST moment of the web
mixing density (or no web structure at all). N2's enhancement R = ⟨A²⟩/⟨A⟩² multiplies
only the 0νββ-insert's spatial SECOND moment — the mode-rate observable. A g-side
constraint is therefore R-blind by construction: it cannot see, let alone exclude, a
rate enhancement at fixed g. (This is N2 §4's point, here verified against every
printed instance.)

**Why class-(3) cannot adjudicate.** The S6 apparatus stakes OCCURRENCE and is
adjudicated by full-funnel exclusion — a yes/no on the standard-mode search, not a rate
number; its one embedded half-life is the wrong mode with no printed conversion
(§4.3); and the corpus's own q-κ seal forbids T½ arithmetic against corpus constants.
Occurrence language, however emphatic, contains no number for R × printed to exceed.

**What would fire the kill.** A printed (archive-external) Majoron-mode 0νββ limit in
any of three equivalent currencies: g_lim directly; a Majoron-mode T₁/₂ limit PLUS the
mode's T₁/₂(g) normalization; or a quantitative detectability line for P-ν4. Any such
number below 8×10⁻³ — five-plus orders above the corpus's entire printed coupling
window [8×10⁻¹⁰, 1.3×10⁻⁸] — fires it at every anchor.

## 6. Honest limits / caveats

1. **In-archive certification only.** O2 certifies the bound is unprinted in this
   archive; it imports no external literature values (campaign rule), so it neither
   fires nor defuses the conditional kill — it certifies the kill's premise ("the
   archive prints NO such bound") at archaeology grade.
2. **Pattern-based coverage.** The sweep is lexical: a bound expressed with none of the
   14 patterns' tokens AND none of the mode tokens within ±5 lines would be missed.
   Mitigations: two independent families by construction (A hits the observable
   vocabulary, B the bound/experiment vocabulary, C the numerics), a candidate
   detector that fails loud, and the corpus's own structural evidence (q-κ seal;
   "occurrence, never a half-life") that no such print was ever intended.
3. **Context filter risk.** B2/B5/B6/C1 hits without a mode token in ±5 lines were not
   catalogued (raw counts logged). A Majoron-mode bound printed > 5 lines from any mode
   token is judged implausible given the corpus's formulaic row style, but is the
   filter's blind spot; the unfiltered patterns (A1–A5, B1, B3, B4, B4b, C2) carry no
   such filter and independently cover the bound vocabulary.
4. **Binary files skipped** (47: PNGs, PDFs, .npz archives — extension list in JSON).
   All are campaign-generated figures/data, not corpus documents.
5. **Classification is rule-based + hand-adjudicated overrides** (7 overrides, all
   printed in `o2_sweep.py` with reasons); the class-1/not-class-1 boundary — the only
   one the verdict depends on — was additionally protected by the fail-loud candidate
   detector, independent of the classifier.

## 7. Deliverables

- `o2_sweep.py` — the sweep driver (scope, patterns, filters, classifier, overrides,
  candidate detector; deterministic; re-run = `python3 o2_sweep.py`).
- `o2_results.json` — protocol block (file counts, pattern log, raw counts, filters,
  supplementary pass) + the complete 370-entry catalog (file, line, layer, patterns,
  class, note, verbatim untruncated quote) + the verdict block.
- This memo. **F-T7-O2 status: ABSENT certified — the archive prints no Majoron-mode
  0νββ rate-level bound; F-T7-N2's verdict (iii) and its conditional kill stand, the
  blocking number confirmed unprinted; adjudication requires an archive-external
  g_lim.**

— end of memo —
