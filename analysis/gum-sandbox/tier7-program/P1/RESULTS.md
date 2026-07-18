# P1 — corpus3 v4.2 delta: the Phase-N/O fold (F-T7-P1 editor record)

**Workstream:** Tier-7 Phase P, P1 (document editor, REVISION_CHARTER_v2
rules; marker ⟦T7⟧ retained). **Date:** 2026-07-18. **Payload authority:**
corpus3/T7_FOLD_PAYLOAD.md, ADDENDUM — Phase N/O payload (Q6–Q9) — every
folded number copied from the payload, none restated from memory.
**Within-model; nothing here bears on nature.**

## What was done

1. **Rename + version headers.** `git mv corpus3/01-GUM-Omega-Paper-v4.1-ext.md
   corpus3/01-GUM-Omega-Paper-v4.2-ext.md`; all **six** version-header
   lines bumped v4.1-ext → v4.2-ext (main title + five PART-N running
   heads; the "(Tier-6/7 fold)" parenthetical kept) — the sole
   permitted deletions under the charter's version-header exception.
2. **Seven new ⟦T7⟧ blocks inserted** (additive; corpus prose style;
   payload-sourced; artifact-cited):
   - **Head matter** — ⟦T7⟧ REVISION NOTE (v4.2-ext) directly after the
     v4.1 note: the Phase-N/O fold (F-T7-N2/N3/O1/O2/O3); charter +
     payload-addendum pointers; adoption is the authors'.
   - **Sec. 0** — dated fourth/fifth-phase public record [Q9].
   - **Sec. VI.C** (after the F-R14 dispute box) — the Q8 block: the
     computed decision matrix, ~104 cells 21/21 checks, exactly ONE
     sound cell (P-acoustic + supertrace), all five totals at
     N_knot = 3, tuning fraction 1/4, the minimum bill computed;
     decision support only.
   - **Sec. VII.I** (after the Majoron dictionary / P-ν4 paragraph) —
     the Q6 block: the honest second-moment rate (R × claimed;
     R = 1.1×10¹⁴–8.9×10¹⁷), ABSENT-CERTIFIED archaeology (658 files,
     zero rate-level bounds), the armed conditional kill
     (g_lim < 8×10⁻³), the F-R9 terminal menu; S1/S6 untouched.
   - **Sec. VIII.F** (after the ⟦T6⟧ proposal #5 box) — the Q7 block:
     kill-bin zeros at field/flow grade (O1: 64,000 points, zero
     violations, margin 0.621 = 22,068× the p95 path error, measure
     ≤ 2.8×10⁻⁷, ~5,400× sharper than L7b's 1.5×10⁻³), the N3
     line-criterion refutation, T4-W5 residual narrowed to the
     certified/analytic form.
   - **Sec. X** — Phase-N/O campaign record [Q9].
   - **Sec. XI** — register update: T4-W5 residual narrowed (line
     route closed by proof; field/flow grade attained;
     certified/analytic form the terminal open); the F-R9 terminal
     menu + the armed external kill; F-R14 carries the computed
     matrix pointer.
3. **README.manifest-v5.0-ext.md** — census filename cell
   v4.1 → v4.2 (the noted version-header-class exception under S-30
   rule 2, matching the v4.1 precedent), a ⟦T7⟧ v4.2 census-delta
   paragraph, and a ⟦T7⟧ version-history continuation line.
4. **09-External-Replication-Record-v2.0.md** — ⟦T7⟧ §9 appended in
   the record's format: scope/method + verdict census, the
   five-workstream credit table (N2/N3/O1/O2/O3), the findings'
   standing, and a §9-binding epistemic notice.

## Verification (gates)

| Gate | Wording (Phase-P addendum) | Measured | Verdict |
|---|---|---|---|
| G1 | Additive-only diff vs git HEAD's v4.1-ext content | `diff <(git show HEAD:…v4.1-ext.md) …v4.2-ext.md`: **6 deleted lines, all six version headers** (main title + 5 PART heads); 20 added lines (6 bumped headers + 14 lines = 7 ⟦T7⟧ blocks with separators). Manifest: 1 deleted line = the census filename cell (the noted exception); Replication Record: **0 deletions**. | PASS |
| G2 | Marker preservation + counts | Paper occurrence counts before→after: ⟦T7⟧ 11→23 (+12, all inside the 7 new blocks), ⟦T6⟧ 31→32 (+1, a mention inside the new head note), ⟦rev⟧ 170→171 (+1, ditto), ⟦H4 resolved⟧ 7→7, ⟦I4 reported⟧ 5→5. Leading-marker block count 8→**15** (the eight prior ⟦T7⟧ blocks verbatim — proven by the G1 diff: no non-header line touched). | PASS |
| G3 | Payload fidelity, ≥ 15 Q6–Q9 numbers verbatim | **39/39 spot-checked strings verbatim** in the paper (34 exact matches against the payload text; 5 identical after payload line-wrap normalization): R range, RSD, Lorenz shares, g_equiv, orders, f_req, KamLAND-Zen numbers, g_lim kill, axial support, covering, 64,000/80,000/0.621/22,068×/2.8×10⁻⁷/5,400×/74.0%, all five O3 totals, 21/21, 1/18, 5/11, tuning 1/4, workflow IDs, three verdicts. | PASS |
| G4 | Cross-file consistency | Census cell, reading-order pointer, and version-history line all name **01-GUM-Omega-Paper-v4.2-ext.md**; Record §9 and the paper's Sec. 0/X blocks carry identical verdict censuses, workflow IDs, and artifact paths; every new block cites its finding (F-T7-N2/N3/O1/O2/O3) and tier7-program/ artifact; the same numbers appear payload-identically in paper and Record §9. | PASS |

## Key numbers (of record, from the payload addendum)

- Q6: R = 1.1×10¹⁴–8.9×10¹⁷; g_equiv = √R·g = 8.5×10⁻³–12.3
  (6.9–10.1 orders above P-ν4); f_req = 490–44,000 TeV; O2 census
  658 files / 14 patterns / 370 hits, zero Majoron-mode rate-level
  bounds; armed kill g_lim < 8×10⁻³.
- Q7: 64,000 detector points, zero violations, 80,000/80,000
  agreement; min margin 0.621 = 22,068× p95 path error; measure
  ≤ 2.8×10⁻⁷ per kill member (~5,400× sharper than 1.5×10⁻³).
- Q8: one sound cell (P-acoustic + supertrace); totals at N_knot = 3
  (naive/supertrace): frame −5/3 / +1/3; covector −7/6 / +5/6; vector
  −119/66 / +13/66; cone-only −17/12 / +7/12; P-acoustic(w = −1/4)
  −1/6 / +11/6; tuning fraction 1/4.
- Edit census: paper 6 header bumps + 7 new ⟦T7⟧ blocks (8→15);
  manifest +2 ⟦T7⟧ paragraphs + 1 cell; Record +1 section (§9).

## Method

Charter v2 rules 1–7 applied with marker ⟦T7⟧ (the Phase-N/O pass);
single-source payload discipline (rule 7); additive-only editing with
the version-header exception; verification by (i) full diff against
`git show HEAD:analysis/gum-sandbox/corpus3/01-GUM-Omega-Paper-v4.1-ext.md`,
(ii) scripted occurrence counts of all five marker families, and
(iii) a scripted 39-string verbatim comparison of paper/Record text
against the payload addendum.

## Caveats

- The five G3 strings not found by exact substring match in the payload
  are hard-wrapped there; after whitespace normalization all 39 are
  byte-identical — no number was restated from memory.
- The manifest's census-row filename-cell edit is the one non-additive
  touch outside the paper's headers; it follows the v4.1 fold's own
  precedent and is recorded as the header-class exception both here and
  in the manifest's ⟦T7⟧ census-delta paragraph.
- The v4.1-ext text (including its eight ⟦T7⟧ blocks and its REVISION
  NOTE, which still self-describes as v4.1) is carried verbatim as a
  dated record — per S-30, prior notes are never reworded.
- Adoption of every folded block is the authors'; grades never rise by
  replication; offers remain offers.

**Within-model; nothing here bears on nature.**
