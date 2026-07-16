# Audit Logs — full agent/workflow execution transcripts

Complete chat + thinking + tool-call logs for every subagent and
workflow run of the GUM replication campaign and physics-core program,
preserved for scientific reproducibility and future audits.

## Contents
- `subagent-transcripts.tar.xz` — all `agent-<id>.jsonl` transcripts
  (41 at snapshot) + `agent-<id>.meta.json` + `workflows/wf_*/`
  (journal.jsonl with each workflow agent's structured return, plus
  per-agent transcripts and the executed workflow scripts).
- `coordinator-session.jsonl.xz` — the coordinating session's own
  transcript (the adjudication layer's full reasoning record).
- `MANIFEST.json` — agent-id → phase map (which transcript belongs to
  which tier/phase/finding), sizes, and ongoing-at-snapshot flags.

## Format
JSONL: one JSON object per line — user/assistant messages (assistant
entries include `thinking` blocks), tool_use calls with full inputs,
and tool_result contents. Read with any JSONL reader; each record
carries timestamps, sessionId, and gitBranch.

## Provenance chain
Findings F-R1–F-R8, the exact results (locked thresholds,
G* = 16√2/9, the oblate compacton), and every crate's gate table can
be traced from the committed RESULTS/ADJUDICATION files back through
these transcripts to the exact tool calls that produced them. The
deterministic artifacts themselves (gates binaries, golden Merkle
roots) re-verify independently of these logs; the logs additionally
preserve the reasoning, dead ends, spec errors, and corrections as
they happened (the print-your-own-defects record in its rawest form).

## Snapshot discipline
This archive is a live-session snapshot; agents flagged ONGOING in the
manifest (the N=192 crossing run, the theory-audit workflow) continue
to append. A final sweep re-archives everything at session end —
compare archive mtimes if in doubt.
