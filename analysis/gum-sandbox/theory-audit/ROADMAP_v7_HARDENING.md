# GUM Program — Roadmap v7: Hardening & Frontier (post-Phase-J)

Successor to ROADMAP_v6 (executed: Phases I and J; ledger F-R1–F-R17;
T-ledger empty; no named opens). What remains is HARDENING: the three
micro-opens the J memos printed honestly, plus the strongest guarantee
an external replication package can carry — a fresh end-to-end re-run
of every command it promises.

## Phase K — four workstreams ▶ EXECUTING

**K1 — the /6 vs /(15π/8) discriminator (J1's remaining split).** The
edge-referenced amplitude convention is recovered AS A CLASS; the two
in-class candidates differ by 1.87% in 𝔭 (0.7σ against ±0.03 — archive
error bars alone cannot split them). Find and execute a discriminating
computation: an independent archive quantity sensitive to the amplitude
convention at better than ~2% (candidates: the I.2/G4 two-knot tail
amplitude, the B.4 band edges themselves, the I.1 gate constants). A
verified null — no such quantity exists in the archive — is a
legitimate terminal state and closes the item as "indistinguishable
within the archive"; print it as such. → k1_* artifacts.

**K2 — the ε < 0.005 probe (J4's stated limit).** Extend the J4 scan
to ε ∈ {0.001, 0.002} (two resolutions each, honest convergence gates;
if the solver cannot converge at a point, print the limit rather than
forcing it). Does the linear law hold to the smallest honestly
solvable ε? → k2_* artifacts extending j4_results.json's dataset.

**K3 — GPU pipeline extension (the "Beyond" on-ramp made concrete).**
Port the sector-measure + E_static gradient pipeline to WGSL f64 and
verify tolerance-band against the fs-gum-statics CPU goldens on
llvmpipe — the full readiness demonstration for real-GPU physics
(scheduled after K1/K2/K4 land; the heavy item).

**K4 — the REPRODUCE gauntlet (package freshness).** Execute every
command in DISCHARGE_PACKAGE/REPRODUCE.md fresh on this host, top to
bottom: record PASS/FAIL, runtime, and any drift from the printed
expectations (goldens must reproduce bit-identically within-backend;
tolerance-band items within band). Deliverable:
DISCHARGE_PACKAGE/FRESHNESS.md — the dated re-run record an external
replicator sees first. Any failure is a finding, not an embarrassment:
print it. → k4 artifacts.

## Beyond (unscheduled)
The corpus's response to corpus2/DISCHARGE_PACKAGE; real-GPU hardware.
