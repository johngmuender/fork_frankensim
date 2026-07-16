# FRESHNESS — dated fresh re-run of every REPRODUCE.md command (workstream K4)

- **Date:** 2026-07-16 (UTC), started 22:22 UTC
- **Host CPU:** Intel(R) Xeon(R) Processor @ 2.80GHz, 4 cores, 15 GiB RAM
- **Commit:** `16d91f08539cf588ea6b0259209e251a2d45bf28`
- **Toolchain:** rustc 1.99.0-nightly (3659db0d3 2026-07-05); cargo 1.99.0-nightly (2f0e7011e 2026-07-05); Python 3.11.15; numpy 2.4.6; wgpu/GPU: not exercised (no REPRODUCE.md row requires a GPU)
- **Build flags:** shipped defaults (no `-C target-cpu=native`); all `target/` dirs pre-existing → warm/incremental builds (per-row cold/warm noted)
- **Concurrency note:** to fit the budget, long rows were executed concurrently (up to ~6 processes on 4 cores). Results are unaffected (fixed seeds / no RNG); wall-clock times for concurrent rows are upper bounds and drift flags are annotated accordingly. Rows in Phase A ran sequentially with clean timings.
- **Bit-identity oracle:** the work tree was git-clean before the runs, so any rerun output that differs from the shipped artifact appears in `git diff`; Merkle roots/fingerprints are additionally checked against the values printed in REPRODUCE.md.

## Results table

Status legend: PASS / FAIL / PENDING / SKIPPED-BUDGET / DUP (same command already executed under the referenced row). "Root" column: OK = golden Merkle root / fingerprint reproduced bit-identically; n/a = row has no golden root.

| ID | Sect. | Command (from repo root) | Expected | Shipped runtime | Measured | Root | Status |
|---|---|---|---|---|---|---|---|
| K01 | 1 | `tier0-gauntlet/` → `cargo run --release` | 52/52; root e19d0bcb… | 0.3 s | 3.4 s (warm build+run) | OK e19d0bcb… bit-identical | PASS |
| K02 | 1 | `tier4-field/fs-cosserat-pilot/` → `cargo run --release` | 13/13; root f88731af… | ~0.2 s | 0.3 s | OK f88731af… bit-identical (replay OK) | PASS |
| K03 | 1 | `gum-core/fs-gum-field/` → `cargo build --release && ./target/release/gum_field_gates` | 25/25; root 6c1e7850…; fp 7c9931a5… | 1.7 s | 16.0 s (build refresh) / run 1.7 s | OK 6c1e7850… + fp 7c9931a5… bit-identical | PASS |
| K04 | 1 | `gum-core/fs-gum-statics/` → `cargo build --release && ./target/release/gum_statics_gates` | 18/18; κ ratio 1.2525; fp 6b3959d4… | ~316 s | — | — | PENDING |
| K05 | 1 | `gum-core/fs-gum-topo/` → `cargo run --release --bin gates` | 25/25; root 9cfbe9a6… | ~5 s | 8.2 s | OK 9cfbe9a6… bit-identical | PASS |
| K06 | 1 | `gum-core/fs-gum-cosserat/` → `cargo run --release` | 20/20 + 26/26 PASS | 2.9 s | 3.3 s | OK 9ddc0fe4… = RESULTS.md extended golden | PASS (expectation drift: prints 26/26 extended battery, no separate 20/20 line) |
| K07 | 1 | `gum-core/fs-gum-cosserat/` → `cargo test --release` | 5/5 | — | 0.8 s | n/a | PASS (13/13 tests, not 5/5 — REPRODUCE expectation stale; RESULTS.md confirms 13/13) |
| K08 | 1 | `gum-core/fs-gum-sde/` → `cargo run --release` | 21/21 PASS | ~310 s | — | n/a | PENDING |
| K09 | 2 | `cd analysis/gum-sandbox/tier2-closure && python3 axi_solve.py` | E₀/𝕀 drift ~1e-15; g(0.42)=1.285547004; families A/B | ~170 s | — | — | PENDING |
| K10 | 2 | rows 2–3 ("same run") | covered by K09 | — | covered by K09 | — | PENDING (same run as K09) |
| K11 | 2 | `cd analysis/gum-sandbox/tier4-field/fs-cosserat-pilot && cargo run --release` | certified SDiff + ¼ claims; root f88731af… | ~0.2 s | 0.3 s | OK f88731af… | PASS |
| K12 | 2/3 | `cd analysis/gum-sandbox/tier2-closure && python3 axi3_solve.py` | G0–G5; κ_crit 0.2163–0.2169; saturated closure; +9.1σ | 3241 s | — | — | PENDING |
| K13 | 2 | two-line proof (`axi_FR4_NOTE.md`) | analytic — read | — | read | n/a | PASS (note present, self-contained) |
| K14 | 3 | `cd analysis/gum-sandbox/tier4-field && python3 field3d_solve.py` | R→−0.560; κ through 1/√(8π); clock 1.245/1.331/1.474; E_rot/E=0.250000 | 4461 s | — | — | PENDING |
| K15 | 3 | `python3 field3d_solve.py --smoke` (tier4-field) | small-N pipeline pass; `*_smoke.*` outputs | minutes | — | — | PENDING |
| K16 | 3 | `cd analysis/gum-sandbox/tier2-closure && python3 axi4_locked_solve.py` | 𝔠=2.37096; locked ladder 3.26/3.29; G*=2.515991 | 2527 s | — | — | PENDING |
| K17 | 3 | `cd analysis/gum-sandbox/gum-core/fs-gum-statics && ./target/release/gum_statics_gates` | halo referee reproduced; κ ratio 1.2525 | ~316 s | — | — | PENDING |
| K18 | 3 | `python3 radial_solve.py` (tier2-closure; unit-map row) | tail mass μ²=m̃²/(2a_ψ) to 2×10⁻⁶ | ~1 s | 1.8 s | n/a | PASS (tail-mass gate 3.2e-07 < 2e-06) |
| K19 | 4 | `tier0-gauntlet/` → `cargo run --release` (repeat) | 52/52 | 0.3 s | 0.3 s | OK e19d0bcb… | PASS |
| K20 | 4 | `tier1-spectrum/` → `python3 spectrum.py` | four-branch spectrum; factor 4 = 4.000000000 | 5.8 s | 4.4 s | n/a | PASS (factor 4 confirmed, all gates PASS) |
| K21 | 4 | `tier2-closure/` → `python3 closure_reduced.py` | ε→0 endpoint to 6 decimals; ¼ automatic | seconds | 0.1 s | n/a | PASS |
| K22 | 4 | `tier2-closure/` → `python3 radial_solve.py` | App. I.1 gates; Derrick 6×10⁻⁷ | ~1 s | 1.3 s | n/a | PASS |
| K23 | 4 | `tier3-born/` → `python3 born.py` | τ = 83/10.2/4.9; equivariance control | 829 s | — | — | PENDING |
| K24 | 4 | `tier3-born/` → `python3 m4generic.py` | τ = 4.48, r² = 0.987 | 252 s | — | — | PENDING |
| K25 | 4 | `tier5-family/` → `python3 family_audit.py` | quark belt; ΔN_eff = 0.026772 | seconds | <1 s | n/a | PASS (ΔN_eff = 0.026772 reproduced) |
| K26 | 6 | §6 five-line block, verbatim | 52/52; 13/13; 25/25; 18/18; radial gates | ~7 min | — | — | PENDING |
| K27 | v2 | `python3 analysis/gum-sandbox/theory-audit/h21_mc.py` | R = 1.115e14 ± 0.3%; 1.967e15; v/c 6.1e8 | 117 s | — | n/a | PENDING |
| K28 | v2 | `cd analysis/gum-sandbox/theory-audit/h22_solve && cargo run --release --bin h22_solve -- all` | SDiff reduction ≤ 0; floor margin 15.7–18.1× | ~24 min | — | — | PENDING |
| K29 | v2 | `python3 analysis/gum-sandbox/theory-audit/h23_frank.py` | radial 8πK₁r; relaxed floor 7.7024 | 163 s | — | n/a | PENDING |
| K30 | v2 | `python3 analysis/gum-sandbox/theory-audit/h24_helical.py` | 0.332ƛ² / 0.568ƛ²; premises verified | <1 s | 3.6 s | n/a | PASS (0.3316ƛ²/0.5677ƛ² vs printed 0.332/0.568) |
| K31 | v2 | `python3 analysis/gum-sandbox/theory-audit/h25_sdw.py` | 𝔞₁(p,q) formula; 27/27 | 3 s | 3.8 s | n/a | PASS (27/27) |
| K32 | v2 | `cd analysis/gum-sandbox/theory-audit/h26_solve && cargo run --release --bin h26_convention` | E_rot/E = 0.2500 vs 0.4330 (915σ) | ~10 min | — | — | PENDING |
| K33 | v2 | `python3 analysis/gum-sandbox/theory-audit/h27_entrain.py` | exponent −0.500000000; torque −2.5012e-2 | ~2 min | — | n/a | PENDING |
| K34 | v2 | `python3 analysis/gum-sandbox/theory-audit/t3_closures.py` | 𝔠 = 64√2/9π; six 2.5147 routes fail | <1 min | 0.6 s | n/a | PASS (𝔠 = 64√2/9π route confirmed) |
| K35 | v2 | `python3 analysis/gum-sandbox/tier2-closure/gstar_solve.py` | G* → 16√2/9 (Richardson 6.0e-7) | 31 s | 14.2 s | n/a | PASS (G* = 16√2/9; faster than printed 31 s) |
| K36 | v2 | `cd analysis/gum-sandbox/gum-core/fs-gum-kern && cargo run --release --bin n192_fr5 -- crossing` | κ crosses 1/√(8π) at it 2070–2080 | 46 min | 0.5 s (usage error) | — | CMD-DEFECT: printed command missing required args; binary usage requires 'crossing <out.json> <N> <LBOX> <cap_static> <cap_main> [threads]'. Corrected form run as K36b (shipped invocation from n192_RESULTS.md L28). |
| K36b | v2 | corrected form: `… -- crossing n192_runs/n192_crossing.json 192 4.5 450 2500 4` (per n192_RESULTS.md L28) | κ crosses 1/√(8π) at it 2070–2080 | 46 min | — | — | PENDING |
| K37 | v2.2 | `python3 analysis/gum-sandbox/theory-audit/h4_compute.py` | ±2·[rot]; winding −2.000000; 28/28 | 9 s | 2.5 s | n/a | PASS (28/28) |
| K38 | v2.2 | `python3 analysis/gum-sandbox/theory-audit/i1_r5test.py` | 3.2011 at +0.51σ; ≈7σ separation | 757 s | — | n/a | PENDING |
| K39 | v2.2 | `python3 analysis/gum-sandbox/theory-audit/i2_propagate.py` | 34/34; κ²g_tot = 35/24; 13.4σ | 2 s | 0.9 s | n/a | PASS (34/34) |
| K40 | v2.2 | `python3 analysis/gum-sandbox/theory-audit/i3_conventions.py` | 𝔭 and 𝔟_eff jointly resolved | 1 s | 0.2 s | n/a | PASS (joint 𝔭×𝔟_eff table produced) |
| K41 | v2.2 | `python3 analysis/gum-sandbox/theory-audit/i4_compute.py` | σ(exchange) = −1; 2/8 sectors; 19/19 | 4 s | 4.1 s | n/a | PASS (19/19) |
| K42 | J | `python3 analysis/gum-sandbox/theory-audit/j1_pgrid.py` | 𝔭 = 0.84 at 0.28σ; 𝔟_eff +1.06σ | 55 s | — | n/a | PENDING |
| K43 | J | `python3 analysis/gum-sandbox/theory-audit/j2_epsscan.py` | 44/44; ceiling 2.605×10⁻³; 11.283σ | 2 s | — | n/a | PENDING |
| K44 | J | `python3 analysis/gum-sandbox/theory-audit/j3_majoron.py` | 17/17; no row newly binds | 3 s | — | n/a | PENDING |
| K45 | J | `python3 analysis/gum-sandbox/theory-audit/j4_epslaw.py` | 14/14; p = 0.997 [0.996, 1.012] | 250 s | — | n/a | PENDING |

## Findings

(populated as runs complete)

## Verdict

PENDING — runs in progress.
