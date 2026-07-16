# REPRODUCE — every load-bearing number, its command, and its guarantee

All paths are relative to the frankensim repository root. The campaign's
artifacts (scripts, JSON results, run logs, figures, adjudications) live
where they were produced, under `analysis/gum-sandbox/`; nothing in this
package is a copy. Everything below is deterministic: fixed seeds or no
RNG at all; the Rust gate suites additionally run themselves twice and
assert a bit-identical Merkle root (replay), so a changed 15th digit is a
loud failure, not drift.

## 0. Environment

- **Python 3** with `numpy` + `matplotlib` only (no scipy anywhere;
  verified against numpy 2.4.6). Single machine, 4 cores; runtimes below
  are from that box.
- **Rust**: the repository's pinned nightly (verified with
  `cargo 1.99.0-nightly (2026-07-05)`). Every crate is a standalone
  `[workspace]` opt-out with path dependencies on the repo's `crates/fs-*`
  (std-only closure) — so a full frankensim checkout is required, but no
  external crates are downloaded and no network is touched.
- Verification status column: **[re-run 2026-07-12]** = executed while
  assembling this package, output matched the recorded golden values;
  **[cite log]** = hour-class solve, not re-run here — its complete
  printed output, JSON, and figures are in the repository and the command
  is verified to exist and parse.

## 1. Certified gate suites (Rust — the fast, fail-closed spine)

| Suite | Command (from the listed directory) | Expected output | Runtime | Status |
|---|---|---|---|---|
| Tier-0 constants gauntlet, **52/52** | `analysis/gum-sandbox/tier0-gauntlet/` → `cargo run --release` | `== 52/52 checks PASS ==`; Merkle root `e19d0bcb3e97e430d9cca8d0823705f7b6b0807bf3601d5f40923c5f010a85f3`; checker modes false/true/false | 0.3 s (binary), ~1 min cold build | **[re-run 2026-07-12]** — root matched |
| Tier-4B certified 3-D diagnostics, **13/13** | `analysis/gum-sandbox/tier4-field/fs-cosserat-pilot/` → `cargo run --release` | 13/13 PASS; golden root `f88731afb676d76570ac719fc5b79ced053a1fa95e17750510e375debc64fd7a`; in-process replay | ~0.2 s | [cite log — `RESULTS.md`/`RESULTS.txt`] |
| fs-gum-field, **25/25** (stored field, sector energies, schemes) | `analysis/gum-sandbox/gum-core/fs-gum-field/` → `cargo build --release && ./target/release/gum_field_gates` | `OVERALL: PASS`; root `6c1e785071a753e40a33c4d8532886d2ca20bd245d7f184f67b5eed9c6c9fb1f`; fingerprint `7c9931a5…`; replay identical | 1.7 s | **[re-run 2026-07-12]** — root + fingerprint matched |
| **fs-gum-statics, 18/18 — the F-R5 halo referee** | `analysis/gum-sandbox/gum-core/fs-gum-statics/` → `cargo build --release && ./target/release/gum_statics_gates` | 18/18 PASS incl. gate G-H (halo referee at N = 48: threshold crossing, halo growth, differential control) and G-K (clock ratio κ(L_clock)/threshold = 1.2525 vs Python 1.245); FD-vs-analytic gradients ≤ 1.2×10⁻⁶; 75 gate numbers bit-identical across the two included runs | ~5 min | **[re-run 2026-07-12]** — see note below |
| fs-gum-topo, 25/25 (Hopf two ways, Z/2 suite) | `analysis/gum-sandbox/gum-core/fs-gum-topo/` → `cargo run --release --bin gates` | 25/25 PASS; root `9cfbe9a6c42a41e1a891413c6763fca1f0d5be656a174c90a25e32f1744b89c1` | ~5 s | [cite log] |
| fs-gum-cosserat, 20/20 + Wχ 26/26 | `analysis/gum-sandbox/gum-core/fs-gum-cosserat/` → `cargo run --release` (also `cargo test --release`, 5/5) | 20/20 + 26/26 PASS; golden root in `RESULTS.md` | 2.9 s | [cite log] |
| fs-gum-sde, 21/21 (incl. Nelson-vs-Bohm) | `analysis/gum-sandbox/gum-core/fs-gum-sde/` → `cargo run --release` | 21/21 PASS | ~310 s | [cite log] |

Every suite bundles its claims as fail-closed `Certified<f64>` into an
EvidencePackage and checks it three ways (deny-all must refuse;
certificate capability must pass; tampered root must fail).

**Note on the statics re-run:** executed during package assembly
(316.1 s, both internal runs): `OVERALL: PASS`, clock ratio
κ(L_clock)/threshold = 1.2525 as recorded, run fingerprint (BLAKE3 over
every gate f64)
`6b3959d4feb3573323a165aee43543ee09dddf3ac6c6e293ea6252480ad5d280`.
If your build box differs, the two-run internal replay is the guarantee
that matters — a PASS on your machine certifies the same claim contents.

## 2. F-R4 — commands behind each load-bearing number

| Number (as used in README §5) | Source artifact | Command | Runtime |
|---|---|---|---|
| SDiff invariance measured: E₀ drift 1.6×10⁻¹⁵, 𝕀 drift 1.2×10⁻¹⁵ at λ = 0.6 (`V3_cross`) | `tier2-closure/axi_results.json`, `axi_run.log` | `cd analysis/gum-sandbox/tier2-closure && python3 axi_solve.py` | ~170 s, no RNG |
| g(0.42) = 1.285547004 vs paper formula 1.285547005 (9 digits); scaling rung 2√(ê₀𝔦₀) = 2.0532877 to 7 digits | same run (gates V4, endA_analytic) | same | — |
| Family A: g ≡ 1, ε→0 at 𝔠 = 2.0533, κ = 0.9354, V = √2; Family B: 𝔠 = 2.279, κ = 0.907 at ε = 0.05 (⟨r1⟩ κ +5.8σ) | same run | same | — |
| SDiff invariance + ¼ invariant as **certified** interval claims (independent Rust engine) | `tier4-field/fs-cosserat-pilot/RESULTS.md` | `cd analysis/gum-sandbox/tier4-field/fs-cosserat-pilot && cargo run --release` | ~0.2 s |
| Compacton marginal spectrum: 12 modes, all dE/dp = 0 (FD floor), all d²E/dp² > 0, all dI/dp ≠ 0 | `tier2-closure/axi3_results.json` §R3, `axi3_run.log` | `cd analysis/gum-sandbox/tier2-closure && python3 axi3_solve.py` | 3241 s, no RNG — [cite log] |
| The two-line proof itself | `tier2-closure/axi_FR4_NOTE.md` | (analytic — read it) | — |

## 3. F-R5 — commands behind each load-bearing number

The single command `python3 axi3_solve.py` (from
`analysis/gum-sandbox/tier2-closure/`; 3240.9 s wall, deterministic, no
RNG; writes `axi3_results.json`, `axi3_run.log`, `axi3_modes.png`)
produces: gates G0–G5 (incl. G0b Bogomolny normalization to 1.1×10⁻¹⁶
and G5 reproducing Step 2's closure to 8×10⁻⁶); the threshold derivation
check (shell probes κ_crit 0.2163–0.2169 vs ideal 1/√(8π) = 0.19947, dI
ratio 4.000); the placement of every Step-2/corpus tuple above threshold;
the direct closures with their `!! BOUNDARY` records; the saturated
closure (κ_paper = 1/√2 exact, internal identity to 2.2×10⁻¹⁶,
G ≥ π/√2 verified to 2×10⁻¹⁶, bracket 𝔠 ∈ [2.82843, 3.20146]); and the
R2 pulls against the ⟨r1⟩ benchmark (+9.1σ direct). **[cite log:
`axi3_run.log` holds every printed table.]**

| Confirmation layer | Command | Key expected outputs | Runtime |
|---|---|---|---|
| **3-D, no symmetry (Tier 4A)** | `cd analysis/gum-sandbox/tier4-field && python3 field3d_solve.py` (fixed seeds 424242/11/12) | R monotone to −0.560 below the axisymmetric reference; κ through 1/√(8π); halo box-limited (1.3× box deeper); 88% of dI priced at 1/(16π); control self-limits at (3/2)·I; clock block: κ(L_clock) = 1.245×/1.331×/1.474× threshold; E_rot/E = 0.250000 at L_clock | 4461 s — [cite log: `field3d_run.log`, `field3d_results.json`] |
| — pipeline smoke test (fast) | `python3 field3d_solve.py --smoke` | small-N pipeline pass, separate `*_smoke.*` outputs | minutes |
| **Direction-locked closure (Tier 4C)** | `cd analysis/gum-sandbox/tier2-closure && python3 axi4_locked_solve.py` | the locked onset family 1/√(2⟨sin²θ⟩_w) with all probe rows; the exact identities: polar channel = 1.00000, s = sin² channel = √(7/12) (⟨sin²θ⟩_w = 6/7), uniform = √3/2; uniform-saturated hedgehog-core closure 𝔠 = **2.37096**; locked cap-ladder runaway to 𝔠 = 3.26/3.29; ring-saturated G\* = 2.515991 | 2527 s, no RNG — [cite log: `axi4_locked_run.log`] |
| **Rust halo referee (independent port)** | `cd analysis/gum-sandbox/gum-core/fs-gum-statics && ./target/release/gum_statics_gates` | halo threshold crossing + growth + differential control reproduced at N = 48; clock ratio 1.2525 | ~5 min — **[re-run 2026-07-12]** |
| Unit-map validation (the √2 the convention question turns on) | `python3 axi_solve.py` gates V2/V4 + `python3 radial_solve.py` (tail mass μ² = m̃²/(2a_ψ) to 2×10⁻⁶) | 7–9 digit agreement on three independent paper quantities | 170 s / ~1 s |

## 4. The replicated-successes side (context numbers)

| Result | Command | Runtime | Status |
|---|---|---|---|
| 52/52 archive arithmetic incl. ⟨r10⟩ bridge (m₃ = 0.0468 eV), F-R1/F-R3 exhibits, WS-A4, NR-D1b | `analysis/gum-sandbox/tier0-gauntlet/` → `cargo run --release` | 0.3 s | **[re-run 2026-07-12]** |
| Four-branch spectrum; factor 4 = 4.000000000; Theorem II.2; M-1 mechanism | `analysis/gum-sandbox/tier1-spectrum/` → `python3 spectrum.py` (seed 20260711, no randomness consumed) | 5.8 s | [cite log] |
| Reduced closure: ε→0 endpoint to 6 decimals; ¼ automatic; F-R2 exhibit | `analysis/gum-sandbox/tier2-closure/` → `python3 closure_reduced.py` | seconds | [cite log: `run.log`] |
| Radial rung: all App. I.1 gates; Derrick 6×10⁻⁷; blind tail-mass confirmation | same dir → `python3 radial_solve.py` | ~1 s | [cite log] |
| Born relaxation: τ = 83/10.2/4.9; equivariance control | `analysis/gum-sandbox/tier3-born/` → `python3 born.py` (seeds 42/12345, N = 20,000) | 829 s | [cite log] |
| Generic-set M = 4 rerun (τ = 4.48, r² = 0.987 — corpus's "M ≳ 4" vindicated) | same dir → `python3 m4generic.py` | 252 s | [cite log] |
| Tier 5a: quark-belt reconstruction, lepton pulls, ΔN_eff = 0.026772; F-R6/F-R7/F-R8 exhibits | `analysis/gum-sandbox/tier5-family/` → `python3 family_audit.py` | seconds | [cite log: `family_results.json`] |

## 5. Determinism guarantees, stated precisely

1. **No RNG at all**: `radial_solve.py`, `axi_solve.py`, `axi3_solve.py`,
   `axi4_locked_solve.py`, `closure_reduced.py`, `family_audit.py`
   (hand-rolled deterministic Nelder–Mead where optimization is needed;
   fixed quadrature domains so objectives are smooth deterministic
   functions of parameters). Reruns reproduce the logs byte-for-byte
   except timings.
2. **Fixed seeds, consumed identically**: `spectrum.py` (20260711),
   `born.py`/`m4generic.py` (42/12345), `field3d_solve.py`
   (424242/11/12 — and the halo seed is applied identically to over-spun
   and control runs, so the two-sided verdict does not depend on it).
   `fs-gum-statics` uses a hand-rolled 64-bit LCG (constants printed in
   its RESULTS) — no platform RNG anywhere.
3. **Bit-replay certification**: every Rust suite runs its full pipeline
   twice in-process and asserts identical BLAKE3 Merkle roots; the
   recorded golden roots above are the cross-machine handshake. Tier-0
   and fs-gum-field were re-executed while assembling this package and
   matched their recorded roots exactly.
4. **Convergence evidence is in-band**: every solver report carries its
   own grid ladder (N = 360 vs 720 re-evaluations; N = 48/64/96 gate
   tables), boundary/cap flags (`interior = NO`, walls-hit columns), and
   clock residuals — the error budget travels with the number, not with
   our say-so.

## 6. What you would run first (suggested 15-minute path)

```
cd analysis/gum-sandbox/tier0-gauntlet        && cargo run --release      # 52/52, root e19d0bcb…
cd ../tier4-field/fs-cosserat-pilot           && cargo run --release      # 13/13, root f88731af…
cd ../../gum-core/fs-gum-field                && cargo build --release && ./target/release/gum_field_gates    # 25/25, root 6c1e7850…
cd ../fs-gum-statics                          && cargo build --release && ./target/release/gum_statics_gates  # 18/18, halo referee (~5 min)
cd ../../tier2-closure                        && python3 radial_solve.py  # ~1 s, App. I.1 gates
```

Then the decisive solves, in increasing cost: `axi_solve.py` (170 s,
F-R4's measurement), `axi4_locked_solve.py` (42 min, the exact
identities), `axi3_solve.py` (54 min, F-R5's derivation + saturation),
`field3d_solve.py` (74 min, the 3-D confirmation). Each writes a JSON
whose fields are named in the corresponding `*_RESULTS.md`.


## v2 additions (theory-audit round, 2026-07-16)

All deterministic; seeds printed in each script. Paths relative to repo
root.

| number | command | runtime | expected |
|---|---|---|---|
| F-R9 R-table | `python3 analysis/gum-sandbox/theory-audit/h21_mc.py` | 117 s | R = 1.115e14 ± 0.3% at (a=52.6 fm, L=1 µm); 1.967e15 at 4.2 µm; rescue v/c 6.1e8 at r0=2 fm |
| F-R4 floor on the collision problem | `cd analysis/gum-sandbox/theory-audit/h22_solve && cargo run --release --bin h22_solve -- all` | ~24 min | SDiff reduction beyond leak ≤ 0; floor margin 15.7–18.1× |
| F-R11 honest floor | `python3 analysis/gum-sandbox/theory-audit/h23_frank.py` | 163 s | radial 8πK₁r; hyperbolic (8π/15)(3K₁+2K₃)r; relaxed floor 7.7024 |
| F-R10 counterexamples | `python3 analysis/gum-sandbox/theory-audit/h24_helical.py` | <1 s | ⟨r_⊥²⟩ = 0.332ƛ² (annulus), 0.568ƛ² (torus knot); all premises verified |
| T-H10 master form | `python3 analysis/gum-sandbox/theory-audit/h25_sdw.py` | 3 s | 𝔞₁(p,q) = −(5p+q)/(12(2p+q)); 27/27 internal checks |
| anchor discriminator | `cd analysis/gum-sandbox/theory-audit/h26_solve && cargo run --release --bin h26_convention` | ~10 min | E_rot/E = 0.2500 vs 0.4330 (915σ) |
| entrainment/torque | `python3 analysis/gum-sandbox/theory-audit/h27_entrain.py` | ~2 min | phase exponent −0.500000000; torque −2.5012e-2 vs analytic −2.5e-2 |
| repaired closure | `python3 analysis/gum-sandbox/theory-audit/t3_closures.py` | <1 min | 𝔠 = 64√2/9π; T3.1 KKT-on-boundary; six 2.5147 routes fail |
| G* exact | `python3 analysis/gum-sandbox/tier2-closure/gstar_solve.py` | 31 s | G* → 16√2/9 (Richardson 6.0e-7) |
| N=192 crossing | `cd analysis/gum-sandbox/gum-core/fs-gum-kern && cargo run --release --bin n192_fr5 -- crossing` | 46 min | κ crosses 1/√(8π) at it 2070–2080 |


## v2.2 additions (synthesis round / Phase I, 2026-07-16)

All deterministic (no RNG); runtimes re-verified on the shipping host.

| number | command | runtime | expected |
|---|---|---|---|
| T-H4 discharge / F-R16 | `python3 analysis/gum-sandbox/theory-audit/h4_compute.py` | 9 s | ∂(gen) = ±2·[rot] on both routes; clutching winding −2.000000; 28/28 |
| I1 ⟨r5⟩ corroboration test | `python3 analysis/gum-sandbox/theory-audit/i1_r5test.py` | 757 s | saturated anchor 64√2/9π = 3.2011 at +0.51σ; branch separation ≈7σ; both analytic floors respected at every ε |
| I2 propagation memo checks | `python3 analysis/gum-sandbox/theory-audit/i2_propagate.py` | 2 s | 34/34; S4′ → κ²g_tot = 35/24, ⟨r1⟩ pull 13.4σ (vs 0.7σ restricted); ⟨r10⟩ bridge ≤ 0.12σ |
| I3 joint convention solve | `python3 analysis/gum-sandbox/theory-audit/i3_conventions.py` | 1 s | one edge-referencing element resolves 𝔭 and 𝔟_eff jointly (B.4/I.2 consistent) |
| I4 exchange lift / F-R17 | `python3 analysis/gum-sandbox/theory-audit/i4_compute.py` | 4 s | σ(exchange) = −1 ≠ +1 = σ(rotation); B=2 sectors 2/8; 19/19 |

### Phase-J additions (the executed obligations, 2026-07-16)

| number | command | runtime | expected |
|---|---|---|---|
| J1 B.4-window 𝔭-grid | `python3 analysis/gum-sandbox/theory-audit/j1_pgrid.py` | 55 s | all solver gates PASS; 𝔭(edge-/6) = 0.84 at benchmark 0.28σ; 𝔟_eff +1.06σ |
| J2 (4.8′)/(4.9′) | `python3 analysis/gum-sandbox/theory-audit/j2_epsscan.py` | 2 s | 44/44; ceiling 2.605×10⁻³; anchor forfeit 11.283σ; dichotomy survives |
| J3 Majoron battery | `python3 analysis/gum-sandbox/theory-audit/j3_majoron.py` | 3 s | 17/17; old-window control exact; no row newly binds |
