# Simulating GUM — Updated Assessment (v2)

**Provenance.** This document is a re-review, performed under the Fable 5 model, of the full context assembled this session: the GUM Session Map (`08-GUM-Session-Map-v1.0.md`), the theory of record (`01-GUM-Omega-Paper-v2.0.1.md`, read end-to-end including AUD-15), the six evidence-grounded frankensim cluster audits, and the two prior assessments written under Opus 4.8 (`GUM_SIMULATION_FEASIBILITY.md`, `EXTERNAL_BASE_CODE_ASSESSMENT.md`). Per the corpus's own annotate-never-rewrite practice, the prior documents stand unedited; this v2 states **what stands, what is corrected, and what is new**.

---

## 0. Delta log — what changed under re-review

| # | Prior position (Opus 4.8) | Updated position (Fable 5) | Kind |
|---|---------------------------|----------------------------|------|
| Δ1 | Goal framed as "simulate the GUM theory"; build the 3-D SU(2) field type first (Phase 0) | Goal reframed as **adversarial replication of the paper's own tagged computations ⟨r1⟩–⟨r11⟩** — the external replication the corpus's own SWP Branch S mandates. The 3-D field is a *late* tier, not the entry point | Reframe |
| Δ2 | Effort dominated by new 3-D field machinery | **Dimensionality inversion:** nearly all of the paper's load-bearing numbers live in 0–2 dimensions (algebraic identities, radial ODEs, axisymmetric 2-D solves, 2-D ensembles). Tiers 0–3 below need no 3-D field, no FEM, no GPU — and run almost entirely on existing frankensim crates | New |
| Δ3 | Stochastic-QM sector (R9) graded "build new… a separate research program" (SDE + Fisher engine) | **Downgraded in difficulty.** The paper's actual numerical claim there (⟨V-SIM-4⟩, Born-rule relaxation) is **deterministic de Broglie–Bohm dynamics — no SDE required**: split-step Schrödinger (fs-fft) + trajectory integration (fs-time) + a ~50-line H̄ functional. The Nelson-diffusion variant remains optional extra work | Correction |
| Δ4 | mumax3 caveat: "dissipative LLG… not symplectic" | Sharpened to the correct structural statement: **LLG is first-order gyroscopic (spins have no inertia); GUM's orientation sector is second-order inertial (`Jφ̈ = …`)**. mumax3 is excellent for *statics* (soliton profiles, DMI chirality, skyrmion/Hopf charges) and for cross-validating Tier-4 diagnostics; its *dynamics* is the wrong equation class regardless of damping | Correction |
| Δ5 | External-code recommendation: a per-sector polyglot pipeline (mumax3 + stark + FEniCS + frankensim) | Softened. Under the replication-first framing, **Tiers 0–3 are essentially pure frankensim** (fs-cheb, fs-la, fs-fft, fs-ascent, fs-adjoint, fs-time, fs-ivl, fs-evidence). External engines matter only from Tier 4. And for the knot/closure core, the closest external lineage is not micromagnetics but the **BPS-Skyrme numerical literature** (Adam–Sánchez-Guillén–Wereszczyński — the paper's own W₆₊₀ *is* the BPS-Skyrme model) and cuSkyrmion-class arrested Newton flow | Revision |
| Δ6 | Numerical hazards not itemized | New section on concrete hazards, the sharpest being the **compacton-edge log divergence** (verified numerically below): the ε>0 boundary layer the paper itself flags at O(ε ln 1/ε) is the stiff spot of the whole Tier-2 program | New |
| Δ7 | — | **Governance homology observed:** frankensim's evidence colors (Verified/Validated/Estimated with anti-laundering composition) and GUM's ledger classes ([DF]/[DW]/[CAL]/[IM]/[CJ] with language locks and grade demotion) are structurally homologous claim lattices. A replication campaign would be the first place the two systems meet — and fs-evidence/fs-checker can carry GUM-replication claims *natively* | New |

Everything in the Opus capability matrix (which crates exist, what they contain, file-level evidence) **stands** — it was grounded in direct code reading and I re-checked the pivotal claims. What changes is the *strategy built on top of it*.

---

## 1. The reframe: replication, not "simulation of the universe"

"Simulate GUM" is ill-posed at whole-theory scope: no simulator evolves a vacuum with grain scale 10⁻²⁶ m, and the theory's contact with nature runs through S1/S2/bench experiments, not through computation. But the Ω-paper makes a *specific, bounded, checkable* set of numerical claims — the tagged releases ⟨r1⟩–⟨r11⟩ and the exact constants of Sections IV/VII — and the corpus's own governance (Sec. XIII, Branch S) declares external replication of its frozen pipelines *mandatory before any stronger claim*. The archive we hold contains the claims and their derivation specs, **but not the release code itself**. Therefore:

> **The most meaningful GUM simulation program is the first independent replication of the corpus's own numerical claims, executed on an independent codebase, under frankensim's certification discipline.**

This reframing does three things. It converts an unbounded project into a bounded one. It aligns the work with the corpus's own rules (a replication that *fails* is a first-class result — arguably the more interesting one). And it produces artifacts frankensim is uniquely built to certify: every replicated number ships as an `Evidence<T>` with interval enclosures, a model card, and a package `fs-checker` can re-verify without running the solver.

**The decision-relevant ordering.** The corpus's centerpiece is the ħ-closure: `ħ = 𝔠Λ√J` with `𝔠 = 128√42/105π = 2.515`, anchored by the benchmark solve ⟨r1⟩ (`𝔠 = 2.37 ± 0.09` at ε = 0.05) and the ε-scan `𝔠(ε) = 𝔠₀[1 − c_g ε^{2/3}]`. If an independent axisymmetric solve does **not** reproduce that chain, the derivation of Planck's constant — the paper's flagship — collapses within-model. So the replication program has a natural kill-shot tier (Tier 2 below), exactly in the corpus's own falsificationist spirit.

---

## 2. The dimensionality ladder (why this is much cheaper than it looked)

Re-reading Section IV and Appendices B/G/K with implementation in mind, the paper's load-bearing numerics sort by dimension:

| Dim | Claims living there | What the computation actually is | frankensim coverage |
|-----|--------------------|----------------------------------|--------------------|
| **0-D (algebra)** | `ê₀/𝔦₀ = 7/4`, `𝔠₀ = 128√42/105π`, `V = √2`, `E_rot/E = ¼`, `κ(g) = √(7/8g)`, spin-selection `w·j(1−j)=1 ⇒ j=½`, `ω_th = √(7/3)ω₀`, `C₆ = 64Λm̃/15π` (Haar), P-F1′ nulls | Closed-form identities + 1-D quadratures. I re-verified the chain end-to-end this session: `2√(g ê₀𝔦₀)` at `g=3/2` reproduces `𝔠₀` to machine precision; `𝔠₀(1−0.42·0.05^{2/3}) = 2.371` matches ⟨r1⟩ | `fs-cheb::integral`, `fs-ivl` certified enclosures — **complete** |
| **1-D (radial ODE)** | Compacton `f₀(r) = 2arccos(r/R*)`, its direct energy integral (= C₆, "two computations one number"), Derrick virial, hedgehog `∫b = 1`, Yukawa tail `𝔭 = 0.84`, rigid/scaling closure rungs | Chebyshev BVP / constrained 1-D minimization | `fs-cheb` (diff_matrix, roots), `fs-ascent`, `fs-ad` — **complete**, with one hazard (§4) |
| **2-D (axisymmetric)** | **⟨r1⟩ benchmark**: non-rigid isorotating closure — `𝔠 = 2.37±0.09`, `κ = 0.802±0.018`, `g* = 1.31`, `λ* = 0.42`, over-spin onset `κ = 1.000±0.004`; ⟨r6⟩ ε-scan `c_g = 0.42±0.04` | 2-D (r,z) constrained saddle over an axisymmetric SU(2) ansatz (two angle fields), with the ¼ and √2 invariants as built-in self-tests | `fs-ascent` + `fs-adjoint` (Sobolev smoothing) + `fs-cheb`; **new code: the SU(2) axisymmetric energy + outgoing-wave handling** |
| **2-D (ensembles)** | ⟨V-SIM-4⟩ Born relaxation: near-exponential H̄ decay for M ≳ 4, power-law tails for small M | Split-step Schrödinger on a 2-D box + de Broglie guidance trajectories + coarse-grained H̄(t). **Deterministic — no SDE** | `fs-fft` (2-D), `fs-time` (RK45), `fs-rand` (initial ensembles), `fs-eproc` (stopping) — **nearly complete** |
| **k-space (spectrum)** | Four branches B1–B4; exact B2 masslessness; Klein–Gordon B3; band-edge mass (7.1) | The linearized action has a closed-form symbol: assemble a small (6×6-class) matrix at each k, sweep, eig. Floquet dressing = 1-D periodic modulation via FFT. **This tier is nearly trivial** — the hard photon claims (masslessness *to all orders*) are symmetry theorems, verified symbolically, not numerically | `fs-la::jacobi_eigh` / `eigen_complex`, `fs-fft` — **complete** |
| **3-D (fields)** | Dynamic stability of a knot with conserved K; ⟨r2⟩/⟨r3⟩ bootstrap servo at field level; blue-fog condensate 𝔪 = 1.9±0.4, sinθ_c, phason f; disclination web; nucleation | Full SU(2)-texture lattice field theory with second-order inertial dynamics | **The genuinely new build** (the Opus "fs-cosserat" plan) — plus GPU realism (§4) |
| **Under-specified** | ⟨r8⟩ family belts (A = 5.6, B = 5.7), ⟨r10⟩ bridge (A_core = 3.05, κ_far = 0.1065), ⟨r11⟩ quark belts | 2-D belt integrals whose integrands must be *reconstructed from prose* (App. K + Course Ch. 18). Spec-underdetermination risk is real — and would itself be a reportable finding ("the frozen pipeline is not reconstructible from the published record") | Cheap to run once pinned; the risk is specification, not computation |

**The headline consequence:** Tiers 0 through 3 — which cover the paper's centerpiece (the ħ-closure chain), its spectrum claims, and its quantum-sector rates — need no 3-D field type, no FEM, no external engine, and no GPU. The Opus plan's Phase 0 ("build the 3-D SU(2) field first") had the pyramid upside down.

---

## 3. Revised roadmap (tiers, deliverables, kill-content)

**Tier 0 — the constants gauntlet** *(days; pure frankensim).*
Independently re-derive every exact constant and identity of Sections IV/VII: the 7/4 ratio, 𝔠₀, the j-family, V = √2, ¼, ω_th, the Haar quadrature for C₆, `∫b = 1` for the hedgehog. This is AUD-15's §V15.1–V15.6 re-executed *outside* the corpus, shipped as an `EvidencePackage` with `fs-ivl` interval enclosures, checkable by `fs-checker`. *Kill-content:* any mismatch is an arithmetic defect in the corpus (none expected — I've spot-verified the chain).

**Tier 1 — spectrum and Floquet bands** *(days–weeks; fs-la + fs-fft).*
Symbol-matrix sweep ω(k) for B1–B4 from the quadratic action (2.2); verify B2 gaplessness to machine precision and the KG form of B3; exhibit the band-edge mass theorem (7.1) on a 1-D periodic modulation. *Kill-content:* a gapped B2 or non-KG B3 from the stated action would falsify the linearization claims — this is where M-1-class objectivity errors would resurface numerically.

**Tier 2 — the closure chain (the kill-shot tier)** *(weeks–months; the flagship).*
(a) 1-D: compacton and hedgehog with the two-way C₆ check and Derrick virial. (b) Algebraic: the rigid rung (κ = √(7/6), to be *rejected* as superradiant), the scaling and shape-exact rungs. (c) 2-D: the axisymmetric non-rigid isorotating solve at ε = 0.05, targeting the full ⟨r1⟩ tuple (𝔠, κ, g*, λ*, onset at κ = 1.000) with ¼ and √2 as internal self-tests, then the ε-scan exponent 2/3. New code: an axisymmetric SU(2) ansatz energy (two angle fields, `Real`-generic per fs-material's pattern), unit-norm handled by parametrization, radiation onset via a resolvent/quasinormal computation with outgoing-wave boundary conditions on the radial spectral grid. *Kill-content:* **failure to reproduce 𝔠 = 2.37 ± 0.09 kills the ħ-derivation within-model** — the single most valuable possible output of the entire program, in either direction.

**Tier 3 — Born-rule relaxation (⟨V-SIM-4⟩)** *(weeks; cheap, high visibility).*
2-D box, M-mode superpositions, split-step Schrödinger via fs-fft, guidance trajectories via fs-time, coarse-grained H̄(t); reproduce near-exponential decay for M ≳ 4 and the small-M pockets. Acceptance test mirrors the corpus's own: reproduce the published pilot-wave relaxation phenomenology (Valentini–Westman / Towler–Russell–Valentini) first. *No SDE integrator needed* (Δ3); the Nelson-diffusion variant (osmotic velocity, Fisher functional) is an optional follow-on and is where `fs-rand`'s Wiener machinery would finally be used.

**Tier 4 — the 3-D texture field** *(months; the genuinely new build).*
The Opus "fs-cosserat" plan survives here, with two amendments: use a **periodic spectral grid** (not FEM — simpler, matches vacuum studies, plays to fs-fft) and **arrested Newton flow** for statics (cribbing cuSkyrmion/BPS-Skyrme practice) before attempting inertial dynamics via a Verlet⊗Lie splitting (fs-time). Diagnostics: `b_P` density and K (adapting the fs-rep-mesh solid-angle kernel), energy partition, charge conservation under evolution. Cross-validate the S²-sector diagnostics (skyrmion/Hopf numbers) against mumax3 on shared configurations — mumax3 as *referee*, not engine (Δ4).

**Tier 5 — frontier sectors** *(open-ended; flagged risks).*
Blue-fog condensate and disclination web (needs large 3-D boxes — the one place GPU is unavoidable; frankensim is deliberately CPU-deterministic, so this either accepts long runtimes under fs-exec tiling or steps outside the determinism envelope); family/bridge/quark belt integrals (spec-reconstruction risk, Δ2 table); nucleation/IBC. These are research projects, not engineering.

---

## 4. Numerical hazards (new in v2)

1. **The compacton edge is log-divergent in the perturbing sectors.** `f₀′ = −2/(R*√(1−(r/R*)²))` blows up at r = R*, and I verified numerically that the W₂ gradient-energy integrand behaves as `2/u` in the edge distance u — so the raw compacton's quadratic-sector energy diverges logarithmically. This is *consistent with the paper* (App. B.5 calls the sharp edge an ε = 0 artifact; App. H prices the boundary layer at O(ε ln 1/ε)) but it means every ε > 0 solve in Tier 2 must resolve a boundary layer: Chebyshev domain-splitting at the compacton radius, or matched Yukawa-tail asymptotics. This is the stiff spot of the whole program.
2. **Sixth-order nonlinearity.** `½Λ²b_P²` puts the *square of a gradient-determinant* in the energy; its EL equations are quasi-linear and stiff. The BPS-Skyrme numerical literature (the paper's refs [8]/[24]) is the established playbook — arrested Newton flow, not naive gradient descent.
3. **Manifold-valued DOFs.** In 1-D/2-D the clean solution is *parametrization* (angle fields), sidestepping frankensim's per-site-manifold gap entirely; only Tier 4 needs quaternion-field machinery with renormalization or Riemannian steps.
4. **Outgoing-wave boundaries.** The superradiance rejection and onset measurement (κ = 1.000 ± 0.004) require radiation boundary conditions (complex scaling or PML on the radial grid) — standard, but not present in any frankensim crate; small new code.
5. **Determinism vs. GPU.** frankensim's whole value proposition is bit-deterministic cross-ISA replay. Tiers 0–3 keep that. Tier 5's blue-fog boxes realistically need GPU, which breaks the golden-hash discipline; the honest options are metric-only goldens with tolerance bands (the fs-flagship-e2e pattern) or accepting CPU cost.

---

## 5. External code, re-weighted

The Opus external assessment's *facts* stand; the *weights* shift under the replication-first framing:

- **frankensim alone suffices for Tiers 0–3.** This is a stronger pro-frankensim conclusion than v1's polyglot-pipeline emphasis. The crates that matter most: `fs-cheb` (the Orr–Sommerfeld pattern is literally the Tier-1/2 template), `fs-la`, `fs-fft`, `fs-ascent`/`fs-adjoint`, `fs-ad`, `fs-ivl`, `fs-time`, `fs-evidence`/`fs-package`/`fs-checker`.
- **BPS-Skyrme codes (cuSkyrmion-class, arrested Newton flow) replace mumax3 as the closest crib for the knot core** — because W₆₊₀ *is* the BPS-Skyrme model. mumax3/Ubermag remain the right tool for the *collective* sector (S² director, DMI chirality, hopfions) and as an independent referee for Tier-4 topological diagnostics.
- **stark_micropolar** remains the best reference for micropolar *elasticity* patterns (SO(3) DOFs + couple-stress in a variational FEM), but note it is quasi-static/implicit with numerical dissipation — a pattern source, not an engine, for a program whose dynamic tiers need conservative integration. It only becomes relevant if the program ever needs body-fitted 3-D micropolar FEM, which Tiers 0–4 do not.
- **hoomd-blue:** unchanged — not relevant.

---

## 6. What simulation can and cannot adjudicate (unchanged, restated precisely)

A complete, successful execution of Tiers 0–4 would establish: *the corpus's numerical claims are independently reproducible from its stated equations* — upgrading the [CAL] rows from "trust the protocol" to "externally replicated," which is precisely the burden the corpus's own AUD-15 verdict says internal convergence cannot discharge. A failure anywhere — especially Tier 2 — would be a within-model falsification of the affected derivation, reportable in the corpus's own obituary format.

What no simulation can do: test GUM against nature. The constants have no dimensional contact with measured reality until Λ, J, m̃ are fixed by experiment; the theory's empirical exposure remains S1 (Σm_ν), S2 (HL-LHC nulls), the bench ratios, and the A3 gate. A replication program should say so on every artifact it ships — which, conveniently, is what `fs-evidence`'s anti-laundering color lattice and GUM's own language locks *both* enforce. The two governance systems are structurally homologous (Δ7): claims carry graded status, composition is conservative, upgrades require external adjudication. An `fs-gum-replication` campaign would be the first artifact certified under both disciplines at once.

---

## 7. Bottom line (v2)

The Opus assessment's crate-level findings stand; its strategy is superseded. **The right GUM simulation program is a tiered, replication-first campaign that starts where the paper's own numbers live — in 0–2 dimensions — and runs almost entirely on existing frankensim crates**, with the 3-D SU(2) texture field deferred to a fourth tier and external engines (BPS-Skyrme flow methods; mumax3 as referee) entering only there. Its centerpiece and kill-shot is the independent re-derivation of the ⟨r1⟩ closure benchmark (`𝔠 = 2.37 ± 0.09`), on which the corpus's derivation of Planck's constant rests. The program is bounded, certifiable end-to-end with fs-evidence/fs-checker, honest about the one thing it cannot do (adjudicate nature), and — fittingly — is exactly the external replication the corpus wrote into its own rules and then could not perform on itself.
