# fs-gum-cosserat — GUM physics core, Phase B3: micropolar/Cosserat continuum

Run date: 2026-07-12. Binary: `gum-cosserat-gates` (`cargo build --release`,
warning-free; `cargo test --release` 5/5). Battery runtime: **2.9 s** (gate
budget 10 min). Gates: **20/20 PASS**, 20/20 claims certified, fs-checker
3/3 modes as expected, two-run bit-identical replay.

- Spec: `analysis/gum-sandbox/gum-core/gap_survey_v3.json`, dimension
  "Micropolar/Cosserat continuum".
- Port source (line-by-line): `analysis/gum-sandbox/tier1-spectrum/spectrum.py`.
- Referees: `tier1-spectrum/REPORT.md` + `tier1-spectrum/theorem_II2_invariance.csv`
  (read at runtime; the battery's k grid IS the CSV's 400-point grid).
- Certified tier: fs-cosserat-pilot pattern, domain `gum-core:cosserat:v1`,
  producer `fs-gum-cosserat`.
- Golden Merkle root:
  `452d780c84039427fa43331f989492038fdf80265d73455124801b643d5c9a5f`
  (identical across the two in-process runs AND across separate invocations).

## Conventions (stated, matching spectrum.py exactly)

- State q̂ = (u_x, u_y, u_z, φ_x, φ_y, φ_z) ∈ C⁶; plane wave e^{i k·x}.
- W₂ = (λ/2)|e_kk|² + μ e_(ij)e_(ij) + μ_c e_[ij]e_[ij] + (α/2)|Γ_kk|²
  + (β/2)Γ_(ij)² + (γ/2)Γ_[ij]² + mass term (Case A: (m_V²/2)|ψ|² with
  ψ = φ − ½curl u; Case B: (m_V²/2)|φ|²) + χ₃ Re(e_[ij]* Γ_[ij]);
  e_ij = i k_i u_j − ε_ijl φ_l, Γ_ij = i k_i φ_j,
  ψ_a = φ_a − (i/2) ε_abc k_b u_c. Each term c·Σ|T|² contributes
  K += 2c R†R (W = ½ q†Kq); K is Hermitized exactly after assembly.
- T = (ρ₀/2)|u̇|² + (J/2)|φ̇|², M = diag(ρ₀×3, J×3);
  ω² M v = K v reduced via M^(−1/2) K M^(−1/2).
- ARBITRARY k ∈ R³ (spectrum.py hardcodes k ∥ z): the z-axis special case
  reproduces spectrum.py's operator matrices entry-for-entry, and the
  spectrum is verified rotation-invariant (G16).
- Eigensolve: 6×6 complex Hermitian H = A + iB via the real-symmetric
  12×12 embedding S = [[A, −B], [B, A]] on the deterministic cyclic-Jacobi
  `jacobi_eigh`; eigenvalues arrive in exact duplicate pairs (worst pair
  residual observed: 0.0); the complex eigenvector is x + iy from any real
  eigenvector (x; y). Values-only oracle: the complex QR `eig` (G17).
- Benchmark moduli: ρ₀ = J = λ = μ = 1, μ_c = 5, α = β = γ = 0.5, m_V = 1;
  fit windows k ≤ 1e-3 (gapless) and k ≤ 0.1 (gapped), as in Tier-1.

## Gate table

| Gate | What it certifies | Measured | Referee/target | Verdict |
|---|---|---|---|---|
| G01 | c_L = √((λ+2μ)/ρ₀) = √3, small-k fit, rel ≤ 1e-12 | 1.732050807569 (rel 0.0) | 1.732050807569 | PASS |
| G02 | raw ω/k at k = 1e-4, rel ≤ 1e-12 | 1.732050807569 (rel 0.0) | √3 | PASS |
| G03 | mode census at k ∈ {1e-3, 0.1, 1, 3}: 1 long-u + 1 long-φ + 4 transverse, weight > 0.99 | all weights 1.000000 | REPORT D1 | PASS |
| G04 | B2 gap ≤ 1e-12, Case A, m_V ∈ {0,1,5,20} | 7.75e-15 (all four) | 0 (Python 7.9e-15) | PASS |
| G05 | B2 speed c² = μ/ρ₀ = 1 within 1e-6 | 0.999999900 | 1 (Python 0.999999900) | PASS |
| G06 | B3 gap ω₀² = (m_V²+4μ_c)/J = 21 | 20.999999822 (rel 8.465e-9) | 21 (Python 20.999999822) | PASS |
| G07 | factor-4 dial: dω₀²/dμ_c, dω₀²/d(m_V²) | 4.000000000 / 1.000000000 | 4 / 1 | PASS |
| G08 | B4 gap 21 + curvature α+β = 1 | 21.000000000 (rel 3.4e-16) / 1.000000000 | 21 / 1 | PASS |
| G09 | Case B c_B2² = μ + μ_c m_V²/(4μ_c+m_V²), no gap opens | 1.238095104 / 3.777777648 / 5.761904760 | 1.238095238 / 3.777777778 / 5.761904762 (≤1e-6) | PASS |
| G10 | r(k→0) = 4μ_c/(4μ_c+m_V²) (Case B), r = 1 (Case A), 9 digits | 0.952380953 / 0.444444445 / 0.047619048; rA = 1.000000000 | 0.952380952 / 0.444444444 / 0.047619048 | PASS |
| G11 | Case A Thm II.2: IR invariance + structural finite-k drift | drift 2.5716e-4 / 3.0396e-3 / 5.2630e-3; dev(k≤0.01) = 1.45e-11; speed dev 2.9e-15 | REPORT 2.572e-4 / 3.040e-3 / 5.263e-3 | PASS |
| G12 | Case B light branch SHIFTS | max rel dev 1.9772 | 1.977 (±1e-3) | PASS |
| G13 | χ₃ = 0.3: light split ∝ k⁵; split maxima | exponent 4.994; light 3.415785e-2; heavy 1.971132e-1 | 5±0.1; 3.415785e-2; 1.971132e-1 (Python 4.994) | PASS |
| G14 | decoupling cross-check (u_z ∪ φ_z ∪ transverse vs full 6×6) | 0.0 | ≤ 1e-11 (Python 3.4e-13) | PASS |
| G15 | global min ω² ≥ −1e-12 | +1.000e-8 | ≥ −1e-12 (Python min 1.0e-8) | PASS |
| G16 | rotation invariance, 8 random k̂ × \|k\| ∈ {1e-3,0.1,1,3} × χ₃ ∈ {0,0.3} | 2.744e-15 (scaled) | ≤ 1e-12 | PASS |
| G17 | complex-QR oracle vs embedding; pair residual | 2.01e-15 re / 2.3e-17 im / pair 0.0 | ≤ 1e-10 | PASS |
| G18 | CSV curves ω_B2(k), 400 k × 2 cases × 4 m_V | worst dev/tol 0.009; rel(ω) k≥0.1: 1.54e-12 | conditioning-aware tol; literal 1e-9 for k ≥ 0.1 | PASS |
| G19 | TIME DOMAIN: branch frequencies | max vs ω_h 1.5e-14; vs ω 4.2e-8; B3(k→0) vs √21 1.79e-7 | 1e-9 / 1e-6 / 1e-6 | PASS |
| G20 | TIME DOMAIN: energy over 2e5 Verlet steps/mode | max dev 2.500e-7; envelope growth 2.7e-14 | ≤ 1e-5 / ≤ 1e-9 | PASS |

Checker modes: deny-all = false (anti-laundering works), cert-capability =
true, tampered root = false — all as expected. Replay: second full battery
run (eigensolves, sweeps, Verlet trajectories, certificates, Merkle
assembly) reproduces the root bit-for-bit.

## Time-domain validation — the repo's first

This is the first TIME-DOMAIN validation of the coupled second-order
micropolar system anywhere in the repository (Tier-1 was frequency-domain
only; the survey lists it as a missing piece). The mass-matrix-aware
symplectic Verlet (kick–drift–kick with q̇ = M⁻¹p) evolves the 12-real
embedding of each branch's complex plane-wave amplitude for 2×10⁵ steps
(h = 1e-3/ω per mode, so ω_h/ω − 1 = (hω)²/24 ≈ 4.2e-8):

| Mode | ω (eigensolve) | ω_h (discrete) | ω measured | vs ω_h | vs ω | max energy dev | envelope growth |
|---|---|---|---|---|---|---|---|
| B2 light transverse, k=1 | 0.947870595456 | 0.947870634951 | 0.947870634951 | 1.5e-14 | 4.2e-8 | 2.500e-7 | 2.0e-14 |
| B1 long. acoustic, k=1 | 1.732050807569 | 1.732050879738 | 1.732050879738 | 1.2e-14 | 4.2e-8 | 2.500e-7 | 0.0 |
| B4 long. twist, k=1 | 4.690415759823 | 4.690415955257 | 4.690415955257 | 4.2e-15 | 4.2e-8 | 2.500e-7 | 0.0 |
| B3 heavy transverse, k=1 | 5.181847289748 | 5.181847505658 | 5.181847505658 | 1.1e-14 | 4.2e-8 | 2.500e-7 | 2.7e-14 |
| B3 at k=1e-3 (KG ring) | — | — | 4.582576513273 | — | vs √21 = 4.582575694956: 1.786e-7 | 2.500e-7 | 0.0 |

The measured frequency equals the DISCRETE Verlet frequency
ω_h = (2/h)·asin(hω/2) to ~1e-14 (the estimator uses the exact three-term
identity s_{n+1} − 2s_n + s_{n−1} = −4 sin²(θ/2) s_n averaged over the
run), and the continuum ω to the expected O(h²ω²/24) ≈ 4.2e-8 — 1e-6-class
as gated. Energy stays inside the bounded Verlet oscillation
(amplitude 2.5e-7·E, consistent with (hω)²/8) with first-to-second-half
envelope growth ≤ 2.7e-14 over 2×10⁵ steps: symplectic, no secular drift.

## Deviations from the letter of the spec (with reasons)

1. **fs-la is vendored, not path-depended.** The survey's build spec calls
   for `fs-la` (jacobi_eigh + eigen_complex oracle). In this checkout fs-la
   is unbuildable from any out-of-workspace crate: its dependency closure
   pulls `fs-exec`, whose Cargo.toml requires
   `asupersync = { path = "../../../asupersync" }` — a sibling repository
   not present next to fork_frankensim. The two kernels this crate needs
   are pure std + fs-math code, so `src/fsla_vendored.rs` carries them
   VERBATIM (function bodies byte-identical, provenance header with exact
   source lines); Cargo.toml documents the one-line swap back once fs-la
   builds from here.
2. **CSV match uses a conditioning-aware tolerance below k = 0.1.** The
   survey's blanket "rel ≤ 1e-9 at all 400 points" is not achievable by any
   reimplementation: the CSV's own ω values carry numpy/LAPACK absolute
   eigenvalue noise ~eps·‖H‖, which at k ~ 1e-4 (ω² ~ 1e-8 against
   ‖H‖ up to ~840 at m_V = 20 Case B) is ~1e-7 RELATIVE — visibly present
   in the CSV itself (Case A columns differ from each other by ~6e-8 at
   row 1 where they are analytically near-identical). Gate G18 therefore
   enforces |ω²_ours − ω²_csv| ≤ 1e-9·ω²_csv + 64·eps·tr(H_t) pointwise
   (worst observed ratio 0.009 — margin ×100) AND the literal 1e-9
   relative on ω for k ≥ 0.1 (worst observed 1.54e-12, margin ×650).
3. **Gate thresholds where the Python pilot's own gates were looser than
   the survey's prose:** G06 gates the B3 gap at rel 1e-6 (spectrum.py's
   own d4 gate) because the measured value 20.999999822 (rel 8.5e-9) is a
   deterministic fit-window bias shared bit-for-bit with the pilot, not
   noise; the note records the 9-digit value. Everything the task listed
   numerically (√3 to 1e-12, B2 gap ≤ 1e-12, factor-4 = 4.000000000,
   Case B c², r to 9 digits, χ₃ exponent ≈ 5, CSV, replay) is gated at the
   stated strength.
4. Not deviations, but worth noting: spectrum.py's Case A pointwise
   invariance (D6 literal < 1e-12) and χ₃ all-k degeneracy (D7 literal)
   FAILED in the Python pilot for structural reasons its REPORT documents;
   this port reproduces the same structural finite-k numbers (G11, G13)
   rather than the impossible literal claims.

## Files

- `Cargo.toml` — standalone crate (empty `[workspace]`, path deps into
  `crates/`, fs-cosserat-pilot pattern).
- `src/lib.rs` — crate docs + exports.
- `src/moduli.rs` — `Moduli`/`MassCase`, Tier-1 benchmark.
- `src/symbol.rs` — arbitrary-k operator matrices, tensor splits, K(k)
  assembly, mass matrix (+ unit tests incl. the e_[ij]e_[ij] = 2|ψ|²
  fold-in identity).
- `src/eigh.rs` — Hermitian eigensolve via the 12×12 real-symmetric
  embedding, pair dedup, eigenvector recovery, generalized M^(−1/2)
  reduction, QR oracle (+ residual tests).
- `src/branches.rs` — sub-block branch extraction, doublet pairing,
  co-motion fidelity, small-k and log-log fits, classification.
- `src/verlet.rs` — mass-matrix-aware symplectic Verlet, plane-wave mode
  evolution, frequency estimator, energy-envelope diagnostics.
- `src/fsla_vendored.rs` — verbatim fs-la `jacobi_eigh` + complex `eig`
  (see Deviations, item 1).
- `src/bin/gates.rs` — the 20-gate referee battery + certified tier +
  replay (`gum-cosserat-gates`).

## Epistemic notice (binding)

Every PASS above certifies a WITHIN-MODEL property of a speculative
theory's linearized dispersion relation. It validates the port of the
quadratic action, the complex-Hermitian eigensolve path, the arbitrary-k
generalization, and the symplectic integrator — never the physics.

---

# Phase E2 addendum — W_χ chiral-coupling module (dimension "Skyrme/chiral")

Run date: 2026-07-12. Same binary (`gum-cosserat-gates`), extended battery:
**26/26 PASS** (the 20 Phase B3 gates above, unchanged and re-passing, plus
G21–G26), 26/26 claims certified, fs-checker 3/3 modes as expected, two-run
in-process bit-identical replay AND identical root across separate
invocations. `cargo build --release` warning-free; `cargo test --release`
13/13. Battery runtime: **3.0 s** (added-gate budget 2 min). The Phase B3
golden root above corresponded to the 20-claim package; the extended
battery's golden Merkle root is
`9ddc0fe40548bf9c97cd505ef30339bbe8b4b8ff537ef6f3ef64377cb5f57b15`.

- Spec: `gap_survey_v3.json`, dimension "Skyrme/chiral", the W_chi
  missing-pieces item (three-coupling density + gradient, free χ's,
  trace-vs-deviatoric switch, Dzyaloshinskii soft-sector functional with
  𝔪 = χ²/(γΔ²) as the only text-pinned referee, hard χ = 0 default).
- Corpus source: `substrate-suite/01-GUM-Omega-Paper-v2.0.1.md` Sec. II.C
  eq (2.2) (the W_χ term), Sec. II.H (the chiral vacuum; the E-H1
  Dzyaloshinskii audit), App. K.5; Substrate Course Ch. 6 (eq 6.1–6.3).
- Referee cross-check: Tier-5a `family_audit.py` C7a
  (sinθ_c = √(1−1/𝔪) = 0.688 → 0.69 at 𝔪 = 1.9).

## New modules

- `src/wchi.rs` — real-space W_χ energy for STORED micropolar fields
  (u, φ) on a periodic cell-centred N³ grid: e_ij = ∂_i u_j − ε_ijk φ_k,
  Γ_ij = ∂_i φ_j, W_χ = χ₁e_kkΓ_ll + χ₂e_(ij)Γ_(ij) + χ₃e_[ij]Γ_[ij];
  density, total energy (h³-weighted, fixed cell-major order), and the
  analytic gradient w.r.t. BOTH fields in two-pass flux form
  (P^e = ∂w/∂e, P^Γ = ∂w/∂Γ, central-difference adjoint = −D, pointwise
  −ε term into ∂/∂φ). Second-order central differences, periodic wrap.
- `src/softsector.rs` — the Dzyaloshinskii soft-sector condensate
  functional of corpus II.H(i): f(θ,q) = sin²θ(−χ_s q + ½γ_s q²) + gap;
  margin 𝔪, pitch q\*, criterion, analytic + numeric tilt (deterministic
  120-iter golden section, Newton/parabolic polish).
- `src/symbol.rs` extension — K(k) now assembles ALL THREE χ couplings
  (χ₃ path refactored into the shared Hermitized cross-term helper with
  the identical loop order — bit-identical results; χ₁/χ₂ added the same
  way), plus the public `deviatoric` projector.
- `src/moduli.rs` — `chi1`, `chi2`, `chi_deviatoric` fields (defaults 0 /
  trace-included; `bench()` unchanged in every existing value).

## Gate table (added rows)

| Gate | What it certifies | Measured | Referee/target | Verdict |
|---|---|---|---|---|
| G21 | W_χ real-space analytic gradient vs central FD, random smooth periodic fields, 8³ grid, χ's singly + mixed, both conventions, rel ≤ 1e-6 | worst rel 4.35e-13 | ≤ 1e-6 (W_χ quadratic ⇒ FD exact to round-off) | PASS |
| G22 | dispersion-path consistency: stored-plane-wave W_χ energy = (V/4)v†K_χ(k_eff)v, same couplings, 24 random (k, polarization) × {χ₁, χ₂, χ₃, mixed} × both conventions | worst scaled dev 3.05e-15 | ≤ 1e-12 (exact discrete identity, k_eff = sin(kh)/h) | PASS |
| G23 | convention ledger: deviatoric(χ₁,χ₂) ≡ trace-included(χ₁−χ₂/3,χ₂) (symbol 4.28e-16, real-space density 2.99e-16); χ = 0 default exactly zero density/energy/gradient; spectrum rotation-invariant with χ₁/χ₂/χ₃ on | 3.27e-15 (rotation) | ≤ 1e-13 / exact zeros / ≤ 1e-12 | PASS |
| G24 | Dzyaloshinskii tilt at the r9 margin 𝔪 = 1.9: minimize f(θ) → sinθ_c = √(1−1/𝔪) = 0.688247201612 (Tier-5a C7a 0.688 → 0.69), stationarity ≤ 1e-12, pitch q\* = χ_s/γ_s ≤ 1e-9, inside 0.69 ± 0.11 | sinθ_c 0.688247201611685 (dev 0.0), f′/Δ² 1.3e-16, pitch rel 4.5e-15 | √(1−1/1.9), Thm H-3 | PASS |
| G25 | 𝔪-scan [1.1, 3], 39 points, three (γ_s, Δ_s) parametrizations (ε-blindness: only 𝔪 matters): numeric tilt vs √(1−1/𝔪) ≤ 1e-10 | worst 3.33e-16 | ≤ 1e-10 | PASS |
| G26 | Dzyaloshinskii threshold at 𝔪 = 1: no tilt below (θ_c ≤ 1e-6 at 𝔪 ∈ {0.3, 0.6, 0.9, 0.999}, f″(0) > 0), tilt = analytic above (𝔪 ∈ {1.001, 1.01, 1.1}, ≤ 1e-8), FD-bisected 𝔪\* = 1.000000000 (±1e-6); literal quadratic-gap reading shares the threshold | 𝔪\* 0.9999999998 | 1 | PASS |

## Convention ledger — what the corpus text PINS vs leaves FREE

PINNED (with location):

1. The three-coupling FORM W_χ = χ₁e_kkΓ_ll + χ₂e_(ij)Γ_(ij) +
   χ₃e_[ij]Γ_[ij], coefficient 1 on each contraction, on
   e_ij = ∂_iu_j − ε_ijkφ_k and Γ_ij = ∂_iφ_j — II.C eq (2.2) ("chiral
   transduction, F2"). Parity-odd (tensor × pseudotensor); lives on the
   MICROPOLAR (u, φ) fields, not the quaternion texture.
2. Units class [χ_i] = E·L⁻¹ (App. A, per the survey entry).
3. χ₃ quenched on B2 at tree level, (ka)²-suppressed regeneration —
   Thm II.3 / App. C.6 (the Tier-1 k⁵ splitting, gate G13, is its
   validated phenomenology).
4. Soft sector (II.H(i) / Course 6.1): f(θ,q) = sin²θ(−χ_sq + ½γ_sq²)
   + ½Δ_s²θ² TO QUADRATIC ORDER in θ; pitch q\* = χ_s/γ_s; criterion
   χ_s² > γ_sΔ_s² (threshold 𝔪 = 1).
5. Thm H2-1 + ⟨r9⟩: 𝔪 = χ̄²/(γ̄Δ̄²) = 1.9 ± 0.4 — the ONLY text-pinned
   NUMBER in the chiral sector.
6. Thm H-3 / K.5: condensate tilt sinθ_c = √(1 − 1/𝔪) = 0.69 ± 0.11.

FREE (module behavior):

1. Values of χ₁, χ₂, χ₃ — runtime parameters, DEFAULT 0. At the default
   the real-space density/energy/gradient are hard-guarded exact zeros
   (G23b) and the symbol path skips every χ accumulation, so the Phase B3
   spectrum battery re-passes bit-for-bit and the campaign's χ-free
   knot-statics functional E_static = t(E2+E4)+E6+E0 is untouched.
2. Trace-included vs deviatoric χ₂ contraction — the text never says
   (survey missing-pieces item). Both implemented behind
   `Moduli::chi_deviatoric` (default trace-included, matching this
   crate's W₂ convention); they are the exact reparametrization
   deviatoric(χ₁, χ₂) ≡ trace-included(χ₁ − χ₂/3, χ₂), gated (G23a) —
   which is WHY the text can afford to leave it unpinned.
3. χ_s, γ_s, Δ_s individually — only 𝔪 and q\* are pinned combinations;
   G25 runs three different (γ_s, Δ_s) parametrizations per margin.
4. Any W_χ ↔ quaternion-knot coupling — none appears in any solved
   statics functional anywhere in the corpus or campaign.
5. **The finite-θ gap completion (documented finding).** The printed
   ½Δ_s²θ² gap taken literally to ALL orders is inconsistent with the
   pinned Thm H-3 tilt: its stationarity is 𝔪 sinθcosθ = θ, giving
   sinθ_c ≈ 0.7924 at 𝔪 = 1.9, not 0.688. Requiring (a) the printed
   chiral term −(χ_s²/2γ_s)sin²θ and (b) sin²θ_c = 1 − 1/𝔪 for EVERY
   𝔪 > 1 forces G′(θ) = Δ_s²tanθ, i.e. determines the gap potential
   UNIQUELY: G(θ) = −Δ_s² ln cosθ = ½Δ_s²θ² + O(θ⁴) — the printed term
   is exactly its quadratic truncation. `softsector` implements BOTH
   readings (`f_tilt` = completed, `f_tilt_quadgap` = literal); the tilt
   gates run on the completed form, the literal form's identical
   threshold is co-gated (G26), and the corpus criterion/pitch hold in
   both.

## Discretization convention (this module's, stated)

Second-order central differences on the periodic grid; the discrete plane
wave is an exact eigenfunction with symbol i·k_eff, k_eff = sin(kh)/h, so
G22's real-space-vs-symbol comparison is an EXACT identity (measured
3e-15), cross-validating the grid energy, the gradient's flux fields, and
the χ₁/χ₂-extended symbol assembly in one shot. Real fields
u = Re(v̂e^{ik·x}) average quadratics to half the complex form and
W = ½q̂†Kq̂ supplies another half: E_grid = (V/4)v̂†K_χv̂.

## Deviations from the letter of the spec (with reasons)

1. **The tilt gate does not minimize the LITERAL printed f(θ).** The spec
   says "implement the functional the paper states" and "verify
   θ_c matches sinθ_c = √(1−1/𝔪)"; as documented above these two demands
   are mutually exclusive under the all-orders literal reading (0.7924 vs
   0.688 at 𝔪 = 1.9), and the corpus itself only ever uses the printed
   density in small-θ expansion while separately printing the tilt as
   Thm H-3. Resolution: implement both readings, gate the tilt on the
   uniquely-determined −Δ²ln cosθ completion (whose quadratic truncation
   IS the printed gap), co-gate the threshold on the literal reading.
2. G21's 1e-6 FD gate passes at 4e-13 because W_χ is exactly quadratic in
   the fields — central FD is exact to round-off; the gate threshold is
   kept at the spec'd 1e-6 (the margin is structural, not tuned).
3. G22's "1e-12-class" is enforced at ≤ 1e-12 scaled by the
   absolute-value quadratic form (V/4)Σ|v||K||v| — the natural scale for
   an indefinite cross term whose value can pass through zero.

## Epistemic notice (binding, unchanged)

Every new PASS certifies a WITHIN-MODEL property: the internal
consistency of the W_χ real-space and symbol implementations, and the
within-model threshold structure of a speculative theory's soft-sector
functional — never the physics.
