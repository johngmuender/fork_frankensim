# Can frankensim actually simulate GUM? — a capability audit and extension plan

**Question.** Given the frankensim workspace (126 `fs-*` crates, a "certified simulation" stack), what — if anything — can be used to *truly* simulate the GUM chiral-micropolar physics as described in the GUM-Ω v2.0.1 paper, and what would have to be extended or built to model the speculative program?

**Method.** I read the workspace manifest and grounded on the ~11 most GUM-critical crates directly (`fs-feec`, `fs-ga`, `fs-solid`, `fs-material`, `fs-time`, `fs-cheb`, `fs-lattice`, `fs-symmetry`, `fs-couple`, `fs-flux`, `fs-topols`), then ran six parallel deep-reads over ~55 crates clustered by capability: (A) continuum/FEEC/elasticity, (B) topology/soliton/geometry, (C) solvers/spectral/numerics, (D) AD/optimization/variational, (E) stochastic/quantum/UQ, (F) certification + e2e campaign templates. Findings below are grounded in specific types/functions with file references.

**One-sentence answer.** frankensim is an unusually strong *numerical, variational, and certification foundation* — its symplectic + Lie-group integrators, discrete exterior calculus, eigen/spectral solvers, adjoint-based minimizers, and evidence/ledger discipline are directly reusable for the **deterministic classical-field sectors** of GUM (soliton profiles, the linear spectrum, the ħ-closure numerics) — but it contains **none of the GUM-specific physics objects** (an SU(2)-valued *field*, micropolar/couple-stress continuum mechanics, the topological-charge density, the Skyrme/chiral energy terms), and it is essentially **irrelevant to GUM's stochastic-QM and gravitation sectors**; a real GUM simulator is therefore a *new physics core built on top of frankensim's numerics*, not a reuse of any existing crate.

---

## 1. What GUM actually requires to be simulated (from the Ω paper)

The paper's physics is a classical field theory of a chiral micropolar (Cosserat) solid, quantized stochastically. A faithful simulator needs, in dependency order:

| # | Requirement | Paper anchor |
|---|-------------|--------------|
| R1 | **Two coupled fields:** displacement `u∈ℝ³` and an independent **micro-rotation texture `Q̃(x)∈SU(2)`** (orientation at every point) | Sec. II.A |
| R2 | **Relative/objective kinematics:** polar decomposition `F=R[u]U`, relative texture `P̃=R̃[u]†Q̃`, relative strain `e_ij`, wryness `Γ_ij=∂_iφ_j`, Maurer–Cartan currents `L_i=Q̃⁻¹∂_iQ̃` | Sec. II.A–B |
| R3 | **The stored energy** `W = W₂ + W_χ + W₄ + W₆₊₀`: micropolar elasticity + locking + wryness (asymmetric stress, **couple-stress**); chiral parity-odd transduction `χ e·Γ`; Skyrme quartic `Tr([L_i,L_j]²)`; sextic topological `½Λ²b²` + potential `𝒱(σ_P)` | Sec. II.C |
| R4 | **Topological charge:** degree `K=∫b d³x`, `b=−(1/24π²)ε_ijk Tr(L_iL_jL_k)` (`π₃(S³)`); Hopf charge (`π₃(S²)`) for the neutrino; disclination winding (`π₁`) for electric charge | Sec. II.B, VII.I, V.D |
| R5 | **Soliton solver:** minimize `W` at fixed `K` to get the Bogomolny compacton `f₀(r)=2arccos(r/R*)`, hedgehog profile, Derrick virial check | Sec. VII.B–C |
| R6 | **Linear spectrum:** dispersion branches `ω(k)` — acoustic B1, gapless photon doublet B2, Klein–Gordon matter B3, pump B4 — plus Floquet–Bloch band structure of the periodic substrate | Sec. II.E, VII.A |
| R7 | **The ħ-closure:** the isorotating knot — inertia integral `𝕀=(16πJ/3)∫r²sin²f dr`, spin-clock constraints `L=ħ/2`, `ħω=E`, invariants `E_rot/E=¼`, `V=√2`, radiative-consistency (reject superradiant) | Sec. IV |
| R8 | **Dynamics:** second-order-in-time symplectic evolution of the coupled `u`/`Q̃` fields; SO(3)/SU(2)-preserving update of the orientation field | Sec. II.D |
| R9 | **Stochastic-QM sector:** Nelson diffusion `dX=b dt+dW`, osmotic velocity `u=ν∇ln ρ`, quantum potential = δ(Fisher info)/δρ, Madelung closure, Born-rule H-theorem relaxation | Sec. III |
| R10 | **Gravitation:** Einstein–Cartan defect geometry, cone-slaving, Gibbs–Duhem vacuum energy | Sec. VI |
| R11 | **Certification:** dimensional units, deterministic replay, evidence colors, golden hashes, stability certificates | Sec. IX–X |

R1–R8 are the *classical-field core* (tractable). R9 is a *different numerical paradigm* (stochastic). R10 is *research-grade differential geometry*. R11 is exactly what frankensim is built for.

---

## 2. Capability matrix — requirement × repository

Legend: **✓ reusable** (exists, use largely as-is) · **◑ extend** (a close relative exists; adapt/generalize) · **✗ build new** (nothing usable).

| Req | Status | What exists / what's missing | Key crates (file refs) |
|-----|--------|------------------------------|------------------------|
| **R1** field pair `u`, `Q̃(x)` | ✗ build new | A *rigid-body* rotor exists; a **continuum orientation field** does not. `fs-solid::rod` is a **1-D Cosserat rod** with per-node unit-quaternion directors updated multiplicatively — the exact pattern, but 1-D. No 3-D `Q:mesh→SU(2)` container anywhere. | `fs-solid/src/rod.rs:37` (quats), `:147` (exp update); `fs-feec` cochains are scalar-per-cell only |
| **R2** objective kinematics, `L=Q⁻¹dQ` | ◑/✗ | `fs-ga` has the *pointwise* algebra (unit quaternions=SU(2), `gp`, `reverse`=inverse, `motor_log`/`exp_bivector` Lie exp/log, half-angle) but **no field-derivative operator** — `dQ` and `L_i` on a texture are absent. Polar decomposition `R[u]` not exposed. | `fs-ga/src/facade.rs:82,107`; `pga.rs:434,457`; **no** Maurer–Cartan of a field |
| **R3** stored energy (micropolar+chiral+Skyrme+sextic) | ◑/✗ | The **AD-energy pattern is ideal**: `fs-material` writes `energy<T:Real>` and gets exact stress/tangent by AD. But the public surface is a *closed enum* of two symmetric isotropic hyperelastic laws in Voigt-6 space — **no asymmetric strain, no wryness argument, no couple-stress conjugate, no chiral/Skyrme/sextic terms, no user closure**. `fs-solid` continuum is symmetric-stress Cauchy. | `fs-material/src/hyper.rs:107` (`energy<T:Real>`), `:138` (AD Piola); Cosserat/couple-stress = explicit no-claim in `fs-solid/CONTRACT.md:255` |
| **R4** topological charge `K`, Hopf, disclination | ✗ build new | **The single sharpest gap.** Every "topology" invariant in the repo is of a *region/shape*, never of a *map*: `fs-topo` exact Betti/Euler of solids, `fs-feec::betti` of complexes, `fs-rep-mesh::winding_exact` = solid-angle **degree of a surface's Gauss map** (van Oosterom–Strackee) — a real degree kernel, but for inside/outside, not for an SU(2) texture. No `b_P`, no Hopf linking, no `π₁` disclination winding. | `fs-rep-mesh/src/winding.rs:36,55` (adaptable kernel); `fs-topo/src/cubical.rs:106`; **no map-degree** |
| **R5** soliton minimizer | ◑ | Strong: `fs-ascent` L-BFGS + trust-region-Newton-Krylov + augmented-Lagrangian (KKT certs); `fs-adjoint` matrix-free IFT adjoint + **Sobolev/H¹ gradient smoothing** (needed to kill grid-noise in a field gradient); `fs-ad` `Real`-generic energy + FD gradient gate; **`fs-topopt` is a worked "minimize a PDE-energy-of-a-field at fixed integral with exact sensitivities."** Missing: **per-site manifold (S²/SU(2)) DOFs** (only whole-variable Riemannian opt), and a **fixed-`K` (Bogomolny) constraint**. | `fs-ascent/src/{lbfgs,trust,auglag,riemann}.rs`; `fs-adjoint/src/{ift,sobolev}.rs`; `fs-topopt/src/lib.rs:37` |
| **R6** linear spectrum / dispersion | ✓ | **Directly covered.** `fs-la` symmetric (`jacobi_eigh`), matrix-free sparse (`lanczos`, `lobpcg`), and complex nonsymmetric (`eigen_complex::eig`). `fs-cheb::orr_sommerfeld` is a **ready-made template**: assemble a collocation operator pair `(A,B)`, reduce the generalized eigenproblem `Aφ=λBφ`, extract sorted growth rates — exactly the B1–B4 branch computation. `fs-fft` 1–3D for Bloch band structure. | `fs-la/src/eigen*.rs`; `fs-cheb/src/orr_sommerfeld.rs:44–114`; `fs-fft` |
| **R7** ħ-closure numerics | ◑ | All *primitives* exist: `fs-cheb` for the radial profile `f(r)` BVP (`diff_matrix`, `dirichlet_laplace_eigs`, build/differentiate/roots), `fs-cheb::integral` for `𝕀` and `C₆` Haar quadratures, `fs-ascent` AL for the constrained isorotation solve, `fs-ivl` for certified roots. The *coupled closure system* (spin+clock constraints selecting `𝔠`, `κ`) is bespoke assembly on these. | `fs-cheb/src/lib.rs:330,341,378`; `fs-ivl/src/newton.rs` |
| **R8** symplectic + Lie-group dynamics | ✓ (primitives) | **A real asset.** `fs-time` has Störmer–Verlet **with a discrete adjoint** (`verlet_step`/`verlet_adjoint`), SO(3) Lie-group exp-map integration (`lie::quat_exp_step`, `rigid_body_step`, norm-drift ~1e-12 over 1e5 steps), generalized-α, IMEX/exponential. These are the *correct* integrator classes for the coupled `u` (symplectic) + `Q̃` (Lie-group) evolution — but they step **externally-assembled** ODE-scale systems; wiring them to an assembled field PDE is new code. | `fs-time/src/{symplectic,lie,galpha,stiff}.rs` |
| **R9** Nelson stochastic QM | ✗ build new | Substrate only: `fs-rand` counter-based Philox + Gaussian draws + Sobol/QMC (Wiener increments assemblable), `fs-eproc` e-processes (H-theorem *stopping* signal), `fs-uq` MC/QMC, `fs-assimilate` Kalman. **No SDE integrator, no Fisher-information functional, no H-theorem diagnostic, no Madelung closure exist anywhere** (verified: no `sde|langevin|leapfrog` in these crates; `fs-uq/kl.rs` is Karhunen–Loève, not KL divergence). | `fs-rand/src/{philox,lib}.rs`; `fs-eproc`; **no dynamics** |
| **R10** gravitation | ✗ out of scope | Nothing: no Einstein–Cartan, no torsion/curvature from defect density, no cosmology. (`fs-ga` PGA is projective rigid-motion geometry, not spacetime geometry.) | — |
| **R11** certification wrapper | ✓ | **Exactly frankensim's purpose.** `fs-evidence` `Evidence<T>`/`Certified<T>` with anti-laundering color lattice; `fs-qty` compile-time + runtime SI dimensions; `fs-package`+`fs-checker` build a claim bundle a *solver-free* checker re-verifies; `fs-ledger` write-time color gate + deterministic replay; `fs-math` bit-identical cross-ISA determinism; `fs-sos` Lyapunov/PSD certificate (limited). e2e templates: `fs-metamat-e2e`/`fs-flutter-e2e` are the canonical "Cert" skeleton. | `fs-evidence/src/lib.rs:705,742,960`; `fs-checker`; `fs-metamat-e2e` |

**Score:** of the eleven requirements — **2 fully reusable** (R6 spectrum, R8 integrators; plus R11 certification), **3 extend-a-relative** (R3, R5, R7), **4 build-new** (R1, R2, R4, R9), **1 out-of-scope** (R10). The reusable set is precisely the *generic numerics and certification*; the build-new set is precisely the *GUM-specific physics*.

---

## 3. What is directly reusable (the foundation is real)

These are usable essentially as-is and are genuinely strong:

- **Discrete exterior calculus** — `fs-feec`: cochains on tet complexes, Whitney P₁Λᵏ forms, discrete Hodge stars, **exact `dd=0`** integer incidence, first-kind Nédélec H(curl) / Raviart–Thomas H(div) vector families (r=1–4), tensor de Rham complex, and a Hodge decomposition into exact⊕coexact⊕harmonic. This is the right structure-preserving discretization backbone — it "kills spurious pressure/EM modes structurally instead of by stabilization folklore," which matters for a gauge-field-like sector.
- **Structure-preserving time integration** — `fs-time`: symplectic Verlet (+ discrete adjoint via revolve checkpointing) and SO(3) Lie-group exp-map integrators. The two ingredients a coupled `u`/`Q̃` wave equation needs, both tested and drift-free.
- **Eigen/spectral/Krylov** — `fs-la` (Jacobi/Lanczos/LOBPCG/complex-QR), `fs-cheb` (collocation + generalized-eigen BVP, the Orr–Sommerfeld pattern), `fs-solver` (CG/MINRES/GMRES/p-MG, transposed solves), `fs-sparse` (SaAMG/ILU0/Chebyshev preconditioners, SpGEMM Galerkin), `fs-fft` (1–3D). Everything a dispersion-branch or band-structure computation needs.
- **Variational minimization + adjoints** — `fs-ascent` (L-BFGS / TR-Newton-Krylov / augmented-Lagrangian with KKT certs / Riemannian L-BFGS on S^{n−1}, SO(3)), `fs-adjoint` (matrix-free IFT adjoint, Sobolev smoothing, checkpointed transient adjoint), `fs-ad` (`Real`-generic AD + FD gradient gate). `fs-topopt` is a directly analogous worked example.
- **Pointwise Lie/Clifford algebra** — `fs-ga`: unit quaternions (=SU(2)), geometric product, grade projection, exact Lie exp/log, `slerp` — deterministic, bit-identical.
- **Certified arithmetic + determinism** — `fs-ivl` (intervals, Taylor models, Krawczyk-certified roots), `fs-math` (cross-ISA bit-identical transcendentals; the whole stack records golden hashes).
- **The certification wrapper** — `fs-evidence`/`fs-qty`/`fs-package`/`fs-checker`/`fs-ledger`/`fs-sos`, with `fs-metamat-e2e`/`fs-flutter-e2e` as the ~150-line campaign template.

**Reusable numerical building blocks the paper itself would want:** the Haar/quadrature machinery for `C₆=64Λm̃/15π` (via `fs-cheb::integral`), the Chebyshev radial profile for `f₀(r)`, the generalized-eigen pattern for the branch spectrum, and the certified-root machinery to pin invariants like `E_rot/E=¼` to an interval.

---

## 4. What must be extended or built (the physics is absent)

### 4.1 The linchpin gap: a manifold-valued *field*
Nothing in the repo represents a field of orientations `Q:mesh→SU(2)` (or a director `n:mesh→S²`). `fs-ga` gives the algebra at a point; `fs-solid::rod` gives it along a 1-D curve; `fs-feec` gives scalar-valued forms. A GUM core must introduce:
- a **manifold-valued nodal field** type (unit quaternion per vertex, or `n∈S²` per vertex), reusing `fs-ga` for the pointwise algebra and `fs-time::lie` for updates;
- the **discrete Maurer–Cartan current** `L_i=Q̃⁻¹∂_iQ̃` on mesh edges (a new operator, but expressible with `fs-feec` incidence + `fs-ga` products);
- the **relative texture** `P̃=R̃[u]†Q̃`, which needs a polar-decomposition `R[u]` of the deformation gradient (small addition; `fs-ga`/`fs-la` can supply it).

Everything else in the physics core hangs off this type.

### 4.2 Micropolar continuum mechanics (extend `fs-solid`/`fs-material`)
- Generalize `fs-solid`'s rod kinematics (SO(3) directors, multiplicative updates, Cosserat Γ/κ strains — `rod.rs`) to a **3-D continuum**: asymmetric strain `e_ij`, wryness `Γ_ij`, **couple-stress** conjugate, and the asymmetric-stress torque balance `Jφ̈=∂m+ε:T`.
- Open `fs-material` to a **user-supplied stored energy** `W(e,Γ,P̃)` generic over `Real` (the AD-energy mechanism already yields exact stress/tangent), including the parity-odd chiral cross terms `χ e·Γ`, the Skyrme quartic `Tr([L_i,L_j]²)`, and the sextic `½Λ²b² + 𝒱(σ_P)`. This is the four-installment potential `W₂+W_χ+W₄+W₆₊₀`.

### 4.3 Topological-charge machinery (build new, adapt `fs-rep-mesh`)
- The **baryon density** `b_P=−(1/24π²)ε_ijk Tr(L_iL_jL_k)` and its integral `K` — the numerical kernel is close to `fs-rep-mesh::winding_exact`'s solid-angle summation, but the map is the SU(2) texture, not surface normals; new code on the R4.1 field type.
- The **Hopf invariant** (linking of `n⁻¹(pt)` preimages) for the neutrino sector; the **disclination `π₁` winding** for electric charge. Both new.

### 4.4 The stochastic-QM sector (build new on `fs-rand`)
- An **SDE/Euler–Maruyama integrator** drawing Wiener increments `√dt·N(0,1)` from `fs-rand::next_normal`; an **osmotic-velocity** and **Fisher-information** functional; an **ensemble H-theorem** relaxation diagnostic using `fs-eproc` as the stopping rule; a **Madelung** `(ρ,S)⟷ψ` closure. None exist; this is a separate solver from the field PDE.

### 4.5 Certification-side extensions
- `fs-sos` currently certifies only **univariate SOS / 2×2 Lyapunov** and has **no SDP solver to *search* for a certificate** — so rigorous stability certification of a nonlinear high-dimensional soliton Hessian is out of reach today; the most one can ship now is `is_psd` on a computed Hessian (a numerical eigenvalue check, not a rigorous regional certificate). Extending `fs-sos` to a real SDP/Lasserre solver is a large, independent effort the workspace already lists as staged.
- Per-site manifold DOFs in `fs-ascent`'s runner (currently metadata + projected-gradient v1) would need a real engine to do Riemannian optimization of the director field cleanly (else fall back to augmented-Lagrangian with `N` pointwise `|φ_i|²=1` constraints).

---

## 5. A concrete, phased plan for a real GUM simulator

Realistic scope: simulate the **deterministic classical-field sectors** and certify them. This reproduces specific *within-model* computations (not the physics claims about nature).

**Phase 0 — new crate `fs-cosserat` (L3), the field type.**
Manifold-valued nodal field (`fs-ga` quaternions), discrete `L_i=Q̃⁻¹∂_iQ̃` on `fs-feec` incidence, relative texture `P̃`, dimensional tags via `fs-qty`. *Deliverable:* build a hedgehog ansatz `P̃=exp(if(r)x̂·σ)` and verify `∫b d³x = 1` (the paper's App. A.1 check) — the first certified number.

**Phase 1 — energy + soliton solver.**
Encode `W₂+W_χ+W₄+W₆₊₀` as a `Real`-generic energy (extend `fs-material`); minimize at fixed `K` with `fs-ascent` L-BFGS + `fs-adjoint` Sobolev-smoothed gradients + augmented-Lagrangian unit-norm constraints; Derrick virial check. *Deliverable:* the BPS compacton `f₀(r)=2arccos(r/R*)`, `C₆=64Λm̃/15π` by two independent quadratures — certified with `fs-evidence` colors + `fs-ivl` bounds. (Template: `fs-topopt` + `fs-metamat-e2e`.)

**Phase 2 — the linear spectrum.**
Linearize about the ground texture; assemble the transverse/longitudinal blocks; solve `ω(k)` with `fs-cheb` (the Orr–Sommerfeld generalized-eigen pattern) / `fs-la`; `fs-fft` for Floquet–Bloch bands. *Deliverable:* the four branches, a **gapless photon doublet** (`m_γ=0` to tolerance) and a gapped Klein–Gordon matter branch — certified as an eigenvalue enclosure.

**Phase 3 — the ħ-closure.**
Radial isorotation solve (`fs-cheb` BVP + `fs-ascent` AL): inertia `𝕀`, the spin+clock constraints, and the invariants `E_rot/E=¼`, `V=√2`; reject the superradiant rung by `fs-time`'s radiation logic. *Deliverable:* the closure table (`𝔠`, `κ`) reproduced and interval-certified.

**Phase 4 — dynamics (optional).**
Wire `fs-time` symplectic (`u`) + Lie-group (`Q̃`) integrators to the assembled coupled system; watch energy/charge conservation as the honest test. *Deliverable:* a stable propagating B2 mode and a static knot, with conserved `K` and bounded energy drift.

**Phase 5 — the campaign crate `fs-gum-e2e` (L4).**
Mirror `fs-metamat-e2e`/`fs-flutter-e2e`: one `run_campaign(...) -> GumReport`; per operating point emit `Color::Verified{lo,hi}` (else `Estimated{dispersion:∞}`); attach a `ModelCard` for the chiral-micropolar closure; bundle claims into an `EvidencePackage` so **`fs-checker`** re-verifies without running the solver. This is what makes it a *frankensim* result rather than a plot.

**Explicitly deferred (research-grade / out of scope):** the Nelson stochastic-QM sector (R9 — a whole separate SDE + Fisher-info engine), gravitation (R10), and the phenomenological spine (families/blue-fog, neutrino, confinement web), which in the paper are analytic/variational arguments with tuned inputs, not lattice simulations.

---

## 6. Honest verdict

- **Directly usable, and genuinely good:** the numerics (eigen/spectral/Krylov/FFT/Chebyshev), the structure-preserving + Lie-group integrators, the AD/adjoint/optimization stack, certified interval arithmetic and cross-ISA determinism, and the evidence/ledger/checker certification discipline. frankensim is, in fact, an *unusually well-suited* substrate for the deterministic parts — better than a from-scratch codebase, because the hard parts (symplectic+SO(3) integration, matrix-free adjoints, generalized-eigen BVPs, exact-`dd=0` DEC, deterministic replay, anti-laundering evidence) already exist and are tested.
- **The physics is absent and must be built:** there is no SU(2) *field*, no micropolar/couple-stress continuum, no topological-charge density, no chiral/Skyrme/sextic energy, and no stochastic-QM dynamics. These are not tweaks; they are the theory. The closest seeds are `fs-solid::rod` (Cosserat kinematics, 1-D static) and `fs-ga`+`fs-time::lie` (pointwise SU(2) + SO(3) stepping).
- **Scope realism:** a first `fs-cosserat`/`fs-gum-e2e` could credibly and *certifiably reproduce specific within-model numbers* — `∫b=1`, `C₆=64Λm̃/15π`, `E_rot/E=¼`, a gapless photon branch, a stable K=1 knot with conserved charge. That is a meaningful, honest deliverable. It would **not** — and no simulator could — validate GUM as a description of nature; the paper's own ledger already grades those claims within-model and self-scores the flagship at ≈2σ.
- **Effort estimate:** Phases 0–3 (the certifiable classical-field core) are a substantial but bounded project — on the order of a few new crates (`fs-cosserat` field type + energy, a topological-charge module, `fs-gum-e2e`) plus targeted extensions to `fs-solid`/`fs-material`/`fs-ascent`, leaning heavily on existing `fs-feec`/`fs-time`/`fs-cheb`/`fs-ascent`/`fs-adjoint`/`fs-evidence`. Phase 4 adds the dynamics; R9/R10 are separate research programs.

**Bottom line:** frankensim can't simulate GUM today — no crate implements the theory — but it is close to an ideal *foundation* to build one on: keep the numerics, the integrators, and the certification wholesale; write the SU(2)-field core, the micropolar energy, and the topological-charge machinery new; and stand up an `fs-gum-e2e` campaign that certifies the reproducible classical-field results the way the workspace certifies everything else.
