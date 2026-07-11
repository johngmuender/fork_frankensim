# External open-source code as a base for a GUM simulator

**Question.** Given the GUM physics gaps that frankensim does not fill — (G1) a micropolar/Cosserat continuum with independent micro-rotation DOF + couple-stress; (G2) an SU(2)/S²-valued orientation *field*; (G3) topological-charge density on a field (`b_P`, Hopf, disclination winding); (G4) Skyrme/chiral energy terms; (G5) a Nelson/stochastic-QM sector — are `stark_micropolar` or `hoomd-blue` better bases, could they be incorporated, and what else exists?

**Short answer.** Neither named repo is a whole-theory base, but they sit on opposite sides of the gap and only one is relevant. **`stark_micropolar` is a genuine hit for G1** (real Cosserat micropolar mechanics with micro-rotation + couple-stress, differentiable, permissive license) and a strong *reference* for the elastic sector. **`hoomd-blue` is the wrong paradigm** (particle MD/MC, not a continuum field theory) and is not a base. The parts neither covers — G2/G3/G4 (the topological-soliton/chiral-field half) — are exactly what the **micromagnetics** ecosystem (`mumax3`/`Ubermag`/OOMMF) and modern **Skyrme-model soliton solvers** (`cuSkyrmion`, `soliton_solver`) already do, and those are the closest existing tools for that half. GUM is explicitly assembled from these very literatures (Cosserat elasticity + chiral liquid crystals/DMI + Skyrme solitons), so their codes are the natural cribs. No single project does the *coupled* theory; a real simulator remains an integration effort, but far less "from scratch" than it looked.

---

## 1. The two named repos

### `InteractiveComputerGraphics/stark_micropolar` — a real micropolar continuum (fits G1)
- **What it is:** the reference implementation of *"Curved Three-Director Cosserat Shells with Strong Coupling"* (Löschner, Fernández-Fernández, Jeske, Bender, CGF 2024), extending the **`stark`** simulation framework. It models **three-director Cosserat shells** — a micropolar continuum with **independent micro-rotational DOF** and **couple-stresses intrinsic to the formulation**, strongly coupled to translation.
- **Method:** FEM (Tri3/6/10, Quad4/9) in an **incremental-potential** formulation, **implicit** strongly-coupled solve, **differentiable** via runtime derivative compilation (`symx`). C++ / Eigen / OpenMP, optional MKL. **Apache-2.0.** Key file: `stark/src/models/deformables/surface/EnergyMicropolarShell.cpp`.
- **Fit vs GUM:**
  - **G1 micropolar/couple-stress: strong.** This is exactly the class of mechanics fs-solid lacks — independent orientation DOF, asymmetric stress, couple-stress, a hyperelastic energy expressed variationally. The differentiable incremental-potential design is also close in spirit to frankensim's AD-energy pattern.
  - **G2 orientation field: partial.** It carries **director frames** per element (an SO(3)-like micro-rotation), which is 90% of what an SU(2) texture field needs — but the target is the shell's director triad, not a bulk `Q̃:ℝ³→SU(2)`.
  - **Limitations:** it is a **shell** code (2-D surfaces in 3-D), not a 3-D **bulk** continuum; graphics/animation-oriented (robustness over spectral accuracy); **no topological charge, no chiral transduction, no Skyrme/sextic energy** (G3/G4 absent). The parent `stark` does 3-D deformables+rigids strongly coupled, so the natural move is "use `stark`'s 3-D FEM + this micropolar energy pattern, extend to a 3-D micropolar bulk with the GUM energy."
- **Verdict:** the **best existing base/reference for GUM's elastic (u + micro-rotation) sector** — worth reading closely and plausibly *incorporating* the micropolar energy assembly. Not a topological-soliton or field-theory engine.

### `glotzerlab/hoomd-blue` — particle MD/MC (wrong paradigm)
- **What it is:** a mature, GPU-accelerated **particle simulator** — molecular dynamics + hard-particle Monte Carlo of anisotropic shapes, Python front end over C++/CUDA, extensible custom potentials. **BSD-3.**
- **Fit vs GUM:** GUM's core is a **continuum field theory** (`u(x)`, `Q̃(x)` on a mesh), not an N-particle system. hoomd-blue **solves no PDE and computes no field topological charge**. Its per-particle orientation/quaternion + rigid-body support and GPU infrastructure are real, but they serve particle ensembles, not a director field.
- **Verdict:** **not a base.** Relevant only if one deliberately reframed GUM as a *knot-particle ensemble* (a coarse caricature, not the field theory), or wanted to borrow GPU/orientation-particle infrastructure. Skip for the field-theoretic core.

---

## 2. The tools that cover the half neither named repo does (G2/G3/G4)

GUM's paper cites **Dzyaloshinskii condensation, blue phases, heliknotons (Ackerman–Smalyukh), and chiral magnets** as its "primitives of discovery." Those are *precisely* the systems micromagnetics and chiral-liquid-crystal codes simulate. That makes the following the closest existing tools for the topological/chiral sector:

### `mumax3` (+ `OOMMF`, `Ubermag`, `Fidimag`) — micromagnetics = the S²-director/topological/chiral sector
- **What it is:** GPU-accelerated micromagnetic simulator. Evolves a **unit-vector field `m:ℝ³→S²`** on a finite-difference grid under Landau–Lifshitz dynamics, with an energy of exchange + **Dzyaloshinskii–Moriya interaction (DMI)** + anisotropy + demagnetization.
- **Why it maps onto GUM astonishingly well:**
  - **G2 director field:** `m∈S²` per cell — exactly the blue-fog director `n̂∈S²` of GUM's collective sector.
  - **G3 topological charge:** mumax3 computes the **skyrmion/winding number** (Berg–Lüscher lattice topological charge) natively, and the **Hopf number** for 3-D hopfions — i.e. GUM's `π₃(S²)` neutrino charge and the `π₂`/`π₃` invariants, already implemented.
  - **G4 chiral energy:** **DMI *is* the Dzyaloshinskii chiral coupling** GUM's `W_χ` and the blue-fog condensation criterion are built on. Simulating skyrmions/hopfions in a chiral magnet is simulating GUM's own analog system.
- **Limitations:** it is a **pure S² sigma model** (fixed `|m|`), not coupled to a displacement field `u` or a micropolar elastic continuum; dynamics are **dissipative LLG** (relaxes to minima; good for *finding* solitons and charges, not the paper's symplectic Hamiltonian evolution or the massless-photon spectrum); it is an `S²` theory, not the full `SU(2)` texture. License GPL-3.0 (mumax3), Go+CUDA; OOMMF is NIST public-domain (C++); **`Ubermag`** wraps OOMMF/mumax3 in Python with a weak-form-ish DSL and is the most *extensible/scriptable* entry point.
- **Verdict:** the **single closest existing tool for GUM's topological-charge + chiral + director-soliton half (G2/G3/G4)** and the **fastest path to a first real result** — a chiral soliton (skyrmion/hopfion) with a *computed* topological charge — because you configure rather than build. It does **not** provide the elastic/spectrum half.

### `cuSkyrmion` and `soliton_solver` — direct Skyrme/nonlinear-field soliton solvers (G3/G4/R5)
- **`cuSkyrmion`** (2026, CUDA/OpenGL): a **3-D Skyrme-model** solver via **arrested Newton flow** — it minimizes the actual `SU(2)` Skyrme energy to static solitons and tracks **baryon number**. This is GUM's `W₄` Skyrme term + `π₃(S³)` degree, done exactly.
- **`soliton_solver`** (2026, Numba-CUDA): a **theory-agnostic** GPU finite-difference solver for 2-D nonlinear field theories with a **modular plug-in "theory registry"** — designed so you *add your own Lagrangian/energy*. That extensibility is the right shape for coding GUM's `W₂+W_χ+W₄+W₆` and watching a topological soliton form.
- **Limitations:** research-grade, niche, single-group; `soliton_solver` is 2-D; both are sigma-model/Skyrme, not the coupled micropolar-elastic theory. But for the **Skyrme energy + topological degree** sub-problem they are near-exact references.
- **Verdict:** the best **algorithmic references (and possible plug-in hosts)** for GUM's `SU(2)` Skyrme/topological sector.

### `FEniCS` / `Firedrake` — the "build any weak form" substrate (G1+G2+G4, most flexible)
- **What they are:** automated FEM from a symbolic weak-form DSL (UFL). You write the energy/weak form; the framework generates and solves it. Firedrake adds strong **FEEC** support; both have published **Cosserat micropolar / strain-gradient elasticity** implementations.
- **Fit vs GUM:** the most flexible base for the **full coupled `u`+`Q̃` micropolar field theory** — you can encode asymmetric strain, wryness, couple-stress, the chiral cross terms, and even a manifold-constrained director field as custom weak forms, and get MPI-parallel solves. Mature, Python, LGPL-class.
- **Limitations:** you supply *all* the physics (no micropolar/Skyrme/topological-charge out of the box); topological charge and the SU(2)-field constraint are hand-coded; not GPU-first. It is a substrate, not a solution — but a *very* capable one.
- **Verdict:** the best choice if the goal is a **single framework that can express the entire coupled action** and you accept writing the physics in UFL.

---

## 3. Mapping: GUM gap → best existing open-source code

| GUM gap | Best existing base | Fit | What's still missing / caveat |
|---------|-------------------|-----|-------------------------------|
| **G1** micropolar/Cosserat continuum, couple-stress, u+micro-rotation | **stark_micropolar / stark**; or **FEniCS/Firedrake** | strong | stark is *shells* (extend to 3-D bulk); Firedrake needs the energy written in UFL |
| **G2** SU(2)/S²-valued orientation *field* | **mumax3/Ubermag** (S² field) ; stark (SO(3) directors) | strong (S²) / partial (SU(2)) | S² ≠ full SU(2) texture; no coupling to `u` in mumax3 |
| **G3** topological-charge density (skyrmion #, Hopf #, winding) | **mumax3** (native); **cuSkyrmion** (baryon #) | strong | disclination `π₁` winding still custom; on a *coupled* field, new code |
| **G4** Skyrme + chiral energy | **cuSkyrmion** (Skyrme), **mumax3** (DMI=chiral), **soliton_solver** (plug-in) | strong | assembling the *combined* `W₂+W_χ+W₄+W₆` with `u` is new integration |
| **G5** Nelson/stochastic-QM (SDE + Fisher info) | **none** (no standard OSS) | absent | custom, as before (build on an SDE integrator + `fs-rand`) |
| dynamics: symplectic + Lie-group | **frankensim `fs-time`** | strong | mumax3/stark dynamics are dissipative/implicit, not symplectic |
| **certification / determinism / units** | **frankensim** (`fs-evidence`/`fs-qty`/`fs-checker`/`fs-ledger`) | unique | none of the external tools have this |

---

## 4. Recommendation

**No single repo is a drop-in base — and that is expected**, because GUM deliberately fuses four literatures (Cosserat elasticity, chiral liquid crystals/DMI, Skyrme solitons, stochastic mechanics). The pragmatic answer is a **per-sector, crib-and-integrate** strategy rather than adopting one base:

1. **Fastest real result (topological/chiral sector):** start with **Ubermag/mumax3**. It already has the S² director field, DMI (the paper's own Dzyaloshinskii chiral coupling), and native skyrmion/Hopf topological-charge diagnostics — you can produce a chiral soliton with a *computed* charge by configuration, not construction. This directly exercises G2/G3/G4 and mirrors GUM's own analog systems (heliknotons, blue phases, chiral magnets).
2. **Elastic (u + micro-rotation) sector:** use **`stark`/`stark_micropolar`** as the reference/base for micropolar mechanics with couple-stress (G1), or **FEniCS/Firedrake** if you want one framework to hold the *entire* coupled weak form.
3. **Skyrme energy + baryon number:** crib **`cuSkyrmion`** (arrested Newton flow, `SU(2)` Skyrme, degree) and the plug-in design of **`soliton_solver`**.
4. **Nelson stochastic-QM (G5):** still custom — no established OSS; build an Euler–Maruyama SDE integrator + Fisher-information functional (on `fs-rand`'s Wiener increments, or NumPy/JAX).
5. **Keep frankensim for what only it does:** the **certification wrapper** (evidence colors, dimensional units, deterministic cross-ISA replay, the solver-free `fs-checker`) and its **symplectic + SO(3) integrators** (`fs-time`) — none of the external tools provide these.

**Honest caveat on "incorporate":** these projects span C++ (stark, OOMMF), Go+CUDA (mumax3), Python (Ubermag, FEniCS/Firedrake, soliton_solver), and Rust (frankensim), under mixed licenses (Apache-2.0, BSD-3, GPL-3.0, LGPL, public-domain). "Incorporate" therefore means **use as a component in a pipeline or crib the algorithms/energy assembly**, not merge codebases. And crucially: each of these simulates a *piece* of GUM's borrowed analogy — none simulates the **coupled** action with the relative texture `P̃=R̃[u]†Q̃` binding elasticity to the SU(2) sector. That coupling, plus G5, is the irreducible new work; but with mumax3 for the chiral/topological half, a micropolar engine (stark or FEniCS) for the elastic half, cuSkyrmion/soliton_solver as Skyrme references, and frankensim for certification, the "build from scratch" surface shrinks dramatically.

### Other projects worth a look (by sector)
- **Chiral liquid crystals / heliknotons:** the Smalyukh-group nematic/`Q`-tensor and blue-phase relaxation codes (the literal experimental system GUM's neutrino/blue-fog sector is modeled on).
- **Nematic `Q`-tensor FEM:** general Landau–de Gennes solvers (some in FEniCS/Firedrake) for the tensor-order-parameter version of the director field.
- **Skyrme/nuclear:** Manton–Sutcliffe-lineage soliton relaxation codes; lattice `O(3)`/`CP¹` sigma-model solvers.
- **Cosserat rods/shells beyond stark:** PyElastica (Cosserat rods) and discrete-elastic-rod libraries (1-D, like fs-solid::rod — lower relevance).
- **Differentiable physics substrates:** JAX-MD / Taichi / Warp if a GPU, autodiff-first custom field solver is preferred over adapting a domain code.

---

## 5. Bottom line
- **`stark_micropolar`: yes, a genuine and useful base/reference for the micropolar-elastic gap (G1)** — permissively licensed, differentiable, the right mechanics; extend from shells to a 3-D micropolar bulk with the GUM energy.
- **`hoomd-blue`: no** — particle MD/MC, wrong paradigm for a continuum field theory.
- **The strongest *additional* find is micromagnetics (`mumax3`/`Ubermag`)**, which already implements GUM's director field, chiral (DMI) coupling, and topological-charge computation — the closest existing tool for G2/G3/G4 and the quickest route to a first certifiable chiral-soliton result — complemented by `cuSkyrmion`/`soliton_solver` for the Skyrme energy and `FEniCS`/`Firedrake` for a single-framework coupled weak form.
- **frankensim remains the certification and deterministic-numerics layer** around whichever physics engine(s) you adopt. The realistic architecture is a **polyglot pipeline**: an external engine (or engines) computes the field/soliton/spectrum; frankensim wraps it in evidence, units, and a solver-free checker.
