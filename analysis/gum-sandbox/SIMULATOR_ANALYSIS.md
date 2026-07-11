# Extreme-Detail Analysis: GUM-SWEET Research Simulator v3.5 + Explanation PDF

**Artifacts**

| File | What it is |
|------|-----------|
| `gum_sweet_research_simulator_v3_5_energy_pinning_good_20251130.py` | 241-line Python script (Colab export) |
| `...v3_5...good_20251130.ipynb` | **Byte-identical** to the `.py` (single code cell) |
| `..._EXPLANATION.pdf` | 27-page "GumCraft Code Components Addendum B" — a florid frame-by-frame narration of the video output, plus screenshots |

This is the **computational/visualization wing** of the same GUM program analyzed in `SUBSTRATE_SUITE_ANALYSIS.md` and `ANALYSIS.md` (branding: "GUM AI Foundation," "GumCraft," "GUM-SWEET TOY," "golden tuning V42.11"). The core question this analysis answers: **what does the code actually compute, and how does that relate to what the PDF says it computes?** The short answer: the code is a competent-but-quirky classical `φ⁶` scalar-field toy; the PDF narrates it as topological particle physics it does not contain.

---

## 1. What the code actually is

A single **real scalar field** `φ(x,y,z)` on a **128³ periodic grid** (box half-width `L=12`), evolved by a classical wave/Klein–Gordon-type equation with a **sextic self-potential**, plus an ad-hoc **"energy pinning"** rescale, rendered as a 2×2 matplotlib dashboard to an `.mp4`.

### 1.1 The field equation
- **Potential** (lines 44, 56): `V(φ) = Aφ² + Bφ⁴ + Cφ⁶` with `A=−0.15, B=0.8, C=0.2` — a generic Landau sextic with a slightly unstable origin (`A<0`).
- **Potential force** `dV/dφ = 2Aφ + 4Bφ³ + 6Cφ⁵` (`force_field_active`, line 45), plus an exponential **barrier** `2·exp(2(φ−4.5))` for `φ>4.5` (lines 47–50) that hard-stops runaway growth.
- **Laplacian** (lines 165–168): standard 7-point stencil `Σ(6 neighbours) − 6φ` via `np.roll` (periodic BCs). **Grid spacing is not included** — the stencil implicitly uses `Δx=1`, whereas the physical spacing is `2L/(N−1)=24/127≈0.19`. So the "Laplacian" is in grid units, decoupled from `L`.
- **Total force** (lines 170–174): `f_elastic = LAMBDA·lap` (`λ=0.89`), `f_shear = MU·0.01·lap` (`μ=42 → 0.42·lap`), `f_potential = −dV/dφ`. Acceleration `= total_force · inertia_scale`, `inertia_scale = 1.230e-3/1.317e-3 ≈ 0.934` (a constant).

### 1.2 The integrator (this part is legitimate)
`update()` runs `SUB_STEPS=20` velocity-Verlet (kick–drift–kick) steps per frame (lines 181–186):
```
acc,lap = get_acceleration(phi); pi += 0.5*acc*DT
phi += pi*DT
acc,lap = get_acceleration(phi); pi += 0.5*acc*DT
```
Velocity-Verlet is a correct symplectic integrator; with `DT=5e-4` and effective wave speed `c²≈(0.89+0.42)·0.934≈1.22`, the CFL number `c·DT/Δx≈5.5e-4 ≪ 1` — heavily over-resolved in time, so the scheme is very stable and would **not** blow up from CFL. The run will settle into a saturated blob rather than diverge.

### 1.3 Initialization
- **"Alpha particle (K=4)"** (lines 73–83): four amplitude-4.5 Gaussians at tetrahedral centres `(±1.2,±1.2,±1.2)` (even sign products).
- **"Canary Neutrino (K=1)"** (lines 86–88): one amplitude-4.2 Gaussian at `(−6,0,0)`.
- `phi = np.clip(phi, 0, 5.2)` (line 90) — a hard amplitude cap.
- **Warm-up** (lines 106–116): 200 steps of a *different*, heavily damped integrator (semi-implicit Euler with `pi *= 0.90` every step), after which `TARGET_ENERGY` is captured from the relaxed state (lines 120–122).

---

## 2. The "energy pinning" — dissected

After every Verlet substep (lines 188–198):
```
E = E_kin + E_g + E_pot
if E > 1.0:
    s = (TARGET_ENERGY / E) ** 0.001
    phi *= s ;  pi *= s
```
This is the file's headline feature ("Eliminates numerical dissipation by pinning Total Energy"). Three things are worth stating plainly:

1. **It is not physics.** A symplectic integrator already conserves energy to high accuracy; genuine conservation is *achieved*, not *imposed*. This code instead multiplicatively rescales the field toward a target number each step — a homeostat that **masks** drift rather than removing its cause.
2. **A single scale factor cannot consistently restore the energy**, because the components scale differently under `φ→sφ, π→sπ`: `E_kin∝s²`, `E_grad∝s²`, but `E_pot` mixes `s²` (quadratic), `s⁴` (quartic) and `s⁶` (sextic) pieces. So "pinning total E" deforms the field in a way that has no correspondence to the dynamics. (The `0.001` exponent makes each nudge tiny, so it is mostly cosmetic — but cosmetic in exactly the direction that produces the "locked" appearance.)
3. **The energy being pinned is not the integrator's conserved Hamiltonian anyway.** `E_g = −0.5 Σ φ·lap` uses coefficient 1, but the dynamics use `1.31·lap` as the gradient force; and `E_kin=0.5 Σ π²` assumes unit mass while the acceleration carries `inertia_scale≈0.934`. The reported "Hamiltonian components" are therefore a *different* quantity from what the equations of motion conserve. The on-screen "perfectly flat, no drift" total is a number held near target by fiat, not a conserved physical energy.

**Net effect:** the visible "lock at φ≈4.2 / stable core-halo" is produced by the **`clip(0,5.2)` cap + the `φ>4.5` exponential barrier + the pinning homeostat**, i.e. the field is squeezed against an amplitude ceiling. It is amplitude clamping dressed as "42 resonance."

---

## 3. Correctness / issues table

| # | Item | Finding |
|---|------|---------|
| 1 | Velocity-Verlet loop | **Correct** symplectic integrator; stable at this `DT`. |
| 2 | Discrete Laplacian | Correct 7-pt stencil, but **ignores `Δx`** (grid units ≠ physical units). |
| 3 | "Elastic" vs "shear" | `f_shear = 0.42·lap` is **mathematically identical** to the elastic term — a single scalar field has no shear tensor. "Shear (MU=42)" is a *label* on `0.42×` the same Laplacian, no distinct physics. |
| 4 | Energy pinning | Non-physical; inconsistent across components; masks rather than conserves (§2). |
| 5 | Energy accounting | Coefficients don't match the dynamics' Hamiltonian → the "conserved" total isn't the dynamical energy. |
| 6 | `inertia_scale`, `RHO_0` | Density-flavoured names for a constant `≈0.934`; two different densities (`1.317e-3`, `1.230e-3`) appear with no derivation. |
| 7 | `N=128` hard-coded (line 62) | Contradicts the header comment "64³ Grid (Scalable to 128)"; 128³ float32 × 600 frames × 20 substeps × ~12 `np.roll`/substep is heavy (many minutes to hours) and memory-hungry. |
| 8 | `meshgrid` default `indexing='xy'` | Axis-0 is Y, not X; the canary probe index (`phi[canary_idx_x, N//2, N//2]`) is therefore slightly off the intended point (harmless for symmetric seeds). |
| 9 | Dead code | `import time` unused; `history_time`/energy lists grow unbounded. |
| 10 | Output | Requires `ffmpeg`/`FFMpegWriter`; wrapped in try/except so a missing codec prints a failure rather than crashing. |
| 11 | Topological charge | **There is none.** A real scalar `φ:ℝ³→ℝ` supports no `π₃` winding. "K=4"/"K=1" are labels on Gaussian bump-counts, not computed invariants. |

**Verdict on the code as software:** it runs, it's numerically stable, and the integrator is sound. As *a physics artifact* it is a generic damped/clamped `φ⁶` scalar blob-former, with two structural defects (pinning is non-physical; the reported energy isn't the conserved one) and a raft of decoratively-named terms that carry no independent physics.

---

## 4. The PDF explanation vs. the code — a large gap

The PDF ("Addendum B") narrates the video output in lavish GUM-theoretic language. Cross-checking each major claim against the actual code:

| PDF claim | Reality in the code |
|-----------|---------------------|
| "42 resonance… ω_torsion = √(2μ/ρ₀) = 42… measured ω=42.004±0.002 universal lock" | **No frequency is ever computed** — no FFT, no `ω`. And the quoted formula gives `√(2·42/1.317e-3) ≈ 252`, not 42. The "42" is simply the input constant `MU=42` (and `φ≈4.2` is the clip/barrier ceiling and the canary seed amplitude). It's a numerological motif ("V42.11", the Hitchhiker's 42), not a measurement. |
| "FCC lattice core / hedgehog texture / spin currents / torsional currents / gyroscopic Magnus G×v" | None exist. A single real scalar has no spin, no torsion, no vector current, no lattice crystal (the cubic *mesh* is the numerical grid, not a physical FCC solid). |
| "Alpha (K=4) breathes 4× stronger than neutrino (K=1) — vibron hierarchy / topological charge" | The alpha seed is **4 superimposed Gaussians**, the canary is **1** — so the central value is trivially larger. "4× because K=4" is 4 bumps vs 1 bump relabeled as topological charge. |
| "Sextic growth detected (120.3), DT modulated to 2.5e-4 / prints every 10 steps / coherence metric / dt→1e-7 at spikes" | **The attached code has fixed `DT=5e-4`, no adaptive timestep, no sextic-growth detector, no coherence metric, and prints every *frame* in a different format.** These logs describe a *different, hypothetical* script the PDF calls "Mods Planned"/"added aggressive observability" — not this file. The "run output" is not this code's output. |
| "20 GeV sub-harmonic (Fermi match) / scale-invariant from N=16 to 1.24M" | Unsupported; there is one `N=128` run, no scan, no energy-scale calibration. |
| "No diffusion, pure determinism, superfluid halo" | Every classical PDE is deterministic; "superfluid" is decoration. |
| "Pinning conserves E, no drift — proving the pinning magic tames the sextic beast" | The total is held near target *by rescaling*, and it isn't the dynamical energy (§2–3). "Conserved with no drift" is what the homeostat is built to display. |

**This is pareidolia**: a rich physical mythology read into a blob that saturates against an amplitude cap. The narration is confident, specific, and — measured against the code — substantially fabricated (most sharply, the run logs and the `ω=42.004±0.002` "measurement").

---

## 5. Relationship to the GUM corpus

Even taken on the program's own terms, **this simulator does not implement GUM.** The rigorous paper's topological/mass sector is built on an **SU(2)-valued micropolar texture** `P̃(x)` with a Skyrme term and a Bogomolny bound, and its sextic `W₆₊₀ = ½Λ²b² + m̃²(1−σ_P) + c₂(1−σ_P)²` is a potential on `σ_P = ½Tr P̃` — *not* a polynomial in a single real scalar. Particles there are topological knots with a genuine `π₃(S³)` degree `K`.

This code has **none of that structure**: one real scalar, a Landau `Aφ²+Bφ⁴+Cφ⁶`, Gaussian lumps with `K` written on them by hand. It is, at most, a *cartoon of the sextic-potential motif*, missing the SU(2) texture, the topology, the couple-stress/chirality, and therefore every feature (knots, spin-½, confinement, the family logic) that the paper's physics actually turns on. So the simulator sits to the corpus roughly as the earlier WebGL sandbox does: an evocative visual with GUM vocabulary layered on generic field dynamics — one more rung below even the "older variant" teaching tool, because here the narration actively over-claims measurements the program's own discipline (ledger tags, "no external sentence exceeds its license") would forbid.

---

## 6. Bottom line

- **As code:** a runnable, numerically stable `φ⁶` scalar-field toy with a correct Verlet core, wrapped in a non-physical "energy-pinning" homeostat and several cosmetically-named terms (a redundant "shear," a constant dressed as density). It will produce a pretty core-halo blob video. It is not wrong to *run*; it is wrong to *interpret* as the PDF does.
- **As physics:** it computes a clamped classical scalar field. There is no topological charge, no resonance measurement, no spin/torsion, no FCC crystal, no scale-invariance study — all of which the PDF asserts. The "42 lock" is `MU=42` + a `φ≈4.2` amplitude ceiling; the "conservation" is imposed by rescaling.
- **Most important single finding:** the PDF's "run output" (adaptive-dt logs, sextic-growth detector, `ω=42.004±0.002`, coherence metric) **does not correspond to the attached code**, which has fixed `DT`, no FFT, and no such instrumentation. The narrative describes a different or imagined program and presents its results as measured facts.

### Constructive notes (if the goal is a *real* toy)
1. Drop energy pinning; a symplectic integrator conserves energy on its own — then plot the *true* Hamiltonian and show the drift honestly (that is the actual test of the scheme).
2. Make the reported energy match the dynamics (same gradient coefficient; include `Δx`; account for `inertia_scale` as an effective mass in `E_kin`).
3. If topology is wanted, the field must have the right target space (e.g. an SU(2)/`O(3)` multiplet) and `K` must be *computed* from a winding integral — otherwise "K=4" is a caption, not a charge.
4. If "resonance at ω" is a claim, *measure* it: record `φ(t)` at a point and FFT it. Right now no frequency is computed anywhere.
5. Keep the vocabulary matched to what the code does — per GUM's own ledger discipline, "no external sentence may exceed its license."
