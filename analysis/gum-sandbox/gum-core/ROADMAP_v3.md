# GUM Physics Core — Roadmap v3 (the engineering program)

Successor to `GUM_SIMULATION_ASSESSMENT_v2_FABLE.md`'s tier ladder: the
replication campaign (Tiers 0–5a) is complete; this roadmap builds the
**fs-gum physics core** on frankensim numerics, closing the five gaps of
`GAP_ANALYSIS_v3.md`. Design rules inherited from the campaign: every
crate ships with machine-checked gates against frozen campaign referee
numbers; deterministic iteration order and bit-identical replay; the
fs-evidence/fs-package/fs-checker certification pattern (4B template);
standalone crates under `analysis/gum-sandbox/gum-core/` (empty
`[workspace]`, path deps) so the fork's workspace is untouched.

## Phase B1 — `fs-gum-field` (G2 + G4's measurement path) ▶ EXECUTING
Stored SU(2) field foundation: `Field3` cell-centered N³ quaternion field
(SoA, ghost rind, vacuum fill, renormalization), central-4th + corner
stencils, Maurer–Cartan currents, sector energies E2/E4/E6/E0 + inertia I
+ degree on stored fields, hedgehog seeding from a radial profile, and
the in-crate 1-D radial solver (graded grid, damped Newton, ε-dial).
Gates: 4B's certified compacton enclosures; field3d's ε=0.05 hedgehog
values (E2 = 12.117459 −0.34%-class FD accuracy documented, E0/I at
1e-8-class, degree 0.9999, virial); Yukawa μ = 1/√(2t); bit-replay.

## Phase B2 — `fs-gum-topo` (G3) ▶ EXECUTING
Director projection n = R(q)ẑ; degree (central + corner + VOS geometric
variant); **Hopf invariant** by both methods (Whitehead with fs-fft
Coulomb-gauge curl-inverse incl. zero-mode/obstruction checks; preimage
linking via marching tetrahedra + Gauss pairs) cross-validating each
other; **disclination Z/2 loop holonomy** + half-integer framed winding +
dual-plaquette defect-line tracing (lines cannot end — built-in check).
Gates: hedgehog degree 1; an analytic Hopf-1 texture; a +1/2 disclination
line; method-vs-method agreement.

## Phase B3 — `fs-gum-cosserat` (G1) ▶ EXECUTING
The micropolar linear-dynamics module: (u, φ) symbol assembly for
arbitrary k (all seven energy terms; Case A objective mass vs Case B —
the M-1 defect — and χ₃), Hermitian 6×6 via real-symmetric 12×12
embedding on fs-la's jacobi_eigh, generalized reduction, branch
classification + co-motion diagnostics, and a mass-matrix-aware
symplectic Verlet giving the first time-domain validation of the coupled
second-order system. Gates: the full Tier-1 referee battery (c_L to
1e-12-class, B2 gap ≤ 1e-12, ω₀² = 21 with factor-4 dial, Case-B
c_B2² triples, r = 4μc/(4μc+m_V²) to 9 digits, χ₃ k⁵ law, rotational
invariance of the arbitrary-k spectrum, energy-bounded long Verlet run).

## Phase B4 — `fs-gum-sde` (G5) ▶ EXECUTING
`WienerStream` on fs-rand Philox with the (trajectory, step, component) →
draw-index contract (order-independent ensembles, per-trajectory replay);
Euler–Maruyama + additive-noise SRA1 (with ΔZ); OU analytic gates
(moment closure; measured strong/weak convergence orders); mode-sum ψ
machinery ported from born.py + the **osmotic velocity**; reflecting
boundaries; the Nelson H-theorem demo: ρ → |ψ|² relaxation with
coarse-grained H̄(t) against Tier-3's qualitative referee (M=16 relaxes
to the noise floor; the equilibrium ensemble stays there).

## Phase E — `fs-gum-statics` (G4's solver half) — SPEC FROZEN, NEXT
The medium job (~5–7 kLOC, weeks-class): analytic variational gradients
for every sector (flux-form (2+4), sextic cofactor, E0, Routhian term),
scheme adjoints, arrested Newton flow with the documented near-BPS guards
(corner objective, degree anchor, Bogomolny wall — parameters frozen in
field3d_RESULTS.md), the W_χ module (χ couplings as free parameters,
default 0; Dzyaloshinskii margin 𝔪 as the text-pinned referee), and the
full gate suite (Derrick virial, clock bisection, FD-vs-analytic gradient
at 1e-5). Deferred past this session's phase-B wave; its complete spec is
in `gap_survey_v3.json` (dimension 4). With Rust speed, the N=96 gate
suite becomes CI-runnable (~10 min vs 75 min Python).

## Phase F — integration & program services (after E)
- e2e certified demo: seed hedgehog → relax (E) → isorotate + clock →
  diagnostics (B1/B2) → EvidencePackage golden root.
- The F-R5 service run: the corpus-facing deliverable — a reproducible,
  certified demonstration package of the halo instability + saturated
  closure, for the corpus's authors to run against their frozen ⟨r1⟩
  code (the discharge path of F-R4/F-R5).
- Out of envelope, unchanged: gravitation sector, blue-fog GPU boxes,
  nucleation (research-class; revisit after the corpus responds).

## Execution record
- Phase A (gap analysis workflow, 5 parallel surveyors): ✅ this document
  + GAP_ANALYSIS_v3.md + gap_survey_v3.json.
- Phases B1–B4: launched as four parallel build agents; each lands with
  its gate table and RESULTS.md; adjudication + commit per phase.
- Phase E/F: specs frozen here; execution follows the B-wave.
