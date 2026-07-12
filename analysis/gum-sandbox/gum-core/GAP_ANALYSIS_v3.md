# GUM Physics Core — Gap Analysis v3 (post-campaign)

Five-dimension parallel survey (full structured output: `gap_survey_v3.json`)
re-auditing the v1/v2 gap list against what exists **after** the replication
campaign. Headline: the original claim ("a real GUM simulator is a new
physics core, not a reuse of any existing crate") remains true for *Rust
crates* — but it is now misleading as a statement about the *repository*:
every physics algorithm the core needs exists in **validated, referee-locked
Python** under `analysis/gum-sandbox/`, and the build has become an
engineering port with frozen specs, not research.

## Verdict per named gap

| # | Gap (as claimed) | Verdict now | What closed it / what remains |
|---|---|---|---|
| G1 | No micropolar/Cosserat continuum | **CONFIRMED** for crates (fs-solid is displacement-only Cauchy; the only "Cosserat" is the 1-D rod). But the full linearized micropolar dynamics exists validated in `tier1-spectrum/spectrum.py` (symbol matrices, all seven energy terms, Case A/B mass, χ₃ chiral, four branches to 9+ digits). Remains: Rust port + complex-Hermitian eigensolve path (fs-la lacks Hermitian-with-vectors) + mass-aware symplectic Verlet (fs-time is unit-mass) + arbitrary-k assembly. ~1.2–1.8 kLOC. |
| G2 | No SU(2)-valued field type | **CONFIRMED with refinement**: fs-ga quats are single-value; fs-cosserat-pilot (Tier 4B) evaluates *analytic* configs with zero storage. The stored-field machinery (SoA layout, ghost rind, renormalization, 4th-order + corner stencils) exists only in `field3d_solve.py`. Remains: the Rust `Field3` type + stencils + deterministic-iteration contract. ~1.2–1.8 kLOC. |
| G3 | No topological-charge density on a field | **TWO-THIRDS FALSE**: degree/b_P is computed twice (4B certified O(h²); field3d on stored fields). Genuinely missing: **Hopf invariant** (both Whitehead/FFT-Coulomb and preimage-linking methods), **disclination winding** (π₁(RP²) = Z/2 loop holonomy + defect-line tracing), and any S²/RP² *director* representation (the Hopf projection n = R(q)ẑ is never computed). ~1.8–2.6 kLOC. |
| G4 | No Skyrme/chiral energy terms | **FALSE for the campaign, TRUE for crates**: E2/E4/E6/E0, arrested Newton flow, near-BPS guards, radial/axisymmetric/3-D engines — all validated Python with frozen referee numbers. The **W_χ micropolar coupling is genuinely absent everywhere** (and its χ₁/χ₂/χ₃ couplings are *unpinned by the corpus text* — must ship as free parameters, default 0, with the Dzyaloshinskii margin 𝔪 as the only text-pinned referee). The Rust port incl. analytic gradients is the one **medium** job: ~5–7 kLOC, weeks-class. |
| G5 | No Nelson/stochastic-QM sector | **CONFIRMED, with one precision**: fs-rand has the right primitives (Philox counter streams, fixed-2-draw normals — random-access, replayable) but no Wiener abstraction, no SDE integrator (zero grep hits for Maruyama/Milstein/Brownian), and the osmotic velocity is computed nowhere (born.py does only the current velocity). Physics note from the survey: Nelson noise is additive-constant, so **Milstein degenerates to Euler–Maruyama exactly** — the higher-order scheme must be an additive-noise SRK (SRA1-class, needs the correlated ΔZ draw). ~1.8–2.6 kLOC. |

Gravitation (Einstein–Cartan defect geometry) and the blue-fog condensate
remain outside this analysis — theorem-level / GPU-scale respectively, per
the v2 assessment's Tier-5 classification.

## The structural conclusion

The campaign inverted the build economics: in v1/v2 the physics core was
the *risk*; now every algorithm, constant, guard parameter, and gate value
is frozen in referee-locked artifacts, and the risk is reduced to
**transcription fidelity + bit-determinism discipline** — exactly what the
4B pilot proved the evidence stack can certify. The port order below
(ROADMAP_v3.md) is chosen so each crate lands with machine-checked gates
against campaign numbers on day one.
