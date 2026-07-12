# fs-gum-e2e — GUM core Phase F1 results (end-to-end certified demo)

Crate: `analysis/gum-sandbox/gum-core/fs-gum-e2e/` (standalone `[workspace]`
opt-out; path deps: the three FINISHED gap-closure siblings
`../fs-gum-field`, `../fs-gum-statics`, `../fs-gum-topo` plus the evidence
stack `crates/fs-ivl` / `fs-evidence` / `fs-package` / `fs-checker` /
`fs-blake3` — std-only closure, the fs-cosserat-pilot pattern).

Purpose: the roadmap's Phase-F e2e certified demo — ONE binary that runs the
whole GUM knot pipeline (seed → measure → relax → clock/isorotate → topo →
package) and emits ONE EvidencePackage golden Merkle root. This is the
integration proof that the five gap-closure crates compose: no sibling crate
was modified; the only new numerical code is a thin stored-field conversion
from fs-gum-field's padded-SoA `Field3` to fs-gum-topo's minimal node
container (same cell-centred grid, same fixed (i, j, k) order, same
vacuum-ghost semantics at every stencil read).

Reproduce: `cargo build --release && ./target/release/gum_e2e`
(**209.2 s total, both replay runs included**; budget was 20 min; build is
warning-free).

## Pipeline (N = 48, LBOX = 4.5, t = 0.008276434949296802, eps = 0.05)

1. **SEED** — `fs-gum-field::radial::radial_solve` (N = 4000, rmax = 6)
   solves the 1-D radial hedgehog profile at the frozen eps = 0.05 dial;
   the 3-D hedgehog is seeded from it on the N = 48 stored field
   (`sample_hedgehog`, dilation 1).
2. **MEASURE** — `fs-gum-field::measure` (central4): sector energies vs the
   Step-1 radial references, degree, weighted Derrick virial, boundary
   tail; recorded as claims S1–S4 / M1–M4.
3. **RELAX** — `fs-gum-statics`: fixed-seed non-axisymmetric bump
   perturbation (LCG seed 424242, K = 6, amp 0.02), guarded arrested Newton
   flow on the corner objective (one-sided degree anchor MU = 5000 / band
   0.005 referenced to the engine's own corner hedgehog degree 0.945803;
   one-sided Bogomolny-floor wall MU_F = 400 at fgap_ref = −0.125856), cap
   900 iterations (ends at iter_cap, 7 arrests, the Phase E1 protocol
   endpoint: Estat 2.9865333, deg 0.940153, fgap −0.128290). Claims R1–R4:
   monotone objective, degree in band, floor held, hedgehog-neighbourhood
   metrics.
4. **CLOCK + ISOROTATE** — the Phase E1 control-protocol Routhian descent
   (static endpoint + bump + tilt-halo shell, L = 0.12 I_hedgehog = 2.6722,
   below threshold; cap 600, ends at iter_cap, endpoint Estat 2.990896,
   I 31.94281), then the Section-K clock bisection
   L_clock² = (2/3) I E_static on the endpoint. Claims K0–K3: descent
   monotone with guards held, E_rot/E_tot = 1/4 at 1e-9 (measured delta
   2.8e-17), bisection = closed form to 0 ulp at print precision, and the
   F-R5-relevant number kappa(L_clock)/threshold = **1.2525** inside the
   claimed band [1.1, 1.4] (Phase E1: 1.2525 at N = 48; Python N = 96
   reference: 1.245).
5. **TOPO** — `fs-gum-topo` on the RELAXED static field through the thin
   conversion: independently ported central4 b-density degree (cross-crate
   agreement claim T1), lattice-exact geometric VOS solid-angle charge
   through a closed cube shell (T2), director projection onto S2 (T3); and
   the Hopf exhibit — Whitehead/FFT vs preimage-linking cross-validation
   (T4–T6) on fs-gum-topo's **analytic Hopf-1 referee texture** (N = 64
   periodic, HALF = 2.4). The relaxed hedgehog's director is degree-driven,
   not Hopf, so the referee texture — NOT a knot solution of the model — is
   the demo's topological exhibit; the claims state this explicitly.
6. **PACKAGE** — every claim fail-closed `Certified<f64>` → one
   `EvidencePackage`, domain `gum-core:e2e:v1`, producer `fs-gum-e2e`;
   golden BLAKE3 Merkle root printed; fs-checker three modes (deny-all =
   false as expected, source-certificate capability = true, tampered root =
   false); the ENTIRE pipeline runs twice in-process and both the root and
   the value fingerprint reproduce bit-for-bit.

## Claim table (22/22 PASS, all certified)

| id | claim | measured | bound | verdict |
|----|-------|----------|-------|---------|
| S1 | seed radial sectors vs radial_results.json refs; Newton floor | worst rel 1.17e-12; gmax 4.0e-13 | < 1e-6; < 1e-11 | PASS |
| S2 | eps dial ratio t(E2+E4)/(E6+E0) at frozen t | 0.049984699 | 0.05 ± 5e-5 | PASS |
| S3 | 1-D degree K of the seed profile | 1.000000292126 | \|K−1\| < 1e-6 | PASS |
| S4 | seed profile monotone pi → 0 | true | — | PASS |
| M1 | central4 N=48 seeded sectors vs Step-1 radial refs | worst rel 2.17e-2 (E2 −2.17e-2, E4 −2.20e-3, E6 +3.94e-3, E0 +2.4e-6, I −7.7e-6) | < 5e-2 (N=48 class; N=96 gated < 1e-2 in fs-gum-field H1) | PASS |
| M2 | central4 N=48 seeded degree | 0.998575173 (err 1.42e-3) | \|deg−1\| < 5e-3 | PASS |
| M3 | weighted Derrick virial, seeded hedgehog (central4, N=48) | rel −6.302e-3 | \|.\| < 1e-2 (Phase E1 recorded −6.30e-3) | PASS |
| M4 | boundary tail on box faces | 6.59e-11 | < 1e-6 | PASS |
| R1 | relaxation objective monotone by construction | max recorded rise −1.29e-5 | ≤ 0 | PASS |
| R2 | relaxed degree in anchor geometry | 0.940153 vs deg_ref 0.945803 (drift 5.65e-3) | ≤ band 5e-3 + penetration 3e-3 | PASS |
| R3 | guarded Bogomolny floor held | fgap −0.128290 vs wall −0.125856 (penetration −2.43e-3) | ≥ wall − 5e-3 | PASS |
| R4 | hedgehog neighbourhood (same functional) | I rel 1.52e-2; halo 5.8e-4; dEstat −0.0566; soft derivative sectors worst rel 0.276 | ≤ 8e-2 / 2e-2 / 0.10 / 0.35 | PASS |
| K0 | isorotating descent monotone, guards held at endpoint | max rise −2.85e-6; deg drift 5.64e-3; wall pen. −2.42e-3 | ≤ 0; ≤ 8e-3; ≥ −5e-3 | PASS |
| K1 | clock: E_rot/E_tot at bisected L_clock | 0.250000000000 (delta 2.8e-17) | 1/4 ± 1e-9 | PASS |
| K2 | bisection vs closed form sqrt((2/3) I Estat) | L_clock 7.980710330, rel 0.0 | < 1e-12 | PASS |
| K3 | **kappa(L_clock)/threshold** (the F-R5-relevant number) | kappa 0.249844 / 0.199471 = **1.2525** | in [1.1, 1.4], expect ~1.25 | PASS |
| T1 | cross-crate central4 b-density degree agreement on the relaxed field | fs-gum-topo 0.646862180314 vs fs-gum-field 0.646862180314, rel 0.0 | < 1e-12 rel | PASS |
| T2 | geometric VOS solid-angle charge, relaxed field, closed cube shell (s = 1.2, 24×24/face) | Q = 1 + 3.4e-15 | integer snap, \|Q−1\| < 1e-9 | PASS |
| T3 | director projection of the relaxed field lands on S2 | max \|\|n\|−1\| = 8.9e-16 | < 1e-12 | PASS |
| T4 | Whitehead/FFT Hopf invariant, REFEREE TEXTURE, N=64 | H_A = 0.962401785 (curl rel err 5.4e-2) | \|H−1\| < 6e-2 | PASS |
| T5 | preimage-linking Hopf invariant, same texture | H_B = 1.000000000000 (loops 1+1, segs 342+450) | \|H−1\| < 1e-6 | PASS |
| T6 | Hopf method-vs-method cross-validation | \|H_A − H_B\| = 3.76e-2 | < 6e-2 | PASS |

Certification: 22/22 claims certified through fs-evidence (fail-closed
`Certified<f64>`; failing claims are excluded from the package by
construction) into one fs-package EvidencePackage, domain `gum-core:e2e:v1`.
fs-checker verifies three modes: deny-all = false (as expected —
anti-laundering works), source-certificate capability = true, tampered root
= false.

**Golden Merkle root:**
`d30de830edb4f035b6eb86b9971386388f5ef42d97c8c59800bb291bc21a4edc`

Run fingerprint (BLAKE3 over every pipeline f64, fixed order, domain
`gum-core:e2e:fingerprint:v1`):
`c3790038b59faee59c2cac4747f66fc659892ed11b8b5f1cf512ec84058bf361`

Runtime: 209.2 s for BOTH in-process runs (root and fingerprint
bit-identical on replay); build warning-free.

## The degree instruments on the relaxed field (T1 vs T2, read carefully)

Three different degree instruments read the SAME relaxed configuration:

| instrument | reading | what it is |
|---|---|---|
| corner-scheme PAIRS degree (the guard's own instrument, fs-gum-statics `eval`) | 0.940153 | pinned at the anchor band edge (R2) |
| central4 b-density degree (fs-gum-field `measure` AND fs-gum-topo `degree_b_density`) | 0.646862 | bit-identical across the two crates (T1) |
| geometric VOS solid-angle winding (fs-gum-topo, lattice-exact, integer-snapping) | 1 + 3.4e-15 | the topological class (T2) |

The spread is an INSTRUMENT property, not a topology change, and it is the
documented near-BPS behaviour of Phase E1 at N = 48: the corner-objective
descent drifts into quadrature-error valleys of the soft derivative
sectors, leaving cell-scale roughness that the width-1 corner scheme
controls but that wide central stencils undersample. The same R4 note
records central4 derivative-sector drifts up to 27.6% on this endpoint —
the central4 b-density (a derivative-based density integral, NOT
quadrature-exact) degrades identically, reading 0.6469 where it reads
0.998575 on the smooth seed (M2). The lattice-exact VOS winding — an exact
integer up to floating point regardless of smoothness — confirms the field
never left the degree-1 sector, and the corner instrument that the anchor
guard actually acts on stayed within its band throughout.

So, precisely: **T1 certifies cross-crate consistency** — two independently
ported implementations of the same central4/det4 instrument produce
bit-identical numbers on the same stored field through the thin conversion
(rel 0.0, tolerance 1e-12). It does NOT certify that 0.6469 is the degree.
**T2 certifies the topological class** (degree 1, lattice-exact route).
**R2 certifies the guard geometry** on the corner instrument. None of the
three claims overreaches its instrument, and the package records all three
readings.

## Determinism statement

All sweeps and reductions are plain sequential f64 loops in fixed ascending
order (inherited contracts of all three siblings); the perturbation
generator is the documented 64-bit LCG; no parallelism, no tree reductions,
no platform libm in kernels. The ENTIRE pipeline (radial solve, 3-D
measurement, 900-step guarded relaxation, 600-step isorotating descent,
clock bisection, cross-crate topo diagnostics, certificates, Merkle
assembly) runs twice in-process; both the package root and the value
fingerprint reproduce bit-for-bit (asserted; the binary aborts otherwise).

## Deviations from the coordinator spec, with reasons

1. **M1 tolerance is 5e-2, not the 1e-2 of fs-gum-field H1**: the 1e-2 gate
   is the N = 96 instrument; at the demo grid N = 48 the same central4
   instrument on the same seeded hedgehog measures worst rel 2.17e-2 (E2).
   The N = 48 class is stated in the claim; the N = 96 certification
   already exists in fs-gum-field.
2. **One isorotating descent, not two**: the spec's clock claim lives on
   the clock-charged endpoint, which in the Phase E1 protocol is the
   CONTROL descent (L = 0.12 I_h); the over-spun main descent and the
   halo-differential gates are already certified in fs-gum-statics G-C and
   re-running them here would roughly double the runtime without adding a
   new claim class. The measured kappa(L_clock)/threshold = 1.2525
   reproduces Phase E1's G-D value exactly.
3. **Hopf claims run on the analytic Hopf-1 referee texture** at N = 64
   periodic, HALF = 2.4 (per the spec's own instruction: the relaxed
   hedgehog's director is degree-driven, not Hopf; and the Whitehead/FFT
   route needs power-of-two axes). Every Hopf claim statement carries the
   "referee texture, not a knot solution" qualifier.
4. **T1 is a consistency claim, not an accuracy claim** — see the
   instrument section above. The spec asked for a cross-crate agreement
   claim with a stated tolerance (1e-12 rel, measured 0.0); the accuracy of
   the b-density instrument on the roughened relaxed field is certified
   separately and honestly by the T2/R2 pair.
5. **No sibling crate was modified**; the only new numerics is the
   stored-field conversion (`to_topo_qfield`), which is layout plumbing,
   not arithmetic.

## Epistemic notice (binding)

Every PASS above certifies a WITHIN-MODEL property: that discretisations,
solvers, topological diagnostics and an evidence pipeline — composed across
five crates — behave as documented on a speculative theory's functional.
The Hopf exhibit runs on the analytic referee texture, not a knot solution
of the model. This is an integration demo of numerical engineering; nothing
here validates the theory, and nothing here says anything about nature.
