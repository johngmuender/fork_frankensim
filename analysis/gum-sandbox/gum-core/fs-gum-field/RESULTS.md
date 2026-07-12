# fs-gum-field — GUM core Phase B1 results

Crate: `analysis/gum-sandbox/gum-core/fs-gum-field/` (standalone `[workspace]`
opt-out, path deps on `crates/fs-*`, std-only dependency closure — the
fs-cosserat-pilot pattern).

Scope: the stored SU(2)/quaternion field type and its **measurement** kernels
— gap-survey-v3 dimension "SU(2) field type" plus the measurement half of
dimension "Skyrme/chiral". Port sources: `tier4-field/field3d_solve.py`
(stored-field stencils, layout, sector energies, gates) and
`tier2-closure/radial_solve.py` (1-D solver). Certification template:
`tier4-field/fs-cosserat-pilot`. Descent machinery (gradients, arrested
Newton flow, penalty guards) is deliberately out of scope; the padded-ghost
SoA layout is kept identical to the Python engine so the adjoint stencils
port 1:1 in the next phase.

Reproduce: `cargo build --release && ./target/release/gum_field_gates`
(~2 s total, both replay runs included; build is warning-free).

## Conventions

Tier-2 unit map (`axi_RESULTS.md`, a6 = a0 = m = 1):

```
q = (q0,q1,q2,q3), |q| = 1, vacuum (1,0,0,0); ghost rind = 2 layers of vacuum
E2 = (1/4pi) INT sum_i |D_i q|^2
E4 = (1/4pi) INT sum_{i<j} (|D_i q|^2 |D_j q|^2 - (D_i q . D_j q)^2)
E6 = pi^3 INT b^2,  b = -(1/2 pi^2) det[q, D_x q, D_y q, D_z q]
E0 = (1/4pi) INT (1 - q0);   I = INT 2 (q1^2 + q2^2);   degree = INT b
E_static = t (E2 + E4) + E6 + E0,   t = 0.008276434949296802 (eps = 0.05)
weighted Derrick virial: t E2 - t E4 - 3 E6 + 3 E0 = 0
```

The b-density uses the first-row-expansion 4x4 determinant of
fs-cosserat-pilot (`sector::det4`), analytically identical to field3d's
6-term PAIRS expansion; the sign is fixed so the f: pi -> 0 hedgehog has
degree +1. Grid: cell-centred N^3, x_i = -half + (i+0.5) h, h = 2 half/N;
HALF = 2.4 for the compacton, LBOX = 4.5 for the eps = 0.05 hedgehog.
Schemes: `Central2` (pilot stencil), `Central4` (c1 = 8/(12h), c2 = 1/(12h);
the field3d measurement engine), `Corner` (compact midpoint scheme at the
(N+1)^3 corners; the field3d descent engine, O(h^2)). E0 and I are
scheme-independent quadrature-exact cell sums.

## Gate table (25/25 PASS)

| id | gate | measured | tolerance | verdict |
|----|------|----------|-----------|---------|
| A1 | compacton degree N=48 (stored field, order-2) | 0.9830385118051, err 1.696e-2; vs pilot rel -3.4e-14 | < 3.0e-2 and pilot match < 1e-9 | PASS |
| A2 | compacton degree N=96 | 0.9956367241585, err 4.363e-3; vs pilot rel +3.0e-14 | < 8.0e-3 and pilot match < 1e-9 | PASS |
| A3 | degree-error N-scaling | measured order p = 1.96 | decreasing, p > 1 | PASS |
| A4 | compacton E0 N=96 vs enclosure 16 sqrt2/15 | 1.508472719913 (drift -2.175e-5); vs pilot rel +6.4e-14 | enclosure +- 2e-3, pilot match | PASS |
| A5 | compacton E6 N=96 vs enclosure (= E0 at BPS) | 1.502173763242 (drift -6.321e-3); vs pilot rel -1.5e-13 | enclosure +- 2e-2, pilot match | PASS |
| A6 | compacton I N=96 vs enclosure 128 pi R*^3/105 | 21.66361366437 (drift -7.298e-4); vs pilot rel -5.1e-14 | enclosure +- 2e-2, pilot match | PASS |
| R1 | radial sectors N=4000 vs radial_results.json | worst rel 1.17e-12; gmax 4.0e-13 | < 1e-6, gmax < 1e-11 | PASS |
| R2 | eps ratio at frozen t | 0.049984699 | 0.05 +- 5e-5 | PASS |
| R3 | 1-D weighted Derrick virial (rel) | +2.385e-6 | < 1e-5 | PASS |
| R4 | Yukawa tail mu_fit vs 1/sqrt(2t) | 7.7726254884 vs 7.7725468626 (rel +1.01e-5); rms 3.0e-6 | < 1e-4 (and recorded mu_fit < 1e-5) | PASS |
| R5 | BPS deficit of E6+E0 over 32 sqrt2/15 | +0.004817757679649 | > 0, matches recorded < 1e-6 | PASS |
| R6 | 1-D degree K | 1.000000292126 | |K-1| < 1e-6, recorded match < 1e-9 | PASS |
| R7 | profile monotone pi -> 0 | true | — | PASS |
| R8 | eps-dial outer loop (N=1000) | t = 0.008276339459 in 3 outers, ratio 0.04998470 | ratio +- 5e-5, t within 1% of frozen | PASS |
| H1 | central4 N=96 sectors vs radial refs | rels: E2 -2.69e-3, E4 -2.71e-4, E6 +4.06e-4, E0 +9.9e-9, I -5.1e-8 | worst < 1e-2 | PASS |
| H2 | central4 N=96 degree | 0.999911418 (err 8.9e-5) | < 5e-3, recorded match < 1e-5 | PASS |
| H3 | central4 N64/N96 convergence ratios | E2 3.62, E4 3.41, E6 4.06; E0/I 5.1e-8 (N96), 8.4e-6 (N64) | ratios in [2.5, 6.0]; E0/I < 1e-7 / 1e-5 | PASS |
| H4 | 3-D weighted Derrick virial (rel), N=96 | -6.615760e-4 | < 2e-3, recorded match < 1e-6 | PASS |
| H5 | boundary tail on box faces | 5.500e-11 | < 1e-6 | PASS |
| H6 | field3d referee rel errors reproduced | worst delta 3.6e-13 | < 1e-5 | PASS |
| C1 | corner N=96 O(h^2) fingerprint | worst rel 2.21e-2, referee delta 2.6e-13 | < 5e-2, referee match < 1e-5 | PASS |
| C2 | corner N64/N96 ratios | E2 2.11, E4 2.18, E6 2.22 | in [1.7, 2.8] | PASS |
| C3 | corner N=96 degree | 0.986011569 | matches recorded 0.986012 < 1e-5 | PASS |
| F1 | renormalize drift + idempotence | first pass 2.16e-2, second 3.3e-16 | > 1e-4, then < 1e-14 | PASS |
| F2 | Maurer-Cartan identity \|a_i\|^2 = \|D_i q\|^2 | worst rel 8.2e-16; max \|Re a_i\| 2.2e-1 (FD tangency diag) | < 1e-12 | PASS |

Certification: 25/25 claims certified through fs-evidence (fail-closed
`Certified<f64>`) into one fs-package EvidencePackage; fs-checker verifies
three modes (deny-all = false as expected, source-certificate capability =
true, tampered root = false), domain `gum-core:field:v1`.

Golden Merkle root:
`6c1e785071a753e40a33c4d8532886d2ca20bd245d7f184f67b5eed9c6c9fb1f`

Run fingerprint (BLAKE3 over every gate f64, fixed order, domain
`gum-core:field:fingerprint:v1`):
`7c9931a50d7539dec295b617720fc09c81371c3dd0fa19bfe82cfabf22f97cfe`

## Determinism statement

All sweeps and reductions are plain sequential f64 loops in fixed ascending
(i, j, k) order — no tree reductions, no parallelism, no allocation-order
dependence feeding any number. Elementary functions in kernels and seeds
come from `fs_math::det` (strict-mode, declared ULP budgets); the only
intrinsics are IEEE `sqrt` and the `cbrt`-for-R* precedent inherited from
fs-cosserat-pilot. The entire pipeline (radial solves, 3-D sweeps, fits,
certificates, Merkle assembly) runs twice in-process; both the package root
and the value fingerprint reproduce bit-for-bit (asserted, gate D).

Notes against the pilot: the pilot evaluates its analytic closure at stencil
points, this crate reads stored neighbours through the vacuum ghost rind;
for the compacton on HALF = 2.4 both paths see identical field values, and
the A-gates reproduce the pilot's certified decimals to < 2e-13 rel (the
residual is the pilot's own 13-digit print rounding plus ULP-level
`x + h` vs `-half + (i+1.5)h` argument differences — not a scheme
difference).

## Epistemic notice (binding)

Every PASS above certifies a WITHIN-MODEL property: that a discretisation,
a port, and an evidence pipeline behave as documented on a speculative
theory's functional. Nothing here validates the theory, and nothing here
says anything about nature.
