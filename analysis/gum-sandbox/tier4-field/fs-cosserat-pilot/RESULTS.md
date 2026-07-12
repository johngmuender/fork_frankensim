# Tier-4B — certified 3-D texture-field diagnostics on the analytic BPS compacton (the "Cosserat pilot")

Reproduce with `cargo run --release` from this directory (standalone crate,
outside the workspace, like `tier0-gauntlet`; runtime ~0.2 s including the
in-process replay). Everything is evaluated on ANALYTIC configurations — no
solver, no RNG: the exact eps = 0 BPS compacton hedgehog and its
volume-preserving family-A pullback. Purpose: demonstrate that frankensim's
evidence machinery (`fs-ivl` -> `fs-evidence` -> `fs-package` -> `fs-checker`)
can certify 3-D texture-field diagnostics with a bit-deterministic golden
Merkle root.

## Conventions (stated; all match tier2-closure/axi_solve.py)

- SU(2) field `q = (q0, q1, q2, q3)`, `|q| = 1`; hedgehog
  `q = (cos f(r), sin f(r) rhat)`.
- Compacton profile `f0(r) = 2 arccos(r/R*)` for `r <= R*`, vacuum outside;
  **R\* = 2^(5/6)** (axi_solve.py's `RSTAR`), the (6+0)-sector radius at
  which E6 = E0 = half the BPS bound. Trig-free closed form used on the grid:
  `q0 = 2 r^2/R*^2 - 1`, `q_i = (2/R*^2) sqrt(R*^2 - r^2) x_i`.
- Closed forms derived here and verified by fs-ivl quadrature (u = r/R*):
  `E0 = 2 R*^3 INT u^2(1-u^2) du = 16 sqrt2/15`,
  `E6 = (64/R*^3) INT u^2(1-u^2) du = 16 sqrt2/15`,
  `I = (64 pi R*^3/3) INT u^4(1-u^2) du = 128 pi R*^3/105` —
  matching axi_solve.py's `E6_COMPACTON = E0 = 16 sqrt2/15`,
  `I_COMPACTON = 128 pi R*^3/105` exactly.
- **D1 method choice**: the central-difference b-density (not the
  face-triangulated solid-angle kernel):
  `b = -(1/2 pi^2) det(q, dq/dx, dq/dy, dq/dz)`, i.e. the epsilon-contraction
  `-(1/12 pi^2) eps_{ijk} eps^{abcd} q_a d_i q_b d_j q_c d_k q_d`, sign fixed
  so the `f: pi -> 0` hedgehog has degree +1. Derivatives are 2nd-order
  central differences of the analytic field (stencil points evaluated in
  closed form — no arrays, no boundary special-casing).
- Grid: cell-centred N^3 over `[-2.4, 2.4]^3` (covers the deformed support
  `e R* = 2.1125` with stencil margin), N = 48 (h = 0.1) and 96 (h = 0.05),
  fixed (i,j,k) summation order, plain f64 — fully deterministic.
- SDiff family-A map (axi_RESULTS.md, family A, d = 1):
  `q_lam(x,y,z) = q_base(x/e, y/e, z/p)`, `e = lam^(-1/3)`, `p = lam^(2/3)`,
  unit volume factor `e*e*p = 1`; spot check at lam = 0.6.
- 1-D certified references: fs-ivl outward-rounded interval Riemann
  enclosures (200 000 cells) of the reduced smooth integrands, intersected
  with the closed-form interval (two routes, one enclosure), plus the exact
  rational identities `INT u^2(1-u^2) = 2/15`, `INT u^4(1-u^2) = 2/35`.
- Certificate canonical string `id|statement|lo_bits|hi_bits`, BLAKE3 domain
  **`gum-tier4:cert:v1`** (policy domain `gum-tier4:policy:v1`), producer
  `fs-cosserat-pilot`. No HashMap anywhere near the hashes; claims enter the
  package in fixed Vec order.

## Claim table (all `Certified<f64>`, fail-closed; 13/13 PASS)

| id | claim | measured value / enclosure | bound | verdict |
|---|---|---|---|---|
| D1a | grid degree B, N = 48 | 0.9830385118051 (err 1.696e-2) | \|B-1\| < 3.0e-2 | **PASS** |
| D1b | grid degree B, N = 96 | 0.9956367241585 (err 4.363e-3) | \|B-1\| < 8.0e-3 | **PASS** |
| D1c | degree-error N-scaling | observed order p = 1.96 | err decreasing, p > 1 | **PASS** |
| D2a | E0 1-D certified enclosure vs 16 sqrt2/15 | 1.508494466531, width 1.6e-15 | routes intersect, width < 1e-4 | **PASS** |
| D2b | grid E0, N = 96, two-way check | 1.508472719913 (drift -2.175e-5; N=48 -1.258e-4) | within enclosure +- 2.0e-3 | **PASS** |
| D3a | E6 1-D certified enclosure vs 128/(15 R*^3) | 1.508494466531, width 1.6e-15 (= E0: BPS) | routes intersect, width < 1e-4 | **PASS** |
| D3b | grid E6, N = 96, two-way check | 1.502173763242 (drift -6.321e-3; N=48 -2.770e-2) | within enclosure +- 2.0e-2 | **PASS** |
| D4a | I 1-D certified enclosure vs 128 pi R*^3/105 | 21.66434346988, width 3.6e-14 | routes intersect, width < 1e-3 | **PASS** |
| D4b | grid I, N = 96, two-way check | 21.66361366437 (drift -7.298e-4; N=48 -4.145e-3) | within enclosure +- 2.0e-2 | **PASS** |
| D5a | SDiff drift \|E0(0.6) - E0(1)\|, N = 96 | 3.271e-5 (N=48: 1.213e-4) | < 2.0e-3 | **PASS** |
| D5b | SDiff drift \|E6(0.6) - E6(1)\|, N = 96 | 1.129e-3 (N=48: 4.455e-3) | < 2.0e-2 | **PASS** |
| D5c | SDiff drift \|I(0.6) - I(1)\|, N = 96 | 8.915e-4 (N=48: 8.018e-3) | < 2.0e-2 | **PASS** |
| D6 | E_rot/E_tot = 1/4 with clock L^2 = (2/3) I E_static | interval [0.25 +- 1.2e-15], width 2.4e-15 | contains 0.25, width < 1e-4, rational identity 3*4 = 12 | **PASS** |

Tolerances were set from the O(h^2) midpoint-rule error model with the
compacton's boundary sqrt-cusp, calibrated once on the N = 48/96 pair; every
measured drift sits below its bound (margins 1.8x–92x) and the raw drifts are
printed in each claim's note so the reader certifies against the numbers, not
the bounds alone.

## N-scaling of the D1 degree error

| N | h | B | \|B - 1\| |
|---|---|---|---|
| 48 | 0.1 | 0.983039 | 1.696e-2 |
| 96 | 0.05 | 0.995637 | 4.363e-3 |

Observed convergence order **p = 1.96** (clean O(h^2) despite the profile's
sqrt-cusp at r = R*, because the b-density itself vanishes like
sqrt(1 - u^2) there). Every other grid diagnostic also contracts by ~4x from
N = 48 to 96 (see the table's parenthetical N=48 drifts); the D5 invariance
drifts are the discretisation error made visible, exactly as the F-R4 lemma
predicts (analytic invariance is exact: E0 and I carry no gradients, and
INT b^2 is invariant because b transforms with unit Jacobian under a
volume-preserving pullback).

## Golden root and checker modes

**Golden Merkle root (BLAKE3, fs-package v6):**

```
f88731afb676d76570ac719fc5b79ced053a1fa95e17750510e375debc64fd7a
```

| mode | capabilities | expected | observed |
|---|---|---|---|
| 1 | deny-all | FAIL (unauthenticated Verified claims refused — anti-laundering) | passed = false ✓ |
| 2 | source-certificate capability (recompute-and-compare against `gum-tier4:cert:v1` content hash) | PASS | passed = true ✓ |
| 3 | tampered root (first byte xor 0xff), same capability | FAIL | passed = false ✓ |

## Replay determinism

- **In-process**: the entire pipeline (grid sums, interval quadratures,
  certificates, Merkle assembly) runs twice inside the binary; the second
  root is asserted bit-identical to the first (`replay determinism: second
  full run root == first`).
- **External**: the binary was executed twice and the full stdout diffed —
  byte-identical, same golden root both times.

## Epistemic notice

Every PASS certifies a WITHIN-MODEL property: that frankensim's evidence
machinery can wrap 3-D texture-field diagnostics of a speculative theory's
analytic configuration in fail-closed certificates with a reproducible golden
root, and that the discretisation converges to the certified 1-D enclosures.
This is certification of diagnostics, not physics validation — nothing here
bears on nature, and nothing here adjudicates the corpus's physical claims
(see tier2-closure adjudications for those).
