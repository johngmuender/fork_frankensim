# fs-gum-topo — GUM core Phase B2: topological-charge diagnostics

Spec: `analysis/gum-sandbox/gum-core/gap_survey_v3.json`, dimension
"Topological-charge diagnostics". Standalone crate (empty `[workspace]`,
path deps into `crates/`), library + gates binary
(`cargo run --release --bin gates`). 25/25 referee checks PASS,
bit-identical two-run replay, total gates wall time ~5 s
(budget 10 min).

## What is implemented

**Field container** (`src/field.rs`). Minimal self-contained
`Field3<const C>`: flat row-major `Vec<[f64; C]>` on an n^3 node grid,
fixed (i, j, k) order everywhere (determinism feeds the certificate
hashes), `Bc::Periodic` or `Bc::Vacuum(ghost)` (field3d_solve.py's
padding), trilinear sampling. `QField = Field3<4>`,
`DirectorField = Field3<3>`. A sibling crate (fs-gum-field) is building
the production stored-field type concurrently; this one is deliberately
minimal so a later integration pass can unify them without touching the
algorithms.

**Hopf map / director projection** (`hopf_project`):
n = R(q) z-hat = (2(q1 q3 + q0 q2), 2(q2 q3 - q0 q1), 1 - 2(q1^2 + q2^2)),
node-wise, vacuum ghosts projected consistently ((1,0,0,0) -> z-hat).

**S3 degree on stored quaternion fields** (`src/degree.rs`).
b = -(1/2 pi^2) det(q, dq/dx, dq/dy, dq/dz), sign fixed so the
f: pi -> 0 hedgehog has degree +1. Three schemes: `Central2` (the
fs-cosserat-pilot density, ported operation-for-operation), `Central4`
and `Corner` (midpoint at (N+1)^3 corners) ported from
field3d_solve.py. Returns (deg, b-field, E6 = pi^3 INT b^2). Plus the
geometric route: `hedgehog_charge_geometric` sums van
Oosterom-Strackee signed spherical-triangle solid angles of the
director over a closed triangulated probe surface (`cube_shell`) —
lattice-exact, integer-snapping.

**Hopf invariant, Method A — Whitehead/FFT** (`src/hopf.rs`).
F_ij = n.(d_i n x d_j n) by central differences; B = (F_yz, F_zx, F_xy)
(hedgehog flux normalization 4 pi); net-flux obstruction check BEFORE
the solve (see deviations); Coulomb-gauge A-hat = i (kappa x B-hat)/|kappa|^2
with the DISCRETE wavevector kappa_m = sin(2 pi m/N)/h, pinned to
exactly 0 at m = 0 and m = N/2 (the f64 sin(pi) ~ 1.2e-16 residue
otherwise turns the zero modes into a catastrophic amplification —
found and fixed against the max-imaginary monitor); A-hat(0) = 0 makes
the answer canonical given zero net flux;
H = (1/(4 pi)^2) INT A.B, fixed loop order. `SIGN_WHITEHEAD = +1`
fixed once against the referee.

**Hopf invariant, Method B — preimage linking** (`src/linking.rs`).
Two regular values y1 (z-hat tilted by fixed irrational-looking angles,
deterministic re-tilt-and-retry schedule for degeneracies) and
y2 = -y1; marching tetrahedra over the 6 Freudenthal/Kuhn tets per cube
(face diagonals consistent across cubes); per tet the frame coordinates
u = n.e1, v = n.e2 are linear, so the preimage segment is unique;
hemisphere check w = n.y > 0 discards the antipodal sheet; segments
oriented by grad u x grad v and chained by EXACT face keys (sorted
global node triples — no floating-point matching); Gauss linking by the
exact polyline segment-pair solid-angle sum (two VOS triangles per
pair). Output is an integer up to floating point. `SIGN_LINK = +1`
fixed once against the same referee.

**Disclination diagnostics** (`src/disclination.rs`). `z2_holonomy`:
sign-gauge transport; the lift closes to sigma n_0 with sigma the
product of the signs of consecutive RAW dots (the running-lift parity
telescopes to exactly this — an earlier flip-event-parity draft was
caught by the winding/holonomy consistency referee G6c and fixed);
CoreHit tolerance 0.1. `planar_winding`: atan2 increments of the lifted
director projected on the loop plane, half-integer iff the Z/2 class is
nontrivial; fails typed on director parallel to the normal.
`disclination_lines`: Z/2 holonomy of every elementary plaquette;
pierced plaquettes' dual edges chained through cubes; structural checks:
even pierced-face count per cube (lines cannot end — the dd = 0
analogue), junction detection, closed-or-boundary-terminated assertion.

**Certification** (`src/certify.rs`). The pilot's fail-closed pattern as
reusable code, domain `gum-core:topo:v1`: Check -> Certified<f64> ->
EvidencePackage claims -> BLAKE3 Merkle root; fs-checker deny-all
(expect fail), source-certificate capability (expect pass), tampered
root (expect fail); the gates binary runs the whole pipeline twice and
asserts a bit-identical root.

## Gate table (25/25 PASS)

| id | claim | measured | bound | verdict |
|---|---|---|---|---|
| G1a | Central2 degree on stored field, N=48, reproduces pilot golden 0.9830385118051 | delta 3.286e-14 | < 1e-12 | PASS |
| G1b | Central2 degree, N=96, golden 0.9956367241585 | delta 2.909e-14 | < 1e-12 | PASS |
| G1c | Central2 error N-scaling | p = 1.96 (1.696e-2 -> 4.363e-3) | decreasing, p > 1 | PASS |
| G1d | Central4 degree, N=96 | \|deg-1\| = 7.914e-5 (N=48: 9.179e-4) | < 5e-3 | PASS |
| G1e | Corner degree, N=96, O(h^2) scaling | \|deg-1\| = 5.092e-3, p = 1.93 | < 2.5e-2, p > 1 | PASS |
| G1f | E6 from b-field vs 16 sqrt2/15 | 1.502173763242 (drift -6.321e-3; = pilot D3b digits) | < 2e-2 | PASS |
| G1g | geometric VOS hedgehog charge, cube shell | Q = 1 + 6.2e-15 | integer snap, < 1e-9 | PASS |
| G2a | hopf_project lands on S2 | max \|\|n\|-1\| = 1.2e-15 | < 1e-12 | PASS |
| G3a | Whitehead H, N=64 | H = 0.962401785 | \|H-1\| < 6e-2 | PASS |
| G3b | Whitehead H, N=128 | H = 0.988917007 | \|H-1\| < 2e-2 | PASS |
| G3c | Whitehead error N-scaling | p = 1.76 (3.760e-2 -> 1.108e-2) | decreasing, p > 1 | PASS |
| G3d | curl_h A vs B (rel L2) | 4.352e-2 (the non-solenoidal O(h^2) part of discrete B) | < 5e-2 | PASS |
| G3e | net-flux obstruction referee (baby-skyrmion slab) | typed NetFlux error | rejected, never a number | PASS |
| G4a | Gauss-linking kernel calibration | Lk(Hopf link) = -1.000000000000, Lk(unlinked) = -4e-18 | \|Lk\| = 1 / 0 | PASS |
| G4b | preimage-linking H, N=64 | 1.000000000000 (loops 1+1, segs 342+450) | \|H-1\| < 1e-6 | PASS |
| G4c | preimage-linking H, N=96 (non-pow2, FFT-free) | 1.000000000000 (segs 518+680) | \|H-1\| < 1e-6 | PASS |
| G4d | preimage-linking H, N=128 | 1.000000000000 (segs 686+908) | \|H-1\| < 1e-6 | PASS |
| G4e | method-vs-method cross-validation, N=128 | \|H_A - H_B\| = 1.108e-2 | < 2e-2 | PASS |
| G5a | uniform control | H_A = 0 exactly, H_B = 0, 0 loops | = 0 | PASS |
| G5b | global O(3) target rotation | \|H(rot) - H\| = 3.3e-16 | < 1e-9 | PASS |
| G5c | parity mirror | H_A = -0.962401785 (drift 3.3e-16), H_B = -1.000000000 | = -H / -1 | PASS |
| G6a | Z/2 holonomy, +1/2 line | NonTrivial at r~0.6 and r~1.8; Trivial off-axis and uniform | exact classes | PASS |
| G6b | planar winding | +0.500000000000 (both radii), -0.5 for -1/2, 0 off-axis | half-integers exact | PASS |
| G6c | winding-holonomy consistency | half-integer iff NonTrivial | exact | PASS |
| G6d | defect-line tracing | one straight boundary-terminated line (48 dual points, xy straightness 0.0), -1/2 ditto, control empty; even-count check held | structural | PASS |

Certification: 25/25 claims certified; Merkle root
`9cfbe9a6c42a41e1a891413c6763fca1f0d5be656a174c90a25e32f1744b89c1`;
checker modes deny-all = false / capability = true / tampered = false
(all as expected); replay root bit-identical.

## Method-vs-method Hopf agreement

| N | Whitehead H_A | preimage-linking H_B | \|H_A - H_B\| |
|---|---|---|---|
| 64 | 0.962401785 | 1.000000000000 | 3.76e-2 |
| 96 | (FFT route needs 2^m axes) | 1.000000000000 | — |
| 128 | 0.988917007 | 1.000000000000 | 1.11e-2 |

H_A converges to the integer at observed order 1.76 (slightly below 2:
the compacton's sqrt-cusp at r = R* enters n = R(q) z-hat quadratically,
unlike the b-density where it self-cancels); H_B returns the exact
integer at every resolution, as the survey predicted for the
preimage/linking construction. Mirrored texture: both methods give -1;
uniform control: both give exactly 0.

## Deviations from the survey spec, with reasons

1. **fs-fft not used.** `crates/fs-fft` pulls `fs-exec ->
   asupersync`, an external sibling repository not present in this
   tree, so the crate cannot build fs-fft standalone. The build spec's
   stated fallback applies: `src/fft.rs`, a minimal internal radix-2
   FFT (deterministic fixed butterfly/axis order, power-of-two only,
   forward unnormalized / inverse 1/n). The Whitehead API keeps the
   structured `NotPow2` rejection so a later swap to fs-fft is drop-in.
2. **Obstruction check restated in lattice-exact form.** The survey's
   mean-B test (tol 1e-10 x max|B|) is unachievable on the
   central-difference B: its box mean carries an O(h^{3/2})
   discretization residue (measured 3.1e-2 at N=32, 1.2e-2 at N=64 on
   the sqrt-cusp referee — z-component only; x/y vanish by symmetry at
   1e-17). The net-flux obstruction is instead measured by the
   Berg-Luscher solid-angle winding of the three coordinate-slice maps,
   an exact integer up to floating point (0 admissible, else typed
   `NetFlux`), which separates topology from noise at machine
   precision; the raw mean is still computed and reported. The
   baby-skyrmion slab referee (net B_z flux) is rejected as required.
3. **Whitehead referee resolutions 64/128 instead of 48/96.** The FFT
   route needs 2^m axes; the survey itself notes this constraint. The
   degree diagnostics keep the pilot's 48/96 pair (and reproduce its
   certified digits to ~3e-14, the ulp-level cost of storing the field
   before differencing); the FFT-free preimage route is additionally
   gated at N=96 to show it has no power-of-two restriction.
4. **Whitehead tolerances calibrated on the measured pair** (6e-2 at
   N=64, 2e-2 at N=128, cross-validation 2e-2), the pilot's procedure:
   bounds set once from the observed error model, raw drifts printed in
   every note so the reader certifies against numbers, not bounds.
   Observed order 1.76 rather than the clean 1.96 of the b-density —
   the cusp does not self-cancel in the director field.
5. **fs-ivl enclosures are used as claim carriers** (point intervals
   plus the pilot's fail-closed Certified path); no new 1-D quadrature
   references were needed here since every referee target is an exact
   integer/half-integer or a pilot-certified golden.

## Known limitations

- **Grid resolution near defect cores.** All loop diagnostics are
  CoreHit-guarded (default tol 0.1 on transported dots): a loop passing
  within ~one cell of a disclination core returns a typed error rather
  than a wrong class; reroute or refine. Similarly the preimage
  extractor treats edge-grazing crossings as degeneracies and re-tilts;
  a field with |grad n| ~ 1/h near the requested preimage can exhaust
  the retry schedule (observed benignly at N=32 on the hopfion).
- The Whitehead route is defined on periodic power-of-two grids with
  zero net flux; compact-support textures on other grids must be
  resampled/padded (or use the preimage route, which is
  grid-agnostic but needs curves to close inside the box —
  torus-wrapping preimages surface as `OpenCurve`, the same obstruction
  class the FFT route rejects as `NetFlux`).
- `disclination_lines` resolves cubes with exactly 2 pierced faces;
  >= 4 (line junctions at the grid scale) return a typed `Degenerate`
  rather than guessing a pairing.
- The corner scheme remains the robust-not-accurate quadrature
  (field3d_RESULTS.md section G): |deg-1| = 5.1e-3 at N=96 on the
  compacton, ~64x the central4 error, with clean O(h^2) contraction.

## Epistemic notice (binding)

Every PASS above certifies a WITHIN-MODEL property of diagnostic
machinery evaluated on a speculative theory's analytic configurations.
It validates the discretizations, the topological algorithms, and the
evidence pipeline — never the physics.

---

# Phase F3 update (2026-07-12): fs-fft is now importable — swap deferred

Deviation 1 above ("fs-fft not used") is no longer forced. Phase F3 made
`crates/fs-fft` (and `crates/fs-la`) buildable from out-of-workspace
crates in this checkout: fs-fft's `fs-exec` dependency is now optional
behind a default-ON `exec` feature, so

```toml
fs-fft = { path = "../../../../crates/fs-fft", default-features = false }
```

resolves and builds without the `../../../asupersync` sibling repository
(the full serial `Fft`/`FftNd`/`RealFft` surface remains; only the
TilePool-tiled `*_pooled` API is compiled out). See
`analysis/gum-sandbox/gum-core/PORTABILITY.md`.

The swap is deliberately NOT performed in F3: `src/fft.rs` is
golden-gated, and fs-fft's Stockham kernels are not guaranteed to
reproduce this file's radix-2 butterfly order bit-for-bit, so a swap
would risk a golden re-freeze for zero functional gain. When a future
phase performs it, keep the structured `NotPow2` rejection at the
Whitehead API boundary (fs-fft's `Fft::new` panics on non-power-of-two
instead of returning a typed error).
