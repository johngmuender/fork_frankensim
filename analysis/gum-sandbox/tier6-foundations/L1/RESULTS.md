# L1 — Contextuality exhibit: the flipped Stern-Gerlach (F-T6-2a-EXEC)

**Goal.** First in-repo execution of the de Broglie-Bohm guidance-sector
demonstration (Bell 1982 / Bricmont 2019) that "measurements don't measure":
with the same initial wave function and the same initial particle position, a
reversed Stern-Gerlach gradient produces the *identical* spatial trajectory
while the outcome LABEL (spin up vs spin down) flips.

Within-model; nothing here bears on nature.

## Gates

| id | spec | measured | verdict |
|----|------|----------|---------|
| G1 | median (z=0) crossings across all 4000 trajectories x 480 steps = 0 | 0 crossings (original run; flipped run also 0) | **PASS** |
| G2 | every trajectory: exit side unchanged under k -> -k, spin label flipped; max \|Z_flip(T)-Z_orig(T)\| at machine/integration precision | all 4000 same side; all 4000 labels flipped (upper side reads "up" at k=+2.5, "down" at k=-2.5); max \|dZ(T)\| = 0.0 exactly; max \|rho_flip - rho_orig\| = 0.0 exactly | **PASS** |
| G3 | upper-exit fraction = 0.5 within 3 sigma binomial, N = 4000 | 1971/4000 = 0.49275, deviation 0.92 sigma (sigma_binom = 0.0079) | **PASS** |
| G4 | dt -> dt/2, re-integrate 100 trajectories: max deviation <= 1e-6 x separation | 5.85e-9 vs threshold 1.20e-5 (margin ~2000x) | **PASS** |

## Key numbers

- Grid: N = 8192 points on periodic z in [-40, 40]; sigma = 1, k = 2.5, dt = 0.005, T = 2.4.
- Packet separation at T: 12.00 sigma_0 (spec >= 8 sigma).
- Ensemble: N = 4000 positions from |phi|^2, seed 20260717.
- G2 invariance: the flip k -> -k maps psi_up <-> psi_down exactly (bitwise, in
  floating point: exp(i(-k)z) == exp(-ikz)); total density rho and total current
  j are component-symmetric, so the velocity field — and hence every trajectory —
  is *bitwise identical* between the two runs. Measured max |Z_flip(T) - Z_orig(T)| = 0.0
  and max |rho_flip(z,T) - rho_orig(z,T)| = 0.0, i.e. the invariance is exact at
  machine precision, not merely within integration tolerance.
- G4 time-convergence margin: 5.85e-9 absolute, ~4.9e-10 of the packet separation.

## Method

Two-component spinor Psi(z,0) = phi(z)(|up> + |down>)/sqrt(2) with phi a unit
Gaussian; the impulsive SG coupling (impulsive limit of H_int = -lambda z sigma_z)
multiplies the components by exp(+/- i k z). Subsequent free evolution (H = p^2/2)
by split-step FFT in half steps of dt/2, so the Bohmian velocity field
v = j_total/rho_total (current computed spectrally) is available at t, t+dt/2 and
t+dt without temporal interpolation; trajectories are advanced by classical RK4
with Catmull-Rom cubic interpolation of v in z. Crossing counting: sign changes of
Z between consecutive steps. Labels: at Z(T), whichever component density
|psi_up|^2 vs |psi_down|^2 dominates locally. The flipped run reuses the identical
Z(0) ensemble. G4 reruns the full evolution and integration at dt/2 for the first
100 trajectories.

## Connection to GUM

This exhibit is the Gleason-Jauch/non-contextuality violation Bell argued any
pilot-wave theory must display: the "spin measurement" outcome is not a readout
of a preexisting spin value but a joint product of the guided position and the
apparatus configuration — reverse the gradient and the same particle, on the
same trajectory, gets the opposite label. GUM Sec. III adopts the guidance
equation verbatim, so it inherits this contextuality wholesale. It is also the
conceptual basis for the corpus's claim (Sec. VIII.F) that its S3/A3
arrival-time gate statistics are apparatus-active: if outcome labels are
constituted by the detection arrangement rather than revealed by it, then
arrival-time statistics at a gate are properties of the particle-plus-apparatus
dynamics, not of the particle alone. Filed as finding exhibit **F-T6-2a-EXEC**.

## Caveats

- 1-D impulsive-kick idealization: the magnet is a momentum kick at t=0, not a
  finite-duration gradient region; Bell's conclusion is insensitive to this, but
  the trajectories inside a real SG field region would differ in detail.
- The exact bitwise invariance under k -> -k (max deviation 0.0) relies on the
  component-swap symmetry of this initial state being exact in floating point;
  a generic asymmetric spinor would show the flip only up to integration
  precision, and the exit *sides* themselves could change (only the symmetric
  case pins every side).
- G1's zero-crossing count is a theorem for the symmetric state (v(0,t) = 0 by
  antisymmetry of j); the numerical count of 0 confirms the integrator respects
  the symmetry but is not an independent physical result.
- Packet separation is measured between component-density centroids; residual
  overlap at T = 2.4 leaves the total density at z = 0 at 1.3e-3 of its peak
  value and does not affect any gate.
- Label assignment uses local component-density dominance at Z(T); the dominance
  ratio is > 10^3 for 99.9% of trajectories and >= 77 for the single worst
  trajectory (the one closest to z = 0, |Z(T)| = 0.88), so every label is
  unambiguous.

## Files

- `l1_flipped_sg.py` — simulation + gates
- `l1_figure.py` — figure generation
- `L1_results.json` — machine-readable results
- `L1_traj.npz` — stored trajectories/densities
- `L1_flipped_sg.png` — figure
