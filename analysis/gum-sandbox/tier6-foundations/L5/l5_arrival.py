#!/usr/bin/env python3
"""
WORKSTREAM L5 -- A3 arrival-time phenomenon on an independent engine.

Spin-dependent Bohmian arrival times in a 2-D hard-wall waveguide
(Das-Durr class, Sci. Rep. 9, 2242 (2019)).

Model (hbar = m = 1):
  psi = phi(x,z,t) chi, chi fixed (no magnetic field).
  Waveguide walls at z=0, z=L (L=1): sine basis via odd extension in z,
  full 2-D FFT on the (Nx x 2*Nz) extended array.  Free Hamiltonian is
  diagonal in that basis, so time evolution is EXACT spectral propagation
  (phase multiplication); dt is only the trajectory integrator step.
  Initial state: truncated Gaussian in x (sigma_x=1, cut at +/-3 sigma,
  renormalized, centered x0=-5, boost e^{i k0 x}, k0=4) x sin(pi z/L).

Guidance law:
  j_conv = Im(phi* grad phi)
  transverse spin (sigma_y eigenstate, M_y = +rho/2):
      j_spin = ( -(1/2) d rho/dz , +(1/2) d rho/dx )
  axial spin (sigma_x / sigma_z eigenstate): j_spin = 0 in-plane.
  v = j / rho.

Trajectories: N=2000 samples from |phi_0|^2 (fixed seed), RK4 with
bilinear interpolation of stored velocity fields at t, t+dt/2, t+dt.
Arrival time = first crossing of x = d (d_near=1, d_far=25); particles
that never cross within T are recorded as non-arrivals.
Both spin cases use the IDENTICAL initial ensemble and share phi.
"""
import json
import os
import time

import numpy as np
from scipy import fft as sfft
from scipy import stats

WORKERS = os.cpu_count() or 1

# ----------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------
L = 1.0
XMIN, XMAX = -15.0, 45.0
NX = 1024
NZ = 128                    # z intervals in [0, L]; extended array 2*NZ
SIGMA_X = 1.0
X0 = -5.0
K0 = 4.0
TRUNC = 3.0 * SIGMA_X
DT = 2.0e-3
T_FINAL = 8.0
N_PART = 2000
SEED = 12345
D_NEAR = 1.0
D_FAR = 25.0
VCLAMP = 500.0              # safety clamp on grid velocities (near-wall 1/z)
RHO_EPS_REL = 1e-14         # relative floor for rho in v = j/rho

OUTDIR = os.path.dirname(os.path.abspath(__file__))


def build_state():
    """Return grids, initial extended array phi0_ext, and spectral k arrays."""
    dx = (XMAX - XMIN) / NX
    xg = XMIN + dx * np.arange(NX)
    dz = L / NZ
    # extended z grid [0, 2L): indices 0..2NZ-1 ; physical rows 0..NZ
    zg_ext = dz * np.arange(2 * NZ)

    fx = np.exp(-(xg - X0) ** 2 / (4.0 * SIGMA_X ** 2)).astype(complex)
    fx[np.abs(xg - X0) > TRUNC] = 0.0          # hard truncation (compact support)
    fx *= np.exp(1j * K0 * xg)

    gz = np.sin(np.pi * zg_ext / L)            # on [0,2L) this IS the odd extension

    phi0 = fx[:, None] * gz[None, :]
    # normalize over the physical half (rows 0..NZ-1; row NZ is the z=L wall, =0)
    nrm = np.sqrt(np.sum(np.abs(phi0[:, :NZ]) ** 2) * dx * dz)
    phi0 /= nrm

    kx = 2.0 * np.pi * sfft.fftfreq(NX, dx)
    kz = 2.0 * np.pi * sfft.fftfreq(2 * NZ, dz)
    KX = kx[:, None]
    KZ = kz[None, :]
    K2 = KX ** 2 + KZ ** 2
    return xg, dx, zg_ext, dz, phi0, KX, KZ, K2


def fields_from_spectrum(F, KX, KZ):
    """phi, dphi/dx, dphi/dz on the extended grid from spectrum F."""
    stack = np.stack([F, 1j * KX * F, 1j * KZ * F])
    out = sfft.ifft2(stack, axes=(-2, -1), workers=WORKERS)
    return out[0], out[1], out[2]


def velocity_fields(F, KX, KZ):
    """Return (vxT, vzT, vxA, vzA) on the physical grid rows 0..NZ (NZ+1 rows),
    plus wall amplitude diagnostic."""
    phi, dpx, dpz = fields_from_spectrum(F, KX, KZ)
    ph = phi[:, : NZ + 1]
    dx_ = dpx[:, : NZ + 1]
    dz_ = dpz[:, : NZ + 1]
    rho = np.abs(ph) ** 2
    jx_c = np.imag(np.conj(ph) * dx_)
    jz_c = np.imag(np.conj(ph) * dz_)
    drho_dx = 2.0 * np.real(np.conj(ph) * dx_)
    drho_dz = 2.0 * np.real(np.conj(ph) * dz_)
    denom = np.maximum(rho, RHO_EPS_REL * rho.max())
    vxA = np.clip(jx_c / denom, -VCLAMP, VCLAMP)
    vzA = np.clip(jz_c / denom, -VCLAMP, VCLAMP)
    vxT = np.clip((jx_c - 0.5 * drho_dz) / denom, -VCLAMP, VCLAMP)
    vzT = np.clip((jz_c + 0.5 * drho_dx) / denom, -VCLAMP, VCLAMP)
    wall = max(np.abs(phi[:, 0]).max(), np.abs(phi[:, NZ]).max())
    return (vxT, vzT, vxA, vzA), wall


def bilinear(V, px, pz, dx, dz):
    """Vectorized bilinear interpolation of grid field V (NX, NZ+1)."""
    gx = np.clip((px - XMIN) / dx, 0.0, NX - 1.000001)
    gz = np.clip(pz / dz, 0.0, NZ - 1e-6)
    i = gx.astype(np.int64)
    j = gz.astype(np.int64)
    fx = gx - i
    fz = gz - j
    return (V[i, j] * (1 - fx) * (1 - fz) + V[i + 1, j] * fx * (1 - fz)
            + V[i, j + 1] * (1 - fx) * fz + V[i + 1, j + 1] * fx * fz)


def sample_ensemble(rng, xg, dx):
    """N samples from |phi_0|^2 (product form: truncated Gaussian x, 2 sin^2 z)."""
    # x: rejection sampling of N(X0, SIGMA_X) truncated at +/- 3 sigma
    xs = np.empty(0)
    while xs.size < N_PART:
        cand = rng.normal(X0, SIGMA_X, size=2 * N_PART)
        cand = cand[np.abs(cand - X0) <= TRUNC]
        xs = np.concatenate([xs, cand])
    xs = xs[:N_PART]
    # z: inverse CDF of pdf 2 sin^2(pi z);  CDF(z) = z - sin(2 pi z)/(2 pi)
    zfine = np.linspace(0.0, 1.0, 200001)
    cdf = zfine - np.sin(2 * np.pi * zfine) / (2 * np.pi)
    u = rng.uniform(0.0, 1.0, size=N_PART)
    zs = np.interp(u, cdf, zfine)
    return np.column_stack([xs, zs])


def run_trajectories(dt, K2, KX, KZ, phi0, dx, dz, P0, t_final=T_FINAL,
                     collect_snapshots=False):
    """Advance transverse (T) and axial (A) ensembles with RK4 in the exact
    spectrally-propagated velocity field.  Returns arrival-time arrays and
    diagnostics."""
    nsteps = int(round(t_final / dt))
    Ehalf = np.exp(-0.5j * K2 * (dt / 2.0))
    F = sfft.fft2(phi0, workers=WORKERS)

    ens = {"T": P0.copy(), "A": P0.copy()}
    active = {"T": np.ones(N_PART, bool), "A": np.ones(N_PART, bool)}
    tau = {s: {d: np.full(N_PART, np.nan) for d in ("near", "far")} for s in "TA"}
    dets = {"near": D_NEAR, "far": D_FAR}

    vf, wall0 = velocity_fields(F, KX, KZ)   # fields at t=0
    wall_max = wall0
    snaps = []
    traj_idx = np.arange(0, N_PART, 25)       # 80 traced particles
    traj = {"t": [], "T": [], "A": []} if collect_snapshots else None
    if collect_snapshots:
        traj["t"].append(0.0)
        traj["T"].append(ens["T"][traj_idx].copy())
        traj["A"].append(ens["A"][traj_idx].copy())
    t = 0.0
    t_start = time.time()
    for n in range(nsteps):
        F = F * Ehalf
        vf_h, w1 = velocity_fields(F, KX, KZ)      # t + dt/2
        F = F * Ehalf
        vf_n, w2 = velocity_fields(F, KX, KZ)      # t + dt
        wall_max = max(wall_max, w1, w2)

        for s, (ivx, ivz) in (("T", (0, 1)), ("A", (2, 3))):
            P = ens[s]
            act = active[s]
            if not act.any():
                continue
            px, pz = P[act, 0], P[act, 1]
            k1x = bilinear(vf[ivx], px, pz, dx, dz)
            k1z = bilinear(vf[ivz], px, pz, dx, dz)
            ax, az = px + 0.5 * dt * k1x, pz + 0.5 * dt * k1z
            k2x = bilinear(vf_h[ivx], ax, az, dx, dz)
            k2z = bilinear(vf_h[ivz], ax, az, dx, dz)
            bx, bz = px + 0.5 * dt * k2x, pz + 0.5 * dt * k2z
            k3x = bilinear(vf_h[ivx], bx, bz, dx, dz)
            k3z = bilinear(vf_h[ivz], bx, bz, dx, dz)
            cx, cz = px + dt * k3x, pz + dt * k3z
            k4x = bilinear(vf_n[ivx], cx, cz, dx, dz)
            k4z = bilinear(vf_n[ivz], cx, cz, dx, dz)
            nx = px + dt / 6.0 * (k1x + 2 * k2x + 2 * k3x + k4x)
            nz = pz + dt / 6.0 * (k1z + 2 * k2z + 2 * k3z + k4z)
            nz = np.clip(nz, 1e-9, L - 1e-9)

            # arrival detection: first upward crossing of each detector line
            idx = np.where(act)[0]
            for dname, dpos in dets.items():
                ta = tau[s][dname]
                cross = np.isnan(ta[idx]) & (px < dpos) & (nx >= dpos)
                if cross.any():
                    ci = idx[cross]
                    frac = (dpos - px[cross]) / (nx[cross] - px[cross])
                    ta[ci] = t + dt * frac
            P[idx, 0] = nx
            P[idx, 1] = nz
            # freeze escapees (left/right margins of the periodic box)
            esc = (nx < XMIN + 0.5) | (nx > XMAX - 0.5)
            if esc.any():
                active[s][idx[esc]] = False
        vf = vf_n
        t = (n + 1) * dt
        if collect_snapshots and (n + 1) % max(1, nsteps // 8) == 0:
            snaps.append((t, ens["T"].copy(), ens["A"].copy()))
        if collect_snapshots and (n + 1) % 20 == 0:
            traj["t"].append(t)
            traj["T"].append(ens["T"][traj_idx].copy())
            traj["A"].append(ens["A"][traj_idx].copy())
    elapsed = time.time() - t_start
    norm_final = float(np.sum(np.abs(sfft.ifft2(F, workers=WORKERS)[:, :NZ]) ** 2)
                       * dx * dz)
    return dict(tau=tau, wall_max=wall_max, norm_final=norm_final,
                active=active, snaps=snaps, traj=traj, traj_idx=traj_idx,
                elapsed=elapsed, ens=ens)


def main():
    xg, dx, zg_ext, dz, phi0, KX, KZ, K2 = build_state()

    # ---------------- G1 diagnostics (exact propagation) ----------------
    F0 = sfft.fft2(phi0, workers=WORKERS)
    dA = dx * dz  # physical-cell measure; extended array double-counts by 2
    norm0 = float(np.sum(np.abs(phi0[:, :NZ]) ** 2) * dA)
    P0spec = np.sum(np.abs(F0) ** 2)
    E0 = float(np.sum(0.5 * K2 * np.abs(F0) ** 2) / P0spec)
    # z-mode purity: power per kz row (mode n lives at kz index n and 2NZ-n)
    rowpow = np.sum(np.abs(F0) ** 2, axis=0)
    mode1 = float(rowpow[1] + rowpow[2 * NZ - 1])
    leak = float((rowpow.sum() - mode1) / rowpow.sum())

    norms, walls, energies = [], [], []
    for tchk in np.linspace(0.0, T_FINAL, 17):
        Ft = F0 * np.exp(-0.5j * K2 * tchk)
        phit = sfft.ifft2(Ft, workers=WORKERS)
        norms.append(float(np.sum(np.abs(phit[:, :NZ]) ** 2) * dA))
        walls.append(float(max(np.abs(phit[:, 0]).max(),
                               np.abs(phit[:, NZ]).max())))
        energies.append(float(np.sum(0.5 * K2 * np.abs(Ft) ** 2)
                              / np.sum(np.abs(Ft) ** 2)))
    norm_dev = float(np.max(np.abs(np.array(norms) - norm0)))
    wall_max_diag = float(np.max(walls))
    e_drift = float(np.max(np.abs(np.array(energies) - E0)) / abs(E0))

    # ---------------- ensembles and main run ----------------
    rng = np.random.default_rng(SEED)
    P0 = sample_ensemble(rng, xg, dx)

    print("main run: dt =", DT)
    main_run = run_trajectories(DT, K2, KX, KZ, phi0, dx, dz, P0,
                                collect_snapshots=True)
    print("  elapsed %.1f s" % main_run["elapsed"])

    print("refinement run: dt =", DT / 2)
    ref_run = run_trajectories(DT / 2, K2, KX, KZ, phi0, dx, dz, P0)
    print("  elapsed %.1f s" % ref_run["elapsed"])

    # ---------------- gate evaluation ----------------
    def arr(run, s, d):
        a = run["tau"][s][d]
        return a[np.isfinite(a)]

    tT_near = arr(main_run, "T", "near")
    tA_near = arr(main_run, "A", "near")
    tT_far = arr(main_run, "T", "far")
    tA_far = arr(main_run, "A", "far")
    tT_near_ref = arr(ref_run, "T", "near")

    # G2: transverse hard cutoff
    tau_max = float(tT_near.max())
    iqr_T = float(np.subtract(*np.percentile(tT_near, [75, 25])))
    gap = T_FINAL - tau_max
    tau_max_ref = float(tT_near_ref.max())
    tau_max_rel = abs(tau_max_ref - tau_max) / tau_max
    g2_pass = (gap >= 2.0 * iqr_T) and (tau_max_rel <= 0.02)

    # G3: axial contrast
    n_axial_late = int(np.sum(tA_near > tau_max))
    frac_axial_late = n_axial_late / max(1, tA_near.size)
    g3_pass = frac_axial_late > 0.05

    # G4: far-field KS
    ks = stats.ks_2samp(tT_far, tA_far)
    g4_pass = ks.pvalue > 0.05

    # G5: near-field mean ratio with bootstrap
    rng_b = np.random.default_rng(777)
    ratio = float(tT_near.mean() / tA_near.mean())
    boots = np.empty(2000)
    for b in range(2000):
        rT = rng_b.choice(tT_near, tT_near.size)
        rA = rng_b.choice(tA_near, tA_near.size)
        boots[b] = rT.mean() / rA.mean()
    ratio_err = float(boots.std(ddof=1))

    g1_pass = (norm_dev < 1e-6) and (wall_max_diag < 1e-10) and (e_drift <= 1e-6)

    results = {
        "workstream": "L5",
        "file_id": "F-T6-L5-EXEC",
        "params": dict(L=L, xmin=XMIN, xmax=XMAX, Nx=NX, Nz=NZ,
                       sigma_x=SIGMA_X, x0=X0, k0=K0, trunc_sigma=3.0,
                       dt=DT, dt_refined=DT / 2, T=T_FINAL, N=N_PART,
                       seed=SEED, d_near=D_NEAR, d_far=D_FAR,
                       vclamp=VCLAMP, rho_eps_rel=RHO_EPS_REL),
        "G1": dict(norm_deviation_max=norm_dev,
                   wall_amplitude_max=wall_max_diag,
                   wall_amplitude_max_traj_run=float(main_run["wall_max"]),
                   energy_drift_rel=e_drift,
                   energy_mean=E0,
                   z_mode_leakage=leak,
                   verdict="PASS" if g1_pass else "FAIL"),
        "G2": dict(tau_max=tau_max, tau_max_refined=tau_max_ref,
                   tau_max_rel_change=tau_max_rel,
                   iqr_transverse_near=iqr_T,
                   empty_gap=gap, required_gap=2.0 * iqr_T,
                   arrivals_in_gap=0,
                   verdict="PASS" if g2_pass else "FAIL"),
        "G3": dict(frac_axial_beyond_tau_max=frac_axial_late,
                   n_axial_beyond=n_axial_late,
                   n_axial_arrivals=int(tA_near.size),
                   axial_last_arrival=float(tA_near.max()),
                   verdict="PASS" if g3_pass else "FAIL"),
        "G4": dict(ks_stat=float(ks.statistic), ks_p=float(ks.pvalue),
                   n_far_T=int(tT_far.size), n_far_A=int(tA_far.size),
                   verdict="PASS" if g4_pass else "FAIL"),
        "G5": dict(mean_tau_T_near=float(tT_near.mean()),
                   mean_tau_A_near=float(tA_near.mean()),
                   ratio_T_over_A=ratio, ratio_bootstrap_err=ratio_err,
                   verdict="REPORT-ONLY"),
        "counts": dict(
            near_T_arrivals=int(tT_near.size), near_A_arrivals=int(tA_near.size),
            far_T_arrivals=int(tT_far.size), far_A_arrivals=int(tA_far.size),
            nonarrival_near_T=int(N_PART - tT_near.size),
            nonarrival_near_A=int(N_PART - tA_near.size),
            nonarrival_far_T=int(N_PART - tT_far.size),
            nonarrival_far_A=int(N_PART - tA_far.size),
            frozen_escapees_T=int(np.sum(~main_run["active"]["T"])),
            frozen_escapees_A=int(np.sum(~main_run["active"]["A"])),
        ),
        "quartiles": dict(
            T_near=[float(q) for q in np.percentile(tT_near, [25, 50, 75])],
            A_near=[float(q) for q in np.percentile(tA_near, [25, 50, 75])],
            T_far=[float(q) for q in np.percentile(tT_far, [25, 50, 75])],
            A_far=[float(q) for q in np.percentile(tA_far, [25, 50, 75])],
        ),
    }

    with open(os.path.join(OUTDIR, "L5_results.json"), "w") as fh:
        json.dump(results, fh, indent=2)

    np.savez_compressed(
        os.path.join(OUTDIR, "L5_raw.npz"),
        tau_T_near=main_run["tau"]["T"]["near"],
        tau_A_near=main_run["tau"]["A"]["near"],
        tau_T_far=main_run["tau"]["T"]["far"],
        tau_A_far=main_run["tau"]["A"]["far"],
        tau_T_near_refined=ref_run["tau"]["T"]["near"],
        P0=P0,
        final_T=main_run["ens"]["T"], final_A=main_run["ens"]["A"],
    )
    # traced trajectories for the figure panel
    if main_run["traj"] is not None:
        tr = main_run["traj"]
        np.savez_compressed(os.path.join(OUTDIR, "L5_traj.npz"),
                            ts=np.array(tr["t"]),
                            T=np.stack(tr["T"]), A=np.stack(tr["A"]),
                            idx=main_run["traj_idx"])

    print(json.dumps({k: results[k] for k in
                      ("G1", "G2", "G3", "G4", "G5", "counts")}, indent=2))


if __name__ == "__main__":
    main()
