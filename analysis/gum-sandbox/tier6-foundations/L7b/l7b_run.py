#!/usr/bin/env python3
"""
WORKSTREAM L7b -- spin-family affinity test (T4-W5 POVM-exclusion testbed,
the A3-relevant leg).  File id: F-T6-L7b.

Physics: spin state chi(n-hat) on the Bloch sphere, NO magnetic field, so
psi = phi(x,z,t) chi with phi spin-independent -- ONE field evolution total.
Pauli guidance current:
    j = Im(phi* grad phi) + (n_y/2) (-d rho/dz, +d rho/dx),
so only the out-of-plane Bloch component n_y enters, and the velocity field
is exactly affine in n_y at the FIELD level:
    v(x,z,t; n_y) = a(x,z,t) + n_y b(x,z,t).
With the factorized exact-in-z state phi = f(x,t) sin(pi z) e^{-i pi^2 t/2}
(validated L5prime engine):
    v_x = Im(f'/f) - n_y * pi * cot(pi z)
    v_z = n_y * Re(f'/f)
which reproduces L5prime's transverse case at n_y = +1 and axial at n_y = 0.

The trajectory map (and hence the first-crossing statistics) is NOT affine
in n_y -- that is what this workstream measures.  A POVM would force
Pi_B(n_y) = A + B n_y per bin (affine on the Bloch sphere restricted by the
n_y-only dynamics); measured nonlinearity in n_y, or an interval of n_y with
Pi_B = 0 next to measured positivity, excludes every single POVM for the
family.

Configuration (inherited from validated L5prime, pre-registered here):
  k0 = 2, C2 smoothstep taper (W=1 for |x-x0|<=2.5 sigma, 0 at 3.5 sigma),
  literal amplitude exp(-(x-x0)^2/(2 sigma^2)), hard-wall waveguide z in
  [0,1] ground sine mode, box x in [-15, 45], reference Nx = 2048,
  dt = 2.5e-4 (<= 5e-4 spec), N = 2000, seed 20260717, d_near = 1, T = 16.
  n_y in {-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1}, identical initial
  ensemble for every n_y.

Refinements: double-Nx run (Nx = 4096, same box) for the G3 zero-bin
verification and tau_max stability; extended-box control ([-15, 105], same
dx) to bound FFT wrap-around contamination of late bins.

hbar = m = 1.
"""
import json
import os
import time

import numpy as np
from scipy import fft as sfft
from scipy import stats

WORKERS = os.cpu_count() or 1

# ----------------------------------------------------------------------
# Parameters (pre-registered)
# ----------------------------------------------------------------------
L = 1.0
SIGMA_X = 1.0
X0 = -5.0
K0 = 2.0
W_IN, W_OUT = 2.5, 3.5
D_NEAR = 1.0
T_FINAL = 16.0
N_PART = 2000
SEED = 20260717
VCLAMP = 500.0
RHO_EPS_REL = 1e-14

XMIN, XMAX = -15.0, 45.0          # pre-registered box (reference)
XMIN_E, XMAX_E = -15.0, 105.0     # extended box (wrap-around control, same dx)

NY_VALUES = np.array([-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0])
N_NY = NY_VALUES.size

N_BINS_T = 40                     # 40 uniform time bins on [0, 16] ...
BIN_EDGES = np.linspace(0.0, T_FINAL, N_BINS_T + 1)
N_BINS = N_BINS_T + 1             # ... plus the non-arrival bin (index 40)
N_BOOT = 2000
BOOT_SEED = 31415
SIGMA_FLOOR = float(np.sqrt((1.0 / N_PART) * (1.0 - 1.0 / N_PART) / N_PART))
# sigma of a single-count bin probability (3.53e-4); floor for zero-count bins

OUTDIR = os.path.dirname(os.path.abspath(__file__))
L5P_RAW = os.path.join(os.path.dirname(OUTDIR), "L5prime", "L5prime_raw.npz")

# reference numbers from L5prime (RESULTS.md / raw file)
L5P_TAU_MAX_T = 5.130950475155492   # extended-box reference, window 16


# ----------------------------------------------------------------------
# Initial profile (identical to L5prime "literal")
# ----------------------------------------------------------------------
def smoothstep(t):
    return t * t * t * (10.0 + t * (-15.0 + 6.0 * t))


def taper(x):
    u = np.abs(np.asarray(x, float) - X0) / SIGMA_X
    W = np.ones_like(u)
    W[u >= W_OUT] = 0.0
    m = (u > W_IN) & (u < W_OUT)
    W[m] = smoothstep((W_OUT - u[m]) / (W_OUT - W_IN))
    return W


def amp(x):
    d2 = (np.asarray(x, float) - X0) ** 2 / SIGMA_X ** 2
    return np.exp(-d2 / 2.0) * taper(x)


def sample_ensemble(seed):
    rng = np.random.default_rng(seed)
    xf = np.linspace(X0 - W_OUT * SIGMA_X, X0 + W_OUT * SIGMA_X, 700001)
    pdf = amp(xf) ** 2
    cdf = np.concatenate([[0.0],
                          np.cumsum(0.5 * (pdf[1:] + pdf[:-1]) * np.diff(xf))])
    cdf /= cdf[-1]
    xs = np.interp(rng.uniform(size=N_PART), cdf, xf)
    zfine = np.linspace(0.0, 1.0, 200001)
    zcdf = zfine - np.sin(2 * np.pi * zfine) / (2 * np.pi)
    zs = np.interp(rng.uniform(size=N_PART), zcdf, zfine)
    return np.column_stack([xs, zs])


# ----------------------------------------------------------------------
# Factorized exact-in-z engine, one field evolution -> 9 spin ensembles
# ----------------------------------------------------------------------
def build_f(Nx, xmin, xmax):
    dx = (xmax - xmin) / Nx
    xg = xmin + dx * np.arange(Nx)
    f = amp(xg).astype(complex)
    f *= np.exp(1j * K0 * xg)
    f /= np.sqrt(np.sum(np.abs(f) ** 2) * dx)
    kx = 2.0 * np.pi * sfft.fftfreq(Nx, dx)
    return xg, dx, f, kx


def lin_interp(V, px, xg0, dx, Nx):
    g = np.clip((px - xg0) / dx, 0.0, Nx - 1.000001)
    i = g.astype(np.int64)
    fr = g - i
    return V[i] * (1 - fr) + V[i + 1] * fr


def run_family(Nx, dt, xmin=XMIN, xmax=XMAX, t_final=T_FINAL, P0=None,
               label=""):
    """Evolve f once; advance all 9 n_y trajectory ensembles simultaneously.
    Velocity: v_x = a_x + n_y b_x with a_x = Im(f'/f) (clamped), b_x =
    -pi cot(pi z); v_z = n_y * b_z with b_z = Re(f'/f) (clamped); v_x
    re-clamped after the cot term, exactly as in L5prime (n_y=1 -> its
    transverse branch, n_y=0 -> its axial branch, bit-for-bit)."""
    xg, dx, f0, kx = build_f(Nx, xmin, xmax)
    F = sfft.fft(f0)
    Eh = np.exp(-0.5j * kx ** 2 * (dt / 2.0))
    nsteps = int(round(t_final / dt))
    norm0 = float(np.sum(np.abs(f0) ** 2) * dx)
    spec_pow = np.abs(F) ** 2
    E1 = float(np.sum(0.5 * kx ** 2 * spec_pow) / np.sum(spec_pow))

    def fields(F):
        f = sfft.ifft(F)
        fp = sfft.ifft(1j * kx * F)
        rho = np.abs(f) ** 2
        den = np.maximum(rho, RHO_EPS_REL * rho.max())
        w = np.conj(f) * fp
        a = np.clip(np.imag(w) / den, -VCLAMP, VCLAMP)
        b = np.clip(np.real(w) / den, -VCLAMP, VCLAMP)
        return a, b

    # flattened ensembles: identical initial points for every n_y
    NF = N_NY * N_PART
    ny_flat = np.repeat(NY_VALUES, N_PART)
    P = np.tile(P0, (N_NY, 1))
    active = np.ones(NF, bool)
    tau = np.full(NF, np.nan)

    def vel(ab, px, pz, ny):
        a, b = ab
        va = lin_interp(a, px, xg[0], dx, Nx)
        vb = lin_interp(b, px, xg[0], dx, Nx)
        vx = np.clip(va - ny * np.pi / np.tan(np.pi * pz / L),
                     -VCLAMP, VCLAMP)
        return vx, ny * vb

    n_check = 17
    check_every = max(1, nsteps // (n_check - 1))
    norms, energies = [], []
    ab = fields(F)
    t_start = time.time()
    for n in range(nsteps):
        t = n * dt
        if n % check_every == 0:
            f = sfft.ifft(F)
            norms.append(float(np.sum(np.abs(f) ** 2) * dx))
            p = np.abs(F) ** 2
            energies.append(float(np.sum(0.5 * kx ** 2 * p) / np.sum(p)))
        F = F * Eh
        ab_h = fields(F)
        F = F * Eh
        ab_n = fields(F)
        idx = np.where(active)[0]
        if idx.size:
            px, pz = P[idx, 0], P[idx, 1]
            ny = ny_flat[idx]
            k1x, k1z = vel(ab, px, pz, ny)
            k2x, k2z = vel(ab_h, px + 0.5 * dt * k1x,
                           np.clip(pz + 0.5 * dt * k1z, 1e-9, L - 1e-9), ny)
            k3x, k3z = vel(ab_h, px + 0.5 * dt * k2x,
                           np.clip(pz + 0.5 * dt * k2z, 1e-9, L - 1e-9), ny)
            k4x, k4z = vel(ab_n, px + dt * k3x,
                           np.clip(pz + dt * k3z, 1e-9, L - 1e-9), ny)
            nx = px + dt / 6.0 * (k1x + 2 * k2x + 2 * k3x + k4x)
            nz = np.clip(pz + dt / 6.0 * (k1z + 2 * k2z + 2 * k3z + k4z),
                         1e-9, L - 1e-9)
            cross = np.isnan(tau[idx]) & (px < D_NEAR) & (nx >= D_NEAR)
            if cross.any():
                ci = idx[cross]
                frac = (D_NEAR - px[cross]) / (nx[cross] - px[cross])
                tau[ci] = t + dt * frac
            P[idx, 0] = nx
            P[idx, 1] = nz
            esc = (nx < xmin + 0.5) | (nx > xmax - 0.5)
            if esc.any():
                active[idx[esc]] = False
        ab = ab_n
        if label and n % 16000 == 0:
            print("  [%s] step %d/%d  t=%.2f  active=%d  elapsed=%.0fs"
                  % (label, n, nsteps, t, active.sum(),
                     time.time() - t_start), flush=True)
    fT = sfft.ifft(F)
    norms.append(float(np.sum(np.abs(fT) ** 2) * dx))
    p = np.abs(F) ** 2
    energies.append(float(np.sum(0.5 * kx ** 2 * p) / np.sum(p)))
    norms = np.array(norms)
    energies = np.array(energies)
    return dict(tau=tau.reshape(N_NY, N_PART),
                active=active.reshape(N_NY, N_PART),
                E1=E1,
                norm_dev_max=float(np.max(np.abs(norms - norm0))),
                energy_drift_rel=float(np.max(np.abs(energies - E1))
                                       / abs(E1)),
                elapsed=time.time() - t_start)


# ----------------------------------------------------------------------
# 2-D wall-amplitude diagnostic (spectral checkpoints; no trajectories)
# ----------------------------------------------------------------------
def diag_2d(Nx=2048, Nz=128, xmin=XMIN, xmax=XMAX, t_final=T_FINAL,
            n_checks=17):
    dx = (xmax - xmin) / Nx
    xg = xmin + dx * np.arange(Nx)
    dz = L / Nz
    zg_ext = dz * np.arange(2 * Nz)
    fx = amp(xg).astype(complex) * np.exp(1j * K0 * xg)
    gz = np.sin(np.pi * zg_ext / L)
    phi0 = fx[:, None] * gz[None, :]
    dA = dx * dz
    phi0 /= np.sqrt(np.sum(np.abs(phi0[:, :Nz]) ** 2) * dA)
    kx = 2.0 * np.pi * sfft.fftfreq(Nx, dx)
    kz = 2.0 * np.pi * sfft.fftfreq(2 * Nz, dz)
    K2 = kx[:, None] ** 2 + kz[None, :] ** 2
    F0 = sfft.fft2(phi0, workers=WORKERS)
    norm0 = float(np.sum(np.abs(phi0[:, :Nz]) ** 2) * dA)
    Pspec = np.sum(np.abs(F0) ** 2)
    E0 = float(np.sum(0.5 * K2 * np.abs(F0) ** 2) / Pspec)
    rowpow = np.sum(np.abs(F0) ** 2, axis=0)
    mode1 = float(rowpow[1] + rowpow[2 * Nz - 1])
    leak = float((rowpow.sum() - mode1) / rowpow.sum())
    norms, walls, energies = [], [], []
    for tchk in np.linspace(0.0, t_final, n_checks):
        Ft = F0 * np.exp(-0.5j * K2 * tchk)
        phit = sfft.ifft2(Ft, workers=WORKERS)
        norms.append(float(np.sum(np.abs(phit[:, :Nz]) ** 2) * dA))
        walls.append(float(max(np.abs(phit[:, 0]).max(),
                               np.abs(phit[:, Nz]).max())))
        energies.append(float(np.sum(0.5 * K2 * np.abs(Ft) ** 2)
                              / np.sum(np.abs(Ft) ** 2)))
    return dict(norm_deviation_max=float(np.max(np.abs(np.array(norms)
                                                       - norm0))),
                wall_amplitude_max=float(np.max(walls)),
                energy_drift_rel=float(np.max(np.abs(np.array(energies)
                                                     - E0)) / abs(E0)),
                energy_mean=E0, z_mode_leakage=leak)


# ----------------------------------------------------------------------
# Binning / fits / bootstrap
# ----------------------------------------------------------------------
def bin_indices(tau_row):
    """Map each of the N_PART trajectories to a bin index 0..40
    (40 = non-arrival within [0, T])."""
    out = np.full(N_PART, N_BINS_T, dtype=np.int64)
    fin = np.isfinite(tau_row) & (tau_row <= T_FINAL)
    out[fin] = np.minimum((tau_row[fin] / (T_FINAL / N_BINS_T)).astype(np.int64),
                          N_BINS_T - 1)
    return out


def counts_matrix(tau):
    """(N_NY, N_BINS) count matrix."""
    C = np.zeros((N_NY, N_BINS), dtype=np.int64)
    for i in range(N_NY):
        C[i] = np.bincount(bin_indices(tau[i]), minlength=N_BINS)
    return C


def bootstrap_sigma(tau, n_boot=N_BOOT, seed=BOOT_SEED):
    """Per-(n_y, bin) bootstrap sigma of the bin probability, with SHARED
    resample indices across n_y (the ensembles share initial points).
    Returns sigma (N_NY, N_BINS) and the full resampled probability tensor
    (n_boot, N_NY, N_BINS) for the quadratic-coefficient bootstrap."""
    rng = np.random.default_rng(seed)
    bidx = np.stack([bin_indices(tau[i]) for i in range(N_NY)])  # (9, N)
    P_boot = np.empty((n_boot, N_NY, N_BINS))
    for r in range(n_boot):
        sel = rng.integers(0, N_PART, N_PART)
        rows = bidx[:, sel]  # (9, N)
        for i in range(N_NY):
            P_boot[r, i] = np.bincount(rows[i], minlength=N_BINS)
    P_boot /= N_PART
    return P_boot.std(axis=0, ddof=1), P_boot


def fit_family(Phat):
    """Least-squares fits per bin over the 9 n_y points.
    Affine: Pi = A + B ny.  Quadratic: Pi = A + B ny + C ny^2.
    Returns dict of coefficient arrays and residual matrices (N_NY, N_BINS)."""
    X1 = np.column_stack([np.ones(N_NY), NY_VALUES])            # (9, 2)
    X2 = np.column_stack([np.ones(N_NY), NY_VALUES, NY_VALUES ** 2])
    M1 = np.linalg.pinv(X1)                                     # (2, 9)
    M2 = np.linalg.pinv(X2)                                     # (3, 9)
    coef1 = M1 @ Phat            # (2, N_BINS)
    coef2 = M2 @ Phat            # (3, N_BINS)
    resid1 = Phat - X1 @ coef1   # (9, N_BINS)
    resid2 = Phat - X2 @ coef2
    return dict(X1=X1, X2=X2, M1=M1, M2=M2, coef1=coef1, coef2=coef2,
                resid1=resid1, resid2=resid2)


def main():
    P0 = sample_ensemble(SEED)

    # ensemble identity check vs L5prime raw
    l5p = np.load(L5P_RAW)
    p0_match = float(np.max(np.abs(P0 - l5p["P0"])))
    print("P0 max abs diff vs L5prime:", p0_match, flush=True)

    print("2-D wall diagnostic ...", flush=True)
    d2 = diag_2d()
    print("  ", json.dumps(d2), flush=True)

    runs = {}
    for label, Nx, dt, xmin, xmax in [
            ("ref",    2048, 2.5e-4, XMIN,   XMAX),    # reference
            ("refine", 4096, 2.5e-4, XMIN,   XMAX),    # double Nx
            ("extbox", 4096, 2.5e-4, XMIN_E, XMAX_E),  # wrap control, same dx
    ]:
        print("run %s: Nx=%d dt=%g box=[%g,%g]" % (label, Nx, dt, xmin, xmax),
              flush=True)
        runs[label] = run_family(Nx, dt, xmin, xmax, P0=P0, label=label)
        r = runs[label]
        print("  elapsed %.1f s  norm_dev %.3e  E_drift %.3e"
              % (r["elapsed"], r["norm_dev_max"], r["energy_drift_rel"]),
              flush=True)
        for i, ny in enumerate(NY_VALUES):
            t = r["tau"][i]
            fin = t[np.isfinite(t) & (t <= T_FINAL)]
            print("    ny=%+.2f  n_arr=%d  tau_max=%.6f  mean=%.4f"
                  % (ny, fin.size, fin.max(), fin.mean()), flush=True)

    ref = runs["ref"]
    tau = ref["tau"]

    # ------------------------------------------------------------------
    # G1: engine + L5prime cross-check
    # ------------------------------------------------------------------
    i_p1 = int(np.where(NY_VALUES == 1.0)[0][0])
    i_0 = int(np.where(NY_VALUES == 0.0)[0][0])
    t_p1 = tau[i_p1]
    t_0 = tau[i_0]
    w_p1 = t_p1[np.isfinite(t_p1) & (t_p1 <= T_FINAL)]
    w_0 = t_0[np.isfinite(t_0) & (t_0 <= T_FINAL)]
    l5T = l5p["tau_T_near"]
    l5A = l5p["tau_A_near"]
    l5T16 = l5T[np.isfinite(l5T) & (l5T <= T_FINAL)]
    l5A16 = l5A[np.isfinite(l5A) & (l5A <= T_FINAL)]
    tau_max_p1 = float(w_p1.max())
    rel_tau_max = abs(tau_max_p1 - L5P_TAU_MAX_T) / L5P_TAU_MAX_T
    ks_T = stats.ks_2samp(w_p1, l5T16)
    ks_A = stats.ks_2samp(w_0, l5A16)
    g1 = dict(
        p0_max_abs_diff_vs_L5prime=p0_match,
        norm_dev_1d=ref["norm_dev_max"],
        energy_drift_rel_1d=ref["energy_drift_rel"],
        energy_1d=ref["E1"],
        diag_2d=d2,
        tau_max_ny1=tau_max_p1,
        tau_max_L5prime=L5P_TAU_MAX_T,
        rel_tau_max_vs_L5prime=float(rel_tau_max),
        ks_ny1_vs_L5prime_T=dict(stat=float(ks_T.statistic),
                                 p=float(ks_T.pvalue),
                                 n=[int(w_p1.size), int(l5T16.size)]),
        ks_ny0_vs_L5prime_A=dict(stat=float(ks_A.statistic),
                                 p=float(ks_A.pvalue),
                                 n=[int(w_0.size), int(l5A16.size)]),
    )
    g1_pass = (ref["norm_dev_max"] < 1e-6
               and d2["wall_amplitude_max"] < 1e-10
               and ref["energy_drift_rel"] <= 1e-6
               and d2["norm_deviation_max"] < 1e-6
               and d2["energy_drift_rel"] <= 1e-6
               and rel_tau_max <= 1e-3
               and ks_T.pvalue > 0.5 and ks_A.pvalue > 0.5)
    g1["verdict"] = "PASS" if g1_pass else "FAIL"
    print("G1:", json.dumps(g1, indent=1), flush=True)

    # ------------------------------------------------------------------
    # G2: per-bin affine fit in n_y + bootstrap residual test
    # ------------------------------------------------------------------
    C = counts_matrix(tau)                       # (9, 41)
    Phat = C / float(N_PART)                     # (9, 41) probabilities
    print("bootstrap ...", flush=True)
    sigma_boot, P_boot = bootstrap_sigma(tau)
    sigma_eff = np.maximum(sigma_boot, SIGMA_FLOOR)
    fits = fit_family(Phat)
    ratio = np.abs(fits["resid1"]) / sigma_eff   # (9, 41)
    max_ratio_per_bin = ratio.max(axis=0)        # (41,)
    n_bins_over_5 = int(np.sum(max_ratio_per_bin > 5.0))
    # populated mid-range check: bins with >= 20 counts for every n_y
    min_count_per_bin = C.min(axis=0)
    populated = np.where(min_count_per_bin >= 20)[0]
    # quadratic diagnostic: C coefficient significance via bootstrap refits
    coef2_boot = np.einsum("cj,rjb->rcb", fits["M2"],
                           P_boot)               # (n_boot, 3, 41)
    sigma_C = coef2_boot[:, 2, :].std(axis=0, ddof=1)
    C_quad = fits["coef2"][2]
    C_signif = C_quad / np.maximum(sigma_C, 1e-300)
    g2 = dict(
        bin_edges=[float(b) for b in BIN_EDGES],
        nonarrival_bin_index=N_BINS_T,
        counts={"%+.2f" % NY_VALUES[i]: [int(c) for c in C[i]]
                for i in range(N_NY)},
        populated_bins_min20=[int(b) for b in populated],
        n_populated_bins_min20=int(populated.size),
        sigma_floor=SIGMA_FLOOR,
        affine_A=[float(v) for v in fits["coef1"][0]],
        affine_B=[float(v) for v in fits["coef1"][1]],
        max_resid_over_sigma_per_bin=[float(v) for v in max_ratio_per_bin],
        n_bins_over_5sigma=n_bins_over_5,
        bins_over_5sigma=[int(b) for b in np.where(max_ratio_per_bin > 5)[0]],
        max_violation_sigma=float(max_ratio_per_bin.max()),
        max_violation_bin=int(np.argmax(max_ratio_per_bin)),
        quad_C=[float(v) for v in C_quad],
        quad_C_signif=[float(v) for v in C_signif],
        max_quad_C_signif_abs=float(np.max(np.abs(C_signif))),
        n_bins_quad_C_over_5=int(np.sum(np.abs(C_signif) > 5.0)),
    )
    g2["verdict"] = "PASS" if n_bins_over_5 >= 3 else "FAIL"
    print("G2: n_bins_over_5sigma=%d max=%.1f sigma  quadC>5sig in %d bins"
          % (n_bins_over_5, g2["max_violation_sigma"],
             g2["n_bins_quad_C_over_5"]), flush=True)

    # ------------------------------------------------------------------
    # G3: affine-vanishing kill
    # ------------------------------------------------------------------
    C_refine = counts_matrix(runs["refine"]["tau"])
    C_ext = counts_matrix(runs["extbox"]["tau"])
    p95_zero = 1.0 - 0.05 ** (1.0 / N_PART)   # exact 95% CL upper bound on p
    # (0 of 2000 -> p < 1.497e-3 at 95% CL; < 1.9e-3 spec satisfied)
    kill_bins = []
    for nstar in (0.75, 0.5):
        outer = np.abs(NY_VALUES) >= nstar - 1e-12
        inner = ~outer
        for b in range(N_BINS_T):
            if np.all(C[outer, b] == 0) and np.any(Phat[inner, b] > 1e-2):
                kill_bins.append(dict(
                    bin=int(b),
                    t_range=[float(BIN_EDGES[b]), float(BIN_EDGES[b + 1])],
                    nstar=float(nstar),
                    zero_ny=[float(v) for v in NY_VALUES[outer]],
                    counts_all=[int(c) for c in C[:, b]],
                    max_inner_p=float(Phat[inner, b].max()),
                    max_inner_ny=float(NY_VALUES[inner][
                        np.argmax(Phat[inner, b])]),
                    zero_confirmed_refine=bool(np.all(C_refine[outer, b] == 0)),
                    zero_confirmed_extbox=bool(np.all(C_ext[outer, b] == 0)),
                ))
        if kill_bins:
            break
    # a kill bin passes if zeros hold at reference AND the double-Nx
    # refinement (wrap control reported alongside)
    confirmed = [k for k in kill_bins
                 if k["zero_confirmed_refine"] and k["zero_confirmed_extbox"]]
    g3 = dict(p_upper_95CL_zero_of_2000=float(p95_zero),
              kill_bins=kill_bins,
              n_kill_bins=len(kill_bins),
              n_confirmed=len(confirmed),
              lemma=("An affine function Pi_B(n_y) = A + B n_y that vanishes "
                     "at >= 2 distinct n_y values (a fortiori on an interval) "
                     "is identically zero; measured Pi_B > 1e-2 at an "
                     "interior n_y therefore contradicts every POVM. "
                     "Quantitatively: affine with Pi_B(+/-1) < 1.5e-3 (95% "
                     "CL) forces Pi_B(0) = (Pi_B(-1)+Pi_B(+1))/2 < 1.5e-3, "
                     "excluded by the measured interior value."))
    g3["verdict"] = "PASS" if len(confirmed) >= 1 else "FAIL"
    print("G3: kill bins %d (confirmed %d)" % (len(kill_bins), len(confirmed)),
          flush=True)

    # ------------------------------------------------------------------
    # G4: tau_max(n_y) profile + gap check per n_y
    # ------------------------------------------------------------------
    prof = []
    for i, ny in enumerate(NY_VALUES):
        t = tau[i]
        fin = t[np.isfinite(t) & (t <= T_FINAL)]
        tr = runs["refine"]["tau"][i]
        finr = tr[np.isfinite(tr) & (tr <= T_FINAL)]
        te = runs["extbox"]["tau"][i]
        fine = te[np.isfinite(te) & (te <= T_FINAL)]
        iqr = float(np.subtract(*np.percentile(fin, [75, 25])))
        tmax = float(fin.max())
        prof.append(dict(ny=float(ny), n_arrivals=int(fin.size),
                         n_nonarrival=int(N_PART - fin.size),
                         tau_max=tmax,
                         tau_max_refine=float(finr.max()),
                         tau_max_extbox=float(fine.max()),
                         n_arrivals_extbox=int(fine.size),
                         mean=float(fin.mean()),
                         iqr=iqr, gap=float(T_FINAL - tmax),
                         gap_ge_2iqr=bool(T_FINAL - tmax >= 2 * iqr),
                         frozen_escapees=int(np.sum(~ref["active"][i]))))
    tmaxs = np.array([p["tau_max"] for p in prof])
    variation = float((tmaxs.max() - tmaxs.min()) / tmaxs.min())
    # monotone in |n_y|: average the +/- branches, check decreasing in |n_y|
    absvals = np.unique(np.abs(NY_VALUES))
    tmax_by_abs = []
    for a in absvals:
        sel = np.abs(NY_VALUES) == a
        tmax_by_abs.append(float(tmaxs[sel].mean()))
    mono = bool(np.all(np.diff(tmax_by_abs) < 0))
    branch_p = tmaxs[NY_VALUES > 0]
    branch_m = tmaxs[NY_VALUES < 0][::-1]
    g4 = dict(profile=prof, variation_rel=variation,
              tau_max_by_absny={"%.2f" % a: v
                               for a, v in zip(absvals, tmax_by_abs)},
              monotone_decreasing_in_absny=mono,
              monotone_branch_pos=bool(np.all(np.diff(branch_p) < 0)),
              monotone_branch_neg=bool(np.all(np.diff(branch_m) < 0)))
    g4["verdict"] = "PASS" if (variation > 0.2 and mono) else \
        ("PARTIAL" if variation > 0.2 else "FAIL")
    print("G4: variation %.2f monotone(|ny|) %s" % (variation, mono),
          flush=True)

    # ------------------------------------------------------------------
    # G5: report-only effect sizes
    # ------------------------------------------------------------------
    means = {"%+.2f" % p["ny"]: p["mean"] for p in prof}
    g5 = dict(max_linearity_violation_sigma=g2["max_violation_sigma"],
              max_violation_bin_t_range=[
                  float(BIN_EDGES[g2["max_violation_bin"]]),
                  float(BIN_EDGES[min(g2["max_violation_bin"] + 1,
                                      N_BINS_T)])]
              if g2["max_violation_bin"] < N_BINS_T else "nonarrival",
              mean_arrival_by_ny=means,
              mean_arrival_range=[float(min(means.values())),
                                  float(max(means.values()))],
              max_quad_C_signif=g2["max_quad_C_signif_abs"],
              verdict="REPORT-ONLY")

    results = dict(
        workstream="L7b",
        file_id="F-T6-L7b",
        params=dict(L=L, sigma_x=SIGMA_X, x0=X0, k0=K0,
                    taper="S2 smoothstep, W=1 for |x-x0|<=2.5 sigma, 0 at 3.5",
                    amplitude="literal exp(-(x-x0)^2/(2 sigma^2)) * W",
                    box=[XMIN, XMAX], box_extended_control=[XMIN_E, XMAX_E],
                    Nx_ref=2048, dt_ref=2.5e-4, T=T_FINAL, N=N_PART,
                    seed=SEED, d_near=D_NEAR, ny_values=[float(v) for v in
                                                         NY_VALUES],
                    n_bins_time=N_BINS_T, n_boot=N_BOOT, boot_seed=BOOT_SEED,
                    vclamp=VCLAMP, rho_eps_rel=RHO_EPS_REL),
        G1=g1, G2=g2, G3=g3, G4=g4, G5=g5,
        run_elapsed={k: runs[k]["elapsed"] for k in runs},
    )
    with open(os.path.join(OUTDIR, "l7b_results.json"), "w") as fh:
        json.dump(results, fh, indent=2)
    np.savez_compressed(
        os.path.join(OUTDIR, "l7b_raw.npz"),
        P0=P0, ny_values=NY_VALUES,
        tau_ref=tau, tau_refine=runs["refine"]["tau"],
        tau_extbox=runs["extbox"]["tau"],
        counts_ref=C, counts_refine=C_refine, counts_extbox=C_ext,
        Phat=Phat, sigma_boot=sigma_boot, sigma_eff=sigma_eff,
        affine_coef=fits["coef1"], quad_coef=fits["coef2"],
        resid_affine=fits["resid1"], ratio=ratio,
        quad_C_signif=C_signif, bin_edges=BIN_EDGES)
    print("verdicts:", {g: results[g]["verdict"]
                        for g in ("G1", "G2", "G3", "G4", "G5")})


if __name__ == "__main__":
    main()
