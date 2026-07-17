#!/usr/bin/env python3
"""
WORKSTREAM L5' -- pre-registered follow-up converting L5's two spec-artifact
PARTIALs (F-T6-L5-EXEC-2).

Same physics as L5 (spin-dependent Bohmian arrival times, 2-D hard-wall
waveguide, Das-Durr class), with the two pre-registered fixes:

  1. C1+-smoothed compact truncation: x-amplitude
         f0(x) = exp(-(x-x0)^2 / (2 sigma_x^2)) * W(x) * e^{i k0 x}
     with W the polynomial smoothstep S2(t) = 6t^5 - 15t^4 + 10t^3 taper:
     W = 1 for |x-x0| <= 2.5 sigma_x, monotone to 0 at |x-x0| = 3.5 sigma_x,
     zero beyond (C2 wave function at the support edge; compact support kept).
     NOTE: the pre-registered amplitude exponent is 1/(2 sigma_x^2) (position
     density std sigma_x/sqrt2), whereas L5's code used 1/(4 sigma_x^2)
     (density std sigma_x).  We follow the pre-registered formula literally
     and run the L5-convention amplitude exp(-d^2/(4 sigma^2)) * W as a
     sensitivity check ("alt" run).
  2. Slower packet: k0 = 2 (was 4).  d_near = 1, d_far = 25 kept.
     Evolve to T = 24; near-field analysis window T_near = 16.

Engines (both inherited from L5, revalidated against each other here):
  * factorized exact-in-z engine (primary; the state is exactly
    f(x,t) sin(pi z) e^{-i pi^2 t/2}, f propagated exactly in 1-D Fourier
    space; guidance: transverse dx/dt = Im(f'/f) - pi cot(pi z),
    dz/dt = Re(f'/f); axial dx/dt = Im(f'/f), dz/dt = 0).
  * full 2-D FFT engine (cross-check; also supplies wall-amplitude G1).

Domain: the pre-registered box x in [-15, 45] cannot contain the packet to
T = 24 at k0 = 2 (front ~ x0 + (k0 + 3 sigma_v) T ~ 93; the FFT box is
periodic, so overrun wraps).  Reference runs therefore use an extended box
[-15, 105] at the SAME dx (Nx = 4096 on the extended box == the
pre-registered Nx = 2048 on [-15, 45]).  A control run on the pre-registered
box quantifies the wrap-around effect; this is a containment fix, not a
physics change.

hbar = m = 1.  Seed 20260717.  N = 2000, identical ensemble for both spin
cases and all resolutions.
"""
import json
import os
import time

import numpy as np
from scipy import fft as sfft
from scipy import stats

WORKERS = os.cpu_count() or 1

# ----------------------------------------------------------------------
# Pre-registered parameters
# ----------------------------------------------------------------------
L = 1.0
SIGMA_X = 1.0
X0 = -5.0
K0 = 2.0                       # pre-registered fix 2 (was 4)
W_IN, W_OUT = 2.5, 3.5         # taper: W=1 for u<=2.5, W=0 at u>=3.5 (u=|x-x0|/sigma)
D_NEAR = 1.0
D_FAR = 25.0
T_EVOLVE = 24.0                # far-field evolution horizon
T_NEAR = 16.0                  # near-field analysis window
N_PART = 2000
SEED = 20260717
VCLAMP = 500.0
RHO_EPS_REL = 1e-14

XMIN_P, XMAX_P = -15.0, 45.0    # pre-registered box (control run)
XMIN_E, XMAX_E = -15.0, 105.0   # extended box (containment; same dx at 2x Nx)

OUTDIR = os.path.dirname(os.path.abspath(__file__))


# ----------------------------------------------------------------------
# Initial profile
# ----------------------------------------------------------------------
def smoothstep(t):
    """S2(t) = 6t^5 - 15t^4 + 10t^3 on [0,1]."""
    return t * t * t * (10.0 + t * (-15.0 + 6.0 * t))


def taper(x):
    u = np.abs(np.asarray(x, float) - X0) / SIGMA_X
    W = np.ones_like(u)
    W[u >= W_OUT] = 0.0
    m = (u > W_IN) & (u < W_OUT)
    W[m] = smoothstep((W_OUT - u[m]) / (W_OUT - W_IN))
    return W


def amp(x, convention="literal"):
    """Real envelope (no boost).  'literal' = pre-registered exp(-d^2/2s^2);
    'alt' = L5 code convention exp(-d^2/4s^2).  Both times the S2 taper."""
    d2 = (np.asarray(x, float) - X0) ** 2 / SIGMA_X ** 2
    g = np.exp(-d2 / 2.0) if convention == "literal" else np.exp(-d2 / 4.0)
    return g * taper(x)


def sample_ensemble(seed, convention="literal"):
    """N samples from |phi_0|^2: x by inverse CDF of amp^2 (analytic, grid-
    independent), z by inverse CDF of 2 sin^2(pi z)."""
    rng = np.random.default_rng(seed)
    xf = np.linspace(X0 - W_OUT * SIGMA_X, X0 + W_OUT * SIGMA_X, 700001)
    pdf = amp(xf, convention) ** 2
    cdf = np.concatenate([[0.0],
                          np.cumsum(0.5 * (pdf[1:] + pdf[:-1]) * np.diff(xf))])
    cdf /= cdf[-1]
    xs = np.interp(rng.uniform(size=N_PART), cdf, xf)
    zfine = np.linspace(0.0, 1.0, 200001)
    zcdf = zfine - np.sin(2 * np.pi * zfine) / (2 * np.pi)
    zs = np.interp(rng.uniform(size=N_PART), zcdf, zfine)
    return np.column_stack([xs, zs])


# ----------------------------------------------------------------------
# Factorized exact-in-z engine (primary; adapted from validated l5_refine.py)
# ----------------------------------------------------------------------
def build_f(Nx, xmin, xmax, convention):
    dx = (xmax - xmin) / Nx
    xg = xmin + dx * np.arange(Nx)
    f = amp(xg, convention).astype(complex)
    f *= np.exp(1j * K0 * xg)
    f /= np.sqrt(np.sum(np.abs(f) ** 2) * dx)
    kx = 2.0 * np.pi * sfft.fftfreq(Nx, dx)
    return xg, dx, f, kx


def lin_interp(V, px, xg0, dx, Nx):
    g = np.clip((px - xg0) / dx, 0.0, Nx - 1.000001)
    i = g.astype(np.int64)
    fr = g - i
    return V[i] * (1 - fr) + V[i + 1] * fr


def run_1d(Nx, dt, P0, xmin=XMIN_E, xmax=XMAX_E, convention="literal",
           t_final=T_EVOLVE):
    xg, dx, f0, kx = build_f(Nx, xmin, xmax, convention)
    F = sfft.fft(f0)
    Eh = np.exp(-0.5j * kx ** 2 * (dt / 2.0))
    nsteps = int(round(t_final / dt))
    N = P0.shape[0]

    # 1-D diagnostics (exact propagation: norm/energy conserved by construction)
    spec_pow = np.abs(F) ** 2
    E1 = float(np.sum(0.5 * kx ** 2 * spec_pow) / np.sum(spec_pow))
    norm0 = float(np.sum(np.abs(f0) ** 2) * dx)

    def fields(F):
        f = sfft.ifft(F)
        fp = sfft.ifft(1j * kx * F)
        rho = np.abs(f) ** 2
        den = np.maximum(rho, RHO_EPS_REL * rho.max())
        w = np.conj(f) * fp
        a = np.clip(np.imag(w) / den, -VCLAMP, VCLAMP)
        b = np.clip(np.real(w) / den, -VCLAMP, VCLAMP)
        return a, b

    ens = {"T": P0.copy(), "A": P0.copy()}
    active = {"T": np.ones(N, bool), "A": np.ones(N, bool)}
    tau = {s: {d: np.full(N, np.nan) for d in ("near", "far")} for s in "TA"}
    dets = {"near": D_NEAR, "far": D_FAR}

    def vel(ab, s, px, pz):
        a, b = ab
        va = lin_interp(a, px, xg[0], dx, Nx)
        if s == "A":
            return va, np.zeros_like(va)
        vb = lin_interp(b, px, xg[0], dx, Nx)
        vx = np.clip(va - np.pi / np.tan(np.pi * pz / L), -VCLAMP, VCLAMP)
        return vx, vb

    ab = fields(F)
    t_start = time.time()
    for n in range(nsteps):
        t = n * dt
        F = F * Eh
        ab_h = fields(F)
        F = F * Eh
        ab_n = fields(F)
        for s in "TA":
            P = ens[s]
            act = active[s]
            if not act.any():
                continue
            idx = np.where(act)[0]
            px, pz = P[idx, 0], P[idx, 1]
            k1x, k1z = vel(ab, s, px, pz)
            k2x, k2z = vel(ab_h, s, px + 0.5 * dt * k1x,
                           np.clip(pz + 0.5 * dt * k1z, 1e-9, L - 1e-9))
            k3x, k3z = vel(ab_h, s, px + 0.5 * dt * k2x,
                           np.clip(pz + 0.5 * dt * k2z, 1e-9, L - 1e-9))
            k4x, k4z = vel(ab_n, s, px + dt * k3x,
                           np.clip(pz + dt * k3z, 1e-9, L - 1e-9))
            nx = px + dt / 6.0 * (k1x + 2 * k2x + 2 * k3x + k4x)
            nz = np.clip(pz + dt / 6.0 * (k1z + 2 * k2z + 2 * k3z + k4z),
                         1e-9, L - 1e-9)
            for dname, dpos in dets.items():
                ta = tau[s][dname]
                cross = np.isnan(ta[idx]) & (px < dpos) & (nx >= dpos)
                if cross.any():
                    ci = idx[cross]
                    frac = (dpos - px[cross]) / (nx[cross] - px[cross])
                    ta[ci] = t + dt * frac
            P[idx, 0] = nx
            P[idx, 1] = nz
            esc = (nx < xmin + 0.5) | (nx > xmax - 0.5)
            if esc.any():
                active[s][idx[esc]] = False
        ab = ab_n
    fT = sfft.ifft(F)
    norm_final = float(np.sum(np.abs(fT) ** 2) * dx)
    return dict(tau=tau, active=active, ens=ens, E1=E1, norm0=norm0,
                norm_final=norm_final, elapsed=time.time() - t_start)


# ----------------------------------------------------------------------
# Full 2-D FFT engine (cross-check; adapted from validated l5_arrival.py)
# ----------------------------------------------------------------------
def run_2d(Nx, Nz, dt, P0, xmin=XMIN_E, xmax=XMAX_E, convention="literal",
           t_final=T_EVOLVE, n_checks=17):
    dx = (xmax - xmin) / Nx
    xg = xmin + dx * np.arange(Nx)
    dz = L / Nz
    zg_ext = dz * np.arange(2 * Nz)

    fx = amp(xg, convention).astype(complex) * np.exp(1j * K0 * xg)
    gz = np.sin(np.pi * zg_ext / L)
    phi0 = fx[:, None] * gz[None, :]
    dA = dx * dz
    phi0 /= np.sqrt(np.sum(np.abs(phi0[:, :Nz]) ** 2) * dA)

    kx = 2.0 * np.pi * sfft.fftfreq(Nx, dx)
    kz = 2.0 * np.pi * sfft.fftfreq(2 * Nz, dz)
    KX, KZ = kx[:, None], kz[None, :]
    K2 = KX ** 2 + KZ ** 2

    F0 = sfft.fft2(phi0, workers=WORKERS)

    # ---- G1 diagnostics ----
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
    diag = dict(norm_deviation_max=float(np.max(np.abs(np.array(norms) - norm0))),
                wall_amplitude_max=float(np.max(walls)),
                energy_drift_rel=float(np.max(np.abs(np.array(energies) - E0))
                                       / abs(E0)),
                energy_mean=E0, z_mode_leakage=leak)

    def velocity_fields(F):
        stack = np.stack([F, 1j * KX * F, 1j * KZ * F])
        out = sfft.ifft2(stack, axes=(-2, -1), workers=WORKERS)
        phi, dpx, dpz = out[0], out[1], out[2]
        ph = phi[:, : Nz + 1]
        dx_ = dpx[:, : Nz + 1]
        dz_ = dpz[:, : Nz + 1]
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
        wall = max(np.abs(phi[:, 0]).max(), np.abs(phi[:, Nz]).max())
        return (vxT, vzT, vxA, vzA), wall

    def bilinear(V, px, pz):
        gx = np.clip((px - xmin) / dx, 0.0, Nx - 1.000001)
        gz = np.clip(pz / dz, 0.0, Nz - 1e-6)
        i = gx.astype(np.int64)
        j = gz.astype(np.int64)
        fx_ = gx - i
        fz_ = gz - j
        return (V[i, j] * (1 - fx_) * (1 - fz_) + V[i + 1, j] * fx_ * (1 - fz_)
                + V[i, j + 1] * (1 - fx_) * fz_ + V[i + 1, j + 1] * fx_ * fz_)

    nsteps = int(round(t_final / dt))
    Ehalf = np.exp(-0.5j * K2 * (dt / 2.0))
    F = F0.copy()
    N = P0.shape[0]
    ens = {"T": P0.copy(), "A": P0.copy()}
    active = {"T": np.ones(N, bool), "A": np.ones(N, bool)}
    tau = {s: {d: np.full(N, np.nan) for d in ("near", "far")} for s in "TA"}
    dets = {"near": D_NEAR, "far": D_FAR}

    vf, wall0 = velocity_fields(F)
    wall_max = wall0
    t_start = time.time()
    for n in range(nsteps):
        t = n * dt
        F = F * Ehalf
        vf_h, w1 = velocity_fields(F)
        F = F * Ehalf
        vf_n, w2 = velocity_fields(F)
        wall_max = max(wall_max, w1, w2)
        for s, (ivx, ivz) in (("T", (0, 1)), ("A", (2, 3))):
            P = ens[s]
            act = active[s]
            if not act.any():
                continue
            idx = np.where(act)[0]
            px, pz = P[idx, 0], P[idx, 1]
            k1x = bilinear(vf[ivx], px, pz)
            k1z = bilinear(vf[ivz], px, pz)
            ax, az = px + 0.5 * dt * k1x, pz + 0.5 * dt * k1z
            k2x = bilinear(vf_h[ivx], ax, az)
            k2z = bilinear(vf_h[ivz], ax, az)
            bx, bz = px + 0.5 * dt * k2x, pz + 0.5 * dt * k2z
            k3x = bilinear(vf_h[ivx], bx, bz)
            k3z = bilinear(vf_h[ivz], bx, bz)
            cx, cz = px + dt * k3x, pz + dt * k3z
            k4x = bilinear(vf_n[ivx], cx, cz)
            k4z = bilinear(vf_n[ivz], cx, cz)
            nx = px + dt / 6.0 * (k1x + 2 * k2x + 2 * k3x + k4x)
            nz = np.clip(pz + dt / 6.0 * (k1z + 2 * k2z + 2 * k3z + k4z),
                         1e-9, L - 1e-9)
            for dname, dpos in dets.items():
                ta = tau[s][dname]
                cross = np.isnan(ta[idx]) & (px < dpos) & (nx >= dpos)
                if cross.any():
                    ci = idx[cross]
                    frac = (dpos - px[cross]) / (nx[cross] - px[cross])
                    ta[ci] = t + dt * frac
            P[idx, 0] = nx
            P[idx, 1] = nz
            esc = (nx < xmin + 0.5) | (nx > xmax - 0.5)
            if esc.any():
                active[s][idx[esc]] = False
        vf = vf_n
    diag["wall_amplitude_max_traj_run"] = float(wall_max)
    return dict(tau=tau, active=active, ens=ens, diag=diag,
                elapsed=time.time() - t_start)


# ----------------------------------------------------------------------
# Summaries / gates
# ----------------------------------------------------------------------
def fin(a):
    return a[np.isfinite(a)]


def summarize(run, label):
    out = {"label": label}
    for s in "TA":
        for d in ("near", "far"):
            a = run["tau"][s][d]
            f = fin(a)
            out["n_%s_%s" % (s, d)] = int(f.size)
            out["max_%s_%s" % (s, d)] = float(f.max()) if f.size else None
    tT = run["tau"]["T"]["near"]
    w = fin(tT)
    w16 = w[w <= T_NEAR]
    out["T_near_max_win16"] = float(w16.max()) if w16.size else None
    out["T_near_n_in_16_24"] = int(np.sum((w > T_NEAR) & (w <= T_EVOLVE)))
    srt = np.sort(w16)
    out["T_near_top5_win16"] = [float(v) for v in srt[-5:]]
    return out


def main():
    P0 = sample_ensemble(SEED, "literal")
    P0_alt = sample_ensemble(SEED, "alt")

    ladder = [
        # label,                Nx,   dt,      box,  convention
        ("Nx2048_dt2e-3",       2048, 2.0e-3,  "E", "literal"),  # 2-D validation match
        ("Nx4096_dt1e-3",       4096, 1.0e-3,  "E", "literal"),
        ("Nx4096_dt5e-4",       4096, 5.0e-4,  "E", "literal"),
        ("Nx4096_dt2.5e-4",     4096, 2.5e-4,  "E", "literal"),  # REFERENCE
        ("Nx8192_dt2.5e-4",     8192, 2.5e-4,  "E", "literal"),  # refine: double Nx
        ("Nx4096_dt1.25e-4",    4096, 1.25e-4, "E", "literal"),  # refine: halve dt
        ("box2048_dt2.5e-4",    2048, 2.5e-4,  "P", "literal"),  # pre-registered box control
        ("alt_Nx4096_dt5e-4",   4096, 5.0e-4,  "E", "alt"),      # amplitude-convention check
    ]
    runs = {}
    for label, Nx, dt, box, conv in ladder:
        xmin, xmax = (XMIN_E, XMAX_E) if box == "E" else (XMIN_P, XMAX_P)
        p0 = P0_alt if conv == "alt" else P0
        print("run", label, flush=True)
        runs[label] = run_1d(Nx, dt, p0, xmin, xmax, conv)
        print("  elapsed %.1f s" % runs[label]["elapsed"],
              json.dumps(summarize(runs[label], label)), flush=True)

    print("run 2-D cross-check Nx2048 Nz128 dt2e-3 (extended box)", flush=True)
    run2d = run_2d(2048, 128, 2.0e-3, P0)
    print("  elapsed %.1f s" % run2d["elapsed"],
          json.dumps(summarize(run2d, "2d")), flush=True)
    print("  2-D diag:", json.dumps(run2d["diag"]), flush=True)

    # ---- cross-engine validation at matched (Nx, dt) ----
    t1d = runs["Nx2048_dt2e-3"]["tau"]["T"]["near"]
    t2d = run2d["tau"]["T"]["near"]
    both = np.isfinite(t1d) & np.isfinite(t2d)
    dd = np.abs(t1d[both] - t2d[both])
    validation = dict(n_common=int(both.sum()),
                      median_abs_diff=float(np.median(dd)),
                      p99_abs_diff=float(np.percentile(dd, 99)),
                      max_abs_diff=float(dd.max()),
                      ks_p=float(stats.ks_2samp(fin(t1d), fin(t2d)).pvalue))

    # ---- gate evaluation on the REFERENCE run ----
    ref = runs["Nx4096_dt2.5e-4"]
    refN = runs["Nx8192_dt2.5e-4"]
    refT = runs["Nx4096_dt1.25e-4"]

    tT_all = ref["tau"]["T"]["near"]
    tA_all = ref["tau"]["A"]["near"]
    tT = fin(tT_all)
    tA = fin(tA_all)
    tT16 = tT[tT <= T_NEAR]
    tA16 = tA[tA <= T_NEAR]
    fT = fin(ref["tau"]["T"]["far"])
    fA = fin(ref["tau"]["A"]["far"])

    # G1 (2-D engine diagnostics + 1-D norm/energy conservation)
    d2 = run2d["diag"]
    norm_dev_1d = abs(ref["norm_final"] - ref["norm0"])
    g1_pass = (d2["norm_deviation_max"] < 1e-6
               and d2["wall_amplitude_max"] < 1e-10
               and d2["energy_drift_rel"] <= 1e-6)

    # G2: transverse near-field hard cutoff (full ensemble, no exclusions)
    tau_max = float(tT16.max())
    n_T_gap = int(np.sum((tT > tau_max) & (tT <= T_NEAR)))       # 0 by construction
    n_T_late = int(np.sum(tT > T_NEAR))                          # any in (16, 24]
    iqr_T = float(np.subtract(*np.percentile(tT16, [75, 25])))
    gap = T_NEAR - tau_max
    tau_max_N = float(max(fin(refN["tau"]["T"]["near"])[
        fin(refN["tau"]["T"]["near"]) <= T_NEAR]))
    tau_max_T = float(max(fin(refT["tau"]["T"]["near"])[
        fin(refT["tau"]["T"]["near"]) <= T_NEAR]))
    rel_N = abs(tau_max_N - tau_max) / tau_max
    rel_T = abs(tau_max_T - tau_max) / tau_max

    # sample-maximum trajectory identity + convergence across the ladder
    imax = int(np.nanargmax(np.where(tT_all <= T_NEAR, tT_all, np.nan)))
    max_traj = dict(index=imax, x0=float(P0[imax, 0]), z0=float(P0[imax, 1]),
                    tau_by_run={lab: (float(runs[lab]["tau"]["T"]["near"][imax])
                                      if np.isfinite(runs[lab]["tau"]["T"]["near"][imax])
                                      else None)
                                for lab, _, _, _, conv in ladder if conv == "literal"},
                    tau_2d=(float(run2d["tau"]["T"]["near"][imax])
                            if np.isfinite(run2d["tau"]["T"]["near"][imax]) else None))
    vals = [v for k, v in max_traj["tau_by_run"].items()
            if v is not None and k != "box2048_dt2.5e-4"]
    max_traj["spread_rel"] = float((max(vals) - min(vals)) / np.mean(vals))

    # full-ensemble per-trajectory convergence (every edge trajectory)
    def pairdiff(a, b):
        both = np.isfinite(a) & np.isfinite(b)
        d = np.abs(a[both] - b[both])
        return dict(n_common=int(both.sum()),
                    n_status_changed=int(np.sum(np.isfinite(a) != np.isfinite(b))),
                    median=float(np.median(d)), p99=float(np.percentile(d, 99)),
                    max=float(d.max()))
    conv_full = dict(
        ref_vs_doubleNx=pairdiff(tT_all, refN["tau"]["T"]["near"]),
        ref_vs_halvedt=pairdiff(tT_all, refT["tau"]["T"]["near"]))

    g2_pass = (n_T_gap == 0 and n_T_late == 0 and gap >= 2.0 * iqr_T
               and rel_N <= 0.02 and rel_T <= 0.02)

    # G3: axial contrast
    n_axial_late = int(np.sum(tA16 > tau_max))
    frac_axial_late = n_axial_late / max(1, tA16.size)
    g3_pass = frac_axial_late > 0.05

    # G4: far-field KS
    ks = stats.ks_2samp(fT, fA)
    ks2d = stats.ks_2samp(fin(run2d["tau"]["T"]["far"]),
                          fin(run2d["tau"]["A"]["far"]))
    g4_pass = ks.pvalue > 0.05

    # G5: mean-arrival ratio with bootstrap (near field, analysis window)
    rng_b = np.random.default_rng(777)
    ratio = float(tT16.mean() / tA16.mean())
    boots = np.empty(2000)
    for b in range(2000):
        boots[b] = (rng_b.choice(tT16, tT16.size).mean()
                    / rng_b.choice(tA16, tA16.size).mean())
    ratio_err = float(boots.std(ddof=1))

    # pre-registered-box control and alt-convention sensitivity
    box = runs["box2048_dt2.5e-4"]
    tTb = fin(box["tau"]["T"]["near"]); tTb16 = tTb[tTb <= T_NEAR]
    tAb = fin(box["tau"]["A"]["near"]); tAb16 = tAb[tAb <= T_NEAR]
    box_ctl = dict(tau_max=float(tTb16.max()),
                   tau_max_rel_vs_ref=float(abs(tTb16.max() - tau_max) / tau_max),
                   frac_axial_late=float(np.mean(tAb16 > tTb16.max())),
                   ks_far_p=float(stats.ks_2samp(fin(box["tau"]["T"]["far"]),
                                                 fin(box["tau"]["A"]["far"])).pvalue),
                   n_T_far=int(fin(box["tau"]["T"]["far"]).size),
                   n_A_far=int(fin(box["tau"]["A"]["far"]).size))
    alt = runs["alt_Nx4096_dt5e-4"]
    tTa = fin(alt["tau"]["T"]["near"]); tTa16 = tTa[tTa <= T_NEAR]
    tAa = fin(alt["tau"]["A"]["near"]); tAa16 = tAa[tAa <= T_NEAR]
    alt_ctl = dict(tau_max=float(tTa16.max()),
                   frac_axial_late=float(np.mean(tAa16 > tTa16.max())),
                   n_axial_late=int(np.sum(tAa16 > tTa16.max())),
                   iqr=float(np.subtract(*np.percentile(tTa16, [75, 25]))),
                   ks_far_p=float(stats.ks_2samp(fin(alt["tau"]["T"]["far"]),
                                                 fin(alt["tau"]["A"]["far"])).pvalue))

    # energy convergence across Nx (evidence the smoothing removed the UV issue)
    E_by_Nx = {}
    for Nx in (2048, 4096, 8192):
        _, _, f0, kx = build_f(Nx, XMIN_E, XMAX_E, "literal")
        Fq = sfft.fft(f0)
        p = np.abs(Fq) ** 2
        E_by_Nx[str(Nx)] = float(np.sum(0.5 * kx ** 2 * p) / np.sum(p))

    results = {
        "workstream": "L5prime",
        "file_id": "F-T6-L5-EXEC-2",
        "params": dict(L=L, sigma_x=SIGMA_X, x0=X0, k0=K0,
                       taper="S2 smoothstep, W=1 for |x-x0|<=2.5 sigma, "
                             "W=0 at 3.5 sigma",
                       amplitude_convention="literal: exp(-(x-x0)^2/(2 sigma^2)) * W",
                       box_preregistered=[XMIN_P, XMAX_P],
                       box_extended=[XMIN_E, XMAX_E],
                       Nx_ref=4096, dx_ref=(XMAX_E - XMIN_E) / 4096,
                       dt_ref=2.5e-4, T_evolve=T_EVOLVE, T_near_window=T_NEAR,
                       N=N_PART, seed=SEED, d_near=D_NEAR, d_far=D_FAR,
                       vclamp=VCLAMP, rho_eps_rel=RHO_EPS_REL),
        "G1": dict(norm_deviation_max=d2["norm_deviation_max"],
                   wall_amplitude_max=d2["wall_amplitude_max"],
                   wall_amplitude_max_traj_run=d2["wall_amplitude_max_traj_run"],
                   energy_drift_rel=d2["energy_drift_rel"],
                   energy_mean_2d=d2["energy_mean"],
                   z_mode_leakage=d2["z_mode_leakage"],
                   norm_dev_1d_ref=float(norm_dev_1d),
                   energy_1d_ref=float(ref["E1"]),
                   energy_1d_by_Nx=E_by_Nx,
                   verdict="PASS" if g1_pass else "FAIL"),
        "G2": dict(tau_max=tau_max,
                   arrivals_in_gap=n_T_gap,
                   transverse_near_arrivals_in_16_24=n_T_late,
                   empty_gap=gap, iqr_transverse_near=iqr_T,
                   required_gap=2.0 * iqr_T,
                   tau_max_doubleNx=tau_max_N, rel_change_doubleNx=rel_N,
                   tau_max_halvedt=tau_max_T, rel_change_halvedt=rel_T,
                   max_trajectory=max_traj,
                   full_ensemble_convergence=conv_full,
                   verdict="PASS" if g2_pass else "FAIL"),
        "G3": dict(frac_axial_beyond_tau_max=frac_axial_late,
                   n_axial_beyond=n_axial_late,
                   n_axial_arrivals_win16=int(tA16.size),
                   axial_last_arrival_win16=float(tA16.max()),
                   axial_last_arrival_T24=float(tA.max()),
                   frac_axial_beyond_tau_max_T24=float(np.mean(tA > tau_max)),
                   verdict="PASS" if g3_pass else "FAIL"),
        "G4": dict(ks_stat=float(ks.statistic), ks_p=float(ks.pvalue),
                   n_far_T=int(fT.size), n_far_A=int(fA.size),
                   ks_p_2d_engine=float(ks2d.pvalue),
                   verdict="PASS" if g4_pass else "FAIL"),
        "G5": dict(mean_tau_T_near=float(tT16.mean()),
                   mean_tau_A_near=float(tA16.mean()),
                   ratio_T_over_A=ratio, ratio_bootstrap_err=ratio_err,
                   verdict="REPORT-ONLY"),
        "counts": dict(
            near_T_win16=int(tT16.size), near_A_win16=int(tA16.size),
            near_T_T24=int(tT.size), near_A_T24=int(tA.size),
            far_T=int(fT.size), far_A=int(fA.size),
            nonarrival_near_T_win16=int(N_PART - tT16.size),
            nonarrival_near_A_win16=int(N_PART - tA16.size),
            nonarrival_near_T_T24=int(N_PART - tT.size),
            nonarrival_near_A_T24=int(N_PART - tA.size),
            nonarrival_far_T=int(N_PART - fT.size),
            nonarrival_far_A=int(N_PART - fA.size),
            frozen_escapees_T=int(np.sum(~ref["active"]["T"])),
            frozen_escapees_A=int(np.sum(~ref["active"]["A"]))),
        "quartiles": dict(
            T_near=[float(q) for q in np.percentile(tT16, [25, 50, 75])],
            A_near=[float(q) for q in np.percentile(tA16, [25, 50, 75])],
            T_far=[float(q) for q in np.percentile(fT, [25, 50, 75])],
            A_far=[float(q) for q in np.percentile(fA, [25, 50, 75])]),
        "cross_engine_validation": validation,
        "ladder": {lab: summarize(runs[lab], lab)
                   for lab, _, _, _, _ in ladder},
        "run_2d": summarize(run2d, "2d_Nx2048_Nz128_dt2e-3"),
        "preregistered_box_control": box_ctl,
        "alt_amplitude_convention": alt_ctl,
        "L5_comparison": dict(
            L5=dict(tau_max=2.28503, gap=5.715, iqr=0.3819,
                    frac_axial_beyond=0.028, ks_far_p=0.307,
                    edge_trajectory="idx 1249 non-convergent (7 resolutions)"),
        ),
    }
    results["L5_comparison"]["L5prime"] = dict(
        tau_max=tau_max, gap=gap, iqr=iqr_T,
        frac_axial_beyond=frac_axial_late, ks_far_p=float(ks.pvalue),
        edge_trajectory="idx %d spread_rel %.2e across ladder"
                        % (imax, max_traj["spread_rel"]))

    with open(os.path.join(OUTDIR, "L5prime_results.json"), "w") as fh:
        json.dump(results, fh, indent=2)

    np.savez_compressed(
        os.path.join(OUTDIR, "L5prime_raw.npz"),
        P0=P0,
        tau_T_near=ref["tau"]["T"]["near"], tau_A_near=ref["tau"]["A"]["near"],
        tau_T_far=ref["tau"]["T"]["far"], tau_A_far=ref["tau"]["A"]["far"],
        tau_T_near_doubleNx=refN["tau"]["T"]["near"],
        tau_T_near_halvedt=refT["tau"]["T"]["near"],
        tau_T_near_2d=run2d["tau"]["T"]["near"],
        tau_T_far_2d=run2d["tau"]["T"]["far"],
        tau_A_far_2d=run2d["tau"]["A"]["far"],
        final_T=ref["ens"]["T"], final_A=ref["ens"]["A"],
        **{"ladder_%s" % lab: runs[lab]["tau"]["T"]["near"]
           for lab, _, _, _, conv in ladder if conv == "literal"})

    print(json.dumps({k: results[k] for k in
                      ("G1", "G2", "G3", "G4", "G5", "counts")}, indent=2))
    print("verdicts:", {g: results[g]["verdict"] for g in
                        ("G1", "G2", "G3", "G4", "G5")})


if __name__ == "__main__":
    main()
