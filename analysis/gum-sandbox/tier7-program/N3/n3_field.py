#!/usr/bin/env python3
"""
WORKSTREAM N3 -- L7 continuum proof attempt of the cutoff zeros via the
field-level crossing criterion (T4-W5 residual).  File id: F-T7-N3.

Theorem (crossing criterion).  If a trajectory of dX/dt = v_x, dZ/dt = v_z
first crosses the line x = d at time t0 (X(t0) = d, X(t) < d for t < t0),
then v_x(d, Z(t0), t0) >= 0: the left difference quotient
(X(t0) - X(t0 - h))/h is >= 0 for all small h > 0, so its limit dX/dt(t0)
= v_x at the crossing point is >= 0.  Contrapositive: if
sup_z v_x(d, z, t) < 0 for ALL t in (tau*, T], no trajectory first-crosses
x = d in (tau*, T]; the continuum first-crossing density vanishes there and
its support within [0, T] is contained in [0, tau*].

This workstream MEASURES the hypothesis on the validated L5prime/L7b engine
(exact-in-z factorized evolution; velocity affine in n_y at field level:
v_x(d, z, t; n_y) = A(t) + n_y * S(z),  A(t) = Im(f'/f)(d, t) clamped as in
the engine, S(z) = -pi cot(pi z); the x-part A is z-independent and the
spin part S is t-independent -- an exact structural property of the
factorized state phi = f(x,t) sin(pi z) e^{-i pi^2 t/2}).

Configuration inherited unchanged: k0 = 2, C2 smoothstep taper, literal
amplitude exp(-(x-x0)^2/(2 sigma^2)), box x in [-15, 45], Nx = 2048
reference, dt = 2.5e-4 (engine-accumulation replication), T = 16,
d = d_near = 1, hbar = m = 1.  Output sampling dt_out = 0.01 (spec
<= 0.01).  Refinements: Nx = 4096 (doubling) and dt = 1.25e-4 (halving).

Also included, as DIAGNOSIS (not a pre-registered gate): a deterministic
covering of the initial support with a 150 x 100 grid of trajectories for
the kill set n_y in {-1, -0.75, +0.75, +1}, to test the L7b kill-bin zeros
at field-like resolution independently of the |phi0|^2 sample.
"""
import json
import os
import time

import numpy as np
from scipy import fft as sfft

# ----------------------------------------------------------------------
# Parameters (inherited verbatim from L7b / L5prime)
# ----------------------------------------------------------------------
L = 1.0
SIGMA_X = 1.0
X0 = -5.0
K0 = 2.0
W_IN, W_OUT = 2.5, 3.5
D_NEAR = 1.0
T_FINAL = 16.0
VCLAMP = 500.0
RHO_EPS_REL = 1e-14
XMIN, XMAX = -15.0, 45.0
NX_REF = 2048
DT_REF = 2.5e-4
DT_OUT = 0.01

NY_GATE = [1.0, 0.75, -0.75, -1.0, 0.25, 0.0]   # gate set + control + axial
NY_COVER = np.array([-1.0, -0.75, 0.75, 1.0])   # kill set, covering diagnostic

NZ_LINE = 1999                                   # z-grid on the line (open)
FLOORS = [1e-1, 1e-2, 1e-3]                      # rho floors: sin^2(pi z) >= eps

TAU_MAX_SAMPLED_P1 = 5.130950                    # L5prime G2 / spec value
TAU_MAX_L7B = {"+1.00": 5.1309462, "+0.75": 5.9993, "-0.75": 6.0263,
               "-1.00": 5.1289, "+0.25": 15.627, "+0.00": 15.836}
KILL_BINS = [(6.4, 6.8), (6.8, 7.2)]

OUTDIR = os.path.dirname(os.path.abspath(__file__))
T6 = os.path.join(os.path.dirname(os.path.dirname(OUTDIR)), "tier6-foundations")
L5P_RAW = os.path.join(T6, "L5prime", "L5prime_raw.npz")
L7B_RAW = os.path.join(T6, "L7b", "l7b_raw.npz")


# ----------------------------------------------------------------------
# Initial profile + engine pieces (verbatim from l7b_run.py)
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


def ab_fields(F, kx, den_floor=RHO_EPS_REL):
    f = sfft.ifft(F)
    fp = sfft.ifft(1j * kx * F)
    rho = np.abs(f) ** 2
    den = np.maximum(rho, den_floor * rho.max())
    w = np.conj(f) * fp
    a = np.clip(np.imag(w) / den, -VCLAMP, VCLAMP)
    b = np.clip(np.real(w) / den, -VCLAMP, VCLAMP)
    return a, b, rho


# ----------------------------------------------------------------------
# Line field A(t) = a(x=d, t) (engine's clamped/interpolated v_x at z=1/2)
# ----------------------------------------------------------------------
def line_field_exact(Nx, tgrid):
    """Exact spectral evaluation at each output time (the propagator is
    exact; no time-stepping error exists for the field itself)."""
    xg, dx, f0, kx = build_f(Nx, XMIN, XMAX)
    F0 = sfft.fft(f0)
    A = np.empty(tgrid.size)
    RHO = np.empty(tgrid.size)
    norm_dev = 0.0
    norm0 = float(np.sum(np.abs(f0) ** 2) * dx)
    p0 = np.abs(F0) ** 2
    E0 = float(np.sum(0.5 * kx ** 2 * p0) / np.sum(p0))
    e_drift = 0.0
    pt = np.array([D_NEAR])
    for j, t in enumerate(tgrid):
        F = F0 * np.exp(-0.5j * kx ** 2 * t)
        a, b, rho = ab_fields(F, kx)
        A[j] = lin_interp(a, pt, xg[0], dx, Nx)[0]
        RHO[j] = lin_interp(rho, pt, xg[0], dx, Nx)[0]
        norm_dev = max(norm_dev, abs(float(np.sum(rho) * dx) - norm0))
        p = np.abs(F) ** 2
        e_drift = max(e_drift, abs(float(np.sum(0.5 * kx ** 2 * p)
                                         / np.sum(p)) - E0) / abs(E0))
    return A, RHO, dict(norm_dev=norm_dev, energy_drift_rel=e_drift, E0=E0)


def line_field_stepped(Nx, dt, tgrid):
    """Replicates the engine's cumulative half-step accumulation of the
    propagator (the only dt-dependence the field has: float rounding)."""
    xg, dx, f0, kx = build_f(Nx, XMIN, XMAX)
    F = sfft.fft(f0)
    Eh = np.exp(-0.5j * kx ** 2 * (dt / 2.0))
    nsteps = int(round(T_FINAL / dt))
    per_out = int(round(DT_OUT / dt))
    A = np.empty(tgrid.size)
    pt = np.array([D_NEAR])
    a, _, _ = ab_fields(F, kx)
    A[0] = lin_interp(a, pt, xg[0], dx, Nx)[0]
    jout = 1
    for n in range(nsteps):
        F = F * Eh
        F = F * Eh
        if (n + 1) % per_out == 0 and jout < tgrid.size:
            a, _, _ = ab_fields(F, kx)
            A[jout] = lin_interp(a, pt, xg[0], dx, Nx)[0]
            jout += 1
    return A


# ----------------------------------------------------------------------
# 2-D wall-amplitude diagnostic (verbatim from l7b_run.py)
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
    F0 = sfft.fft2(phi0)
    norm0 = float(np.sum(np.abs(phi0[:, :Nz]) ** 2) * dA)
    Pspec = np.sum(np.abs(F0) ** 2)
    E0 = float(np.sum(0.5 * K2 * np.abs(F0) ** 2) / Pspec)
    rowpow = np.sum(np.abs(F0) ** 2, axis=0)
    mode1 = float(rowpow[1] + rowpow[2 * Nz - 1])
    leak = float((rowpow.sum() - mode1) / rowpow.sum())
    norms, walls, energies = [], [], []
    for tchk in np.linspace(0.0, t_final, n_checks):
        Ft = F0 * np.exp(-0.5j * K2 * tchk)
        phit = sfft.ifft2(Ft)
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
# Criterion evaluation:  v_x(d, z, t; n_y) = A(t) + n_y S(z), clamped
# ----------------------------------------------------------------------
def criterion(A, tgrid):
    zg = np.arange(1, NZ_LINE + 1) / (NZ_LINE + 1.0)
    S = -np.pi / np.tan(np.pi * zg / L)          # t-independent spin part
    sin2 = np.sin(np.pi * zg) ** 2               # rho(d,z,t) = 2|f|^2 sin^2
    out = {}
    heat_p1 = None
    for ny in NY_GATE:
        V = np.clip(A[:, None] + ny * S[None, :], -VCLAMP, VCLAMP)
        if ny == 1.0:
            heat_p1 = V.astype(np.float32)
        regions = {"full": np.ones(zg.size, bool)}
        for eps in FLOORS:
            regions["floor%g" % eps] = sin2 >= eps
        entry = {}
        for rname, mask in regions.items():
            M = V[:, mask].max(axis=1)
            entry[rname] = tau_star_pack(M, tgrid)
        # mid-channel restriction (most favorable single point: S(1/2)=0)
        entry["mid"] = tau_star_pack(A.copy(), tgrid)
        out["%+.2f" % ny] = entry
    return out, heat_p1, zg


def tau_star_pack(M, tgrid):
    nonneg = M >= 0.0
    if nonneg.any():
        tau = float(tgrid[nonneg][-1])
    else:
        tau = 0.0
    exists = bool(tau < T_FINAL - 1e-12)
    pack = dict(tau_star=tau, exists_below_T=exists,
                frac_t_nonneg=float(nonneg.mean()),
                M_at_T=float(M[-1]),
                M_max=float(M.max()), M_min=float(M.min()))
    # margin on (tau*+0.1, T]
    w = tgrid > tau + 0.1
    if exists and w.any():
        delta = -M[w]
        pack["delta_min_after_tau01"] = float(delta.min())
        pack["delta_median_after_tau01"] = float(np.median(delta))
    else:
        pack["delta_min_after_tau01"] = None
    # positivity through the kill bins and at sampled tau_max
    for (b0, b1) in KILL_BINS:
        m = (tgrid > b0) & (tgrid <= b1)
        pack["M_min_bin_%g_%g" % (b0, b1)] = float(M[m].min())
    return pack


# ----------------------------------------------------------------------
# Deterministic covering diagnostic (kill n_y set, grid initial data)
# ----------------------------------------------------------------------
def covering_run(Nx, dt, nx_grid=150, nz_grid=100, label="cover"):
    xg, dx, f0, kx = build_f(Nx, XMIN, XMAX)
    F = sfft.fft(f0)
    Eh = np.exp(-0.5j * kx ** 2 * (dt / 2.0))
    nsteps = int(round(T_FINAL / dt))

    # cell-centre covering of the compact support x in [x0-3.5s, x0+3.5s]
    xs = (X0 - W_OUT * SIGMA_X) + (np.arange(nx_grid) + 0.5) \
        * (2 * W_OUT * SIGMA_X / nx_grid)
    zs = (np.arange(nz_grid) + 0.5) / nz_grid
    XX, ZZ = np.meshgrid(xs, zs, indexing="ij")
    P0 = np.column_stack([XX.ravel(), ZZ.ravel()])
    wgt = (amp(P0[:, 0]) ** 2) * 2.0 * np.sin(np.pi * P0[:, 1]) ** 2
    wgt = wgt / wgt.sum()
    npts = P0.shape[0]

    NFAM = NY_COVER.size
    ny_flat = np.repeat(NY_COVER, npts)
    P = np.tile(P0, (NFAM, 1))
    active = np.ones(NFAM * npts, bool)
    tau = np.full(NFAM * npts, np.nan)

    def fields(F):
        a, b, _ = ab_fields(F, kx)
        return a, b

    def vel(ab, px, pz, ny):
        a, b = ab
        va = lin_interp(a, px, xg[0], dx, Nx)
        vb = lin_interp(b, px, xg[0], dx, Nx)
        vx = np.clip(va - ny * np.pi / np.tan(np.pi * pz / L),
                     -VCLAMP, VCLAMP)
        return vx, ny * vb

    ab = fields(F)
    t_start = time.time()
    for n in range(nsteps):
        t = n * dt
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
            esc = (nx < XMIN + 0.5) | (nx > XMAX - 0.5)
            if esc.any():
                active[idx[esc]] = False
        ab = ab_n
        if n % 8000 == 0:
            print("  [%s] step %d/%d t=%.2f active=%d elapsed=%.0fs"
                  % (label, n, nsteps, t, active.sum(),
                     time.time() - t_start), flush=True)
    tau = tau.reshape(NFAM, npts)
    active = active.reshape(NFAM, npts)
    summary = {}
    for i, ny in enumerate(NY_COVER):
        t_i = tau[i]
        crossed = np.isfinite(t_i)
        tc = t_i[crossed]
        wc = wgt[crossed]
        smax = TAU_MAX_SAMPLED_P1 if abs(ny) == 1.0 else \
            max(TAU_MAX_L7B["+0.75"], TAU_MAX_L7B["-0.75"])
        late = crossed & (t_i > smax + 0.01)
        kb = crossed & (t_i > KILL_BINS[0][0]) & (t_i <= KILL_BINS[1][1])
        nonc = ~crossed
        summary["%+.2f" % ny] = dict(
            n_grid=int(npts), n_crossed=int(crossed.sum()),
            tau_max_grid=float(tc.max()),
            weight_crossed=float(wc.sum()),
            n_late_beyond_sampled=int(late.sum()),
            weight_late_beyond_sampled=float(wgt[late].sum()),
            late_list=[dict(x0=float(P0[j, 0]), z0=float(P0[j, 1]),
                            tau=float(t_i[j]), weight=float(wgt[j]))
                       for j in np.where(late)[0][:20]],
            n_in_kill_bins=int(kb.sum()),
            weight_in_kill_bins=float(wgt[kb].sum()),
            n_noncross=int(nonc.sum()),
            weight_noncross=float(wgt[nonc].sum()),
            noncross_z_range=[float(P0[nonc, 1].min()),
                              float(P0[nonc, 1].max())] if nonc.any() else None,
            noncross_frozen=int((~active[i] & nonc).sum()),
            tau_p999_weighted=float(weighted_quantile(tc, wc, 0.999))
            if tc.size else None,
        )
    return dict(summary=summary, elapsed=time.time() - t_start), tau, P0, wgt


def weighted_quantile(v, w, q):
    o = np.argsort(v)
    cw = np.cumsum(w[o]) / w.sum()
    return v[o][np.searchsorted(cw, q)]


# ----------------------------------------------------------------------
def main():
    t0 = time.time()
    tgrid = np.round(np.arange(0.0, T_FINAL + 1e-9, DT_OUT), 10)

    print("line field, exact, Nx=2048 ...", flush=True)
    A_ref, RHO_ref, diag1 = line_field_exact(NX_REF, tgrid)
    print("line field, exact, Nx=4096 ...", flush=True)
    A_2x, RHO_2x, diag2 = line_field_exact(2 * NX_REF, tgrid)
    print("line field, stepped, dt=2.5e-4 ...", flush=True)
    A_st = line_field_stepped(NX_REF, DT_REF, tgrid)
    print("line field, stepped, dt=1.25e-4 ...", flush=True)
    A_st2 = line_field_stepped(NX_REF, DT_REF / 2.0, tgrid)

    b_Nx = np.abs(A_ref - A_2x)
    b_dt = np.abs(A_st - A_st2)
    b_mode = np.abs(A_ref - A_st)
    B_tot = b_Nx + b_dt

    def bstats(b, mask):
        v = b[mask]
        return dict(median=float(np.median(v)), p95=float(np.percentile(v, 95)),
                    max=float(v.max()))

    m_after = tgrid > TAU_MAX_SAMPLED_P1
    m_kill = (tgrid > KILL_BINS[0][0]) & (tgrid <= KILL_BINS[1][1])
    bounds = dict(
        Nx_doubling=dict(after_sampled_taumax=bstats(b_Nx, m_after),
                         kill_bins=bstats(b_Nx, m_kill),
                         full=bstats(b_Nx, np.ones_like(m_after))),
        dt_halving=dict(after_sampled_taumax=bstats(b_dt, m_after),
                        kill_bins=bstats(b_dt, m_kill),
                        full=bstats(b_dt, np.ones_like(m_after))),
        exact_vs_stepped=dict(after_sampled_taumax=bstats(b_mode, m_after),
                              kill_bins=bstats(b_mode, m_kill)),
        total_pointwise=dict(after_sampled_taumax=bstats(B_tot, m_after),
                             kill_bins=bstats(B_tot, m_kill)),
    )

    print("2-D wall diagnostic ...", flush=True)
    d2 = diag_2d()

    print("criterion evaluation ...", flush=True)
    crit_ref, heat_p1, zg_line = criterion(A_ref, tgrid)
    crit_2x, _, _ = criterion(A_2x, tgrid)
    crit_dt2, _, _ = criterion(A_st2, tgrid)

    # tau* stability (mid restriction, the only non-degenerate variant)
    stab = {}
    for key in crit_ref:
        e = {}
        for rname in ("full", "floor0.01", "mid"):
            tr = crit_ref[key][rname]["tau_star"]
            tN = crit_2x[key][rname]["tau_star"]
            tT = crit_dt2[key][rname]["tau_star"]
            e[rname] = dict(ref=tr, doubleNx=tN, halvedt=tT,
                            rel_Nx=float(abs(tN - tr) / max(tr, 1e-12)),
                            rel_dt=float(abs(tT - tr) / max(tr, 1e-12)))
        stab[key] = e

    # consistency vs sampled ensembles
    l5p = np.load(L5P_RAW)
    l7b = np.load(L7B_RAW)
    tau_A = l5p["tau_A_near"]
    tau_A = tau_A[np.isfinite(tau_A) & (tau_A <= T_FINAL)]
    ny_vals = l7b["ny_values"]
    i_p1 = int(np.where(ny_vals == 1.0)[0][0])
    tau_p1 = l7b["tau_ref"][i_p1]
    tau_p1 = tau_p1[np.isfinite(tau_p1) & (tau_p1 <= T_FINAL)]
    A_at_axial = np.interp(tau_A, tgrid, A_ref)
    late_ax = tau_A[tau_A > TAU_MAX_SAMPLED_P1]
    A_at_late_ax = np.interp(late_ax, tgrid, A_ref)
    consistency = dict(
        sampled_tau_max_p1_l7b=float(tau_p1.max()),
        sampled_tau_max_p1_spec=TAU_MAX_SAMPLED_P1,
        sampled_axial_last_arrival=float(tau_A.max()),
        n_axial_beyond_taumax=int(late_ax.size),
        A_at_axial_arrivals=dict(min=float(A_at_axial.min()),
                                 frac_ge_m005=float((A_at_axial >= -0.05).mean())),
        A_at_late_axial=dict(min=float(A_at_late_ax.min()),
                             median=float(np.median(A_at_late_ax)),
                             frac_ge_0=float((A_at_late_ax >= 0).mean())),
        tau_star_mid_ref=crit_ref["+0.00"]["mid"]["tau_star"],
        axial_support_inside_tau_star_mid=bool(
            float(tau_A.max()) <= crit_ref["+0.00"]["mid"]["tau_star"] + DT_OUT),
    )

    print("covering diagnostic (60k trajectories) ...", flush=True)
    cover, tau_cover, P0_cover, wgt_cover = covering_run(NX_REF, DT_REF)

    # ------------------------------------------------------------------
    # Gates
    # ------------------------------------------------------------------
    g1 = dict(diag_2d=d2, line_norm_dev=diag1["norm_dev"],
              line_energy_drift_rel=diag1["energy_drift_rel"],
              energy_1d=diag1["E0"])
    g1["verdict"] = "PASS" if (d2["norm_deviation_max"] < 1e-6
                               and d2["wall_amplitude_max"] < 1e-10
                               and d2["energy_drift_rel"] <= 1e-6
                               and diag1["norm_dev"] < 1e-6
                               and diag1["energy_drift_rel"] <= 1e-6) else "FAIL"

    p1 = crit_ref["+1.00"]
    g2 = dict(tau_star_full=p1["full"]["tau_star"],
              exists_below_T_full=p1["full"]["exists_below_T"],
              tau_star_floors={r: p1[r]["tau_star"]
                               for r in p1 if r.startswith("floor")},
              tau_star_mid=p1["mid"]["tau_star"],
              exists_below_T_mid=p1["mid"]["exists_below_T"],
              stability=stab["+1.00"],
              frac_t_nonneg_full=p1["full"]["frac_t_nonneg"],
              frac_t_nonneg_mid=p1["mid"]["frac_t_nonneg"],
              note=("M(t) = max_z v_x(d,z,t) >= 0 at EVERY output time on "
                    "the full line and on every rho-floored region: the "
                    "spin term -n_y pi cot(pi z) has a positive branch at "
                    "one wall (max over a floored region = A(t) + "
                    "|n_y| pi cot(pi z_floor)), and even the most "
                    "favorable single-point restriction z = 1/2 (where the "
                    "spin term vanishes and v_x = A(t) for every n_y) "
                    "stays nonnegative until t = %.2f, far beyond the "
                    "sampled cutoff 5.131.  No tau* < T exists in any "
                    "variant that the theorem's hypothesis accepts."
                    % p1["mid"]["tau_star"]))
    g2["verdict"] = "FAIL"

    g3 = dict(tau_star_used=p1["full"]["tau_star"],
              sampled_tau_max=TAU_MAX_SAMPLED_P1,
              consistent=bool(p1["full"]["tau_star"] >= TAU_MAX_SAMPLED_P1),
              gap=float(p1["full"]["tau_star"] - TAU_MAX_SAMPLED_P1),
              gap_mid=float(p1["mid"]["tau_star"] - TAU_MAX_SAMPLED_P1),
              explanation=("No contradiction (tau* >= sampled tau_max in "
                           "every variant), but the gap ~10.7-10.9 is NOT "
                           "the sampled-max-below-continuum-edge gap the "
                           "gate anticipated: tau* is pinned near T by the "
                           "mid-channel identity v_x(d,1/2,t) = A(t) (the "
                           "axial line velocity, positive through the kill "
                           "bins and recurring to 15.87) and by the "
                           "positive wall branch of the spin term.  tau* "
                           "is therefore not an estimate of the continuum "
                           "support edge; the criterion is vacuous above "
                           "the sampled cutoff."))
    g3["verdict"] = "PARTIAL"

    g4 = dict(window=None,
              delta_min=None,
              disc_bound_after_sampled=bounds["total_pointwise"][
                  "after_sampled_taumax"],
              disc_bound_kill_bins=bounds["total_pointwise"]["kill_bins"],
              A_min_kill_bins=float(A_ref[m_kill].min()),
              note=("The margin window (tau*+0.1, T] is empty (tau* = T on "
                    "the full line and on all floored regions).  Where the "
                    "gate wanted delta(t) = -M(t) > 0, the measurement "
                    "gives M(t) > 0: in the kill bins M >= A(t) >= %.3f "
                    "(mid-channel alone), i.e. the WRONG SIGN by ~1e4 x "
                    "the discretization bound (%.1e median).  The "
                    "hypothesis of the theorem is not merely unverified -- "
                    "it is refuted at proof-grade-modulo-discretization."
                    % (float(A_ref[m_kill].min()),
                       bounds["total_pointwise"]["kill_bins"]["median"])))
    g4["verdict"] = "FAIL"

    p025 = crit_ref["+0.25"]
    g5 = dict(tau_star_full_p025=p025["full"]["tau_star"],
              exists_below_T_p025=p025["full"]["exists_below_T"],
              tau_star_mid_p025=p025["mid"]["tau_star"],
              same_for_p1=dict(tau_star_full=p1["full"]["tau_star"],
                               tau_star_mid=p1["mid"]["tau_star"]),
              note=("Literal clause holds: n_y = +0.25 shows no tau* below "
                    "T (M(t) recurs >= 0 to late times).  But n_y = +1 "
                    "behaves IDENTICALLY under the criterion (also no "
                    "tau* < T; mid-restriction gives the same 15.87 for "
                    "both, since v_x(d,1/2,t) = A(t) is n_y-independent): "
                    "the criterion does not discriminate the cutoff family "
                    "member from the no-cutoff one, which was this gate's "
                    "purpose."))
    g5["verdict"] = "PARTIAL"

    # axial validity demonstration (diagnostic, not a gate)
    ax = crit_ref["+0.00"]
    axial_demo = dict(
        tau_star=ax["mid"]["tau_star"],
        stability=stab["+0.00"]["mid"],
        sampled_axial_last=float(tau_A.max()),
        gap=float(ax["mid"]["tau_star"] - float(tau_A.max())),
        delta_min_after=ax["mid"]["delta_min_after_tau01"],
        note=("For n_y = 0 the trajectories keep z fixed and v_x is "
              "z-independent, so sup_z v_x = A(t) and the criterion "
              "APPLIES VALIDLY: axial continuum first-crossing support "
              "within [0,16] is contained in [0, tau*_axial] -- a true "
              "field statement, sharp to the sampled last arrival.  The "
              "criterion works exactly where the spin term is absent."))

    results = dict(
        workstream="N3", file_id="F-T7-N3",
        params=dict(L=L, sigma_x=SIGMA_X, x0=X0, k0=K0,
                    taper="S2 smoothstep, W=1 for |x-x0|<=2.5 sigma, 0 at 3.5",
                    amplitude="literal exp(-(x-x0)^2/(2 sigma^2)) * W",
                    box=[XMIN, XMAX], Nx_ref=NX_REF, dt_ref=DT_REF,
                    dt_out=DT_OUT, T=T_FINAL, d=D_NEAR, nz_line=NZ_LINE,
                    floors=FLOORS, ny_gate=NY_GATE,
                    ny_cover=[float(v) for v in NY_COVER],
                    vclamp=VCLAMP, rho_eps_rel=RHO_EPS_REL),
        field_structure=dict(
            identity=("v_x(d,z,t;n_y) = A(t) + n_y S(z) with "
                      "A(t) = Im(f'/f)(d,t) [z-independent] and "
                      "S(z) = -pi cot(pi z) [t-independent]; exact for the "
                      "factorized state.  S has a positive branch (z->1 for "
                      "n_y>0, z->0 for n_y<0), unbounded until the engine "
                      "clamp; S(1/2) = 0 for every n_y."),
            A_stats=dict(frac_t_nonneg=float((A_ref >= 0).mean()),
                         last_t_nonneg=float(tgrid[A_ref >= 0][-1]),
                         A_min_kill_bins=float(A_ref[m_kill].min()),
                         A_max_kill_bins=float(A_ref[m_kill].max()),
                         A_at_T=float(A_ref[-1]),
                         rho_line_at_T=float(RHO_ref[-1]))),
        discretization_bounds=bounds,
        criterion=crit_ref,
        tau_star_stability=stab,
        consistency=consistency,
        axial_validity_demo=axial_demo,
        covering_diagnostic=cover,
        G1=g1, G2=g2, G3=g3, G4=g4, G5=g5,
        elapsed_total=time.time() - t0,
    )
    with open(os.path.join(OUTDIR, "n3_results.json"), "w") as fh:
        json.dump(results, fh, indent=2)
    np.savez_compressed(
        os.path.join(OUTDIR, "n3_raw.npz"),
        tgrid=tgrid, A_ref=A_ref, A_2x=A_2x, A_st=A_st, A_st2=A_st2,
        RHO_ref=RHO_ref, heat_p1=heat_p1, zg_line=zg_line,
        b_Nx=b_Nx, b_dt=b_dt,
        tau_cover=tau_cover, P0_cover=P0_cover, wgt_cover=wgt_cover,
        tau_axial_sampled=tau_A)
    print("verdicts:", {g: results[g]["verdict"]
                        for g in ("G1", "G2", "G3", "G4", "G5")})
    print("total elapsed %.1f s" % results["elapsed_total"])


if __name__ == "__main__":
    main()
