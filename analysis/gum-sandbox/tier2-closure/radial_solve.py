#!/usr/bin/env python3
"""
Tier-2b STEP 1 of the GUM replication program:
1-D radial (hedgehog) FIELD-LEVEL solve at eps = 0.05.

Hedgehog ansatz P(x) = exp(i f(r) x_hat.sigma), f(0)=pi, f(inf)=0.

Sector energies (radial integrals), total E = a2 E2 + a4 E4 + a6 E6 + a0 E0:
  E2 = INT (f'^2 r^2 + 2 sin^2 f) dr
  E4 = INT sin^2 f (2 f'^2 + sin^2 f / r^2) dr
  E6 = INT sin^4 f f'^2 / r^2 dr
  E0 = INT m^2 (1 - cos f) r^2 dr
Fixed: a6 = a0 = 1, m = 1.  Dial: a2 = a4 = t, tuned so that at the solved
profile  ratio = t (E2+E4) / (E6+E0) = 0.05 +- 0.001.

Method (documented):
  * Midpoint (cell-centered) discretization of the energy on a graded radial
    grid (fine near r=0, near the compacton edge r ~ R* = 2^(5/6), and through
    the tail region).  f_0 = pi and f_N = 0 are hard Dirichlet constraints.
  * Minimization: (a) diagonally preconditioned gradient flow with adaptive
    step (energy-monotone) to reach the basin, then (b) damped (Levenberg)
    Newton on the exact gradient of the discretized energy.  The energy is a
    sum of per-cell terms coupling nearest neighbours only, so the exact
    discrete Hessian is symmetric tridiagonal; each Newton step is a Thomas
    (tridiagonal) solve.  Newton is run to max|grad| < 1e-13 (machine floor),
    which makes the per-step energy change << 1e-12.
  * Outer loop: fixed-point iteration on t (ratio is ~linear in t).
  * Grid-refinement check: everything recomputed at N and 2N (same t).

Derived identities used below (derivations in radial_RESULTS.md):
  * Derrick virial (r -> lambda r):  E2 -> lam E2, E4 -> lam^-1 E4,
    E6 -> lam^-3 E6, E0 -> lam^3 E0  =>  t E2 - t E4 - 3 E6 + 3 E0 = 0.
  * Linearized EL (small f):  t (f'' + 2 f'/r - 2 f/r^2) = (m^2/2) f,
    i.e. mu_true = m / sqrt(2 t).  (The task sheet's quoted formula
    mu = m sqrt(a0/a2) misses the 1/2 from (1-cos f) ~ f^2/2; both reported.)
  * BPS (Bogomolny) bound: E6 + E0 >= 2 sqrt(a6 a0) m sqrt(2) * 16/15
    = 32 sqrt(2)/15  (m = a6 = a0 = 1), saturated by the compacton
    f0 = 2 arccos(r/R*), R* = 2^(5/6).

Deterministic, no RNG. Runtime ~1-2 min. Outputs:
  radial_profile.png, radial_results.json, radial_RESULTS.md is written by hand
  from the printed table (this script prints everything it computes).
"""

import json
import time

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# constants / conventions
# ----------------------------------------------------------------------
M = 1.0                      # potential mass m
A6 = 1.0                     # sextic coefficient a6
A0 = 1.0                     # potential coefficient a0
EPS_TARGET = 0.05            # target sector ratio
EPS_TOL = 5.0e-5             # inner tolerance on the ratio (spec: +-1e-3)
RSTAR = 2.0 ** (5.0 / 6.0)   # BPS compacton radius for a6=a0=m=1
BPS_BOUND = 32.0 * np.sqrt(2.0) / 15.0   # Bogomolny bound on E6+E0

OUTDIR = "/home/user/fork_frankensim/analysis/gum-sandbox/tier2-closure"


# ----------------------------------------------------------------------
# grid
# ----------------------------------------------------------------------
def make_grid(N, rmax):
    """Graded grid on [0, rmax]: N cells, N+1 nodes.

    Node density w(r) is boosted near r=0 (core), near the former compacton
    edge r ~ R* (the known stiff spot), and mildly through the Yukawa-tail
    window; deterministic inverse-CDF construction."""
    s = np.linspace(0.0, 1.0, 200001)
    rr = s * rmax
    w = (0.55
         + 2.0 * np.exp(-(rr / 0.5) ** 2)                    # core
         + 3.0 * np.exp(-((rr - RSTAR) / 0.45) ** 2)          # edge region
         + 0.8 * np.exp(-((rr - (RSTAR + 1.2)) / 0.9) ** 2))  # tail window
    cw = np.concatenate([[0.0], np.cumsum(0.5 * (w[1:] + w[:-1]) * np.diff(rr))])
    cw /= cw[-1]
    r = np.interp(np.linspace(0.0, 1.0, N + 1), cw, rr)
    r[0] = 0.0
    r[-1] = rmax
    return r


def init_profile(r):
    """Smoothed compacton: f = 2 arccos(x), x = u/(1+u^8)^(1/8), u = r/R*.
    Matches the compacton in the interior, smooth power tail outside."""
    u = r / RSTAR
    x = u / (1.0 + u ** 8) ** 0.125
    f = 2.0 * np.arccos(np.clip(x, 0.0, 1.0))
    f[0] = np.pi
    f[-1] = 0.0
    return f


# ----------------------------------------------------------------------
# discretized energy: midpoint rule, exact gradient + tridiagonal Hessian
# ----------------------------------------------------------------------
def assemble(r, f, t, need_hess=True):
    """Return E, (E2,E4,E6,E0), grad (interior nodes), Hessian tri-diagonals.

    Per cell c:  e_c = h_c * L(r_m, f_m, d)  with midpoint values
    f_m = (f_i+f_{i+1})/2, d = (f_{i+1}-f_i)/h_c.
    L = t (d^2 r^2 + 2 s^2) + t s^2 (2 d^2 + s^2/r^2) + s^4 d^2/r^2
        + m^2 (1-cos f) r^2,   s = sin f_m.
    """
    h = np.diff(r)
    rm = 0.5 * (r[1:] + r[:-1])
    rm2 = rm * rm
    fm = 0.5 * (f[1:] + f[:-1])
    d = np.diff(f) / h
    d2 = d * d
    s = np.sin(fm)
    c = np.cos(fm)
    s2 = s * s
    s3 = s2 * s
    s4 = s2 * s2
    sin2f = 2.0 * s * c
    cos2f = c * c - s2

    L2 = d2 * rm2 + 2.0 * s2
    L4 = s2 * (2.0 * d2 + s2 / rm2)
    L6 = s4 * d2 / rm2
    L0 = (M * M) * (1.0 - np.cos(fm)) * rm2
    E2 = float(np.sum(h * L2))
    E4 = float(np.sum(h * L4))
    E6 = float(np.sum(h * L6))
    E0 = float(np.sum(h * L0))
    E = t * (E2 + E4) + A6 * E6 + A0 * E0

    # dL/dd and dL/df at midpoints (weights included: a2=a4=t, a6=a0=1)
    L_d = 2.0 * t * d * rm2 + 4.0 * t * s2 * d + 2.0 * s4 * d / rm2
    L_f = (2.0 * t * sin2f * (1.0 + d2)
           + 4.0 * (t + d2) * s3 * c / rm2
           + (M * M) * s * rm2)

    # node gradient: cell c contributes gL to node c, gR to node c+1
    gL = 0.5 * h * L_f - L_d
    gR = 0.5 * h * L_f + L_d
    g = gR[:-1] + gL[1:]          # interior nodes i = 1..N-1

    if not need_hess:
        return E, (E2, E4, E6, E0), g, None, None

    L_dd = 2.0 * t * rm2 + 4.0 * t * s2 + 2.0 * s4 / rm2
    L_fd = 4.0 * t * d * sin2f + 8.0 * s3 * c * d / rm2
    L_ff = (4.0 * t * cos2f * (1.0 + d2)
            + 4.0 * (t + d2) * (3.0 * s2 * c * c - s4) / rm2
            + (M * M) * c * rm2)

    HLL = 0.25 * h * L_ff - L_fd + L_dd / h
    HRR = 0.25 * h * L_ff + L_fd + L_dd / h
    HLR = 0.25 * h * L_ff - L_dd / h
    diag = HRR[:-1] + HLL[1:]     # interior nodes
    off = HLR[1:-1]               # coupling (i, i+1), i = 1..N-2
    return E, (E2, E4, E6, E0), g, diag, off


def thomas(diag, off, b):
    """Solve symmetric tridiagonal system (diag; off) x = b (Thomas algo)."""
    n = diag.size
    dd = diag.tolist()
    ee = off.tolist()
    bb = b.tolist()
    for i in range(1, n):
        mfac = ee[i - 1] / dd[i - 1]
        dd[i] -= mfac * ee[i - 1]
        bb[i] -= mfac * bb[i - 1]
    x = [0.0] * n
    x[-1] = bb[-1] / dd[-1]
    for i in range(n - 2, -1, -1):
        x[i] = (bb[i] - ee[i] * x[i + 1]) / dd[i]
    return np.asarray(x)


# ----------------------------------------------------------------------
# minimizers
# ----------------------------------------------------------------------
def flow(r, f, t, gmax_target=1e-2, maxsteps=6000):
    """Diagonally preconditioned, energy-monotone gradient flow."""
    E, _, g, diag, _ = assemble(r, f, t)
    dt = 0.5
    for _ in range(maxsteps):
        if np.max(np.abs(g)) < gmax_target:
            break
        prec = np.maximum(diag, 1e-3 * np.max(diag))
        fn = f.copy()
        fn[1:-1] = f[1:-1] - dt * g / prec
        np.clip(fn, 0.0, np.pi, out=fn)
        En, _, gn, diagn, _ = assemble(r, fn, t)
        if En < E:
            f, E, g, diag = fn, En, gn, diagn
            dt = min(dt * 1.1, 5.0)
        else:
            dt *= 0.5
            if dt < 1e-9:
                break
    return f


def newton(r, f, t, gtol=1e-13, maxit=300):
    """Damped (Levenberg) Newton on grad(E)=0; exact tridiagonal Hessian.

    Accepts a step only if max|grad| decreases; damping lambda adapts.
    Returns f, final sector energies, achieved max|grad|, last |dE| per step,
    and iteration count."""
    E, sect, g, diag, off = assemble(r, f, t)
    lam = 1e-3
    dE_last = np.inf
    it = 0
    while it < maxit:
        gmax = float(np.max(np.abs(g)))
        if gmax < gtol:
            break
        scale = np.abs(diag) + 1e-30
        accepted = False
        while lam < 1e14:
            step = thomas(diag + lam * scale, off, -g)
            if np.all(np.isfinite(step)):
                fn = f.copy()
                fn[1:-1] += step
                En, sectn, gn, diagn, offn = assemble(r, fn, t)
                gmaxn = float(np.max(np.abs(gn)))
                if np.isfinite(En) and gmaxn < gmax:
                    dE_last = abs(E - En)
                    f, E, sect, g, diag, off = fn, En, sectn, gn, diagn, offn
                    lam = max(lam * 0.25, 1e-14)
                    accepted = True
                    break
            lam *= 10.0
        if not accepted:
            break   # stalled at machine floor
        it += 1
    gmax = float(np.max(np.abs(g)))
    return f, E, sect, g, gmax, dE_last, it


def solve_profile(r, f0, t):
    f = flow(r, f0, t, gmax_target=1e-2)
    return newton(r, f, t)


# ----------------------------------------------------------------------
# diagnostics
# ----------------------------------------------------------------------
def degree(r, f):
    """K = -(2/pi) INT sin^2 f f' dr  (signed; midpoint quadrature)."""
    fm = 0.5 * (f[1:] + f[:-1])
    return float(-(2.0 / np.pi) * np.sum(np.sin(fm) ** 2 * np.diff(f)))


def tail_fit(r, f, mu_lo, mu_hi):
    """Fit f ~ A (1/(mu r)^2 + 1/(mu r)) e^{-mu r} on the window
    f in [1e-6, 1e-2] (log-space least squares; golden section over mu,
    analytic solve for log A)."""
    mask = (f > 1e-6) & (f < 1e-2) & (r > RSTAR)
    rw = r[mask]
    y = np.log(f[mask])

    def misfit(mu):
        x = mu * rw
        model = -x + np.log(1.0 / x + 1.0 / (x * x))
        lgA = float(np.mean(y - model))
        res = y - model - lgA
        return float(np.sqrt(np.mean(res ** 2))), lgA

    invphi = (np.sqrt(5.0) - 1.0) / 2.0
    a, b = mu_lo, mu_hi
    c = b - invphi * (b - a)
    d = a + invphi * (b - a)
    fc, _ = misfit(c)
    fd, _ = misfit(d)
    for _ in range(120):
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - invphi * (b - a)
            fc, _ = misfit(c)
        else:
            a, c, fc = c, d, fd
            d = a + invphi * (b - a)
            fd, _ = misfit(d)
    mu = 0.5 * (a + b)
    rms, lgA = misfit(mu)
    return mu, float(np.exp(lgA)), rms, int(mask.sum()), (float(rw[0]), float(rw[-1]))


def diagnostics(r, f, t, gmax, dE_last, n_newton):
    E, (E2, E4, E6, E0), g, _, _ = assemble(r, f, t)
    ratio = t * (E2 + E4) / (A6 * E6 + A0 * E0)
    K = degree(r, f)
    vir = t * E2 - t * E4 - 3.0 * A6 * E6 + 3.0 * A0 * E0
    vir_rel = abs(vir) / E
    # continuum EL residual: discrete gradient / nodal quadrature weight
    w = 0.5 * (r[2:] - r[:-2])
    el = g / w
    el_max = float(np.max(np.abs(el)))
    mu_true = M / np.sqrt(2.0 * t)     # derived linearization
    mu_spec = M / np.sqrt(t)           # task-sheet formula (no 1/2)
    mu_fit, Afit, rms, npts, win = tail_fit(r, f, 0.3 * mu_true, 2.5 * mu_true)
    # tail sanity
    i_chk = np.searchsorted(r, 0.97 * r[-1])
    f_far = float(np.max(np.abs(f[i_chk:])))
    mono = bool(np.all(np.diff(f) <= 1e-12))
    return {
        "N_cells": int(len(r) - 1),
        "Rmax": float(r[-1]),
        "t": float(t),
        "E2": E2, "E4": E4, "E6": E6, "E0": E0,
        "E_total": float(E),
        "eps_ratio": float(ratio),
        "degree_K": K,
        "virial": float(vir),
        "virial_rel": float(vir_rel),
        "EL_residual_max": el_max,
        "grad_max": gmax,
        "dE_last_step": float(dE_last),
        "newton_iters": n_newton,
        "mu_fit": float(mu_fit),
        "mu_analytic_derived": float(mu_true),
        "mu_task_sheet": float(mu_spec),
        "A_fit": float(Afit),
        "tail_rms_logresid": float(rms),
        "tail_npts": npts,
        "tail_window_r": win,
        "A_over_mu": float(Afit / mu_fit),
        "A_over_mu2": float(Afit / mu_fit ** 2),
        "A_times_mu": float(Afit * mu_fit),
        "A_times_mu2": float(Afit * mu_fit ** 2),
        # order-1 candidates (need e^{-mu R*} to tame the k_1 normalization)
        "A_exp_muRstar": float(Afit * np.exp(-mu_fit * RSTAR)),
        "A_exp_over_mu": float(Afit * np.exp(-mu_fit * RSTAR) / mu_fit),
        "A_exp_over_muRstar": float(Afit * np.exp(-mu_fit * RSTAR)
                                    / (mu_fit * RSTAR)),
        "A_exp_over_2pi": float(Afit * np.exp(-mu_fit * RSTAR)
                                / (2.0 * np.pi)),
        "f_fit_at_Rstar": float(Afit * np.exp(-mu_fit * RSTAR)
                                * (1.0 / (mu_fit * RSTAR)
                                   + 1.0 / (mu_fit * RSTAR) ** 2)),
        "BPS_bound": float(BPS_BOUND),
        "BPS_deficit": float((A6 * E6 + A0 * E0) / BPS_BOUND - 1.0),
        "f_beyond_0.97Rmax_max": f_far,
        "f_monotone_decreasing": mono,
    }


# ----------------------------------------------------------------------
# main
# ----------------------------------------------------------------------
def main():
    t_start = time.time()
    N1 = 4000
    rmax = 6.0

    # --- outer eps-dial loop at N1 ------------------------------------
    r1 = make_grid(N1, rmax)
    f1 = init_profile(r1)
    t = 0.008   # initial guess for a2=a4
    history = []
    for outer in range(15):
        f1, E, sect, g, gmax, dE_last, nit = solve_profile(r1, f1, t)
        E2, E4, E6, E0 = sect
        ratio = t * (E2 + E4) / (A6 * E6 + A0 * E0)
        history.append((t, ratio))
        print(f"[outer {outer}] t = {t:.10f}  ratio = {ratio:.7f}  "
              f"gmax = {gmax:.2e}  newton_its = {nit}")
        if abs(ratio - EPS_TARGET) < EPS_TOL:
            break
        t *= EPS_TARGET / ratio

    # --- Rmax check ----------------------------------------------------
    i_chk = np.searchsorted(r1, 0.97 * rmax)
    assert np.max(np.abs(f1[i_chk:])) < 1e-8, "increase Rmax"

    D1 = diagnostics(r1, f1, t, gmax, dE_last, nit)

    # --- grid refinement: 2N, same t ------------------------------------
    N2 = 2 * N1
    r2 = make_grid(N2, rmax)
    f2 = np.interp(r2, r1, f1)
    f2[0] = np.pi
    f2[-1] = 0.0
    f2, E, sect, g, gmax2, dE_last2, nit2 = solve_profile(r2, f2, t)
    D2 = diagnostics(r2, f2, t, gmax2, dE_last2, nit2)

    # --- report ----------------------------------------------------------
    keys = ["t", "E2", "E4", "E6", "E0", "E_total", "eps_ratio", "degree_K",
            "virial", "virial_rel", "EL_residual_max", "grad_max",
            "dE_last_step", "mu_fit", "mu_analytic_derived", "mu_task_sheet",
            "A_fit", "A_over_mu", "A_over_mu2", "A_times_mu", "A_times_mu2",
            "A_exp_muRstar", "A_exp_over_mu", "A_exp_over_muRstar",
            "A_exp_over_2pi", "f_fit_at_Rstar",
            "tail_rms_logresid", "BPS_bound", "BPS_deficit",
            "f_beyond_0.97Rmax_max"]
    print("\n===== results (N = %d vs 2N = %d cells) =====" % (N1, N2))
    print(f"{'quantity':26s} {'N':>22s} {'2N':>22s} {'|delta|':>12s}")
    for k in keys:
        v1, v2 = D1[k], D2[k]
        print(f"{k:26s} {v1:22.12g} {v2:22.12g} {abs(v2 - v1):12.3e}")
    print("degree K (2N), 6 dp: %.6f" % D2["degree_K"])
    print("monotone decreasing: N=%s 2N=%s" %
          (D1["f_monotone_decreasing"], D2["f_monotone_decreasing"]))
    print("runtime: %.1f s" % (time.time() - t_start))

    out = {"N_run": D1, "2N_run": D2, "eps_dial_history": history,
           "RSTAR": float(RSTAR), "conventions":
           {"a6": A6, "a0": A0, "m": M, "a2=a4": t,
            "virial_identity": "t*E2 - t*E4 - 3*E6 + 3*E0 = 0",
            "mu_derived": "mu = m/sqrt(2 a2)",
            "BPS_bound": "E6+E0 >= 32*sqrt(2)/15 (a6=a0=m=1)"}}
    with open(OUTDIR + "/radial_results.json", "w") as fh:
        json.dump(out, fh, indent=2)

    # --- figure ----------------------------------------------------------
    make_figure(r2, f2, D2)
    print("wrote radial_results.json, radial_profile.png")


def make_figure(r, f, D):
    col_sol = "#2a78d6"   # series 1 (blue): eps = 0.05 solution
    col_cmp = "#e34948"   # series 6 (red, dashed): eps = 0 compacton
    col_fit = "#4a3aa7"   # violet, dotted: tail fit
    ink = "#333333"
    grid_c = "#dddddd"

    fig, ax = plt.subplots(figsize=(8.0, 5.4), dpi=160)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    rc = np.linspace(0.0, RSTAR, 800)
    fc = 2.0 * np.arccos(rc / RSTAR)
    ax.plot(rc, fc, ls="--", lw=2.0, color=col_cmp,
            label=r"$\epsilon=0$ compacton $2\arccos(r/R_*)$, $R_*=2^{5/6}$")
    ax.plot([RSTAR, 4.0], [0.0, 0.0], ls="--", lw=2.0, color=col_cmp)
    ax.plot(r, f, lw=2.0, color=col_sol,
            label=rf"$\epsilon=0.05$ solution ($t=a_2=a_4={D['t']:.5f}$)")
    ax.axvline(RSTAR, color=grid_c, lw=1.0, zorder=0)
    ax.text(RSTAR + 0.03, 0.35, r"$R_*$", color="#888888", fontsize=10)

    ax.set_xlim(0.0, 4.0)
    ax.set_ylim(-0.08, np.pi + 0.08)
    ax.set_xlabel(r"$r$  (units of $\tilde\Lambda^{-1}$, $m=a_6=a_0=1$)",
                  color=ink)
    ax.set_ylabel(r"$f(r)$", color=ink)
    ax.set_yticks([0, np.pi / 2, np.pi], ["0", r"$\pi/2$", r"$\pi$"])
    ax.grid(color=grid_c, lw=0.6)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(colors=ink)
    ax.set_title("Tier-2b hedgehog profile, field-level solve at "
                 r"$\epsilon = 0.05$  ($B=1$)", color=ink, fontsize=12)
    ax.legend(frameon=False, loc="upper right", fontsize=9)

    # log inset: tail + fit
    axi = ax.inset_axes([0.47, 0.28, 0.5, 0.42])
    mask = (f > 1e-9) & (r > 1.4)
    axi.semilogy(r[mask], f[mask], lw=1.8, color=col_sol)
    mu, A = D["mu_fit"], D["A_fit"]
    rr = np.linspace(*D["tail_window_r"], 300)
    x = mu * rr
    axi.semilogy(rr, A * (1.0 / x + 1.0 / x ** 2) * np.exp(-x),
                 ls=":", lw=2.2, color=col_fit)
    axi.text(0.05, 0.08,
             rf"fit $A\,k_1(\mu r)$: $\mu={mu:.3f}$"
             "\n" rf"($m\sqrt{{a_0/2a_2}}={D['mu_analytic_derived']:.3f}$)",
             transform=axi.transAxes, fontsize=8, color=col_fit)
    axi.axhspan(1e-6, 1e-2, color="#f2f2f2", zorder=0)
    axi.set_ylim(1e-9, 3)
    axi.set_xlim(1.4, 4.2)
    axi.set_ylabel(r"$f(r)$ (log)", fontsize=8, color=ink)
    axi.tick_params(labelsize=8, colors=ink)
    axi.grid(color=grid_c, lw=0.5)
    for sp in ("top", "right"):
        axi.spines[sp].set_visible(False)
    axi.set_title("tail, fit window shaded", fontsize=8, color=ink)

    fig.tight_layout()
    fig.savefig(OUTDIR + "/radial_profile.png", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
