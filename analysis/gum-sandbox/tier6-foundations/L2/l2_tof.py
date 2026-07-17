#!/usr/bin/env python3
"""
WORKSTREAM L2 -- Momentum time-of-flight: dBB exact + Nelson leg.

Bricmont Appendix-1 momentum ToF for the 1D Gaussian ground-state packet
under free evolution, done twice:
  (a) deterministic de Broglie-Bohm (dX/dt = v),
  (b) Nelson stochastic mechanics (dX = (v+u)dt + sqrt(2 nu) dW),
using the closed-form free-packet solution (no PDE solve).

Units: hbar = m = 1.  Psi(x,0) = pi^{-1/4} exp(-x^2/2).
  Psi(x,t) = (1+it)^{-1/2} pi^{-1/4} exp(-x^2/(2(1+it)))
  rho(x,t) = pi^{-1/2} (1+t^2)^{-1/2} exp(-x^2/(1+t^2))
  v(x,t)   = t x / (1+t^2)            (current velocity)
  u(x,t)   = -2 nu x / (1+t^2)        (osmotic velocity, u = nu d/dx ln rho)
ToF estimator: p = X(T)/T, target |FT Psi|^2 = pi^{-1/2} exp(-p^2),
i.e. Gaussian(0, sigma = 1/sqrt(2)).
"""
import json
import numpy as np
from scipy import stats

RNG_SEED = 20260717
T_FINAL = 100.0
N_PART = 20000
N_TRAJ_G1 = 1000
NU_BASE = 0.5
NU_LIST = [0.1, 0.5, 1.0]
SIGMA_P = 1.0 / np.sqrt(2.0)

OUTDIR = "/home/user/fork_frankensim/analysis/gum-sandbox/tier6-foundations/L2"

# ----------------------------------------------------------------------
# closed forms
def v_field(x, t):
    return t * x / (1.0 + t * t)

def u_field(x, t, nu):
    return -2.0 * nu * x / (1.0 + t * t)

def rho(x, t):
    s2 = 1.0 + t * t
    return np.exp(-x * x / s2) / np.sqrt(np.pi * s2)

# ----------------------------------------------------------------------
# time grids
def dbb_time_grid():
    """dt = 1e-3 for t < 10, dt = 1e-2 for 10 <= t <= 100 (RK4)."""
    steps = [(1e-3, 10000), (1e-2, 9000)]
    return steps

def nelson_time_grid():
    """dt = 1e-3 for t < 10, then geometric growth (x1.02) capped at 0.05."""
    ts = []
    dts = []
    t = 0.0
    dt = 1e-3
    n_early = 10000
    for _ in range(n_early):
        ts.append(t)
        dts.append(dt)
        t += dt
    t = 10.0  # exact
    while t < T_FINAL - 1e-12:
        dt = min(dt * 1.02, 0.05, T_FINAL - t)
        ts.append(t)
        dts.append(dt)
        t += dt
    return np.array(ts), np.array(dts)

# ----------------------------------------------------------------------
# dBB leg: RK4 on dX/dt = v(X,t)
def integrate_dbb(x0):
    x = x0.copy()
    t = 0.0
    for dt, nsteps in dbb_time_grid():
        for _ in range(nsteps):
            k1 = v_field(x, t)
            k2 = v_field(x + 0.5 * dt * k1, t + 0.5 * dt)
            k3 = v_field(x + 0.5 * dt * k2, t + 0.5 * dt)
            k4 = v_field(x + dt * k3, t + dt)
            x = x + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
            t += dt
        t = round(t, 10)
    return x, t

# ----------------------------------------------------------------------
# Nelson leg: Euler-Maruyama on dX = (v+u)dt + sqrt(2 nu) dW
def integrate_nelson(x0, nu, rng, snap_times=(1.0, 10.0, 100.0)):
    ts, dts = nelson_time_grid()
    x = x0.copy()
    snaps = {}
    snap_set = sorted(snap_times)
    sqrt2nu = np.sqrt(2.0 * nu)
    for t, dt in zip(ts, dts):
        # record snapshot if we are exactly at a requested time
        while snap_set and abs(t - snap_set[0]) < 1e-9:
            snaps[snap_set.pop(0)] = x.copy()
        drift = v_field(x, t) + u_field(x, t, nu)
        x = x + drift * dt + sqrt2nu * np.sqrt(dt) * rng.standard_normal(x.size)
        t_end = t + dt
    while snap_set and abs(t_end - snap_set[0]) < 1e-9:
        snaps[snap_set.pop(0)] = x.copy()
    assert not snap_set, f"missed snapshots: {snap_set}"
    return x, snaps

# ----------------------------------------------------------------------
# G3 helper: binned density vs analytic, with statistical floor
def equivariance_check(x, t, n_bins=61):
    s = np.sqrt((1.0 + t * t) / 2.0)          # std of rho(.,t)
    edges = np.linspace(-4.5 * s, 4.5 * s, n_bins + 1)
    h = edges[1] - edges[0]
    counts, _ = np.histogram(x, bins=edges)
    rho_hat = counts / (x.size * h)
    centers = 0.5 * (edges[:-1] + edges[1:])
    rho_an = rho(centers, t)
    mask = rho_an > 1e-3 * rho_an.max()
    # expected counts per bin and the ~4/sqrt(N_per_bin) relative floor
    n_per_bin = x.size * rho_an * h
    floor_abs = 4.0 * np.sqrt(rho_an / (x.size * h))   # = 4 rho / sqrt(N_per_bin)
    diff = np.abs(rho_hat - rho_an)
    sup = float(np.max(diff[mask]))
    ratio = float(np.max(diff[mask] / floor_abs[mask]))
    return {
        "t": t,
        "n_bins_used": int(mask.sum()),
        "bin_width": float(h),
        "sup_norm": sup,
        "floor_at_argmax": float(floor_abs[mask][np.argmax(diff[mask])]),
        "min_expected_counts_in_mask": float(n_per_bin[mask].min()),
        "max_ratio_diff_over_floor": ratio,
        "pass": bool(ratio <= 1.0),
        "_plot": (centers, rho_hat, rho_an, mask),
    }

# ----------------------------------------------------------------------
def main():
    rng = np.random.default_rng(RNG_SEED)
    results = {"seed": RNG_SEED, "T": T_FINAL, "N": N_PART,
               "sigma_p_target": SIGMA_P}

    # initial conditions X(0) ~ rho(.,0) = N(0, 1/2)
    x0 = rng.normal(0.0, np.sqrt(0.5), N_PART)

    # ---------------- dBB leg ----------------
    xT, t_end = integrate_dbb(x0)
    assert abs(t_end - T_FINAL) < 1e-8, t_end
    stretch = np.sqrt(1.0 + T_FINAL**2)
    x_exact = x0 * stretch

    # G1: trajectory law over first 1000 trajectories
    err = xT[:N_TRAJ_G1] - x_exact[:N_TRAJ_G1]
    rms = float(np.sqrt(np.mean(err**2)))
    g1_tol = 1e-6 * stretch
    results["G1"] = {"rms_error": rms, "tol": float(g1_tol),
                     "max_abs_error": float(np.max(np.abs(err))),
                     "pass": bool(rms <= g1_tol)}

    # G2: dBB ToF KS vs N(0, 1/sqrt(2))
    p_dbb = xT / T_FINAL
    ks2 = stats.kstest(p_dbb, "norm", args=(0.0, SIGMA_P))
    results["G2"] = {"KS_D": float(ks2.statistic), "p_value": float(ks2.pvalue),
                     "sample_std": float(p_dbb.std(ddof=1)),
                     "pass": bool(ks2.pvalue > 0.1)}

    # note: instantaneous dBB momentum at t=0 (Psi real) is exactly zero
    results["dbb_v_at_t0_max"] = float(np.max(np.abs(v_field(x0, 0.0))))

    # ---------------- Nelson leg (baseline nu = 1/2) ----------------
    rng_n = np.random.default_rng(RNG_SEED + 1)
    x0_n = rng_n.normal(0.0, np.sqrt(0.5), N_PART)
    xT_n, snaps = integrate_nelson(x0_n, NU_BASE, rng_n)

    # G3: equivariance at t = 1, 10, 100
    g3 = {}
    g3_plots = {}
    for t_snap in (1.0, 10.0, 100.0):
        x_snap = snaps[t_snap] if t_snap in snaps else xT_n
        chk = equivariance_check(x_snap, t_snap)
        g3_plots[t_snap] = chk.pop("_plot")
        g3[str(t_snap)] = chk
    results["G3"] = {"snapshots": g3,
                     "pass": bool(all(v["pass"] for v in g3.values()))}

    # G4: Nelson ToF KS vs N(0, 1/sqrt(2))
    p_nel = xT_n / T_FINAL
    ks4 = stats.kstest(p_nel, "norm", args=(0.0, SIGMA_P))
    results["G4"] = {"nu": NU_BASE, "KS_D": float(ks4.statistic),
                     "p_value": float(ks4.pvalue),
                     "sample_std": float(p_nel.std(ddof=1)),
                     "pass": bool(ks4.pvalue > 0.1)}

    # G5: nu-dial -- Nelson ToF for nu in {0.1, 0.5, 1.0}
    p_by_nu = {}
    one_sample = {}
    for i, nu in enumerate(NU_LIST):
        if nu == NU_BASE:
            p_by_nu[nu] = p_nel
            one_sample[nu] = {"KS_D": float(ks4.statistic),
                              "p_value": float(ks4.pvalue)}
            continue
        rng_i = np.random.default_rng(RNG_SEED + 10 + i)
        x0_i = rng_i.normal(0.0, np.sqrt(0.5), N_PART)
        xT_i, _ = integrate_nelson(x0_i, nu, rng_i, snap_times=(100.0,))
        p_by_nu[nu] = xT_i / T_FINAL
        ks_i = stats.kstest(p_by_nu[nu], "norm", args=(0.0, SIGMA_P))
        one_sample[nu] = {"KS_D": float(ks_i.statistic),
                          "p_value": float(ks_i.pvalue)}
    pairwise = {}
    for a in range(len(NU_LIST)):
        for b in range(a + 1, len(NU_LIST)):
            na, nb = NU_LIST[a], NU_LIST[b]
            ks = stats.ks_2samp(p_by_nu[na], p_by_nu[nb])
            pairwise[f"{na}_vs_{nb}"] = {"KS_D": float(ks.statistic),
                                         "p_value": float(ks.pvalue)}
    results["G5"] = {"one_sample_vs_target": {str(k): v for k, v in one_sample.items()},
                     "pairwise": pairwise,
                     "pass": bool(all(v["p_value"] > 0.05 for v in pairwise.values()))}

    # ---------------- save ----------------
    with open(f"{OUTDIR}/L2_results.json", "w") as f:
        json.dump(results, f, indent=2)

    np.savez_compressed(f"{OUTDIR}/l2_samples.npz",
                        p_dbb=p_dbb, p_nel=p_nel,
                        p_nu01=p_by_nu[0.1], p_nu10=p_by_nu[1.0],
                        snap1=snaps[1.0], snap10=snaps[10.0])

    make_figure(p_dbb, p_nel, p_by_nu, g3_plots)

    for k in ("G1", "G2", "G3", "G4", "G5"):
        print(k, "PASS" if results[k]["pass"] else "FAIL")
    print(json.dumps({k: results[k] for k in ("G1", "G2", "G4")}, indent=2))
    print("G3 ratios:", {k: v["max_ratio_diff_over_floor"]
                         for k, v in results["G3"]["snapshots"].items()})
    print("G5 pairwise:", results["G5"]["pairwise"])

# ----------------------------------------------------------------------
def make_figure(p_dbb, p_nel, p_by_nu, g3_plots):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # palette (dataviz reference, light mode)
    BLUE, GREEN, MAGENTA, YELLOW = "#2a78d6", "#008300", "#e87ba4", "#eda100"
    INK, MUTED = "#0b0b0b", "#52514e"
    SURF = "#fcfcfb"

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.2), facecolor=SURF)
    for ax in axes.flat:
        ax.set_facecolor(SURF)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.grid(True, color="#e6e5e1", lw=0.8)
        ax.set_axisbelow(True)
        ax.tick_params(colors=MUTED, labelsize=9)

    pgrid = np.linspace(-3.2, 3.2, 400)
    target = np.exp(-pgrid**2) / np.sqrt(np.pi)

    # (a) dBB ToF
    ax = axes[0, 0]
    ax.hist(p_dbb, bins=80, density=True, color=BLUE, alpha=0.55,
            edgecolor=SURF, linewidth=0.3, label="dBB  p = X(T)/T")
    ax.plot(pgrid, target, color=INK, lw=2, label=r"$|\tilde\Psi(p)|^2$")
    ax.set_title("(a) dBB time-of-flight, T = 100  (initial velocity exactly 0)",
                 fontsize=10, color=INK, loc="left")
    ax.set_xlabel("p", color=MUTED); ax.set_ylabel("density", color=MUTED)
    ax.legend(frameon=False, fontsize=9)

    # (b) Nelson ToF, nu = 1/2
    ax = axes[0, 1]
    ax.hist(p_nel, bins=80, density=True, color=GREEN, alpha=0.55,
            edgecolor=SURF, linewidth=0.3, label=r"Nelson $\nu=1/2$")
    ax.plot(pgrid, target, color=INK, lw=2, label=r"$|\tilde\Psi(p)|^2$")
    ax.set_title(r"(b) Nelson time-of-flight, $\nu = 1/2$, T = 100",
                 fontsize=10, color=INK, loc="left")
    ax.set_xlabel("p", color=MUTED); ax.set_ylabel("density", color=MUTED)
    ax.legend(frameon=False, fontsize=9)

    # (c) equivariance
    ax = axes[1, 0]
    cols = {1.0: BLUE, 10.0: GREEN, 100.0: MAGENTA}
    for t_snap in (1.0, 10.0, 100.0):
        centers, rho_hat, rho_an, mask = g3_plots[t_snap]
        s = np.sqrt((1 + t_snap**2) / 2)
        ax.step(centers / s, rho_hat * s, where="mid", color=cols[t_snap],
                lw=1.6, label=f"binned, t = {t_snap:g}")
        ax.plot(centers / s, rho_an * s, color=INK, lw=1.0, ls="--")
    ax.set_title(r"(c) Nelson equivariance: binned vs $\rho(x,t)$ (scaled units)",
                 fontsize=10, color=INK, loc="left")
    ax.set_xlabel(r"$x/\sigma(t)$", color=MUTED)
    ax.set_ylabel(r"$\sigma\,\rho$", color=MUTED)
    ax.legend(frameon=False, fontsize=9)

    # (d) nu-dial CDFs
    ax = axes[1, 1]
    cols_nu = {0.1: BLUE, 0.5: GREEN, 1.0: YELLOW}
    for nu in (0.1, 0.5, 1.0):
        q = np.sort(p_by_nu[nu])
        ax.plot(q, np.arange(1, q.size + 1) / q.size, color=cols_nu[nu],
                lw=1.8, label=rf"$\nu = {nu:g}$")
    ax.plot(pgrid, stats.norm.cdf(pgrid, 0, 1 / np.sqrt(2)), color=INK,
            lw=1.0, ls="--", label=r"$\mathcal{N}(0,\,1/\sqrt{2})$ CDF")
    ax.set_xlim(-3.2, 3.2)
    ax.set_title(r"(d) $\nu$-dial: ToF empirical CDFs coincide (diffusion-blind)",
                 fontsize=10, color=INK, loc="left")
    ax.set_xlabel("p", color=MUTED); ax.set_ylabel("CDF", color=MUTED)
    ax.legend(frameon=False, fontsize=9, loc="lower right")

    fig.suptitle("L2 — Momentum time-of-flight: dBB exact + Nelson leg "
                 "(within-model)", fontsize=12, color=INK, x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(f"{OUTDIR}/L2_tof.png", dpi=160)
    plt.close(fig)

if __name__ == "__main__":
    main()
