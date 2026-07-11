#!/usr/bin/env python3
"""
Tier-3 follow-up: Born-rule relaxation for M = 4 with a NON-DEGENERATE mode set.

The pilot's M = 4 case used the first-2x2 modes, whose spectrum
E in {1, 2.5, 2.5, 4} has a degeneracy and commensurate gaps (few distinct
beat frequencies) -- the corpus's "small-M near-integrable" caveat case, and
it did not relax near-exponentially. This rerun replaces it with a generic set

    S = {(1,2), (2,3), (3,1), (1,4)},  E in {2.5, 6.5, 5.0, 8.5}

(all distinct; all pairwise gaps distinct), keeping the protocol otherwise
identical: N = 20,000, phases seed 42, sampling seed 12345 (+M offset, i.e.
the SAME initial ensemble as the original M = 4 run), rho_0 = |phi_11|^2,
RK4 dt = 2e-3, t_final = 4*pi, H at 32x32 and 16x16.

Outputs: m4generic_hdata.csv, m4generic_decay.png, m4generic_results.json
"""

import json
import os
import time

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from born import (CGS, DT, L, N_PART, SEED_PHASES, SEED_SAMPLE, T_FINAL,
                  ModeSystem, fit_tau, run_ensemble, sample_phi11)

OUTDIR = os.path.dirname(os.path.abspath(__file__))
MODES = [(1, 2), (2, 3), (3, 1), (1, 4)]


class GenericModeSystem(ModeSystem):
    """ModeSystem over an arbitrary list of (m, n) modes (equal weights,
    fixed random phases, seed 42), reusing the separable-bilinear velocity
    and grid evaluators of the base class via coeff_matrix()."""

    def __init__(self, modes):
        mm = np.array([m for m, _ in modes])
        nn = np.array([n for _, n in modes])
        self.kmax = int(max(mm.max(), nn.max()))
        self.mm = mm.astype(float)
        self.nn = nn.astype(float)
        self.im = mm - 1
        self.jn = nn - 1
        self.M = len(modes)
        self.E = 0.5 * (self.mm ** 2 + self.nn ** 2)
        rng = np.random.default_rng(SEED_PHASES)
        self.theta = rng.uniform(0.0, 2.0 * np.pi, self.M)
        self.ks = np.arange(1, self.kmax + 1, dtype=float)

    def coeff_matrix(self, t):
        C = np.zeros((self.kmax, self.kmax), dtype=complex)
        C[self.im, self.jn] = self.coeffs(t)
        return C


def main():
    t0 = time.time()
    sysg = GenericModeSystem(MODES)
    assert len(set(sysg.E.tolist())) == sysg.M, "energies must be distinct"
    print(f"generic M=4 set {MODES}, E = {sorted(sysg.E.tolist())}", flush=True)

    # identical initial ensemble to the original (degenerate) M = 4 run
    rng = np.random.default_rng(SEED_SAMPLE + sysg.M)
    x, y = sample_phi11(N_PART, rng)
    ts, hcg, _ = run_ensemble(sysg, x, y, "M=4 generic")

    results = {"modes": MODES, "E": sysg.E.tolist(), "N": N_PART, "dt": DT,
               "t_final": float(T_FINAL), "seed_phases": SEED_PHASES,
               "seed_sample": SEED_SAMPLE, "runs": {}}
    for cg in CGS:
        tau, h0fit, t_end, r2 = fit_tau(ts, hcg[cg])
        results["runs"][f"cg{cg}"] = {
            "H_initial": float(hcg[cg][0]), "H_final": float(hcg[cg][-1]),
            "tau": float(tau), "H0_fit": h0fit,
            "fit_window_t_end": t_end, "fit_r2": r2}
    results["noise_floor"] = {f"cg{cg}": cg * cg / (2.0 * N_PART) for cg in CGS}

    # ------------------------------------------------------------------
    # m4generic_hdata.csv
    # ------------------------------------------------------------------
    arr = np.column_stack([ts] + [hcg[cg] for cg in CGS])
    np.savetxt(os.path.join(OUTDIR, "m4generic_hdata.csv"), arr,
               delimiter=",", header="t,H32_M4generic,H16_M4generic",
               comments="", fmt="%.8g")

    # ------------------------------------------------------------------
    # overlay figure: degenerate M4 vs generic M4 vs old M9 (from hdata.csv)
    # ------------------------------------------------------------------
    old = np.genfromtxt(os.path.join(OUTDIR, "hdata.csv"), delimiter=",",
                        names=True)
    series = {  # label -> (color, {cg: curve})
        "M = 4, degenerate (first 2$\\times$2)":
            ("#2a78d6", {32: old["H32_M4"], 16: old["H16_M4"]}),
        "M = 4, generic set":
            ("#1baf7a", {32: hcg[32], 16: hcg[16]}),
        "M = 9 (first 3$\\times$3)":
            ("#eda100", {32: old["H32_M9"], 16: old["H16_M9"]}),
    }
    ink, muted = "#0b0b0b", "#52514e"
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6), sharey=True)
    for ax, cg in zip(axes, CGS):
        for lab, (c, curves) in series.items():
            ax.plot(ts, curves[cg], color=c, lw=1.8, label=lab)
        # dotted exponential fit for the generic run
        rr = results["runs"][f"cg{cg}"]
        tt = np.linspace(0, rr["fit_window_t_end"], 50)
        ax.plot(tt, rr["H0_fit"] * np.exp(-tt / rr["tau"]), color="#1baf7a",
                lw=1.0, ls=":", alpha=0.9)
        floor = results["noise_floor"][f"cg{cg}"]
        ax.axhline(floor, color=muted, lw=1.0, ls="-.",
                   label=f"noise floor $\\approx$ {floor:.3g}")
        ax.axvline(T_FINAL, color=muted, lw=1.0, ls=":")
        ymin, ymax = ax.get_ylim()
        ax.annotate("$\\psi$ revival\n$t = 4\\pi$", xy=(T_FINAL, ymax),
                    xytext=(T_FINAL - 0.15, ymax * 0.55), ha="right",
                    fontsize=8, color=muted)
        ax.set_yscale("log")
        ax.set_xlabel("t", color=ink)
        ax.set_title(f"{cg}$\\times${cg} coarse-graining", fontsize=11,
                     color=ink)
        ax.grid(True, which="major", color="#DDDDDD", lw=0.6)
        ax.tick_params(colors=muted)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axes[0].set_ylabel(r"$\bar{H}(t)=\sum \bar{P}\,\ln(\bar{P}/\bar{Q})$",
                       color=ink)
    axes[0].legend(frameon=False, fontsize=9, loc="lower left")
    fig.suptitle("Does a non-degenerate spectrum restore relaxation at M = 4? "
                 f"(same $\\rho_0$ and seeds; N = {N_PART:,}; "
                 "dotted: exponential fit to generic run)",
                 fontsize=11.5, color=ink)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(os.path.join(OUTDIR, "m4generic_decay.png"), dpi=160)
    plt.close(fig)

    results["runtime_s"] = time.time() - t0
    with open(os.path.join(OUTDIR, "m4generic_results.json"), "w") as f:
        json.dump(results, f, indent=2)
    print(json.dumps(results, indent=2))
    print(f"TOTAL RUNTIME: {results['runtime_s']:.1f} s", flush=True)


if __name__ == "__main__":
    main()
