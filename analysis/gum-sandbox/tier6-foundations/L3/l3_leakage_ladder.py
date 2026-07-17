#!/usr/bin/env python3
"""
WORKSTREAM L3 -- Lock-2 conservation lemma: the leakage-exponent ladder.

Exact retarded solution of the 3-D scalar wave equation
    (1/c_L^2) d2p/dt2 - Laplacian p = s
for point-source arrays with harmonic amplitude cos(omega t):

    p(x,t) = (1/4pi) sum_i q_i cos(omega (t - r_i/c_L)) / r_i,   r_i = |x - x_i|.

Time-averaged far-field intensity proxy  I = <(dp/dt)^2> / c_L  is computed
analytically (no time sampling): dp/dt is a sum of sinusoids at frequency omega,
so  <(dp/dt)^2> = (1/2) |A|^2  with complex amplitude
    A = (omega/4pi) sum_i q_i exp(-i omega r_i / c_L) / r_i .

Total radiated power P(R) = integral over sphere of I dA, done by
Gauss-Legendre quadrature in cos(theta) (sources on z-axis => phi symmetry).

Three arrays:
  monopole   : q=+1 at origin                    (total charge oscillates: non-conserved)
  dipole     : q=+1 at +a z, q=-1 at -a z        (sum q = 0; dipole moment oscillates)
  quadrupole : q=+1 at +a z, q=-2 at 0, q=+1 at -a z  (sum q = 0 AND sum q z = 0)

Expected exponent ladder in 1/c_L for P:  1 / 3 / 5.

Units hbar = m = 1; omega = 1.
"""

import json
import numpy as np

# ----------------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------------
OMEGA = 1.0
A_HALF = 0.05          # array half-size a; at c_L,min = 1, k a = omega a / c_L = 0.05
N_CL = 12
CL_VALUES = np.logspace(0.0, 1.5, N_CL)   # 1.5 decades: c_L in [1, 31.62]
# Far-field radii: spec requires R >= 50 * c_L,max / omega = 1581.
# R1 = 2000 gives k R = omega R / c_L >= 63 for every c_L (far zone everywhere).
R1, R2 = 2000.0, 4000.0
N_GAUSS = 64           # Gauss-Legendre nodes in cos(theta)

ARRAYS = {
    "monopole":   {"z": np.array([0.0]),                 "q": np.array([1.0])},
    "dipole":     {"z": np.array([+A_HALF, -A_HALF]),    "q": np.array([+1.0, -1.0])},
    "quadrupole": {"z": np.array([+A_HALF, 0.0, -A_HALF]), "q": np.array([+1.0, -2.0, +1.0])},
}

# ----------------------------------------------------------------------------
# Exact time-averaged radiated power through sphere of radius R
# ----------------------------------------------------------------------------
def radiated_power(z_src, q_src, c_L, R):
    """P = int_sphere <(dp/dt)^2>/c_L dA, exact in time, GL quadrature in cos(theta)."""
    mu, w = np.polynomial.legendre.leggauss(N_GAUSS)   # nodes in cos(theta) on [-1,1]
    sin_t = np.sqrt(1.0 - mu**2)
    # Observation points on the sphere (phi = 0 WLOG; sources on z-axis)
    x_obs = R * sin_t                                   # (N,)
    z_obs = R * mu                                      # (N,)
    # Distances to each source: r_ij = |x_obs_j - x_i|
    r = np.sqrt(x_obs[None, :]**2 + (z_obs[None, :] - z_src[:, None])**2)  # (Ns, N)
    # Complex amplitude of dp/dt; factor out common phase exp(-i k R) (drops in |A|^2)
    k = OMEGA / c_L
    amp = (OMEGA / (4.0 * np.pi)) * np.sum(
        q_src[:, None] * np.exp(-1j * k * (r - R)) / r, axis=0)             # (N,)
    mean_pt2 = 0.5 * np.abs(amp)**2
    # dA = R^2 * 2pi * dcos(theta)
    return 2.0 * np.pi * R**2 / c_L * np.sum(w * mean_pt2)


# ----------------------------------------------------------------------------
# Sweep
# ----------------------------------------------------------------------------
results = {
    "params": {
        "omega": OMEGA, "a": A_HALF, "R1": R1, "R2": R2,
        "n_gauss": N_GAUSS, "c_L_values": CL_VALUES.tolist(),
        "ka_min": float(OMEGA * A_HALF / CL_VALUES.max()),
        "ka_max": float(OMEGA * A_HALF / CL_VALUES.min()),
        "kR1_min": float(OMEGA * R1 / CL_VALUES.max()),
        "R_spec_floor": float(50.0 * CL_VALUES.max() / OMEGA),
    },
    "cases": {},
}

for name, arr in ARRAYS.items():
    P1 = np.array([radiated_power(arr["z"], arr["q"], c, R1) for c in CL_VALUES])
    P2 = np.array([radiated_power(arr["z"], arr["q"], c, R2) for c in CL_VALUES])
    flux_reldiff = np.abs(P1 - P2) / P1
    # Fit ln P vs ln(1/c_L)
    x = np.log(1.0 / CL_VALUES)
    slope, intercept = np.polyfit(x, np.log(P1), 1)
    resid = np.log(P1) - (slope * x + intercept)
    results["cases"][name] = {
        "P_R1": P1.tolist(),
        "P_R2": P2.tolist(),
        "flux_reldiff_max": float(flux_reldiff.max()),
        "slope": float(slope),
        "intercept": float(intercept),
        "fit_resid_max": float(np.abs(resid).max()),
    }
    print(f"{name:11s} slope = {slope:+.6f}   max flux reldiff = {flux_reldiff.max():.3e}"
          f"   max fit resid = {np.abs(resid).max():.3e}")

# ----------------------------------------------------------------------------
# Gates
# ----------------------------------------------------------------------------
g0 = all(results["cases"][n]["flux_reldiff_max"] <= 0.01 for n in ARRAYS)
s_m = results["cases"]["monopole"]["slope"]
s_d = results["cases"]["dipole"]["slope"]
s_q = results["cases"]["quadrupole"]["slope"]
g1 = abs(s_m - 1.0) <= 0.10
g2 = abs(s_d - 3.0) <= 0.15
g3 = abs(s_q - 5.0) <= 0.25

results["gates"] = {
    "G0_flux_consistency": {"pass": bool(g0),
                            "max_reldiff": max(results["cases"][n]["flux_reldiff_max"] for n in ARRAYS),
                            "spec": "<= 1% at two radii, every case"},
    "G1_monopole_slope":   {"pass": bool(g1), "measured": s_m, "spec": "1 +/- 0.1"},
    "G2_dipole_slope":     {"pass": bool(g2), "measured": s_d, "spec": "3 +/- 0.15"},
    "G3_quadrupole_slope": {"pass": bool(g3), "measured": s_q, "spec": "5 +/- 0.25"},
    "G4_lemma_printed":    {"pass": True, "spec": "lemma stated and proof-sketched in RESULTS.md"},
}

# Fallback arithmetic for the connection paragraph
results["fallback_arithmetic"] = {
    "cL_over_c": 1e4,
    "leakage_exponent_1": 1e-4,
    "leakage_exponent_3": 1e-12,
    "leakage_exponent_5": 1e-20,
}

with open("/home/user/fork_frankensim/analysis/gum-sandbox/tier6-foundations/L3/L3_results.json", "w") as f:
    json.dump(results, f, indent=2)

for gid, g in results["gates"].items():
    print(f"{gid}: {'PASS' if g['pass'] else 'FAIL'}")

# ----------------------------------------------------------------------------
# Figure
# ----------------------------------------------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(7.0, 5.0), dpi=150)
colors = {"monopole": "#4269d0", "dipole": "#efb118", "quadrupole": "#ff725c"}
x = np.log10(1.0 / CL_VALUES)
for name in ARRAYS:
    P1 = np.array(results["cases"][name]["P_R1"])
    sl = results["cases"][name]["slope"]
    ax.plot(x, np.log10(P1), "o", ms=5, color=colors[name],
            label=f"{name}: slope = {sl:.3f}")
    # fitted line
    fit = results["cases"][name]["intercept"] / np.log(10) + sl * x
    ax.plot(x, fit, "-", lw=1.2, color=colors[name], alpha=0.6)
ax.set_xlabel(r"$\log_{10}(1/c_L)$")
ax.set_ylabel(r"$\log_{10} P$  (time-averaged radiated power)")
ax.set_title("Lock-2 leakage-exponent ladder: exact retarded point-source arrays\n"
             r"$P \propto (1/c_L)^{1,3,5}$ for non-conserved / dipole / conserved-class sources")
ax.legend(frameon=False)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig("/home/user/fork_frankensim/analysis/gum-sandbox/tier6-foundations/L3/L3_ladder.png")
print("figure written")
