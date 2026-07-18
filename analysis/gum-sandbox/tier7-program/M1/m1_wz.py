#!/usr/bin/env python3
"""
M1 -- WS4-M4c toy: w(z) from the Gamma(H) relaxation profile.

Tier 7 program extensions (ROADMAP_v9_PROGRAM.md, item M1; F-R12 anchor).
Within-model; nothing here bears on nature.  Units: H0 = 1 (hbar = m = 1
elsewhere in the campaign; neither enters).  Flat FRW test-field background:
E^2(x) = Om0 x^3 + (1-Om0), x = 1+z, Om0 = 0.3.

Model pieces
------------
(a) Algebraic tracker  rho_L = alpha H^2 M*^2            -> w = -1 + Om(z) exactly.
(b) Memory form        rho_L = alpha H^2 M*^2 + beta Hdot M*^2
    With Hdot = -(3/2) Om(z) H^2 (flat matter+Lambda), b = -(3/2) beta:
        rho_L = H^2 M*^2 (alpha + b Om),
        W(z) = (w+1)/Om = 1 + b (1-Om) / (alpha + b Om).
(c) Relaxation ODE     d(dq)/dt = -Gamma dq + S H,  Gamma = const (fast
    relaxation / tracking regime; prompt-binding choice).  rho_L = dq^2/2.
    Quasi-static attractor dq = (S/Gamma) H  =>  rho ~ H^2.
    Gradient expansion: dq = (S/Gamma)[H - Hdot/Gamma + ...] =>
    effective beta/alpha = -2 H*/Gamma < 0 (sign forced).
(d) Gamma(H) profile family: Gamma = g0 E^n  =>  dq_qs ~ H^(1-n)
    =>  w = -1 + (1-n) Om(z),  i.e.  W = 1-n.

All w from the continuity equation: w = -1 + (1/3) dln rho / dln(1+z).

Gates G1-G5 measured below; results to m1_results.json.
"""

import json
import numpy as np
from scipy.integrate import solve_ivp, cumulative_trapezoid

OM0 = 0.3
OL0 = 1.0 - OM0

# ----------------------------------------------------------------- background
def E2(x):
    return OM0 * x**3 + OL0

def Om(x):
    return OM0 * x**3 / E2(x)

def E2_a(a):
    return E2(1.0 / a)

def Om_a(a):
    return Om(1.0 / a)

# 4th-order central first derivative on a uniform grid
def d4(y, h):
    d = np.empty_like(y)
    d[2:-2] = (-y[4:] + 8*y[3:-1] - 8*y[1:-3] + y[:-4]) / (12*h)
    # 2nd-order one-sided at edges (edges excluded from gate windows)
    d[0] = (-3*y[0] + 4*y[1] - y[2]) / (2*h)
    d[1] = (y[2] - y[0]) / (2*h)
    d[-2] = (y[-1] - y[-3]) / (2*h)
    d[-1] = (3*y[-1] - 4*y[-2] + y[-3]) / (2*h)
    return d

results = {}

# ============================================================ G2: exact tracker
# Analytic identity (two lines, printed in RESULTS.md):
#   (1) H^2 = H0^2 [Om0 x^3 + OL0]  =>  dln H^2/dln x = 3 Om0 x^3 / E^2 = 3 Om(z).
#   (2) rho ~ H^2, continuity w = -1 + (1/3) dln rho/dln x  =>  w = -1 + Om(z).
# Numerical check: 4th-order finite differences of ln rho on a ln x grid.
u_lo, u_hi = -0.05, np.log(11.0) + 0.05     # pad so z in [0,10] is interior
N_grid = 8001
u = np.linspace(u_lo, u_hi, N_grid)
h_u = u[1] - u[0]
x_grid = np.exp(u)
lnrho_tracker = np.log(E2(x_grid))          # rho ~ H^2 (alpha, M*^2 cancel)
w_fd = -1.0 + d4(lnrho_tracker, h_u) / 3.0
w_exact = -1.0 + Om(x_grid)
mask_z = (x_grid >= 1.0) & (x_grid <= 11.0)
g2_maxerr = float(np.max(np.abs(w_fd[mask_z] - w_exact[mask_z])))
results["G2_tracker_identity"] = {
    "grid_points": N_grid,
    "z_range": [0.0, 10.0],
    "max_abs_dev_w_vs_identity": g2_maxerr,
    "spec": 1e-6,
}

# ==================================================== memory forms: W(z) closed
def W_memory(x, beta, alpha=1.0):
    b = -1.5 * beta
    om = Om(x)
    return 1.0 + b * (1.0 - om) / (alpha + b * om)

def w_memory_a(a, beta, alpha=1.0):
    x = 1.0 / a
    return -1.0 + W_memory(x, beta, alpha) * Om(x)

# FD cross-check of the closed form (memory rho = E^2*(alpha + b*Om))
def check_memory_fd(beta):
    b = -1.5 * beta
    lnrho = np.log(E2(x_grid) * (1.0 + b * Om(x_grid)))
    w_num = -1.0 + d4(lnrho, h_u) / 3.0
    w_cf = -1.0 + W_memory(x_grid, beta) * Om(x_grid)
    return float(np.max(np.abs(w_num[mask_z] - w_cf[mask_z])))

mem_fd_err = {str(b): check_memory_fd(b) for b in (-0.3, -1.0)}

# ------------------------------------------------------- G4 convention testing
printed = {"-0.3": (1.28, 1.03), "-1.0": (1.72, 1.05)}
g4 = {"printed_pairs": printed, "memory_closedform_fd_crosscheck": mem_fd_err}

def Wz(beta, z):
    return float(W_memory(np.array([1.0 + z]), beta)[0])

# Reading 1 (task-spec): pair = (W(z=0.5), W(z=0))
r1 = {}
for b in (-0.3, -1.0):
    W05, W0 = Wz(b, 0.5), Wz(b, 0.0)
    p = printed[str(b)]
    r1[str(b)] = {
        "W(z=0.5)": W05, "W(z=0)": W0,
        "dev_first_pct": 100 * abs(W05 - p[0]) / p[0],
        "dev_second_pct": 100 * abs(W0 - p[1]) / p[1],
    }
g4["reading1_spec_z05_to_z0"] = r1

# Reading 2 (recovered): pair = (W(z=0), W(z=2))
r2 = {}
for b in (-0.3, -1.0):
    W0, W2 = Wz(b, 0.0), Wz(b, 2.0)
    p = printed[str(b)]
    r2[str(b)] = {
        "W(z=0)": W0, "W(z=2)": W2,
        "dev_first_pct": 100 * abs(W0 - p[0]) / p[0],
        "dev_second_pct": 100 * abs(W2 - p[1]) / p[1],
    }
g4["reading2_recovered_z0_to_z2"] = r2

# Scan: which z2 best matches the second printed number (uniqueness of z=2)?
zs = np.linspace(0.0, 6.0, 6001)
best_z = {}
for b in (-0.3, -1.0):
    Wv = W_memory(1.0 + zs, b)
    tgt = printed[str(b)][1]
    best_z[str(b)] = float(zs[np.argmin(np.abs(Wv - tgt))])
g4["z_where_W_equals_second_printed"] = best_z
results["G4_memory_recovery"] = g4

# ================================================================ G3: CPL fits
def cpl_fit(w_of_a, a_lo, a_hi, n=400):
    a = np.linspace(a_lo, a_hi, n)
    wv = w_of_a(a)
    A = np.vstack([np.ones_like(a), 1.0 - a]).T
    coef, *_ = np.linalg.lstsq(A, wv, rcond=None)
    return float(coef[0]), float(coef[1])          # w0, wa

def wa_tangent(w_of_a, h=1e-4):
    """wa = -dw/da at a=1 via a 5-point one-sided stencil.  Samples run along
    decreasing a (a_k = 1 - k h); the standard forward stencil for g'(0) is
    (-25g0+48g1-36g2+16g3-3g4)/12 = -h dw/da, so wa = -dw/da = that / h.
    Sanity-checked against the analytic pure-tracker value 3 Om0 (1-Om0)."""
    a = 1.0 - h * np.arange(5)
    g = w_of_a(a)
    return float((-25*g[0] + 48*g[1] - 36*g[2] + 16*g[3] - 3*g[4]) / (12*h))

def w_pure(a, wcal=1.0):
    return -1.0 + wcal * Om_a(a)

conventions = {
    "A_tangent_a1": ("derivative at a=1 (corpus formula's own definition)", None),
    "B_window_0.6_0.8": ("LS fit a in [0.6,0.8] (literal 'around a=0.7')", (0.6, 0.8)),
    "C_window_0.5_1.0": ("LS fit a in [0.5,1.0] (z in [0,1])", (0.5, 1.0)),
    "D_window_0.7_1.0": ("LS fit a in [0.7,1.0]", (0.7, 1.0)),
}

g3 = {"corpus_formula_wa_at_Wcal1": 3 * OM0 * (1 - OM0)}
tab = {}
for key, (desc, win) in conventions.items():
    if win is None:
        wa = wa_tangent(lambda a: w_pure(a))
    else:
        _, wa = cpl_fit(lambda a: w_pure(a), *win)
    tab[key] = {"description": desc, "wa": wa,
                "dev_from_0.63_pct": 100 * abs(wa - 0.63) / 0.63}
g3["pure_tracker_conventions"] = tab

# family sign scan (algebraic members; dynamics members appended later)
family = []
for wcal in (0.25, 0.5, 1.0, 1.5):
    family.append((f"tracker_Wcal={wcal}", lambda a, c=wcal: w_pure(a, c)))
for b in (-0.1, -0.3, -0.5, -1.0, -2.0):
    family.append((f"memory_beta={b}", lambda a, bb=b: w_memory_a(a, bb)))
for n_prof in (-0.5, 0.0, 0.5):
    family.append((f"GammaProfile_n={n_prof}",
                   lambda a, nn=n_prof: -1.0 + (1 - nn) * Om_a(a)))

sign_scan = {}
for name, wf in family:
    wa_A = wa_tangent(wf)
    _, wa_C = cpl_fit(wf, 0.5, 1.0)
    sign_scan[name] = {"wa_tangent": wa_A, "wa_window_C": wa_C}
g3["family_sign_scan"] = sign_scan

# memory-form tangent wa (for the corpus's "+0.36 at beta=-0.3" print)
g3["memory_wa"] = {
    str(b): {"wa_tangent": wa_tangent(lambda a, bb=b: w_memory_a(a, bb)),
             "wa_window_C": cpl_fit(lambda a, bb=b: w_memory_a(a, bb), 0.5, 1.0)[1]}
    for b in (-0.3, -1.0)
}
results["G3_CPL"] = g3

# ======================================================= dynamics leg (ODE, G1)
# d(dq)/dN = S - (Gamma/E) dq,  N = ln a, integrate z=3000 -> 0 with dq(ini)=0.
S = 1.0
Z_START = 3000.0
N_i = -np.log(1.0 + Z_START)

def run_dq(Gamma, n_prof=0.0, rtol=1e-12, atol=1e-14, dq0=0.0, N_start=None,
           method="DOP853", t_eval=None):
    if N_start is None:
        N_start = N_i
    def rhs(N, y):
        a = np.exp(N)
        E = np.sqrt(E2_a(a))
        return [S - (Gamma * E**n_prof / E) * y[0]]
    sol = solve_ivp(rhs, (N_start, 0.0), [dq0], method=method,
                    rtol=rtol, atol=atol, dense_output=True, t_eval=t_eval)
    assert sol.success
    return sol

def w_dyn_of(sol, Gamma, n_prof=0.0):
    """w on a grid from the ODE solution, using the exact RHS (no FD noise):
       w = -1 - (2/3) (d dq/dN)/dq."""
    def wf(a):
        N = np.log(a)
        dq = sol.sol(N)[0]
        E = np.sqrt(E2_a(a))
        ddq = S - (Gamma * E**n_prof / E) * dq
        return -1.0 - (2.0 / 3.0) * ddq / dq
    return wf

dyn = {}
GAMMAS = (50.0, 200.0, 800.0)
z_eval = np.linspace(0.0, 10.0, 2001)
a_eval = 1.0 / (1.0 + z_eval)
dev0, dev10 = [], []
sols = {}
for G in GAMMAS:
    sol = run_dq(G)
    sols[G] = sol
    wf = w_dyn_of(sol, G)
    wv = wf(a_eval)
    omv = Om_a(a_eval)
    Wdyn = (wv + 1.0) / omv
    rho = 0.5 * sol.sol(np.log(a_eval))[0] ** 2
    dyn[str(G)] = {
        "W_dyn(z=0)": float(Wdyn[0]),
        "W_dyn(z=0.5)": float(Wdyn[np.argmin(np.abs(z_eval - 0.5))]),
        "min_W_dyn_z0_10": float(np.min(Wdyn)),
        "max_abs_w_dev_from_identity_z0_10": float(np.max(np.abs(wv - (-1 + omv)))),
        "w_dev_at_z0": float(wv[0] - (-1 + omv[0])),
        "w_dev_at_z10": float(wv[-1] - (-1 + omv[-1])),
        "w_dev_at_z3": float(wv[np.argmin(np.abs(z_eval - 3.0))]
                             - (-1 + omv[np.argmin(np.abs(z_eval - 3.0))])),
        "z_of_max_dev": float(z_eval[np.argmax(np.abs(wv - (-1 + omv)))]),
        "z_tracking_entry_Gamma_eq_H": float(((G**2 - OL0) / OM0) ** (1/3) - 1),
        "min_rho": float(np.min(rho)),
        "beta_eff_over_alpha_at_z0_predicted": -2.0 / G,
    }
    dev0.append(abs(wv[0] - (-1 + omv[0])))
    dev10.append(abs(wv[-1] - (-1 + omv[-1])))
    # append to the G3 family sign scan
    wa_A = wa_tangent(wf)
    _, wa_C = cpl_fit(wf, 0.5, 1.0)
    g3["family_sign_scan"][f"dynamics_Gamma={G}"] = {
        "wa_tangent": wa_A, "wa_window_C": wa_C}

lg = np.log(np.array(GAMMAS))
dyn["lag_scaling_logslope_z0"] = float(np.polyfit(lg, np.log(dev0), 1)[0])
dyn["lag_scaling_logslope_z10"] = float(np.polyfit(lg, np.log(dev10), 1)[0])

# first-order lag prediction (gradient expansion, derived in RESULTS.md):
#   dq = (S/Gamma)[H - Hdot/Gamma + O(Gamma^-2)]
#   =>  W(z) - 1 = 3 (H/Gamma) (1 - Om/2)  =>  W(0) - 1 = (3 - 1.5*Om0)/Gamma.
for G in GAMMAS:
    pred = (3.0 - 1.5 * OM0) / G
    meas = dyn[str(G)]["W_dyn(z=0)"] - 1.0
    dyn[str(G)]["W0_minus1_predicted_first_order"] = pred
    dyn[str(G)]["W0_minus1_measured"] = meas
    dyn[str(G)]["W0_meas_over_pred"] = meas / pred

# n=0.5 profile dynamics check: quasi-static Wcal = 1-n = 0.5
sol_n = run_dq(400.0, n_prof=0.5)
wf_n = w_dyn_of(sol_n, 400.0, n_prof=0.5)
wv_n = wf_n(a_eval)
Wcal_n_meas = float((wv_n[0] + 1.0) / Om_a(1.0))
dyn["profile_n0.5_Wcal_measured_z0"] = Wcal_n_meas
dyn["profile_n0.5_Wcal_predicted"] = 0.5

# attractor approach rate (linear ODE => difference of two solutions decays at
# exactly Gamma): Gamma=50, displace at z=10 by x1.5, fit ln|Delta| vs t.
G_rate = 50.0
N10 = -np.log(11.0)
base = sols[G_rate]
dq10 = base.sol(N10)[0]
pert = run_dq(G_rate, dq0=1.5 * dq10, N_start=N10)
N_f = np.linspace(N10, 0.0, 4001)
a_f = np.exp(N_f)
E_f = np.sqrt(E2_a(a_f))
t_f = cumulative_trapezoid(1.0 / E_f, N_f, initial=0.0)   # cosmic time (H0=1)
delta = np.abs(pert.sol(N_f)[0] - base.sol(N_f)[0])
rel = delta / np.abs(base.sol(N_f)[0])
mfit = (G_rate * t_f > 1.0) & (G_rate * t_f < 15.0) & (rel > 1e-11)
rate = -np.polyfit(t_f[mfit], np.log(delta[mfit]), 1)[0]
dyn["approach_rate"] = {
    "Gamma": G_rate, "measured_rate": float(rate),
    "rate_over_Gamma_minus_1": float(rate / G_rate - 1.0),
    "fit_window_Gamma_t": [1.0, 15.0],
}
results["dynamics"] = dyn

# --------------------------------------------------- G1 numerics + positivity
# (i) integrator convergence: DOP853 rtol 1e-12 vs Radau rtol 1e-10 reference
sol_a = sols[200.0]
sol_b = run_dq(200.0, rtol=1e-10, atol=1e-13, method="Radau")
N_chk = np.linspace(N_i * 0.999, 0.0, 512)
qa, qb = sol_a.sol(N_chk)[0], sol_b.sol(N_chk)[0]
ode_reldiff = float(np.max(np.abs(qa - qb) / np.maximum(np.abs(qa), 1e-300)))
# (ii) FD-vs-analytic w on the algebraic grid (from G2)
# (iii) positivity across every form
min_rho_dyn = min(dyn[str(G)]["min_rho"] for G in GAMMAS)
b_grid = np.array([-0.1, -0.3, -0.5, -1.0, -2.0])
min_f_mem = float(min(np.min(1.0 + (-1.5 * b) * Om(x_grid[mask_z])) for b in b_grid))
results["G1_numerics_positivity"] = {
    "ode_cross_integrator_reldiff": ode_reldiff,
    "fd_vs_analytic_w_maxerr": g2_maxerr,
    "min_rho_dynamics_over_forms": min_rho_dyn,
    "min_memory_shape_factor_1_plus_bOm": min_f_mem,
    "tracker_rho_positive": True,
    "spec": 1e-8,
}

# ===================================================== G5: Branch-A operational
# tracker limit: (w+1)/Om = 1 exactly (same measurement as G2, restated as W)
W_fd = (w_fd[mask_z] + 1.0) / Om(x_grid[mask_z])
g5 = {"tracker_max_abs_W_minus_1": float(np.max(np.abs(W_fd - 1.0)))}
# deviation shape with memory: (w+1) - Om = b Om(1-Om)/(1+b Om)
for b in (-0.3, -1.0):
    bb = -1.5 * b
    zz = np.linspace(0, 10, 2001)
    omz = Om(1.0 + zz)
    devs = bb * omz * (1 - omz) / (1 + bb * omz)
    g5[f"memory_beta={b}"] = {
        "deviation_form": "(w+1)-Om = b*Om*(1-Om)/(1+b*Om), b=-1.5*beta",
        "peak_z": float(zz[np.argmax(devs)]),
        "peak_value": float(np.max(devs)),
        "dev_at_z0": float(devs[0]),
        "dev_at_z10": float(devs[-1]),
        "W_monotone_decreasing_in_z": bool(np.all(np.diff(
            1 + bb * (1 - omz) / (1 + bb * omz)) <= 0)),
        "W_ge_1_everywhere": bool(np.all(1 + bb * (1 - omz) / (1 + bb * omz) >= 1.0)),
    }
results["G5_branchA"] = g5

# ============================================= self-consistency (one iteration)
# rho_de,0(x) = 0.7 * E0^2 * f(x)/f(1);  E1^2 = 0.3 x^3 + rho_de,0;
# recompute rho_de,1 from E1 (same functional form), renormalize, re-derive w.
def self_consistency(beta):
    b = -1.5 * beta
    E20 = E2(x_grid)
    om0 = Om(x_grid)
    f0 = (1.0 + b * om0)
    rho0 = OL0 * E20 * f0 / (1.0 + b * OM0)      # normalized: rho_de(z=0) = 0.7
    E21 = OM0 * x_grid**3 + rho0                 # one Friedmann iteration
    # rho_de,1 from E1: alpha H^2 + beta Hdot with Hdot = -(1/2) dE^2/dlnx
    dE21 = d4(E21, h_u)
    rho1 = E21 + beta * (-0.5 * dE21)
    i0 = np.argmin(np.abs(x_grid - 1.0))
    rho1 = OL0 * rho1 / rho1[i0]
    w0 = -1.0 + d4(np.log(np.maximum(rho0, 1e-300)), h_u) / 3.0
    w1 = -1.0 + d4(np.log(np.maximum(rho1, 1e-300)), h_u) / 3.0
    i05 = np.argmin(np.abs(x_grid - 1.5))
    return {
        "w_iter0_z0": float(w0[i0]), "w_iter1_z0": float(w1[i0]),
        "w_iter0_z05": float(w0[i05]), "w_iter1_z05": float(w1[i05]),
        "delta_w_z0": float(w1[i0] - w0[i0]),
        "delta_w_z05": float(w1[i05] - w0[i05]),
    }

sc = {"pure_tracker": self_consistency(0.0),
      "beta=-0.3": self_consistency(-0.3),
      "beta=-1.0": self_consistency(-1.0),
      "note": ("pure alpha H^2 has NO self-consistent accelerating fixed point "
               "(iteration flows toward EdS; the corpus's own printed pathology); "
               "the test-field/LCDM-background reading is what the identity is "
               "a statement about")}
results["self_consistency_one_iteration"] = sc

# ------------------------------------------------------------------ save JSON
with open("/home/user/fork_frankensim/analysis/gum-sandbox/tier7-program/M1/"
          "m1_results.json", "w") as f:
    json.dump(results, f, indent=2)

# ------------------------------------------------------------------- printout
print("== G2 tracker identity: max|w-(-1+Om)| =", g2_maxerr)
print("== G4 reading1 (z=0.5,0):", json.dumps(r1, indent=1))
print("== G4 reading2 (z=0,2):", json.dumps(r2, indent=1))
print("== G4 best-z for 2nd printed:", best_z)
print("== G3 conventions:", json.dumps(tab, indent=1))
print("== G3 memory wa:", json.dumps(g3["memory_wa"], indent=1))
print("== G3 min wa over family:",
      min(min(v["wa_tangent"], v["wa_window_C"])
          for v in g3["family_sign_scan"].values()))
print("== dynamics:", json.dumps(dyn, indent=1))
print("== G1:", json.dumps(results["G1_numerics_positivity"], indent=1))
print("== G5:", json.dumps(g5, indent=1))
print("== self-consistency:", json.dumps(sc, indent=1))
