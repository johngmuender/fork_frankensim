#!/usr/bin/env python3
"""
L7a - POVM-exclusion testbed, quadraticity/polarization test (1-D spatial family, exact).

Tier 6 (foundations), Phase L7, workstream L7a of the GUM replication campaign.
File as F-T6-L7a.

Units: hbar = m = 1.

Model
-----
Free Gaussian components (closed form, no PDE solve):
    psi_j(x,0) = (pi sigma^2)^{-1/4} exp(-(x-x0)^2/(2 sigma^2)) exp(i k_j x)
evolves under H = p^2/2 to
    psi_j(x,t) = (pi sigma^2)^{-1/4} sqrt(sigma^2/w) *
                 exp( -(x - x0 - k_j t)^2 / (2 w) + i k_j x - i k_j^2 t / 2 ),
    w = sigma^2 + i t.
(Derived by Fourier transform of the initial Gaussian and exact quadratic
integration of the free propagator; verified below by norm conservation.)

Family: psi_{theta,phi} = [cos(theta) psi_1 + sin(theta) e^{i phi} psi_2] / n,
        n^2 = 1 + sin(2 theta) Re(e^{i phi} <psi_1|psi_2>).
The overlap <psi_1|psi_2> = exp(-sigma^2 (k2-k1)^2 / 4) e^{i (k2-k1) x0} is NOT
zero, so each member is normalized explicitly (the POVM form for unnormalized
psi is a ratio of quadratics; normalizing restores an exact quadratic-form test:
data(theta,phi) * n^2 must be a 4-parameter sesquilinear form on span{psi_1,psi_2}).

Exact Bohmian first-arrival CDF (1-D no-crossing theorem):
    CDF_psi(t) = max_{s<=t} N_psi(s),   N_psi(s) = int_d^inf |psi_s|^2 dx  (d=0).

POVM (quadratic-form) fit per time bin:
    F(theta,phi) = a cos^2(theta) + b sin^2(theta)
                 + sin(2 theta) (c cos(phi) + d sin(phi)),
fitted to n^2 * data by linear least squares (equivalently: basis functions
divided by n^2 fitted to data directly - we do the latter).
"""

import json
import time as walltime

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import erfc

# ----------------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------------
SIG = 1.0
X0 = -8.0
K1, K2 = 1.2, 2.4
D = 0.0                      # detector position
T_HOR = 12.0                 # horizon
DT = 0.005                   # time grid step (<= 0.005 per spec)
NT = int(round(T_HOR / DT)) + 1
TGRID = np.linspace(0.0, T_HOR, NT)

NBIN = 48                    # uniform time bins over [0, T]
BIN_EDGES = np.linspace(0.0, T_HOR, NBIN + 1)
EDGE_IDX = np.round(BIN_EDGES / DT).astype(int)
assert np.allclose(TGRID[EDGE_IDX], BIN_EDGES)

THETAS = np.arange(7) * np.pi / 12.0          # 0 .. pi/2, 7 values
PHIS = np.arange(8) * np.pi / 4.0             # 0 .. 7pi/4, 8 values

NTRAJ = 4000
DT_TRAJ = 0.002

RNG_DIR = "/home/user/fork_frankensim/analysis/gum-sandbox/tier6-foundations/L7a"

# ----------------------------------------------------------------------------
# Closed-form free Gaussian component and derivative
# ----------------------------------------------------------------------------
def psi_comp(x, t, k, x0, sig=SIG):
    """Normalized free-evolved boosted Gaussian, closed form."""
    w = sig**2 + 1j * t
    pref = (np.pi * sig**2) ** (-0.25) * np.sqrt(sig**2 / w)
    return pref * np.exp(-(x - x0 - k * t) ** 2 / (2.0 * w)
                         + 1j * k * x - 0.5j * k * k * t)


def logderiv_comp(x, t, k, x0, sig=SIG):
    """d/dx log psi_j, closed form."""
    w = sig**2 + 1j * t
    return -(x - x0 - k * t) / w + 1j * k


# ----------------------------------------------------------------------------
# Quadrature: composite Gauss-Legendre panels (spectrally accurate per panel)
# ----------------------------------------------------------------------------
def panel_nodes(a, b, h=0.5, npt=12):
    m = int(round((b - a) / h))
    xg, wg = leggauss(npt)
    edges = np.linspace(a, b, m + 1)
    mid = 0.5 * (edges[:-1] + edges[1:])
    half = 0.5 * (edges[1:] - edges[:-1])
    nodes = (mid[:, None] + half[:, None] * xg[None, :]).ravel()
    weights = (half[:, None] * wg[None, :]).ravel()
    return nodes, weights


XR, WR = panel_nodes(0.0, 160.0)      # right half-line (detector side), truncated
XL, WL = panel_nodes(-160.0, 0.0)     # left half-line


def pair_integrals(k_a, x0_a, k_b, x0_b):
    """
    For components a, b return, on TGRID:
      right-half integrals Iaa_R, Ibb_R, Iab_R  (Iab = int psi_a^* psi_b),
      full-line integrals  Iaa_F, Ibb_F, Iab_F  (norm / overlap checks),
    computed by the panel Gauss-Legendre rule.
    """
    Iaa_R = np.empty(NT); Ibb_R = np.empty(NT)
    Iab_R = np.empty(NT, complex)
    Iaa_F = np.empty(NT); Ibb_F = np.empty(NT)
    Iab_F = np.empty(NT, complex)
    for i, t in enumerate(TGRID):
        paR = psi_comp(XR, t, k_a, x0_a); pbR = psi_comp(XR, t, k_b, x0_b)
        paL = psi_comp(XL, t, k_a, x0_a); pbL = psi_comp(XL, t, k_b, x0_b)
        aaR = np.dot(WR, np.abs(paR) ** 2); bbR = np.dot(WR, np.abs(pbR) ** 2)
        abR = np.dot(WR, np.conj(paR) * pbR)
        aaL = np.dot(WL, np.abs(paL) ** 2); bbL = np.dot(WL, np.abs(pbL) ** 2)
        abL = np.dot(WL, np.conj(paL) * pbL)
        Iaa_R[i] = aaR; Ibb_R[i] = bbR; Iab_R[i] = abR
        Iaa_F[i] = aaR + aaL; Ibb_F[i] = bbR + bbL; Iab_F[i] = abR + abL
    return Iaa_R, Ibb_R, Iab_R, Iaa_F, Ibb_F, Iab_F


# ----------------------------------------------------------------------------
# Family machinery
# ----------------------------------------------------------------------------
def family_grid():
    TH, PH = np.meshgrid(THETAS, PHIS, indexing="ij")
    return TH.ravel(), PH.ravel()          # 56 preparations


def family_N(Iaa_R, Ibb_R, Iab_R, overlap):
    """N_psi(t) for every (theta,phi), normalized states. Shape (56, NT)."""
    th, ph = family_grid()
    c2 = np.cos(th) ** 2
    s2 = np.sin(th) ** 2
    s2t = np.sin(2.0 * th)
    n2 = 1.0 + s2t * (np.cos(ph) * overlap.real - np.sin(ph) * overlap.imag)
    quad = (c2[:, None] * Iaa_R[None, :]
            + s2[:, None] * Ibb_R[None, :]
            + s2t[:, None] * (np.cos(ph)[:, None] * Iab_R.real[None, :]
                              - np.sin(ph)[:, None] * Iab_R.imag[None, :]))
    return quad / n2[:, None], n2


def design_matrix(overlap):
    """56 x 4 design for the POVM form on NORMALIZED states (basis / n^2)."""
    th, ph = family_grid()
    n2 = 1.0 + np.sin(2 * th) * (np.cos(ph) * overlap.real
                                 - np.sin(ph) * overlap.imag)
    B = np.column_stack([np.cos(th) ** 2,
                         np.sin(th) ** 2,
                         np.sin(2 * th) * np.cos(ph),
                         np.sin(2 * th) * np.sin(ph)])
    return B / n2[:, None]


def fit_povm_form(data_mat, design):
    """
    data_mat: (56, nbins). Least-squares fit of the POVM form per bin.
    Returns per-bin max abs residual (nbins,), residual matrix, coefs.
    """
    coef, *_ = np.linalg.lstsq(design, data_mat, rcond=None)
    resid = data_mat - design @ coef
    return np.max(np.abs(resid), axis=0), resid, coef


def cdf_and_bins(Nmat):
    """Running-max CDF and per-bin arrival probabilities (+ non-arrival)."""
    CDF = np.maximum.accumulate(Nmat, axis=1)
    edges = CDF[:, EDGE_IDX]                              # (56, 49)
    incr = np.diff(edges, axis=1)                         # (56, 48)
    nonarr = 1.0 - CDF[:, -1]                             # (56,)
    return CDF, np.column_stack([incr, nonarr])           # (56, 49)


# ----------------------------------------------------------------------------
# Bohmian trajectory fan (G1)
# ----------------------------------------------------------------------------
def velocity(x, t, c1, c2, par1, par2):
    p1 = psi_comp(x, t, *par1); p2 = psi_comp(x, t, *par2)
    l1 = logderiv_comp(x, t, *par1); l2 = logderiv_comp(x, t, *par2)
    num = c1 * p1 * l1 + c2 * p2 * l2
    den = c1 * p1 + c2 * p2
    return np.imag(num / den)


def run_fan(theta, phi, overlap, par1, par2, n_traj=NTRAJ, dt=DT_TRAJ):
    n2 = 1.0 + np.sin(2 * theta) * (np.cos(phi) * overlap.real
                                    - np.sin(phi) * overlap.imag)
    c1 = np.cos(theta) / np.sqrt(n2)
    c2 = np.sin(theta) * np.exp(1j * phi) / np.sqrt(n2)

    # quantile sampling from |psi_0|^2
    xg = np.linspace(X0 - 6.5, X0 + 6.5, 260001)
    rho = np.abs(c1 * psi_comp(xg, 0.0, *par1)
                 + c2 * psi_comp(xg, 0.0, *par2)) ** 2
    cdf0 = np.concatenate([[0.0], np.cumsum(0.5 * (rho[1:] + rho[:-1])
                                            * np.diff(xg))])
    cdf0 /= cdf0[-1]
    q = (np.arange(n_traj) + 0.5) / n_traj
    x = np.interp(q, cdf0, xg)

    nstep = int(round(T_HOR / dt))
    crossed = np.zeros(n_traj, bool)
    tcross = np.full(n_traj, np.inf)
    min_gap = np.inf
    n_inversion_steps = 0
    for istep in range(nstep):
        t = istep * dt
        x_old = x
        k1v = velocity(x, t, c1, c2, par1, par2)
        k2v = velocity(x + 0.5 * dt * k1v, t + 0.5 * dt, c1, c2, par1, par2)
        k3v = velocity(x + 0.5 * dt * k2v, t + 0.5 * dt, c1, c2, par1, par2)
        k4v = velocity(x + dt * k3v, t + dt, c1, c2, par1, par2)
        x = x + dt / 6.0 * (k1v + 2 * k2v + 2 * k3v + k4v)
        gap = np.min(np.diff(x))
        min_gap = min(min_gap, gap)
        if gap < 0:
            n_inversion_steps += 1
        newly = (~crossed) & (x >= D)
        frac = np.clip((D - x_old[newly]) / (x[newly] - x_old[newly]), 0, 1)
        tcross[newly] = t + frac * dt
        crossed |= newly
    tc_sorted = np.sort(tcross[np.isfinite(tcross)])
    ecdf = np.searchsorted(tc_sorted, TGRID, side="right") / n_traj
    return ecdf, min_gap, n_inversion_steps, crossed.sum()


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    t0 = walltime.time()
    res = {}

    # ---- components, overlap ----
    par1, par2 = (K1, X0), (K2, X0)
    ov_analytic = np.exp(-SIG**2 * (K2 - K1) ** 2 / 4.0) \
        * np.exp(1j * (K2 - K1) * X0)
    Iaa_R, Ibb_R, Iab_R, Iaa_F, Ibb_F, Iab_F = pair_integrals(*par1, *par2)

    norm_dev = max(np.max(np.abs(Iaa_F - 1.0)), np.max(np.abs(Ibb_F - 1.0)))
    ov_dev_analytic = np.max(np.abs(Iab_F - ov_analytic))
    ov_dev_const = np.max(np.abs(Iab_F - Iab_F[0]))
    overlap = Iab_F[0]

    # right-integral vs 1 - left-integral consistency = norm closure above
    res["overlap_analytic"] = [ov_analytic.real, ov_analytic.imag]
    res["overlap_numeric_t0"] = [overlap.real, overlap.imag]
    res["overlap_abs"] = float(np.abs(overlap))
    res["overlap_max_dev_from_analytic"] = float(ov_dev_analytic)
    res["overlap_max_dev_from_const_in_t"] = float(ov_dev_const)
    res["component_norm_max_dev"] = float(norm_dev)

    # ---- family presence probability, CDF, bins ----
    Nmat, n2_fam = family_N(Iaa_R, Ibb_R, Iab_R, overlap)
    CDF, bins_bohm = cdf_and_bins(Nmat)
    backflow_per_prep = np.max(np.maximum.accumulate(Nmat, axis=1) - Nmat,
                               axis=1)
    star = int(np.argmax(backflow_per_prep))
    th_f, ph_f = family_grid()
    res["max_backflow_over_family"] = float(backflow_per_prep[star])
    res["star_prep"] = {"index": star, "theta": float(th_f[star]),
                        "phi": float(ph_f[star]),
                        "theta_over_pi": float(th_f[star] / np.pi),
                        "phi_over_pi": float(ph_f[star] / np.pi),
                        "backflow": float(backflow_per_prep[star])}
    res["n_preps"] = len(th_f)
    res["n2_family_min_max"] = [float(n2_fam.min()), float(n2_fam.max())]

    # backflow-carrying bins: CDF increment differs from N increment
    N_edges = Nmat[:, EDGE_IDX]
    N_incr = np.diff(N_edges, axis=1)
    C_incr = bins_bohm[:, :NBIN]
    bin_backflow_dev = np.max(np.abs(C_incr - N_incr), axis=0)   # (48,)
    backflow_bins = np.where(bin_backflow_dev > 1e-10)[0]
    res["bin_backflow_dev"] = bin_backflow_dev.tolist()
    res["backflow_bins"] = backflow_bins.tolist()

    # ---- POVM-form fits ----
    design = design_matrix(overlap)
    # G2 comparator floor: presence probability N at the 49 bin edges
    resid_floor_perbin, _, _ = fit_povm_form(N_edges, design)
    floor = float(np.max(resid_floor_perbin))
    res["G2_floor_perbin_max"] = resid_floor_perbin.tolist()
    res["G2_floor"] = floor

    # Bohmian per-bin fit (48 increments + non-arrival)
    resid_bohm_perbin, resid_bohm, _ = fit_povm_form(bins_bohm, design)
    res["bohm_resid_perbin"] = resid_bohm_perbin.tolist()
    res["bohm_resid_global_max"] = float(np.max(resid_bohm_perbin))

    # G3 evaluation
    g3_thresh = max(1e3 * floor, 1e-4)
    g3_bins = [int(i) for i in backflow_bins
               if resid_bohm_perbin[i] > g3_thresh]
    res["G3_threshold"] = g3_thresh
    res["G3_passing_backflow_bins"] = g3_bins
    res["G3_max_resid_in_backflow_bins"] = float(
        np.max(resid_bohm_perbin[backflow_bins])) if len(backflow_bins) else 0.0

    # ---- G1: trajectory fan at the star preparation ----
    ecdf, min_gap, n_inv, n_crossed = run_fan(th_f[star], ph_f[star],
                                              overlap, par1, par2)
    dkw = np.sqrt(np.log(2.0 / 0.05) / (2.0 * NTRAJ))
    sup = float(np.max(np.abs(ecdf - CDF[star])))
    res["G1_fan"] = {"n_traj": NTRAJ, "dt_traj": DT_TRAJ,
                     "sup_norm_ecdf_vs_runmax": sup,
                     "dkw_95_band": float(dkw),
                     "dkw_2x_band": float(2 * dkw),
                     "min_ordering_gap": float(min_gap),
                     "n_steps_with_inversion": int(n_inv),
                     "n_crossed": int(n_crossed),
                     "runmax_CDF_at_T": float(CDF[star, -1])}

    # ---- G4: no-backflow control family (same k, different centers) ----
    K_C = 1.8
    x0a, x0b = -8.0, -11.0
    adjusted = False
    for attempt in range(2):
        par_a, par_b = (K_C, x0a), (K_C, x0b)
        cIaa_R, cIbb_R, cIab_R, cIaa_F, cIbb_F, cIab_F = \
            pair_integrals(*par_a, *par_b)
        c_norm_dev = max(np.max(np.abs(cIaa_F - 1.0)),
                         np.max(np.abs(cIbb_F - 1.0)))
        c_overlap = cIab_F[0]
        cNmat, _ = family_N(cIaa_R, cIbb_R, cIab_R, c_overlap)
        c_backflow = float(np.max(np.maximum.accumulate(cNmat, axis=1)
                                  - cNmat))
        if c_backflow <= 1e-10:
            break
        if attempt == 0:
            x0b = -12.0
            adjusted = True
    _, bins_ctrl = cdf_and_bins(cNmat)
    c_design = design_matrix(c_overlap)
    resid_ctrl_perbin, _, _ = fit_povm_form(bins_ctrl, c_design)
    res["G4"] = {"k_control": K_C, "x0a": x0a, "x0b": x0b,
                 "separation_adjusted": adjusted,
                 "control_overlap_abs": float(np.abs(c_overlap)),
                 "control_norm_max_dev": float(c_norm_dev),
                 "control_max_backflow": c_backflow,
                 "control_resid_perbin_max": float(np.max(resid_ctrl_perbin)),
                 "control_resid_perbin": resid_ctrl_perbin.tolist(),
                 "ten_x_floor": 10.0 * floor}

    # ---- gate verdicts ----
    g1_pass = (norm_dev <= 1e-10 and sup <= 2 * dkw and n_inv == 0
               and backflow_per_prep[star] > 0)
    g2_pass = floor <= 1e-6
    g3_pass = len(g3_bins) > 0
    g4_pass = (c_backflow <= 1e-10
               and float(np.max(resid_ctrl_perbin)) <= 10.0 * floor)
    res["gates"] = {"G1": bool(g1_pass), "G2": bool(g2_pass),
                    "G3": bool(g3_pass), "G4": bool(g4_pass)}
    res["params"] = {"sigma": SIG, "x0": X0, "k1": K1, "k2": K2,
                     "detector": D, "T": T_HOR, "dt": DT, "nbins": NBIN,
                     "thetas": THETAS.tolist(), "phis": PHIS.tolist()}
    res["walltime_s"] = walltime.time() - t0

    with open(f"{RNG_DIR}/L7a_results.json", "w") as f:
        json.dump(res, f, indent=2)

    # save arrays for the figure
    np.savez(f"{RNG_DIR}/l7a_arrays.npz",
             TGRID=TGRID, Nstar=Nmat[star], CDFstar=CDF[star], ecdf=ecdf,
             resid_bohm_perbin=resid_bohm_perbin,
             resid_floor_perbin=resid_floor_perbin,
             resid_ctrl_perbin=resid_ctrl_perbin,
             bin_backflow_dev=bin_backflow_dev,
             backflow_bins=backflow_bins, star=star,
             theta_star=th_f[star], phi_star=ph_f[star],
             Nmat=Nmat, dkw=dkw)

    # console summary
    print(json.dumps({k: v for k, v in res.items()
                      if k not in ("bohm_resid_perbin", "G2_floor_perbin_max",
                                   "bin_backflow_dev")}, indent=2))


if __name__ == "__main__":
    main()
