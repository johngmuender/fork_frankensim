#!/usr/bin/env python3
"""
L7a diagnostics: settle whether the absence/presence of backflow in the two
families is exact physics or a grid artifact.

Current matrix: for psi = c1 psi_a + c2 psi_b, the probability current at the
detector is j(0,t) = c^dagger J(t) c with
    J_jl(t) = [psi_j^*(0,t) dx_psi_l(0,t) - (dx_psi_j(0,t))^* psi_l(0,t)] / (2i)
Hermitian. lambda_min(J(t)) >= 0 for all t  <=>  no state in span{psi_a,psi_b}
ever has negative detector current  =>  zero backflow for the WHOLE continuum
family (any theta, phi). lambda_min < 0 in a window => negative current states
exist there.
"""
import json
import numpy as np

import importlib.util
spec = importlib.util.spec_from_file_location(
    "l7a", "/home/user/fork_frankensim/analysis/gum-sandbox/tier6-foundations/L7a/l7a_povm_exclusion.py")
l7a = importlib.util.module_from_spec(spec)
spec.loader.exec_module.__self__ if False else None
# import module without running main()
import sys
sys.modules["l7a"] = l7a
spec.loader.exec_module(l7a)

psi_comp = l7a.psi_comp
logderiv_comp = l7a.logderiv_comp
D = l7a.D


def current_matrix_lambda_min(par_a, par_b, tgrid):
    """lambda_min of the 2x2 Hermitian current matrix at x=D, for each t."""
    lam = np.empty(len(tgrid))
    for i, t in enumerate(tgrid):
        ps = np.array([psi_comp(D, t, *par_a), psi_comp(D, t, *par_b)])
        dps = np.array([psi_comp(D, t, *par_a) * logderiv_comp(D, t, *par_a),
                        psi_comp(D, t, *par_b) * logderiv_comp(D, t, *par_b)])
        M = np.conj(ps)[:, None] * dps[None, :]
        J = (M - M.conj().T) / 2j
        lam[i] = np.linalg.eigvalsh(J)[0]
    return lam


def family_backflow_scan(par_a, par_b, overlap, n_th=121, n_ph=128):
    """Max backflow over a fine continuum-approximating (theta,phi) grid."""
    ths = np.linspace(0, np.pi / 2, n_th)
    phs = np.linspace(0, 2 * np.pi, n_ph, endpoint=False)
    # reuse the pair-integral machinery at the registered dt
    global_pars = (par_a[0], par_a[1], par_b[0], par_b[1])
    Iaa_R, Ibb_R, Iab_R, *_ = l7a.pair_integrals(*global_pars)
    TH, PH = np.meshgrid(ths, phs, indexing="ij")
    th = TH.ravel(); ph = PH.ravel()
    c2 = np.cos(th) ** 2; s2 = np.sin(th) ** 2; s2t = np.sin(2 * th)
    n2 = 1.0 + s2t * (np.cos(ph) * overlap.real - np.sin(ph) * overlap.imag)
    quad = (c2[:, None] * Iaa_R[None, :] + s2[:, None] * Ibb_R[None, :]
            + s2t[:, None] * (np.cos(ph)[:, None] * Iab_R.real[None, :]
                              - np.sin(ph)[:, None] * Iab_R.imag[None, :]))
    N = quad / n2[:, None]
    bf = np.max(np.maximum.accumulate(N, axis=1) - N, axis=1)
    imax = int(np.argmax(bf))
    return float(bf.max()), float(th[imax]), float(ph[imax])


def main():
    out = {}
    tfine = np.arange(0.0, 12.0 + 1e-12, 0.0005)

    # ---- main family: co-located, k1 != k2 ----
    par1, par2 = (l7a.K1, l7a.X0), (l7a.K2, l7a.X0)
    ov_main = np.exp(-0.25 * (l7a.K2 - l7a.K1) ** 2) * np.exp(
        1j * (l7a.K2 - l7a.K1) * l7a.X0)
    lam_main = current_matrix_lambda_min(par1, par2, tfine)
    out["main_lambda_min_over_t"] = float(lam_main.min())
    out["main_lambda_min_argmin_t"] = float(tfine[np.argmin(lam_main)])
    bf, th_b, ph_b = family_backflow_scan(par1, par2, ov_main)
    out["main_fine_scan_max_backflow"] = bf
    out["main_fine_scan_argmax"] = {"theta": th_b, "phi": ph_b}

    # local wavenumber mismatch at detector vs time (diagnosis)
    tt = np.array([2.0, 3.33, 5.0, 6.67, 10.0])
    out["local_dk_at_detector"] = {
        "t": tt.tolist(),
        "dk_loc": ((l7a.K2 - l7a.K1) / (1 + tt ** 2)).tolist(),
        "note": "phase-gradient mismatch of the two components at x=0: "
                "dk_loc = (k2-k1) sigma^4/(sigma^4+t^2); backflow needs "
                "amplitude-weighted mismatch to beat the mean local k ~ 1.6"}

    # ---- control family: same k, displaced centers; separation scan ----
    KC = 1.8
    seps = [1.0, 1.5, 2.0, 3.0, 4.0]
    ctrl = {}
    for s in seps:
        pa, pb = (KC, -8.0), (KC, -8.0 - s)
        ovc = np.exp(-s ** 2 / 4.0)  # same k: overlap real, exp(-sep^2/(4 sig^2))
        lam = current_matrix_lambda_min(pa, pb, tfine)
        bf, th_b, ph_b = family_backflow_scan(pa, pb, ovc + 0j)
        ctrl[str(s)] = {"lambda_min": float(lam.min()),
                        "lambda_min_t": float(tfine[np.argmin(lam)]),
                        "fine_scan_max_backflow": bf,
                        "argmax_theta": th_b, "argmax_phi": ph_b}
    out["control_separation_scan"] = ctrl

    # ---- control family (sep=3, registered): max-backflow grid prep,
    #      trajectory fan validation there ----
    pa, pb = (KC, -8.0), (KC, -11.0)
    ovc = np.exp(-9.0 / 4.0) + 0j
    Iaa_R, Ibb_R, Iab_R, Iaa_F, Ibb_F, Iab_F = l7a.pair_integrals(
        KC, -8.0, KC, -11.0)
    ovc_num = Iab_F[0]
    cN, _ = l7a.family_N(Iaa_R, Ibb_R, Iab_R, ovc_num)
    bf_pp = np.max(np.maximum.accumulate(cN, axis=1) - cN, axis=1)
    cstar = int(np.argmax(bf_pp))
    th_f, ph_f = l7a.family_grid()
    out["ctrl_grid_star"] = {"index": cstar,
                             "theta_over_pi": float(th_f[cstar] / np.pi),
                             "phi_over_pi": float(ph_f[cstar] / np.pi),
                             "backflow": float(bf_pp[cstar])}
    CDFc = np.maximum.accumulate(cN, axis=1)
    ecdf, min_gap, n_inv, n_crossed = l7a.run_fan(
        th_f[cstar], ph_f[cstar], ovc_num, pa, pb, dt=0.001)
    dkw = np.sqrt(np.log(2 / 0.05) / (2 * l7a.NTRAJ))
    out["ctrl_fan"] = {
        "dt_traj": 0.001,
        "sup_norm_ecdf_vs_runmax": float(np.max(np.abs(ecdf - CDFc[cstar]))),
        "dkw_2x_band": float(2 * dkw),
        "min_ordering_gap": float(min_gap),
        "n_steps_with_inversion": int(n_inv),
        "n_crossed": int(n_crossed)}

    # continuity check dN/dt vs j(0,t) for the control star prep
    n2c = 1 + np.sin(2 * th_f[cstar]) * (np.cos(ph_f[cstar]) * ovc_num.real
                                         - np.sin(ph_f[cstar]) * ovc_num.imag)
    c1 = np.cos(th_f[cstar]) / np.sqrt(n2c)
    c2 = np.sin(th_f[cstar]) * np.exp(1j * ph_f[cstar]) / np.sqrt(n2c)
    jt = np.empty(l7a.NT)
    for i, t in enumerate(l7a.TGRID):
        p = c1 * psi_comp(D, t, *pa) + c2 * psi_comp(D, t, *pb)
        dp = (c1 * psi_comp(D, t, *pa) * logderiv_comp(D, t, *pa)
              + c2 * psi_comp(D, t, *pb) * logderiv_comp(D, t, *pb))
        jt[i] = np.imag(np.conj(p) * dp)
    dNdt = np.gradient(cN[cstar], l7a.TGRID)
    out["ctrl_continuity_max_dev"] = float(
        np.max(np.abs(dNdt[2:-2] - jt[2:-2])))
    out["ctrl_min_current"] = float(jt.min())

    with open("/home/user/fork_frankensim/analysis/gum-sandbox/"
              "tier6-foundations/L7a/l7a_diagnostics.json", "w") as f:
        json.dump(out, f, indent=2)
    np.savez("/home/user/fork_frankensim/analysis/gum-sandbox/"
             "tier6-foundations/L7a/l7a_diag_arrays.npz",
             tfine=tfine, lam_main=lam_main,
             ctrl_jt=jt, ctrl_Nstar=cN[cstar], ctrl_CDFstar=CDFc[cstar],
             ctrl_ecdf=ecdf, cstar=cstar)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
