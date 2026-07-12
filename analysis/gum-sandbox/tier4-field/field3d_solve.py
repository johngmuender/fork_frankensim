#!/usr/bin/env python3
"""
Tier-4A of the GUM replication campaign: fully 3-D Cartesian, symmetry-free
adjudication of the Step-3 fixed-L Routhian claims (F-R5) at eps = 0.05.

Reproduce with `python3 field3d_solve.py` (deterministic: the only RNG uses
fixed seeds; writes field3d_results.json, field3d_decay.png, field3d_run.log).
`python3 field3d_solve.py --smoke` runs a tiny-N smoke test (separate outputs).

Conventions (Tier-2 unit map, axi_RESULTS.md, a6 = a0 = m = 1):
  field       q = (q0,q1,q2,q3), |q| = 1, vacuum q = (1,0,0,0)
  E2 = (1/4pi) INT sum_i |a_i|^2      -> FD: sum_i |D_i q|^2
  E4 = (1/4pi) INT sum_{i<j} |a_i x a_j|^2
       -> FD: sum_{i<j} (|D_i|^2 |D_j|^2 - (D_i.D_j)^2)   (D_i q _|_ q)
  E6 = pi^3 INT b^2,  b = sgn/(2 pi^2) det[q, D_x q, D_y q, D_z q]
       -> density det^2/(4pi); degree = INT b (sgn fixed so hedgehog deg = +1)
  E0 = (1/4pi) INT (1 - q0)
  I  = INT 2 (q1^2 + q2^2)            (isorotation inertia, 3rd iso-axis)
  E_static = t (E2 + E4) + E6 + E0,  t = 0.0082764349 (eps = 0.05 dial)
  R(q; L)  = E_static + L^2 / (2 I)   (fixed-L Routhian)
  halo threshold (Step 3, F-R5): kappa_ours = L/I > 1/sqrt(8 pi) = 0.19947
  weighted Derrick virial: t E2 - t E4 - 3 E6 + 3 E0 = 0
  clock: L^2 = (2/3) I E_static  ->  E_rot/E_tot = 1/4

Grid/discretization: cell-centred N^3 Cartesian box, half-width LBOX = 4.5
(~2.5 R*, R* = 2^(5/6); Yukawa tail mu = 7.77 -> boundary tail < 1e-10 of max).
4th-order central differences with a fixed-vacuum ghost rind (2 layers); all
integrals are plain h^3 cell sums.  1-D hedgehog profile re-solved here with
the Step-1 method (graded-grid midpoint energy, gradient flow + damped Newton
with exact tridiagonal Hessian), t frozen -- no eps re-dial.

Minimizer: arrested Newton flow (velocity field on q, acceleration =
-(tangent-projected) functional gradient, velocity zeroed and the step
reverted whenever the objective increases; adaptive dt), unit norm enforced
by renormalization after every update, velocity re-projected to the tangent
space of the updated field.

Stages:
  G  gates: hedgehog sectors at N = 48/64/96 vs step-1 radial references,
     degree, virial, FD-vs-analytic gradient check, boundary-tail check
  A  on-grid axisymmetric family-A reference: min over dilation d of
     R(hedgehog(x/d); L) -- the same-discretization stand-in for Step-2 solA
  S  static minimization from hedgehog + fixed-seed perturbation
  M  fixed-L Routhian descent at L = 8.4979 (Step-2 solA charge)
  C  control descent at L = 4.0 (below threshold at the hedgehog)
  B  box-size sanity: stage-M rerun at 1.3x box (same N, coarser h)
  K  clock: L_clock from L^2 = (2/3) I E_static on the stage-C solution
     (bisection at frozen fields = single fixed-point pass), kappa(L_clock)
     vs the halo threshold.

Python 3 + numpy only; matplotlib (Agg) for the figure.
"""

import json
import sys
import time

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# constants / conventions
# ----------------------------------------------------------------------
FOURPI = 4.0 * np.pi
T_FROZEN = 0.008276434949296802      # eps = 0.05 dial (Tier-2 step 1)
RSTAR = 2.0 ** (5.0 / 6.0)           # compacton radius, a6 = a0 = m = 1
KC = 1.0 / np.sqrt(8.0 * np.pi)      # halo threshold kappa_ours (F-R5)
L_MAIN = 8.4979                      # Step-2 solA isospin (task value)
L_CTRL = 4.0                         # control: kappa(hedgehog) < KC
HALO_RADIUS = 1.5 * RSTAR            # halo diagnostic radius

# Step-1 radial references (2-D/1-D solves, eps = 0.05 hedgehog)
REF = dict(E2=12.117459, E4=6.191126, E6=1.523936, E0=1.507588, I=22.26826)
# Step-2 solA (spectral axisymmetric stationary point, from axi_results.json)
SOLA = dict(L=8.497910448038807, R=4.513717738366562, Estat=3.3852883016361837,
            I=31.99778365944103, kappa=0.2655780956107401, d=1.1284383902516937)

OUTDIR = "/home/user/fork_frankensim/analysis/gum-sandbox/tier4-field"

SMOKE = "--smoke" in sys.argv
TAG = "_smoke" if SMOKE else ""

# grid / iteration budget (SMOKE shrinks everything; correctness unchanged)
N_MAIN = 32 if SMOKE else 96
N_GATE_LO = 24 if SMOKE else 64
N_GATE_XLO = 16 if SMOKE else 48
LBOX = 4.5
LBOX_BIG = 4.5 * 1.3
MAXIT_STATIC = 60 if SMOKE else 1200
MAXIT_MAIN = 80 if SMOKE else 2400
MAXIT_CTRL = 60 if SMOKE else 1400
MAXIT_BIG = 60 if SMOKE else 1600
INSTR_EVERY = 10 if SMOKE else 25
PRINT_EVERY = 20 if SMOKE else 200
GTOL_RATIO = 1.0e-6
SEED_PERT = 424242

LOGFH = None


def log(msg=""):
    print(msg)
    if LOGFH is not None:
        LOGFH.write(msg + "\n")
        LOGFH.flush()


# ----------------------------------------------------------------------
# 1-D radial hedgehog profile at frozen t (Step-1 method, condensed)
# ----------------------------------------------------------------------
def radial_grid(N, rmax):
    s = np.linspace(0.0, 1.0, 200001)
    rr = s * rmax
    w = (0.55 + 2.0 * np.exp(-(rr / 0.5) ** 2)
         + 3.0 * np.exp(-((rr - RSTAR) / 0.45) ** 2)
         + 0.8 * np.exp(-((rr - (RSTAR + 1.2)) / 0.9) ** 2))
    cw = np.concatenate([[0.0], np.cumsum(0.5 * (w[1:] + w[:-1]) * np.diff(rr))])
    cw /= cw[-1]
    r = np.interp(np.linspace(0.0, 1.0, N + 1), cw, rr)
    r[0] = 0.0
    r[-1] = rmax
    return r


def radial_assemble(r, f, t):
    h = np.diff(r)
    rm = 0.5 * (r[1:] + r[:-1])
    rm2 = rm * rm
    fm = 0.5 * (f[1:] + f[:-1])
    d = np.diff(f) / h
    d2 = d * d
    s = np.sin(fm)
    c = np.cos(fm)
    s2, s3, s4 = s * s, s * s * s, (s * s) ** 2
    sin2f = 2.0 * s * c
    cos2f = c * c - s2
    E2 = float(np.sum(h * (d2 * rm2 + 2.0 * s2)))
    E4 = float(np.sum(h * (s2 * (2.0 * d2 + s2 / rm2))))
    E6 = float(np.sum(h * (s4 * d2 / rm2)))
    E0 = float(np.sum(h * ((1.0 - np.cos(fm)) * rm2)))
    E = t * (E2 + E4) + E6 + E0
    L_d = 2.0 * t * d * rm2 + 4.0 * t * s2 * d + 2.0 * s4 * d / rm2
    L_f = 2.0 * t * sin2f * (1.0 + d2) + 4.0 * (t + d2) * s3 * c / rm2 + s * rm2
    gL = 0.5 * h * L_f - L_d
    gR = 0.5 * h * L_f + L_d
    g = gR[:-1] + gL[1:]
    L_dd = 2.0 * t * rm2 + 4.0 * t * s2 + 2.0 * s4 / rm2
    L_fd = 4.0 * t * d * sin2f + 8.0 * s3 * c * d / rm2
    L_ff = (4.0 * t * cos2f * (1.0 + d2)
            + 4.0 * (t + d2) * (3.0 * s2 * c * c - s4) / rm2 + c * rm2)
    HLL = 0.25 * h * L_ff - L_fd + L_dd / h
    HRR = 0.25 * h * L_ff + L_fd + L_dd / h
    HLR = 0.25 * h * L_ff - L_dd / h
    diag = HRR[:-1] + HLL[1:]
    off = HLR[1:-1]
    return E, (E2, E4, E6, E0), g, diag, off


def radial_thomas(diag, off, b):
    n = diag.size
    dd = diag.copy()
    bb = b.copy()
    for i in range(1, n):
        m = off[i - 1] / dd[i - 1]
        dd[i] -= m * off[i - 1]
        bb[i] -= m * bb[i - 1]
    x = np.zeros(n)
    x[-1] = bb[-1] / dd[-1]
    for i in range(n - 2, -1, -1):
        x[i] = (bb[i] - off[i] * x[i + 1]) / dd[i]
    return x


def radial_solve(t, N=4000, rmax=6.0):
    r = radial_grid(N, rmax)
    u = r / RSTAR
    f = 2.0 * np.arccos(np.clip(u / (1.0 + u ** 8) ** 0.125, 0.0, 1.0))
    f[0] = np.pi
    f[-1] = 0.0
    E, _, g, diag, _ = radial_assemble(r, f, t)
    dt = 0.5
    for _ in range(6000):                       # preconditioned flow
        if np.max(np.abs(g)) < 1e-2:
            break
        prec = np.maximum(diag, 1e-3 * np.max(diag))
        fn = f.copy()
        fn[1:-1] = f[1:-1] - dt * g / prec
        np.clip(fn, 0.0, np.pi, out=fn)
        En, _, gn, diagn, _ = radial_assemble(r, fn, t)
        if En < E:
            f, E, g, diag = fn, En, gn, diagn
            dt = min(dt * 1.1, 5.0)
        else:
            dt *= 0.5
            if dt < 1e-9:
                break
    lam = 1e-3
    for _ in range(300):                        # damped Newton to the floor
        E, sect, g, diag, off = radial_assemble(r, f, t)
        gmax = float(np.max(np.abs(g)))
        if gmax < 1e-12:
            break
        scale = np.abs(diag) + 1e-30
        ok = False
        while lam < 1e14:
            step = radial_thomas(diag + lam * scale, off, -g)
            fn = f.copy()
            fn[1:-1] += step
            En, _, gn, _, _ = radial_assemble(r, fn, t)
            if np.isfinite(En) and np.max(np.abs(gn)) < gmax:
                f, lam, ok = fn, max(lam * 0.25, 1e-14), True
                break
            lam *= 10.0
        if not ok:
            break
    E, sect, g, _, _ = radial_assemble(r, f, t)
    return r, f, sect, float(np.max(np.abs(g)))


# ----------------------------------------------------------------------
# 3-D Cartesian engine: 4th-order central differences, vacuum ghosts
# ----------------------------------------------------------------------
# pair expansion of det[q, Dx, Dy, Dz]:
#   det = sum_p s_p (q_a Dx_b - q_b Dx_a)(Dy_c Dz_d - Dy_d Dz_c)
PAIRS = [(0, 1, 2, 3, 1.0), (0, 2, 1, 3, -1.0), (0, 3, 1, 2, 1.0),
         (1, 2, 0, 3, 1.0), (1, 3, 0, 2, -1.0), (2, 3, 0, 1, 1.0)]


class Engine:
    def __init__(self, N, Lbox, t=T_FROZEN):
        self.N = N
        self.Lbox = Lbox
        self.t = t
        self.h = 2.0 * Lbox / N
        self.h3 = self.h ** 3
        self.x1 = (np.arange(N) + 0.5) * self.h - Lbox
        x = self.x1
        self.r = np.sqrt(x[:, None, None] ** 2 + x[None, :, None] ** 2
                         + x[None, None, :] ** 2)
        self.sgn6 = -1.0    # measured on the hedgehog: raw degree < 0
        sh = (N, N, N)
        self.qp = np.zeros((4, N + 4, N + 4, N + 4))
        self.qp[0] = 1.0                       # vacuum ghosts
        self.D = np.empty((3, 4) + sh)         # D_i q
        self.P = np.empty((3, 4) + sh)         # fluxes d(dens)/d(D_i q)
        self.g = np.empty((4,) + sh)
        self.M6 = np.empty((6,) + sh)
        self.Pc6 = np.empty((6,) + sh)
        self.det = np.empty(sh)
        self.w6 = np.empty(sh)
        self.nn = np.empty((3,) + sh)          # |D_i|^2
        self.dd = np.empty((3,) + sh)          # D_x.D_y, D_x.D_z, D_y.D_z
        self.T = [np.empty(sh) for _ in range(4)]

    # -- fields --------------------------------------------------------
    def hedgehog(self, rnodes, fnodes, d=1.0):
        """Family-A dilation q_d(x) = q_base(x/d) of the radial profile."""
        f = np.interp(self.r / d, rnodes, fnodes)
        sf_r = np.sin(f) / self.r
        q = np.empty((4, self.N, self.N, self.N))
        q[0] = np.cos(f)
        q[1] = sf_r * self.x1[:, None, None]
        q[2] = sf_r * self.x1[None, :, None]
        q[3] = sf_r * self.x1[None, None, :]
        return q

    def perturbation(self, seed=SEED_PERT, K=6, amp=0.02):
        """Fixed-seed smooth non-axisymmetric bump field (4 components)."""
        rng = np.random.default_rng(seed)
        x = self.x1
        dq = np.zeros((4, self.N, self.N, self.N))
        for _ in range(K):
            c = rng.uniform(-1.25 * RSTAR, 1.25 * RSTAR, size=3)
            w = rng.uniform(0.5, 0.9)
            A = rng.uniform(-amp, amp, size=4)
            gau = np.exp(-((x[:, None, None] - c[0]) ** 2
                           + (x[None, :, None] - c[1]) ** 2
                           + (x[None, None, :] - c[2]) ** 2) / (2.0 * w * w))
            dq += A[:, None, None, None] * gau
        return dq

    @staticmethod
    def normalize(q):
        nrm = np.sqrt(q[0] ** 2 + q[1] ** 2 + q[2] ** 2 + q[3] ** 2)
        q /= nrm
        return q

    # -- 4th-order derivatives ------------------------------------------
    def _view(self, a, axis, off):
        N = self.N
        sl = [slice(2, 2 + N)] * 3
        sl[axis] = slice(2 + off, 2 + off + N)
        return self.qp[a][tuple(sl)]

    def energy(self, q, L=None, need_grad=True):
        """Sector energies + (optionally) tangent-projected gradient of
        E_static (L None) or of R = E_static + L^2/(2I).  The returned
        gradient array is an internal buffer: copy it before the next call."""
        N, t, h3 = self.N, self.t, self.h3
        c1 = 8.0 / (12.0 * self.h)
        c2 = 1.0 / (12.0 * self.h)
        self.qp[:, 2:-2, 2:-2, 2:-2] = q
        D, T0 = self.D, self.T[0]
        for ax in range(3):
            for a in range(4):
                Dv = D[ax, a]
                np.subtract(self._view(a, ax, 1), self._view(a, ax, -1), out=Dv)
                Dv *= c1
                np.subtract(self._view(a, ax, 2), self._view(a, ax, -2), out=T0)
                T0 *= c2
                Dv -= T0
        nn, dd = self.nn, self.dd
        for ax in range(3):
            v, o = D[ax], nn[ax]
            np.multiply(v[0], v[0], out=o)
            for a in range(1, 4):
                np.multiply(v[a], v[a], out=T0)
                o += T0
        for k, (i, j) in enumerate(((0, 1), (0, 2), (1, 2))):
            o = dd[k]
            np.multiply(D[i, 0], D[j, 0], out=o)
            for a in range(1, 4):
                np.multiply(D[i, a], D[j, a], out=T0)
                o += T0
        E2 = (float(np.sum(nn[0])) + float(np.sum(nn[1]))
              + float(np.sum(nn[2]))) * h3 / FOURPI
        S4 = 0.0
        for k, (i, j) in enumerate(((0, 1), (0, 2), (1, 2))):
            np.multiply(nn[i], nn[j], out=T0)
            S4 += float(np.sum(T0))
            np.multiply(dd[k], dd[k], out=T0)
            S4 -= float(np.sum(T0))
        E4 = S4 * h3 / FOURPI
        det, M6, Pc6 = self.det, self.M6, self.Pc6
        Dx, Dy, Dz = D[0], D[1], D[2]
        for p, (a, b, c, d, sg) in enumerate(PAIRS):
            Mp, Pc = M6[p], Pc6[p]
            np.multiply(q[a], Dx[b], out=Mp)
            np.multiply(q[b], Dx[a], out=T0)
            Mp -= T0
            np.multiply(Dy[c], Dz[d], out=Pc)
            np.multiply(Dy[d], Dz[c], out=T0)
            Pc -= T0
            np.multiply(Mp, Pc, out=T0)
            if p == 0:
                np.copyto(det, T0)
            elif sg > 0:
                det += T0
            else:
                det -= T0
        np.multiply(det, det, out=T0)
        E6 = float(np.sum(T0)) * h3 / FOURPI
        E0 = (q[0].size - float(np.sum(q[0]))) * h3 / FOURPI
        deg = self.sgn6 * float(np.sum(det)) * h3 / (2.0 * np.pi ** 2)
        np.multiply(q[1], q[1], out=T0)
        Ival = float(np.sum(T0))
        np.multiply(q[2], q[2], out=T0)
        Ival = (Ival + float(np.sum(T0))) * 2.0 * h3
        Estat = t * (E2 + E4) + E6 + E0
        out = dict(E2=E2, E4=E4, E6=E6, E0=E0, I=Ival, deg=deg, Estat=Estat)
        out["R"] = Estat if L is None else Estat + L * L / (2.0 * Ival)
        if not need_grad:
            return out, None
        # ---- gradient --------------------------------------------------
        g, P, T1, T2, T3 = self.g, self.P, self.T[1], self.T[2], self.T[3]
        g[0].fill(-1.0 / FOURPI)               # e0 density part
        g[1].fill(0.0)
        g[2].fill(0.0)
        g[3].fill(0.0)
        if L is not None:                      # -(L^2/2I^2) dI/dq (density)
            fac = -(L * L) / (2.0 * Ival * Ival) * 4.0
            np.multiply(q[1], fac, out=T0)
            g[1] += T0
            np.multiply(q[2], fac, out=T0)
            g[2] += T0
        # (2+4)-sector fluxes, weight t/4pi
        w24 = 2.0 * t / FOURPI
        DOT = {(0, 1): 0, (0, 2): 1, (1, 2): 2}
        for ax, (o1, o2) in ((0, (1, 2)), (1, (0, 2)), (2, (0, 1))):
            np.add(nn[o1], nn[o2], out=T1)
            T1 += 1.0                          # "+1" carries the E2 flux
            d1 = dd[DOT[tuple(sorted((ax, o1)))]]
            d2 = dd[DOT[tuple(sorted((ax, o2)))]]
            for a in range(4):
                Pv = P[ax, a]
                np.multiply(T1, D[ax, a], out=Pv)
                np.multiply(d1, D[o1, a], out=T0)
                Pv -= T0
                np.multiply(d2, D[o2, a], out=T0)
                Pv -= T0
                Pv *= w24
        # sextic: d(det^2/4pi)/d det = det/2pi
        np.multiply(det, 1.0 / (2.0 * np.pi), out=self.w6)
        w6 = self.w6
        for p, (a, b, c, d, sg) in enumerate(PAIRS):
            np.multiply(w6, Pc6[p], out=T2)
            if sg < 0:
                T2 *= -1.0
            np.multiply(T2, Dx[b], out=T0)
            g[a] += T0
            np.multiply(T2, Dx[a], out=T0)
            g[b] -= T0
            np.multiply(T2, q[a], out=T0)
            P[0, b] += T0
            np.multiply(T2, q[b], out=T0)
            P[0, a] -= T0
            np.multiply(w6, M6[p], out=T3)
            if sg < 0:
                T3 *= -1.0
            np.multiply(T3, Dz[d], out=T0)
            P[1, c] += T0
            np.multiply(T3, Dz[c], out=T0)
            P[1, d] -= T0
            np.multiply(T3, Dy[c], out=T0)
            P[2, d] += T0
            np.multiply(T3, Dy[d], out=T0)
            P[2, c] -= T0
        # divergence: adjoint of the 4th-order stencil (ghost flux = 0)
        for a in range(4):
            gv = g[a]
            Pv = P[0, a]
            gv[1:, :, :] += c1 * Pv[:-1, :, :]
            gv[:-1, :, :] -= c1 * Pv[1:, :, :]
            gv[:-2, :, :] += c2 * Pv[2:, :, :]
            gv[2:, :, :] -= c2 * Pv[:-2, :, :]
            Pv = P[1, a]
            gv[:, 1:, :] += c1 * Pv[:, :-1, :]
            gv[:, :-1, :] -= c1 * Pv[:, 1:, :]
            gv[:, :-2, :] += c2 * Pv[:, 2:, :]
            gv[:, 2:, :] -= c2 * Pv[:, :-2, :]
            Pv = P[2, a]
            gv[:, :, 1:] += c1 * Pv[:, :, :-1]
            gv[:, :, :-1] -= c1 * Pv[:, :, 1:]
            gv[:, :, :-2] += c2 * Pv[:, :, 2:]
            gv[:, :, 2:] -= c2 * Pv[:, :, :-2]
        g *= h3
        # tangent projection
        np.multiply(g[0], q[0], out=T0)
        for a in range(1, 4):
            np.multiply(g[a], q[a], out=T1)
            T0 += T1
        for a in range(4):
            np.multiply(T0, q[a], out=T1)
            g[a] -= T1
        return out, g

    # -- diagnostics ----------------------------------------------------
    def halo_fraction(self, q):
        """I-fraction outside r > 1.5 R* around the (1-q0) centroid."""
        w = 1.0 - q[0]
        W = float(np.sum(w))
        x = self.x1
        cx = float(np.sum(w.sum(axis=(1, 2)) * x)) / W
        cy = float(np.sum(w.sum(axis=(0, 2)) * x)) / W
        cz = float(np.sum(w.sum(axis=(0, 1)) * x)) / W
        r2c = ((x[:, None, None] - cx) ** 2 + (x[None, :, None] - cy) ** 2
               + (x[None, None, :] - cz) ** 2)
        dI = q[1] ** 2 + q[2] ** 2
        tot = float(np.sum(dI))
        halo = float(np.sum(dI[r2c > HALO_RADIUS ** 2]))
        return halo / tot, (cx, cy, cz)


def virial(out, t=T_FROZEN):
    return t * out["E2"] - t * out["E4"] - 3.0 * out["E6"] + 3.0 * out["E0"]


# ----------------------------------------------------------------------
# arrested Newton flow
# ----------------------------------------------------------------------
def anf(eng, q, L, maxit, label, dt0=0.02, dt_max=0.60,
        gtol_ratio=GTOL_RATIO, instr_every=INSTR_EVERY):
    """Arrested Newton flow on E_static (L None) or R(q; L).  Returns the
    final field, final sector dict, the instrumented time series, and a
    status string.  Objective is monotone by construction (reject+arrest)."""
    h3 = eng.h3
    out, gbuf = eng.energy(q, L=L, need_grad=True)
    g = gbuf.copy()
    obj = out["R"]
    gn0 = float(np.sqrt(np.sum(g * g))) / h3
    v = np.zeros_like(q)
    qb = q.copy()
    dt = dt0
    arrests = 0
    series = []
    status = "iter_cap"

    def record(it, gn, dtv):
        hf, cen = eng.halo_fraction(q)
        series.append(dict(
            it=it, R=obj, Estat=out["Estat"], E2=out["E2"], E4=out["E4"],
            E6=out["E6"], E0=out["E0"], I=out["I"],
            kappa=(0.0 if L is None else L / out["I"]),
            deg=out["deg"], halo=hf, gnorm=gn, dt=dtv, arrests=arrests,
            cen=[round(c, 5) for c in cen]))
        return series[-1]

    rec = record(0, gn0, dt)
    log(f"  [{label}] it=0  obj={obj:.7f}  I={out['I']:.4f} "
        f"kappa={rec['kappa']:.5f} deg={out['deg']:.5f} halo={rec['halo']:.4f} "
        f"gnorm={gn0:.3e}")
    it = 0
    while it < maxit:
        it += 1
        v -= (dt / h3) * g
        q += dt * v
        eng.normalize(q)
        out_new, gbuf = eng.energy(q, L=L, need_grad=True)
        if out_new["R"] > obj:                 # arrest: revert, kill velocity
            q[...] = qb
            v.fill(0.0)
            dt *= 0.6
            arrests += 1
            if dt < 1e-7:
                status = "stall_dt"
                break
            # g of the accepted point is still valid (stored copy)
        else:
            out = out_new
            obj = out["R"]
            np.copyto(g, gbuf)
            np.copyto(qb, q)
            # re-project velocity onto the new tangent space
            dot = v[0] * q[0] + v[1] * q[1] + v[2] * q[2] + v[3] * q[3]
            for a in range(4):
                v[a] -= dot * q[a]
            dt = min(dt * 1.005, dt_max)
        gn = float(np.sqrt(np.sum(g * g))) / h3
        if it % instr_every == 0 or it == maxit:
            rec = record(it, gn, dt)
            if it % PRINT_EVERY == 0 or it == maxit:
                log(f"  [{label}] it={it}  obj={obj:.7f}  I={out['I']:.4f} "
                    f"kappa={rec['kappa']:.5f} deg={out['deg']:.5f} "
                    f"halo={rec['halo']:.4f} gnorm={gn:.3e} dt={dt:.2e} "
                    f"arr={arrests}")
        if gn < gtol_ratio * gn0:
            status = "gtol"
            record(it, gn, dt)
            break
    q[...] = qb
    gn = float(np.sqrt(np.sum(g * g))) / h3
    log(f"  [{label}] DONE ({status}) it={it} obj={obj:.7f} "
        f"gnorm/gnorm0={gn / gn0:.3e} arrests={arrests}")
    return q, out, series, status, dict(gn0=gn0, gn=gn, iters=it,
                                        arrests=arrests, status=status)


# ----------------------------------------------------------------------
# on-grid axisymmetric family-A reference: min_d R(hedgehog(x/d); L)
# ----------------------------------------------------------------------
def axi_reference(eng, rn, fn, L, lo=-0.10, hi=0.45, iters=48):
    def Rofd(ld):
        q = eng.hedgehog(rn, fn, d=np.exp(ld))
        out, _ = eng.energy(q, L=L, need_grad=False)
        return out["R"], out

    invphi = (np.sqrt(5.0) - 1.0) / 2.0
    a, b = lo, hi
    c = b - invphi * (b - a)
    d = a + invphi * (b - a)
    fc, _ = Rofd(c)
    fd, _ = Rofd(d)
    for _ in range(iters):
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - invphi * (b - a)
            fc, _ = Rofd(c)
        else:
            a, c, fc = c, d, fd
            d = a + invphi * (b - a)
            fd, _ = Rofd(d)
    ld = 0.5 * (a + b)
    R, out = Rofd(ld)
    return float(np.exp(ld)), R, out


# ----------------------------------------------------------------------
# gates
# ----------------------------------------------------------------------
def run_gates(rn, fn):
    log("\n===== Section G: validation gates =====")
    gates = {}
    hh = {}
    for N in (N_GATE_XLO, N_GATE_LO, N_MAIN):
        eng = Engine(N, LBOX)
        q = eng.hedgehog(rn, fn)
        t0 = time.time()
        out, _ = eng.energy(q, need_grad=False)
        dt = time.time() - t0
        hh[N] = out
        rel = {k: out[k] / REF[k] - 1.0 for k in REF}
        vir = virial(out)
        log(f"  N={N} (h={eng.h:.5f}, eval {dt:.2f}s): "
            + " ".join(f"{k}={out[k]:.6f}({rel[k]:+.2e})" for k in REF))
        log(f"        degree={out['deg']:+.6f}  virial={vir:+.3e} "
            f"(rel {vir / out['Estat']:+.3e})")
        gates[f"N{N}"] = dict(
            {k: out[k] for k in ("E2", "E4", "E6", "E0", "I")},
            deg=out["deg"], rel={k: rel[k] for k in REF},
            virial=vir, virial_rel=vir / out["Estat"])
    conv = {k: abs(gates[f"N{N_GATE_LO}"]["rel"][k])
            / max(abs(gates[f"N{N_MAIN}"]["rel"][k]), 1e-30) for k in REF}
    log("  convergence ratio |err(N=%d)|/|err(N=%d)|: " % (N_GATE_LO, N_MAIN)
        + " ".join(f"{k}={conv[k]:.2f}" for k in REF)
        + f"   [O(h^2) would be {(N_MAIN / N_GATE_LO) ** 2:.2f}]")
    worst = max(abs(v) for v in gates[f"N{N_MAIN}"]["rel"].values())
    dg = abs(abs(gates[f"N{N_MAIN}"]["deg"]) - 1.0)
    log(f"  gate: worst sector rel err at N={N_MAIN}: {worst:.2e} "
        f"(tol 1e-2) -> {'PASS' if worst < 1e-2 else 'FAIL'}")
    log(f"  gate: |degree - 1| = {dg:.2e} (tol 5e-3) -> "
        f"{'PASS' if dg < 5e-3 else 'FAIL'}")
    gates["convergence_ratio"] = conv
    gates["worst_rel_Nmain"] = worst
    gates["degree_err_Nmain"] = dg
    gates["pass_sectors"] = bool(worst < 1e-2)
    gates["pass_degree"] = bool(dg < 5e-3)

    # boundary-tail check on the initial hedgehog
    engM = Engine(N_MAIN, LBOX)
    qM = engM.hedgehog(rn, fn)
    fb = np.sqrt(qM[1] ** 2 + qM[2] ** 2 + qM[3] ** 2)
    tail = max(float(fb[0].max()), float(fb[-1].max()),
               float(fb[:, 0].max()), float(fb[:, -1].max()),
               float(fb[:, :, 0].max()), float(fb[:, :, -1].max())) / np.pi
    log(f"  gate: boundary tail |sin f|_max / pi = {tail:.3e} (tol 1e-6) -> "
        f"{'PASS' if tail < 1e-6 else 'FAIL'}")
    gates["boundary_tail"] = tail
    gates["pass_tail"] = bool(tail < 1e-6)

    # FD-vs-analytic gradient check (small grid, perturbed hedgehog)
    Ng = 20
    eng = Engine(Ng, LBOX)
    q = eng.hedgehog(rn, fn)
    q = eng.normalize(q + eng.perturbation(seed=11, amp=0.05))
    rng = np.random.default_rng(12)
    worst_gc = 0.0
    for L in (None, L_MAIN):
        out, g = eng.energy(q, L=L, need_grad=True)
        g = g.copy()
        for _ in range(3):
            u = rng.standard_normal(q.shape)
            u -= (u[0] * q[0] + u[1] * q[1] + u[2] * q[2] + u[3] * q[3]) * q
            u /= np.sqrt(np.sum(u * u))
            eps = 1e-5
            qp = eng.normalize(q + eps * u)
            qm = eng.normalize(q - eps * u)
            op, _ = eng.energy(qp, L=L, need_grad=False)
            om, _ = eng.energy(qm, L=L, need_grad=False)
            fd = (op["R"] - om["R"]) / (2.0 * eps)
            an = float(np.sum(g * u))
            worst_gc = max(worst_gc, abs(fd - an) / max(abs(fd), 1e-30))
    log(f"  gate: FD-vs-analytic gradient (E_static and Routhian, 6 dirs): "
        f"worst rel {worst_gc:.2e} (tol 1e-5) -> "
        f"{'PASS' if worst_gc < 1e-5 else 'FAIL'}")
    gates["gradcheck_worst"] = worst_gc
    gates["pass_gradcheck"] = bool(worst_gc < 1e-5)
    return gates, hh


# ----------------------------------------------------------------------
# figure
# ----------------------------------------------------------------------
def make_figure(serM, serB, serC, refs, path):
    surface, ink, ink2 = "#fcfcfb", "#0b0b0b", "#52514e"
    gridc, basec = "#e1e0d9", "#c3c2b7"
    col_main, col_ctrl, col_big = "#2a78d6", "#1baf7a", "#e34948"

    def cols(ser, k):
        return [s["it"] for s in ser], [s[k] for s in ser]

    fig, axes = plt.subplots(3, 1, figsize=(8.4, 9.6), dpi=160, sharex=True)
    fig.patch.set_facecolor(surface)
    panels = [("R", "fixed-$L$ Routhian  $R = E_{stat} + L^2/2I$"),
              ("kappa", r"$\kappa_{ours} = L/I$"),
              ("halo", r"halo fraction of $I$  ($r > 1.5\,R_*$)")]
    for ax, (key, ttl) in zip(axes, panels):
        ax.set_facecolor(surface)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(basec)
        ax.grid(color=gridc, lw=0.6)
        ax.tick_params(colors=ink2, labelsize=9)
        ax.plot(*cols(serM, key), color=col_main, lw=1.8,
                label=f"$L={L_MAIN}$, box $\\pm{LBOX}$")
        ax.plot(*cols(serB, key), color=col_big, lw=1.8, ls="--",
                label=f"$L={L_MAIN}$, box $\\pm{LBOX_BIG:.2f}$ (1.3x)")
        ax.plot(*cols(serC, key), color=col_ctrl, lw=1.8,
                label=f"$L={L_CTRL}$ control")
        ax.set_title(ttl, color=ink, fontsize=11, loc="left")
    axes[0].axhline(refs["R_axi_grid_main"], color=col_main, lw=1.0, ls=":",
                    alpha=0.8)
    axes[0].axhline(refs["R_axi_grid_big"], color=col_big, lw=1.0, ls=":",
                    alpha=0.8)
    axes[0].axhline(SOLA["R"], color=ink2, lw=1.0, ls=":")
    axes[0].text(0.99, SOLA["R"], "  solA $R$ (spectral) = %.4f" % SOLA["R"],
                 color=ink2, fontsize=8, va="bottom", ha="right",
                 transform=axes[0].get_yaxis_transform())
    axes[0].text(0.99, refs["R_axi_grid_main"],
                 "axisymmetric family-A ref, this grid = %.4f  "
                 % refs["R_axi_grid_main"], color=col_main, fontsize=8,
                 va="top", ha="right",
                 transform=axes[0].get_yaxis_transform())
    axes[1].axhline(KC, color=ink2, lw=1.2)
    axes[1].text(0.01, KC, r"  halo threshold $1/\sqrt{8\pi}=0.19947$",
                 color=ink2, fontsize=8, va="bottom",
                 transform=axes[1].get_yaxis_transform())
    axes[1].axhline(SOLA["kappa"], color=ink2, lw=1.0, ls=":")
    axes[1].text(0.01, SOLA["kappa"], "  solA $\\kappa$ = 0.26558",
                 color=ink2, fontsize=8, va="bottom",
                 transform=axes[1].get_yaxis_transform())
    axes[2].set_xlabel("arrested-Newton-flow iteration", color=ink2,
                       fontsize=10)
    axes[0].legend(frameon=False, fontsize=9, loc="center right",
                   labelcolor=ink2)
    fig.suptitle("Tier-4A: symmetry-free fixed-$L$ Routhian descent, "
                 f"3-D Cartesian $N={N_MAIN}$, $\\epsilon=0.05$",
                 color=ink, fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(path, facecolor=surface)
    plt.close(fig)


# ----------------------------------------------------------------------
# main
# ----------------------------------------------------------------------
def main():
    global LOGFH
    LOGFH = open(f"{OUTDIR}/field3d_run{TAG}.log", "w")
    t_start = time.time()
    log(f"Tier-4A field3d_solve.py  (smoke={SMOKE})  N={N_MAIN} "
        f"LBOX={LBOX} h={2 * LBOX / N_MAIN:.5f}  t={T_FROZEN:.10f}")
    log(f"L_main={L_MAIN} (solA)  L_ctrl={L_CTRL}  threshold KC={KC:.6f}")

    # ---- 1-D profile ---------------------------------------------------
    rn, fn, sect, gmax = radial_solve(T_FROZEN)
    ratio = T_FROZEN * (sect[0] + sect[1]) / (sect[2] + sect[3])
    log(f"\nradial profile (t frozen): E2={sect[0]:.6f} E4={sect[1]:.6f} "
        f"E6={sect[2]:.6f} E0={sect[3]:.6f}  eps_ratio={ratio:.6f} "
        f"gmax={gmax:.1e}")

    # ---- gates ----------------------------------------------------------
    gates, hh = run_gates(rn, fn)

    # ---- engines --------------------------------------------------------
    eng = Engine(N_MAIN, LBOX)
    engB = Engine(N_MAIN, LBOX_BIG)

    # ---- on-grid axisymmetric family-A references -----------------------
    log("\n===== Section A: on-grid axisymmetric family-A references =====")
    dA, RA, outA = axi_reference(eng, rn, fn, L_MAIN)
    log(f"  main box, L={L_MAIN}: d*={dA:.5f}  R_A={RA:.6f} "
        f"(spectral solA: d={SOLA['d']:.5f} R={SOLA['R']:.6f})  "
        f"I={outA['I']:.4f} kappa={L_MAIN / outA['I']:.5f}")
    dAb, RAb, outAb = axi_reference(engB, rn, fn, L_MAIN)
    log(f"  big  box, L={L_MAIN}: d*={dAb:.5f}  R_A={RAb:.6f}  "
        f"I={outAb['I']:.4f} kappa={L_MAIN / outAb['I']:.5f}")
    dAc, RAc, outAc = axi_reference(eng, rn, fn, L_CTRL)
    log(f"  main box, L={L_CTRL}: d*={dAc:.5f}  R_A={RAc:.6f}  "
        f"I={outAc['I']:.4f} kappa={L_CTRL / outAc['I']:.5f}")
    axiref = dict(main=dict(d=dA, R=RA, I=outA["I"], Estat=outA["Estat"]),
                  big=dict(d=dAb, R=RAb, I=outAb["I"], Estat=outAb["Estat"]),
                  ctrl=dict(d=dAc, R=RAc, I=outAc["I"], Estat=outAc["Estat"]),
                  solA_spectral=SOLA)

    # ---- Section S: static minimization ---------------------------------
    log("\n===== Section S: static minimization (perturbed hedgehog) =====")
    q = eng.hedgehog(rn, fn)
    out_h, _ = eng.energy(q, need_grad=False)
    q = eng.normalize(q + eng.perturbation())
    out_p, _ = eng.energy(q, need_grad=False)
    log(f"  hedgehog E_static={out_h['Estat']:.7f}; perturbed "
        f"E_static={out_p['Estat']:.7f} (dE=+{out_p['Estat'] - out_h['Estat']:.5f})")
    q, outS, serS, stS, msS = anf(eng, q, None, MAXIT_STATIC, "static")
    hfS, cenS = eng.halo_fraction(q)
    virS = virial(outS)
    log(f"  static solution: Estat={outS['Estat']:.7f} "
        f"(hedgehog {out_h['Estat']:.7f}, delta {outS['Estat'] - out_h['Estat']:+.2e})")
    log("  sectors: " + " ".join(
        f"{k}={outS[k]:.6f}(vs hh {outS[k] / out_h[k] - 1:+.2e})"
        for k in ("E2", "E4", "E6", "E0", "I")))
    log(f"  degree={outS['deg']:.6f}  virial={virS:+.3e} "
        f"(rel {virS / outS['Estat']:+.3e})  centroid={cenS}")
    q_static = q.copy()
    static_res = dict(E_hedgehog=out_h["Estat"], E_perturbed=out_p["Estat"],
                      final={k: outS[k] for k in outS}, virial=virS,
                      virial_rel=virS / outS["Estat"], centroid=list(cenS),
                      meta=msS, series=serS)

    # ---- Section M: main Routhian descent at L = 8.4979 -----------------
    log(f"\n===== Section M: Routhian descent, L={L_MAIN} "
        f"(kappa_hedgehog={L_MAIN / outS['I']:.4f} > KC={KC:.5f}) =====")
    q = eng.normalize(q_static.copy() + eng.perturbation())
    q, outM, serM, stM, msM = anf(eng, q, L_MAIN, MAXIT_MAIN, "L=8.4979")
    log(f"  final: R={outM['R']:.6f} vs on-grid axi ref {RA:.6f} "
        f"(dR={outM['R'] - RA:+.6f}) vs spectral solA {SOLA['R']:.6f}")
    log(f"  I: {serM[0]['I']:.4f} -> {outM['I']:.4f};  kappa: "
        f"{L_MAIN / serM[0]['I']:.5f} -> {L_MAIN / outM['I']:.5f} "
        f"(threshold {KC:.5f});  halo: {serM[0]['halo']:.4f} -> "
        f"{serM[-1]['halo']:.4f};  deg={outM['deg']:.5f}")

    # ---- Section C: control descent at L = 4.0 ---------------------------
    log(f"\n===== Section C: control descent, L={L_CTRL} "
        f"(kappa_hedgehog={L_CTRL / outS['I']:.4f} < KC) =====")
    q = eng.normalize(q_static.copy() + eng.perturbation())
    q, outC, serC, stC, msC = anf(eng, q, L_CTRL, MAXIT_CTRL, "L=4.0")
    log(f"  final: R={outC['R']:.6f} vs on-grid axi ref {RAc:.6f} "
        f"(dR={outC['R'] - RAc:+.6f})")
    log(f"  I: {serC[0]['I']:.4f} -> {outC['I']:.4f};  kappa -> "
        f"{L_CTRL / outC['I']:.5f};  halo -> {serC[-1]['halo']:.4f};  "
        f"deg={outC['deg']:.5f}")
    q_ctrl = q.copy()

    # ---- Section B: box-size sanity (1.3x box, same N) -------------------
    log(f"\n===== Section B: box-size sanity, L={L_MAIN}, box "
        f"+-{LBOX_BIG:.2f}, h={engB.h:.5f} =====")
    qB = engB.hedgehog(rn, fn)
    outBh, _ = engB.energy(qB, need_grad=False)
    log(f"  big-box hedgehog: Estat={outBh['Estat']:.6f} "
        f"deg={outBh['deg']:.5f} (start = perturbed hedgehog, not re-relaxed)")
    qB = engB.normalize(qB + engB.perturbation())
    qB, outB, serB, stB, msB = anf(engB, qB, L_MAIN, MAXIT_BIG, "bigbox")
    log(f"  final: R={outB['R']:.6f} vs big-box axi ref {RAb:.6f} "
        f"(dR={outB['R'] - RAb:+.6f});  main-box dR was "
        f"{outM['R'] - RA:+.6f}")
    log(f"  kappa -> {L_MAIN / outB['I']:.5f};  halo -> "
        f"{serB[-1]['halo']:.4f};  deg={outB['deg']:.5f}")

    # ---- Section K: clock -------------------------------------------------
    log("\n===== Section K: clock invariant on the control solution =====")
    Estat_c, I_c = outC["Estat"], outC["I"]
    erot4 = (L_CTRL ** 2 / (2 * I_c)) / (Estat_c + L_CTRL ** 2 / (2 * I_c))
    # bisection on phi(L) = L^2 - (2/3) I Estat at frozen fields
    lo, hi = 0.0, 20.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mid * mid - (2.0 / 3.0) * I_c * Estat_c > 0:
            hi = mid
        else:
            lo = mid
    L_clock = 0.5 * (lo + hi)
    kap_clock = L_clock / I_c
    erotc = (L_clock ** 2 / (2 * I_c)) / (Estat_c + L_clock ** 2 / (2 * I_c))
    log(f"  control solution: Estat={Estat_c:.6f} I={I_c:.5f}")
    log(f"  E_rot/E_tot at L={L_CTRL}: {erot4:.6f}  (not 1/4: clock not "
        f"imposed)")
    log(f"  L_clock (bisection, frozen fields) = {L_clock:.6f}  "
        f"[closed form sqrt((2/3) I Estat) = "
        f"{np.sqrt(2.0 / 3.0 * I_c * Estat_c):.6f}]")
    log(f"  E_rot/E_tot at L_clock: {erotc:.6f}  (= 1/4 identically when the "
        f"clock holds)")
    log(f"  kappa_ours(L_clock) = {kap_clock:.6f} vs threshold {KC:.6f} -> "
        f"clock-charged solution is "
        f"{'ABOVE (over-spun regime)' if kap_clock > KC else 'below'} "
        f"the halo threshold (ratio {kap_clock / KC:.4f})")
    clock = dict(Estat=Estat_c, I=I_c, erot_frac_L4=erot4, L_clock=L_clock,
                 kappa_clock=kap_clock, erot_frac_clock=erotc,
                 above_threshold=bool(kap_clock > KC),
                 ratio_to_threshold=kap_clock / KC)

    # ---- outputs ----------------------------------------------------------
    results = dict(
        meta=dict(N=N_MAIN, LBOX=LBOX, LBOX_BIG=LBOX_BIG,
                  h=2 * LBOX / N_MAIN, t=T_FROZEN, L_main=L_MAIN,
                  L_ctrl=L_CTRL, KC=KC, seed=SEED_PERT, smoke=SMOKE,
                  halo_radius=HALO_RADIUS,
                  maxit=dict(static=MAXIT_STATIC, main=MAXIT_MAIN,
                             ctrl=MAXIT_CTRL, big=MAXIT_BIG),
                  runtime_s=None),
        radial=dict(E2=sect[0], E4=sect[1], E6=sect[2], E0=sect[3],
                    eps_ratio=ratio, gmax=gmax),
        gates=gates,
        axiref=axiref,
        static=static_res,
        descent_main=dict(L=L_MAIN, status=stM, meta=msM,
                          final={k: outM[k] for k in outM},
                          R_axi_grid=RA, dR_vs_axi=outM["R"] - RA,
                          series=serM),
        descent_ctrl=dict(L=L_CTRL, status=stC, meta=msC,
                          final={k: outC[k] for k in outC},
                          R_axi_grid=RAc, dR_vs_axi=outC["R"] - RAc,
                          series=serC),
        descent_big=dict(L=L_MAIN, status=stB, meta=msB,
                         final={k: outB[k] for k in outB},
                         R_axi_grid=RAb, dR_vs_axi=outB["R"] - RAb,
                         series=serB),
        clock=clock,
    )
    results["meta"]["runtime_s"] = time.time() - t_start
    with open(f"{OUTDIR}/field3d_results{TAG}.json", "w") as fh:
        json.dump(results, fh, indent=1)
    make_figure(serM, serB, serC,
                dict(R_axi_grid_main=RA, R_axi_grid_big=RAb),
                f"{OUTDIR}/field3d_decay{TAG}.png")
    log(f"\nwrote field3d_results{TAG}.json, field3d_decay{TAG}.png, "
        f"field3d_run{TAG}.log;  runtime {time.time() - t_start:.1f} s")
    LOGFH.close()


if __name__ == "__main__":
    main()
