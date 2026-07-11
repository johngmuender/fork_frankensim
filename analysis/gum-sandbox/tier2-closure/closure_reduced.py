#!/usr/bin/env python3
"""Tier 2a - the REDUCED closure loop (GUM replication program).

Replicates the internal consistency of the <r1> benchmark tuple
(c, kappa, g*, lambda*, V, E_rot/E at eps = 0.05) using the paper's own
reduced structure (GUM-Omega Sec. IV.C/IV.J, App. G.5 as corrected by
AUD-15 F-A15-7):

  - shape moduli enter ONLY through the inertia enhancement g(lambda)V
    (the oblate-spheroid curve g(lambda) = (3/2)[1/b^2 - (1-b^2)artanh(b)/b^3],
    b^2 = 1 - lambda^2);
  - E_field(V) = (e0/2)(1/V + V) on the BPS locus;
  - E_rot = c^2/(8 i0 g V) at spin L = c/2;
  - the eps-suppressed sectors are modeled as a penalty E_eps = eps*c_p/lambda^p
    (p = 1 reproduces the corpus's geometric exponent 2/3; p = 4/3, the naive
    slab-gradient estimate, gives 3/5 - both are run for sensitivity).

Exactly ONE constant (c_p) is calibrated, from ONE datum (lambda* = 0.42 at
eps = 0.05). Everything else - g*, c, kappa, V, the eps-scan law - is then
PREDICTED and compared against the corpus's claims.

Closure system (units Lambda*m~ for energy, Lambda*sqrt(J) for action):
  spin:    kappa = c/(2 i0 g V)
  V-stat:  V^2 = 1 + u/g,          u = c^2/(4 e0 i0)
  clock:   E_field + E_eps + E_rot = c*kappa   =>  u = g(1 + beta*V),
           V^2 - beta*V - 2 = 0,   beta = eps*c_p/(e0*lambda^p)*p_factor
  lam-stat: c^2(-g')/(8 i0 g^2 V) = eps*c_p*p/lambda^(p+1)

E_rot/E_tot = 1/4 is automatic (pure kinematics of spin+clock) - as the
paper's Cor. IV.2 claims.
"""
import math

E0 = 64.0 / (15.0 * math.pi)     # e0 = 1.3581...
I0 = 256.0 / (105.0 * math.pi)   # i0 = 0.7761...
C0 = 128.0 * math.sqrt(42.0) / (105.0 * math.pi)  # 2.5147...


def g_oblate(lam):
    """AUD-15 F-A15-7 corrected oblate branch; g(1)=1, g(0+)=3/2."""
    if lam >= 1.0:
        return 1.0
    b2 = 1.0 - lam * lam
    b = math.sqrt(b2)
    return 1.5 * (1.0 / b2 - (1.0 - b2) * math.atanh(b) / (b2 * b))


def gprime(lam, h=1e-7):
    return (g_oblate(lam + h) - g_oblate(lam - h)) / (2 * h)


def solve_at(eps, c_p, p=1.0):
    """Solve the reduced closure at eps; returns dict or None."""
    def residual(lam):
        g = g_oblate(lam)
        gp = gprime(lam)
        beta = eps * c_p / (E0 * lam ** p)
        V = 0.5 * (beta + math.sqrt(beta * beta + 8.0))
        u = g * (1.0 + beta * V)
        c2 = 4.0 * E0 * I0 * u
        lhs = c2 * (-gp) / (8.0 * I0 * g * g * V)
        rhs = eps * c_p * p / lam ** (p + 1.0)
        return lhs - rhs

    # bisection on lam in (lo, hi): residual > 0 near sphere (rotor wants
    # pancake), < 0 near pancake (penalty wins)
    lo, hi = 1e-4, 0.999
    flo, fhi = residual(lo), residual(hi)
    if flo * fhi > 0:
        return None
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        fm = residual(mid)
        if flo * fm <= 0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
    lam = 0.5 * (lo + hi)
    g = g_oblate(lam)
    beta = eps * c_p / (E0 * lam ** p)
    V = 0.5 * (beta + math.sqrt(beta * beta + 8.0))
    u = g * (1.0 + beta * V)
    c = 2.0 * math.sqrt(E0 * I0 * u)
    kappa = c / (2.0 * I0 * g * V)
    e_rot = c * c / (8.0 * I0 * g * V)
    e_field = 0.5 * E0 * (1.0 / V + V)
    e_eps = eps * c_p / lam ** p
    e_tot = e_field + e_eps + e_rot
    return dict(lam=lam, g=g, V=V, c=c, kappa=kappa,
                erot_frac=e_rot / e_tot, clock_check=c * kappa - e_tot)


def calibrate_cp(eps, lam_target, p=1.0):
    """Find c_p such that lambda*(eps) = lam_target (bisection on log c_p)."""
    def lam_of(cp):
        r = solve_at(eps, cp, p)
        return r["lam"] if r else float("nan")
    lo, hi = 1e-4, 1e4
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        lm = lam_of(mid)
        # larger c_p (stronger penalty) pushes lambda up toward sphere
        if lm < lam_target:
            lo = mid
        else:
            hi = mid
    return math.sqrt(lo * hi)


def main():
    print("== TIER 2a: REDUCED CLOSURE LOOP ==\n")

    # --- eps -> 0 limit must reproduce the exact family -------------------
    r0 = solve_at(1e-9, 1.0, 1.0)
    print(f"eps->0 limit:  c = {r0['c']:.6f} (exact 2.514754), "
          f"V = {r0['V']:.6f} (sqrt2 = 1.414214), kappa = {r0['kappa']:.6f} "
          f"(sqrt(7/12) = 0.763763), E_rot/E = {r0['erot_frac']:.6f}")

    # --- g(lambda) curve spot checks (AUD-15 V15.2) ------------------------
    print(f"\ng(0.37) = {g_oblate(0.37):.4f}, g(0.42) = {g_oblate(0.42):.4f}, "
          f"g(0.47) = {g_oblate(0.47):.4f}  [corpus: 1.28-1.32 spans measured "
          f"g* = 1.31 +/- 0.04]")

    for p, pname in [(1.0, "1/lambda (corpus-exponent form)"),
                     (4.0 / 3.0, "lambda^-4/3 (slab-gradient estimate)")]:
        print(f"\n----- penalty form E_eps = eps*c_p*{pname} -----")
        cp = calibrate_cp(0.05, 0.42, p)
        r = solve_at(0.05, cp, p)
        print(f"calibrated c_p = {cp:.4f} from the ONE datum lambda*(0.05) = 0.42")
        print(f"PREDICTIONS at eps = 0.05 vs <r1> benchmark:")
        print(f"  g*     = {r['g']:.4f}    vs 1.31  +/- 0.04")
        print(f"  c      = {r['c']:.4f}    vs 2.37  +/- 0.09")
        print(f"  kappa  = {r['kappa']:.4f}    vs 0.802 +/- 0.018")
        print(f"  V      = {r['V']:.4f}    vs 1.409 +/- 0.010  (sqrt2 = 1.4142)")
        print(f"  E_rot/E= {r['erot_frac']:.6f}  vs 0.2500 +/- 0.0002")
        print(f"  clock residual = {r['clock_check']:.2e} (must be ~0)")

        # --- eps-scan: fit  c(eps) = C0 (1 - c_g eps^nu) --------------------
        eps_list = [0.001 * (10 ** (i / 12.0)) for i in range(25)]  # 1e-3..1e-1
        rows = []
        for e in eps_list:
            rr = solve_at(e, cp, p)
            if rr and rr["c"] < C0:
                rows.append((e, rr["c"], rr["lam"]))
        # log-log fit of deficit d = 1 - c/C0 vs eps
        xs = [math.log(e) for e, c, _ in rows]
        ys = [math.log(1.0 - c / C0) for e, c, _ in rows]
        n = len(xs)
        sx, sy = sum(xs), sum(ys)
        sxx = sum(x * x for x in xs)
        sxy = sum(x * y for x, y in zip(xs, ys))
        nu = (n * sxy - sx * sy) / (n * sxx - sx * sx)
        lnA = (sy - nu * sx) / n
        cg_fit = math.exp(lnA)
        print(f"  eps-scan fit:  1 - c/c0 = {cg_fit:.3f} * eps^{nu:.3f}")
        print(f"     corpus law: 1 - c/c0 = c_g eps^(2/3), c_g = 0.42 +/- 0.04 "
              f"(with O(eps^(1/3) ln eps) corrections acknowledged)")
        # local c_g at the benchmark point assuming exponent 2/3 exactly
        e_b, c_b, _ = min(rows, key=lambda t: abs(t[0] - 0.05))
        cg_local = (1.0 - c_b / C0) / (e_b ** (2.0 / 3.0))
        print(f"  local c_g at eps ~ 0.05 (exponent pinned to 2/3): {cg_local:.3f}")

        with open(f"epsscan_p{p:.2f}.csv", "w") as f:
            f.write("eps,c,lambda\n")
            for e, c, lam in rows:
                f.write(f"{e:.6e},{c:.8f},{lam:.6f}\n")

    print("\nEPISTEMIC NOTICE: reduced-model replication (one calibrated "
          "constant).\nThe full 2-D field-theoretic solve remains future work; "
          "PASS here verifies the\ncorpus's INTERNAL consistency loop "
          "(lambda*, g*, c, kappa, V), not the field theory.")


if __name__ == "__main__":
    main()
