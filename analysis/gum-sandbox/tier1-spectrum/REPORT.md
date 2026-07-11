# Tier-1 linear-spectrum pilot - GUM replication program

Run date: 2026-07-11. Script: `spectrum.py` (deterministic; seed 20260711 set, no randomness consumed). Runtime: 5.8 s. numpy-only solver (M diagonal, standard Hermitian eigenproblem on M^-1/2 K M^-1/2; scipy not installed).

Benchmark moduli: rho0=1, J=1, lambda=1, mu=1, mu_c=5, alpha=0.5, beta=0.5, gamma=0.5, m_V=1. k grid: 400 log-spaced points in [1e-4, 3].

## Global sanity

- Most negative omega^2 encountered anywhere (all runs, all k): 1.000e-08  (tolerance -1e-12) -> **PASS**
- Exact-decoupling cross-check (union of u_z / phi_z / transverse sub-blocks vs full 6x6 spectrum): max |Delta omega^2| = 3.411e-13 -> PASS
- Hermiticity of K including chi3 term: max |K - K^dag| = 0.000e+00

## D1 - Mode count and classification

6 branches at every k. Classification of the full 6x6 eigenvectors (dominant-sector weight in parentheses):

- k = 0.001: omega^2=1e-06 [transverse, 1.000000], omega^2=1e-06 [transverse, 1.000000], omega^2=3e-06 [long-u (u_z), 1.000000], omega^2=21 [long-phi (phi_z), 1.000000], omega^2=21 [transverse, 1.000000], omega^2=21 [transverse, 1.000000]
- k = 0.1: omega^2=0.009988 [transverse, 1.000000], omega^2=0.009988 [transverse, 1.000000], omega^2=0.03 [long-u (u_z), 1.000000], omega^2=21.01 [long-phi (phi_z), 1.000000], omega^2=21.06 [transverse, 1.000000], omega^2=21.06 [transverse, 1.000000]
- k = 1: omega^2=0.8985 [transverse, 1.000000], omega^2=0.8985 [transverse, 1.000000], omega^2=3 [long-u (u_z), 1.000000], omega^2=22 [long-phi (phi_z), 1.000000], omega^2=26.85 [transverse, 1.000000], omega^2=26.85 [transverse, 1.000000]
- k = 3: omega^2=5.823 [transverse, 1.000000], omega^2=5.823 [transverse, 1.000000], omega^2=27 [long-u (u_z), 1.000000], omega^2=30 [long-phi (phi_z), 1.000000], omega^2=75.93 [transverse, 1.000000], omega^2=75.93 [transverse, 1.000000]

Counts at every probed k: 1 longitudinal-u (B1), 1 longitudinal-phi (B4), 4 transverse (B2 light doublet + B3 heavy doublet) -> **PASS**

Labels: B1 = longitudinal acoustic (u_z), B2 = light transverse doublet, B3 = heavy transverse doublet, B4 = longitudinal twist (phi_z). Note B3 and B4 are degenerate at k = 0 (same gap), so at k -> 0 the gapped level is triply degenerate; they separate at O(k^2).

## D2 - B1 acoustic speed (Case A, m_V = 1)

- measured c_L = omega/k (k->0 fit over k <= 0.001): 1.732050807569
- raw omega/k at k = 1e-4: 1.732050807569
- expected c_L = sqrt((lambda+2mu)/rho0) = sqrt(3) = 1.732050807569
- relative error: 1.282e-16 -> **PASS**

## D3 - B2 gaplessness and speed (Case A, m_V = 1)

- extrapolated gap omega^2(k->0) from linear fit in k^2: 7.862e-15 (expected ~0) 
- small-k speed c_B2 = 0.999999950162 (compare sqrt(mu/rho0) = 1.000000000000)
- intra-doublet splitting of B2 across all k (chi3 = 0): 2.842e-14
- verdict: **PASS**

## D4 - B3 Klein-Gordon fit (Case A, m_V = 1)

- fit omega^2 = omega0^2 + c_psi^2 k^2 over k <= 0.1:
  - measured omega0^2 = 20.999999822
  - measured c_psi^2 = 5.750872825  (c_psi = 2.398098)
- candidate formula (m_V^2 + 4 mu_c)/J = 21.000000000 -> relative error 8.465e-09 -> **PASS**
- measured coefficient scan (gap vs moduli): d(omega0^2)/d(mu_c) = 4.000000000 (candidate: 4), d(omega0^2)/d(m_V^2) = 1.000000000 (candidate: 1). The measured relation is omega0^2 = (4 mu_c + m_V^2)/J with the factor 4 confirmed under this paper's conventions.
- intra-doublet splitting of B3 across all k: 1.421e-14

## D5 - B4 gap (Case A, m_V = 1)

- measured omega0^2 = 21.000000000, small-k curvature c^2 = 1.000000000 (analytic slope alpha+beta = 1)
- equals the B3 gap (4 mu_c + m_V^2)/J = 21 -> **PASS**

## D6 - HEADLINE: Theorem II.2 numerically (m_V in {0, 1, 5, 20})

Light-doublet dispersion omega_B2(k) for both cases; full table in `theorem_II2_invariance.csv`.

### Case A (objective mass on psi)

| m_V | gap omega^2(0) | c_B2^2 (k->0) | max_k |omega/omega(m_V=0) - 1| | r(k=1e-4) | r(k=1e-2) |
|---|---|---|---|---|---|
| 0 | 7.74e-15 | 0.999999900 | 0.000e+00 | 1.000000000 | 1.000002469 |
| 1 | 7.86e-15 | 0.999999900 | 2.572e-04 | 1.000000000 | 1.000002352 |
| 5 | 8.02e-15 | 0.999999900 | 3.040e-03 | 1.000000000 | 1.000001097 |
| 20 | 7.51e-15 | 0.999999895 | 5.263e-03 | 1.000000000 | 1.000000118 |

- max relative deviation of omega_B2 across the m_V sweep, all k in [1e-4, 3]: **5.263e-03**
- same restricted to k <= 0.01: 1.730e-06
- gap invariance: gap = 0 for every m_V to 8.0e-15 -> PASS
- IR speed invariance: max |c_B2(m_V) - c_B2(0)| = 2.659e-09

- **Literal spec check (max relative deviation < 1e-12): FAIL** (measured 5.263e-03).

  **Honest deviation note.** Pointwise invariance of the full omega_B2(k) curve under m_V does NOT hold at finite k in this implementation, and the residual is structural, not numerical noise. In Case A, m_V^2 enters the energy only through |psi|^2, i.e. only in the combination A_c = 2 mu_c + m_V^2/2 multiplying the relative-rotation term (the identity e_[ij] e_[ij] = 2|psi|^2 makes the mu_c and m_V terms the same operator). The light branch has a small but nonzero psi admixture at finite k (psi ~ k^3/(4 A_c) x u), so its frequency acquires an O(k^4/A_c) dependence on m_V: analytically, at small k, omega^2 = (mu/rho0) k^2 [1 - ((beta+gamma)/8 + J mu/(2 rho0)) k^2 / A_c + ...] for this benchmark, and the curve flows toward the psi = 0 (locked MacCullagh) limit omega^2 = k^2 (8 mu + (beta+gamma) k^2)/(8 rho0 + 2 J k^2) as m_V -> inf. The measured max deviation (~5e-3 in omega, at k = 3, between m_V = 0 and 20) matches this closed form. What IS invariant to machine precision: the gap (identically 0 for all m_V) and the k -> 0 dispersion (speed and co-motion). So the *infrared* content of Theorem II.2 is reproduced exactly; the <1e-12 pointwise-in-k claim is not, and the deviation decays as k^2 relative (k^4 absolute) toward the IR.

### Case B (M-1 defect: mass on phi)

| m_V | gap omega^2(0) | c_B2^2 (k->0) | analytic c^2 = mu + mu_c m_V^2/(4mu_c+m_V^2) | max_k |omega/omega(m_V=0) - 1| | r(k=1e-4) | r(k=1e-2) |
|---|---|---|---|---|---|---|
| 0 | 7.74e-15 | 0.999999900 | 1.000000000 | 0.000e+00 | 1.000000000 | 1.000002469 |
| 1 | 1.02e-14 | 1.238095105 | 1.238095238 | 1.127e-01 | 0.952380953 | 0.952384259 |
| 5 | 9.73e-15 | 3.777777648 | 3.777777778 | 1.048e+00 | 0.444444445 | 0.444447642 |
| 20 | -1.23e-15 | 5.761904765 | 5.761904762 | 1.977e+00 | 0.047619048 | 0.047619107 |

- **The Case B light branch SHIFTS, as required: max relative deviation 1.977e+00** (i.e. up to 198% in omega at m_V = 20).
- HOW it shifts: no gap opens (B2 stays exactly gapless in Case B too, because at k = 0 the displacement u drops out of every energy term, leaving three zero modes); instead the light-branch SPEED changes: c_B2^2 rises from mu = 1 at m_V = 0 toward mu + mu_c = 6 as m_V -> inf, following the measured values above (they match the analytic 2x2 result c^2 = mu + mu_c m_V^2/(4 mu_c + m_V^2) to fit precision).
- Co-motion fidelity r(k->0) = |phi|/((k/2)|u|): Case A locks r = 1 for every m_V (the MacCullagh wave keeps its rotational content); Case B suppresses it as r = 4 mu_c/(4 mu_c + m_V^2): measured r(0) = 1.000000, r(1) = 0.952381, r(5) = 0.444444, r(20) = 0.047619 vs analytic 1.000000, 0.952381, 0.444444, 0.047619. This is the numerical content of 'the M-1 defect kills the photon's rotational content'.

## D7 - Stretch: tree achirality with chi3 = 0.3 (Case A, m_V = 1)

- chiral term implemented as W_chi3 = chi3 Re(e_[ij]* Gamma_[ij]), i.e. K += chi3 (R_ea^dag R_Ga + R_Ga^dag R_ea); Hermiticity residual 0.0e+00.
- max splitting of the two light transverse (circular-polarization) branches over all k: **max |omega_+ - omega_-| = 3.415785e-02**
- splitting at sample k: k=0.01: 3.514e-11, k=0.1: 3.458e-07, k=1: 2.399e-03, k=3: 3.416e-02
- max heavy-doublet (B3) splitting: 1.971132e-01
- **Literal spec check (splitting at machine precision): FAIL** (measured 3.416e-02).

  **Honest deviation note.** The two circular polarizations of the LIGHT doublet are NOT exactly degenerate at finite k under this term; the splitting is real but strongly IR-suppressed: a log-log fit of Delta(omega^2) vs k over k in [0.02, 0.2] gives exponent 4.994 (prefactor exp(-4.95) = 0.007), i.e. Delta(omega^2) ~ k^5. Analytically this is because W_chi3 = -chi3 Re(psi* . curl phi) vanishes exactly on the locked (psi = 0) wave; the light branch only violates psi = 0 at O(k^3/A_c), giving Delta(omega^2) ~ chi3 k^5/A_c — zero at tree level in the IR sense (achirality is exact as k -> 0 and exact on the constrained MacCullagh wave), but not 'machine precision at every k'. The HEAVY doublet splits at O(chi3 k), which is the expected chiral response of the psi-carrying branch. If the paper's Theorem meant exact all-k degeneracy, this implementation does not reproduce it; if it meant the physical (IR / on-shell photon) statement, it does.

## Deliverables

- `spectrum.py` - this script (single file, reproducible)
- `dispersion_caseA.png`, `dispersion_caseB.png` - all 6 branches, log-log inset
- `theorem_II2_invariance.csv` - omega_B2(k) for both cases x m_V in {0,1,5,20}
- `REPORT.md` - this report

## Verdict summary

| Diagnostic | Verdict |
|---|---|
| Global: all omega^2 >= -1e-12 | PASS (min 1.0e-08) |
| D1 mode count / classification | PASS |
| D2 c_L = sqrt(3) | PASS |
| D3 B2 gapless, c_B2 = 1 | PASS |
| D4 B3 gap = (m_V^2 + 4 mu_c)/J | PASS (factor 4 confirmed) |
| D5 B4 gap | PASS |
| D6 Case A invariance (literal < 1e-12) | FAIL (measured 5.3e-03; IR statement exact) |
| D6 Case A gap/speed/r invariance (IR) | PASS |
| D6 Case B shifts (speed, no gap) | PASS (up to 198% at m_V=20) |
| D7 light-doublet achirality (literal machine-eps) | FAIL (measured 3.4e-02, ~k^5 suppressed, exact as k->0) |
