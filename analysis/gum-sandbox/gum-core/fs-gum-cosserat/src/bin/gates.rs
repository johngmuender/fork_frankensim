//! GUM physics core, Phase B3 — the Tier-1 referee battery for
//! fs-gum-cosserat (gates binary).
//!
//! Every gate below is bound to a number the validated Python pilot
//! (analysis/gum-sandbox/tier1-spectrum/spectrum.py + REPORT.md +
//! theorem_II2_invariance.csv) already measured, plus the two pieces that
//! are NEW in this crate: arbitrary-k rotation invariance and the
//! time-domain symplectic validation (the repo's first).
//!
//! Certified tier (fs-cosserat-pilot pattern, domain "gum-core:cosserat:v1"):
//! passing gates become fail-closed Certified claims in an fs-package
//! EvidencePackage; fs-checker verifies the Merkle root in three modes; the
//! whole battery runs twice in-process and must reproduce the root
//! bit-for-bit (jacobi_eigh and every reduction here are deterministic).
//!
//! Epistemic notice (binding): every PASS certifies a WITHIN-MODEL property
//! of a speculative theory's linearized dispersion — it validates the port,
//! the eigensolve path, and the integrator, never the physics.

use fs_blake3::ContentHash;
use fs_evidence::{Evidence, ProvenanceHash};
use fs_gum_cosserat::branches::{branches, classify_full, fit_loglog, fit_w2};
use fs_gum_cosserat::eigh::{oracle_check, solve};
use fs_gum_cosserat::moduli::{MassCase, Moduli};
use fs_gum_cosserat::softsector::SoftSector;
use fs_gum_cosserat::symbol::{mass_matrix, stiffness};
use fs_gum_cosserat::verlet::evolve_mode;
use fs_gum_cosserat::wchi::{energy_displaced, grad_dot, wchi_energy, wchi_grad, MicroField};
use fs_math::c64::C64;
use fs_package::origin::{
    SourceCertificateRequest, SourceCertificateVerifier, VerificationCapabilities,
    VerificationDecision,
};
use fs_package::{Claim, EvidencePackage, Provenance};

const EPS64: f64 = f64::EPSILON;
const SQRT3: f64 = 1.732_050_807_568_877_2;
const K_FIT_ACOUSTIC: f64 = 1e-3;
const K_FIT_GAPPED: f64 = 0.1;
const MV_LIST: [f64; 4] = [0.0, 1.0, 5.0, 20.0];
const TD_STEPS: usize = 200_000;

// ---------------------------------------------------------------------------
// Certificate policy (fs-cosserat-pilot pattern, new domain)
// ---------------------------------------------------------------------------

fn certificate_hash(id: &str, statement: &str, lo: f64, hi: f64) -> ContentHash {
    let canon = format!("{id}|{statement}|{:016x}|{:016x}", lo.to_bits(), hi.to_bits());
    fs_blake3::hash_domain("gum-core:cosserat:v1", canon.as_bytes())
}

struct GatesCertVerifier;

impl SourceCertificateVerifier for GatesCertVerifier {
    fn verify(&self, req: &SourceCertificateRequest<'_>) -> VerificationDecision {
        let policy =
            fs_blake3::hash_domain("gum-core:cosserat:policy:v1", b"recompute-and-compare");
        let expect = certificate_hash(req.claim_id, req.statement, req.lo, req.hi);
        if req.producer == "fs-gum-cosserat" && req.certificate_hash == expect {
            VerificationDecision::accept(policy)
        } else {
            VerificationDecision::reject(policy)
        }
    }
}

// ---------------------------------------------------------------------------
// Gate bookkeeping
// ---------------------------------------------------------------------------

struct Gate {
    id: &'static str,
    statement: String,
    measured: f64,
    target: f64,
    pass: bool,
    note: String,
}

fn gate(
    id: &'static str,
    statement: &str,
    measured: f64,
    target: f64,
    pass: bool,
    note: String,
) -> Gate {
    Gate { id, statement: statement.to_string(), measured, target, pass, note }
}

// ---------------------------------------------------------------------------
// Deterministic direction sampler (splitmix-style LCG; no external RNG)
// ---------------------------------------------------------------------------

fn next_unit(state: &mut u64) -> f64 {
    *state = state.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
    ((*state >> 11) as f64) / ((1u64 << 53) as f64)
}

fn rand_dir(state: &mut u64) -> [f64; 3] {
    loop {
        let x = 2.0 * next_unit(state) - 1.0;
        let y = 2.0 * next_unit(state) - 1.0;
        let z = 2.0 * next_unit(state) - 1.0;
        let n2 = x * x + y * y + z * z;
        if n2 > 0.04 && n2 <= 1.0 {
            let n = n2.sqrt();
            return [x / n, y / n, z / n];
        }
    }
}

// ---------------------------------------------------------------------------
// CSV referee data
// ---------------------------------------------------------------------------

struct CsvData {
    kk: Vec<f64>,
    /// omega_B2 columns in order (A,0), (A,1), (A,5), (A,20), (B,0), (B,1), (B,5), (B,20).
    omega: [Vec<f64>; 8],
}

fn read_csv() -> CsvData {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tier1-spectrum/theorem_II2_invariance.csv");
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("cannot read {}: {e}", path.display()));
    let mut kk = Vec::with_capacity(400);
    let mut omega: [Vec<f64>; 8] = Default::default();
    for (li, line) in text.lines().enumerate() {
        if li == 0 || line.trim().is_empty() {
            continue;
        }
        let mut it = line.split(',');
        kk.push(it.next().expect("k col").trim().parse::<f64>().expect("k parse"));
        for col in &mut omega {
            col.push(it.next().expect("omega col").trim().parse::<f64>().expect("omega parse"));
        }
    }
    assert_eq!(kk.len(), 400, "expected the 400-point Tier-1 grid");
    CsvData { kk, omega }
}

// ---------------------------------------------------------------------------
// Sweep bookkeeping (spectrum.py Diagnostic 6)
// ---------------------------------------------------------------------------

struct Sweep {
    om: Vec<f64>,
    w2mean: Vec<f64>,
    gap: f64,
    c2: f64,
    r: Vec<f64>,
}

fn with_case(case: MassCase, m_v: f64) -> Moduli {
    let mut m = Moduli::bench();
    m.case = case;
    m.m_v = m_v;
    m
}

struct Globals {
    min_w2: f64,
    xval: f64,
    pair_resid: f64,
}

impl Globals {
    fn absorb(&mut self, b: &fs_gum_cosserat::branches::Branches) {
        self.min_w2 = self.min_w2.min(b.min_w2);
        self.xval = self.xval.max(b.xval);
        self.pair_resid = self.pair_resid.max(b.pair_resid);
    }
}

fn run_sweep(m: &Moduli, kk: &[f64], g: &mut Globals) -> Sweep {
    let b = branches(m, kk);
    g.absorb(&b);
    let w2mean: Vec<f64> = b.b2.iter().map(|p| 0.5 * (p[0] + p[1])).collect();
    let om: Vec<f64> = w2mean.iter().map(|&w| w.max(0.0).sqrt()).collect();
    let (gap, c2) = fit_w2(kk, &w2mean, K_FIT_ACOUSTIC);
    Sweep { om, w2mean, gap, c2, r: b.r_light }
}

// ---------------------------------------------------------------------------
// One full deterministic battery pass
// ---------------------------------------------------------------------------

struct RunOutput {
    gates: Vec<Gate>,
    certified: usize,
    root_hex: String,
    mode_deny: bool,
    mode_cert: bool,
    mode_tamper: bool,
    report: Vec<String>,
}

#[allow(clippy::too_many_lines)]
fn run(csv: &CsvData) -> RunOutput {
    let kk = &csv.kk;
    let mut gates: Vec<Gate> = Vec::new();
    let mut report: Vec<String> = Vec::new();
    let mut glob = Globals { min_w2: f64::INFINITY, xval: 0.0, pair_resid: 0.0 };

    let bench = Moduli::bench();
    let md = mass_matrix(&bench);

    // ================= benchmark Case A run =================
    let b_a = branches(&bench, kk);
    glob.absorb(&b_a);
    let w2_b2: Vec<f64> = b_a.b2.iter().map(|p| 0.5 * (p[0] + p[1])).collect();
    let w2_b3: Vec<f64> = b_a.b3.iter().map(|p| 0.5 * (p[0] + p[1])).collect();
    let b2_intra = b_a.b2.iter().map(|p| (p[1] - p[0]).abs()).fold(0.0, f64::max);
    let b3_intra = b_a.b3.iter().map(|p| (p[1] - p[0]).abs()).fold(0.0, f64::max);

    // ---- G01/G02: B1 longitudinal acoustic speed c_L = sqrt(3) ----
    let (_g1, s1) = fit_w2(kk, &b_a.b1, K_FIT_ACOUSTIC);
    let c_l = s1.sqrt();
    let c_l_raw = b_a.b1[0].max(0.0).sqrt() / kk[0];
    let rel = ((c_l - SQRT3) / SQRT3).abs();
    gates.push(gate(
        "G01",
        "B1 acoustic speed c_L = sqrt((lambda+2mu)/rho0) = sqrt(3) to 1e-12 (small-k fit)",
        c_l,
        SQRT3,
        rel < 1e-12,
        format!("c_L = {c_l:.12}, rel err {rel:.3e}"),
    ));
    let rel_raw = ((c_l_raw - SQRT3) / SQRT3).abs();
    gates.push(gate(
        "G02",
        "raw omega/k at k = 1e-4 equals sqrt(3) to 1e-12",
        c_l_raw,
        SQRT3,
        rel_raw < 1e-12,
        format!("omega/k = {c_l_raw:.12}, rel err {rel_raw:.3e}"),
    ));

    // ---- G03: mode census / classification (D1) ----
    let mut census_ok = true;
    let mut census_note = String::new();
    for kc in [1e-3, 0.1, 1.0, 3.0] {
        let rows = classify_full(kc, [0.0, 0.0, 1.0], &bench);
        let n_u = rows.iter().filter(|r| r.1 == "long-u").count();
        let n_p = rows.iter().filter(|r| r.1 == "long-phi").count();
        let n_t = rows.iter().filter(|r| r.1 == "transverse").count();
        let wmin = rows.iter().map(|r| r.2).fold(1.0, f64::min);
        census_ok &= n_u == 1 && n_p == 1 && n_t == 4 && wmin > 0.99;
        census_note.push_str(&format!("k={kc}: {n_u}u/{n_p}phi/{n_t}t wmin={wmin:.6}; "));
    }
    gates.push(gate(
        "G03",
        "mode census at k in (1e-3, 0.1, 1, 3): 1 long-u + 1 long-phi + 4 transverse, dominant weight > 0.99",
        if census_ok { 1.0 } else { 0.0 },
        1.0,
        census_ok,
        census_note,
    ));

    // ---- G05: B2 gaplessness and speed (bench) ----
    let (g2, s2) = fit_w2(kk, &w2_b2, K_FIT_ACOUSTIC);
    gates.push(gate(
        "G05",
        "B2 light-doublet speed c_B2^2 = mu/rho0 = 1 within 1e-6 (bench, small-k fit)",
        s2,
        1.0,
        (s2 - 1.0).abs() < 1e-6,
        format!("c_B2^2 = {s2:.9}, gap = {g2:.3e}, intra-doublet split (chi3=0) {b2_intra:.3e}"),
    ));

    // ---- G06/G07: B3 Klein-Gordon gap and the factor-4 dial ----
    let (g3, s3) = fit_w2(kk, &w2_b3, K_FIT_GAPPED);
    let w0sq = (bench.m_v * bench.m_v + 4.0 * bench.mu_c) / bench.j;
    let rel3 = ((g3 - w0sq) / w0sq).abs();
    gates.push(gate(
        "G06",
        "B3 Klein-Gordon gap omega0^2 = (m_V^2 + 4 mu_c)/J = 21 within rel 1e-6",
        g3,
        w0sq,
        rel3 < 1e-6,
        format!("omega0^2 = {g3:.9}, rel err {rel3:.3e}, c_psi^2 = {s3:.9}"),
    ));
    let gap_b3 = |mu_c: f64, m_v: f64, g: &mut Globals| -> f64 {
        let mut p = Moduli::bench();
        p.mu_c = mu_c;
        p.m_v = m_v;
        let grid: Vec<f64> = kk.iter().copied().filter(|&k| k <= K_FIT_GAPPED).collect();
        let b = branches(&p, &grid);
        g.absorb(&b);
        let mean: Vec<f64> = b.b3.iter().map(|q| 0.5 * (q[0] + q[1])).collect();
        let (gg, _) = fit_w2(&grid, &mean, K_FIT_GAPPED);
        gg
    };
    let g_50 = gap_b3(5.0, 0.0, &mut glob);
    let g_60 = gap_b3(6.0, 0.0, &mut glob);
    let g_52 = gap_b3(5.0, 2.0, &mut glob);
    let coef_mu_c = g_60 - g_50;
    let coef_mv2 = (g_52 - g_50) / 4.0;
    gates.push(gate(
        "G07",
        "factor-4 dial: d(omega0^2)/d(mu_c) = 4 and d(omega0^2)/d(m_V^2) = 1 within 1e-6",
        coef_mu_c,
        4.0,
        (coef_mu_c - 4.0).abs() < 1e-6 && (coef_mv2 - 1.0).abs() < 1e-6,
        format!("d/d(mu_c) = {coef_mu_c:.9}, d/d(m_V^2) = {coef_mv2:.9}"),
    ));

    // ---- G08: B4 longitudinal twist gap + curvature ----
    let (g4, s4) = fit_w2(kk, &b_a.b4, K_FIT_GAPPED);
    let rel4 = ((g4 - w0sq) / w0sq).abs();
    let curv = (s4 - (bench.alpha + bench.beta)).abs();
    gates.push(gate(
        "G08",
        "B4 twist gap = (m_V^2 + 4 mu_c)/J = 21 (rel 1e-6) with small-k curvature alpha + beta = 1 (1e-6)",
        g4,
        w0sq,
        rel4 < 1e-6 && curv < 1e-6,
        format!("gap = {g4:.9} (rel {rel4:.3e}), curvature = {s4:.9}"),
    ));

    // ================= D6 sweep: both cases x m_V ================
    let mut sweeps_a: Vec<Sweep> = Vec::new();
    let mut sweeps_b: Vec<Sweep> = Vec::new();
    for &mv in &MV_LIST {
        sweeps_a.push(run_sweep(&with_case(MassCase::A, mv), kk, &mut glob));
        sweeps_b.push(run_sweep(&with_case(MassCase::B, mv), kk, &mut glob));
    }

    // ---- G04: B2 gap <= 1e-12, Case A, all m_V ----
    let max_gap_a = sweeps_a.iter().map(|s| s.gap.abs()).fold(0.0, f64::max);
    gates.push(gate(
        "G04",
        "B2 stays exactly gapless in Case A: extrapolated gap below 1e-12 for m_V in (0, 1, 5, 20)",
        max_gap_a,
        0.0,
        max_gap_a <= 1e-12,
        format!(
            "gaps = {:?}",
            sweeps_a.iter().map(|s| format!("{:.2e}", s.gap)).collect::<Vec<_>>()
        ),
    ));

    // ---- G11: Case A IR invariance + structural finite-k drift ----
    let dev_case = |sw: &[Sweep], mask_kmax: f64| -> Vec<f64> {
        let base = &sw[0].om;
        sw.iter()
            .map(|s| {
                let mut d = 0.0f64;
                for i in 0..kk.len() {
                    if kk[i] <= mask_kmax {
                        let den = if base[i] > 0.0 { base[i] } else { 1.0 };
                        d = d.max(((s.om[i] - base[i]) / den).abs());
                    }
                }
                d
            })
            .collect()
    };
    let dev_a_full = dev_case(&sweeps_a, f64::INFINITY);
    let dev_a_small = dev_case(&sweeps_a, 0.01);
    let dev_b_full = dev_case(&sweeps_b, f64::INFINITY);
    let speed_dev_a = sweeps_a
        .iter()
        .map(|s| (s.c2.max(0.0).sqrt() - sweeps_a[0].c2.max(0.0).sqrt()).abs())
        .fold(0.0, f64::max);
    let drift_targets = [0.0, 2.572e-4, 3.040e-3, 5.263e-3]; // REPORT.md D6/A table
    let mut drift_ok = true;
    let mut drift_note = String::new();
    for (i, s) in dev_a_full.iter().enumerate() {
        drift_ok &= (s - drift_targets[i]).abs() < 1e-6;
        drift_note.push_str(&format!("mV={}: {s:.4e}; ", MV_LIST[i]));
    }
    let ir_dev = dev_a_small.iter().fold(0.0f64, |a, &b| a.max(b));
    gates.push(gate(
        "G11",
        "Case A Theorem II.2: IR invariance (dev(k<=0.01) <= 3e-6, speed dev <= 1e-8) and the structural finite-k drift matches REPORT.md (2.572e-4, 3.040e-3, 5.263e-3 within 1e-6)",
        dev_a_full[3],
        5.263e-3,
        ir_dev <= 3e-6 && speed_dev_a <= 1e-8 && drift_ok,
        format!("{drift_note}dev(k<=0.01) = {ir_dev:.3e}, speed dev = {speed_dev_a:.3e}"),
    ));

    // ---- G09: Case B speeds ----
    let c2_analytic = |mv: f64| -> f64 {
        let (mu, mu_c) = (bench.mu, bench.mu_c);
        if 4.0 * mu_c + mv * mv > 0.0 {
            mu + mu_c - 4.0 * mu_c * mu_c / (4.0 * mu_c + mv * mv)
        } else {
            mu
        }
    };
    let mut cb_ok = true;
    let mut cb_note = String::new();
    for (i, &mv) in MV_LIST.iter().enumerate() {
        let want = c2_analytic(mv);
        let got = sweeps_b[i].c2;
        cb_ok &= (got - want).abs() < 1e-6 && sweeps_b[i].gap.abs() <= 1e-12;
        cb_note.push_str(&format!("mV={mv}: c^2 = {got:.9} (analytic {want:.9}); "));
    }
    gates.push(gate(
        "G09",
        "Case B M-1 defect: no gap opens (gap <= 1e-12) and c_B2^2 = mu + mu_c m_V^2/(4mu_c + m_V^2) = (1.238095238, 3.777777778, 5.761904762) within 1e-6",
        sweeps_b[3].c2,
        c2_analytic(20.0),
        cb_ok,
        cb_note,
    ));

    // ---- G12: Case B branch shifts ----
    let dev_b = dev_b_full.iter().fold(0.0f64, |a, &b| a.max(b));
    gates.push(gate(
        "G12",
        "Case B light branch SHIFTS: max relative deviation across the m_V sweep = 1.977 within 1e-3",
        dev_b,
        1.977,
        (dev_b - 1.977).abs() < 1e-3,
        format!(
            "per-mV devs = {:?}",
            dev_b_full.iter().map(|d| format!("{d:.4e}")).collect::<Vec<_>>()
        ),
    ));

    // ---- G10: co-motion fidelity r(k->0) ----
    let mut r_ok = true;
    let mut r_note = String::new();
    for (i, &mv) in MV_LIST.iter().enumerate() {
        let r_a = sweeps_a[i].r[0];
        let r_b = sweeps_b[i].r[0];
        let want_b = 4.0 * bench.mu_c / (4.0 * bench.mu_c + mv * mv);
        r_ok &= (r_a - 1.0).abs() < 5e-9 && (r_b - want_b).abs() < 5e-9;
        r_note.push_str(&format!("mV={mv}: rA = {r_a:.9}, rB = {r_b:.9} (want {want_b:.9}); "));
    }
    gates.push(gate(
        "G10",
        "co-motion fidelity r(k->0): Case A locks r = 1; Case B gives r = 4mu_c/(4mu_c + m_V^2) — both to 9 digits (5e-9)",
        sweeps_b[3].r[0],
        4.0 * bench.mu_c / (4.0 * bench.mu_c + 400.0),
        r_ok,
        r_note,
    ));

    // ================= D7: chi3 chirality ================
    let mut m_chi = Moduli::bench();
    m_chi.chi3 = 0.3;
    let b_chi = branches(&m_chi, kk);
    glob.absorb(&b_chi);
    let om_split: Vec<f64> = b_chi
        .b2
        .iter()
        .map(|p| (p[1].max(0.0).sqrt() - p[0].max(0.0).sqrt()).abs())
        .collect();
    let max_split_light = om_split.iter().fold(0.0f64, |a, &b| a.max(b));
    let w2_split: Vec<f64> = b_chi.b2.iter().map(|p| (p[1] - p[0]).abs()).collect();
    let (p_exp, _) = fit_loglog(kk, &w2_split, |i| {
        kk[i] >= 0.02 && kk[i] <= 0.2 && w2_split[i] > 0.0
    });
    let max_split_heavy = b_chi
        .b3
        .iter()
        .map(|p| (p[1].max(0.0).sqrt() - p[0].max(0.0).sqrt()).abs())
        .fold(0.0f64, f64::max);
    gates.push(gate(
        "G13",
        "chi3 = 0.3 chirality: light-doublet splitting scales as k^5 (fit exponent within 0.1 of 5), max light split = 3.415785e-2 and max heavy split = 1.971132e-1 within 1e-6; b2 intra-doublet split at chi3 = 0 below 1e-12",
        p_exp,
        5.0,
        (p_exp - 5.0).abs() < 0.1
            && (max_split_light - 3.415785e-2).abs() < 1e-6
            && (max_split_heavy - 1.971132e-1).abs() < 1e-6
            && b2_intra <= 1e-12,
        format!(
            "exponent = {p_exp:.3}, max light split = {max_split_light:.6e}, max heavy split = {max_split_heavy:.6e}, b3 intra (chi3=0) = {b3_intra:.2e}"
        ),
    ));

    // ================= rotation invariance (NEW: arbitrary k) ================
    let mut rot_dev = 0.0f64;
    let mut seed: u64 = 20_260_711;
    for m in [&bench, &m_chi] {
        for kmag in [1e-3, 0.1, 1.0, 3.0] {
            let kz = stiffness([0.0, 0.0, kmag], m);
            let wz = solve(&kz, &md, None);
            let scale = wz.w2[5].abs().max(1.0);
            for _ in 0..8 {
                let d = rand_dir(&mut seed);
                let kg = stiffness([kmag * d[0], kmag * d[1], kmag * d[2]], m);
                let wg = solve(&kg, &md, None);
                for i in 0..6 {
                    rot_dev = rot_dev.max((wg.w2[i] - wz.w2[i]).abs() / scale);
                }
            }
        }
    }
    gates.push(gate(
        "G16",
        "rotation invariance: 8 random k directions x |k| in (1e-3, 0.1, 1, 3), chi3 in (0, 0.3) — eigenvalues match the z-axis case to 1e-12 (scaled)",
        rot_dev,
        0.0,
        rot_dev <= 1e-12,
        format!("max scaled eigenvalue deviation = {rot_dev:.3e}"),
    ));

    // ================= oracle + embedding pairing ================
    let mut orc_dev = 0.0f64;
    let mut orc_im = 0.0f64;
    let mut seed2: u64 = 987_654_321;
    for m in [&bench, &m_chi] {
        for kmag in [1e-3, 1.0, 3.0] {
            for dir in [[0.0, 0.0, 1.0], rand_dir(&mut seed2)] {
                let kmat = stiffness([kmag * dir[0], kmag * dir[1], kmag * dir[2]], m);
                let modes = solve(&kmat, &md, None);
                let scale = modes.w2[5].abs().max(1.0);
                let (dre, dim) = oracle_check(&kmat, &md, &modes.w2);
                orc_dev = orc_dev.max(dre / scale);
                orc_im = orc_im.max(dim / scale);
            }
        }
    }
    gates.push(gate(
        "G17",
        "values-only oracle: fs-la complex QR eig on M^(-1/2)KM^(-1/2) matches the real-symmetric-embedding eigenvalues to 1e-10 (scaled); worst embedding pair residual below 1e-11 (scaled)",
        orc_dev,
        0.0,
        orc_dev <= 1e-10 && orc_im <= 1e-10 && glob.pair_resid <= 1e-11 * 830.0,
        format!(
            "oracle re dev = {orc_dev:.3e}, im = {orc_im:.3e}, pair residual = {:.3e}",
            glob.pair_resid
        ),
    ));

    // ================= CSV referee curves ================
    // Floating-point honesty: the CSV's own omega values carry the absolute
    // eigensolver noise of numpy's LAPACK (~eps * ||H||), which at k ~ 1e-4
    // is ~1e-7 RELATIVE on the light branch (omega^2 ~ 1e-8 vs ||H|| up to
    // ~840 at m_V = 20 Case B). A blanket 1e-9 relative match is therefore
    // not achievable by ANY reimplementation at small k; the gate uses the
    // conditioning-aware tolerance |w2_ours - w2_csv| <= 1e-9 w2_csv +
    // 64 eps tr(H_t), plus the strict literal 1e-9 relative on omega for
    // k >= 0.1 where conditioning permits it.
    let combos: [(MassCase, usize); 8] = [
        (MassCase::A, 0),
        (MassCase::A, 1),
        (MassCase::A, 2),
        (MassCase::A, 3),
        (MassCase::B, 0),
        (MassCase::B, 1),
        (MassCase::B, 2),
        (MassCase::B, 3),
    ];
    let mut csv_ratio = 0.0f64; // worst dev / tolerance
    let mut csv_strict = 0.0f64; // worst rel dev on omega, k >= 0.1
    let mut csv_relmax = 0.0f64; // worst rel dev on omega, all k
    for (c, (case, imv)) in combos.iter().enumerate() {
        let sw = match case {
            MassCase::A => &sweeps_a[*imv],
            MassCase::B => &sweeps_b[*imv],
        };
        let m = with_case(*case, MV_LIST[*imv]);
        for i in 0..kk.len() {
            let w2_ref = csv.omega[c][i] * csv.omega[c][i];
            let w2_ours = sw.w2mean[i].max(0.0);
            // trace of the transverse block of M^(-1/2) K M^(-1/2)
            let kmat = stiffness([0.0, 0.0, kk[i]], &m);
            let tr: f64 = [0usize, 1, 3, 4]
                .iter()
                .map(|&a| kmat[a * 6 + a].re / md[a])
                .sum();
            let tol = 1e-9 * w2_ref + 64.0 * EPS64 * tr;
            let dev = (w2_ours - w2_ref).abs();
            csv_ratio = csv_ratio.max(dev / tol);
            let om_rel = (w2_ours.sqrt() - csv.omega[c][i]).abs() / csv.omega[c][i];
            csv_relmax = csv_relmax.max(om_rel);
            if kk[i] >= 0.1 {
                csv_strict = csv_strict.max(om_rel);
            }
        }
    }
    gates.push(gate(
        "G18",
        "CSV referee: omega_B2(k) reproduces theorem_II2_invariance.csv (400 k x 2 cases x 4 m_V) within the conditioning-aware tolerance 1e-9 w2 + 64 eps tr(H); literal 1e-9 relative holds for k >= 0.1",
        csv_strict,
        0.0,
        csv_ratio <= 1.0 && csv_strict <= 1e-9,
        format!(
            "worst dev/tol = {csv_ratio:.3}, max rel(omega) k>=0.1 = {csv_strict:.3e}, all k = {csv_relmax:.3e}"
        ),
    ));

    // ================= decoupling + global sanity ================
    gates.push(gate(
        "G14",
        "exact-decoupling cross-check: sub-block union (u_z, phi_z, transverse) matches the full 6x6 spectrum, max dev(omega^2) below 1e-11",
        glob.xval,
        0.0,
        glob.xval <= 1e-11,
        format!("max dev = {:.3e}", glob.xval),
    ));
    gates.push(gate(
        "G15",
        "global sanity: most negative omega^2 anywhere >= -1e-12",
        glob.min_w2,
        0.0,
        glob.min_w2 >= -1e-12,
        format!("min omega^2 = {:.3e}", glob.min_w2),
    ));

    // ================= time domain (NEW: the repo's first) ================
    let kmat1 = stiffness([0.0, 0.0, 1.0], &bench);
    let full1 = solve(&kmat1, &md, None);
    let mut i_b1 = usize::MAX;
    let mut i_b4 = usize::MAX;
    let mut trans: Vec<usize> = Vec::new();
    for n in 0..6 {
        let v = &full1.vecs[n];
        let tot: f64 = v.iter().map(|z| z.norm_sq()).sum();
        if v[2].norm_sq() / tot > 0.99 {
            i_b1 = n;
        } else if v[5].norm_sq() / tot > 0.99 {
            i_b4 = n;
        } else {
            trans.push(n);
        }
    }
    let (i_b2, i_b3) = (trans[0], *trans.last().expect("transverse modes"));
    let mut td_freq_h = 0.0f64;
    let mut td_freq_c = 0.0f64;
    let mut td_edev = 0.0f64;
    let mut td_drift = 0.0f64;
    let mut td_note = String::new();
    for (label, idx) in
        [("B2", i_b2), ("B1", i_b1), ("B4", i_b4), ("B3", i_b3)]
    {
        let omega = full1.w2[idx].max(0.0).sqrt();
        let run = evolve_mode(&kmat1, &md, &full1.vecs[idx], full1.w2[idx], TD_STEPS, 1e-3 / omega);
        let dh = ((run.omega_meas - run.omega_h) / run.omega_h).abs();
        let dc = ((run.omega_meas - run.omega) / run.omega).abs();
        td_freq_h = td_freq_h.max(dh);
        td_freq_c = td_freq_c.max(dc);
        td_edev = td_edev.max(run.e_dev_max);
        td_drift = td_drift.max(run.e_drift);
        td_note.push_str(&format!(
            "{label}(k=1): omega = {:.12} meas = {:.12} (vs omega_h {dh:.1e}, vs omega {dc:.1e}, Edev {:.1e}); ",
            run.omega, run.omega_meas, run.e_dev_max
        ));
        report.push(format!(
            "time-domain {label} k=1: omega_eig = {:.12}, omega_h = {:.12}, omega_meas = {:.12}, energy dev max = {:.3e}, envelope growth = {:.3e} ({} steps, h = {:.3e})",
            run.omega, run.omega_h, run.omega_meas, run.e_dev_max, run.e_drift, run.steps, run.h
        ));
    }
    // B3 ringing at k -> 0: the Klein-Gordon gap frequency sqrt(21).
    let kmat_s = stiffness([0.0, 0.0, 1e-3], &bench);
    let full_s = solve(&kmat_s, &md, None);
    let run_s = evolve_mode(
        &kmat_s,
        &md,
        &full_s.vecs[5],
        full_s.w2[5],
        TD_STEPS,
        1e-3 / full_s.w2[5].sqrt(),
    );
    let want = w0sq.sqrt();
    let d_gap = ((run_s.omega_meas - want) / want).abs();
    td_edev = td_edev.max(run_s.e_dev_max);
    td_drift = td_drift.max(run_s.e_drift);
    report.push(format!(
        "time-domain B3 k=1e-3: omega_meas = {:.12} vs sqrt(21) = {want:.12} (rel {d_gap:.3e}), energy dev max = {:.3e}, envelope growth = {:.3e}",
        run_s.omega_meas, run_s.e_dev_max, run_s.e_drift
    ));
    gates.push(gate(
        "G19",
        "time domain: each branch's plane wave (B1, B2, B3, B4 at k = 1) rings at its eigensolve frequency — measured omega matches the discrete omega_h to 1e-9 and the continuum omega to 1e-6; B3 at k = 1e-3 rings at sqrt((m_V^2 + 4mu_c)/J) = sqrt(21) to 1e-6",
        td_freq_c,
        0.0,
        td_freq_h <= 1e-9 && td_freq_c <= 1e-6 && d_gap <= 1e-6,
        format!("{td_note}B3(k=1e-3) vs sqrt(21): {d_gap:.3e}"),
    ));
    gates.push(gate(
        "G20",
        "time domain: symplectic energy conservation over 2e5 Verlet steps per mode — max relative energy deviation below 1e-5 (bounded oscillation), first-to-second-half envelope growth (secular-drift witness) below 1e-9",
        td_edev,
        0.0,
        td_edev <= 1e-5 && td_drift <= 1e-9,
        format!("max energy dev = {td_edev:.3e}, max envelope growth = {td_drift:.3e}"),
    ));

    // ================= PHASE E2: W_chi chiral-coupling module ================
    // Spec: gap_survey_v3.json "Skyrme/chiral" (W_chi missing-pieces item) +
    // corpus II.C eq (2.2) and II.H. chi1/chi2/chi3 are FREE parameters
    // (default 0 — the chi = 0 default is what keeps every gate above and
    // the knot-statics functional unchanged); the trace-vs-deviatoric chi2
    // convention is unpinned by the text and carried as a switch.

    let chi_only = |chi1: f64, chi2: f64, chi3: f64, dev: bool| -> Moduli {
        let mut m = Moduli::bench();
        m.lam = 0.0;
        m.mu = 0.0;
        m.mu_c = 0.0;
        m.alpha = 0.0;
        m.beta = 0.0;
        m.gamma = 0.0;
        m.m_v = 0.0;
        m.chi1 = chi1;
        m.chi2 = chi2;
        m.chi3 = chi3;
        m.chi_deviatoric = dev;
        m
    };
    let smooth_field = |seed0: u64, n: usize| -> MicroField {
        let mut f = MicroField::zeros(n, 2.0 * std::f64::consts::PI);
        let mut s = seed0;
        for comp in 0..6 {
            for _ in 0..4 {
                let m = [
                    (next_unit(&mut s) * 7.0) as i64 - 3,
                    (next_unit(&mut s) * 7.0) as i64 - 3,
                    (next_unit(&mut s) * 7.0) as i64 - 3,
                ];
                let amp = 2.0 * next_unit(&mut s) - 1.0;
                let ph = 2.0 * std::f64::consts::PI * next_unit(&mut s);
                f.add_mode(comp, m, amp, ph);
            }
        }
        f
    };

    // ---- G21: FD-vs-analytic gradient of the real-space W_chi energy ----
    let mut fd_worst = 0.0f64;
    let mut fd_note = String::new();
    let chi_sets = [
        (0.4, 0.0, 0.0, false),
        (0.0, -0.3, 0.0, false),
        (0.0, 0.0, 0.25, false),
        (0.4, -0.3, 0.25, false),
        (0.4, -0.3, 0.25, true),
    ];
    for (ci, &(c1, c2, c3, dv)) in chi_sets.iter().enumerate() {
        let m = chi_only(c1, c2, c3, dv);
        let f = smooth_field(41_000 + ci as u64, 8);
        let d = smooth_field(52_000 + ci as u64, 8);
        let (gu, gphi) = wchi_grad(&f, &m);
        let dot = grad_dot(&gu, &gphi, &d);
        let t = 1e-3; // W_chi is exactly quadratic: central FD exact to round-off
        let fd = (energy_displaced(&f, &d, t, &m) - energy_displaced(&f, &d, -t, &m)) / (2.0 * t);
        let rel = (fd - dot).abs() / dot.abs().max(fd.abs()).max(1e-30);
        fd_worst = fd_worst.max(rel);
        fd_note.push_str(&format!("chi=({c1},{c2},{c3},dev={dv}): rel {rel:.2e}; "));
    }
    gates.push(gate(
        "G21",
        "W_chi real-space analytic gradient vs central FD on random smooth periodic (u, phi) fields (8^3 grid, chi1/chi2/chi3 singly and mixed, both trace conventions): worst relative deviation below 1e-6",
        fd_worst,
        0.0,
        fd_worst <= 1e-6,
        fd_note,
    ));

    // ---- G22: plane-wave consistency real-space W_chi vs symbol K_chi ----
    // The discrete plane wave is an exact eigenfunction of the central
    // difference with symbol i*k_eff, k_eff = sin(k h)/h, so
    // E_grid = (V/4) v^dag K_chi(k_eff) v is an EXACT identity (1e-12-class).
    let mut pw_worst = 0.0f64;
    let mut seed_pw: u64 = 77_100_100;
    let pw_sets = [
        (0.7, -0.4, 0.3, false),
        (0.7, -0.4, 0.3, true),
        (0.9, 0.0, 0.0, false),
        (0.0, 0.8, 0.0, false),
        (0.0, 0.8, 0.0, true),
        (0.0, 0.0, 0.6, false),
    ];
    let mut field_pw = MicroField::zeros(8, 2.0 * std::f64::consts::PI);
    for &(c1, c2, c3, dv) in &pw_sets {
        let m = chi_only(c1, c2, c3, dv);
        for _ in 0..4 {
            // integer mode in [-3,3]^3, not all zero (2m never aliases to 0 mod 8)
            let mut mode = [0i64; 3];
            loop {
                for md_ in &mut mode {
                    *md_ = (next_unit(&mut seed_pw) * 7.0) as i64 - 3;
                }
                if mode != [0, 0, 0] {
                    break;
                }
            }
            let mut v = [C64::ZERO; 6];
            for vc in &mut v {
                *vc = C64::new(
                    2.0 * next_unit(&mut seed_pw) - 1.0,
                    2.0 * next_unit(&mut seed_pw) - 1.0,
                );
            }
            field_pw.set_plane_wave(mode, &v);
            let (_, e_real) = wchi_energy(&field_pw, &m);
            let kmat = stiffness(field_pw.k_eff(mode), &m);
            let mut quad = 0.0f64;
            let mut scale = 0.0f64;
            for a in 0..6 {
                for b in 0..6 {
                    quad += (v[a].conj() * kmat[a * 6 + b] * v[b]).re;
                    scale += v[a].abs() * kmat[a * 6 + b].abs() * v[b].abs();
                }
            }
            let vol = field_pw.lbox().powi(3);
            let dev = (e_real - 0.25 * vol * quad).abs() / (0.25 * vol * scale).max(1e-30);
            pw_worst = pw_worst.max(dev);
        }
    }
    gates.push(gate(
        "G22",
        "dispersion-path consistency: real-space W_chi energy of a stored plane wave equals the symbol quadratic form (V/4) v^dag K_chi(k_eff) v (same couplings, k_eff = sin(kh)/h exact for central differences) — 24 random (k, polarization) samples across chi1/chi2/chi3 and both conventions, scaled deviation below 1e-12",
        pw_worst,
        0.0,
        pw_worst <= 1e-12,
        format!("worst scaled deviation = {pw_worst:.3e} over 24 samples"),
    ));

    // ---- G23: convention switch identity, chi = 0 default, rotation invariance ----
    // (a) deviatoric(chi1, chi2) == trace-included(chi1 - chi2/3, chi2) exactly
    let mut conv_dev = 0.0f64;
    let mut seed_cv: u64 = 33_444_555;
    for kmag in [0.3, 1.7] {
        for _ in 0..3 {
            let d = rand_dir(&mut seed_cv);
            let kv = [kmag * d[0], kmag * d[1], kmag * d[2]];
            let kd = stiffness(kv, &chi_only(0.2, 0.9, 0.15, true));
            let kt = stiffness(kv, &chi_only(0.2 - 0.9 / 3.0, 0.9, 0.15, false));
            let scale = kt.iter().map(|z| z.abs()).fold(1.0, f64::max);
            for i in 0..36 {
                conv_dev = conv_dev.max((kd[i] - kt[i]).abs() / scale);
            }
        }
    }
    // ... and the same identity for the real-space DENSITY, pointwise
    // (the total of a random field is a near-cancelling sum — the density
    // comparison is the stronger statement anyway)
    let f_cv = smooth_field(61_000, 8);
    let (dens_dev, _) = wchi_energy(&f_cv, &chi_only(0.2, 0.9, 0.15, true));
    let (dens_tr, _) = wchi_energy(&f_cv, &chi_only(0.2 - 0.9 / 3.0, 0.9, 0.15, false));
    let dens_scale = dens_tr.iter().map(|x| x.abs()).fold(1e-30, f64::max);
    let conv_rs = dens_dev
        .iter()
        .zip(dens_tr.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0f64, f64::max)
        / dens_scale;
    // (b) chi = 0 default: density, total, and gradient EXACTLY zero
    let m0 = Moduli::bench();
    let (dens0, tot0) = wchi_energy(&f_cv, &m0);
    let (gu0, gphi0) = wchi_grad(&f_cv, &m0);
    let zero_ok = tot0 == 0.0
        && dens0.iter().all(|&x| x == 0.0)
        && gu0.iter().chain(gphi0.iter()).all(|v| v.iter().all(|&x| x == 0.0));
    // (c) rotation invariance of the full spectrum with chi1/chi2/chi3 on
    let mut m_full = Moduli::bench();
    m_full.chi1 = 0.2;
    m_full.chi2 = 0.15;
    m_full.chi3 = 0.3;
    let mut rot_chi = 0.0f64;
    let mut seed_rc: u64 = 91_929_394;
    for kmag in [0.1, 1.0, 3.0] {
        let kz = stiffness([0.0, 0.0, kmag], &m_full);
        let wz = solve(&kz, &md, None);
        let scale = wz.w2[5].abs().max(1.0);
        for _ in 0..4 {
            let d = rand_dir(&mut seed_rc);
            let kg = stiffness([kmag * d[0], kmag * d[1], kmag * d[2]], &m_full);
            let wg = solve(&kg, &md, None);
            for i in 0..6 {
                rot_chi = rot_chi.max((wg.w2[i] - wz.w2[i]).abs() / scale);
            }
        }
    }
    gates.push(gate(
        "G23",
        "W_chi conventions: (a) deviatoric(chi1, chi2) equals trace-included(chi1 - chi2/3, chi2) exactly (symbol below 1e-13 scaled, real-space density pointwise below 1e-12 scaled — the switch is a reparametrization, which is WHY the text can leave it unpinned); (b) chi = 0 default gives exactly zero density/energy/gradient (knot statics untouched); (c) full spectrum with chi1/chi2/chi3 on is rotation-invariant to 1e-12 (scaled)",
        conv_dev.max(rot_chi),
        0.0,
        conv_dev <= 1e-13 && conv_rs <= 1e-12 && zero_ok && rot_chi <= 1e-12,
        format!(
            "symbol conv dev = {conv_dev:.3e}, real-space conv dev = {conv_rs:.3e}, chi=0 exact zeros = {zero_ok}, rotation dev = {rot_chi:.3e}"
        ),
    ));

    // ---- G24: Dzyaloshinskii tilt at the r9 margin m = 1.9 ----
    let sin_h3 = |mm: f64| (1.0 - 1.0 / mm).sqrt();
    let s19 = SoftSector::with_margin(1.9, 0.7, 1.1);
    let th19 = s19.theta_c_numeric();
    let sin_num = th19.sin();
    let sin_want = sin_h3(1.9); // 0.6882472016116852 -> corpus prints 0.69 +/- 0.11
    let tilt_dev = (sin_num - sin_want).abs();
    let stat_resid = s19.df_tilt(th19).abs() / (s19.delta_s * s19.delta_s);
    let q_dev = ((s19.pitch_numeric(0.8) - s19.pitch_star()) / s19.pitch_star()).abs();
    let in_band = (sin_num - 0.69).abs() <= 0.11;
    gates.push(gate(
        "G24",
        "Dzyaloshinskii soft sector at the ONLY text-pinned margin m = 1.9 (r9, Thm H2-1): numerical minimization of f(theta) gives sin(theta_c) = sqrt(1 - 1/m) = 0.688247201612 (Thm H-3; Tier-5a C7a's 0.688 -> 0.69) within 1e-12, stationarity residual below 1e-12, numerical pitch equals q* = chi_s/gamma_s within 1e-9, and the tilt sits inside the printed 0.69 +/- 0.11",
        sin_num,
        sin_want,
        tilt_dev <= 1e-12 && stat_resid <= 1e-12 && q_dev <= 1e-9 && in_band,
        format!(
            "sin theta_c = {sin_num:.15} (analytic {sin_want:.15}, dev {tilt_dev:.2e}), f'(theta_c)/Delta^2 = {stat_resid:.2e}, pitch rel dev = {q_dev:.2e}, in 0.69+/-0.11 band = {in_band}"
        ),
    ));

    // ---- G25: margin scan m in [1.1, 3] against the Thm H-3 tilt ----
    let mut scan_worst = 0.0f64;
    let gd_combos = [(0.7, 1.1), (1.3, 0.9), (2.0, 0.5)];
    for i in 0..39usize {
        let mm = 0.05f64.mul_add(i as f64, 1.1);
        let (gs, ds) = gd_combos[i % 3];
        let s = SoftSector::with_margin(mm, gs, ds);
        let dev = (s.theta_c_numeric().sin() - sin_h3(mm)).abs();
        scan_worst = scan_worst.max(dev);
    }
    gates.push(gate(
        "G25",
        "margin scan m in [1.1, 3] (39 points, three (gamma_s, Delta_s) parametrizations — only the DIMENSIONLESS margin matters, Thm H2-1): numerical tilt matches sin(theta_c) = sqrt(1 - 1/m) within 1e-10 everywhere",
        scan_worst,
        0.0,
        scan_worst <= 1e-10,
        format!("worst |sin theta_num - sqrt(1 - 1/m)| = {scan_worst:.3e}"),
    ));

    // ---- G26: the Dzyaloshinskii threshold at m = 1 ----
    let mut below_max_tilt = 0.0f64;
    let mut below_ok = true;
    for mm in [0.3, 0.6, 0.9, 0.999] {
        let s = SoftSector::with_margin(mm, 1.3, 0.9);
        let th = s.theta_c_numeric();
        below_max_tilt = below_max_tilt.max(th);
        below_ok &= !s.condenses() && th <= 1e-6 && s.d2f_tilt(0.0) > 0.0;
    }
    let mut above_ok = true;
    let mut above_note = String::new();
    for mm in [1.001, 1.01, 1.1] {
        let s = SoftSector::with_margin(mm, 1.3, 0.9);
        let dev = (s.theta_c_numeric().sin() - sin_h3(mm)).abs();
        above_ok &= s.condenses() && dev <= 1e-8 && s.d2f_tilt(0.0) < 0.0;
        above_note.push_str(&format!("m={mm}: dev {dev:.1e}; "));
    }
    // bisect the f''(0) sign change (numeric FD, no analytic shortcut)
    let fpp0 = |mm: f64| -> f64 {
        let s = SoftSector::with_margin(mm, 1.3, 0.9);
        let d = 1e-4;
        2.0 * (s.f_tilt(d) - s.f_tilt(0.0)) / (d * d)
    };
    let (mut lo, mut hi) = (0.5f64, 2.0f64);
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if fpp0(mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let m_star = 0.5 * (lo + hi);
    // the literal printed quadratic gap shares the same threshold
    let quad_below = SoftSector::with_margin(0.9, 1.0, 1.0);
    let quad_above = SoftSector::with_margin(1.9, 1.0, 1.0);
    let quad_ok = (0..14).all(|i| quad_below.f_tilt_quadgap(0.1 * (i + 1) as f64) > 0.0)
        && quad_above.f_tilt_quadgap(0.5) < 0.0;
    gates.push(gate(
        "G26",
        "Dzyaloshinskii threshold structure: NO tilt for m <= 1 (theta_c below 1e-6 at m in (0.3, 0.6, 0.9, 0.999), f''(0) > 0), tilt matching sqrt(1 - 1/m) within 1e-8 just above (m in (1.001, 1.01, 1.1), f''(0) < 0); numeric-FD bisection of the f''(0) sign change lands at m* = 1 within 1e-6; the literal quadratic-gap reading shares the same threshold (co-gated)",
        m_star,
        1.0,
        below_ok && above_ok && (m_star - 1.0).abs() <= 1e-6 && quad_ok,
        format!(
            "max tilt below threshold = {below_max_tilt:.2e}; {above_note}bisected m* = {m_star:.9}; quadgap threshold ok = {quad_ok}"
        ),
    ));

    // sort gates by id for a stable table and stable Merkle input order
    gates.sort_by(|a, b| a.id.cmp(b.id));

    // ============== EVIDENCE + PACKAGE + INDEPENDENT CHECK ==============
    let mut pkg = EvidencePackage::new(Provenance::new(
        "fs-gum-cosserat-v0.1.0-phaseB3-tier1-referee-battery",
        "frankensim-fs-la-fs-math-fs-ivl-fs-evidence-fs-package-fs-checker",
    ));
    let mut certified = 0usize;
    for g in &gates {
        if !g.pass || !g.measured.is_finite() {
            continue;
        }
        let prov = ProvenanceHash::of_bytes(
            format!("{}|{}|{:016x}", g.id, g.statement, g.measured.to_bits()).as_bytes(),
        );
        let ev = Evidence::enclosed(g.measured, g.measured, g.measured, prov);
        match ev.certified() {
            Ok(_cert) => {
                certified += 1;
                pkg = pkg.with_claim(Claim::from_certificate(
                    g.id,
                    &g.statement,
                    g.measured,
                    g.measured,
                    "fs-gum-cosserat",
                    certificate_hash(g.id, &g.statement, g.measured, g.measured).to_hex(),
                ));
            }
            Err(e) => println!("  [fs-evidence refused {}: {e:?}]", g.id),
        }
    }
    let root = pkg.try_merkle_root().expect("package root");
    let root_hex = root.to_hex();

    let mode_deny = fs_checker::check_against_root(&pkg, root).passed();
    let verifier = GatesCertVerifier;
    let caps = VerificationCapabilities::deny_all().with_source_certificates(&verifier);
    let mode_cert = fs_checker::check_with_capabilities(&pkg, Some(root), None, &caps).passed();
    let mut wrong = root;
    wrong.0[0] ^= 0xff;
    let mode_tamper =
        fs_checker::check_with_capabilities(&pkg, Some(wrong), None, &caps).passed();

    RunOutput { gates, certified, root_hex, mode_deny, mode_cert, mode_tamper, report }
}

fn main() {
    let t0 = std::time::Instant::now();
    println!("== fs-gum-cosserat GATES: GUM Phase B3 micropolar/Cosserat referee battery ==");
    println!(
        "port source: tier1-spectrum/spectrum.py; referees: REPORT.md + theorem_II2_invariance.csv"
    );
    println!(
        "benchmark: rho0=J=lambda=mu=1, mu_c=5, alpha=beta=gamma=0.5, m_V=1; 400-point k grid [1e-4, 3]\n"
    );

    let csv = read_csv();
    let first = run(&csv);

    let mut passes = 0usize;
    for g in &first.gates {
        if g.pass {
            passes += 1;
        }
        println!(
            "{:<4} {}  measured {:.12e}  target {:.6e}",
            g.id,
            if g.pass { "PASS" } else { "FAIL" },
            g.measured,
            g.target
        );
        println!("      {}", g.statement);
        if !g.note.is_empty() {
            println!("      note: {}", g.note);
        }
    }
    println!("\n== {passes}/{} gates PASS ==\n", first.gates.len());

    for line in &first.report {
        println!("{line}");
    }

    println!(
        "\nfs-evidence: {}/{} passing gates certified (fail-closed; failing gates are excluded by construction).",
        first.certified, first.gates.len()
    );
    println!("fs-package Merkle root: {}", first.root_hex);
    println!(
        "fs-checker mode 1 (deny-all):        passed = {} (expected false — anti-laundering)",
        first.mode_deny
    );
    println!("fs-checker mode 2 (cert capability): passed = {} (expected true)", first.mode_cert);
    println!("fs-checker mode 3 (tampered root):   passed = {} (expected false)", first.mode_tamper);

    // Replay determinism: the ENTIRE battery (eigensolves, sweeps, Verlet
    // trajectories, certificates, Merkle assembly) runs a second time and
    // must reproduce the root bit-for-bit.
    let second = run(&csv);
    let replay_ok = second.root_hex == first.root_hex;
    println!(
        "replay determinism: second full run root {} first ({})",
        if replay_ok { "==" } else { "!=" },
        second.root_hex
    );

    let all_pass = passes == first.gates.len()
        && first.certified == first.gates.len()
        && !first.mode_deny
        && first.mode_cert
        && !first.mode_tamper
        && replay_ok;
    assert!(replay_ok, "replay determinism violated");
    println!(
        "\nOVERALL: {}   (runtime {:.1} s)",
        if all_pass { "PASS (all gates, all checker modes, replay)" } else { "FAIL" },
        t0.elapsed().as_secs_f64()
    );

    println!(
        "\nEPISTEMIC NOTICE: every PASS certifies a WITHIN-MODEL property of a\n\
         speculative theory's linearized dispersion. It validates the port, the\n\
         eigensolve path, and the integrator — never the physics."
    );

    if !all_pass {
        std::process::exit(1);
    }
}
