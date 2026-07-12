//! GUM core Phase F1 — `fs-gum-e2e`: the end-to-end certified demo.
//!
//! ONE binary that runs the whole GUM knot pipeline across the three
//! finished gap-closure crates and emits ONE EvidencePackage golden root —
//! the integration proof that the gap-closure crates compose:
//!
//!   1. SEED     `fs-gum-field::radial` solves the 1-D radial profile at the
//!               frozen eps = 0.05 dial t = 0.0082764349; the 3-D hedgehog
//!               is seeded from it on the N = 48, LBOX = 4.5 stored field.
//!   2. MEASURE  `fs-gum-field::measure` (central4): sector energies,
//!               degree, weighted Derrick virial, boundary tail — claims.
//!   3. RELAX    `fs-gum-statics`: fixed-seed bump perturbation (LCG
//!               424242), guarded arrested Newton flow (one-sided degree
//!               anchor MU = 5000 / band 0.005, one-sided Bogomolny-floor
//!               wall MU_F = 400 / band 0.01) to stall/cap; claims: monotone
//!               objective, degree in band, floor held, hedgehog-
//!               neighbourhood metrics.
//!   4. CLOCK + ISOROTATE  a fixed-L Routhian descent (the Phase-E1 control
//!               protocol, L = 0.12 I_hedgehog, below threshold), then the
//!               Section-K clock bisection L_clock^2 = (2/3) I E_static on
//!               its endpoint; claims: E_rot/E_tot = 1/4 at 1e-9, bisection
//!               vs closed form, and the F-R5-relevant number
//!               kappa(L_clock)/threshold in the stated band [1.1, 1.4]
//!               (expect ~1.25; Phase E1 measured 1.2525 at N = 48, the
//!               Python N = 96 reference 1.245).
//!   5. TOPO     `fs-gum-topo` on the RELAXED field through a thin
//!               stored-field conversion: independent central4 degree
//!               (cross-crate agreement claim), geometric VOS solid-angle
//!               charge (integer snap), director projection onto S2; plus
//!               the Hopf exhibit — Whitehead/FFT vs preimage-linking cross-
//!               validation on fs-gum-topo's ANALYTIC HOPF-1 REFEREE TEXTURE
//!               (the relaxed hedgehog's director is degree-driven, not
//!               Hopf; the referee texture, NOT a knot solution, is the
//!               demo's topological exhibit, stated as such in the claims).
//!   6. PACKAGE  every claim fail-closed Certified<f64> -> one
//!               EvidencePackage (domain `gum-core:e2e:v1`) -> golden BLAKE3
//!               Merkle root; fs-checker three modes (deny-all expect-fail,
//!               source-certificate capability expect-pass, tampered root
//!               expect-fail); the ENTIRE pipeline runs twice in-process and
//!               the root must be bit-identical.
//!
//! Determinism contract (inherited from all three siblings): plain
//! sequential f64 loops in fixed ascending order, deterministic LCG
//! perturbations, no parallelism, no platform libm in kernels.
//!
//! Epistemic notice (binding): every PASS certifies a WITHIN-MODEL property
//! of discretisations, solvers, diagnostics and an evidence pipeline
//! composed across crates — numerical engineering on a speculative theory's
//! functional. Nothing here validates the theory, and nothing here says
//! anything about nature.

#![forbid(unsafe_code)]

use fs_blake3::ContentHash;
use fs_evidence::{Evidence, ProvenanceHash};
use fs_gum_field::radial::{degree as radial_degree, monotone_decreasing, radial_solve};
use fs_gum_field::{bps_floor, measure, virial, Field3, Scheme, Sectors, T_FROZEN};
use fs_gum_statics::diag::{
    add_scaled, bump_field, clock_bisect, erot_frac, halo_fraction, halo_seed_field,
    restore_cells, save_cells,
};
use fs_gum_statics::{anf, eval, kappa_threshold, AnfParams, AnfResult, Opts, FLOOR_BAND};
use fs_gum_topo::degree::{
    cube_shell, degree_b_density, hedgehog_charge_geometric, Scheme as TopoScheme,
};
use fs_gum_topo::hopf::hopf_whitehead;
use fs_gum_topo::linking::hopf_preimage_linking;
use fs_gum_topo::referee::{hopfion_director, vector_director};
use fs_gum_topo::{hopf_project, Bc, QField};
use fs_ivl::Interval;
use fs_package::origin::{
    SourceCertificateRequest, SourceCertificateVerifier, VerificationCapabilities,
    VerificationDecision,
};
use fs_package::{Claim, EvidencePackage, Provenance};

use std::time::Instant;

// ---------------------------------------------------------------------------
// frozen pipeline parameters (Phase E1 protocol, N = 48 demo resolution)
// ---------------------------------------------------------------------------

/// Demo grid: N = 48 cells per axis (the Phase E1 relaxation grid).
const N_RUN: usize = 48;
/// eps = 0.05 hedgehog box half-width (field3d LBOX).
const LBOX: f64 = 4.5;
/// Campaign perturbation seed (Phase E1 G-B/G-C protocol).
const SEED_PERT: u64 = 424242;
/// Static relaxation iteration cap (Phase E1 G-B).
const MAXIT_STATIC: usize = 900;
/// Isorotating (Routhian control) descent iteration cap (Phase E1 G-C).
const MAXIT_ISO: usize = 600;
/// Isorotation spin target kappa_0 = L / I_hedgehog (the Phase E1 control:
/// below threshold 1/sqrt(8 pi), so the descent dilates without shedding).
const KAPPA_ISO: f64 = 0.12;
/// Hopf-referee grid (power-of-two for the Whitehead/FFT route).
const N_HOPF: usize = 64;
/// Hopf-referee half-box (the compacton/pilot box).
const HALF_HOPF: f64 = 2.4;
/// Geometric-charge probe: cube-shell half-width and quads per face.
const SHELL_S: f64 = 1.2;
const SHELL_M: usize = 24;

/// radial_results.json N_run references (N = 4000, t frozen) — the same
/// referee constants certified by fs-gum-field gate R1.
const RAD_REF: [f64; 4] =
    [12.117459252278454, 6.191126358808831, 1.523935580428207, 1.50758847423602];

/// Step-1 radial references (6-dp, the field3d REF dict) for the 3-D
/// hedgehog sector rel errors.
const REF6: [(&str, f64); 5] = [
    ("E2", 12.117459),
    ("E4", 6.191126),
    ("E6", 1.523936),
    ("E0", 1.507588),
    ("I", 22.26826),
];

// ---------------------------------------------------------------------------
// certificate plumbing (sibling pattern, domain gum-core:e2e:v1)
// ---------------------------------------------------------------------------

/// The certification domain of this demo's claims.
const DOMAIN: &str = "gum-core:e2e:v1";

fn certificate_hash(id: &str, statement: &str, lo: f64, hi: f64) -> ContentHash {
    let canon = format!("{id}|{statement}|{:016x}|{:016x}", lo.to_bits(), hi.to_bits());
    fs_blake3::hash_domain(DOMAIN, canon.as_bytes())
}

struct E2eCertVerifier;

impl SourceCertificateVerifier for E2eCertVerifier {
    fn verify(&self, req: &SourceCertificateRequest<'_>) -> VerificationDecision {
        let policy = fs_blake3::hash_domain(DOMAIN, b"policy:recompute-and-compare");
        let expect = certificate_hash(req.claim_id, req.statement, req.lo, req.hi);
        if req.producer == "fs-gum-e2e" && req.certificate_hash == expect {
            VerificationDecision::accept(policy)
        } else {
            VerificationDecision::reject(policy)
        }
    }
}

// ---------------------------------------------------------------------------
// helpers
// ---------------------------------------------------------------------------

struct Check {
    id: &'static str,
    statement: String,
    enclosure: Interval,
    target: f64,
    pass: bool,
    note: String,
}

fn check(
    id: &'static str,
    statement: &str,
    measured: f64,
    target: f64,
    pass: bool,
    note: String,
) -> Check {
    // A non-finite measurement can never be enclosed: record a point-zero
    // sentinel with pass = false so the claim is excluded fail-closed.
    let (enclosure, pass) = if measured.is_finite() {
        (Interval::point(measured), pass)
    } else {
        (Interval::point(0.0), false)
    };
    Check { id, statement: statement.to_string(), enclosure, target, pass, note }
}

fn rel(a: f64, b: f64) -> f64 {
    a / b - 1.0
}

fn sector_of(s: &Sectors, key: &str) -> f64 {
    match key {
        "E2" => s.e2,
        "E4" => s.e4,
        "E6" => s.e6,
        "E0" => s.e0,
        "I" => s.i,
        _ => unreachable!("unknown sector key"),
    }
}

/// Thin conversion from fs-gum-field's padded-SoA stored field to
/// fs-gum-topo's minimal node container: same cell-centred grid
/// (node i at -half + (i + 0.5) h), same fixed (i, j, k) order, same
/// vacuum-ghost boundary semantics (fs-gum-field's 2-layer rind and
/// fs-gum-topo's `Bc::Vacuum` ghost resolve to the identical (1,0,0,0)
/// value at every out-of-range read of the width-2 central4 stencil).
fn to_topo_qfield(f: &Field3) -> QField {
    let n = f.n();
    let mut data = Vec::with_capacity(n * n * n);
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                data.push(f.get(i, j, k));
            }
        }
    }
    let h = f.h();
    let o = -f.half() + 0.5 * h;
    QField::from_parts([n; 3], h, [o; 3], Bc::Vacuum([1.0, 0.0, 0.0, 0.0]), data)
}

struct RunOutput {
    checks: Vec<Check>,
    certified: usize,
    root_hex: String,
    fingerprint_hex: String,
    mode_deny: bool,
    mode_cert: bool,
    mode_tamper: bool,
}

// ---------------------------------------------------------------------------
// one full deterministic pass of the pipeline
// ---------------------------------------------------------------------------

#[allow(clippy::too_many_lines)]
fn run(verbose: bool) -> RunOutput {
    let t0 = Instant::now();
    let mut fp: Vec<f64> = Vec::new(); // determinism fingerprint values
    let mut checks: Vec<Check> = Vec::new();
    let say = |msg: &str| {
        if verbose {
            println!("{msg}");
        }
    };

    // ================= 1. SEED: 1-D radial profile -> 3-D hedgehog ==========
    let rp = radial_solve(T_FROZEN, 4000, 6.0);
    let [re2, re4, re6, re0] = rp.sectors;
    let ratio = T_FROZEN * (re2 + re4) / (re6 + re0);
    let k1d = radial_degree(&rp.r, &rp.f);
    let mono = monotone_decreasing(&rp.f);
    let worst_rad = rp
        .sectors
        .iter()
        .zip(RAD_REF.iter())
        .map(|(&a, &b)| rel(a, b).abs())
        .fold(0.0_f64, f64::max);
    fp.extend_from_slice(&[re2, re4, re6, re0, ratio, k1d, rp.gmax]);
    say(&format!(
        "[1 SEED] radial solve N=4000: E2={re2:.6} E4={re4:.6} E6={re6:.6} E0={re0:.6} gmax={:.1e} ({:.1}s)",
        rp.gmax,
        t0.elapsed().as_secs_f64()
    ));
    checks.push(check(
        "S1",
        "seed radial solve (N=4000, t frozen) reproduces the radial_results.json sector references to < 1e-6 rel with Newton floor gmax < 1e-11",
        worst_rad,
        0.0,
        worst_rad < 1.0e-6 && rp.gmax < 1.0e-11,
        format!("worst sector rel {worst_rad:.2e}; gmax = {:.1e}", rp.gmax),
    ));
    checks.push(check(
        "S2",
        "eps dial: ratio t (E2+E4)/(E6+E0) at the frozen t equals 0.05 within 5e-5",
        ratio,
        0.05,
        (ratio - 0.05).abs() < 5.0e-5,
        format!("ratio = {ratio:.9}"),
    ));
    checks.push(check(
        "S3",
        "1-D degree K of the seed profile equals 1 to < 1e-6",
        k1d,
        1.0,
        (k1d - 1.0).abs() < 1.0e-6,
        format!("K = {k1d:.12}"),
    ));
    checks.push(check(
        "S4",
        "seed profile is monotone non-increasing pi -> 0 (no overshoot)",
        if mono { 1.0 } else { 0.0 },
        1.0,
        mono,
        String::new(),
    ));
    let mut f = Field3::new(N_RUN, LBOX);
    f.sample_hedgehog(&rp.r, &rp.f, 1.0);

    // ================= 2. MEASURE: central4 sectors on the seed =============
    let sec = measure(&f, Scheme::Central4, T_FROZEN);
    let rels: Vec<f64> = REF6.iter().map(|&(k, v)| rel(sector_of(&sec, k), v)).collect();
    let worst_sec = rels.iter().fold(0.0_f64, |m, &v| m.max(v.abs()));
    let vir_rel = virial(&sec, T_FROZEN) / sec.estat;
    let tail = f.boundary_tail();
    fp.extend_from_slice(&[sec.e2, sec.e4, sec.e6, sec.e0, sec.i, sec.deg, vir_rel, tail]);
    say(&format!(
        "[2 MEASURE] central4 N={N_RUN}: Estat={:.7} deg={:.6} I={:.5} virial rel={vir_rel:+.2e} ({:.1}s)",
        sec.estat,
        sec.deg,
        sec.i,
        t0.elapsed().as_secs_f64()
    ));
    checks.push(check(
        "M1",
        "seeded hedgehog central4 N=48 sectors vs the Step-1 radial references: worst rel error < 5e-2 (the N=48 discretisation class; the same instrument is gated < 1e-2 at N=96 in fs-gum-field H1)",
        worst_sec,
        0.0,
        worst_sec < 5.0e-2,
        REF6.iter()
            .zip(rels.iter())
            .map(|(&(k, _), &e)| format!("{k} {e:+.2e}"))
            .collect::<Vec<_>>()
            .join("  "),
    ));
    checks.push(check(
        "M2",
        "seeded hedgehog central4 N=48 degree within 5e-3 of 1 (fs-gum-topo G1d recorded 9.2e-4 at this resolution)",
        sec.deg,
        1.0,
        (sec.deg - 1.0).abs() < 5.0e-3,
        format!("deg = {:.9}, err {:.2e}", sec.deg, (sec.deg - 1.0).abs()),
    ));
    checks.push(check(
        "M3",
        "weighted Derrick virial on the seeded hedgehog (central4, N=48): rel magnitude < 1e-2 (Phase E1 recorded -6.30e-3 on this grid; the N=96 referee is -6.6e-4)",
        vir_rel,
        0.0,
        vir_rel.abs() < 1.0e-2,
        format!("virial rel = {vir_rel:+.3e}"),
    ));
    checks.push(check(
        "M4",
        "boundary tail max vector amplitude over pi on the box faces < 1e-6 (decay-to-vacuum ghost policy)",
        tail,
        0.0,
        tail < 1.0e-6,
        format!("tail = {tail:.3e}"),
    ));

    // ================= 3. RELAX: guarded arrested Newton flow ================
    let (out_h, _) = eval(&f, &Opts::estatic(Scheme::Corner), false);
    let deg_ref = out_h.deg;
    let fgap_ref = out_h.e6 + out_h.e0 - bps_floor() - FLOOR_BAND;
    let opts_static = Opts::estatic(Scheme::Corner).with_guards(deg_ref, fgap_ref);
    let bump = bump_field(N_RUN, LBOX, SEED_PERT, 6, 0.02);
    add_scaled(&mut f, &bump, 1.0);
    let _ = f.renormalize();
    let mut prm = AnfParams::new(MAXIT_STATIC);
    prm.print_every = if verbose { 300 } else { 0 };
    say(&format!(
        "[3 RELAX] corner hedgehog ref: Estat={:.7} deg_ref={:.6} wall={:+.6}; perturbed (LCG {SEED_PERT}), ANF cap {MAXIT_STATIC}",
        out_h.estat, deg_ref, fgap_ref
    ));
    let res_s: AnfResult = anf(&mut f, &opts_static, &prm, "static");
    let out_s = res_s.out;
    let (halo_s, _) = halo_fraction(&f);
    let mut q_static = vec![0.0_f64; 4 * N_RUN * N_RUN * N_RUN];
    save_cells(&f, &mut q_static);
    let mut max_rise_s = f64::NEG_INFINITY;
    for w in res_s.series.windows(2) {
        max_rise_s = max_rise_s.max(w[1].obj - w[0].obj);
    }
    fp.extend_from_slice(&[
        deg_ref, fgap_ref, out_s.estat, out_s.deg, out_s.i, out_s.floor_gap, out_s.epen,
        out_s.efpen, halo_s, max_rise_s, res_s.iters as f64, res_s.arrests as f64,
    ]);
    say(&format!(
        "[3 RELAX] done ({}, {} iters, {} arrests): Estat={:.7} deg={:.6} fgap={:+.6} halo={:.1e} ({:.1}s)",
        res_s.status.as_str(),
        res_s.iters,
        res_s.arrests,
        out_s.estat,
        out_s.deg,
        out_s.floor_gap,
        halo_s,
        t0.elapsed().as_secs_f64()
    ));
    checks.push(check(
        "R1",
        "relaxation objective is monotone non-increasing by construction (max recorded rise <= 0 across the instrumented series)",
        max_rise_s,
        0.0,
        max_rise_s <= 0.0,
        format!("max recorded rise {max_rise_s:+.2e}"),
    ));
    let deg_drift = (out_s.deg - deg_ref).abs();
    checks.push(check(
        "R2",
        "relaxed degree within the one-sided anchor geometry: drift from deg_ref <= band 5e-3 + penalty penetration 3e-3",
        out_s.deg,
        deg_ref,
        deg_drift <= 8.0e-3,
        format!("deg {:.6} vs deg_ref {deg_ref:.6}, drift {deg_drift:.2e}", out_s.deg),
    ));
    let fgap_pen = out_s.floor_gap - fgap_ref;
    checks.push(check(
        "R3",
        "relaxed energy holds the guarded Bogomolny floor: floor gap >= wall - 5e-3",
        out_s.floor_gap,
        fgap_ref,
        out_s.floor_gap >= fgap_ref - 5.0e-3,
        format!("fgap {:+.6} vs wall {fgap_ref:+.6} (penetration {fgap_pen:+.2e})", out_s.floor_gap),
    ));
    let sect_rels = [
        rel(out_s.e2, out_h.e2).abs(),
        rel(out_s.e4, out_h.e4).abs(),
        rel(out_s.e6, out_h.e6).abs(),
        rel(out_s.e0, out_h.e0).abs(),
    ];
    let worst_drift = sect_rels.iter().fold(0.0_f64, |m, &v| m.max(v));
    let i_rel = rel(out_s.i, out_h.i).abs();
    let de_stat = out_s.estat - out_h.estat;
    checks.push(check(
        "R4",
        "relaxed solution stays in the hedgehog neighbourhood (same functional, Phase E1 calibrated tolerances): quadrature-exact I rel <= 8e-2, halo <= 2e-2, energy drift magnitude <= 0.10, soft derivative sectors worst rel <= 0.35",
        worst_drift,
        0.0,
        worst_drift <= 0.35 && i_rel <= 0.08 && de_stat.abs() <= 0.10 && halo_s <= 0.02,
        format!(
            "I rel {i_rel:.2e}, halo {halo_s:.1e}, dEstat {de_stat:+.5}, worst sector rel {worst_drift:.2e}"
        ),
    ));

    // ================= 4. CLOCK + ISOROTATE ==================================
    let l_iso = KAPPA_ISO * out_h.i;
    restore_cells(&mut f, &q_static);
    add_scaled(&mut f, &bump, 1.0);
    let hs = halo_seed_field(N_RUN, LBOX, 0.05, 3.2, 0.6);
    add_scaled(&mut f, &hs, 1.0);
    let _ = f.renormalize();
    let opts_iso = Opts::routhian(Scheme::Corner, l_iso).with_guards(deg_ref, fgap_ref);
    let mut prm_i = AnfParams::new(MAXIT_ISO);
    prm_i.print_every = if verbose { 300 } else { 0 };
    say(&format!(
        "[4 ISOROTATE] Routhian descent L={l_iso:.4} (kappa_0 ~ {KAPPA_ISO}, below threshold {:.5}), cap {MAXIT_ISO}",
        kappa_threshold()
    ));
    let res_i: AnfResult = anf(&mut f, &opts_iso, &prm_i, "isorot");
    let out_i = res_i.out;
    let mut max_rise_i = f64::NEG_INFINITY;
    for w in res_i.series.windows(2) {
        max_rise_i = max_rise_i.max(w[1].obj - w[0].obj);
    }
    let iso_deg_drift = (out_i.deg - deg_ref).abs();
    let iso_fgap_pen = out_i.floor_gap - fgap_ref;
    checks.push(check(
        "K0",
        "isorotating descent: objective monotone by construction and both guards held at the endpoint (degree drift <= 8e-3, floor gap >= wall - 5e-3)",
        max_rise_i,
        0.0,
        max_rise_i <= 0.0 && iso_deg_drift <= 8.0e-3 && out_i.floor_gap >= fgap_ref - 5.0e-3,
        format!(
            "max rise {max_rise_i:+.2e}; deg drift {iso_deg_drift:.2e}; wall penetration {iso_fgap_pen:+.2e}"
        ),
    ));
    let l_clock = clock_bisect(out_i.i, out_i.estat);
    let l_closed = ((2.0 / 3.0) * out_i.i * out_i.estat).sqrt();
    let erot = erot_frac(l_clock, out_i.i, out_i.estat);
    let kap_clock = l_clock / out_i.i;
    let kratio = kap_clock / kappa_threshold();
    fp.extend_from_slice(&[
        l_iso, out_i.estat, out_i.i, out_i.deg, out_i.floor_gap, max_rise_i, l_clock, erot,
        kap_clock, kratio,
    ]);
    say(&format!(
        "[4 CLOCK] endpoint Estat={:.6} I={:.5} -> L_clock={l_clock:.6}, E_rot/E={erot:.9}, kappa(L_clock)={kap_clock:.6} = {kratio:.4}x threshold ({:.1}s)",
        out_i.estat,
        out_i.i,
        t0.elapsed().as_secs_f64()
    ));
    checks.push(check(
        "K1",
        "clock: rotational energy fraction E_rot/E_tot at the bisected L_clock equals 1/4 within 1e-9",
        erot,
        0.25,
        (erot - 0.25).abs() <= 1.0e-9,
        format!("E_rot/E_tot = {erot:.12} (delta {:.1e})", (erot - 0.25).abs()),
    ));
    checks.push(check(
        "K2",
        "clock bisection (200 steps on [0, 20]) matches the closed form sqrt((2/3) I E_static) to < 1e-12 rel",
        l_clock,
        l_closed,
        rel(l_clock, l_closed).abs() <= 1.0e-12,
        format!("L_clock = {l_clock:.9} vs closed form {l_closed:.9} (rel {:.1e})", rel(l_clock, l_closed).abs()),
    ));
    checks.push(check(
        "K3",
        "the F-R5-relevant number: kappa(L_clock)/threshold on the isorotating endpoint lies in the stated band [1.1, 1.4] (expect ~1.25; the clock-charged configuration sits in the over-spun regime)",
        kratio,
        1.25,
        (1.1..=1.4).contains(&kratio),
        format!(
            "kappa(L_clock) = {kap_clock:.6}, threshold {:.6}, ratio {kratio:.4} (Phase E1: 1.2525 at N=48; Python N=96: 1.245)",
            kappa_threshold()
        ),
    ));

    // ================= 5. TOPO: cross-crate diagnostics ======================
    restore_cells(&mut f, &q_static); // topo runs on the RELAXED static field
    let sec_relaxed = measure(&f, Scheme::Central4, T_FROZEN);
    let qt = to_topo_qfield(&f);
    let dt = degree_b_density(&qt, TopoScheme::Central4);
    let cross_rel = rel(dt.deg, sec_relaxed.deg).abs();
    fp.extend_from_slice(&[sec_relaxed.deg, dt.deg, cross_rel]);
    checks.push(check(
        "T1",
        "cross-crate degree agreement on the relaxed field: fs-gum-topo's independently ported central4 b-density degree matches fs-gum-field's measure() degree to < 1e-12 rel (identical stencil and det convention, independent code paths through the thin field conversion)",
        dt.deg,
        sec_relaxed.deg,
        cross_rel < 1.0e-12,
        format!(
            "fs-gum-topo deg = {:.12}, fs-gum-field deg = {:.12}, rel {cross_rel:.2e}",
            dt.deg, sec_relaxed.deg
        ),
    ));
    let vdir = vector_director(&qt, [0.0, 0.0, 1.0]);
    let shell = cube_shell([0.0; 3], SHELL_S, SHELL_M);
    let qgeo = hedgehog_charge_geometric(&vdir, &shell);
    fp.push(qgeo);
    checks.push(check(
        "T2",
        "geometric VOS solid-angle hedgehog charge of the relaxed field's vector director through a closed cube shell snaps to the integer 1 (< 1e-9) — the lattice-exact route, independent of any b-density",
        qgeo,
        1.0,
        qgeo.round() == 1.0 && (qgeo - 1.0).abs() < 1.0e-9,
        format!("Q_raw = {qgeo:.12} (offset {:.2e}); shell s = {SHELL_S}, {SHELL_M}x{SHELL_M} quads/face", (qgeo - 1.0).abs()),
    ));
    let ndir = hopf_project(&qt);
    let mut worst_s2 = 0.0_f64;
    for v in &ndir.data {
        let d = ((v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt() - 1.0).abs();
        if d > worst_s2 {
            worst_s2 = d;
        }
    }
    fp.push(worst_s2);
    checks.push(check(
        "T3",
        "director projection of the relaxed field lands on S2: max norm deviation < 1e-12 (Hopf map n = R(q) z-hat on the renormalised stored field)",
        worst_s2,
        0.0,
        worst_s2 < 1.0e-12,
        format!("max deviation {worst_s2:.2e}"),
    ));
    say(&format!(
        "[5 TOPO] cross-crate degree rel {cross_rel:.1e}, VOS Q = {qgeo:.9}, S2 dev {worst_s2:.1e} ({:.1}s)",
        t0.elapsed().as_secs_f64()
    ));
    // Hopf exhibit — fs-gum-topo's ANALYTIC HOPF-1 REFEREE TEXTURE (the
    // Hopf projection of the exact eps = 0 compacton, H = 1 by
    // construction). The relaxed hedgehog's own director is degree-driven,
    // not Hopf, so the referee texture is the demo's topological exhibit —
    // it is NOT a knot solution of the model.
    let href = hopfion_director(N_HOPF, HALF_HOPF, Bc::Periodic);
    let (h_white, white_note) = match hopf_whitehead(&href) {
        Ok(w) => (
            w.hopf,
            format!("H_A = {:.9}, curl rel err {:.2e}, max Im {:.1e}", w.hopf, w.curl_rel_err, w.max_imag),
        ),
        Err(e) => (f64::NAN, format!("ERROR: {e}")),
    };
    let (h_link, link_note) = match hopf_preimage_linking(&href) {
        Ok(o) => (
            o.hopf,
            format!(
                "H_B = {:.12}, loops {}+{}, segments {}+{}",
                o.hopf, o.loops[0], o.loops[1], o.segments[0], o.segments[1]
            ),
        ),
        Err(e) => (f64::NAN, format!("ERROR: {e}")),
    };
    fp.extend_from_slice(&[h_white, h_link]);
    say(&format!(
        "[5 TOPO] Hopf referee texture N={N_HOPF}: Whitehead {h_white:.9}, preimage linking {h_link:.9} ({:.1}s)",
        t0.elapsed().as_secs_f64()
    ));
    checks.push(check(
        "T4",
        "Hopf exhibit (ANALYTIC HOPF-1 REFEREE TEXTURE, not a knot solution): Whitehead/FFT Hopf invariant at N=64 within 6e-2 of 1 (fs-gum-topo G3a calibration)",
        h_white,
        1.0,
        (h_white - 1.0).abs() < 6.0e-2,
        white_note,
    ));
    checks.push(check(
        "T5",
        "Hopf exhibit (same referee texture): preimage-linking Hopf invariant at N=64 equals the exact integer 1 within 1e-6",
        h_link,
        1.0,
        (h_link - 1.0).abs() < 1.0e-6,
        link_note,
    ));
    checks.push(check(
        "T6",
        "Hopf method-vs-method cross-validation on the referee texture: Whitehead and preimage-linking agree within 6e-2 (the N=64 Whitehead discretisation error dominates)",
        (h_white - h_link).abs(),
        0.0,
        (h_white - h_link).abs() < 6.0e-2,
        format!("difference {:.3e}", (h_white - h_link).abs()),
    ));

    // ================= 6. PACKAGE: evidence + checker ========================
    let mut pkg = EvidencePackage::new(Provenance::new(
        "fs-gum-e2e-v0.1.0-gum-core-phaseF1-e2e-certified-demo",
        "frankensim-fs-gum-field-fs-gum-statics-fs-gum-topo-fs-ivl-fs-evidence-fs-package-fs-checker-fs-blake3-std-only",
    ));
    let mut certified = 0usize;
    for c in &checks {
        if !c.pass || !c.enclosure.lo().is_finite() || !c.enclosure.hi().is_finite() {
            continue;
        }
        let prov = ProvenanceHash::of_bytes(
            format!(
                "{}|{}|{:016x}|{:016x}",
                c.id,
                c.statement,
                c.enclosure.lo().to_bits(),
                c.enclosure.hi().to_bits()
            )
            .as_bytes(),
        );
        let ev =
            Evidence::enclosed(c.enclosure.midpoint(), c.enclosure.lo(), c.enclosure.hi(), prov);
        if ev.certified().is_ok() {
            certified += 1;
            pkg = pkg.with_claim(Claim::from_certificate(
                c.id,
                &c.statement,
                c.enclosure.lo(),
                c.enclosure.hi(),
                "fs-gum-e2e",
                certificate_hash(c.id, &c.statement, c.enclosure.lo(), c.enclosure.hi()).to_hex(),
            ));
        }
    }
    let root = pkg.try_merkle_root().expect("package root");
    let root_hex = root.to_hex();

    let mode_deny = fs_checker::check_against_root(&pkg, root).passed();
    let verifier = E2eCertVerifier;
    let caps = VerificationCapabilities::deny_all().with_source_certificates(&verifier);
    let mode_cert = fs_checker::check_with_capabilities(&pkg, Some(root), None, &caps).passed();
    let mut wrong = root;
    wrong.0[0] ^= 0xff;
    let mode_tamper =
        fs_checker::check_with_capabilities(&pkg, Some(wrong), None, &caps).passed();

    let mut bytes = Vec::with_capacity(fp.len() * 8);
    for v in &fp {
        bytes.extend_from_slice(&v.to_bits().to_le_bytes());
    }
    let fingerprint_hex = fs_blake3::hash_domain("gum-core:e2e:fingerprint:v1", &bytes).to_hex();

    RunOutput { checks, certified, root_hex, fingerprint_hex, mode_deny, mode_cert, mode_tamper }
}

fn main() {
    println!("== GUM CORE PHASE F1: fs-gum-e2e certified end-to-end demo ==");
    println!(
        "pipeline: seed (fs-gum-field radial) -> measure (fs-gum-field) -> relax (fs-gum-statics\n\
         guarded ANF) -> clock + isorotate (fs-gum-statics) -> topo (fs-gum-topo) -> package\n\
         N = {N_RUN}, LBOX = {LBOX}, t = {T_FROZEN:.15} (eps = 0.05), domain {DOMAIN}\n"
    );

    let t0 = Instant::now();
    let first = run(true);
    let mut passes = 0usize;
    println!();
    for c in &first.checks {
        if c.pass {
            passes += 1;
        }
        println!(
            "{:<3} {}  [{:.12e}, {:.12e}]  target {:.6e}  {}",
            c.id,
            if c.pass { "PASS" } else { "FAIL" },
            c.enclosure.lo(),
            c.enclosure.hi(),
            c.target,
            c.statement
        );
        if !c.note.is_empty() {
            println!("       note: {}", c.note);
        }
    }
    println!("\n== {passes}/{} claims PASS ==\n", first.checks.len());

    println!(
        "fs-evidence: {}/{} claims certified (fail-closed Certified<f64>; failing claims are \
         excluded from the package by construction).",
        first.certified,
        first.checks.len()
    );
    println!("fs-package golden Merkle root: {}", first.root_hex);
    println!("run fingerprint (BLAKE3 over all pipeline numbers): {}", first.fingerprint_hex);
    println!(
        "fs-checker mode 1 (deny-all):        passed = {} (expected false — unauthenticated \
         Verified claims are refused; anti-laundering works)",
        first.mode_deny
    );
    println!("fs-checker mode 2 (cert capability): passed = {} (expected true)", first.mode_cert);
    println!("fs-checker mode 3 (tampered root):   passed = {} (expected false)", first.mode_tamper);

    // In-process double run: the ENTIRE pipeline (radial solve, relaxation,
    // isorotating descent, clock, topo diagnostics, certificates, Merkle
    // assembly) runs a second time; root and fingerprint must be
    // bit-identical.
    println!("\n[replay] second full pipeline run (silent)...");
    let second = run(false);
    let replay_ok =
        second.root_hex == first.root_hex && second.fingerprint_hex == first.fingerprint_hex;
    println!(
        "replay determinism: second-run root {} first, fingerprint {} first",
        if second.root_hex == first.root_hex { "==" } else { "!=" },
        if second.fingerprint_hex == first.fingerprint_hex { "==" } else { "!=" },
    );
    println!("total runtime: {:.1} s (both runs)", t0.elapsed().as_secs_f64());

    let all_pass = passes == first.checks.len()
        && first.certified == first.checks.len()
        && !first.mode_deny
        && first.mode_cert
        && !first.mode_tamper
        && replay_ok;
    assert!(replay_ok, "replay determinism violated");
    println!(
        "\nOVERALL: {}",
        if all_pass { "PASS (all claims certified, all checker modes, bit-identical replay)" } else { "FAIL" }
    );

    println!(
        "\nEPISTEMIC NOTICE: every PASS certifies a WITHIN-MODEL property of discretisations,\n\
         solvers, diagnostics and an evidence pipeline composed across crates — numerical\n\
         engineering on a speculative theory's functional. The Hopf exhibit runs on the\n\
         analytic referee texture, not a knot solution. Nothing here validates the theory,\n\
         and nothing here says anything about nature."
    );
    if !all_pass {
        std::process::exit(1);
    }
}
