//! GUM core Phase B1 gate binary for `fs-gum-field`.
//!
//! Gates (referees frozen by earlier tiers; every number is machine-checked
//! here, certified through fs-evidence into an fs-package bundle, verified
//! by fs-checker in three modes, and replayed for bit determinism):
//!
//!   A  eps = 0 compacton, stored field, 2nd-order stencil, HALF = 2.4:
//!      deg/E0/E6/I vs fs-cosserat-pilot's certified enclosures
//!      (E0 = E6 = 16 sqrt2/15, I = 128 pi R*^3/105, R* = 2^(5/6)) and the
//!      pilot's own certified grid decimals; measured degree order O(h^2).
//!   R  1-D radial hedgehog solve at frozen t = 0.0082764349 (N = 4000,
//!      rmax = 6): sectors vs radial_results.json, eps ratio, Derrick
//!      virial, degree, Yukawa tail mu_fit vs 1/sqrt(2t), BPS deficit,
//!      monotonicity, plus the eps-dial outer loop at N = 1000.
//!   H  3-D stored hedgehog at LBOX = 4.5, central4 measurement engine,
//!      N = 48/64/96: sector values vs the Step-1 radial references
//!      (E2 = 12.117459, E4 = 6.191126, E6 = 1.523936, E0 = 1.507588,
//!      I = 22.26826), degree, N64/N96 convergence ratios, weighted Derrick
//!      virial, boundary tail, and the recorded field3d rel-error referees.
//!   C  compact corner scheme at N = 64/96: the documented O(h^2)
//!      fingerprint (worst rel ~2.2e-2, ratios 2.1-2.2, deg 0.98601).
//!   F  field mechanics: renormalisation drift + idempotence on a
//!      deterministically perturbed hedgehog; Maurer-Cartan identity
//!      |a_i|^2 = |D_i q|^2 and the scalar-part diagnostic.
//!   D  determinism: the ENTIRE pipeline runs twice; the BLAKE3 fingerprint
//!      of every gate number and the package Merkle root must be
//!      bit-identical.
//!
//! Epistemic notice (binding): every PASS certifies a WITHIN-MODEL property
//! of a discretisation and its port. This validates numerical engineering
//! on a speculative theory's functional — never anything about nature.

use fs_blake3::ContentHash;
use fs_evidence::{Evidence, ProvenanceHash};
use fs_gum_field::radial::{
    degree, eps_dial, monotone_decreasing, radial_solve, tail_fit, RadialProfile,
};
use fs_gum_field::{
    bps_floor, mc_currents, measure, stencil, virial, Field3, Scheme, Sectors, T_FROZEN,
};
use fs_ivl::Interval;
use fs_math::det;
use fs_package::origin::{
    SourceCertificateRequest, SourceCertificateVerifier, VerificationCapabilities,
    VerificationDecision,
};
use fs_package::{Claim, EvidencePackage, Provenance};

use core::f64::consts::PI;
use std::time::Instant;

// ---------------------------------------------------------------------------
// referee numbers (frozen by earlier tiers; sources in comments)
// ---------------------------------------------------------------------------

/// Step-1 radial references used by field3d's rel errors (6-dp constants,
/// exactly the Python REF dict).
const REF6: [(&str, f64); 5] = [
    ("E2", 12.117459),
    ("E4", 6.191126),
    ("E6", 1.523936),
    ("E0", 1.507588),
    ("I", 22.26826),
];

/// radial_results.json N_run (N = 4000, t frozen), full precision.
const RAD_REF: [f64; 4] = [12.117459252278454, 6.191126358808831, 1.523935580428207, 1.50758847423602];
const RAD_DEG_REF: f64 = 1.0000002921260709;
const RAD_MU_FIT_REF: f64 = 7.772625488362893;
const RAD_BPS_DEFICIT_REF: f64 = 0.004817757679630885;

/// fs-cosserat-pilot certified grid decimals (RESULTS.md, 13/13, golden
/// root f88731af...).
const PILOT_DEG48: f64 = 0.9830385118051;
const PILOT_DEG96: f64 = 0.9956367241585;
const PILOT_E0_96: f64 = 1.508472719913;
const PILOT_E6_96: f64 = 1.502173763242;
const PILOT_I_96: f64 = 21.66361366437;

/// field3d_results.json gates: central4 N=96/N=64 rel errors (E2, E4, E6)
/// and degrees; corner N=96/N=64 the same (the O(h^2) fingerprint).
const C4_96_REL: [f64; 3] = [-0.0026901517531369157, -0.00027140175595197213, 0.0004057739877256772];
const C4_96_DEG: f64 = 0.99991141809055;
const C4_96_VIRIAL_REL: f64 = -0.0006615759772494911;
const CO_96_REL: [f64; 3] = [-0.01946776571657083, -0.014456361235233306, -0.022138143230497564];
const CO_96_DEG: f64 = 0.9860115688934769;

/// Pilot tolerance constants (grid-error certification bounds).
const TOL_DEG_N48: f64 = 3.0e-2;
const TOL_DEG_N96: f64 = 8.0e-3;
const TOL_E0_N96: f64 = 2.0e-3;
const TOL_E6_N96: f64 = 2.0e-2;
const TOL_I_N96: f64 = 2.0e-2;

/// Compacton half-box (pilot) and eps = 0.05 box half-width (field3d).
const HALF_COMPACTON: f64 = 2.4;
const LBOX: f64 = 4.5;

// ---------------------------------------------------------------------------
// certificate plumbing (fs-cosserat-pilot pattern, domain gum-core:field:v1)
// ---------------------------------------------------------------------------

fn certificate_hash(id: &str, statement: &str, lo: f64, hi: f64) -> ContentHash {
    let canon = format!("{id}|{statement}|{:016x}|{:016x}", lo.to_bits(), hi.to_bits());
    fs_blake3::hash_domain("gum-core:field:v1", canon.as_bytes())
}

struct FieldCertVerifier;

impl SourceCertificateVerifier for FieldCertVerifier {
    fn verify(&self, req: &SourceCertificateRequest<'_>) -> VerificationDecision {
        let policy = fs_blake3::hash_domain("gum-core:field:policy:v1", b"recompute-and-compare");
        let expect = certificate_hash(req.claim_id, req.statement, req.lo, req.hi);
        if req.producer == "fs-gum-field" && req.certificate_hash == expect {
            VerificationDecision::accept(policy)
        } else {
            VerificationDecision::reject(policy)
        }
    }
}

// ---------------------------------------------------------------------------
// small helpers
// ---------------------------------------------------------------------------

fn ipi() -> Interval {
    Interval::new(PI, PI.next_up())
}

fn ii(x: i64) -> Interval {
    Interval::point(x as f64)
}

fn ipt(x: f64) -> Interval {
    Interval::point(x)
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
    enclosure: Interval,
    target: f64,
    pass: bool,
    note: String,
) -> Check {
    Check { id, statement: statement.to_string(), enclosure, target, pass, note }
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
// one full deterministic pass
// ---------------------------------------------------------------------------

#[allow(clippy::too_many_lines)]
fn run(verbose: bool) -> RunOutput {
    let mut fp: Vec<f64> = Vec::new(); // determinism fingerprint values
    let mut checks: Vec<Check> = Vec::new();
    let t0 = Instant::now();
    let say = |msg: &str, verbose: bool| {
        if verbose {
            println!("{msg}");
        }
    };

    // ===================== A: eps = 0 compacton vs pilot =====================
    let mut fc48 = Field3::new(48, HALF_COMPACTON);
    fc48.sample_compacton();
    let mut fc96 = Field3::new(96, HALF_COMPACTON);
    fc96.sample_compacton();
    let a48 = measure(&fc48, Scheme::Central2, 0.0);
    let a96 = measure(&fc96, Scheme::Central2, 0.0);
    fp.extend_from_slice(&[a48.deg, a48.e0, a48.e6, a48.i, a96.deg, a96.e0, a96.e6, a96.i]);
    say(&format!("[A] compacton central2 done ({:.1}s)", t0.elapsed().as_secs_f64()), verbose);

    let sqrt2 = ii(2).sqrt();
    let rstar3 = ii(32).sqrt(); // R*^3 = 2^(5/2), exact enclosure
    let e0_ref = ii(16) * sqrt2 / ii(15); // certified closed form 16 sqrt2/15
    let i_ref = ii(128) * ipi() * rstar3 / ii(105); // 128 pi R*^3 / 105

    let err48 = (a48.deg - 1.0).abs();
    let err96 = (a96.deg - 1.0).abs();
    let order_p = (err48 / err96).ln() / 2.0_f64.ln();
    checks.push(check(
        "A1",
        "stored-field compacton degree, N=48, order-2 stencil: |B-1| < 3.0e-2 and reproduces the pilot certified decimal to < 1e-9 rel",
        ipt(a48.deg),
        1.0,
        err48 < TOL_DEG_N48 && rel(a48.deg, PILOT_DEG48).abs() < 1.0e-9,
        format!("|B-1| = {err48:.3e}; vs pilot 0.9830385118051 rel {:+.2e}", rel(a48.deg, PILOT_DEG48)),
    ));
    checks.push(check(
        "A2",
        "stored-field compacton degree, N=96: |B-1| < 8.0e-3 and reproduces the pilot certified decimal to < 1e-9 rel",
        ipt(a96.deg),
        1.0,
        err96 < TOL_DEG_N96 && rel(a96.deg, PILOT_DEG96).abs() < 1.0e-9,
        format!("|B-1| = {err96:.3e}; vs pilot 0.9956367241585 rel {:+.2e}", rel(a96.deg, PILOT_DEG96)),
    ));
    checks.push(check(
        "A3",
        "compacton degree-error N-scaling: error strictly decreases from N=48 to N=96 with measured order p > 1 (O(h^2) expected)",
        ipt(order_p),
        2.0,
        err96 < err48 && order_p > 1.0,
        format!("measured order p = {order_p:.2} (err48 {err48:.3e} -> err96 {err96:.3e})"),
    ));
    let e0_wide = Interval::new(e0_ref.lo() - TOL_E0_N96, e0_ref.hi() + TOL_E0_N96);
    checks.push(check(
        "A4",
        "compacton grid E0(N=96) within the certified enclosure 16 sqrt2/15 widened by 2.0e-3, and matches the pilot decimal to < 1e-9 rel",
        ipt(a96.e0),
        e0_ref.midpoint(),
        e0_wide.contains(a96.e0) && rel(a96.e0, PILOT_E0_96).abs() < 1.0e-9,
        format!("grid - ref = {:+.3e}; vs pilot 1.508472719913 rel {:+.2e}", a96.e0 - e0_ref.midpoint(), rel(a96.e0, PILOT_E0_96)),
    ));
    let e6_wide = Interval::new(e0_ref.lo() - TOL_E6_N96, e0_ref.hi() + TOL_E6_N96);
    checks.push(check(
        "A5",
        "compacton grid E6(N=96) within the certified enclosure (E6 = E0 at the BPS point) widened by 2.0e-2, and matches the pilot decimal to < 1e-9 rel",
        ipt(a96.e6),
        e0_ref.midpoint(),
        e6_wide.contains(a96.e6) && rel(a96.e6, PILOT_E6_96).abs() < 1.0e-9,
        format!("grid - ref = {:+.3e}; vs pilot 1.502173763242 rel {:+.2e}", a96.e6 - e0_ref.midpoint(), rel(a96.e6, PILOT_E6_96)),
    ));
    let i_wide = Interval::new(i_ref.lo() - TOL_I_N96, i_ref.hi() + TOL_I_N96);
    checks.push(check(
        "A6",
        "compacton grid I(N=96) within the certified enclosure 128 pi R*^3/105 widened by 2.0e-2, and matches the pilot decimal to < 1e-9 rel",
        ipt(a96.i),
        i_ref.midpoint(),
        i_wide.contains(a96.i) && rel(a96.i, PILOT_I_96).abs() < 1.0e-9,
        format!("grid - ref = {:+.3e}; vs pilot 21.66361366437 rel {:+.2e}", a96.i - i_ref.midpoint(), rel(a96.i, PILOT_I_96)),
    ));

    // ===================== R: 1-D radial solve at frozen t ===================
    let prof: RadialProfile = radial_solve(T_FROZEN, 4000, 6.0);
    let [re2, re4, re6, re0] = prof.sectors;
    let etot = T_FROZEN * (re2 + re4) + re6 + re0;
    let ratio = T_FROZEN * (re2 + re4) / (re6 + re0);
    let rvir = T_FROZEN * re2 - T_FROZEN * re4 - 3.0 * re6 + 3.0 * re0;
    let rvir_rel = rvir / etot;
    let rdeg = degree(&prof.r, &prof.f);
    let mu_true = 1.0 / (2.0 * T_FROZEN).sqrt();
    let tf = tail_fit(&prof.r, &prof.f, 0.3 * mu_true, 2.5 * mu_true);
    let deficit = (re6 + re0) / bps_floor() - 1.0;
    let mono = monotone_decreasing(&prof.f);
    fp.extend_from_slice(&[re2, re4, re6, re0, ratio, rvir, rdeg, tf.mu, tf.a, tf.rms, deficit, prof.gmax]);
    say(&format!("[R] radial solve done, gmax = {:.1e} ({:.1}s)", prof.gmax, t0.elapsed().as_secs_f64()), verbose);

    let worst_rad = prof
        .sectors
        .iter()
        .zip(RAD_REF.iter())
        .map(|(&a, &b)| rel(a, b).abs())
        .fold(0.0_f64, f64::max);
    checks.push(check(
        "R1",
        "radial sectors (N=4000, t frozen) reproduce radial_results.json N_run E2/E4/E6/E0 to < 1e-6 rel; Newton reaches |g|_max < 1e-11",
        ipt(worst_rad),
        0.0,
        worst_rad < 1.0e-6 && prof.gmax < 1.0e-11,
        format!("worst sector rel {worst_rad:.2e}; E2 = {re2:.9}, E4 = {re4:.9}, E6 = {re6:.9}, E0 = {re0:.9}; gmax = {:.1e}", prof.gmax),
    ));
    checks.push(check(
        "R2",
        "eps dial value: ratio t (E2+E4)/(E6+E0) at the frozen t is 0.05 within 5e-5 (the Tier-2 dial tolerance)",
        ipt(ratio),
        0.05,
        (ratio - 0.05).abs() < 5.0e-5,
        format!("ratio = {ratio:.9} (ref 0.04998469914518915)"),
    ));
    checks.push(check(
        "R3",
        "weighted Derrick virial on the 1-D solution: |t E2 - t E4 - 3 E6 + 3 E0| / E_total < 1e-5",
        ipt(rvir_rel),
        0.0,
        rvir_rel.abs() < 1.0e-5,
        format!("virial rel = {rvir_rel:+.3e} (ref +2.38e-6); E_total = {etot:.10}"),
    ));
    checks.push(check(
        "R4",
        "Yukawa tail: fitted mu equals the derived linearisation mu = 1/sqrt(2t) to < 1e-4 rel and the recorded mu_fit to < 1e-5 rel",
        ipt(tf.mu),
        mu_true,
        rel(tf.mu, mu_true).abs() < 1.0e-4 && rel(tf.mu, RAD_MU_FIT_REF).abs() < 1.0e-5,
        format!("mu_fit = {:.10} vs 1/sqrt(2t) = {mu_true:.10} (rel {:+.2e}); tail rms = {:.2e}, {} pts", tf.mu, rel(tf.mu, mu_true), tf.rms, tf.npts),
    ));
    checks.push(check(
        "R5",
        "Bogomolny: E6 + E0 exceeds the continuum floor 32 sqrt2/15 by the recorded +0.48% deficit (match < 1e-6)",
        ipt(deficit),
        RAD_BPS_DEFICIT_REF,
        deficit > 0.0 && (deficit - RAD_BPS_DEFICIT_REF).abs() < 1.0e-6,
        format!("deficit = {deficit:.9} (ref 0.004817757679630885)"),
    ));
    checks.push(check(
        "R6",
        "1-D degree K = -(2/pi) INT sin^2 f f' dr equals 1 to < 1e-6 and the recorded value to < 1e-9",
        ipt(rdeg),
        1.0,
        (rdeg - 1.0).abs() < 1.0e-6 && (rdeg - RAD_DEG_REF).abs() < 1.0e-9,
        format!("K = {rdeg:.12}"),
    ));
    checks.push(check(
        "R7",
        "solved profile is monotone non-increasing (pi -> 0, no overshoot)",
        ipt(if mono { 1.0 } else { 0.0 }),
        1.0,
        mono,
        String::new(),
    ));
    let (t_dial, dial_prof, dial_hist) = eps_dial(0.008, 0.05, 5.0e-5, 1000, 6.0, 15);
    let dial_ratio = dial_hist.last().map_or(0.0, |&(_, r)| r);
    fp.extend_from_slice(&[t_dial, dial_ratio, dial_prof.gmax]);
    say(&format!("[R] eps dial done: t = {t_dial:.12} in {} outers ({:.1}s)", dial_hist.len(), t0.elapsed().as_secs_f64()), verbose);
    checks.push(check(
        "R8",
        "eps-dial outer loop (N=1000): converges to ratio 0.05 within 5e-5 and lands within 1% of the frozen N=4000 dial t",
        ipt(t_dial),
        T_FROZEN,
        (dial_ratio - 0.05).abs() < 5.0e-5 && rel(t_dial, T_FROZEN).abs() < 1.0e-2,
        format!("t_dial = {t_dial:.12} (rel to frozen {:+.2e}), ratio = {dial_ratio:.8}, {} outer iterations", rel(t_dial, T_FROZEN), dial_hist.len()),
    ));

    // ===================== H: 3-D hedgehog, central4 =========================
    let seed = |n: usize| -> Field3 {
        let mut f = Field3::new(n, LBOX);
        f.sample_hedgehog(&prof.r, &prof.f, 1.0);
        f
    };
    let (h48, h64, h96) = (seed(48), seed(64), seed(96));
    let c4_48 = measure(&h48, Scheme::Central4, T_FROZEN);
    let c4_64 = measure(&h64, Scheme::Central4, T_FROZEN);
    let c4_96 = measure(&h96, Scheme::Central4, T_FROZEN);
    say(&format!("[H] central4 N=48/64/96 done ({:.1}s)", t0.elapsed().as_secs_f64()), verbose);
    for s in [&c4_48, &c4_64, &c4_96] {
        fp.extend_from_slice(&[s.e2, s.e4, s.e6, s.e0, s.i, s.deg]);
    }
    let rel_of = |s: &Sectors| -> Vec<f64> { REF6.iter().map(|&(k, v)| rel(sector_of(s, k), v)).collect() };
    let (r64, r96) = (rel_of(&c4_64), rel_of(&c4_96));
    let worst96 = r96.iter().fold(0.0_f64, |m, &v| m.max(v.abs()));
    checks.push(check(
        "H1",
        "central4 N=96 hedgehog sectors vs the Step-1 radial references: worst rel error < 1e-2 (field3d Section-G gate)",
        ipt(worst96),
        0.0,
        worst96 < 1.0e-2,
        REF6.iter().zip(r96.iter()).map(|(&(k, _), &e)| format!("{k} {e:+.2e}")).collect::<Vec<_>>().join("  "),
    ));
    let dg = (c4_96.deg.abs() - 1.0).abs();
    checks.push(check(
        "H2",
        "central4 N=96 degree: |deg - 1| < 5e-3 and matches the recorded field3d value 0.999911 to < 1e-5",
        ipt(c4_96.deg),
        1.0,
        dg < 5.0e-3 && (c4_96.deg - C4_96_DEG).abs() < 1.0e-5,
        format!("deg = {:.9}, |deg-1| = {dg:.2e}", c4_96.deg),
    ));
    let ratios: Vec<f64> = (0..3).map(|s| r64[s].abs() / r96[s].abs().max(1.0e-30)).collect();
    let e0i_96 = r96[3].abs().max(r96[4].abs());
    let e0i_64 = r64[3].abs().max(r64[4].abs());
    checks.push(check(
        "H3",
        "central4 N64/N96 convergence: derivative-sector (E2,E4,E6) error ratios in [2.5, 6.0] (recorded 3.4-4.1); quadrature-exact E0 and I < 1e-7 rel at N=96, < 1e-5 at N=64 (recorded 5.1e-8 / 8.4e-6)",
        ipt(ratios.iter().copied().fold(f64::INFINITY, f64::min)),
        3.5,
        ratios.iter().all(|&r| (2.5..=6.0).contains(&r)) && e0i_96 < 1.0e-7 && e0i_64 < 1.0e-5,
        format!("ratios E2 {:.2}, E4 {:.2}, E6 {:.2}; E0/I worst rel {e0i_96:.1e} (N=96), {e0i_64:.1e} (N=64)", ratios[0], ratios[1], ratios[2]),
    ));
    let v96 = virial(&c4_96, T_FROZEN);
    let v96_rel = v96 / c4_96.estat;
    checks.push(check(
        "H4",
        "weighted Derrick virial on the central4 N=96 hedgehog: |t E2 - t E4 - 3 E6 + 3 E0|/E_static < 2e-3 and matches the recorded -6.6e-4 to < 1e-6",
        ipt(v96_rel),
        0.0,
        v96_rel.abs() < 2.0e-3 && (v96_rel - C4_96_VIRIAL_REL).abs() < 1.0e-6,
        format!("virial rel = {v96_rel:+.6e} (recorded {C4_96_VIRIAL_REL:+.6e})"),
    ));
    let tail = h64.boundary_tail();
    fp.push(tail);
    checks.push(check(
        "H5",
        "boundary tail |q_vec|_max/pi on the N=64 box faces < 1e-6 (decay-to-vacuum ghost policy; recorded 5.5e-11)",
        ipt(tail),
        0.0,
        tail < 1.0e-6,
        format!("tail = {tail:.3e}"),
    ));
    let match96 = r96
        .iter()
        .take(3)
        .zip(C4_96_REL.iter())
        .map(|(&a, &b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    checks.push(check(
        "H6",
        "central4 N=96 rel errors reproduce the recorded field3d referee rel errors (E2/E4/E6) to < 1e-5 abs",
        ipt(match96),
        0.0,
        match96 < 1.0e-5,
        format!("worst |rel - rel_ref| = {match96:.2e}"),
    ));

    // ===================== C: corner scheme fingerprint ======================
    let co_64 = measure(&h64, Scheme::Corner, T_FROZEN);
    let co_96 = measure(&h96, Scheme::Corner, T_FROZEN);
    say(&format!("[C] corner N=64/96 done ({:.1}s)", t0.elapsed().as_secs_f64()), verbose);
    for s in [&co_64, &co_96] {
        fp.extend_from_slice(&[s.e2, s.e4, s.e6, s.e0, s.i, s.deg]);
    }
    let (cr64, cr96) = (rel_of(&co_64), rel_of(&co_96));
    let cworst = cr96.iter().take(3).fold(0.0_f64, |m, &v| m.max(v.abs()));
    let ce0i = cr96[3].abs().max(cr96[4].abs());
    let cmatch = cr96
        .iter()
        .take(3)
        .zip(CO_96_REL.iter())
        .map(|(&a, &b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    checks.push(check(
        "C1",
        "corner N=96: worst derivative-sector rel error < 5e-2 with the recorded ~2.2e-2 O(h^2) fingerprint reproduced to < 1e-5; E0 and I quadrature-exact (< 1e-6)",
        ipt(cworst),
        0.0,
        cworst < 5.0e-2 && cmatch < 1.0e-5 && ce0i < 1.0e-6,
        format!("worst rel {cworst:.2e}; |rel - rel_ref| worst {cmatch:.2e}; E0/I worst {ce0i:.1e}"),
    ));
    let cratios: Vec<f64> = (0..3).map(|s| cr64[s].abs() / cr96[s].abs().max(1.0e-30)).collect();
    checks.push(check(
        "C2",
        "corner N64/N96 convergence ratios (E2,E4,E6) in [1.7, 2.8] — the clean O(h^2) fingerprint (recorded 2.11-2.22)",
        ipt(cratios.iter().copied().fold(f64::INFINITY, f64::min)),
        2.25,
        cratios.iter().all(|&r| (1.7..=2.8).contains(&r)),
        format!("ratios E2 {:.2}, E4 {:.2}, E6 {:.2}", cratios[0], cratios[1], cratios[2]),
    ));
    checks.push(check(
        "C3",
        "corner N=96 degree matches the recorded 0.986012 to < 1e-5 (the documented O(h^2) under-measurement, reported not hidden)",
        ipt(co_96.deg),
        CO_96_DEG,
        (co_96.deg - CO_96_DEG).abs() < 1.0e-5,
        format!("deg = {:.9} (|deg-1| = {:.2e})", co_96.deg, (co_96.deg - 1.0).abs()),
    ));

    // ===================== F: field mechanics ================================
    let mut pf = Field3::new(32, LBOX);
    pf.sample_hedgehog(&prof.r, &prof.f, 1.0);
    // deterministic trig perturbation (no RNG; amplitude 0.02)
    for i in 0..32 {
        for j in 0..32 {
            for k in 0..32 {
                let (x, y, z) = (pf.x(i), pf.x(j), pf.x(k));
                let mut q = pf.get(i, j, k);
                for (a, qa) in q.iter_mut().enumerate() {
                    let ph = 0.7 * a as f64;
                    *qa += 0.02 * det::sin(1.7 * x + ph) * det::cos(1.3 * y - ph) * det::sin(0.9 * z + 0.1);
                }
                pf.set(i, j, k, q);
            }
        }
    }
    let drift1 = pf.renormalize();
    let drift2 = pf.renormalize();
    fp.extend_from_slice(&[drift1, drift2]);
    checks.push(check(
        "F1",
        "renormalize(): deterministic perturbed hedgehog has norm drift > 1e-4; renormalisation is idempotent (second-pass drift < 1e-14)",
        ipt(drift2),
        0.0,
        drift1 > 1.0e-4 && drift2 < 1.0e-14,
        format!("first-pass drift {drift1:.3e}, second-pass {drift2:.3e}"),
    ));
    let mut mc_worst = 0.0_f64;
    let mut sc_worst = 0.0_f64;
    for i in 0..48 {
        for j in 0..48 {
            for k in 0..48 {
                let q = h48.get(i, j, k);
                let d = stencil::derivs4(&h48, i, j, k);
                let a = mc_currents(q, &d);
                for ax in 0..3 {
                    let na = a[ax][0] * a[ax][0] + a[ax][1] * a[ax][1] + a[ax][2] * a[ax][2] + a[ax][3] * a[ax][3];
                    let nd = d[ax][0] * d[ax][0] + d[ax][1] * d[ax][1] + d[ax][2] * d[ax][2] + d[ax][3] * d[ax][3];
                    let e = (na - nd).abs() / nd.max(1.0e-30);
                    if e > mc_worst {
                        mc_worst = e;
                    }
                    let s = a[ax][0].abs();
                    if s > sc_worst {
                        sc_worst = s;
                    }
                }
            }
        }
    }
    fp.extend_from_slice(&[mc_worst, sc_worst]);
    checks.push(check(
        "F2",
        "Maurer-Cartan currents a_i = conj(q) D_i q on the N=48 hedgehog: |a_i|^2 = |D_i q|^2 to < 1e-12 rel (unit-norm identity)",
        ipt(mc_worst),
        0.0,
        mc_worst < 1.0e-12,
        format!("worst identity rel {mc_worst:.2e}; scalar-part diagnostic max |Re a_i| = {sc_worst:.2e} (continuum value 0; measures FD tangency error)"),
    ));

    // ============== evidence + package + three checker modes ================
    let mut pkg = EvidencePackage::new(Provenance::new(
        "fs-gum-field-v0.1.0-phaseB1-gates",
        "frankensim-fs-math-fs-ga-fs-ivl-fs-evidence-fs-package-fs-checker-std-only",
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
        match ev.certified() {
            Ok(_cert) => {
                certified += 1;
                pkg = pkg.with_claim(Claim::from_certificate(
                    c.id,
                    &c.statement,
                    c.enclosure.lo(),
                    c.enclosure.hi(),
                    "fs-gum-field",
                    certificate_hash(c.id, &c.statement, c.enclosure.lo(), c.enclosure.hi())
                        .to_hex(),
                ));
            }
            Err(e) => println!("  [fs-evidence refused {}: {e:?}]", c.id),
        }
    }
    let root = pkg.try_merkle_root().expect("package root");
    let root_hex = root.to_hex();

    let mode_deny = fs_checker::check_against_root(&pkg, root).passed();
    let verifier = FieldCertVerifier;
    let caps = VerificationCapabilities::deny_all().with_source_certificates(&verifier);
    let mode_cert = fs_checker::check_with_capabilities(&pkg, Some(root), None, &caps).passed();
    let mut wrong = root;
    wrong.0[0] ^= 0xff;
    let mode_tamper =
        fs_checker::check_with_capabilities(&pkg, Some(wrong), None, &caps).passed();

    // determinism fingerprint over every gate number, fixed order
    let mut bytes = Vec::with_capacity(fp.len() * 8);
    for v in &fp {
        bytes.extend_from_slice(&v.to_bits().to_le_bytes());
    }
    let fingerprint_hex = fs_blake3::hash_domain("gum-core:field:fingerprint:v1", &bytes).to_hex();

    RunOutput { checks, certified, root_hex, fingerprint_hex, mode_deny, mode_cert, mode_tamper }
}

fn main() {
    println!("== GUM CORE PHASE B1: fs-gum-field gates ==");
    println!(
        "stored quaternion field, SoA + 2-layer vacuum ghosts; schemes central2/central4/corner;\n\
         compacton box [-{HALF_COMPACTON}, {HALF_COMPACTON}]^3, hedgehog box [-{LBOX}, {LBOX}]^3, t = {T_FROZEN:.15}\n"
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
            "{:<4} {}  [{:.12e}, {:.12e}]  target {:.6e}  {}",
            c.id,
            if c.pass { "PASS" } else { "FAIL" },
            c.enclosure.lo(),
            c.enclosure.hi(),
            c.target,
            c.statement
        );
        if !c.note.is_empty() {
            println!("        note: {}", c.note);
        }
    }
    println!("\n== {passes}/{} gates PASS ==\n", first.checks.len());

    println!(
        "fs-evidence: {}/{} claims certified (fail-closed Certified<f64>; failing gates are \
         excluded from the package by construction).",
        first.certified,
        first.checks.len()
    );
    println!("fs-package Merkle root: {}", first.root_hex);
    println!("run fingerprint (BLAKE3 over all gate numbers): {}", first.fingerprint_hex);
    println!(
        "fs-checker mode 1 (deny-all):        passed = {} (expected false — unauthenticated \
         Verified claims are refused; anti-laundering works)",
        first.mode_deny
    );
    println!("fs-checker mode 2 (cert capability): passed = {} (expected true)", first.mode_cert);
    println!("fs-checker mode 3 (tampered root):   passed = {} (expected false)", first.mode_tamper);

    // Replay determinism: the ENTIRE pipeline (radial solves, 3-D sweeps,
    // fits, certificates, Merkle assembly) runs a second time and must
    // reproduce both the fingerprint and the root bit-for-bit.
    let second = run(false);
    let replay_ok =
        second.root_hex == first.root_hex && second.fingerprint_hex == first.fingerprint_hex;
    println!(
        "replay determinism: second full run root {} first, fingerprint {} first",
        if second.root_hex == first.root_hex { "==" } else { "!=" },
        if second.fingerprint_hex == first.fingerprint_hex { "==" } else { "!=" },
    );
    println!("total gate runtime: {:.1} s", t0.elapsed().as_secs_f64());

    let all_pass = passes == first.checks.len()
        && first.certified == first.checks.len()
        && !first.mode_deny
        && first.mode_cert
        && !first.mode_tamper
        && replay_ok;
    assert!(replay_ok, "replay determinism violated");
    println!(
        "\nOVERALL: {}",
        if all_pass { "PASS (all gates, all checker modes, replay)" } else { "FAIL" }
    );

    println!(
        "\nEPISTEMIC NOTICE: every PASS certifies a WITHIN-MODEL property of a discretisation\n\
         and its port — numerical engineering on a speculative theory's functional. Nothing\n\
         here says anything about nature."
    );
}
