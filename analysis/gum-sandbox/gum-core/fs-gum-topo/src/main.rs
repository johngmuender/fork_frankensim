//! GUM core Phase B2 gates: topological-charge diagnostics referee
//! table (fs-gum-topo). Prints every referee check, wraps the passing
//! ones in the fail-closed Certified/EvidencePackage pipeline (domain
//! "gum-core:topo:v1"), verifies the package in three fs-checker modes,
//! and replays the ENTIRE pipeline a second time asserting a
//! bit-identical Merkle root.
//!
//! Epistemic notice (binding): every PASS certifies a WITHIN-MODEL
//! property of diagnostic machinery on a speculative theory's analytic
//! configurations — never the physics.

use fs_gum_topo::certify::{certify, check, CertOutcome, Check};
use fs_gum_topo::degree::{cube_shell, degree_b_density, hedgehog_charge_geometric, Scheme};
use fs_gum_topo::disclination::{
    disclination_lines, planar_winding, z2_holonomy, GridLoop, Z2, CORE_TOL,
};
use fs_gum_topo::hopf::hopf_whitehead;
use fs_gum_topo::linking::{gauss_linking, hopf_preimage_linking};
use fs_gum_topo::referee::{
    baby_skyrmion_slab, compacton_qfield, disclination_field, hopfion_director, linked_circles,
    mirror_z, rotate_target, uniform_director, unlinked_circles, vector_director,
};
use fs_gum_topo::{hopf_project, Bc, TopoError};
use fs_ivl::Interval;

use core::f64::consts::SQRT_2;

/// Domain half-width — the pilot's grid, so the Central2 digits are
/// directly comparable to the certified goldens.
const HALF: f64 = 2.4;

/// Pilot goldens (fs-cosserat-pilot RESULTS.md D1a/D1b, certified).
const GOLD_B48: f64 = 0.9830385118051;
const GOLD_B96: f64 = 0.9956367241585;
const E6_CLOSED: f64 = 16.0 * SQRT_2 / 15.0;

fn ipt(x: f64) -> Interval {
    Interval::point(x)
}

struct RunOutput {
    checks: Vec<Check>,
    cert: CertOutcome,
}

#[allow(clippy::too_many_lines)]
fn run() -> RunOutput {
    let mut checks: Vec<Check> = Vec::new();

    // ---------------- G1 — S3 degree on stored fields ----------------
    let q48 = compacton_qfield(48, HALF);
    let q96 = compacton_qfield(96, HALF);
    let c2_48 = degree_b_density(&q48, Scheme::Central2);
    let c2_96 = degree_b_density(&q96, Scheme::Central2);
    let c4_48 = degree_b_density(&q48, Scheme::Central4);
    let c4_96 = degree_b_density(&q96, Scheme::Central4);
    let co_48 = degree_b_density(&q48, Scheme::Corner);
    let co_96 = degree_b_density(&q96, Scheme::Corner);

    let d48 = (c2_48.deg - GOLD_B48).abs();
    let d96 = (c2_96.deg - GOLD_B96).abs();
    checks.push(check(
        "G1a",
        "stored-field Central2 degree reproduces the pilot's certified digits at N=48",
        ipt(c2_48.deg),
        GOLD_B48,
        d48 < 1.0e-12,
        format!("B(48) = {:.13} (|delta golden| = {d48:.3e})", c2_48.deg),
    ));
    checks.push(check(
        "G1b",
        "stored-field Central2 degree reproduces the pilot's certified digits at N=96",
        ipt(c2_96.deg),
        GOLD_B96,
        d96 < 1.0e-12,
        format!("B(96) = {:.13} (|delta golden| = {d96:.3e})", c2_96.deg),
    ));
    let e48 = (c2_48.deg - 1.0).abs();
    let e96 = (c2_96.deg - 1.0).abs();
    let p_c2 = (e48 / e96).ln() / 2.0_f64.ln();
    checks.push(check(
        "G1c",
        "Central2 degree error N-scaling: strictly decreasing 48 -> 96 with order > 1",
        ipt(p_c2),
        2.0,
        e96 < e48 && p_c2 > 1.0,
        format!("observed order p = {p_c2:.2} (err48 {e48:.3e} -> err96 {e96:.3e})"),
    ));
    let e96c4 = (c4_96.deg - 1.0).abs();
    checks.push(check(
        "G1d",
        "Central4 stored-field degree gate |deg - 1| < 5e-3 at N=96 (field3d_solve.py gate)",
        ipt(c4_96.deg),
        1.0,
        e96c4 < 5.0e-3,
        format!("|deg-1| = {e96c4:.3e} (N=48: {:.3e})", (c4_48.deg - 1.0).abs()),
    ));
    let e48co = (co_48.deg - 1.0).abs();
    let e96co = (co_96.deg - 1.0).abs();
    let p_co = (e48co / e96co).ln() / 2.0_f64.ln();
    checks.push(check(
        "G1e",
        "Corner (midpoint) stored-field degree: |deg - 1| < 2.5e-2 at N=96 with O(h^2) N-scaling",
        ipt(co_96.deg),
        1.0,
        e96co < 2.5e-2 && e96co < e48co && p_co > 1.0,
        format!("|deg-1| = {e96co:.3e} (N=48: {e48co:.3e}), order p = {p_co:.2}"),
    ));
    let e6_err = (c2_96.e6 - E6_CLOSED).abs();
    checks.push(check(
        "G1f",
        "E6 = pi^3 INT b^2 from the returned b-field matches 16 sqrt2/15 within 2e-2 (pilot D3)",
        ipt(c2_96.e6),
        E6_CLOSED,
        e6_err < 2.0e-2,
        format!("E6(96) = {:.12} (drift {:+.3e}; pilot grid value 1.502173763242)", c2_96.e6, c2_96.e6 - E6_CLOSED),
    ));
    let vdir = vector_director(&q48, [0.0, 0.0, 1.0]);
    let shell = cube_shell([0.0; 3], 0.8, 16);
    let qgeo = hedgehog_charge_geometric(&vdir, &shell);
    let qsnap = qgeo.round();
    checks.push(check(
        "G1g",
        "geometric VOS solid-angle hedgehog charge on a cube shell snaps to the integer 1",
        ipt(qgeo),
        1.0,
        qsnap == 1.0 && (qgeo - 1.0).abs() < 1.0e-9,
        format!("Q_raw = {qgeo:.12e} offset from 1 by {:.3e} (lattice-exact snap)", (qgeo - 1.0).abs()),
    ));

    // ---------------- G2 — Hopf map / director projection ----------------
    let n96 = hopf_project(&q96);
    let mut worst = 0.0_f64;
    for v in &n96.data {
        let d = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt() - 1.0;
        if d.abs() > worst {
            worst = d.abs();
        }
    }
    checks.push(check(
        "G2a",
        "hopf_project(q) lands on S2: max | |n| - 1 | < 1e-12 over the N=96 compacton",
        ipt(worst),
        0.0,
        worst < 1.0e-12,
        format!("max | |n|-1 | = {worst:.3e}"),
    ));

    // ---------------- G3 — Hopf invariant, Whitehead/FFT route ----------------
    let h64 = hopfion_director(64, HALF, Bc::Periodic);
    let h128 = hopfion_director(128, HALF, Bc::Periodic);
    let wa64 = hopf_whitehead(&h64);
    let wa128 = hopf_whitehead(&h128);
    let (w64h, w64note) = match &wa64 {
        Ok(w) => (w.hopf, format!("H(64) = {:.9}, curl rel err {:.3e}, max|Im| {:.3e}", w.hopf, w.curl_rel_err, w.max_imag)),
        Err(e) => (f64::NAN, format!("ERROR: {e}")),
    };
    let (w128h, w128note) = match &wa128 {
        Ok(w) => (w.hopf, format!("H(128) = {:.9}, curl rel err {:.3e}, mean|B| obstruction {:.3e} of max|B| {:.3e}", w.hopf, w.curl_rel_err, w.mean_b.iter().fold(0.0f64, |a, m| a.max(m.abs())), w.max_b)),
        Err(e) => (f64::NAN, format!("ERROR: {e}")),
    };
    let we64 = (w64h - 1.0).abs();
    let we128 = (w128h - 1.0).abs();
    checks.push(check(
        "G3a",
        "Whitehead-integral Hopf invariant of the analytic Hopf-1 texture, N=64: |H - 1| < 6e-2",
        ipt(w64h),
        1.0,
        we64 < 6.0e-2,
        w64note,
    ));
    checks.push(check(
        "G3b",
        "Whitehead-integral Hopf invariant, N=128: |H - 1| < 2e-2",
        ipt(w128h),
        1.0,
        we128 < 2.0e-2,
        w128note,
    ));
    let p_w = (we64 / we128).ln() / 2.0_f64.ln();
    checks.push(check(
        "G3c",
        "Whitehead error N-scaling: strictly decreasing 64 -> 128 with order > 1",
        ipt(p_w),
        2.0,
        we128 < we64 && p_w > 1.0,
        format!("observed order p = {p_w:.2} (err64 {we64:.3e} -> err128 {we128:.3e})"),
    ));
    let curl128 = wa128.as_ref().map(|w| w.curl_rel_err).unwrap_or(f64::NAN);
    checks.push(check(
        "G3d",
        "discrete-wavevector curl-inverse: curl_h A reproduces the solenoidal B (rel L2 < 5e-2)",
        ipt(curl128),
        0.0,
        curl128 < 5.0e-2,
        format!("rel L2 residual {curl128:.3e} = non-solenoidal O(h^2) part of the discrete B"),
    ));
    let slab = baby_skyrmion_slab(32, HALF);
    let obstructed = matches!(hopf_whitehead(&slab), Err(TopoError::NetFlux { .. }));
    checks.push(check(
        "G3e",
        "net-flux obstruction: the baby-skyrmion slab (net B_z flux on T^3) is REJECTED, never a number",
        ipt(if obstructed { 1.0 } else { 0.0 }),
        1.0,
        obstructed,
        "typed NetFlux error returned before any solve".to_string(),
    ));

    // ---------------- G4 — Hopf invariant, preimage-linking route ----------------
    let (lc1, lc2) = linked_circles(64);
    let (uc1, uc2) = unlinked_circles(64);
    let lk_linked = gauss_linking(&lc1, &lc2);
    let lk_unlinked = gauss_linking(&uc1, &uc2);
    checks.push(check(
        "G4a",
        "polyline Gauss-linking kernel: |Lk| = 1 for the Hopf-link circles, 0 for the unlinked pair",
        ipt(lk_linked),
        1.0,
        (lk_linked.abs() - 1.0).abs() < 1.0e-9 && lk_unlinked.abs() < 1.0e-12,
        format!("Lk(linked) = {lk_linked:.12e}, Lk(unlinked) = {lk_unlinked:.3e}"),
    ));
    let pl64 = hopf_preimage_linking(&h64);
    let pl96 = hopf_preimage_linking(&hopfion_director(96, HALF, Bc::Periodic));
    let pl128 = hopf_preimage_linking(&h128);
    let fmt_pl = |r: &Result<fs_gum_topo::linking::LinkingOut, TopoError>| -> (f64, String) {
        match r {
            Ok(o) => (
                o.hopf,
                format!(
                    "H = {:.12e}, loops {}+{}, segments {}+{}, attempt {}",
                    o.hopf, o.loops[0], o.loops[1], o.segments[0], o.segments[1], o.attempt
                ),
            ),
            Err(e) => (f64::NAN, format!("ERROR: {e}")),
        }
    };
    let (h_pl64, n_pl64) = fmt_pl(&pl64);
    let (h_pl96, n_pl96) = fmt_pl(&pl96);
    let (h_pl128, n_pl128) = fmt_pl(&pl128);
    checks.push(check(
        "G4b",
        "preimage-linking Hopf invariant, N=64: the exact integer 1 (|H - 1| < 1e-6)",
        ipt(h_pl64),
        1.0,
        (h_pl64 - 1.0).abs() < 1.0e-6,
        n_pl64,
    ));
    checks.push(check(
        "G4c",
        "preimage-linking Hopf invariant, N=96 (non-power-of-two, FFT-free route): integer 1",
        ipt(h_pl96),
        1.0,
        (h_pl96 - 1.0).abs() < 1.0e-6,
        n_pl96,
    ));
    checks.push(check(
        "G4d",
        "preimage-linking Hopf invariant, N=128: integer 1",
        ipt(h_pl128),
        1.0,
        (h_pl128 - 1.0).abs() < 1.0e-6,
        n_pl128,
    ));
    let cross = (w128h - h_pl128).abs();
    checks.push(check(
        "G4e",
        "method-vs-method cross-validation at N=128: |H_Whitehead - H_linking| < 2e-2",
        ipt(cross),
        0.0,
        cross < 2.0e-2,
        format!("H_A(128) = {w128h:.9}, H_B(128) = {h_pl128:.9}, |diff| = {cross:.3e}"),
    ));

    // ---------------- G5 — controls ----------------
    let uni = uniform_director(32, HALF, Bc::Periodic);
    let w_uni = hopf_whitehead(&uni).map(|w| w.hopf).unwrap_or(f64::NAN);
    let pl_uni = hopf_preimage_linking(&uni);
    let (h_uni_b, uni_loops) = match &pl_uni {
        Ok(o) => (o.hopf, o.loops[0] + o.loops[1]),
        Err(_) => (f64::NAN, usize::MAX),
    };
    checks.push(check(
        "G5a",
        "H = 0 control: uniform director gives exactly 0 in BOTH methods (no preimage curves)",
        ipt(w_uni),
        0.0,
        w_uni == 0.0 && h_uni_b == 0.0 && uni_loops == 0,
        format!("H_A = {w_uni:.1e}, H_B = {h_uni_b:.1e}, preimage loops = {uni_loops}"),
    ));
    let rot = rotate_target(&h64, [1.0, 2.0, 3.0], 1.2345);
    let w_rot = hopf_whitehead(&rot).map(|w| w.hopf).unwrap_or(f64::NAN);
    let rot_drift = (w_rot - w64h).abs();
    checks.push(check(
        "G5b",
        "isometry invariance: globally O(3)-rotated Hopf-1 texture leaves H unchanged (< 1e-9)",
        ipt(w_rot),
        1.0,
        rot_drift < 1.0e-9,
        format!("|H(rot) - H| = {rot_drift:.3e}"),
    ));
    let mir = mirror_z(&h64);
    let w_mir = hopf_whitehead(&mir).map(|w| w.hopf).unwrap_or(f64::NAN);
    let mir_drift = (w_mir + w64h).abs();
    let pl_mir = hopf_preimage_linking(&mir);
    let h_mir_b = pl_mir.as_ref().map(|o| o.hopf).unwrap_or(f64::NAN);
    checks.push(check(
        "G5c",
        "parity: the z-mirrored texture measures H = -1 in BOTH methods",
        ipt(w_mir),
        -1.0,
        mir_drift < 1.0e-9 && (h_mir_b + 1.0).abs() < 1.0e-6,
        format!("H_A(mirror) = {w_mir:.9} (|H_A(mirror)+H_A| = {mir_drift:.3e}), H_B(mirror) = {h_mir_b:.9}"),
    ));

    // ---------------- G6 — disclination diagnostics ----------------
    let dp = disclination_field(48, HALF, 1.0);
    let dm = disclination_field(48, HALF, -1.0);
    let ring_small = GridLoop::rectangle_xy(18, 29, 18, 29, 24);
    let ring_big = GridLoop::rectangle_xy(6, 41, 6, 41, 24);
    let off_axis = GridLoop::rectangle_xy(30, 40, 30, 40, 24);
    let hol = |f: &fs_gum_topo::DirectorField, lp: &GridLoop| z2_holonomy(f, lp, CORE_TOL);
    let z_small = hol(&dp, &ring_small);
    let z_big = hol(&dp, &ring_big);
    let z_off = hol(&dp, &off_axis);
    let z_uni = hol(&uniform_director(48, HALF, Bc::Vacuum([0.0, 0.0, 1.0])), &ring_small);
    let g6a_pass = z_small == Ok(Z2::NonTrivial)
        && z_big == Ok(Z2::NonTrivial)
        && z_off == Ok(Z2::Trivial)
        && z_uni == Ok(Z2::Trivial);
    checks.push(check(
        "G6a",
        "Z/2 holonomy: nontrivial on every z-axis-encircling loop of the +1/2 line, trivial off-axis and on the defect-free control",
        ipt(if g6a_pass { 1.0 } else { 0.0 }),
        1.0,
        g6a_pass,
        format!("encircling r~0.6: {z_small:?}, encircling r~1.8: {z_big:?}, off-axis: {z_off:?}, uniform: {z_uni:?}"),
    ));
    let nu = [0.0, 0.0, 1.0];
    let wp = planar_winding(&dp, &ring_small, nu, 0.1).unwrap_or(f64::NAN);
    let wp_big = planar_winding(&dp, &ring_big, nu, 0.1).unwrap_or(f64::NAN);
    let wm = planar_winding(&dm, &ring_small, nu, 0.1).unwrap_or(f64::NAN);
    let woff = planar_winding(&dp, &off_axis, nu, 0.1).unwrap_or(f64::NAN);
    let g6b_pass = (wp - 0.5).abs() < 1.0e-9
        && (wp_big - 0.5).abs() < 1.0e-9
        && (wm + 0.5).abs() < 1.0e-9
        && woff.abs() < 1.0e-9;
    checks.push(check(
        "G6b",
        "half-integer planar winding: +0.5 around the +1/2 line at both radii, -0.5 for the -1/2 variant, 0 off-axis",
        ipt(wp),
        0.5,
        g6b_pass,
        format!("k(+1/2, r~0.6) = {wp:.12}, k(+1/2, r~1.8) = {wp_big:.12}, k(-1/2) = {wm:.12}, k(off) = {woff:.1e}"),
    ));
    let frac = |k: f64| (k - k.round()).abs();
    let g6c_pass = ((frac(wp) - 0.5).abs() < 1.0e-9) == (z_small == Ok(Z2::NonTrivial))
        && ((frac(woff) - 0.5).abs() < 1.0e-9) == (z_off == Ok(Z2::NonTrivial));
    checks.push(check(
        "G6c",
        "consistency: planar winding is half-integer exactly where the Z/2 holonomy is nontrivial",
        ipt(if g6c_pass { 1.0 } else { 0.0 }),
        1.0,
        g6c_pass,
        "checked on the encircling and off-axis loops".to_string(),
    ));
    let lines_p = disclination_lines(&dp, CORE_TOL);
    let lines_m = disclination_lines(&dm, CORE_TOL);
    let lines_u = disclination_lines(&uniform_director(48, HALF, Bc::Vacuum([0.0, 0.0, 1.0])), CORE_TOL);
    let describe = |r: &Result<Vec<fs_gum_topo::disclination::DefectLine>, TopoError>| -> String {
        match r {
            Ok(ls) => format!(
                "{} line(s): {:?}",
                ls.len(),
                ls.iter().map(|l| (l.points.len(), l.closed)).collect::<Vec<_>>()
            ),
            Err(e) => format!("ERROR: {e}"),
        }
    };
    let straight = |r: &Result<Vec<fs_gum_topo::disclination::DefectLine>, TopoError>| -> f64 {
        match r {
            Ok(ls) if ls.len() == 1 => {
                let p0 = ls[0].points[0];
                ls[0]
                    .points
                    .iter()
                    .map(|p| (p[0] - p0[0]).abs().max((p[1] - p0[1]).abs()))
                    .fold(0.0, f64::max)
            }
            _ => f64::NAN,
        }
    };
    let st_p = straight(&lines_p);
    let g6d_pass = matches!(&lines_p, Ok(ls) if ls.len() == 1 && !ls[0].closed && ls[0].points.len() == 48)
        && st_p < 1.0e-12
        && matches!(&lines_m, Ok(ls) if ls.len() == 1 && !ls[0].closed)
        && matches!(&lines_u, Ok(ls) if ls.is_empty());
    checks.push(check(
        "G6d",
        "defect-line tracing: exactly one straight boundary-terminated line on the z-axis (+1/2 and -1/2), none for the control; even pierced-face count per cell (lines cannot end) held everywhere",
        ipt(if g6d_pass { 1.0 } else { 0.0 }),
        1.0,
        g6d_pass,
        format!(
            "+1/2: {} (xy straightness {st_p:.1e}); -1/2: {}; control: {}",
            describe(&lines_p),
            describe(&lines_m),
            describe(&lines_u)
        ),
    ));

    let cert = certify(&checks);
    RunOutput { checks, cert }
}

fn main() {
    println!("== GUM CORE PHASE B2 (fs-gum-topo): topological-charge diagnostics gates ==");
    println!(
        "compacton hedgehog f0(r) = 2 arccos(r/R*), R* = 2^(5/6), grid [-{HALF}, {HALF}]^3; \
         Hopf referees on N = 64/96/128, degree on the pilot's N = 48/96 pair\n"
    );

    let t0 = std::time::Instant::now();
    let first = run();
    let mut passes = 0usize;
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
    println!("\n== {passes}/{} diagnostics PASS ==\n", first.checks.len());

    println!(
        "fs-evidence: {}/{} claims certified (fail-closed Certified<f64>; \
         failing checks are excluded from the package by construction).",
        first.cert.certified, first.cert.total
    );
    println!("fs-package Merkle root: {}", first.cert.root_hex);
    println!(
        "fs-checker mode 1 (deny-all):        passed = {} (expected false — anti-laundering)",
        first.cert.mode_deny
    );
    println!(
        "fs-checker mode 2 (cert capability): passed = {} (expected true)",
        first.cert.mode_cert
    );
    println!(
        "fs-checker mode 3 (tampered root):   passed = {} (expected false)",
        first.cert.mode_tamper
    );

    // Replay determinism: the ENTIRE pipeline (grids, FFTs, marching
    // tetrahedra, Gauss sums, certificates, Merkle assembly) runs a
    // second time and must reproduce the root bit-for-bit.
    let second = run();
    let replay_ok = second.cert.root_hex == first.cert.root_hex;
    println!(
        "replay determinism: second full run root {} first ({})",
        if replay_ok { "==" } else { "!=" },
        second.cert.root_hex
    );
    println!("total wall time: {:.1} s", t0.elapsed().as_secs_f64());

    let all_pass = passes == first.checks.len()
        && first.cert.certified == first.checks.len()
        && !first.cert.mode_deny
        && first.cert.mode_cert
        && !first.cert.mode_tamper
        && replay_ok;
    assert!(replay_ok, "replay determinism violated");
    println!(
        "\nOVERALL: {}",
        if all_pass { "PASS (all diagnostics, all checker modes, replay)" } else { "FAIL" }
    );

    println!(
        "\nEPISTEMIC NOTICE: every PASS certifies a WITHIN-MODEL diagnostic of a\n\
         speculative theory's analytic configuration. This validates frankensim's\n\
         evidence machinery and the discretisation — it is not physics validation."
    );
}
