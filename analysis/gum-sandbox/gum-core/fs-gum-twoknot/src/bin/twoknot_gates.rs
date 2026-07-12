//! Pre-campaign gates for fs-gum-twoknot.
//!
//! T1 (twin bit-identity, LOAD-BEARING): on a cubic N = 24 grid the
//!    anisotropic engine must reproduce fs-gum-statics `eval` BITWISE —
//!    every Out field and every gradient component — for the bare corner
//!    and central4 objectives AND for the guard-forced objective (deg_ref
//!    = 1 puts the degree below the band, fgap_ref = +10 penetrates the
//!    wall everywhere: both one-sided branches active — the field3d
//!    gradcheck protocol).
//! T2 (twin ANF bit-identity): 25 guarded ANF iterations from the same
//!    perturbed seed must agree bitwise in (obj, estat, deg, gnorm) at
//!    every recorded step and in every final field component.
//! T3 (degree diagnostic): the two-knot product-ansatz seed carries
//!    topological degree ~ 2 in every channel (|deg - 2| < 0.05 at the
//!    gate grid), and the single-knot seed degree ~ 1.
//! T4 (replay): T1's fingerprint values reproduce bit-for-bit on a second
//!    in-process run.

use fs_gum_field::radial::radial_solve;
use fs_gum_field::{Field3, Scheme, T_FROZEN};
use fs_gum_statics::diag::{add_scaled, bump_field, save_cells};
use fs_gum_statics::{anf, AnfParams, Opts};
use fs_gum_twoknot::{
    anf_a, eval_a, seed_single_knot, seed_two_knot, AnfParamsA, Channel, FieldA, OptsA, H_FROZEN,
};

const N_GATE: usize = 24;
const LBOX: f64 = 4.5;

/// Copy the real cells of a cubic Field3 into an equal-size FieldA.
fn to_field_a(f: &Field3) -> FieldA {
    let n = f.n();
    let h = f.h();
    let mut fa = FieldA::new(n, n, n, h);
    assert_eq!(fa.h().to_bits(), h.to_bits(), "h must match bitwise");
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                fa.set(i, j, k, f.get(i, j, k));
            }
        }
    }
    fa
}

fn bits_eq(a: f64, b: f64) -> bool {
    a.to_bits() == b.to_bits()
}

#[allow(clippy::too_many_lines)]
fn run_gates() -> (bool, Vec<f64>) {
    let mut all_pass = true;
    let mut fingerprint: Vec<f64> = Vec::new();
    let mut gate = |id: &str, msg: &str, ok: bool| {
        println!("[{}] {} — {}", if ok { "PASS" } else { "FAIL" }, id, msg);
        all_pass &= ok;
    };

    // shared seed configuration: hedgehog + LCG bump (statics G-A protocol)
    let rp = radial_solve(T_FROZEN, 4000, 6.0);
    let mut f3 = Field3::new(N_GATE, LBOX);
    f3.sample_hedgehog(&rp.r, &rp.f, 1.0);
    let bump = bump_field(N_GATE, LBOX, 424242, 6, 0.02);
    add_scaled(&mut f3, &bump, 1.0);
    let _ = f3.renormalize();
    let fa = to_field_a(&f3);

    // ---- T1: eval twin bit-identity ---------------------------------------
    for (label, scheme) in [("corner", Scheme::Corner), ("central4", Scheme::Central4)] {
        for (glabel, guards) in [("bare", false), ("guards-forced", true)] {
            let mut o3 = Opts::estatic(scheme);
            let mut oa = OptsA::estatic(scheme);
            if guards {
                o3 = o3.with_guards(1.0, 10.0);
                oa = oa.with_guards(1.0, 1.0, 10.0);
            }
            let (out3, g3) = fs_gum_statics::eval(&f3, &o3, true);
            let (outa, ga) = eval_a(&fa, &oa, true);
            let fields = [
                (out3.e2, outa.e2),
                (out3.e4, outa.e4),
                (out3.e6, outa.e6),
                (out3.e0, outa.e0),
                (out3.i, outa.i),
                (out3.deg, outa.deg),
                (out3.estat, outa.estat),
                (out3.floor_gap, outa.floor_gap),
                (out3.epen, outa.epen),
                (out3.efpen, outa.efpen),
                (out3.obj, outa.obj),
            ];
            let mut ok = fields.iter().all(|&(a, b)| bits_eq(a, b));
            let g3 = g3.expect("grad");
            let ga = ga.expect("grad");
            ok &= g3.data.len() == ga.data.len();
            let mut gdiff = 0usize;
            for (a, b) in g3.data.iter().zip(ga.data.iter()) {
                if !bits_eq(*a, *b) {
                    gdiff += 1;
                }
            }
            ok &= gdiff == 0;
            fingerprint.push(outa.obj);
            fingerprint.push(outa.deg);
            gate(
                "T1",
                &format!(
                    "eval twin {label}/{glabel}: Out bitwise {}, grad mismatches {gdiff}/{}",
                    if fields.iter().all(|&(a, b)| bits_eq(a, b)) { "yes" } else { "NO" },
                    ga.data.len()
                ),
                ok,
            );
        }
    }

    // ---- T2: ANF twin bit-identity (25 guarded iterations) ----------------
    {
        let (outh, _) = fs_gum_statics::eval(&f3, &Opts::estatic(Scheme::Corner), false);
        let deg_ref = outh.deg;
        let fgap_ref = outh.floor_gap - fs_gum_statics::FLOOR_BAND;
        // NOTE floor_gap from bare Opts has deg_ref = 1.0; recompute as the
        // statics gate does: fgap = e6 + e0 - FLOOR (deg/deg_ref = 1).
        let fgap_h = outh.e6 + outh.e0 - fs_gum_field::bps_floor();
        let fgap_ref = {
            let _ = fgap_ref;
            fgap_h - fs_gum_statics::FLOOR_BAND
        };
        let o3 = Opts::estatic(Scheme::Corner).with_guards(deg_ref, fgap_ref);
        let oa = OptsA::estatic(Scheme::Corner).with_guards(deg_ref, deg_ref, fgap_ref);
        let mut f3r = fs_gum_statics::diag::clone_field(&f3);
        let mut far = to_field_a(&f3);
        let mut p3 = AnfParams::new(25);
        p3.instr_every = 5;
        let mut pa = AnfParamsA::new(25);
        pa.instr_every = 5;
        let r3 = anf(&mut f3r, &o3, &p3, "twin3");
        let ra = anf_a(&mut far, &oa, &pa, "twinA");
        let mut ok = r3.series.len() == ra.series.len();
        for (s3, sa) in r3.series.iter().zip(ra.series.iter()) {
            ok &= s3.it == sa.it
                && bits_eq(s3.obj, sa.obj)
                && bits_eq(s3.estat, sa.estat)
                && bits_eq(s3.deg, sa.deg)
                && bits_eq(s3.gnorm, sa.gnorm)
                && bits_eq(s3.dt, sa.dt)
                && s3.arrests == sa.arrests;
        }
        let n = N_GATE;
        let mut buf3 = vec![0.0_f64; 4 * n * n * n];
        save_cells(&f3r, &mut buf3);
        let mut fdiff = 0usize;
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    let qa = far.get(i, j, k);
                    let q3 = f3r.get(i, j, k);
                    for a in 0..4 {
                        if !bits_eq(qa[a], q3[a]) {
                            fdiff += 1;
                        }
                    }
                }
            }
        }
        ok &= fdiff == 0;
        fingerprint.push(ra.out.obj);
        fingerprint.push(ra.out.estat);
        fingerprint.push(ra.out.deg);
        gate(
            "T2",
            &format!(
                "ANF twin 25 guarded iters: series bitwise {}, field mismatches {fdiff}, final obj {:.12}",
                if r3.series.len() == ra.series.len() { "checked" } else { "LEN-MISMATCH" },
                ra.out.obj
            ),
            ok,
        );
    }

    // ---- T3: product-ansatz degree diagnostic ------------------------------
    {
        // reduced gate grid at the frozen h (transverse tail slightly
        // clipped: degree tolerance 0.05 absorbs it)
        let (nx, nt) = (120usize, 64usize);
        let s = 18.0 * H_FROZEN; // d = 3.375 = 2 * 18 h
        for ch in [Channel::Attract, Channel::Repulse, Channel::Align] {
            let mut fa2 = FieldA::new(nx, nt, nt, H_FROZEN);
            seed_two_knot(&mut fa2, &rp, s, ch);
            let drift = fa2.renormalize();
            let (out, _) = eval_a(&fa2, &OptsA::estatic(Scheme::Corner), false);
            let ok = (out.deg - 2.0).abs() < 0.05 && drift < 1.0e-12;
            fingerprint.push(out.deg);
            gate(
                "T3",
                &format!(
                    "product-ansatz degree, channel {}: deg = {:.6} (target 2, tol 0.05), unit-norm drift {drift:.2e}",
                    ch.as_str(),
                    out.deg
                ),
                ok,
            );
        }
        let mut fa1 = FieldA::new(nx, nt, nt, H_FROZEN);
        seed_single_knot(&mut fa1, &rp);
        let (out1, _) = eval_a(&fa1, &OptsA::estatic(Scheme::Corner), false);
        fingerprint.push(out1.deg);
        gate(
            "T3",
            &format!("single-knot seed degree: {:.6} (target 1, tol 0.05)", out1.deg),
            (out1.deg - 1.0).abs() < 0.05,
        );
    }

    (all_pass, fingerprint)
}

fn main() {
    println!("fs-gum-twoknot gates (twin bit-identity + degree diagnostics)");
    let t0 = std::time::Instant::now();
    let (pass1, fp1) = run_gates();
    let (pass2, fp2) = run_gates();
    let mut replay = fp1.len() == fp2.len();
    for (a, b) in fp1.iter().zip(fp2.iter()) {
        replay &= a.to_bits() == b.to_bits();
    }
    println!(
        "[{}] T4 — two-run in-process replay: {} fingerprint f64s bitwise {}",
        if replay { "PASS" } else { "FAIL" },
        fp1.len(),
        if replay { "equal" } else { "UNEQUAL" }
    );
    let all = pass1 && pass2 && replay;
    println!(
        "gates: {} ({:.1}s)",
        if all { "ALL PASS" } else { "FAILURES PRESENT" },
        t0.elapsed().as_secs_f64()
    );
    std::process::exit(i32::from(!all));
}
