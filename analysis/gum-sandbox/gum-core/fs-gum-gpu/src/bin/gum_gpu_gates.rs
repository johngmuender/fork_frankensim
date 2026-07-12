//! GUM core Phase G2 gate binary for `fs-gum-gpu` — the wgpu f64 backend
//! correctness rig on lavapipe (software Vulkan).
//!
//! Gates (all machine-checked; band values justified in RESULTS.md):
//!
//!   E1  environment: a Vulkan adapter with SHADER_F64 exists and a device
//!       is created (adapter name/driver printed — the environment proof).
//!   E2  fill_ghosts: upload a hedgehog field with deliberately corrupted
//!       ghost cells, run the device kernel, download — ghosts restored to
//!       vacuum and real cells untouched BIT-EXACTLY (pure stores, so this
//!       one is legitimately a bit gate).
//!   K1  forward sums: GPU central4 forward sweep (workgroup partials +
//!       fixed-order CPU combine) vs fs-gum-field's `measure()` on the
//!       eps = 0.05 hedgehog at N = 48 — per-sector relative bands
//!       (tolerance-band class: FMA contraction + reduction reorder).
//!   K2  same-device replay: the forward partials buffer and the gradient
//!       buffer are BIT-IDENTICAL across two dispatches of the same state.
//!   K3  gradient field (bare E_static weights) vs the CPU engine's
//!       analytic gradient at N = 48: L2 and Linf relative bands.
//!   K3b gradient field with BOTH one-sided guards forced active (anchor +
//!       floor wall) vs CPU — exercises the w6s/cstw/e0c guard dressing.
//!   K4  ANF smoke at N = 32: 50 iterations of guarded descent entirely
//!       on-device except the arrest logic (objective readback per
//!       iteration) vs the same 50 iterations of `fs_gum_statics::anf` on
//!       CPU — GPU objective monotone non-increasing, identical
//!       arrest/acceptance pattern, per-checkpoint and endpoint bands.
//!
//! Epistemic notice (binding): a passing gate certifies the GPU port of a
//! discretisation and its backend plumbing on a software Vulkan rig —
//! never anything about nature, and not yet anything about real-GPU
//! drivers (which may place results elsewhere inside the same bands).

use fs_gum_field::radial::radial_solve;
use fs_gum_field::{measure, Field3, Scheme, GHOST, T_FROZEN};
use fs_gum_gpu::{anf_gpu, GpuEngine};
use fs_gum_statics::diag::{add_scaled, bump_field, clone_field};
use fs_gum_statics::{anf, eval, AnfParams, Opts, FLOOR_BAND};

use std::time::Instant;

const LBOX: f64 = 4.5;
const N_SUM: usize = 48; // forward/gradient gate grid
const N_ANF: usize = 32; // ANF smoke grid
const SEED_PERT: u64 = 424242; // campaign perturbation seed (logical id)
const ANF_ITERS: usize = 50;

// Tolerance bands (derivation in RESULTS.md):
/// Sector sums: n ~ 1.1e5 well-conditioned terms; reduction reorder +
/// per-op FMA contraction bound ~ n * eps-class 2.5e-11 worst-case, but the
/// blocked reduction cancels most of it; measured 1e-14-class.  Band 1e-12.
const BAND_SUM_REL: f64 = 1.0e-12;
/// Gradient: per-entry fixed-length (53-term class) resummation, no global
/// growth; measured 1e-13-class in L2.  Bands 1e-11 (L2) / 1e-10 (Linf).
const BAND_GRAD_L2: f64 = 1.0e-11;
const BAND_GRAD_LINF: f64 = 1.0e-10;
/// ANF endpoints after 50 identical-decision iterations: per-iteration
/// 1e-13-class relative state drift, linear-in-iterations growth class;
/// measured 1e-11-class.  Bands 1e-9 relative (energies), 1e-10 abs (deg).
const BAND_ANF_REL: f64 = 1.0e-9;
const BAND_ANF_DEG: f64 = 1.0e-10;

struct Gates {
    rows: Vec<(String, String, bool)>,
}

impl Gates {
    fn new() -> Self {
        Gates { rows: Vec::new() }
    }
    fn gate(&mut self, id: &str, desc: String, pass: bool) {
        println!("  [{}] {}  {}", if pass { "PASS" } else { "FAIL" }, id, desc);
        self.rows.push((id.to_string(), desc, pass));
    }
    fn all_pass(&self) -> bool {
        self.rows.iter().all(|r| r.2)
    }
}

fn rel(a: f64, b: f64) -> f64 {
    (a / b - 1.0).abs()
}

fn bits_equal(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len() && a.iter().zip(b.iter()).all(|(x, y)| x.to_bits() == y.to_bits())
}

/// Raw padded per-component copy of a field (4 * P^3, component-major).
fn raw_of(f: &Field3) -> Vec<f64> {
    let p = f.padded();
    let mut raw = Vec::with_capacity(4 * p * p * p);
    for c in 0..4 {
        for ip in 0..p {
            for jp in 0..p {
                for kp in 0..p {
                    raw.push(f.get_pad(c, ip, jp, kp));
                }
            }
        }
    }
    raw
}

#[allow(clippy::too_many_lines)]
fn main() {
    let t_all = Instant::now();
    let mut g = Gates::new();
    println!(
        "fs-gum-gpu gates (Phase G2) — N_sum={N_SUM}, N_anf={N_ANF}, LBOX={LBOX}, t={T_FROZEN:.15}"
    );

    // ---- radial profile (frozen Step-1 method) ---------------------------
    let t0 = Instant::now();
    let rp = radial_solve(T_FROZEN, 4000, 6.0);
    println!("radial profile ready ({:.1}s)", t0.elapsed().as_secs_f64());

    // ---- E1: environment --------------------------------------------------
    let t0 = Instant::now();
    let mut e48 = match GpuEngine::new(N_SUM, LBOX) {
        Ok(e) => e,
        Err(err) => {
            println!("  [FAIL] E1 device init: {err}");
            std::process::exit(1);
        }
    };
    println!("adapter: {}", e48.adapter_desc());
    g.gate(
        "E1",
        format!(
            "Vulkan adapter with SHADER_F64, device + 8 f64 pipelines created ({:.1}s): {}",
            t0.elapsed().as_secs_f64(),
            e48.adapter_desc()
        ),
        true,
    );

    // ---- E2: fill_ghosts (bit gate — pure stores) --------------------------
    let mut fh = Field3::new(N_SUM, LBOX);
    fh.sample_hedgehog(&rp.r, &rp.f, 1.0);
    let clean = raw_of(&fh);
    let p = fh.padded();
    let p3 = p * p * p;
    let mut corrupted = clean.clone();
    let real = GHOST..N_SUM + GHOST;
    for c in 0..4 {
        for ip in 0..p {
            for jp in 0..p {
                for kp in 0..p {
                    if !(real.contains(&ip) && real.contains(&jp) && real.contains(&kp)) {
                        corrupted[c * p3 + (ip * p + jp) * p + kp] =
                            -7.5 + c as f64 + 0.001 * ip as f64;
                    }
                }
            }
        }
    }
    e48.upload_raw(&corrupted);
    e48.fill_ghosts();
    let restored = e48.download_raw();
    g.gate(
        "E2",
        format!(
            "fill_ghosts restores corrupted rind to vacuum, real cells untouched (bit-exact over 4*P^3 = {} values)",
            4 * p3
        ),
        bits_equal(&restored, &clean),
    );

    // ---- K1: forward sums vs measure() ------------------------------------
    let opts = Opts::estatic(Scheme::Central4);
    let t0 = Instant::now();
    let p1 = e48.forward_partials();
    let t_fwd = t0.elapsed().as_secs_f64();
    let out = e48.out_from_partials(&opts, &p1);
    let s = measure(&fh, Scheme::Central4, T_FROZEN);
    let mut worst = 0.0_f64;
    for (name, a, b) in [
        ("E2", out.e2, s.e2),
        ("E4", out.e4, s.e4),
        ("E6", out.e6, s.e6),
        ("E0", out.e0, s.e0),
        ("I", out.i, s.i),
        ("deg", out.deg, s.deg),
    ] {
        let r = rel(a, b);
        println!("    {name}: gpu {a:.15e} cpu {b:.15e} rel {r:.2e}");
        worst = worst.max(r);
    }
    g.gate(
        "K1",
        format!(
            "GPU forward sums vs fs-gum-field measure(), eps=0.05 hedgehog N={N_SUM} central4: worst rel {worst:.2e} <= {BAND_SUM_REL:.0e} (forward dispatch {t_fwd:.2}s, partials {} KB)",
            p1.len() * 8 / 1024
        ),
        worst <= BAND_SUM_REL,
    );

    // ---- K2: same-device replay (bit gate) ---------------------------------
    let p2 = e48.forward_partials();
    let fwd_replay = bits_equal(&p1, &p2);
    let t0 = Instant::now();
    e48.gradient(&opts, &out);
    let g1 = e48.download_grad();
    let t_grad = t0.elapsed().as_secs_f64();
    e48.gradient(&opts, &out);
    let g2 = e48.download_grad();
    let grad_replay = bits_equal(&g1, &g2);
    g.gate(
        "K2",
        format!(
            "same-device replay bit-identical: forward partials ({} f64) {fwd_replay}, gradient buffer ({} f64) {grad_replay} (flux+gather dispatch {t_grad:.2}s)",
            p1.len(),
            g1.len()
        ),
        fwd_replay && grad_replay,
    );

    // ---- K3: gradient vs CPU (bare E_static) -------------------------------
    let (_, gc) = eval(&fh, &opts, true);
    let gc = gc.expect("cpu gradient").data;
    let (mut d2, mut c2, mut dinf, mut cinf) = (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
    for (a, b) in g1.iter().zip(gc.iter()) {
        let d = a - b;
        d2 += d * d;
        c2 += b * b;
        dinf = dinf.max(d.abs());
        cinf = cinf.max(b.abs());
    }
    let l2_rel = d2.sqrt() / c2.sqrt();
    let linf_rel = dinf / cinf;
    g.gate(
        "K3",
        format!(
            "GPU gradient (flux+gather) vs CPU analytic gradient, E_static N={N_SUM}: L2 rel {l2_rel:.2e} <= {BAND_GRAD_L2:.0e}, Linf rel {linf_rel:.2e} <= {BAND_GRAD_LINF:.0e}"
        ),
        l2_rel <= BAND_GRAD_L2 && linf_rel <= BAND_GRAD_LINF,
    );

    // ---- K3b: gradient with both guards forced active ----------------------
    // force BOTH guards: deg_ref above the measured degree (anchor active),
    // then fgap_ref above the floor gap AS EVALUATED UNDER THAT deg_ref
    // (floor_gap depends on deg_ref, so it must be probed first).
    let probe = Opts::estatic(Scheme::Central4).with_guards(out.deg + 0.05, 0.0);
    let fgap_probe = e48.out_from_partials(&probe, &p1).floor_gap;
    let opts_g = Opts::estatic(Scheme::Central4).with_guards(out.deg + 0.05, fgap_probe + 0.05);
    let out_g = e48.out_from_partials(&opts_g, &p1);
    assert!(out_g.epen > 0.0 && out_g.efpen > 0.0, "guards not active in K3b");
    e48.gradient(&opts_g, &out_g);
    let gg = e48.download_grad();
    let (_, gcg) = eval(&fh, &opts_g, true);
    let gcg = gcg.expect("cpu guarded gradient").data;
    let (mut d2, mut c2, mut dinf, mut cinf) = (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
    for (a, b) in gg.iter().zip(gcg.iter()) {
        let d = a - b;
        d2 += d * d;
        c2 += b * b;
        dinf = dinf.max(d.abs());
        cinf = cinf.max(b.abs());
    }
    let l2g = d2.sqrt() / c2.sqrt();
    let linfg = dinf / cinf;
    g.gate(
        "K3b",
        format!(
            "GPU gradient with anchor+wall guards ACTIVE (epen {:.2e}, efpen {:.2e}) vs CPU: L2 rel {l2g:.2e} <= {BAND_GRAD_L2:.0e}, Linf rel {linfg:.2e} <= {BAND_GRAD_LINF:.0e}",
            out_g.epen, out_g.efpen
        ),
        l2g <= BAND_GRAD_L2 && linfg <= BAND_GRAD_LINF,
    );
    drop(e48);

    // ---- K4: ANF smoke at N = 32 -------------------------------------------
    let t0 = Instant::now();
    let mut fh32 = Field3::new(N_ANF, LBOX);
    fh32.sample_hedgehog(&rp.r, &rp.f, 1.0);
    let (outh, _) = eval(&fh32, &Opts::estatic(Scheme::Central4), false);
    let opts_anf =
        Opts::estatic(Scheme::Central4).with_guards(outh.deg, outh.floor_gap - FLOOR_BAND);
    let mut fseed = clone_field(&fh32);
    let bump = bump_field(N_ANF, LBOX, SEED_PERT, 6, 0.02);
    add_scaled(&mut fseed, &bump, 1.0);
    let _ = fseed.renormalize();

    // CPU reference: identical 50-iteration protocol
    let mut fcpu = clone_field(&fseed);
    let pcpu = AnfParams {
        maxit: ANF_ITERS,
        dt0: 0.01,
        dt_max: 0.05,
        gtol_ratio: 0.0, // disabled: fixed iteration count on both backends
        instr_every: 10,
        print_every: 0,
    };
    let rc = anf(&mut fcpu, &opts_anf, &pcpu, "cpu");
    let t_cpu = t0.elapsed().as_secs_f64();

    let t0 = Instant::now();
    let mut e32 = match GpuEngine::new(N_ANF, LBOX) {
        Ok(e) => e,
        Err(err) => {
            println!("  [FAIL] K4 device init: {err}");
            std::process::exit(1);
        }
    };
    e32.upload(&fseed);
    e32.fill_ghosts();
    let rg = anf_gpu(&mut e32, &opts_anf, ANF_ITERS, 0.01, 0.05);
    let t_gpu = t0.elapsed().as_secs_f64();
    println!(
        "  ANF 50it N={N_ANF}: cpu {t_cpu:.1}s, gpu {t_gpu:.1}s (lavapipe; throughput irrelevant)"
    );

    let mono = rg.objs.windows(2).all(|w| w[1] <= w[0]);
    g.gate(
        "K4a",
        format!(
            "GPU ANF objective monotone non-increasing over {} accepted states (obj {:.7} -> {:.7})",
            rg.objs.len(),
            rg.objs[0],
            rg.objs[rg.objs.len() - 1]
        ),
        mono,
    );
    g.gate(
        "K4b",
        format!(
            "identical descent decisions: arrests gpu {} == cpu {}, iters {} == {}, status {} == {}",
            rg.arrests,
            rc.arrests,
            rg.iters,
            rc.iters,
            rg.status,
            rc.status.as_str()
        ),
        rg.arrests == rc.arrests && rg.iters == rc.iters && rg.status == rc.status.as_str(),
    );
    // per-checkpoint objective agreement (CPU series records every 10 its)
    let mut worst_ckpt = 0.0_f64;
    for rec in &rc.series {
        if rec.it < rg.objs.len() {
            worst_ckpt = worst_ckpt.max(rel(rg.objs[rec.it], rec.obj));
        }
    }
    g.gate(
        "K4c",
        format!(
            "objective checkpoints (its 0,10,...,50): worst rel {worst_ckpt:.2e} <= {BAND_ANF_REL:.0e}"
        ),
        worst_ckpt <= BAND_ANF_REL,
    );
    let mut worst_end = 0.0_f64;
    for (name, a, b) in [
        ("obj", rg.out.obj, rc.out.obj),
        ("estat", rg.out.estat, rc.out.estat),
        ("E2", rg.out.e2, rc.out.e2),
        ("E4", rg.out.e4, rc.out.e4),
        ("E6", rg.out.e6, rc.out.e6),
        ("E0", rg.out.e0, rc.out.e0),
        ("I", rg.out.i, rc.out.i),
    ] {
        let r = rel(a, b);
        println!("    endpoint {name}: gpu {a:.12e} cpu {b:.12e} rel {r:.2e}");
        worst_end = worst_end.max(r);
    }
    let ddeg = (rg.out.deg - rc.out.deg).abs();
    g.gate(
        "K4d",
        format!(
            "endpoint bands: worst energy/I rel {worst_end:.2e} <= {BAND_ANF_REL:.0e}, |deg gpu - cpu| {ddeg:.2e} <= {BAND_ANF_DEG:.0e}"
        ),
        worst_end <= BAND_ANF_REL && ddeg <= BAND_ANF_DEG,
    );

    // ---- summary -----------------------------------------------------------
    println!("\ngate summary:");
    for (id, _desc, pass) in &g.rows {
        println!("  {}  {}", if *pass { "PASS" } else { "FAIL" }, id);
    }
    println!("total wall time {:.1}s", t_all.elapsed().as_secs_f64());
    if !g.all_pass() {
        std::process::exit(1);
    }
    println!("ALL GATES PASS");
}
