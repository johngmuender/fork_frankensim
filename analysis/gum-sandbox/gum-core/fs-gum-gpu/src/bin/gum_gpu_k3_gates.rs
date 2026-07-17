//! GUM core Workstream K3 gate binary for `fs-gum-gpu` — the FULL-pipeline
//! WGSL f64 verification battery on llvmpipe (software Vulkan), extending
//! the Phase G2 suite (`gum_gpu_gates`) per ROADMAP_v7_HARDENING.md.
//!
//! Battery (all machine-checked; bands justified in K3_RESULTS.md):
//!
//!   E1     environment: Vulkan adapter with SHADER_F64, device created
//!          (adapter identity printed — the evidence line).
//!   S{N}   sector measure GPU-vs-CPU: the five sector integrals
//!          (E2, E4, E6, E0, I) plus deg on the eps = 0.05 hedgehog at
//!          N = 24, 32, 48, 64, GPU central4 forward sweep (fixed-order
//!          workgroup partials + fixed-order CPU combine) vs
//!          fs-gum-field `measure()` — per-sector relative band.
//!   G{N}   E_static gradient GPU-vs-CPU (gather form, T_FROZEN weights):
//!          max-abs and RMS relative difference over ALL sites/components
//!          vs the CPU analytic gradient of `fs_gum_statics::eval`.
//!   D{N}   GPU determinism: two dispatches of the same device state are
//!          BIT-IDENTICAL on every output buffer (forward partials and
//!          gradient buffer compared bitwise).
//!   F{N}   golden fingerprints: FNV-1a-64 over the LE bytes of the
//!          forward-partials and gradient buffers, recorded as
//!          within-backend goldens (the replay copy must hash equal).
//!   DD     directional derivative on BOTH backends at N = 32: central FD
//!          of E_static along seeded random tangent directions
//!          (eps = 1e-5, the statics G-A protocol) vs <grad, dir>, with
//!          the FD evaluated by each backend's own forward path; plus the
//!          cross-backend <g_gpu, u> vs <g_cpu, u> agreement.
//!   AD     arrested-descent smoke at N = 24: the descent loop stepping
//!          with the GPU gradient and arresting on the GPU objective;
//!          every ACCEPTED state downloaded and re-evaluated by the CPU
//!          `eval` — the CPU-evaluated objective must be monotone
//!          non-increasing and strictly decrease overall.
//!
//! Wall-time honesty: per-N CPU and GPU timings are RECORDED but are not
//! evidence — llvmpipe is software Vulkan on this 4-core container; the
//! rig demonstrates correctness/readiness, never GPU performance.
//!
//! Determinism: all seeds fixed; physics numbers of two runs of this
//! binary agree bit-for-bit on a given stack (timing fields exempt).
//! Output `k3_results.json` is stage-flushed next to Cargo.toml.
//!
//! Epistemic notice (binding): every PASS certifies a WITHIN-MODEL
//! property of a discretisation, its adjoint, and a backend port on a
//! software Vulkan rig — never anything about nature.

use fs_gum_field::radial::radial_solve;
use fs_gum_field::{measure, Field3, Scheme, T_FROZEN};
use fs_gum_gpu::GpuEngine;
use fs_gum_statics::diag::{
    add_scaled, bump_field, clone_field, dot_flat, normalize_l2, tangent_project,
};
use fs_gum_statics::{eval, Opts, FLOOR_BAND};

use std::fmt::Write as _;
use std::time::Instant;

const LBOX: f64 = 4.5;
/// Battery grids (64 included: llvmpipe wall-time allowed it; see timings).
const NS: [usize; 4] = [24, 32, 48, 64];
const N_DD: usize = 32; // directional-derivative grid
const N_AD: usize = 24; // arrested-descent smoke grid
const AD_ITERS: usize = 40;
const AD_DT0: f64 = 0.01;
const AD_DT_MAX: f64 = 0.05;
const EPS_FD: f64 = 1.0e-5; // statics G-A FD protocol
const SEED_DD_BASE: u64 = 11; // statics gradcheck base-perturbation seed id
const SEED_AD: u64 = 424242; // campaign perturbation seed (logical id)

// Tolerance bands (derivations in K3_RESULTS.md §2/§3/§4):
/// Sector sums: reduction-reorder + FMA-contraction class, n <= 2.6e5
/// well-conditioned terms; G2 band kept, re-verified per N here.
const BAND_SUM_REL: f64 = 1.0e-12;
/// Gradient: per-entry fixed-length (~57-term) resummation, no n-growth.
const BAND_GRAD_RMS: f64 = 1.0e-11;
const BAND_GRAD_MAX: f64 = 1.0e-10;
/// FD-vs-analytic residual at eps = 1e-5 (truncation-dominated; the
/// statics G-A tolerance).
const BAND_DD: f64 = 1.0e-5;
/// Cross-backend <g, u> agreement (fixed-order CPU dot over 4 N^3 entries
/// that differ at the 1e-15-relative class each).
const BAND_DD_CROSS: f64 = 1.0e-9;

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

/// FNV-1a 64-bit over the little-endian bytes of an f64 slice — the
/// dependency-free fingerprint of a GPU output buffer.
fn fnv1a64(v: &[f64]) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325_u64;
    for x in v {
        for b in x.to_le_bytes() {
            h ^= u64::from(b);
            h = h.wrapping_mul(0x0000_0100_0000_01b3);
        }
    }
    h
}

/// Stage-flushed JSON writer: rewrites the whole file from accumulated
/// completed fragments after every section (crash discipline).
struct Json {
    path: std::path::PathBuf,
    frags: Vec<(String, String)>, // (key, already-serialised value)
}

impl Json {
    fn new() -> Self {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("k3_results.json");
        Json { path, frags: Vec::new() }
    }
    fn put(&mut self, key: &str, value: String) {
        self.frags.retain(|(k, _)| k != key);
        self.frags.push((key.to_string(), value));
        self.flush();
    }
    fn flush(&self) {
        let mut s = String::from("{\n");
        for (i, (k, v)) in self.frags.iter().enumerate() {
            let _ = write!(s, "  \"{k}\": {v}");
            s.push_str(if i + 1 < self.frags.len() { ",\n" } else { "\n" });
        }
        s.push_str("}\n");
        std::fs::write(&self.path, s).expect("write k3_results.json");
    }
}

fn jstr(s: &str) -> String {
    format!("\"{}\"", s.replace('\\', "\\\\").replace('"', "\\\""))
}

/// One battery grid's record (accumulated, then serialised).
struct NRec {
    n: usize,
    sector_rows: Vec<(String, f64, f64, f64)>, // name, gpu, cpu, rel
    worst_sector: f64,
    grad_rms_rel: f64,
    grad_max_rel: f64,
    grad_len: usize,
    det_fwd: bool,
    det_grad: bool,
    fp_partials: u64,
    fp_grad: u64,
    t_cpu_measure: f64,
    t_cpu_eval_grad: f64,
    t_gpu_init: f64,
    t_gpu_forward: f64,
    t_gpu_grad: f64,
}

fn nrec_json(r: &NRec) -> String {
    let mut s = String::from("{\n");
    let _ = write!(s, "    \"n\": {},\n    \"sectors\": {{", r.n);
    for (i, (name, a, b, rr)) in r.sector_rows.iter().enumerate() {
        let _ = write!(
            s,
            "\n      \"{name}\": {{\"gpu\": {a:e}, \"cpu\": {b:e}, \"rel\": {rr:e}}}{}",
            if i + 1 < r.sector_rows.len() { "," } else { "" }
        );
    }
    let _ = write!(
        s,
        "\n    }},\n    \"worst_sector_rel\": {:e},\n    \"grad_rms_rel\": {:e},\n    \"grad_maxabs_rel\": {:e},\n    \"grad_entries\": {},\n    \"replay_bit_identical\": {{\"forward_partials\": {}, \"gradient\": {}}},\n    \"fingerprint_fnv1a64\": {{\"forward_partials\": \"{:016x}\", \"gradient\": \"{:016x}\"}},\n    \"timings_s\": {{\"cpu_measure\": {:.3}, \"cpu_eval_grad\": {:.3}, \"gpu_init\": {:.3}, \"gpu_forward_readback\": {:.3}, \"gpu_flux_gather_readback\": {:.3}}}\n  }}",
        r.worst_sector,
        r.grad_rms_rel,
        r.grad_max_rel,
        r.grad_len,
        r.det_fwd,
        r.det_grad,
        r.fp_partials,
        r.fp_grad,
        r.t_cpu_measure,
        r.t_cpu_eval_grad,
        r.t_gpu_init,
        r.t_gpu_forward,
        r.t_gpu_grad
    );
    s
}

#[allow(clippy::too_many_lines)]
fn main() {
    let t_all = Instant::now();
    let mut g = Gates::new();
    let mut js = Json::new();
    println!(
        "fs-gum-gpu K3 gates — grids {NS:?}, DD N={N_DD}, AD N={N_AD}, LBOX={LBOX}, t={T_FROZEN:.15}"
    );
    js.put("workstream", jstr("K3 full-pipeline WGSL f64 battery (llvmpipe)"));
    js.put(
        "protocol",
        format!(
            "{{\"grids\": [24, 32, 48, 64], \"lbox\": {LBOX}, \"t_frozen\": {T_FROZEN:e}, \"scheme\": \"central4\", \"eps_fd\": {EPS_FD:e}, \"bands\": {{\"sector_rel\": {BAND_SUM_REL:e}, \"grad_rms_rel\": {BAND_GRAD_RMS:e}, \"grad_maxabs_rel\": {BAND_GRAD_MAX:e}, \"dirder_rel\": {BAND_DD:e}, \"dirder_cross_rel\": {BAND_DD_CROSS:e}}}}}"
        ),
    );

    // ---- radial profile (frozen Step-1 method) ---------------------------
    let t0 = Instant::now();
    let rp = radial_solve(T_FROZEN, 4000, 6.0);
    println!("radial profile ready ({:.1}s)", t0.elapsed().as_secs_f64());

    let opts = Opts::estatic(Scheme::Central4);

    // ---- E1 + per-N battery ----------------------------------------------
    let mut first = true;
    let mut nrecs: Vec<String> = Vec::new();
    for n in NS {
        println!("---- N = {n} ----");
        // CPU reference
        let mut fh = Field3::new(n, LBOX);
        fh.sample_hedgehog(&rp.r, &rp.f, 1.0);
        let t0 = Instant::now();
        let s = measure(&fh, Scheme::Central4, T_FROZEN);
        let t_cpu_measure = t0.elapsed().as_secs_f64();
        let t0 = Instant::now();
        let (out_cpu, gc) = eval(&fh, &opts, true);
        let t_cpu_eval_grad = t0.elapsed().as_secs_f64();
        let gc = gc.expect("cpu gradient").data;

        // GPU engine
        let t0 = Instant::now();
        let mut eng = match GpuEngine::new(n, LBOX) {
            Ok(e) => e,
            Err(err) => {
                println!("  [FAIL] E1 device init (N={n}): {err}");
                g.gate("E1", format!("device init N={n}: {err}"), false);
                break;
            }
        };
        let t_gpu_init = t0.elapsed().as_secs_f64();
        if first {
            println!("adapter: {}", eng.adapter_desc());
            js.put("adapter", jstr(&eng.adapter_desc()));
            g.gate(
                "E1",
                format!("Vulkan adapter with SHADER_F64, device created: {}", eng.adapter_desc()),
                true,
            );
            first = false;
        }
        eng.upload(&fh);
        eng.fill_ghosts();

        // S{n}: sector measure
        let t0 = Instant::now();
        let p1 = eng.forward_partials();
        let t_gpu_forward = t0.elapsed().as_secs_f64();
        let out_gpu = eng.out_from_partials(&opts, &p1);
        let mut sector_rows = Vec::new();
        let mut worst_sector = 0.0_f64;
        for (name, a, b) in [
            ("E2", out_gpu.e2, s.e2),
            ("E4", out_gpu.e4, s.e4),
            ("E6", out_gpu.e6, s.e6),
            ("E0", out_gpu.e0, s.e0),
            ("I", out_gpu.i, s.i),
            ("deg", out_gpu.deg, s.deg),
        ] {
            let r = rel(a, b);
            println!("    {name}: gpu {a:.15e} cpu {b:.15e} rel {r:.2e}");
            worst_sector = worst_sector.max(r);
            sector_rows.push((name.to_string(), a, b, r));
        }
        // cross-check the two CPU forward paths agree (measure vs eval)
        assert!(rel(out_cpu.e2, s.e2) < 1.0e-12, "cpu measure/eval mismatch");
        g.gate(
            &format!("S{n}"),
            format!(
                "sector measure GPU vs CPU measure(), hedgehog N={n} central4: worst rel {worst_sector:.2e} <= {BAND_SUM_REL:.0e}"
            ),
            worst_sector <= BAND_SUM_REL,
        );

        // G{n}: gradient
        let t0 = Instant::now();
        eng.gradient(&opts, &out_gpu);
        let g1 = eng.download_grad();
        let t_gpu_grad = t0.elapsed().as_secs_f64();
        let (mut d2, mut c2, mut dinf, mut cinf) = (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
        for (a, b) in g1.iter().zip(gc.iter()) {
            let d = a - b;
            d2 += d * d;
            c2 += b * b;
            dinf = dinf.max(d.abs());
            cinf = cinf.max(b.abs());
        }
        let grad_rms_rel = (d2 / g1.len() as f64).sqrt() / (c2 / g1.len() as f64).sqrt();
        let grad_max_rel = dinf / cinf;
        g.gate(
            &format!("G{n}"),
            format!(
                "E_static gradient (gather) GPU vs CPU analytic, N={n}, {} entries: RMS rel {grad_rms_rel:.2e} <= {BAND_GRAD_RMS:.0e}, max-abs rel {grad_max_rel:.2e} <= {BAND_GRAD_MAX:.0e}",
                g1.len()
            ),
            grad_rms_rel <= BAND_GRAD_RMS && grad_max_rel <= BAND_GRAD_MAX,
        );

        // D{n}: determinism (two dispatches, bitwise)
        let p2 = eng.forward_partials();
        let det_fwd = bits_equal(&p1, &p2);
        eng.gradient(&opts, &out_gpu);
        let g2 = eng.download_grad();
        let det_grad = bits_equal(&g1, &g2);
        g.gate(
            &format!("D{n}"),
            format!(
                "GPU determinism N={n}: replay bit-identical — forward partials ({} f64) {det_fwd}, gradient ({} f64) {det_grad}",
                p1.len(),
                g1.len()
            ),
            det_fwd && det_grad,
        );

        // F{n}: golden fingerprints (replay copy must hash equal)
        let fp_partials = fnv1a64(&p1);
        let fp_grad = fnv1a64(&g1);
        let fp_ok = fnv1a64(&p2) == fp_partials && fnv1a64(&g2) == fp_grad;
        g.gate(
            &format!("F{n}"),
            format!(
                "golden fingerprints N={n} (FNV-1a-64, within-backend): partials {fp_partials:016x}, gradient {fp_grad:016x}, replay hashes equal {fp_ok}"
            ),
            fp_ok,
        );

        println!(
            "    timings N={n} [s]: cpu measure {t_cpu_measure:.3}, cpu eval+grad {t_cpu_eval_grad:.3}; gpu init {t_gpu_init:.3}, forward+readback {t_gpu_forward:.3}, flux+gather+readback {t_gpu_grad:.3}  (llvmpipe: NOT a performance claim)"
        );
        let rec = NRec {
            n,
            sector_rows,
            worst_sector,
            grad_rms_rel,
            grad_max_rel,
            grad_len: g1.len(),
            det_fwd,
            det_grad,
            fp_partials,
            fp_grad,
            t_cpu_measure,
            t_cpu_eval_grad,
            t_gpu_init,
            t_gpu_forward,
            t_gpu_grad,
        };
        nrecs.push(nrec_json(&rec));
        js.put("battery", format!("[{}\n  ]", nrecs.join(", ")));
    }

    // ---- DD: directional derivative on both backends at N = 32 ------------
    println!("---- DD (N = {N_DD}) ----");
    let mut fbase = Field3::new(N_DD, LBOX);
    fbase.sample_hedgehog(&rp.r, &rp.f, 1.0);
    let bump = bump_field(N_DD, LBOX, SEED_DD_BASE, 6, 0.05);
    add_scaled(&mut fbase, &bump, 1.0);
    let _ = fbase.renormalize();
    let (out_b, gcb) = eval(&fbase, &opts, true);
    let gcb = gcb.expect("cpu gradient").data;
    let mut edd = GpuEngine::new(N_DD, LBOX).expect("DD engine");
    edd.upload(&fbase);
    edd.fill_ghosts();
    let pb = edd.forward_partials();
    let out_bg = edd.out_from_partials(&opts, &pb);
    edd.gradient(&opts, &out_bg);
    let ggb = edd.download_grad();
    let mut worst_cpu = 0.0_f64;
    let mut worst_gpu = 0.0_f64;
    let mut worst_cross = 0.0_f64;
    let mut dd_rows = Vec::new();
    for dir in 0..3u64 {
        let seed = 7000 + dir;
        let mut u = bump_field(N_DD, LBOX, seed, 6, 1.0);
        tangent_project(&fbase, &mut u);
        let _ = normalize_l2(&mut u);
        let mut fp = clone_field(&fbase);
        add_scaled(&mut fp, &u, EPS_FD);
        let _ = fp.renormalize();
        let mut fm = clone_field(&fbase);
        add_scaled(&mut fm, &u, -EPS_FD);
        let _ = fm.renormalize();
        // CPU FD
        let (op, _) = eval(&fp, &opts, false);
        let (om, _) = eval(&fm, &opts, false);
        let fd_cpu = (op.obj - om.obj) / (2.0 * EPS_FD);
        // GPU FD (each backend's own forward path)
        edd.upload(&fp);
        edd.fill_ghosts();
        let ojp = edd.forward_out(&opts).obj;
        edd.upload(&fm);
        edd.fill_ghosts();
        let ojm = edd.forward_out(&opts).obj;
        let fd_gpu = (ojp - ojm) / (2.0 * EPS_FD);
        let an_cpu = dot_flat(&gcb, &u);
        let an_gpu = dot_flat(&ggb, &u);
        let r_cpu = (fd_cpu - an_cpu).abs() / fd_cpu.abs().max(1.0e-30);
        let r_gpu = (fd_gpu - an_gpu).abs() / fd_gpu.abs().max(1.0e-30);
        let r_cross = (an_gpu - an_cpu).abs() / an_cpu.abs().max(1.0e-30);
        println!(
            "    dir seed {seed}: fd_cpu {fd_cpu:+.12e} an_cpu {an_cpu:+.12e} res {r_cpu:.2e} | fd_gpu {fd_gpu:+.12e} an_gpu {an_gpu:+.12e} res {r_gpu:.2e} | cross {r_cross:.2e}"
        );
        worst_cpu = worst_cpu.max(r_cpu);
        worst_gpu = worst_gpu.max(r_gpu);
        worst_cross = worst_cross.max(r_cross);
        dd_rows.push(format!(
            "{{\"seed\": {seed}, \"fd_cpu\": {fd_cpu:e}, \"an_cpu\": {an_cpu:e}, \"res_cpu\": {r_cpu:e}, \"fd_gpu\": {fd_gpu:e}, \"an_gpu\": {an_gpu:e}, \"res_gpu\": {r_gpu:e}, \"cross_rel\": {r_cross:e}}}"
        ));
    }
    // restore the base field on device (documented state for any follow-on)
    edd.upload(&fbase);
    edd.fill_ghosts();
    drop(edd);
    js.put(
        "dirder",
        format!(
            "{{\"n\": {N_DD}, \"eps_fd\": {EPS_FD:e}, \"base_obj_cpu\": {:e}, \"base_obj_gpu\": {:e}, \"dirs\": [{}], \"worst_res_cpu\": {worst_cpu:e}, \"worst_res_gpu\": {worst_gpu:e}, \"worst_cross_rel\": {worst_cross:e}}}",
            out_b.obj,
            out_bg.obj,
            dd_rows.join(", ")
        ),
    );
    g.gate(
        "DD",
        format!(
            "directional derivative N={N_DD}, 3 seeded tangent dirs, eps={EPS_FD:.0e}: FD-vs-<g,u> worst res cpu {worst_cpu:.2e} / gpu {worst_gpu:.2e} <= {BAND_DD:.0e}; cross-backend <g,u> rel {worst_cross:.2e} <= {BAND_DD_CROSS:.0e}"
        ),
        worst_cpu <= BAND_DD && worst_gpu <= BAND_DD && worst_cross <= BAND_DD_CROSS,
    );

    // ---- AD: arrested-descent smoke with the GPU gradient ------------------
    println!("---- AD (N = {N_AD}, {AD_ITERS} iterations) ----");
    let mut fh24 = Field3::new(N_AD, LBOX);
    fh24.sample_hedgehog(&rp.r, &rp.f, 1.0);
    let (outh, _) = eval(&fh24, &opts, false);
    let opts_ad = Opts::estatic(Scheme::Central4).with_guards(outh.deg, outh.floor_gap - FLOOR_BAND);
    let mut fseed = clone_field(&fh24);
    let bump = bump_field(N_AD, LBOX, SEED_AD, 6, 0.02);
    add_scaled(&mut fseed, &bump, 1.0);
    let _ = fseed.renormalize();

    let t0 = Instant::now();
    let mut ead = GpuEngine::new(N_AD, LBOX).expect("AD engine");
    ead.upload(&fseed);
    ead.fill_ghosts();
    let h3 = ead.h() * ead.h() * ead.h();
    let mut ftmp = Field3::new(N_AD, LBOX); // ghosts vacuum, matches device policy
    let mut cpu_objs: Vec<f64> = Vec::new();
    let mut gpu_objs: Vec<f64> = Vec::new();
    let cpu_eval_of = |e: &GpuEngine, f: &mut Field3| -> f64 {
        e.download_field(f);
        eval(f, &opts_ad, false).0.obj
    };
    let mut out = ead.forward_out(&opts_ad);
    ead.gradient(&opts_ad, &out);
    ead.snapshot_save();
    ead.v_zero();
    cpu_objs.push(cpu_eval_of(&ead, &mut ftmp));
    gpu_objs.push(out.obj);
    let mut dt = AD_DT0;
    let mut arrests = 0usize;
    let mut accepted = 0usize;
    let mut stalled = false;
    for _it in 0..AD_ITERS {
        ead.v_axpy(dt / h3);
        ead.q_step(dt);
        ead.renormalize();
        let out_new = ead.forward_out(&opts_ad);
        if out_new.obj > out.obj {
            ead.snapshot_restore();
            ead.v_zero();
            dt *= 0.6;
            arrests += 1;
            if dt < 1.0e-7 {
                stalled = true;
                break;
            }
        } else {
            out = out_new;
            ead.gradient(&opts_ad, &out);
            ead.snapshot_save();
            ead.v_project();
            dt = (dt * 1.01).min(AD_DT_MAX);
            accepted += 1;
            cpu_objs.push(cpu_eval_of(&ead, &mut ftmp));
            gpu_objs.push(out.obj);
        }
    }
    let t_ad = t0.elapsed().as_secs_f64();
    let mono_cpu = cpu_objs.windows(2).all(|w| w[1] <= w[0]);
    let mono_gpu = gpu_objs.windows(2).all(|w| w[1] <= w[0]);
    let mut worst_obj_rel = 0.0_f64;
    for (a, b) in gpu_objs.iter().zip(cpu_objs.iter()) {
        worst_obj_rel = worst_obj_rel.max(rel(*a, *b));
    }
    let net = cpu_objs[0] - cpu_objs[cpu_objs.len() - 1];
    println!(
        "    accepted {accepted}, arrests {arrests}, stalled {stalled}; CPU-eval obj {:.9} -> {:.9} (net -{net:.3e}); worst GPU-vs-CPU obj rel {worst_obj_rel:.2e}  ({t_ad:.1}s)",
        cpu_objs[0],
        cpu_objs[cpu_objs.len() - 1]
    );
    js.put(
        "descent_smoke",
        format!(
            "{{\"n\": {N_AD}, \"iters\": {AD_ITERS}, \"dt0\": {AD_DT0}, \"dt_max\": {AD_DT_MAX}, \"seed\": {SEED_AD}, \"accepted\": {accepted}, \"arrests\": {arrests}, \"stalled\": {stalled}, \"cpu_obj_first\": {:e}, \"cpu_obj_last\": {:e}, \"cpu_monotone\": {mono_cpu}, \"gpu_monotone\": {mono_gpu}, \"worst_obj_cross_rel\": {worst_obj_rel:e}, \"wall_s\": {t_ad:.1}}}",
            cpu_objs[0],
            cpu_objs[cpu_objs.len() - 1]
        ),
    );
    g.gate(
        "AD-a",
        format!(
            "arrested descent (GPU gradient) N={N_AD}: CPU-evaluated objective monotone non-increasing over {} accepted states ({:.9} -> {:.9})",
            cpu_objs.len(),
            cpu_objs[0],
            cpu_objs[cpu_objs.len() - 1]
        ),
        mono_cpu,
    );
    g.gate(
        "AD-b",
        format!(
            "descent made real progress: {accepted} accepted steps (>= 10), net CPU-eval decrease {net:.3e} > 0, GPU objective also monotone {mono_gpu}"
        ),
        accepted >= 10 && net > 0.0 && mono_gpu,
    );

    // ---- summary -----------------------------------------------------------
    println!("\ngate summary:");
    for (id, _desc, pass) in &g.rows {
        println!("  {}  {}", if *pass { "PASS" } else { "FAIL" }, id);
    }
    let wall = t_all.elapsed().as_secs_f64();
    println!("total wall time {wall:.1}s");
    let gate_rows: Vec<String> = g
        .rows
        .iter()
        .map(|(id, desc, pass)| {
            format!("{{\"id\": {}, \"pass\": {pass}, \"desc\": {}}}", jstr(id), jstr(desc))
        })
        .collect();
    js.put("gates", format!("[\n    {}\n  ]", gate_rows.join(",\n    ")));
    js.put("all_pass", format!("{}", g.all_pass()));
    js.put("total_wall_s", format!("{wall:.1}"));
    js.put(
        "walltime_notice",
        jstr("llvmpipe is software Vulkan on a 4-core container: every GPU timing here is a correctness/readiness record, NOT a GPU performance claim."),
    );
    if !g.all_pass() {
        std::process::exit(1);
    }
    println!("ALL GATES PASS");
}
