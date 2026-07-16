//! Phase H3.2 — T-H5 adjudication computation (ROADMAP_v5 Phase H3).
//!
//! THE CLAIM UNDER TEST (corpus IV.D, T2 §d.3 DEFECT-CANDIDATE): "the rigid
//! spinning knot [kappa_paper = sqrt(7/6) = 1.080] is superradiant and is not
//! a solution.  Self-consistency forces centrifugal deformation until
//! kappa < 1: the in-gap character of matter is a derived necessity."  The
//! REJECTION of the rigid rung is SOUND (T2 §d.2, convention-robust); what is
//! asserted-not-derived is the CONVERGENCE story — that the deformation flow
//! stops once in-gap, at/near the corpus's benchmark (kappa = 0.802 +- 0.018
//! at eps = 0.05).  T3.1 predicts instead that the honest flow passes
//! kappa = 1 without any distinguished behaviour and lands on the halo
//! saturation locus kappa_paper = 1/sqrt(2) (the constrained optimum IS the
//! saturation point, KKT-marginal; h26 already saw the descent-relaxed clock
//! ROOT drift to 0.6997 ~ 1/sqrt(2) with budget).
//!
//! THE COMPUTATION: the corpus's narrative, executed literally and honestly —
//! a quasi-static selection flow alternating
//!   (a) M_DESC guarded ANF descent steps of the fixed-L Routhian
//!       R(q; L) = E_static(q) + L^2/(2 I(q))  (energy-decreasing deformation
//!       at fixed charge/spin: the ANF objective is monotone non-increasing
//!       BY CONSTRUCTION within each block), with
//!   (b) the corpus's own clock self-consistency L^2 = (2/3) I E_static
//!       (Section-K form of E_tot = c*omega at j = 1/2), re-anchoring L to
//!       the deformed state,
//! starting FROM the rigid-rung state: the relaxed spherical static solution
//! (the engine's rigid/undeformed knot at N = 48) spun at the rigid-rung
//! clock rate, L_rigid = sqrt(7/6) * I_static / (2 sqrt(pi)) via the
//! validated unit map kappa_paper = 2 sqrt(pi) L / I.
//!
//! Instrumented per macro-cycle: kappa_paper (block and post-clock-update),
//! E_rot/E, halo fraction, degree, the bare Routhian at block start/end
//! (=> the honest descent rate dR/diter at fixed L), guard penalties,
//! gradient ratio, arrests.  The kappa = 1 crossing and the endpoint trend
//! are computed downstream (h32_analyze.py).
//!
//! CONTROLS (the F-R5-class rescue test): the identical flow with the halo
//! channels suppressed by a crude in-gap projection — after every descent
//! block, every cell farther than r_c from the current (1 - q0) centroid is
//! RESET to the static-anchor field (vacuum-class out there) and the field
//! renormalised, before the clock update reads (I, E_static).  Two radii:
//!   (i)  r_c = 1.5 R*  — the campaign's halo diagnostic radius (far-field
//!        suppression as literally stated by F-R5's halo referee);
//!   (ii) r_c = 1.25 R* — a compact-support constraint just above the
//!        corpus's own in-family support (dilation V = sqrt(2) has support
//!        radius V^(1/3) R* = 1.12 R*; the benchmark V = 1.409 likewise),
//!        i.e. the tightest reading under which the corpus's claimed
//!        deformations (dilation + oblate g-dial) remain admissible.
//! These are deliberately crude: they are the minimal unstated constraints
//! that could make the corpus's "deform until kappa < 1 then stop"
//! narrative come out — if a constrained flow stabilises near the corpus's
//! benchmark, the narrative is RESCUABLE by a constraint the corpus never
//! states (same class as the F-R5 lock); if not, not even that.  Per-cycle
//! I-radial diagnostics (I-fraction beyond 1.0/1.25/1.5 R*, I-rms radius)
//! locate where the inertia growth actually lives.
//!
//! Protocol (the campaign's frozen N = 48 gate protocol, h26 anchor stage
//! verbatim): LBOX = 4.5, eps = 0.05 (t = T_FROZEN), corner scheme, guards
//! MU = 5000 at the engine hedgehog degree / MU_F = 400 at hedgehog gap
//! - 0.01; static stage = hedgehog + bump(424242, K = 6, amp 0.02), 900
//! guarded ANF iterations.  Unit maps (SESSION_HANDOFF §3): kappa_ours =
//! L/I, kappa_paper = 2 sqrt(pi) kappa_ours, c_paper = 2L/sqrt(2 pi^3);
//! halo threshold kappa_paper = 1/sqrt(2); corpus onset 1.000 +- 0.004;
//! corpus benchmark kappa = 0.802 +- 0.018.
//!
//! Usage: h32_selection <outdir> [m_desc] [ncyc] [threads]
//! (defaults 25 400 4; total descent budget = m_desc * ncyc per run).
//! Same binary + args => bit-identical JSON (fs-gum-kern bit contract; the
//! first three free-flow cycles are re-run and asserted bitwise identical).
//!
//! Epistemic notice (binding): within-model numerical engineering on a
//! speculative theory's functional.  Measured facts about the campaign's
//! own certified functional and the corpus's argument structure; nothing
//! here says anything about nature.  Adjudication is the coordinator's.

use fs_gum_field::radial::radial_solve;
use fs_gum_field::{Field3, Scheme, T_FROZEN};
use fs_gum_statics::diag::{add_scaled, aidx, bump_field, halo_fraction, restore_cells, save_cells};
use fs_gum_statics::{halo_radius, kappa_threshold, AnfParams, Opts, FLOOR_BAND};

use std::fmt::Write as _;
use std::time::Instant;

const LBOX: f64 = 4.5;
const N: usize = 48;
const MAXIT_STATIC: usize = 900;
const SEED_PERT: u64 = 424242;
const PI: f64 = std::f64::consts::PI;

/// One macro-cycle of the quasi-static selection flow.
#[derive(Clone, Copy, Debug)]
struct Cyc {
    cyc: usize,
    /// Spin held fixed during this block (set by the previous clock update).
    l: f64,
    /// Measured on the state that feeds the clock update (post-projection
    /// in the constrained runs).
    i: f64,
    /// I before the in-gap projection (== i when free): i_preproj - i is
    /// the inertia the projection cuts per cycle.
    i_preproj: f64,
    estat: f64,
    /// Bare Routhian at block start (this L) and block end — dR/diter at
    /// fixed L is (r_end - r0)/iters, guaranteed <= 0 within a block.
    r0: f64,
    r_end: f64,
    /// Bare Routhian before the in-gap projection (== r_end when free).
    r_preproj: f64,
    obj: f64,
    kappa_ours: f64,
    kappa_paper: f64,
    erot_frac: f64,
    /// Halo fraction after descent, before any projection.
    halo_pre: f64,
    /// Halo fraction after projection (== halo_pre when free).
    halo_post: f64,
    /// I-radial diagnostics on the clock-feeding state: I-fraction beyond
    /// 1.0 / 1.25 / 1.5 R* of the centroid, and the I-rms radius (in R*).
    ifrac_r100: f64,
    ifrac_r125: f64,
    ifrac_r150: f64,
    r_rms_over_rstar: f64,
    deg: f64,
    floor_gap: f64,
    epen: f64,
    efpen: f64,
    gn_ratio: f64,
    arrests: usize,
    iters: usize,
    status: &'static str,
    /// Cells reset by the in-gap projection this cycle (0 when free).
    ncells_proj: usize,
    /// Clock update output: L_next^2 = (2/3) I E_static.
    l_next: f64,
    kappa_paper_next: f64,
    c_paper_next: f64,
}

fn je(v: f64) -> String {
    if v.is_finite() {
        format!("{v:.17e}")
    } else {
        format!("\"{v}\"")
    }
}

fn cyc_json(c: &Cyc) -> String {
    format!(
        "{{\"cyc\": {}, \"L\": {}, \"I\": {}, \"I_preproj\": {}, \"Estat\": {}, \"R0\": {}, \"R_end\": {}, \"R_preproj\": {}, \"obj\": {}, \"kappa_ours\": {}, \"kappa_paper\": {}, \"Erot_over_Etot\": {}, \"halo_pre\": {}, \"halo_post\": {}, \"ifrac_r100\": {}, \"ifrac_r125\": {}, \"ifrac_r150\": {}, \"r_rms_over_rstar\": {}, \"deg\": {}, \"floor_gap\": {}, \"epen\": {}, \"efpen\": {}, \"gn_ratio\": {}, \"arrests\": {}, \"iters\": {}, \"status\": \"{}\", \"ncells_proj\": {}, \"L_next\": {}, \"kappa_paper_next\": {}, \"c_paper_next\": {}}}",
        c.cyc, je(c.l), je(c.i), je(c.i_preproj), je(c.estat), je(c.r0), je(c.r_end),
        je(c.r_preproj), je(c.obj), je(c.kappa_ours), je(c.kappa_paper), je(c.erot_frac),
        je(c.halo_pre), je(c.halo_post), je(c.ifrac_r100), je(c.ifrac_r125), je(c.ifrac_r150),
        je(c.r_rms_over_rstar), je(c.deg), je(c.floor_gap), je(c.epen), je(c.efpen),
        je(c.gn_ratio), c.arrests, c.iters, c.status, c.ncells_proj, je(c.l_next),
        je(c.kappa_paper_next), je(c.c_paper_next)
    )
}

/// I-radial diagnostics about `center`: (I-fraction beyond 1.0 R*, beyond
/// 1.25 R*, beyond 1.5 R*, I-rms radius / R*).  Same I-density convention
/// as `halo_fraction` (di = q1^2 + q2^2 per cell; fractions are h^3-free).
fn i_radial_diag(f: &Field3, center: [f64; 3]) -> (f64, f64, f64, f64) {
    let n = f.n();
    let rstar = fs_gum_field::rstar();
    let (r1, r2, r3) = (rstar, 1.25 * rstar, 1.5 * rstar);
    let (r1sq, r2sq, r3sq) = (r1 * r1, r2 * r2, r3 * r3);
    let (mut tot, mut o1, mut o2, mut o3, mut mr2) = (0.0_f64, 0.0, 0.0, 0.0, 0.0);
    for i in 0..n {
        let dx = f.x(i) - center[0];
        for j in 0..n {
            let dy = f.x(j) - center[1];
            for k in 0..n {
                let dz = f.x(k) - center[2];
                let q = f.get(i, j, k);
                let di = q[1] * q[1] + q[2] * q[2];
                let r2c = dx * dx + dy * dy + dz * dz;
                tot += di;
                mr2 += di * r2c;
                if r2c > r1sq {
                    o1 += di;
                }
                if r2c > r2sq {
                    o2 += di;
                }
                if r2c > r3sq {
                    o3 += di;
                }
            }
        }
    }
    (o1 / tot, o2 / tot, o3 / tot, (mr2 / tot).sqrt() / rstar)
}

/// Crude in-gap projection: reset every real cell farther than `rc` from
/// `center` to the static-anchor snapshot values, then renormalise.
/// Returns the number of cells reset.
fn project_far(f: &mut Field3, anchor: &[f64], center: [f64; 3], rc: f64) -> usize {
    let n = f.n();
    let rc2 = rc * rc;
    let mut cnt = 0usize;
    for i in 0..n {
        let dx = f.x(i) - center[0];
        for j in 0..n {
            let dy = f.x(j) - center[1];
            for k in 0..n {
                let dz = f.x(k) - center[2];
                if dx * dx + dy * dy + dz * dz > rc2 {
                    let q = [
                        anchor[aidx(n, 0, i, j, k)],
                        anchor[aidx(n, 1, i, j, k)],
                        anchor[aidx(n, 2, i, j, k)],
                        anchor[aidx(n, 3, i, j, k)],
                    ];
                    f.set(i, j, k, q);
                    cnt += 1;
                }
            }
        }
    }
    let _ = f.renormalize();
    cnt
}

struct FlowCtx<'a> {
    q_static: &'a [f64],
    deg_ref: f64,
    fgap_ref: f64,
    threads: usize,
}

/// The quasi-static selection flow.  Restores the static anchor, spins it
/// at l0, then alternates m_desc guarded Routhian descent steps at fixed L
/// with the clock update L <- sqrt((2/3) I E_static).  `proj_rc = Some(rc)`
/// applies the crude in-gap projection at radius rc after every descent
/// block.  Returns the seed measurement and the per-cycle trajectory; the
/// field is left at the endpoint.
fn run_flow(
    ctx: &FlowCtx,
    f: &mut Field3,
    l0: f64,
    m_desc: usize,
    ncyc: usize,
    proj_rc: Option<f64>,
    label: &str,
    verbose: bool,
) -> (Cyc, Vec<Cyc>) {
    let s2pi = 2.0 * PI.sqrt();
    let c_paper_of = |l: f64| 2.0 * l / (2.0 * PI * PI * PI).sqrt();
    restore_cells(f, ctx.q_static);
    let mut l = l0;

    // seed measurement (cycle 0: the rigid-rung state, no descent yet)
    let opts0 = Opts::routhian(Scheme::Corner, l).with_guards(ctx.deg_ref, ctx.fgap_ref);
    let (out0, _) = fs_gum_kern::eval(f, &opts0, false, ctx.threads);
    let (halo0, cen0) = halo_fraction(f);
    let (if0a, if0b, if0c, rr0) = i_radial_diag(f, cen0);
    let erot0 = l * l / (2.0 * out0.i);
    let seed = Cyc {
        cyc: 0,
        l,
        i: out0.i,
        i_preproj: out0.i,
        estat: out0.estat,
        r0: out0.r,
        r_end: out0.r,
        r_preproj: out0.r,
        obj: out0.obj,
        kappa_ours: l / out0.i,
        kappa_paper: s2pi * l / out0.i,
        erot_frac: erot0 / (out0.estat + erot0),
        halo_pre: halo0,
        halo_post: halo0,
        ifrac_r100: if0a,
        ifrac_r125: if0b,
        ifrac_r150: if0c,
        r_rms_over_rstar: rr0,
        deg: out0.deg,
        floor_gap: out0.floor_gap,
        epen: out0.epen,
        efpen: out0.efpen,
        gn_ratio: 1.0,
        arrests: 0,
        iters: 0,
        status: "seed",
        ncells_proj: 0,
        l_next: l,
        kappa_paper_next: s2pi * l / out0.i,
        c_paper_next: c_paper_of(l),
    };
    if verbose {
        println!(
            "[{label}] seed: L={:.6} I={:.5} Estat={:.6} kappa_paper={:.6} (rigid rung sqrt(7/6)={:.6}) halo={:.4}",
            l,
            out0.i,
            out0.estat,
            seed.kappa_paper,
            (7.0_f64 / 6.0).sqrt(),
            halo0
        );
    }

    let mut cycles = Vec::with_capacity(ncyc);
    for cyc in 1..=ncyc {
        // (a) bounded guarded descent at fixed L (energy-decreasing
        //     deformation at fixed charge/spin)
        let opts = Opts::routhian(Scheme::Corner, l).with_guards(ctx.deg_ref, ctx.fgap_ref);
        let prm = AnfParams::new(m_desc);
        let res = fs_gum_kern::anf(f, &opts, &prm, label, ctx.threads);
        let r0 = res.series.first().map_or(f64::NAN, |rec| rec.r);
        let (halo_pre, cen) = halo_fraction(f);
        let r_preproj = res.out.r;
        let i_preproj = res.out.i;

        // (controls only) the crude in-gap projection
        let (out, halo_post, cen_post, ncells_proj) = if let Some(rc) = proj_rc {
            let ncp = project_far(f, ctx.q_static, cen, rc);
            let (o2, _) = fs_gum_kern::eval(f, &opts, false, ctx.threads);
            let (hp, cp) = halo_fraction(f);
            (o2, hp, cp, ncp)
        } else {
            (res.out, halo_pre, cen, 0)
        };
        let (ifa, ifb, ifc, rr) = i_radial_diag(f, cen_post);

        // (b) the corpus's own clock self-consistency
        let l_next = ((2.0 / 3.0) * out.i * out.estat).sqrt();
        let erot = l * l / (2.0 * out.i);
        let c = Cyc {
            cyc,
            l,
            i: out.i,
            i_preproj,
            estat: out.estat,
            r0,
            r_end: out.r,
            r_preproj,
            obj: out.obj,
            kappa_ours: l / out.i,
            kappa_paper: s2pi * l / out.i,
            erot_frac: erot / (out.estat + erot),
            halo_pre,
            halo_post,
            ifrac_r100: ifa,
            ifrac_r125: ifb,
            ifrac_r150: ifc,
            r_rms_over_rstar: rr,
            deg: out.deg,
            floor_gap: out.floor_gap,
            epen: out.epen,
            efpen: out.efpen,
            gn_ratio: res.gn / res.gn0,
            arrests: res.arrests,
            iters: res.iters,
            status: res.status.as_str(),
            ncells_proj,
            l_next,
            kappa_paper_next: s2pi * l_next / out.i,
            c_paper_next: c_paper_of(l_next),
        };
        if verbose && (cyc % 10 == 0 || cyc == 1) {
            println!(
                "[{label}] cyc={cyc:3} L={:.5} -> I={:.4} Estat={:.5} kappa_paper={:.5} (next {:.5}) halo={:.4} if(1,1.25,1.5R*)=({:.3},{:.3},{:.4}) r_rms={:.3}R*{} dR/dit={:+.2e} deg={:.5} arr={}",
                c.l,
                c.i,
                c.estat,
                c.kappa_paper,
                c.kappa_paper_next,
                c.halo_pre,
                c.ifrac_r100,
                c.ifrac_r125,
                c.ifrac_r150,
                c.r_rms_over_rstar,
                if proj_rc.is_some() {
                    format!(" dI_cut={:+.3e}", c.i_preproj - c.i)
                } else {
                    String::new()
                },
                (c.r_end - c.r0) / (c.iters.max(1) as f64),
                c.deg,
                c.arrests
            );
        }
        cycles.push(c);
        l = l_next;
    }
    (seed, cycles)
}

fn main() {
    let t0 = Instant::now();
    let args: Vec<String> = std::env::args().collect();
    let outdir = args.get(1).cloned().unwrap_or_else(|| ".".to_string());
    let m_desc: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(25);
    let ncyc: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(400);
    let threads: usize = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(4);

    println!(
        "h32_selection: N={N} LBOX={LBOX} t={T_FROZEN:.15} m_desc={m_desc} ncyc={ncyc} (budget {} iters/run) threads={threads}",
        m_desc * ncyc
    );

    // ---- the rigid rung in the corpus's own algebra (record) --------------
    let e0hat = 64.0 / (15.0 * PI);
    let i0hat = 256.0 / (105.0 * PI);
    let kappa_rigid_paper = (2.0 * e0hat / (3.0 * i0hat)).sqrt(); // sqrt(7/6)
    println!(
        "rigid rung (paper algebra): kappa^2 = 2*e0/(3*i0) = 7/6 exactly; kappa = {kappa_rigid_paper:.9}; c = 2*i0*kappa = {:.6}",
        2.0 * i0hat * kappa_rigid_paper
    );

    // ---- radial profile + hedgehog references (gates G-B / h26 protocol) --
    let rp = radial_solve(T_FROZEN, 4000, 6.0);
    let mut f = Field3::new(N, LBOX);
    f.sample_hedgehog(&rp.r, &rp.f, 1.0);
    let (out_h, _) = fs_gum_kern::eval(&f, &Opts::estatic(Scheme::Corner), false, threads);
    let deg_ref = out_h.deg;
    let fgap_h = out_h.e6 + out_h.e0 - fs_gum_field::bps_floor();
    let fgap_ref = fgap_h - FLOOR_BAND;
    println!(
        "hedgehog (corner N={N}): Estat={:.7} deg={:.6} I={:.5} -> deg_ref={deg_ref:.6} fgap_ref={fgap_ref:+.6}",
        out_h.estat, out_h.deg, out_h.i
    );

    // ---- static stage (the engine's rigid/undeformed knot) ----------------
    let bump = bump_field(N, LBOX, SEED_PERT, 6, 0.02);
    add_scaled(&mut f, &bump, 1.0);
    let _ = f.renormalize();
    let opts_static = Opts::estatic(Scheme::Corner).with_guards(deg_ref, fgap_ref);
    let prm_s = AnfParams::new(MAXIT_STATIC);
    let res_s = fs_gum_kern::anf(&mut f, &opts_static, &prm_s, "static", threads);
    let mut q_static = vec![0.0_f64; 4 * N * N * N];
    save_cells(&f, &mut q_static);
    let (i_s, e_s) = (res_s.out.i, res_s.out.estat);
    println!(
        "static anchor: status={} iters={} arrests={} Estat={:.7} I={:.5} deg={:.6} ({:.1}s)",
        res_s.status.as_str(),
        res_s.iters,
        res_s.arrests,
        e_s,
        i_s,
        res_s.out.deg,
        t0.elapsed().as_secs_f64()
    );

    // ---- the rigid-rung spin in engine units (unit map, exact) ------------
    // kappa_paper = 2 sqrt(pi) L / I  =>  L_rigid = sqrt(7/6) I_static / (2 sqrt(pi))
    let l_rigid = kappa_rigid_paper * i_s / (2.0 * PI.sqrt());
    // for the record: the clock value at the frozen anchor (h26's L0)
    let l_clock_static = ((2.0 / 3.0) * i_s * e_s).sqrt();
    println!(
        "L_rigid = sqrt(7/6) I_s/(2 sqrt(pi)) = {l_rigid:.6}   (clock at frozen anchor: L = {l_clock_static:.6}, kappa_paper = {:.6})",
        2.0 * PI.sqrt() * l_clock_static / i_s
    );

    let ctx = FlowCtx { q_static: &q_static, deg_ref, fgap_ref, threads };
    let rstar = fs_gum_field::rstar();
    let rc_halo = halo_radius(); // 1.5 R*
    let rc_supp = 1.25 * rstar;

    // ---- run 1: the free (honest) selection flow ---------------------------
    println!("\n=== FREE FLOW (honest energy-decreasing deformation + clock) ===");
    let (seed_free, traj_free) =
        run_flow(&ctx, &mut f, l_rigid, m_desc, ncyc, None, "free", true);
    let last_free = *traj_free.last().expect("free trajectory nonempty");
    println!(
        "free endpoint (cyc {ncyc}): kappa_paper={:.6} (next {:.6}) vs 1/sqrt(2)={:.6}, benchmark 0.802; halo={:.4}; I/I_s={:.4} ({:.1}s)",
        last_free.kappa_paper,
        last_free.kappa_paper_next,
        1.0 / 2.0_f64.sqrt(),
        last_free.halo_pre,
        last_free.i / i_s,
        t0.elapsed().as_secs_f64()
    );

    // ---- run 2: control at the halo diagnostic radius (1.5 R*) ------------
    println!("\n=== CONTROL 1 (far-field re-seeded at r > 1.5 R* each cycle) ===");
    let (seed_ch, traj_ch) =
        run_flow(&ctx, &mut f, l_rigid, m_desc, ncyc, Some(rc_halo), "ctl_halo", true);
    let last_ch = *traj_ch.last().expect("control trajectory nonempty");
    println!(
        "control-1.5R* endpoint (cyc {ncyc}): kappa_paper={:.6} (next {:.6}) vs benchmark 0.802 +- 0.018; halo_pre={:.4} post={:.4}; I/I_s={:.4} ({:.1}s)",
        last_ch.kappa_paper,
        last_ch.kappa_paper_next,
        last_ch.halo_pre,
        last_ch.halo_post,
        last_ch.i / i_s,
        t0.elapsed().as_secs_f64()
    );

    // ---- run 3: control at the compact-support radius (1.25 R*) -----------
    println!("\n=== CONTROL 2 (compact support: re-seeded at r > 1.25 R* each cycle) ===");
    let (seed_cs, traj_cs) =
        run_flow(&ctx, &mut f, l_rigid, m_desc, ncyc, Some(rc_supp), "ctl_supp", true);
    let last_cs = *traj_cs.last().expect("control trajectory nonempty");
    println!(
        "control-1.25R* endpoint (cyc {ncyc}): kappa_paper={:.6} (next {:.6}) vs benchmark 0.802 +- 0.018, exact-tier 0.847; I/I_s={:.4}; dI_cut={:+.3e} ({:.1}s)",
        last_cs.kappa_paper,
        last_cs.kappa_paper_next,
        last_cs.i / i_s,
        last_cs.i_preproj - last_cs.i,
        t0.elapsed().as_secs_f64()
    );

    // ---- determinism replay: first 3 free cycles, bitwise ------------------
    let (seed_r, traj_r) = run_flow(&ctx, &mut f, l_rigid, m_desc, 3, None, "replay", false);
    let a = traj_free[2];
    let b = traj_r[2];
    let replay_ok = seed_r.i.to_bits() == seed_free.i.to_bits()
        && a.i.to_bits() == b.i.to_bits()
        && a.estat.to_bits() == b.estat.to_bits()
        && a.deg.to_bits() == b.deg.to_bits()
        && a.l_next.to_bits() == b.l_next.to_bits();
    println!("\nreplay of free cycles 1-3: bitwise (I, Estat, deg, L_next) identical = {replay_ok}");
    assert!(replay_ok, "determinism replay failed");

    // ---- JSON ---------------------------------------------------------------
    let traj_json = |traj: &[Cyc]| -> String {
        let mut s = String::new();
        for (i, c) in traj.iter().enumerate() {
            let _ = write!(s, "{}\n  {}", if i == 0 { "" } else { "," }, cyc_json(c));
        }
        s
    };
    let run_json = |seed: &Cyc, traj: &[Cyc], rc: Option<f64>, desc: &str| -> String {
        format!(
            "{{\"projection_radius\": {}, \"projection_desc\": \"{desc}\", \"seed\": {}, \"cycles\": [{}\n]}}",
            rc.map_or("null".to_string(), je),
            cyc_json(seed),
            traj_json(traj)
        )
    };
    let body = format!(
        "{{\n\"phase\": \"H3.2 T-H5 selection-flow adjudication (quasi-static: guarded Routhian descent at fixed L alternated with the clock update L^2 = (2/3) I E_static, from the rigid rung)\",\n\"date\": \"2026-07-16\",\n\"protocol\": {{\"N\": {N}, \"LBOX\": {LBOX}, \"t_frozen\": {}, \"eps\": 0.05, \"scheme\": \"corner\", \"guards\": {{\"deg_ref\": {}, \"fgap_ref\": {}, \"mu_deg\": 5000.0, \"mu_floor\": 400.0}}, \"static\": {{\"seed\": {SEED_PERT}, \"kbumps\": 6, \"amp\": 0.02, \"maxit\": {MAXIT_STATIC}, \"I\": {}, \"Estat\": {}}}, \"m_desc\": {m_desc}, \"ncyc\": {ncyc}, \"descent_budget_per_run\": {}, \"threads\": {threads}, \"kernel\": \"fs-gum-kern tiled deterministic layer (bit-identical at any thread count)\", \"rstar\": {}, \"projection_rule\": \"cells with r > r_c of the current (1 - q0) centroid reset to the static-anchor snapshot, then renormalize; applied after each descent block, before the clock update reads (I, E_static)\"}},\n\"rigid_rung\": {{\"kappa_paper\": {}, \"kappa_paper_sq\": \"7/6 exactly\", \"L_rigid_engine\": {}, \"unit_map\": \"kappa_paper = 2 sqrt(pi) L / I; L_rigid = sqrt(7/6) I_static/(2 sqrt(pi))\", \"L_clock_at_frozen_anchor\": {}, \"c_paper_of_L\": \"2 L / sqrt(2 pi^3)\"}},\n\"references\": {{\"halo_threshold_paper\": {}, \"halo_threshold_ours\": {}, \"corpus_onset\": [1.000, 0.004], \"corpus_benchmark_kappa\": [0.802, 0.018], \"h26_exact_tier_kappa\": 0.847028, \"h26_descent_root_kappa\": 0.699745, \"saturation_locus\": {}, \"corpus_deep_bps_kappa\": {}}},\n\"free\": {},\n\"control_halo_1p5rstar\": {},\n\"control_support_1p25rstar\": {},\n\"determinism\": {{\"replay_bitwise_identical\": {replay_ok}}},\n\"runtime_seconds\": {}\n}}\n",
        je(T_FROZEN), je(deg_ref), je(fgap_ref), je(i_s), je(e_s), m_desc * ncyc, je(rstar),
        je(kappa_rigid_paper), je(l_rigid), je(l_clock_static),
        je(1.0 / 2.0_f64.sqrt()), je(kappa_threshold()),
        je(1.0 / 2.0_f64.sqrt()), je((7.0_f64 / 12.0).sqrt()),
        run_json(&seed_free, &traj_free, None, "none (honest flow)"),
        run_json(
            &seed_ch,
            &traj_ch,
            Some(rc_halo),
            "far-field re-seed at 1.5 R* (the campaign halo diagnostic radius)"
        ),
        run_json(
            &seed_cs,
            &traj_cs,
            Some(rc_supp),
            "compact-support re-seed at 1.25 R* (just above the corpus in-family support V^(1/3) R* ~ 1.12 R*)"
        ),
        je(t0.elapsed().as_secs_f64())
    );
    let path = format!("{outdir}/h32_results.json");
    std::fs::write(&path, &body).expect("write json");
    println!("wrote {path}");
    println!("total {:.1}s", t0.elapsed().as_secs_f64());
}
