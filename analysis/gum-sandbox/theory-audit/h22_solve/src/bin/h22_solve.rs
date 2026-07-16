//! Phase H2.2 — NR-D2 P4 two-texture reconciliation test on the campaign's
//! own certified engines (fs-gum-field measurement path, fs-gum-statics
//! guard/ANF protocol; seeding pattern from fs-gum-twoknot).
//!
//! THE CONSTRUCTION (spec reconstruction, printed in full — this is the
//! mapping risk):  NR-D2 §2 posits two overlapping double-twist textures
//! with mismatched local anisotropy axes; the engine has no cholesteric
//! sector, so the corpus geometry is interpreted through the hedgehog-class
//! machinery exactly as the twoknot campaign did for App I.2:
//!
//!   * texture patch      <-> B = 1 hedgehog knot (size R* <-> pitch p)
//!   * mismatched axes    <-> knot 2 isorotated by theta = pi about the
//!                            separation axis x-hat ("repulse" channel).
//!                            By hedgehog equivariance q(Rx) = a q(x) a*
//!                            this internal-axis mismatch IS a relative
//!                            rotation of the second texture's frame — the
//!                            most SDiff-flattering channel (the one whose
//!                            reconciliation looks most like a spatial
//!                            rearrangement).
//!   * reconciled state   <-> identical geometry, identity relative
//!                            rotation ("align" channel)
//!   * overlap            <-> separation d = 2 m h with d < R* .. 2R*
//!   * soft manifold      <-> the compositional SDiff orbit (F-R4's
//!                            reading: Q -> Q o T, det DT = 1)
//!   * eps-lifting        <-> the t-weighted kinetic sectors: at the frozen
//!                            dial eps = t(E2+E4)/(E6+E0) = 0.05, the ONLY
//!                            sectors that vary on the SDiff orbit
//!
//! Grid: cubic Field3 n = 96, half = 4.5 (the certified LBOX referee
//! geometry, frozen h = 0.09375).  Frozen statics guards (MU = 5000,
//! band 0.005; MU_F = 400, seed gap - 0.01) with COMMON references taken
//! across both channels' seeds so the wall can never manufacture the
//! mismatch cost (deviation from the per-seed twoknot convention,
//! documented; endpoint penalties are recorded so an active wall is
//! visible).
//!
//! Modes:
//!   h22_solve single <maxit> <outdir>              reference knot, ANF
//!   h22_solve pair   <channel> <m> <maxit> <outdir> full guarded ANF
//!   h22_solve sdiff  <channel> <m> <maxit> <outdir> SDiff-projected descent
//!   h22_solve twist  <channel> <m> <nalpha> <outdir> exact twist-family scan
//!
//! SDIFF-PROJECTED DESCENT (mode `sdiff`) — the documented choice: a
//! projection of the descent step onto divergence-free displacement flows,
//! applied COMPOSITIONALLY with analytic reseeding (the F-R4 family-A
//! pullback pattern: the configuration is always q_seed(X(x)), never a
//! field resampled from itself, so no cumulative interpolation diffusion
//! can fake a value-sector change):
//!   1. w_c(x) = sum_a g_a(x) d_c q_a(x)   (g = the engine's exact guarded
//!      corner-scheme gradient; d_c q by the certified 4th-order stencil);
//!      minus-w is the steepest advective descent direction, since under
//!      q -> q(x - tau u) the objective changes by -tau <u, w>.
//!   2. u = w - grad phi, lap phi = div w   (Leray projection with the
//!      CONSISTENT wide Laplacian div_central o grad_central so that
//!      div_central u = 0 to CG tolerance; Dirichlet phi = 0, plain CG,
//!      warm-started, sequential deterministic reductions; residuals and
//!      the det(grad X) distortion of the composed map are recorded).
//!   3. departure-map composition X <- X(x - tau u/|u|_max) (trilinear on
//!      the smooth displacement D = X - x, D = 0 outside), analytic reseed
//!      q = q_seed(X), renormalize, arrest on objective increase OR on a
//!      frozen-quantity budget violation (|E6+E0 - seed| <= 0.05,
//!      |deg - seed| <= 0.01: exact SDiff keeps both constant, so the
//!      budget enforces the constraint without biasing the objective;
//!      accept -> tau *= 1.2 capped at 0.02, reject -> tau *= 0.5,
//!      stall 1e-5).
//! Volume preservation is exact for the continuum flow and approximate on
//! the lattice; the drift of the SDiff-frozen quantities (E6 + E0, degree)
//! is measured and reported per iteration — that drift IS the numerical
//! error bar of the frozen-value-lemma check.
//!
//! EXACT TWIST FAMILY (mode `twist`): T_alpha(x, y, z) = (x, R(alpha
//! chi(x)) (y,z)) with chi a smoothstep ramp 0 -> 1 across [-s, +s] —
//! det DT = 1 EXACTLY (planar rotation with x-dependent angle).  This is
//! the corpus's "smooth interpolating reorientation" executed literally as
//! an SDiff element: at alpha = pi the map locally undoes the pi mismatch
//! at knot 2 (the rotation axis passes through both centres).  Energies of
//! q_seed o T_alpha are measured over the alpha scan; E6 + E0 must stay
//! flat (quadrature-level) while E2/E4 pay the twist.
//!
//! Determinism: no RNG anywhere; analytic seeds; plain sequential f64
//! loops in fixed order; CG is sequential with a fixed tolerance/cap.
//! Same binary + args => bit-identical JSON.
//!
//! Epistemic notice (binding): within-model numerical engineering on a
//! speculative theory's functional.  Measured facts only; adjudication is
//! the coordinator's.

use fs_gum_field::radial::{interp, radial_solve, RadialProfile};
use fs_gum_field::{bps_floor, qconj, qmul, rstar, stencil, Field3, Scheme, T_FROZEN};
use fs_gum_statics::diag::{aidx, restore_cells, save_cells};
use fs_gum_statics::{anf, eval, AnfParams, Opts, Out, FLOOR_BAND};
use fs_math::det;

use std::fmt::Write as _;
use std::time::Instant;

// ---- frozen protocol constants ---------------------------------------------
/// Cubic grid cells per axis (the certified transverse referee count).
const N: usize = 96;
/// Box half-width (the certified LBOX = 4.5 referee geometry -> h = 0.09375).
const HALF: f64 = 4.5;
/// Radial-profile protocol (the fs-gum-statics gate values).
const RADIAL_N: usize = 4000;
const RADIAL_RMAX: f64 = 6.0;

// ---- channels (fs-gum-twoknot seed.rs, verbatim semantics) -----------------
#[derive(Clone, Copy, PartialEq, Eq)]
enum Channel {
    /// theta = pi about z-hat (perpendicular to separation).
    Attract,
    /// theta = pi about x-hat (parallel to separation): the mismatch channel.
    Repulse,
    /// Identity: the reconciled/aligned reference.
    Align,
}

impl Channel {
    fn parse(s: &str) -> Option<Self> {
        match s {
            "attract" => Some(Channel::Attract),
            "repulse" => Some(Channel::Repulse),
            "align" => Some(Channel::Align),
            _ => None,
        }
    }
    fn as_str(&self) -> &'static str {
        match self {
            Channel::Attract => "attract",
            Channel::Repulse => "repulse",
            Channel::Align => "align",
        }
    }
    /// Isorotation unit quaternion a (q2' = a q2 conj(a)).
    fn quat(&self) -> [f64; 4] {
        match self {
            Channel::Attract => [0.0, 0.0, 0.0, 1.0],
            Channel::Repulse => [0.0, 1.0, 0.0, 0.0],
            Channel::Align => [1.0, 0.0, 0.0, 0.0],
        }
    }
}

// ---- seeding (twoknot product ansatz on the cubic Field3) ------------------
/// Hedgehog quaternion at (x, y, z) relative to the knot centre
/// (twoknot seed.rs `hedgehog_q`, verbatim).
#[inline]
fn hedgehog_q(prof: &RadialProfile, x: f64, y: f64, z: f64) -> [f64; 4] {
    let r = (x * x + y * y + z * z).sqrt();
    let f = interp(r, &prof.r, &prof.f);
    let q0 = det::cos(f);
    let sfr = if r > 1.0e-300 { det::sin(f) / r } else { 0.0 };
    [q0, sfr * x, sfr * y, sfr * z]
}

/// Two-knot product-ansatz value at an ARBITRARY point (needed for the
/// pullback modes): q = q1(x + s) * (a q2(x - s) a*), centres at -+ s x-hat.
#[inline]
fn seed_q(prof: &RadialProfile, s: f64, a: [f64; 4], p: [f64; 3]) -> [f64; 4] {
    let q1 = hedgehog_q(prof, p[0] + s, p[1], p[2]);
    let q2 = hedgehog_q(prof, p[0] - s, p[1], p[2]);
    let ac = qconj(a);
    qmul(q1, qmul(qmul(a, q2), ac))
}

/// Seed the pair on the grid (twoknot `seed_two_knot`, cubic Field3).
fn seed_pair(f: &mut Field3, prof: &RadialProfile, s: f64, ch: Channel) {
    let a = ch.quat();
    for i in 0..N {
        let x = f.x(i);
        for j in 0..N {
            let y = f.x(j);
            for k in 0..N {
                let z = f.x(k);
                f.set(i, j, k, seed_q(prof, s, a, [x, y, z]));
            }
        }
    }
}

/// Seed a single centred knot (the reference run, same grid class).
fn seed_single(f: &mut Field3, prof: &RadialProfile) {
    for i in 0..N {
        let x = f.x(i);
        for j in 0..N {
            let y = f.x(j);
            for k in 0..N {
                let z = f.x(k);
                f.set(i, j, k, hedgehog_q(prof, x, y, z));
            }
        }
    }
}

/// Effective separation diagnostic (twoknot `effective_separation`).
fn effective_separation(f: &Field3) -> f64 {
    let (mut wp, mut xp, mut wm, mut xm) = (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
    for i in 0..N {
        let x = f.x(i);
        for j in 0..N {
            for k in 0..N {
                let q = f.get(i, j, k);
                let w = 1.0 - q[0];
                if x > 0.0 {
                    wp += w;
                    xp += w * x;
                } else {
                    wm += w;
                    xm += w * x;
                }
            }
        }
    }
    xp / wp - xm / wm
}

/// L2 distance (h^3-weighted) between the stored field and the align-channel
/// seed at the same separation: a rough "how reconciled is this state"
/// diagnostic.
fn dist_to_align_seed(f: &Field3, prof: &RadialProfile, s: f64) -> f64 {
    let h3 = f.h() * f.h() * f.h();
    let a = Channel::Align.quat();
    let mut ss = 0.0_f64;
    for i in 0..N {
        let x = f.x(i);
        for j in 0..N {
            let y = f.x(j);
            for k in 0..N {
                let z = f.x(k);
                let qa = seed_q(prof, s, a, [x, y, z]);
                let q = f.get(i, j, k);
                for c in 0..4 {
                    let d = q[c] - qa[c];
                    ss += d * d;
                }
            }
        }
    }
    (ss * h3).sqrt()
}

// ---- guards (common references across channels; see module docs) ----------
struct GuardRefs {
    deg_ref: f64,
    fgap_ref: f64,
    /// Per-channel seed floor gaps in the statics convention (logged).
    fgap_align_seed: f64,
    fgap_chan_seed: f64,
}

/// Compute COMMON guard references from the align seed and the requested
/// channel's seed (statics floor-gap convention: fgap = E6 + E0 -
/// bps_floor * deg / deg_ref).  deg_ref = min of the two seed degrees
/// (the one-sided degree anchor, twoknot protocol).
///
/// WALL DEVIATION (documented): the twoknot campaign referenced the
/// Bogomolny wall to the SEED's own floor gap - 0.01.  The strongly
/// overlapping mismatch seeds here carry ~5 units of E6 + E0 excess above
/// the continuum B = 2 bound, and a seed-referenced wall freezes the value
/// sectors even in the FULL relaxation — i.e. it would manufacture exactly
/// the frozen-value result this test is supposed to measure.  The wall is
/// therefore referenced to the CONTINUUM per-unit-degree Bogomolny bound:
/// it fires only when E6 + E0 < 2 * (32 sqrt2 / 15) - band at the
/// anchor-pinned degree (fgap_ref = bps_floor - FLOOR_BAND in the statics
/// floor-gap convention, since floor_gap = E6 + E0 - bps_floor at
/// deg = deg_ref).  This keeps the wall's one purpose — blocking the
/// near-BPS lattice collapse below the continuum bound — while leaving all
/// legitimate value-sector relaxation unguarded.  Endpoint penalties are
/// recorded; an active wall at the endpoint is a flagged caveat.
fn common_guards(prof: &RadialProfile, s: f64, ch: Channel) -> GuardRefs {
    let mut f = Field3::new(N, HALF);
    seed_pair(&mut f, prof, s, Channel::Align);
    f.renormalize();
    let (oa, _) = eval(&f, &Opts::estatic(Scheme::Corner), false);
    seed_pair(&mut f, prof, s, ch);
    f.renormalize();
    let (oc, _) = eval(&f, &Opts::estatic(Scheme::Corner), false);
    let deg_ref = oa.deg.min(oc.deg);
    let ga = oa.e6 + oa.e0 - bps_floor() * oa.deg / deg_ref;
    let gc = oc.e6 + oc.e0 - bps_floor() * oc.deg / deg_ref;
    GuardRefs {
        deg_ref,
        // wall at E6 + E0 = 2 * bps_floor (continuum B = 2 bound), i.e.
        // floor_gap = 2*floor - floor*(deg/deg_ref ~ 1) = bps_floor
        fgap_ref: bps_floor() - FLOOR_BAND,
        fgap_align_seed: ga,
        fgap_chan_seed: gc,
    }
}

// ---- JSON helpers (twoknot_run pattern) ------------------------------------
fn je(v: f64) -> String {
    if v.is_finite() {
        format!("{v:e}")
    } else {
        "null".to_string()
    }
}

fn out_json(o: &Out) -> String {
    format!(
        "{{\"e2\": {}, \"e4\": {}, \"e6\": {}, \"e0\": {}, \"i\": {}, \"deg\": {}, \"estat\": {}, \"floor_gap\": {}, \"epen\": {}, \"efpen\": {}, \"obj\": {}}}",
        je(o.e2), je(o.e4), je(o.e6), je(o.e0), je(o.i), je(o.deg), je(o.estat),
        je(o.floor_gap), je(o.epen), je(o.efpen), je(o.obj)
    )
}

fn proto_json(mode: &str, chan: &str, m: usize, s: f64, extra: &str) -> String {
    format!(
        "{{\"mode\": \"{mode}\", \"channel\": \"{chan}\", \"m\": {m}, \"d\": {}, \"x_sep\": {}, \"n\": {N}, \"half\": {}, \"h\": {}, \"t\": {}, \"rstar\": {}, \"scheme\": \"corner\", \"radial_n\": {RADIAL_N}, \"radial_rmax\": {}{extra}}}",
        je(2.0 * s), je(2.0 * s / rstar()), je(HALF), je(2.0 * HALF / N as f64),
        je(T_FROZEN), je(rstar()), je(RADIAL_RMAX)
    )
}

fn write_json(outdir: &str, fname: &str, body: String) {
    std::fs::create_dir_all(outdir).expect("outdir");
    let path = format!("{outdir}/{fname}");
    std::fs::write(&path, body).expect("write json");
    println!("  wrote {path}");
}

// ---- mode: single -----------------------------------------------------------
fn run_single(maxit: usize, outdir: &str) {
    let t_all = Instant::now();
    let rp = radial_solve(T_FROZEN, RADIAL_N, RADIAL_RMAX);
    let mut f = Field3::new(N, HALF);
    seed_single(&mut f, &rp);
    let drift = f.renormalize();
    let (seed_c, _) = eval(&f, &Opts::estatic(Scheme::Corner), false);
    let (seed_c4, _) = eval(&f, &Opts::estatic(Scheme::Central4), false);
    let seed_tail = f.boundary_tail();
    let deg_ref = seed_c.deg;
    let fgap_ref = seed_c.e6 + seed_c.e0 - bps_floor() - FLOOR_BAND;
    let opts = Opts::estatic(Scheme::Corner).with_guards(deg_ref, fgap_ref);
    println!(
        "single: seed Estat={:.9} deg={:.6} fgap_ref={:+.6} tail={:.2e}",
        seed_c.estat, seed_c.deg, fgap_ref, seed_tail
    );
    let mut prm = AnfParams::new(maxit);
    prm.instr_every = 10;
    prm.print_every = 50;
    let t0 = Instant::now();
    let res = anf(&mut f, &opts, &prm, "single");
    let secs = t0.elapsed().as_secs_f64();
    let (fin_c4, _) = eval(&f, &Opts::estatic(Scheme::Central4), false);
    let body = format!(
        "{{\n\"protocol\": {},\n\"guards\": {{\"deg_ref\": {}, \"fgap_ref\": {}, \"mu_deg\": {}, \"mu_floor\": {}}},\n\"seed\": {{\"corner\": {}, \"central4\": {}, \"boundary_tail\": {}, \"norm_drift\": {}}},\n\"final\": {{\"corner\": {}, \"central4\": {}, \"boundary_tail\": {}}},\n\"anf\": {{\"status\": \"{}\", \"iters\": {}, \"arrests\": {}, \"gn0\": {}, \"gn\": {}, \"seconds\": {}}}\n}}\n",
        proto_json("single", "single", 0, 0.0, &format!(", \"maxit\": {maxit}")),
        je(deg_ref), je(fgap_ref), je(fs_gum_statics::MU_DEG), je(fs_gum_statics::MU_FLOOR),
        out_json(&seed_c), out_json(&seed_c4), je(seed_tail), je(drift),
        out_json(&res.out), out_json(&fin_c4), je(f.boundary_tail()),
        res.status.as_str(), res.iters, res.arrests, je(res.gn0), je(res.gn), je(secs)
    );
    write_json(outdir, "h22_single.json", body);
    println!("  total {:.1}s", t_all.elapsed().as_secs_f64());
}

// ---- mode: pair (full guarded ANF) ------------------------------------------
fn run_pair(ch: Channel, m: usize, maxit: usize, outdir: &str) {
    let t_all = Instant::now();
    let h = 2.0 * HALF / N as f64;
    let s = m as f64 * h;
    let rp = radial_solve(T_FROZEN, RADIAL_N, RADIAL_RMAX);
    let gr = common_guards(&rp, s, ch);
    let mut f = Field3::new(N, HALF);
    seed_pair(&mut f, &rp, s, ch);
    let drift = f.renormalize();
    let (seed_c, _) = eval(&f, &Opts::estatic(Scheme::Corner), false);
    let (seed_c4, _) = eval(&f, &Opts::estatic(Scheme::Central4), false);
    let seed_dist = dist_to_align_seed(&f, &rp, s);
    let seed_deff = effective_separation(&f);
    let seed_tail = f.boundary_tail();
    let opts = Opts::estatic(Scheme::Corner).with_guards(gr.deg_ref, gr.fgap_ref);
    println!(
        "pair {}: m={m} d={:.6} x={:.4}  seed Estat={:.9} deg={:.6} deg_ref={:.6} fgap_ref={:+.6}",
        ch.as_str(), 2.0 * s, 2.0 * s / rstar(), seed_c.estat, seed_c.deg, gr.deg_ref, gr.fgap_ref
    );
    let mut prm = AnfParams::new(maxit);
    prm.instr_every = 10;
    prm.print_every = 50;
    let t0 = Instant::now();
    let res = anf(&mut f, &opts, &prm, ch.as_str());
    let secs = t0.elapsed().as_secs_f64();
    let (fin_c4, _) = eval(&f, &Opts::estatic(Scheme::Central4), false);
    let fin_dist = dist_to_align_seed(&f, &rp, s);
    let mut sj = String::from("[");
    for (i, r) in res.series.iter().enumerate() {
        if i > 0 {
            sj.push_str(", ");
        }
        let _ = write!(
            sj,
            "{{\"it\": {}, \"obj\": {}, \"estat\": {}, \"e2\": {}, \"e4\": {}, \"e6\": {}, \"e0\": {}, \"deg\": {}, \"floor_gap\": {}, \"epen\": {}, \"efpen\": {}, \"gnorm\": {}, \"dt\": {}, \"arrests\": {}}}",
            r.it, je(r.obj), je(r.estat), je(r.e2), je(r.e4), je(r.e6), je(r.e0), je(r.deg),
            je(r.floor_gap), je(r.epen), je(r.efpen), je(r.gnorm), je(r.dt), r.arrests
        );
    }
    sj.push(']');
    let body = format!(
        "{{\n\"protocol\": {},\n\"guards\": {{\"deg_ref\": {}, \"fgap_ref\": {}, \"fgap_align_seed\": {}, \"fgap_chan_seed\": {}, \"mu_deg\": {}, \"mu_floor\": {}, \"convention\": \"deg_ref = common min across align+channel seeds; wall at the continuum per-unit-degree B=2 Bogomolny bound (see common_guards docs)\"}},\n\"seed\": {{\"corner\": {}, \"central4\": {}, \"d_eff\": {}, \"dist_to_align_seed\": {}, \"boundary_tail\": {}, \"norm_drift\": {}}},\n\"final\": {{\"corner\": {}, \"central4\": {}, \"d_eff\": {}, \"dist_to_align_seed\": {}, \"boundary_tail\": {}}},\n\"anf\": {{\"status\": \"{}\", \"iters\": {}, \"arrests\": {}, \"gn0\": {}, \"gn\": {}, \"seconds\": {}}},\n\"series\": {sj}\n}}\n",
        proto_json("pair", ch.as_str(), m, s, &format!(", \"maxit\": {maxit}")),
        je(gr.deg_ref), je(gr.fgap_ref), je(gr.fgap_align_seed), je(gr.fgap_chan_seed),
        je(fs_gum_statics::MU_DEG), je(fs_gum_statics::MU_FLOOR),
        out_json(&seed_c), out_json(&seed_c4), je(seed_deff),
        je(seed_dist), je(seed_tail), je(drift),
        out_json(&res.out), out_json(&fin_c4), je(effective_separation(&f)), je(fin_dist),
        je(f.boundary_tail()),
        res.status.as_str(), res.iters, res.arrests, je(res.gn0), je(res.gn), je(secs)
    );
    write_json(outdir, &format!("h22_pair_{}_m{m:02}.json", ch.as_str()), body);
    println!("  total {:.1}s", t_all.elapsed().as_secs_f64());
}

// ---- SDiff-projected descent machinery --------------------------------------
#[inline]
fn c3(i: usize, j: usize, k: usize) -> usize {
    (i * N + j) * N + k
}

/// CG solve of -lap_wide phi = b, where lap_wide = div_central o
/// grad_central (the CONSISTENT operator: neighbours at +-2 per axis,
/// spacing 2h, Dirichlet 0 outside), warm start.  Consistency matters:
/// with this operator the projected u = w - grad_central phi satisfies
/// div_central u = 0 to CG tolerance — with the compact 7-point Laplacian
/// it does not, and the descent exploits the leftover compressive modes
/// (observed as an E6+E0 leak in the first sdiff attempt).
/// Returns (iters, final |r|/|b|).
fn cg_poisson(b: &[f64], phi: &mut [f64], h: f64, tol: f64, cap: usize) -> (usize, f64) {
    let n3 = N * N * N;
    let ih2 = 1.0 / (4.0 * h * h);
    let apply = |x: &[f64], out: &mut [f64]| {
        for i in 0..N {
            for j in 0..N {
                for k in 0..N {
                    let c = x[c3(i, j, k)];
                    let mut s = 6.0 * c;
                    if i >= 2 {
                        s -= x[c3(i - 2, j, k)];
                    }
                    if i + 2 < N {
                        s -= x[c3(i + 2, j, k)];
                    }
                    if j >= 2 {
                        s -= x[c3(i, j - 2, k)];
                    }
                    if j + 2 < N {
                        s -= x[c3(i, j + 2, k)];
                    }
                    if k >= 2 {
                        s -= x[c3(i, j, k - 2)];
                    }
                    if k + 2 < N {
                        s -= x[c3(i, j, k + 2)];
                    }
                    out[c3(i, j, k)] = s * ih2;
                }
            }
        }
    };
    let dot = |a: &[f64], b: &[f64]| -> f64 {
        let mut s = 0.0_f64;
        for (x, y) in a.iter().zip(b.iter()) {
            s += x * y;
        }
        s
    };
    let mut r = vec![0.0_f64; n3];
    let mut ap = vec![0.0_f64; n3];
    apply(phi, &mut ap);
    for i in 0..n3 {
        r[i] = b[i] - ap[i];
    }
    let bnorm = dot(b, b).sqrt().max(1.0e-300);
    let mut p = r.clone();
    let mut rr = dot(&r, &r);
    let mut it = 0usize;
    while it < cap && rr.sqrt() / bnorm > tol {
        apply(&p, &mut ap);
        let alpha = rr / dot(&p, &ap).max(1.0e-300);
        for i in 0..n3 {
            phi[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        let rr_new = dot(&r, &r);
        let beta = rr_new / rr;
        rr = rr_new;
        for i in 0..n3 {
            p[i] = r[i] + beta * p[i];
        }
        it += 1;
    }
    (it, rr.sqrt() / bnorm)
}

/// Trilinear interpolation of the displacement field D (3 * N^3, component-
/// major c3 layout) at physical point p; D = 0 outside the real cells.
fn interp_disp(d: &[f64], h: f64, p: [f64; 3]) -> [f64; 3] {
    let mut out = [0.0_f64; 3];
    // grid coordinates: cell centre i at -HALF + (i + 0.5) h
    let gx = (p[0] + HALF) / h - 0.5;
    let gy = (p[1] + HALF) / h - 0.5;
    let gz = (p[2] + HALF) / h - 0.5;
    let (i0, j0, k0) = (gx.floor() as i64, gy.floor() as i64, gz.floor() as i64);
    let (fx, fy, fz) = (gx - i0 as f64, gy - j0 as f64, gz - k0 as f64);
    let get = |c: usize, i: i64, j: i64, k: i64| -> f64 {
        if i < 0 || j < 0 || k < 0 || i >= N as i64 || j >= N as i64 || k >= N as i64 {
            0.0
        } else {
            d[c * N * N * N + c3(i as usize, j as usize, k as usize)]
        }
    };
    for (c, o) in out.iter_mut().enumerate() {
        let mut acc = 0.0_f64;
        for (di, wi) in [(0i64, 1.0 - fx), (1, fx)] {
            for (dj, wj) in [(0i64, 1.0 - fy), (1, fy)] {
                for (dk, wk) in [(0i64, 1.0 - fz), (1, fz)] {
                    acc += wi * wj * wk * get(c, i0 + di, j0 + dj, k0 + dk);
                }
            }
        }
        *o = acc;
    }
    out
}

/// Resample the field from the analytic seed through the departure map
/// X(x) = x + D(x) (family-A pullback pattern: always seed o map, never
/// field o map).  Returns max |norm - 1| before renormalisation.
fn resample_from_map(
    f: &mut Field3,
    prof: &RadialProfile,
    s: f64,
    a: [f64; 4],
    d: &[f64],
) -> f64 {
    let n3 = N * N * N;
    for i in 0..N {
        let x = f.x(i);
        for j in 0..N {
            let y = f.x(j);
            for k in 0..N {
                let z = f.x(k);
                let c = c3(i, j, k);
                let px = [x + d[c], y + d[n3 + c], z + d[2 * n3 + c]];
                f.set(i, j, k, seed_q(prof, s, a, px));
            }
        }
    }
    f.renormalize()
}

/// Volume-distortion diagnostic of the departure map X = x + D:
/// J = I + grad D (2nd-order central, D = 0 outside real cells);
/// returns (max |det J - 1|, mean |det J - 1|).
fn detj_stats(d: &[f64], h: f64) -> (f64, f64) {
    let n3 = N * N * N;
    let i2h = 1.0 / (2.0 * h);
    let get = |c: usize, i: i64, j: i64, k: i64| -> f64 {
        if i < 0 || j < 0 || k < 0 || i >= N as i64 || j >= N as i64 || k >= N as i64 {
            0.0
        } else {
            d[c * n3 + c3(i as usize, j as usize, k as usize)]
        }
    };
    let (mut mx, mut mean) = (0.0_f64, 0.0_f64);
    for i in 0..N as i64 {
        for j in 0..N as i64 {
            for k in 0..N as i64 {
                let mut jm = [[0.0_f64; 3]; 3];
                for (c, row) in jm.iter_mut().enumerate() {
                    row[0] = (get(c, i + 1, j, k) - get(c, i - 1, j, k)) * i2h;
                    row[1] = (get(c, i, j + 1, k) - get(c, i, j - 1, k)) * i2h;
                    row[2] = (get(c, i, j, k + 1) - get(c, i, j, k - 1)) * i2h;
                    row[c] += 1.0;
                }
                let det = jm[0][0] * (jm[1][1] * jm[2][2] - jm[1][2] * jm[2][1])
                    - jm[0][1] * (jm[1][0] * jm[2][2] - jm[1][2] * jm[2][0])
                    + jm[0][2] * (jm[1][0] * jm[2][1] - jm[1][1] * jm[2][0]);
                let e = (det - 1.0).abs();
                mean += e;
                if e > mx {
                    mx = e;
                }
            }
        }
    }
    (mx, mean / n3 as f64)
}

#[allow(clippy::too_many_lines)]
fn run_sdiff(ch: Channel, m: usize, maxit: usize, outdir: &str) {
    let t_all = Instant::now();
    let h = 2.0 * HALF / N as f64;
    let h3 = h * h * h;
    let s = m as f64 * h;
    let a = ch.quat();
    let n3 = N * N * N;
    let rp = radial_solve(T_FROZEN, RADIAL_N, RADIAL_RMAX);
    let gr = common_guards(&rp, s, ch);
    let opts = Opts::estatic(Scheme::Corner).with_guards(gr.deg_ref, gr.fgap_ref);

    // seed = identity map
    let mut dmap = vec![0.0_f64; 3 * n3];
    let mut f = Field3::new(N, HALF);
    seed_pair(&mut f, &rp, s, ch);
    f.renormalize();
    let (seed_c, _) = eval(&f, &Opts::estatic(Scheme::Corner), false);
    let (seed_c4, _) = eval(&f, &Opts::estatic(Scheme::Central4), false);
    let seed_value = seed_c.e6 + seed_c.e0;
    println!(
        "sdiff {}: m={m} d={:.6}  seed Estat={:.9} E6+E0={:.9} deg={:.6}",
        ch.as_str(), 2.0 * s, seed_c.estat, seed_value, seed_c.deg
    );

    // Frozen-quantity budget: exact SDiff keeps E6+E0 and deg constant, so
    // the numerical flow is confined to that surface within these budgets
    // (constraint enforcement, not objective bias; hitting the budget is an
    // arrest, and a stall against it is itself a reported result).
    const VALUE_BUDGET: f64 = 0.05;
    const DEG_BUDGET: f64 = 0.01;
    let (mut out, mut gopt) = eval(&f, &opts, true);
    let mut obj = out.obj;
    let mut tau = 0.01_f64;
    let mut arrests = 0usize;
    let mut budget_hits = 0usize;
    let mut phi = vec![0.0_f64; n3];
    let mut qb = vec![0.0_f64; 4 * n3];
    save_cells(&f, &mut qb);
    let mut series = String::from("[");
    let mut nrec = 0usize;
    let mut status = "iter_cap".to_string();
    let mut accepted = 0usize;
    let mut last_cg = (0usize, 0.0_f64);
    let mut last_wnorm = 0.0_f64;
    let mut last_unorm = 0.0_f64;
    let mut it = 0usize;
    let t0 = Instant::now();
    while it < maxit {
        it += 1;
        let g = gopt.as_ref().expect("grad");
        // w_c = sum_a g_a d_c q_a  (advective steepest direction is -w)
        let mut w = vec![0.0_f64; 3 * n3];
        for i in 0..N {
            for j in 0..N {
                for k in 0..N {
                    let d4 = stencil::derivs4(&f, i, j, k);
                    let c = c3(i, j, k);
                    for ax in 0..3 {
                        let mut acc = 0.0_f64;
                        for cc in 0..4 {
                            acc += g.data[aidx(N, cc, i, j, k)] * d4[ax][cc];
                        }
                        w[ax * n3 + c] = acc;
                    }
                }
            }
        }
        // div w (2nd-order central, w = 0 outside real cells)
        let mut b = vec![0.0_f64; n3];
        let i2h = 1.0 / (2.0 * h);
        let wat = |ax: usize, i: i64, j: i64, k: i64| -> f64 {
            if i < 0 || j < 0 || k < 0 || i >= N as i64 || j >= N as i64 || k >= N as i64 {
                0.0
            } else {
                w[ax * n3 + c3(i as usize, j as usize, k as usize)]
            }
        };
        for i in 0..N as i64 {
            for j in 0..N as i64 {
                for k in 0..N as i64 {
                    // -lap phi = -div w  (so that u = w - grad phi is div-free)
                    b[c3(i as usize, j as usize, k as usize)] = -((wat(0, i + 1, j, k)
                        - wat(0, i - 1, j, k))
                        + (wat(1, i, j + 1, k) - wat(1, i, j - 1, k))
                        + (wat(2, i, j, k + 1) - wat(2, i, j, k - 1)))
                        * i2h;
                }
            }
        }
        last_cg = cg_poisson(&b, &mut phi, h, 1.0e-8, 800);
        // u = w - grad phi
        let mut umax = 0.0_f64;
        let mut wn = 0.0_f64;
        let mut un = 0.0_f64;
        let mut u = w; // reuse buffer, overwrite in place
        let pat = |i: i64, j: i64, k: i64| -> f64 {
            if i < 0 || j < 0 || k < 0 || i >= N as i64 || j >= N as i64 || k >= N as i64 {
                0.0
            } else {
                phi[c3(i as usize, j as usize, k as usize)]
            }
        };
        for i in 0..N as i64 {
            for j in 0..N as i64 {
                for k in 0..N as i64 {
                    let c = c3(i as usize, j as usize, k as usize);
                    let gp = [
                        (pat(i + 1, j, k) - pat(i - 1, j, k)) * i2h,
                        (pat(i, j + 1, k) - pat(i, j - 1, k)) * i2h,
                        (pat(i, j, k + 1) - pat(i, j, k - 1)) * i2h,
                    ];
                    for ax in 0..3 {
                        let wv = u[ax * n3 + c];
                        wn += wv * wv;
                        let uv = wv - gp[ax];
                        u[ax * n3 + c] = uv;
                        un += uv * uv;
                        let av = uv.abs();
                        if av > umax {
                            umax = av;
                        }
                    }
                }
            }
        }
        last_wnorm = wn.sqrt();
        last_unorm = un.sqrt();
        if umax <= 0.0 {
            status = "u_zero".to_string();
            break;
        }
        // trial: D_new(x) = -tau u(x)/umax + D_old(x - tau u(x)/umax)
        let inv = 1.0 / umax;
        let mut dnew = vec![0.0_f64; 3 * n3];
        for i in 0..N {
            let x = -HALF + (i as f64 + 0.5) * h;
            for j in 0..N {
                let y = -HALF + (j as f64 + 0.5) * h;
                for k in 0..N {
                    let z = -HALF + (k as f64 + 0.5) * h;
                    let c = c3(i, j, k);
                    let step = [
                        tau * u[c] * inv,
                        tau * u[n3 + c] * inv,
                        tau * u[2 * n3 + c] * inv,
                    ];
                    let pq = [x - step[0], y - step[1], z - step[2]];
                    let dold = interp_disp(&dmap, h, pq);
                    dnew[c] = -step[0] + dold[0];
                    dnew[n3 + c] = -step[1] + dold[1];
                    dnew[2 * n3 + c] = -step[2] + dold[2];
                }
            }
        }
        let _dr = resample_from_map(&mut f, &rp, s, a, &dnew);
        let (out_new, gnew) = eval(&f, &opts, true);
        let within_budget = (out_new.e6 + out_new.e0 - seed_value).abs() <= VALUE_BUDGET
            && (out_new.deg - seed_c.deg).abs() <= DEG_BUDGET;
        if out_new.obj < obj && within_budget {
            out = out_new;
            obj = out.obj;
            gopt = gnew;
            dmap = dnew;
            save_cells(&f, &mut qb);
            accepted += 1;
            tau = (tau * 1.2).min(0.02);
        } else {
            restore_cells(&mut f, &qb);
            tau *= 0.5;
            arrests += 1;
            if !within_budget {
                budget_hits += 1;
            }
            if tau < 1.0e-5 {
                status = "stall_tau".to_string();
                // record final row below, then stop
            }
        }
        if nrec > 0 {
            series.push_str(", ");
        }
        let (djmax, djmean) = detj_stats(&dmap, h);
        let _ = write!(
            series,
            "{{\"it\": {it}, \"obj\": {}, \"estat\": {}, \"e2\": {}, \"e4\": {}, \"e6\": {}, \"e0\": {}, \"deg\": {}, \"value_drift\": {}, \"deg_drift\": {}, \"epen\": {}, \"efpen\": {}, \"tau\": {}, \"arrests\": {arrests}, \"budget_hits\": {budget_hits}, \"detj_max\": {}, \"detj_mean\": {}, \"cg_iters\": {}, \"cg_res\": {}, \"wnorm\": {}, \"unorm\": {}}}",
            je(out.obj), je(out.estat), je(out.e2), je(out.e4), je(out.e6), je(out.e0),
            je(out.deg), je(out.e6 + out.e0 - seed_value), je(out.deg - seed_c.deg),
            je(out.epen), je(out.efpen), je(tau), je(djmax), je(djmean),
            last_cg.0, je(last_cg.1), je(last_wnorm), je(last_unorm)
        );
        nrec += 1;
        if it % 10 == 0 || it == 1 {
            println!(
                "  [sdiff] it={it} Estat={:.9} E6+E0 drift={:+.3e} deg drift={:+.3e} tau={:.2e} acc={accepted} arr={arrests} bud={budget_hits} detJ=({:.1e}, {:.1e}) cg=({}, {:.1e})",
                out.estat, out.e6 + out.e0 - seed_value, out.deg - seed_c.deg, tau,
                djmax, djmean, last_cg.0, last_cg.1
            );
        }
        if status == "stall_tau" {
            break;
        }
    }
    series.push(']');
    let secs = t0.elapsed().as_secs_f64();
    // final measurements on the last ACCEPTED configuration
    restore_cells(&mut f, &qb);
    let (fin_c, _) = eval(&f, &Opts::estatic(Scheme::Corner), false);
    let (fin_c4, _) = eval(&f, &Opts::estatic(Scheme::Central4), false);
    let fin_dist = dist_to_align_seed(&f, &rp, s);
    // map diagnostics: max |D|, L2 |D|
    let (mut dmax, mut dl2) = (0.0_f64, 0.0_f64);
    for c in 0..n3 {
        let dd = [dmap[c], dmap[n3 + c], dmap[2 * n3 + c]];
        let m2 = dd[0] * dd[0] + dd[1] * dd[1] + dd[2] * dd[2];
        dl2 += m2;
        if m2.sqrt() > dmax {
            dmax = m2.sqrt();
        }
    }
    dl2 = (dl2 * h3).sqrt();
    println!(
        "  [sdiff] DONE ({status}) it={it} acc={accepted} Estat={:.9} (seed {:.9})  E6+E0 drift={:+.3e}  deg drift={:+.3e}  ({:.1}s)",
        fin_c.estat, seed_c.estat, fin_c.e6 + fin_c.e0 - seed_value, fin_c.deg - seed_c.deg, secs
    );
    let (djmax_fin, djmean_fin) = detj_stats(&dmap, h);
    let body = format!(
        "{{\n\"protocol\": {},\n\"guards\": {{\"deg_ref\": {}, \"fgap_ref\": {}, \"mu_deg\": {}, \"mu_floor\": {}, \"convention\": \"deg_ref = common min across align+channel seeds; wall at the continuum per-unit-degree B=2 Bogomolny bound (see common_guards docs)\"}},\n\"seed\": {{\"corner\": {}, \"central4\": {}}},\n\"final\": {{\"corner\": {}, \"central4\": {}, \"dist_to_align_seed\": {}, \"boundary_tail\": {}, \"map_dmax\": {}, \"map_dl2\": {}, \"detj_max\": {}, \"detj_mean\": {}}},\n\"descent\": {{\"status\": \"{status}\", \"iters\": {it}, \"accepted\": {accepted}, \"arrests\": {arrests}, \"budget_hits\": {budget_hits}, \"seconds\": {}, \"value_drift_final\": {}, \"deg_drift_final\": {}}},\n\"series\": {series}\n}}\n",
        proto_json("sdiff", ch.as_str(), m, s,
            &format!(", \"maxit\": {maxit}, \"cg_tol\": 1e-8, \"cg_cap\": 800, \"tau0\": 0.01, \"tau_max\": 0.02, \"value_budget\": 0.05, \"deg_budget\": 0.01, \"laplacian\": \"wide div_central o grad_central (consistent projector)\"")),
        je(gr.deg_ref), je(gr.fgap_ref), je(fs_gum_statics::MU_DEG), je(fs_gum_statics::MU_FLOOR),
        out_json(&seed_c), out_json(&seed_c4),
        out_json(&fin_c), out_json(&fin_c4), je(fin_dist), je(f.boundary_tail()),
        je(dmax), je(dl2), je(djmax_fin), je(djmean_fin),
        je(secs), je(fin_c.e6 + fin_c.e0 - seed_value), je(fin_c.deg - seed_c.deg)
    );
    write_json(outdir, &format!("h22_sdiff_{}_m{m:02}.json", ch.as_str()), body);
    println!("  total {:.1}s", t_all.elapsed().as_secs_f64());
}

// ---- mode: twist (exact volume-preserving interpolating reorientation) -----
fn run_twist(ch: Channel, m: usize, nalpha: usize, outdir: &str) {
    let t_all = Instant::now();
    let h = 2.0 * HALF / N as f64;
    let s = m as f64 * h;
    let a = ch.quat();
    let rp = radial_solve(T_FROZEN, RADIAL_N, RADIAL_RMAX);
    let mut f = Field3::new(N, HALF);
    // smoothstep ramp chi: 0 at x <= -s, 1 at x >= +s
    let chi = |x: f64| -> f64 {
        if x <= -s {
            0.0
        } else if x >= s {
            1.0
        } else {
            let u = (x + s) / (2.0 * s);
            u * u * (3.0 - 2.0 * u)
        }
    };
    let mut rows = String::from("[");
    for ia in 0..nalpha {
        let alpha = core::f64::consts::PI * ia as f64 / (nalpha - 1) as f64;
        for i in 0..N {
            let x = f.x(i);
            let th = alpha * chi(x);
            let (ct, st) = (det::cos(th), det::sin(th));
            for j in 0..N {
                let y = f.x(j);
                for k in 0..N {
                    let z = f.x(k);
                    // T_alpha: rotate (y, z) about the x-axis by theta(x)
                    let yp = ct * y - st * z;
                    let zp = st * y + ct * z;
                    f.set(i, j, k, seed_q(&rp, s, a, [x, yp, zp]));
                }
            }
        }
        f.renormalize();
        let (oc, _) = eval(&f, &Opts::estatic(Scheme::Corner), false);
        let dist = dist_to_align_seed(&f, &rp, s);
        if ia > 0 {
            rows.push_str(", ");
        }
        let _ = write!(
            rows,
            "{{\"alpha\": {}, \"corner\": {}, \"dist_to_align_seed\": {}}}",
            je(alpha),
            out_json(&oc),
            je(dist)
        );
        println!(
            "  twist alpha={alpha:.4}  Estat={:.9} E2={:.6} E4={:.6} E6={:.9} E0={:.9} deg={:.6} dist={:.4}",
            oc.estat, oc.e2, oc.e4, oc.e6, oc.e0, oc.deg, dist
        );
    }
    rows.push(']');
    let body = format!(
        "{{\n\"protocol\": {},\n\"scan\": {rows}\n}}\n",
        proto_json("twist", ch.as_str(), m, s,
            &format!(", \"nalpha\": {nalpha}, \"ramp\": \"smoothstep on [-s, s]\", \"map\": \"x-axis twist, det DT = 1 exact\""))
    );
    write_json(outdir, &format!("h22_twist_{}_m{m:02}.json", ch.as_str()), body);
    println!("  total {:.1}s", t_all.elapsed().as_secs_f64());
}

// ---- mode: seedscan (channel-vs-align seed energies across separations) ----
fn run_seedscan(m_lo: usize, m_hi: usize, outdir: &str) {
    let h = 2.0 * HALF / N as f64;
    let rp = radial_solve(T_FROZEN, RADIAL_N, RADIAL_RMAX);
    let mut f = Field3::new(N, HALF);
    let mut rows = String::from("[");
    for m in m_lo..=m_hi {
        let s = m as f64 * h;
        seed_pair(&mut f, &rp, s, Channel::Align);
        f.renormalize();
        let (oa, _) = eval(&f, &Opts::estatic(Scheme::Corner), false);
        seed_pair(&mut f, &rp, s, Channel::Repulse);
        f.renormalize();
        let (oc, _) = eval(&f, &Opts::estatic(Scheme::Corner), false);
        if m > m_lo {
            rows.push_str(", ");
        }
        let _ = write!(
            rows,
            "{{\"m\": {m}, \"d\": {}, \"x_sep\": {}, \"align\": {}, \"repulse\": {}}}",
            je(2.0 * s),
            je(2.0 * s / rstar()),
            out_json(&oa),
            out_json(&oc)
        );
        println!(
            "  m={m:2} d={:.4} x={:.3}  E(align)={:.6} E(repulse)={:.6}  dE(rep-al)={:+.6}  dValue={:+.6}",
            2.0 * s,
            2.0 * s / rstar(),
            oa.estat,
            oc.estat,
            oc.estat - oa.estat,
            (oc.e6 + oc.e0) - (oa.e6 + oa.e0)
        );
    }
    rows.push(']');
    let body = format!(
        "{{\n\"protocol\": {},\n\"scan\": {rows}\n}}\n",
        proto_json("seedscan", "align+repulse", m_lo, m_lo as f64 * h,
            &format!(", \"m_lo\": {m_lo}, \"m_hi\": {m_hi}"))
    );
    write_json(outdir, &format!("h22_seedscan_m{m_lo:02}_m{m_hi:02}.json"), body);
}

// ---- main -------------------------------------------------------------------
fn main() {
    let args: Vec<String> = std::env::args().collect();
    let usage = "usage: h22_solve single <maxit> <outdir> | pair <channel> <m> <maxit> <outdir> | sdiff <channel> <m> <maxit> <outdir> | twist <channel> <m> <nalpha> <outdir>";
    if args.len() < 2 {
        eprintln!("{usage}");
        std::process::exit(2);
    }
    match args[1].as_str() {
        "single" if args.len() == 4 => {
            run_single(args[2].parse().expect("maxit"), &args[3]);
        }
        "seedscan" if args.len() == 5 => {
            run_seedscan(args[2].parse().expect("m_lo"), args[3].parse().expect("m_hi"), &args[4]);
        }
        "pair" | "sdiff" | "twist" if args.len() == 6 => {
            let ch = Channel::parse(&args[2]).expect("channel");
            let m: usize = args[3].parse().expect("m");
            let k: usize = args[4].parse().expect("maxit/nalpha");
            match args[1].as_str() {
                "pair" => run_pair(ch, m, k, &args[5]),
                "sdiff" => run_sdiff(ch, m, k, &args[5]),
                _ => run_twist(ch, m, k, &args[5]),
            }
        }
        _ => {
            eprintln!("{usage}");
            std::process::exit(2);
        }
    }
}
