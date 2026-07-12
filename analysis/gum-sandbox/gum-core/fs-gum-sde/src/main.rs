//! fs-gum-sde gate + demo binary (Phase B4).
//!
//! Runs, in order:
//! - G1: WienerStream contract gates (shared ΔW bits, per-trajectory
//!   isolation replay, ensemble order-independence, in-process replay).
//! - G2: analytic-SDE gates — OU moment closure at finite t (N = 10⁵),
//!   measured strong convergence orders (EM on additive OU ≈ 1.0, SRA1 ≈
//!   1.5, EM on multiplicative GBM ≈ 0.5) and weak orders via dt-halving
//!   ladders with FIXED Wiener paths (coarse increments aggregated exactly
//!   from the finest level, including the ΔZ composition rule).
//! - G3: Nelson stationary state — ψ = φ₁₁, ensemble born at |φ₁₁|²:
//!   H̄ must sit at the finite-N noise floor (the sharpest factor-of-2
//!   detector in the crate).
//! - G4: the H-theorem demo — N = 20,000 Nelson trajectories, M = 16
//!   (born.py phases, seed 42; sampling streams seeded 12345), ρ₀ = |φ₁₁|²,
//!   coarse-grained H̄(t) at 32×32 and 16×16 to t = 4π, plus the
//!   equilibrium-born control. Referee: Tier-3 qualitative gate
//!   (near-exponential decay to the ~0.026 floor; control flat at floor).
//!
//! Outputs (deterministic — no wall-clock content): `hdata.csv`,
//! `results.json` next to the crate's Cargo.toml. Exit code 1 on any gate
//! failure.

use std::fmt::Write as _;
use std::time::Instant;

use fs_gum_sde::hfunc::{coarse_q, fit_tau, h_bar, histogram_p, noise_floor, TauFit};
use fs_gum_sde::modes::ModeSystem;
use fs_gum_sde::nelson::{nelson_macro_step, sample_phi11, sample_psi2, CoeffCache};
use fs_gum_sde::sde::{em_step_with, sra1_step_with, EulerMaruyama, Scratch, Sde, Sra1, Stepper};
use fs_gum_sde::wiener::WienerStream;
use fs_math::det;

// ---------------------------------------------------------------------------
// Parameters (born.py values where they exist).
// ---------------------------------------------------------------------------
const SEED_PHASES: u64 = 42; // evolution-noise seed (born.py phase seed)
const SEED_SAMPLE: u64 = 12345; // ensemble-sampling seed (born.py)
const N_PART: usize = 20_000;
const DT: f64 = 2e-3;
const DT_OUT: f64 = 0.1;
const NFINE: usize = 256;
const CGS: [usize; 2] = [32, 16];
const T_FINAL: f64 = 4.0 * std::f64::consts::PI;
/// Tier-3 deterministic (de Broglie–Bohm) referee: results.json M=16 cg32.
const BOHM_TAU_CG32: f64 = 4.852765879922653;

// Kernel registry (one kernel id per independent noise consumer).
const K_NELSON_M16: u32 = 1;
const K_NELSON_CTRL: u32 = 2;
const K_NELSON_G3: u32 = 3;
const K_INIT_PHI11: u32 = 10;
const K_INIT_PSI2: u32 = 11;
const K_INIT_G3: u32 = 12;
const K_OU_MOM_EM: u32 = 20;
const K_OU_MOM_SRA1: u32 = 21;
const K_OU_STRONG: u32 = 22;
const K_GBM_STRONG: u32 = 23;
const K_OU_WEAK: u32 = 24;
const K_CUBIC_STRONG: u32 = 25;
const K_G1: u32 = 90;

// ---------------------------------------------------------------------------
// Gate bookkeeping.
// ---------------------------------------------------------------------------
struct Gate {
    id: &'static str,
    desc: &'static str,
    value: String,
    target: String,
    pass: bool,
}

fn push_gate(gates: &mut Vec<Gate>, id: &'static str, desc: &'static str, value: String, target: String, pass: bool) {
    println!("  [{}] {} — value {} vs target {} => {}", id, desc, value, target, if pass { "PASS" } else { "FAIL" });
    gates.push(Gate { id, desc, value, target, pass });
}

// ---------------------------------------------------------------------------
// Analytic SDEs for the gates.
// ---------------------------------------------------------------------------
/// Ornstein–Uhlenbeck dX = −θX dt + σ dW (additive).
struct Ou {
    theta: f64,
    sigma: f64,
}
impl Sde for Ou {
    fn dim(&self) -> usize {
        1
    }
    fn drift(&self, _t: f64, x: &[f64], out: &mut [f64]) {
        out[0] = -self.theta * x[0];
    }
    fn diffusion(&self, _t: f64, _x: &[f64], out: &mut [f64]) {
        out[0] = self.sigma;
    }
}

/// Additive-noise double well dX = (X − X³) dt + σ dW — the NONLINEAR-drift
/// referee for SRA1's generic strong order 1.5. On LINEAR drift (OU) SRA1
/// superconverges to strong order 2.0: the h^1.5-limiting Itô–Taylor
/// residuals are the zero-mean fluctuations of terms carrying f''(x)
/// (e.g. ½g²f'' via (ΔZ/h)² in the H2 stage), which vanish identically for
/// linear f. A nonlinear drift (f'' = −6x ≠ 0) restores the generic rate.
struct CubicWell {
    sigma: f64,
}
impl Sde for CubicWell {
    fn dim(&self) -> usize {
        1
    }
    fn drift(&self, _t: f64, x: &[f64], out: &mut [f64]) {
        out[0] = x[0] - x[0] * x[0] * x[0];
    }
    fn diffusion(&self, _t: f64, _x: &[f64], out: &mut [f64]) {
        out[0] = self.sigma;
    }
}

/// Geometric Brownian motion dX = μX dt + σX dW (multiplicative — the
/// EM-strong-0.5 referee; SRA1 is NOT applied to it).
struct Gbm {
    mu: f64,
    sigma: f64,
}
impl Sde for Gbm {
    fn dim(&self) -> usize {
        1
    }
    fn drift(&self, _t: f64, x: &[f64], out: &mut [f64]) {
        out[0] = self.mu * x[0];
    }
    fn diffusion(&self, _t: f64, x: &[f64], out: &mut [f64]) {
        out[0] = self.sigma * x[0];
    }
}

// ---------------------------------------------------------------------------
// G2 machinery: moments and dt-halving ladders on fixed Wiener paths.
// ---------------------------------------------------------------------------
fn ou_moments(
    stepper: &impl Stepper,
    ou: &Ou,
    kernel: u32,
    npaths: usize,
    nsteps: u64,
    t_end: f64,
    x0: f64,
) -> (f64, f64) {
    let ws = WienerStream::new(SEED_PHASES, kernel, 1);
    let h = t_end / nsteps as f64;
    let mut scratch = Scratch::new(1);
    let (mut sum, mut sumsq) = (0.0f64, 0.0f64);
    for p in 0..npaths {
        let mut x = [x0];
        for s in 0..nsteps {
            stepper.step(ou, &ws, p as u32, s, s as f64 * h, h, &mut x, &mut scratch);
        }
        sum += x[0];
        sumsq += x[0] * x[0];
    }
    let n = npaths as f64;
    let mean = sum / n;
    (mean, sumsq / n - mean * mean)
}

/// Endpoints of an SDE integrated at several dt-halving levels on the SAME
/// Wiener path per trajectory: the finest-level increments (ΔW, ΔZ) are
/// generated once through the WienerStream addressing, and coarser levels
/// aggregate them exactly (ΔW by summation; ΔZ by the composition rule
/// ΔZ_coarse = Σ_j [ΔZ_j + h_f · (W(t_j) − W(t_start))]).
/// Returns (endpoints[level][path], total W(T) per path).
#[allow(clippy::too_many_arguments)]
fn ladder<S: Sde>(
    sde: &S,
    use_sra1: bool,
    kernel: u32,
    npaths: usize,
    k_fine: u32,
    levels: &[u32],
    t_end: f64,
    x0: f64,
) -> (Vec<Vec<f64>>, Vec<f64>) {
    let ws = WienerStream::new(SEED_PHASES, kernel, 1);
    let nf = 1usize << k_fine;
    let hf = t_end / nf as f64;
    let mut endpoints = vec![vec![0.0f64; npaths]; levels.len()];
    let mut w_tot = vec![0.0f64; npaths];
    let mut dwf = vec![0.0f64; nf];
    let mut dzf = vec![0.0f64; nf];
    let mut scratch = Scratch::new(1);
    let (mut dwb, mut dzb) = ([0.0f64], [0.0f64]);
    for p in 0..npaths {
        let mut wsum = 0.0f64;
        for i in 0..nf {
            ws.dw_dz(p as u32, i as u64, hf, &mut dwb, &mut dzb);
            dwf[i] = dwb[0];
            dzf[i] = dzb[0];
            wsum += dwb[0];
        }
        w_tot[p] = wsum;
        for (li, &k) in levels.iter().enumerate() {
            let nsteps = 1usize << k;
            let stride = 1usize << (k_fine - k);
            let h = t_end / nsteps as f64;
            let mut x = [x0];
            for m in 0..nsteps {
                let base = m * stride;
                let (mut dz, mut wacc) = (0.0f64, 0.0f64);
                for i in 0..stride {
                    dz += dzf[base + i] + hf * wacc;
                    wacc += dwf[base + i];
                }
                let dw = wacc;
                let t = m as f64 * h;
                if use_sra1 {
                    sra1_step_with(sde, t, h, &mut x, &[dw], &[dz], &mut scratch);
                } else {
                    em_step_with(sde, t, h, &mut x, &[dw], &mut scratch);
                }
            }
            endpoints[li][p] = x[0];
        }
    }
    (endpoints, w_tot)
}

/// Least-squares slope of ln(e) vs ln(h).
fn loglog_slope(hs: &[f64], es: &[f64]) -> f64 {
    let n = hs.len() as f64;
    let (mut sx, mut sy, mut sxx, mut sxy) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
    for (&h, &e) in hs.iter().zip(es) {
        let (lx, ly) = (det::ln(h), det::ln(e));
        sx += lx;
        sy += ly;
        sxx += lx * lx;
        sxy += lx * ly;
    }
    (n * sxy - sx * sy) / (n * sxx - sx * sx)
}

/// Strong error ladder from successive-level differences (coupled paths):
/// e_k = E|X_k(T) − X_{k+1}(T)| ∝ h_k^p.
fn strong_slope(levels: &[u32], endpoints: &[Vec<f64>], t_end: f64) -> (Vec<f64>, Vec<f64>, f64) {
    let npaths = endpoints[0].len() as f64;
    let mut hs = Vec::new();
    let mut es = Vec::new();
    for li in 0..levels.len() - 1 {
        let mut acc = 0.0f64;
        for p in 0..endpoints[li].len() {
            acc += (endpoints[li][p] - endpoints[li + 1][p]).abs();
        }
        hs.push(t_end / f64::from(1u32 << levels[li]));
        es.push(acc / npaths);
    }
    let slope = loglog_slope(&hs, &es);
    (hs, es, slope)
}

/// Weak error ladder from successive-level mean differences (common random
/// numbers): e_k = |E[X_k(T) − X_{k+1}(T)]| ∝ h_k^p.
fn weak_slope(
    levels: &[u32],
    endpoints: &[Vec<f64>],
    t_end: f64,
    n_use: usize,
) -> (Vec<f64>, Vec<f64>, f64) {
    let npaths = endpoints[0].len() as f64;
    let mut hs = Vec::new();
    let mut es = Vec::new();
    for li in 0..(levels.len() - 1).min(n_use) {
        let mut acc = 0.0f64;
        for p in 0..endpoints[li].len() {
            acc += endpoints[li][p] - endpoints[li + 1][p];
        }
        hs.push(t_end / f64::from(1u32 << levels[li]));
        es.push((acc / npaths).abs());
    }
    let slope = loglog_slope(&hs, &es);
    (hs, es, slope)
}

// ---------------------------------------------------------------------------
// Nelson ensemble driver (born.py run_ensemble structure).
// ---------------------------------------------------------------------------
struct HCurves {
    ts: Vec<f64>,
    h: [Vec<f64>; 2], // per CGS entry
    /// Fine early-time series (every [`FINE_EVERY`] macro steps while
    /// t ≤ fine_until): the Nelson relaxation turns out to be far faster
    /// than born.py's DT_OUT = 0.1 diagnostic cadence can resolve, so the
    /// τ fits use this series.
    fine_ts: Vec<f64>,
    fine_h: [Vec<f64>; 2],
}

/// Fine-sampling cadence: every 5 macro steps = 0.01 time units.
const FINE_EVERY: u64 = 5;

fn record_h(
    modes: &ModeSystem,
    t: f64,
    xs: &[f64],
    ys: &[f64],
    ts: &mut Vec<f64>,
    h: &mut [Vec<f64>; 2],
) {
    let qf = modes.psi2_grid(t, NFINE);
    ts.push(t);
    for (ci, &cg) in CGS.iter().enumerate() {
        let p = histogram_p(xs, ys, cg);
        let q = coarse_q(&qf, NFINE, cg);
        h[ci].push(h_bar(&p, &q));
    }
}

fn run_nelson_ensemble(
    modes: &ModeSystem,
    ws: &WienerStream,
    xs: &mut [f64],
    ys: &mut [f64],
    t_final: f64,
    fine_until: f64,
    label: &str,
) -> HCurves {
    let t0 = Instant::now();
    let n = xs.len();
    let n_full = (t_final / DT_OUT).floor() as usize;
    let steps_per_out = (DT_OUT / DT).round() as usize;
    let mut out = HCurves {
        ts: Vec::new(),
        h: [Vec::new(), Vec::new()],
        fine_ts: Vec::new(),
        fine_h: [Vec::new(), Vec::new()],
    };
    let mut cache = CoeffCache::new();
    let mut step: u64 = 0;
    let mut t = 0.0f64;
    record_h(modes, t, xs, ys, &mut out.ts, &mut out.h);
    if fine_until > 0.0 {
        record_h(modes, t, xs, ys, &mut out.fine_ts, &mut out.fine_h);
    }
    for block in 1..=n_full {
        for _ in 0..steps_per_out {
            cache.clear();
            for i in 0..n {
                nelson_macro_step(modes, ws, &mut cache, i as u32, step, t, DT, &mut xs[i], &mut ys[i]);
            }
            step += 1;
            t += DT;
            if t <= fine_until && step % FINE_EVERY == 0 {
                record_h(modes, t, xs, ys, &mut out.fine_ts, &mut out.fine_h);
            }
        }
        record_h(modes, t, xs, ys, &mut out.ts, &mut out.h);
        if block % 25 == 0 {
            println!(
                "    [{label}] block {block}/{n_full}  t = {t:.3}  H32 = {:.4}  ({:.1} s)",
                out.h[0].last().unwrap(),
                t0.elapsed().as_secs_f64()
            );
        }
    }
    let rem = t_final - t;
    if rem > 1e-12 {
        let n_extra = ((rem / DT).ceil() as usize).max(1);
        let h = rem / n_extra as f64;
        for _ in 0..n_extra {
            cache.clear();
            for i in 0..n {
                nelson_macro_step(modes, ws, &mut cache, i as u32, step, t, h, &mut xs[i], &mut ys[i]);
            }
            step += 1;
            t += h;
        }
        t = t_final;
        record_h(modes, t, xs, ys, &mut out);
    }
    println!(
        "    [{label}] done in {:.1} s (final H32 = {:.4})",
        t0.elapsed().as_secs_f64(),
        out.h[0].last().unwrap()
    );
    out
}

// ---------------------------------------------------------------------------
// FNV-1a 64 over f64 bit patterns (determinism fingerprint).
// ---------------------------------------------------------------------------
fn fnv1a_f64(mut acc: u64, xs: &[f64]) -> u64 {
    for &x in xs {
        for b in x.to_bits().to_le_bytes() {
            acc ^= u64::from(b);
            acc = acc.wrapping_mul(0x0000_0100_0000_01B3);
        }
    }
    acc
}
const FNV_OFFSET: u64 = 0xCBF2_9CE4_8422_2325;

// ---------------------------------------------------------------------------
// Main.
// ---------------------------------------------------------------------------
#[allow(clippy::too_many_lines)]
fn main() {
    let wall = Instant::now();
    let mut gates: Vec<Gate> = Vec::new();
    let outdir = env!("CARGO_MANIFEST_DIR");

    // =======================================================================
    // G1 — WienerStream contract.
    // =======================================================================
    println!("== G1: WienerStream contract ==");
    {
        // G1.a — EM and SRA1 read the SAME ΔW bits at every address.
        let ws = WienerStream::new(SEED_PHASES, K_G1, 2);
        let (mut a, mut w, mut z) = ([0.0f64; 2], [0.0f64; 2], [0.0f64; 2]);
        let mut same = true;
        for (traj, slot) in [(0u32, 0u64), (7, 123), (19_999, 6283 * 32 + 31)] {
            ws.dw(traj, slot, DT, &mut a);
            ws.dw_dz(traj, slot, DT, &mut w, &mut z);
            same &= a[0].to_bits() == w[0].to_bits() && a[1].to_bits() == w[1].to_bits();
        }
        push_gate(&mut gates, "G1.a", "EM/SRA1 share dW bits", format!("{same}"), "true".into(), same);

        // G1.b/G1.c — isolation replay + order permutation on a short
        // Nelson ensemble (M = 16 drift, real substepping).
        let modes = ModeSystem::born_m16();
        let wsn = WienerStream::new(SEED_PHASES, K_G1 + 1, 2);
        let nn = 64usize;
        let steps = 100u64;
        let mut xs: Vec<f64> = (0..nn).map(|i| sample_phi11(SEED_SAMPLE, K_G1 + 2, i as u32).0).collect();
        let mut ys: Vec<f64> = (0..nn).map(|i| sample_phi11(SEED_SAMPLE, K_G1 + 2, i as u32).1).collect();
        let (x0, y0) = (xs.clone(), ys.clone());
        let mut cache = CoeffCache::new();
        for step in 0..steps {
            let t = step as f64 * DT;
            cache.clear();
            for i in 0..nn {
                nelson_macro_step(&modes, &wsn, &mut cache, i as u32, step, t, DT, &mut xs[i], &mut ys[i]);
            }
        }
        // Isolation replay of trajectory 41 (fresh cache, no ensemble).
        let (mut rx, mut ry) = (x0[41], y0[41]);
        let mut fresh = CoeffCache::new();
        for step in 0..steps {
            let t = step as f64 * DT;
            fresh.clear();
            nelson_macro_step(&modes, &wsn, &mut fresh, 41, step, t, DT, &mut rx, &mut ry);
        }
        let iso = rx.to_bits() == xs[41].to_bits() && ry.to_bits() == ys[41].to_bits();
        push_gate(&mut gates, "G1.b", "trajectory isolation replay (bitwise)", format!("{iso}"), "true".into(), iso);

        let mut xr = x0;
        let mut yr = y0;
        let mut cache2 = CoeffCache::new();
        for step in 0..steps {
            let t = step as f64 * DT;
            cache2.clear();
            for i in (0..nn).rev() {
                nelson_macro_step(&modes, &wsn, &mut cache2, i as u32, step, t, DT, &mut xr[i], &mut yr[i]);
            }
        }
        let perm = (0..nn).all(|i| xr[i].to_bits() == xs[i].to_bits() && yr[i].to_bits() == ys[i].to_bits());
        push_gate(&mut gates, "G1.c", "ensemble order-independence (reverse sweep, bitwise)", format!("{perm}"), "true".into(), perm);
    }

    // =======================================================================
    // G2 — analytic-SDE gates.
    // =======================================================================
    println!("== G2: OU moment closure (N = 1e5) ==");
    let ou = Ou {
        theta: 1.0,
        sigma: 0.5,
    };
    let (x0, t_end) = (1.0f64, 1.0f64);
    let mean_exact = x0 * det::exp(-ou.theta * t_end);
    let var_exact = ou.sigma * ou.sigma * (1.0 - det::exp(-2.0 * ou.theta * t_end)) / (2.0 * ou.theta);
    let (em_mean, em_var) = ou_moments(&EulerMaruyama, &ou, K_OU_MOM_EM, 100_000, 256, t_end, x0);
    let (sra_mean, sra_var) = ou_moments(&Sra1, &ou, K_OU_MOM_SRA1, 100_000, 256, t_end, x0);
    // In-process bit-identical replay of the EM moment ensemble.
    let (em_mean2, em_var2) = ou_moments(&EulerMaruyama, &ou, K_OU_MOM_EM, 100_000, 256, t_end, x0);
    let replay = em_mean.to_bits() == em_mean2.to_bits() && em_var.to_bits() == em_var2.to_bits();
    push_gate(&mut gates, "G1.d", "OU ensemble in-process replay (bitwise)", format!("{replay}"), "true".into(), replay);

    let em_me = (em_mean - mean_exact).abs();
    let em_ve = (em_var / var_exact - 1.0).abs();
    let sra_me = (sra_mean - mean_exact).abs();
    let sra_ve = (sra_var / var_exact - 1.0).abs();
    push_gate(&mut gates, "G2.a", "EM OU |mean err| (exact 0.367879)", format!("{em_me:.2e}"), "< 5e-3".into(), em_me < 5e-3);
    push_gate(&mut gates, "G2.b", "EM OU |var rel err| (exact 0.108030)", format!("{em_ve:.2e}"), "< 2e-2".into(), em_ve < 2e-2);
    push_gate(&mut gates, "G2.c", "SRA1 OU |mean err|", format!("{sra_me:.2e}"), "< 5e-3".into(), sra_me < 5e-3);
    push_gate(&mut gates, "G2.d", "SRA1 OU |var rel err|", format!("{sra_ve:.2e}"), "< 2e-2".into(), sra_ve < 2e-2);

    println!("== G2: convergence orders (fixed Wiener paths, dt-halving) ==");
    let strong_levels: Vec<u32> = (4..=10).collect();
    let (ep_em, _) = ladder(&ou, false, K_OU_STRONG, 2000, 12, &strong_levels, t_end, x0);
    let (ep_sra, _) = ladder(&ou, true, K_OU_STRONG, 2000, 12, &strong_levels, t_end, x0);
    let (_, es_em, p_em) = strong_slope(&strong_levels, &ep_em, t_end);
    let (_, es_sra, p_sra) = strong_slope(&strong_levels, &ep_sra, t_end);
    println!("    EM   strong errors: {es_em:?}");
    println!("    SRA1 strong errors: {es_sra:?}");
    push_gate(&mut gates, "G2.e", "strong order EM on OU (additive => 1.0)", format!("{p_em:.3}"), "[0.85, 1.15]".into(), (0.85..=1.15).contains(&p_em));
    push_gate(&mut gates, "G2.f", "strong order SRA1 on OU (LINEAR: superconverges => 2.0)", format!("{p_sra:.3}"), "[1.85, 2.15]".into(), (1.85..=2.15).contains(&p_sra));

    // Nonlinear additive drift: the generic SRA1 strong order 1.5.
    let cubic = CubicWell { sigma: 0.5 };
    let (cep_em, _) = ladder(&cubic, false, K_CUBIC_STRONG, 2000, 12, &strong_levels, t_end, x0);
    let (cep_sra, _) = ladder(&cubic, true, K_CUBIC_STRONG, 2000, 12, &strong_levels, t_end, x0);
    let (_, ces_em, cp_em) = strong_slope(&strong_levels, &cep_em, t_end);
    let (_, ces_sra, cp_sra) = strong_slope(&strong_levels, &cep_sra, t_end);
    println!("    EM   cubic-well strong errors: {ces_em:?}");
    println!("    SRA1 cubic-well strong errors: {ces_sra:?}");
    push_gate(&mut gates, "G2.g", "strong order EM on cubic well (additive => 1.0)", format!("{cp_em:.3}"), "[0.85, 1.15]".into(), (0.85..=1.15).contains(&cp_em));
    push_gate(&mut gates, "G2.h", "strong order SRA1 on cubic well (additive => 1.5)", format!("{cp_sra:.3}"), "[1.35, 1.65]".into(), (1.35..=1.65).contains(&cp_sra));

    let gbm = Gbm { mu: 1.0, sigma: 1.0 };
    let gbm_levels: Vec<u32> = (4..=9).collect();
    let (ep_gbm, w_tot) = ladder(&gbm, false, K_GBM_STRONG, 2000, 9, &gbm_levels, t_end, x0);
    let mut hs_g = Vec::new();
    let mut es_g = Vec::new();
    for (li, &k) in gbm_levels.iter().enumerate() {
        let mut acc = 0.0f64;
        for p in 0..ep_gbm[li].len() {
            let exact = x0 * det::exp((gbm.mu - 0.5 * gbm.sigma * gbm.sigma) * t_end + gbm.sigma * w_tot[p]);
            acc += (ep_gbm[li][p] - exact).abs();
        }
        hs_g.push(t_end / f64::from(1u32 << k));
        es_g.push(acc / ep_gbm[li].len() as f64);
    }
    let p_gbm = loglog_slope(&hs_g, &es_g);
    println!("    EM   GBM strong errors (vs exact path solution): {es_g:?}");
    push_gate(&mut gates, "G2.i", "strong order EM on GBM (multiplicative => 0.5)", format!("{p_gbm:.3}"), "[0.35, 0.65]".into(), (0.35..=0.65).contains(&p_gbm));

    let weak_levels: Vec<u32> = (1..=7).collect();
    let (wep_em, _) = ladder(&ou, false, K_OU_WEAK, 20_000, 7, &weak_levels, t_end, x0);
    let (wep_sra, _) = ladder(&ou, true, K_OU_WEAK, 20_000, 7, &weak_levels, t_end, x0);
    let (_, we_em, q_em) = weak_slope(&weak_levels, &wep_em, t_end, 5);
    let (_, we_sra, q_sra) = weak_slope(&weak_levels, &wep_sra, t_end, 4);
    println!("    EM   weak errors (CRN successive diffs): {we_em:?}");
    println!("    SRA1 weak errors (CRN successive diffs): {we_sra:?}");
    push_gate(&mut gates, "G2.j", "weak order EM on OU (=> 1)", format!("{q_em:.3}"), "[0.8, 1.3]".into(), (0.8..=1.3).contains(&q_em));
    push_gate(&mut gates, "G2.k", "weak order SRA1 on OU (=> 2)", format!("{q_sra:.3}"), "[1.6, 2.5]".into(), (1.6..=2.5).contains(&q_sra));

    // =======================================================================
    // G3 — Nelson stationary state (psi = phi_11): H̄ pinned at the floor.
    // =======================================================================
    println!("== G3: Nelson stationary state phi_11 (N = {N_PART}, T = 2 pi) ==");
    let g3_curves;
    {
        let modes = ModeSystem::phi11();
        let ws = WienerStream::new(SEED_PHASES, K_NELSON_G3, 2);
        let mut xs = vec![0.0f64; N_PART];
        let mut ys = vec![0.0f64; N_PART];
        for i in 0..N_PART {
            let (x, y) = sample_phi11(SEED_SAMPLE, K_INIT_G3, i as u32);
            xs[i] = x;
            ys[i] = y;
        }
        g3_curves = run_nelson_ensemble(&modes, &ws, &mut xs, &mut ys, 2.0 * std::f64::consts::PI, "G3 phi11");
    }
    let floor32 = noise_floor(32, N_PART);
    let g3_max = g3_curves.h[0].iter().cloned().fold(f64::MIN, f64::max);
    let g3_mean = g3_curves.h[0].iter().sum::<f64>() / g3_curves.h[0].len() as f64;
    push_gate(&mut gates, "G3.a", "stationary phi11: max H32 (floor 0.0256)", format!("{g3_max:.4}"), "< 2x floor (0.0512)".into(), g3_max < 2.0 * floor32);
    push_gate(&mut gates, "G3.b", "stationary phi11: mean H32", format!("{g3_mean:.4}"), "< 1.5x floor (0.0384)".into(), g3_mean < 1.5 * floor32);

    // =======================================================================
    // G4 — the H-theorem demo (M = 16 + equilibrium control, T = 4 pi).
    // =======================================================================
    println!("== G4: Nelson H-theorem demo (M = 16, N = {N_PART}, T = 4 pi) ==");
    let modes16 = ModeSystem::born_m16();
    let main_curves;
    {
        let ws = WienerStream::new(SEED_PHASES, K_NELSON_M16, 2);
        let mut xs = vec![0.0f64; N_PART];
        let mut ys = vec![0.0f64; N_PART];
        for i in 0..N_PART {
            let (x, y) = sample_phi11(SEED_SAMPLE, K_INIT_PHI11, i as u32);
            xs[i] = x;
            ys[i] = y;
        }
        main_curves = run_nelson_ensemble(&modes16, &ws, &mut xs, &mut ys, T_FINAL, "M=16 nelson");
    }
    let ctrl_curves;
    {
        let ws = WienerStream::new(SEED_PHASES, K_NELSON_CTRL, 2);
        let c0 = modes16.coeffs(0.0);
        let bound = modes16.psi2_max_bound(0.0);
        let mut xs = vec![0.0f64; N_PART];
        let mut ys = vec![0.0f64; N_PART];
        for i in 0..N_PART {
            let (x, y) = sample_psi2(&modes16, &c0, bound, SEED_SAMPLE, K_INIT_PSI2, i as u32);
            xs[i] = x;
            ys[i] = y;
        }
        ctrl_curves = run_nelson_ensemble(&modes16, &ws, &mut xs, &mut ys, T_FINAL, "M=16 control");
    }

    let fits: Vec<TauFit> = (0..2).map(|ci| fit_tau(&main_curves.ts, &main_curves.h[ci])).collect();
    let h32_final = *main_curves.h[0].last().unwrap();
    let ctrl_max = ctrl_curves.h[0].iter().cloned().fold(f64::MIN, f64::max);
    let ctrl_mean = ctrl_curves.h[0].iter().sum::<f64>() / ctrl_curves.h[0].len() as f64;

    push_gate(&mut gates, "G4.a", "M=16 cg32 exponential-fit r^2", format!("{:.4}", fits[0].r2), ">= 0.95".into(), fits[0].r2 >= 0.95);
    push_gate(&mut gates, "G4.b", "M=16 H32(4 pi) (floor 0.0256)", format!("{h32_final:.4}"), "< 2x floor (0.0512)".into(), h32_final < 2.0 * floor32);
    println!("  control H32: mean {ctrl_mean:.4}, max {ctrl_max:.4} (floor {floor32:.4})");
    push_gate(&mut gates, "G4.c", "control max H32 (equilibrium stays at floor)", format!("{ctrl_max:.4}"), "< 2x floor (0.0512)".into(), ctrl_max < 2.0 * floor32);
    let tau32 = fits[0].tau;
    push_gate(&mut gates, "G4.d", "Nelson tau (cg32) vs Bohm 4.853 (expect <=)", format!("{tau32:.3}"), "(0, 6.07]".into(), tau32 > 0.0 && tau32 <= 1.25 * BOHM_TAU_CG32);

    println!(
        "  Nelson vs de Broglie-Bohm relaxation (cg32): tau_nelson = {:.3}, tau_bohm = {:.3}, ratio {:.3}",
        tau32,
        BOHM_TAU_CG32,
        tau32 / BOHM_TAU_CG32
    );

    // =======================================================================
    // Outputs.
    // =======================================================================
    let mut hash = FNV_OFFSET;
    hash = fnv1a_f64(hash, &main_curves.h[0]);
    hash = fnv1a_f64(hash, &main_curves.h[1]);
    hash = fnv1a_f64(hash, &ctrl_curves.h[0]);
    hash = fnv1a_f64(hash, &ctrl_curves.h[1]);
    hash = fnv1a_f64(hash, &g3_curves.h[0]);
    let fingerprint = format!("{hash:016x}");
    println!("  determinism fingerprint (FNV-1a 64 over all H curves): {fingerprint}");

    // hdata.csv (t, H32/H16 for nelson + control).
    let mut csv = String::from("t,H32_nelson16,H16_nelson16,H32_ctrl,H16_ctrl\n");
    for i in 0..main_curves.ts.len() {
        let _ = writeln!(
            csv,
            "{},{},{},{},{}",
            main_curves.ts[i], main_curves.h[0][i], main_curves.h[1][i], ctrl_curves.h[0][i], ctrl_curves.h[1][i]
        );
    }
    std::fs::write(format!("{outdir}/hdata.csv"), csv).expect("write hdata.csv");

    // results.json (deterministic content only).
    let mut js = String::from("{\n");
    let _ = writeln!(js, "  \"n_part\": {N_PART},");
    let _ = writeln!(js, "  \"dt\": {DT},");
    let _ = writeln!(js, "  \"t_final\": {T_FINAL},");
    let _ = writeln!(js, "  \"seed_phases\": {SEED_PHASES},");
    let _ = writeln!(js, "  \"seed_sample\": {SEED_SAMPLE},");
    let _ = writeln!(js, "  \"noise_floor\": {{\"cg32\": {}, \"cg16\": {}}},", noise_floor(32, N_PART), noise_floor(16, N_PART));
    let _ = writeln!(js, "  \"ou\": {{");
    let _ = writeln!(js, "    \"mean_exact\": {mean_exact}, \"var_exact\": {var_exact},");
    let _ = writeln!(js, "    \"em\": {{\"mean\": {em_mean}, \"var\": {em_var}}},");
    let _ = writeln!(js, "    \"sra1\": {{\"mean\": {sra_mean}, \"var\": {sra_var}}}");
    let _ = writeln!(js, "  }},");
    let _ = writeln!(js, "  \"orders\": {{");
    let _ = writeln!(js, "    \"strong_em_ou\": {p_em}, \"strong_sra1_ou\": {p_sra},");
    let _ = writeln!(js, "    \"strong_em_cubic\": {cp_em}, \"strong_sra1_cubic\": {cp_sra}, \"strong_em_gbm\": {p_gbm},");
    let _ = writeln!(js, "    \"weak_em_ou\": {q_em}, \"weak_sra1_ou\": {q_sra}");
    let _ = writeln!(js, "  }},");
    let _ = writeln!(js, "  \"g3_phi11\": {{\"H32_max\": {g3_max}, \"H32_mean\": {g3_mean}}},");
    for (name, curves, fit) in [("nelson_m16", &main_curves, Some(&fits)), ("control", &ctrl_curves, None)] {
        let _ = writeln!(js, "  \"{name}\": {{");
        for (ci, cg) in CGS.iter().enumerate() {
            let h = &curves.h[ci];
            let hmax = h.iter().cloned().fold(f64::MIN, f64::max);
            let hmean = h.iter().sum::<f64>() / h.len() as f64;
            let _ = write!(
                js,
                "    \"cg{cg}\": {{\"H_initial\": {}, \"H_final\": {}, \"H_max\": {hmax}, \"H_mean\": {hmean}",
                h[0],
                h.last().unwrap()
            );
            if let Some(f) = fit {
                let _ = write!(
                    js,
                    ", \"tau\": {}, \"H0_fit\": {}, \"fit_window_t_end\": {}, \"fit_r2\": {}",
                    f[ci].tau, f[ci].h0_fit, f[ci].t_end, f[ci].r2
                );
            }
            let _ = writeln!(js, "}}{}", if ci == 0 { "," } else { "" });
        }
        let _ = writeln!(js, "  }},");
    }
    let _ = writeln!(js, "  \"bohm_tau_cg32\": {BOHM_TAU_CG32},");
    let _ = writeln!(js, "  \"tau_ratio_nelson_over_bohm\": {},", tau32 / BOHM_TAU_CG32);
    let _ = writeln!(js, "  \"determinism_fnv64\": \"{fingerprint}\",");
    let _ = writeln!(js, "  \"gates\": [");
    for (i, g) in gates.iter().enumerate() {
        let _ = writeln!(
            js,
            "    {{\"id\": \"{}\", \"desc\": \"{}\", \"value\": \"{}\", \"target\": \"{}\", \"pass\": {}}}{}",
            g.id,
            g.desc,
            g.value,
            g.target,
            g.pass,
            if i + 1 < gates.len() { "," } else { "" }
        );
    }
    let _ = writeln!(js, "  ]");
    let _ = writeln!(js, "}}");
    std::fs::write(format!("{outdir}/results.json"), js).expect("write results.json");

    // =======================================================================
    // Gate table.
    // =======================================================================
    println!("\n== GATE TABLE ==");
    println!("{:<6} {:<52} {:>12} {:>22} {:>6}", "id", "gate", "value", "target", "pass");
    let mut all_pass = true;
    for g in &gates {
        all_pass &= g.pass;
        println!("{:<6} {:<52} {:>12} {:>22} {:>6}", g.id, g.desc, g.value, g.target, if g.pass { "PASS" } else { "FAIL" });
    }
    println!(
        "\n{} gates, {} passed; total wall time {:.1} s",
        gates.len(),
        gates.iter().filter(|g| g.pass).count(),
        wall.elapsed().as_secs_f64()
    );
    if !all_pass {
        std::process::exit(1);
    }
}
