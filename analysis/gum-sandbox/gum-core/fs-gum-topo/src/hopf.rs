//! Hopf invariant, Method A: the Whitehead integral with an FFT
//! Coulomb-gauge curl-inverse on a periodic grid.
//!
//! Pipeline (survey build_spec, Algorithm 2A):
//!   1. F_ij = n . (d_i n x d_j n) with O(h^2) central differences;
//!      B_k = (1/2) eps_kij F_ij, so B = (F_yz, F_zx, F_xy).
//!      Normalization: the flux of B through a sphere around a unit
//!      hedgehog is 4 pi.
//!   2. Obstruction check BEFORE the solve: the net flux of B through
//!      each torus 2-cycle must vanish, or the Hopf invariant is
//!      undefined on T^3 and a typed error is returned, never a number.
//!      The survey's raw mean-B test (tol 1e-10 x max|B|) cannot work
//!      on the central-difference B, whose box mean carries an
//!      O(h^{3/2}) discretization residue (measured ~1e-2 on the
//!      sqrt-cusp referee); the flux is therefore measured with the
//!      LATTICE-EXACT Berg-Luscher solid-angle winding of the three
//!      coordinate-slice maps — an integer up to floating point, which
//!      separates topology from noise at machine precision. The raw
//!      mean is still computed and reported.
//!   3. Coulomb-gauge vector potential via the crate's minimal radix-2
//!      FFT (crates/fs-fft is not usable standalone — see src/fft.rs;
//!      axes must be powers of two, structured rejection otherwise):
//!      A-hat(k) = i (kappa x B-hat) / |kappa|^2 with the DISCRETE
//!      wavevector kappa_m = sin(2 pi m / N) / h — the symbol of the
//!      central difference, so the discrete curl of A reproduces B to
//!      machine precision on the solenoidal part instead of O(h^2).
//!      A-hat = 0 at kappa = 0 (the k = 0 zero mode AND the Nyquist
//!      combinations where every sine vanishes): with the step-2 zero
//!      mean enforced, this makes the answer canonical (the A -> A +
//!      const gauge freedom contributes const . INT B = 0).
//!   4. H = SIGN_WHITEHEAD * (1/(4 pi)^2) INT A . B d^3x, midpoint sum,
//!      fixed (i, j, k) loop order; the sign is fixed once against the
//!      analytic Hopf-1 referee.

use crate::fft::{C64, Fft3};
use crate::field::{Bc, DirectorField, Field3};
use crate::{v3, TopoError, SIGN_WHITEHEAD};
use core::f64::consts::PI;

/// Result of the Whitehead route.
#[derive(Debug, Clone, Copy)]
pub struct WhiteheadOut {
    /// The Hopf invariant estimate (sign-fixed).
    pub hopf: f64,
    /// The raw integral (1/(4 pi)^2) INT A . B before the sign fix.
    pub raw: f64,
    /// Lattice-exact slice winding integers (net flux / 4 pi) through
    /// the (yz, zx, xy) 2-cycles — all zero for an admissible field.
    pub flux: [f64; 3],
    /// Component-wise box mean of the central-difference B (the
    /// discretization-residue monitor).
    pub mean_b: [f64; 3],
    /// max |B| over the box.
    pub max_b: f64,
    /// Relative L2 error of curl_h A against B (Nyquist content and the
    /// non-solenoidal O(h^2) part of the discrete B live here).
    pub curl_rel_err: f64,
    /// Largest imaginary residue left by the inverse FFT (sanity).
    pub max_imag: f64,
}

fn is_pow2(x: usize) -> bool {
    x >= 1 && (x & (x - 1)) == 0
}

/// Central-difference director derivatives at a node (periodic get).
#[inline]
fn dn(f: &DirectorField, i: isize, j: isize, k: isize, inv2h: f64) -> [[f64; 3]; 3] {
    let dx = v3::scale(v3::sub(f.get(i + 1, j, k), f.get(i - 1, j, k)), inv2h);
    let dy = v3::scale(v3::sub(f.get(i, j + 1, k), f.get(i, j - 1, k)), inv2h);
    let dz = v3::scale(v3::sub(f.get(i, j, k + 1), f.get(i, j, k - 1)), inv2h);
    [dx, dy, dz]
}

/// Lattice-exact geometric winding (net flux / 4 pi) of the director
/// map restricted to the coordinate slice `fixed_axis = 0`: sum of the
/// signed spherical-triangle areas of the plaquette corners — an exact
/// integer up to floating point (the Berg-Luscher kernel).
fn slice_winding(f: &DirectorField, fixed_axis: usize) -> f64 {
    let bx = (fixed_axis + 1) % 3;
    let cx = (fixed_axis + 2) % 3;
    let nb = f.n[bx];
    let nc = f.n[cx];
    let read = |b: usize, c: usize| -> [f64; 3] {
        let mut v = [0usize; 3];
        v[bx] = b % nb;
        v[cx] = c % nc;
        crate::v3::normalize(f.at(v[0], v[1], v[2]))
    };
    let mut omega = 0.0;
    for b in 0..nb {
        for c in 0..nc {
            let n00 = read(b, c);
            let n10 = read(b + 1, c);
            let n11 = read(b + 1, c + 1);
            let n01 = read(b, c + 1);
            omega += crate::solid_angle_origin(n00, n10, n11)
                + crate::solid_angle_origin(n00, n11, n01);
        }
    }
    omega / (4.0 * PI)
}

/// The Faraday/curvature field B of a director field:
/// B = (n.(dy n x dz n), n.(dz n x dx n), n.(dx n x dy n)).
fn b_field(f: &DirectorField) -> Vec<[f64; 3]> {
    let [n0, n1, n2] = f.n;
    let inv2h = 1.0 / (2.0 * f.h);
    let mut out = Vec::with_capacity(n0 * n1 * n2);
    for i in 0..n0 {
        for j in 0..n1 {
            for k in 0..n2 {
                let nv = f.at(i, j, k);
                let d = dn(f, i as isize, j as isize, k as isize, inv2h);
                out.push([
                    v3::dot(nv, v3::cross(d[1], d[2])),
                    v3::dot(nv, v3::cross(d[2], d[0])),
                    v3::dot(nv, v3::cross(d[0], d[1])),
                ]);
            }
        }
    }
    out
}

/// The Whitehead-integral Hopf invariant of a periodic director field.
///
/// # Errors
/// [`TopoError::NotPeriodic`] for non-periodic BCs,
/// [`TopoError::NotPow2`] for non-power-of-two axes, and
/// [`TopoError::NetFlux`] when the torus obstruction is nonzero.
pub fn hopf_whitehead(f: &DirectorField) -> Result<WhiteheadOut, TopoError> {
    if f.bc != Bc::Periodic {
        return Err(TopoError::NotPeriodic);
    }
    let [n0, n1, n2] = f.n;
    if !(is_pow2(n0) && is_pow2(n1) && is_pow2(n2)) {
        return Err(TopoError::NotPow2 { dims: f.n });
    }
    let total = n0 * n1 * n2;
    let h = f.h;

    // Step 1: B at every node, fixed loop order.
    let b = b_field(f);

    // Step 2: obstruction check before any solve.
    let mut mean = [0.0_f64; 3];
    let mut max_b = 0.0_f64;
    for v in &b {
        for c in 0..3 {
            mean[c] += v[c];
        }
        let m = v3::norm(*v);
        if m > max_b {
            max_b = m;
        }
    }
    for c in &mut mean {
        *c /= total as f64;
    }
    let flux = [
        slice_winding(f, 0),
        slice_winding(f, 1),
        slice_winding(f, 2),
    ];
    if flux.iter().any(|q| q.abs() > 0.5) {
        return Err(TopoError::NetFlux { flux, mean, max_b });
    }
    if flux.iter().any(|q| q.abs() > 1.0e-6) {
        // Not an integer at all: a slice map grazes a degeneracy.
        return Err(TopoError::Degenerate(format!(
            "slice winding not integer-snapped: {flux:?}"
        )));
    }
    if max_b == 0.0 {
        // Uniform (or curvature-free) director: H = 0 exactly.
        return Ok(WhiteheadOut {
            hopf: 0.0,
            raw: 0.0,
            flux,
            mean_b: mean,
            max_b,
            curl_rel_err: 0.0,
            max_imag: 0.0,
        });
    }

    // Step 3: Coulomb-gauge A via the discrete-wavevector curl-inverse.
    let plan = Fft3::new([n0, n1, n2]);
    let mut spec: [Vec<C64>; 3] = [
        Vec::with_capacity(total),
        Vec::with_capacity(total),
        Vec::with_capacity(total),
    ];
    for c in 0..3 {
        spec[c].extend(b.iter().map(|v| C64::new(v[c], 0.0)));
        plan.forward(&mut spec[c]);
    }
    let tau = 2.0 * PI;
    let kax: Vec<f64> = (0..n0).map(|m| (tau * m as f64 / n0 as f64).sin() / h).collect();
    let kay: Vec<f64> = (0..n1).map(|m| (tau * m as f64 / n1 as f64).sin() / h).collect();
    let kaz: Vec<f64> = (0..n2).map(|m| (tau * m as f64 / n2 as f64).sin() / h).collect();
    for i in 0..n0 {
        for j in 0..n1 {
            for k in 0..n2 {
                let idx = (i * n1 + j) * n2 + k;
                let kv = [kax[i], kay[j], kaz[k]];
                let k2 = kv[0] * kv[0] + kv[1] * kv[1] + kv[2] * kv[2];
                if k2 < 1.0e-30 {
                    // Zero mode / all-Nyquist combinations: A-hat = 0.
                    for s in &mut spec {
                        s[idx] = C64::new(0.0, 0.0);
                    }
                    continue;
                }
                let bh = [spec[0][idx], spec[1][idx], spec[2][idx]];
                // A-hat = i (kappa x B-hat) / |kappa|^2.
                let cxr = kv[1] * bh[2].re - kv[2] * bh[1].re;
                let cxi = kv[1] * bh[2].im - kv[2] * bh[1].im;
                let cyr = kv[2] * bh[0].re - kv[0] * bh[2].re;
                let cyi = kv[2] * bh[0].im - kv[0] * bh[2].im;
                let czr = kv[0] * bh[1].re - kv[1] * bh[0].re;
                let czi = kv[0] * bh[1].im - kv[1] * bh[0].im;
                spec[0][idx] = C64::new(-cxi / k2, cxr / k2);
                spec[1][idx] = C64::new(-cyi / k2, cyr / k2);
                spec[2][idx] = C64::new(-czi / k2, czr / k2);
            }
        }
    }
    let mut max_imag = 0.0_f64;
    let mut adata: Vec<[f64; 3]> = vec![[0.0; 3]; total];
    for c in 0..3 {
        plan.inverse(&mut spec[c]);
        for (av, sv) in adata.iter_mut().zip(spec[c].iter()) {
            av[c] = sv.re;
            let im = sv.im.abs();
            if im > max_imag {
                max_imag = im;
            }
        }
    }

    // Curl-reproduction diagnostic: curl_h A vs B, relative L2.
    let afield = Field3::from_parts(f.n, h, f.origin, Bc::Periodic, adata);
    let inv2h = 1.0 / (2.0 * h);
    let (mut num, mut den) = (0.0_f64, 0.0_f64);
    for i in 0..n0 {
        let ii = i as isize;
        for j in 0..n1 {
            let jj = j as isize;
            for k in 0..n2 {
                let kk = k as isize;
                let d = dn(&afield, ii, jj, kk, inv2h);
                // curl A = (dA_z/dy - dA_y/dz, dA_x/dz - dA_z/dx,
                //           dA_y/dx - dA_x/dy).
                let curl = [
                    d[1][2] - d[2][1],
                    d[2][0] - d[0][2],
                    d[0][1] - d[1][0],
                ];
                let bv = b[afield.idx(i, j, k)];
                let e = v3::sub(curl, bv);
                num += v3::dot(e, e);
                den += v3::dot(bv, bv);
            }
        }
    }
    let curl_rel_err = (num / den).sqrt();

    // Step 4: H = sign * (1/(4 pi)^2) INT A . B, midpoint sum.
    let mut s = 0.0_f64;
    for (av, bv) in afield.data.iter().zip(b.iter()) {
        s += v3::dot(*av, *bv);
    }
    let raw = s * h * h * h / (16.0 * PI * PI);
    Ok(WhiteheadOut {
        hopf: SIGN_WHITEHEAD * raw,
        raw,
        flux,
        mean_b: mean,
        max_b,
        curl_rel_err,
        max_imag,
    })
}
