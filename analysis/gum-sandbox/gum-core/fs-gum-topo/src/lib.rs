//! fs-gum-topo — GUM core Phase B2: topological-charge diagnostics on
//! STORED fields (gap_survey_v3.json, dimension "Topological-charge
//! diagnostics").
//!
//! What lives here, all on a minimal self-contained grid container
//! ([`Field3`], flat `Vec<[f64; C]>`, periodic or vacuum-padded BCs):
//!
//!   1. the Hopf map / director projection n = R(q) z-hat (S3 -> S2)
//!      on stored quaternion fields ([`hopf_project`]);
//!   2. S3 degree b_P on stored quaternion fields: the pilot's
//!      central-difference b-density (ported bit-compatibly from
//!      fs-cosserat-pilot), the central4 and corner (midpoint) schemes
//!      ported from field3d_solve.py, and a geometric van Oosterom-
//!      Strackee solid-angle (integer-snapping) hedgehog charge;
//!   3. the Hopf invariant of S2-valued director fields, BOTH ways:
//!      (a) the Whitehead integral with an FFT Coulomb-gauge
//!      curl-inverse (discrete wavevector, zero-mode handling, net-flux
//!      obstruction check) and (b) preimage linking (marching-tetrahedra
//!      curve extraction + exact polyline Gauss-linking sum);
//!   4. disclination diagnostics for RP2 director fields: Z/2 loop
//!      holonomy by sign-gauge transport, half-integer planar winding
//!      for framed loops, and dual-plaquette defect-line tracing with
//!      the lines-cannot-end structural check.
//!
//! Sign conventions are fixed ONCE against the analytic referees (the
//! survey's rule): the compacton hedgehog has degree +1 and its Hopf
//! projection has H = +1 in both methods; see [`SIGN_WHITEHEAD`] and
//! [`SIGN_LINK`].
//!
//! Epistemic notice (binding): every PASS in the gates binary certifies
//! a WITHIN-MODEL property of diagnostic machinery on a speculative
//! theory's analytic configurations. This validates the discretisation
//! and the evidence pipeline — never the physics.

#![forbid(unsafe_code)]

pub mod certify;
pub mod degree;
pub mod disclination;
pub mod field;
pub mod hopf;
pub mod linking;
pub mod referee;

pub use field::{hopf_project, Bc, DirectorField, Field3, QField};

use core::fmt;

/// Overall sign of the Whitehead integral, fixed once against the
/// analytic Hopf-1 referee (hopf_project of the degree-+1 compacton
/// hedgehog must give H = +1). With this factor the parity-mirrored
/// texture measures -1, as it must.
pub const SIGN_WHITEHEAD: f64 = 1.0;

/// Overall sign of the preimage-linking route, fixed once against the
/// same Hopf-1 referee (both methods must agree on +1).
pub const SIGN_LINK: f64 = 1.0;

/// Typed failures of the topological diagnostics. Obstructed or
/// degenerate configurations return one of these — never a number.
#[derive(Debug, Clone, PartialEq)]
pub enum TopoError {
    /// The mean of B over the periodic box is nonzero: net flux through
    /// a torus 2-cycle, so the Hopf invariant is undefined on T^3.
    NetFlux {
        /// Component-wise box mean of B.
        mean: [f64; 3],
        /// max |B| over the box (the scale the tolerance is relative to).
        max_b: f64,
    },
    /// The FFT route needs power-of-two axes (fs-fft constraint).
    NotPow2 {
        /// The offending grid dimensions.
        dims: [usize; 3],
    },
    /// The Whitehead route is defined here on periodic grids only.
    NotPeriodic,
    /// Loop transport hit a defect core: |n_a . n_b| below tolerance.
    CoreHit {
        /// Index of the loop step at which transport degenerated.
        step: usize,
        /// The offending |dot| value.
        dot: f64,
    },
    /// A geometric degeneracy (non-regular value, edge-grazing crossing,
    /// junction, ...). The message says which; preimage extraction
    /// retries with a re-tilted regular value before giving up.
    Degenerate(String),
    /// A preimage curve failed to close inside the marching region —
    /// either it wraps a torus cycle (same obstruction as `NetFlux`) or
    /// it leaves through the boundary layer.
    OpenCurve,
}

impl fmt::Display for TopoError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            TopoError::NetFlux { mean, max_b } => write!(
                f,
                "net B-flux obstruction: mean B = [{:.3e}, {:.3e}, {:.3e}], max|B| = {:.3e}",
                mean[0], mean[1], mean[2], max_b
            ),
            TopoError::NotPow2 { dims } => {
                write!(f, "FFT axes must be powers of two, got {dims:?}")
            }
            TopoError::NotPeriodic => write!(f, "Whitehead route requires periodic BCs"),
            TopoError::CoreHit { step, dot } => {
                write!(f, "defect core hit at loop step {step}: |dot| = {dot:.3e}")
            }
            TopoError::Degenerate(msg) => write!(f, "degeneracy: {msg}"),
            TopoError::OpenCurve => write!(f, "preimage curve failed to close"),
        }
    }
}

/// Minimal 3-vector helpers on `[f64; 3]` (kept local so the crate is
/// self-contained; fs-geom's Point3/Vec3 can absorb these at the later
/// integration pass).
pub mod v3 {
    /// Dot product.
    #[inline]
    #[must_use]
    pub fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
        a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    }

    /// Cross product.
    #[inline]
    #[must_use]
    pub fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    }

    /// a - b.
    #[inline]
    #[must_use]
    pub fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
    }

    /// a + b.
    #[inline]
    #[must_use]
    pub fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
    }

    /// s * a.
    #[inline]
    #[must_use]
    pub fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
        [a[0] * s, a[1] * s, a[2] * s]
    }

    /// Euclidean norm.
    #[inline]
    #[must_use]
    pub fn norm(a: [f64; 3]) -> f64 {
        dot(a, a).sqrt()
    }

    /// a / |a| (caller guarantees |a| > 0).
    #[inline]
    #[must_use]
    pub fn normalize(a: [f64; 3]) -> [f64; 3] {
        scale(a, 1.0 / norm(a))
    }
}

/// Signed solid angle of the (spherical) triangle spanned by directions
/// `a, b, c` seen from the origin — the van Oosterom–Strackee kernel
/// (the 8-line kernel copied from crates/fs-rep-mesh/src/winding.rs, as
/// the survey recommends, to avoid the mesh dependency). Exact up to
/// floating point; the vectors need not be unit.
#[must_use]
pub fn solid_angle_origin(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let (la, lb, lc) = (v3::norm(a), v3::norm(b), v3::norm(c));
    let det = v3::dot(a, v3::cross(b, c));
    let denom = la * lb * lc
        + v3::dot(a, b) * lc
        + v3::dot(b, c) * la
        + v3::dot(c, a) * lb;
    2.0 * det.atan2(denom)
}

/// 3x3 determinant of rows a, b, c (shared by degree and linking code).
#[inline]
#[must_use]
pub fn det3(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
}

/// 4x4 determinant with rows (q, dq/dx, dq/dy, dq/dz), first-row
/// expansion — bit-identical port of fs-cosserat-pilot's `det4` (the
/// same operation order feeds the golden-digit reproduction gate).
#[must_use]
pub fn det4(r0: &[f64; 4], r1: &[f64; 4], r2: &[f64; 4], r3: &[f64; 4]) -> f64 {
    let minor = |r: &[f64; 4], skip: usize| -> [f64; 3] {
        let mut m = [0.0; 3];
        let mut j = 0;
        for (c, &v) in r.iter().enumerate() {
            if c != skip {
                m[j] = v;
                j += 1;
            }
        }
        m
    };
    r0[0] * det3(minor(r1, 0), minor(r2, 0), minor(r3, 0))
        - r0[1] * det3(minor(r1, 1), minor(r2, 1), minor(r3, 1))
        + r0[2] * det3(minor(r1, 2), minor(r2, 2), minor(r3, 2))
        - r0[3] * det3(minor(r1, 3), minor(r2, 3), minor(r3, 3))
}
