//! GUM core Phase G1 — `fs-gum-kern`: the deterministic parallel kernel
//! layer for the N^3 sweeps of `fs-gum-field` / `fs-gum-statics`
//! (ROADMAP_v4_SCALE.md Phase G1; design measured and fixed by
//! scale_survey_v4.json DIMENSIONs 1-2).
//!
//! ONE parallel primitive: a domain-decomposed sweep over evaluation
//! points — per-tile serial accumulation in a thread-count-INDEPENDENT
//! i-slab decomposition, combined by a fixed pairwise tree — plus the
//! gradient adjoint inverted from scatter to GATHER.  Serial and threaded
//! execution are bit-identical BY CONSTRUCTION: every rounding decision
//! is a pure function of (extent, [`TILE_I`]) and never of thread count,
//! scheduling, or steal order (gate K-D asserts hash equality at 1/2/3/4
//! threads).
//!
//! THE BIT CONTRACT (the fs-la `KC`/`GEMM_BIT_SEMANTICS_VERSION`
//! precedent — these constants legitimately change bits; retuning any of
//! them is a deliberate, documented golden bump):
//!
//!  * [`TILE_I`] = 4 — i-slab thickness in planes; tile count
//!    `ceil(m / TILE_I)` is a pure function of the sweep extent `m`;
//!  * the combine tree — largest-power-of-two-strictly-below split,
//!    ascending tile order (vendored fs-exec shape rule);
//!  * per-tile accumulation — plain f64, ascending (i, j, k), the exact
//!    per-point statements of the old engine;
//!  * the gather order — per output cell: axis x, y, z with eval-point
//!    offsets (-1, +1, -2, +2), then own gq (central4); the 8 corners in
//!    ascending (ei, ej, ek) with the fused flux+gq corner expression
//!    (corner scheme); then cell terms, h^3, tangent projection;
//!  * `norm_flat` per-tile order — (a, i, j, k) ascending.
//!
//! All of it is pinned by [`GUM_KERN_BIT_SEMANTICS`], which domain-tags
//! the frozen goldens in the gates binary.
//!
//! THE DELIBERATE GOLDEN BUMP (documented per the fs-la KC-contract
//! precedent): the tiled canonical order provably differs from the old
//! flat sequential order (measured at the ~12th significant digit —
//! survey DIMENSION 2, prototype `tiledet`; re-measured by gate K-A),
//! and the gather adjoint re-orders each cell's contribution sum.  The
//! gates binary cross-checks the new serial path against the OLD
//! `fs-gum-statics` engine ONCE at machine-eps class, states the measured
//! deltas, then freezes NEW bit goldens for the tiled order.  The old
//! crates and their goldens are untouched.
//!
//! Threading: `std::thread::scope` + an `AtomicUsize` tile dispenser with
//! cache-padded per-tile slots (the fs-la band-dispenser pattern;
//! `workers = threads.min(tiles)`); no rayon, no fs-exec (unbuildable
//! sibling — its `pairwise_fold`/`Compensated` are VENDORED in
//! [`reduce`] with a documented swap-back, the B3/F3 precedent).  The
//! crate is safe Rust throughout: tile disjointness is established by
//! slice splitting, not unsafe.
//!
//! Epistemic notice (binding, inherited): everything here is
//! within-model numerical engineering on a speculative theory's
//! functional.  A passing gate certifies the kernel layer — determinism,
//! equivalence, speedup — never anything about nature.

#![forbid(unsafe_code)]

pub mod anf;
pub mod engine;
pub mod reduce;
pub mod tile;

pub use anf::anf;
pub use engine::{axpy_sub, eval, eval_ws, norm_flat, renormalize, step_q, tangent_project, Ws};
pub use reduce::{pairwise_fold, Compensated, Reduce};
pub use tile::{sweep_reduce, tile_count, tile_range, tiled_run, CachePadded};

/// i-slab tile thickness in planes — PART OF THE BIT CONTRACT (see
/// module docs).  4 planes: >= the widest stencil reach (central4 reads
/// +-2), fine-grained enough for a 4-way dispenser at N = 48
/// (13 corner tiles), and the value the survey prototype validated.
pub const TILE_I: usize = 4;

/// The bit-semantics version string: pins TILE_I, the combine-tree shape,
/// the per-tile accumulation order, and the gather order.  Any change to
/// any of them bumps this string and forces a deliberate golden re-freeze
/// (the fs-la GEMM_BIT_SEMANTICS_VERSION discipline).  Frozen goldens are
/// hashed under this domain.
pub const GUM_KERN_BIT_SEMANTICS: &str = "gum-kern-bits-v1: TILE_I=4 i-slabs; \
    largest-pow2-below pairwise combine in ascending tile order; per-tile plain-f64 \
    ascending (i,j,k); gather adjoint per cell: axes x,y,z offsets (-1,+1,-2,+2) then gq \
    (central4) / corners ascending (ei,ej,ek) fused flux+gq (corner), then cell terms, \
    h^3, tangent projection; norm_flat per-tile (a,i,j,k)";
