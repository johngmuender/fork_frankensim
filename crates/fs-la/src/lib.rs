//! fs-la — Dense linear algebra: GEMM, batched small dense, factorizations, eigensolvers.
//!
//! Layer: L1. See CONTRACT.md for invariants, error model, determinism
//! class, cancellation behavior, and no-claim boundaries. This crate is part
//! of the FrankenSim workspace; the layer dependency direction is enforced by
//! `cargo run -p xtask -- check-all`.

pub mod batched;
pub mod batched_f32;
pub mod eigen;
pub mod eigen_complex;
pub mod factor;
pub mod gemm;
pub mod mixed;
pub mod rand_nla;

pub use gemm::{
    GEMM_DEPGRAPH_RECEIPT_DOMAIN, GEMM_IMPLEMENTATION_VERSION, GEMM_MAX_FMAS_BETWEEN_POLLS,
    GEMM_PANEL_RUN_DOMAIN, GemmMemoryEnvelope, GemmMemoryReport, Trans, gemm_execution_tier,
    gemm_f32, gemm_f64, gemm_f64_op, gemm_f64_parallel, gemm_f64_parallel_with, gemm_mixed,
    gemm_tuning_is_effective,
};

// Executor-coupled surface (the `exec` feature, ON by default): the
// cancellation-aware pooled GEMM API and the build-identity constants whose
// fingerprint inputs include the ../asupersync sibling checkout. Compiled out
// under `default-features = false` so fs-la stays buildable from a checkout
// without that sibling repo (Phase F3 portability).
#[cfg(feature = "exec")]
pub use gemm::{
    GEMM_BUILD_FINGERPRINT, GEMM_DEPGRAPH_RECEIPT, GEMM_DEPGRAPH_RECEIPT_DIGEST,
    GEMM_GRAPH_EVIDENCE, GEMM_GRAPH_EVIDENCE_KIND, GemmCancelled, GemmGraphEvidence,
    GemmGraphEvidenceClass, GemmRunError, GemmRunReport, gemm_build_identity,
    gemm_f64_parallel_with_cancel, gemm_f64_parallel_with_pool,
    gemm_f64_parallel_with_pool_budgeted, gemm_f64_parallel_with_pool_declared,
    gemm_graph_evidence, gemm_panel_run_id,
};

/// Crate version, re-exported for provenance stamping (the Five Explicits'
/// "versions" pillar reaches down to individual crates).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#[cfg(test)]
mod tests {
    #[test]
    fn version_is_stamped() {
        assert!(!super::VERSION.is_empty());
    }
}
