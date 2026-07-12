//! VENDORED deterministic reduction core — `pairwise_fold` + `Compensated`
//! from `crates/fs-exec/src/reduce.rs` (the workspace's deterministic
//! reduction library), trimmed to the ~100 lines this crate needs.
//!
//! WHY VENDORED (the B3 `jacobi_eigh` / F3 `exec`-feature precedent,
//! recorded again in GAP_ANALYSIS_v4_SCALE.md finding 7): fs-exec is
//! unbuildable in this checkout — its manifest path-deps on the absent
//! `../../../asupersync` sibling, so ANY transitive dep on it poisons the
//! build.  SWAP-BACK CONTRACT: if fs-exec (or a dep-free `fs-reduce`
//! split of its reduce module) becomes buildable here, delete this file
//! and re-export `fs_exec::reduce::{pairwise_fold, Compensated}` — the
//! shapes below are copied verbatim, so the swap is bit-neutral.
//!
//! THE SHAPE RULE (part of the bit contract): a fold over `n` ordered
//! items splits at the largest power of two strictly below `n` and
//! recurses — a pure function of `n`, never of scheduling, arrival order,
//! or thread count.  Every rounding decision in a reduction built on this
//! tree is replayable from the item count alone.

/// Merge half of the vendored `fs_exec::kernel::Reduce` trait: an
/// identity element and an order-sensitive binary merge.  Merges see
/// operands in ascending logical (tile) order.
pub trait Reduce {
    /// The identity partial.
    fn identity() -> Self;
    /// Merge two partials (left operand is the lower logical index).
    #[must_use]
    fn merge(self, other: Self) -> Self;
}

/// Fold ordered items up the fixed-shape pairwise tree.  The tree depends
/// only on `items.len()`, so the result is bit-identical for a given
/// input sequence regardless of how the inputs were produced (fs-exec
/// reduce.rs, verbatim).
#[must_use]
pub fn pairwise_fold<T: Reduce>(mut items: Vec<T>) -> T {
    match items.len() {
        0 => T::identity(),
        1 => items.pop().expect("len checked"),
        n => {
            let split = largest_pow2_below(n);
            let right = items.split_off(split);
            pairwise_fold(items).merge(pairwise_fold(right))
        }
    }
}

/// Largest power of two strictly below `n` (n >= 2).
fn largest_pow2_below(n: usize) -> usize {
    debug_assert!(n >= 2);
    let p = usize::BITS - (n - 1).leading_zeros() - 1;
    1 << p
}

/// A Neumaier-compensated partial sum (fs-exec reduce.rs, verbatim):
/// `Reduce`-composable, so compensated global sums ride the same fixed
/// tree as everything else.  The gum kernels use PLAIN per-tile f64
/// accumulation to match the current physics semantics (survey DIMENSION
/// 2 recommendation); this type is offered for future kernels whose
/// physics wants compensation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Compensated {
    /// Running sum.
    pub sum: f64,
    /// Running compensation (lost low-order bits).
    pub comp: f64,
}

impl Compensated {
    /// The zero partial.
    #[must_use]
    pub const fn zero() -> Self {
        Compensated { sum: 0.0, comp: 0.0 }
    }

    /// Accumulate one term (Neumaier's variant: compensation also captures
    /// the case where the term dominates the running sum).
    #[must_use]
    pub fn accumulate(self, x: f64) -> Self {
        let t = self.sum + x;
        let comp = if self.sum.abs() >= x.abs() {
            self.comp + ((self.sum - t) + x)
        } else {
            self.comp + ((x - t) + self.sum)
        };
        Compensated { sum: t, comp }
    }

    /// The compensated total.
    #[must_use]
    pub fn value(self) -> f64 {
        self.sum + self.comp
    }
}

impl Reduce for Compensated {
    fn identity() -> Self {
        Compensated::zero()
    }

    fn merge(self, other: Self) -> Self {
        // Two-sum of the partial sums; compensations add exactly enough
        // that `value()` of the merge equals the compensated total.
        let s = self.sum + other.sum;
        let bb = s - self.sum;
        let err = (self.sum - (s - bb)) + (other.sum - bb);
        Compensated { sum: s, comp: self.comp + other.comp + err }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Concatenation partial: exposes the visit order of the tree.
    struct Concat(Vec<u64>);

    impl Reduce for Concat {
        fn identity() -> Self {
            Concat(Vec::new())
        }
        fn merge(mut self, mut other: Self) -> Self {
            self.0.append(&mut other.0);
            self
        }
    }

    #[test]
    fn tree_shape_is_a_pure_function_of_length() {
        // Concatenation exposes the visit order: it must be ascending for
        // every length, i.e. the tree only regroups, never reorders.
        for n in 0..40usize {
            let items: Vec<Concat> = (0..n as u64).map(|i| Concat(vec![i])).collect();
            let folded = pairwise_fold(items);
            assert!(folded.0.iter().copied().eq(0..n as u64), "n={n}");
        }
        assert_eq!(largest_pow2_below(2), 1);
        assert_eq!(largest_pow2_below(3), 2);
        assert_eq!(largest_pow2_below(8), 4);
        assert_eq!(largest_pow2_below(9), 8);
    }

    #[test]
    fn compensated_partials_merge_to_the_compensated_total() {
        let xs: Vec<f64> = (0..1000).map(|i| 1.0 / f64::from(i + 1)).collect();
        let whole = xs.iter().fold(Compensated::zero(), |c, &x| c.accumulate(x));
        let a = xs[..500].iter().fold(Compensated::zero(), |c, &x| c.accumulate(x));
        let b = xs[500..].iter().fold(Compensated::zero(), |c, &x| c.accumulate(x));
        let merged = a.merge(b).value();
        assert!((whole.value() - merged).abs() <= 1e-12 * whole.value().abs());
    }
}
