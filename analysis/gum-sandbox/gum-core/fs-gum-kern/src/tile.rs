//! The one parallel primitive: a tile-decomposed sweep whose serial and
//! threaded executions are bit-identical BY CONSTRUCTION.
//!
//! CANONICAL TILE DECOMPOSITION (part of the bit contract, the fs-la KC
//! precedent): an extent of `m` planes along the slow (i) axis is cut
//! into i-slabs of [`TILE_I`](crate::TILE_I) = 4 planes; the tile count
//! `ceil(m / TILE_I)` is a pure function of `m`, NEVER of thread count.
//! Per tile the caller accumulates serially in ascending (i, j, k) order;
//! tile partials are combined up the fixed largest-pow2-below pairwise
//! tree ([`crate::reduce::pairwise_fold`]) in ascending tile order.
//! Scheduling and steal order are provably irrelevant: every rounding
//! decision is a function of (m, TILE_I) alone.
//!
//! THREADED PATH (the fs-la band-dispenser pattern, no fs-exec needed):
//! `std::thread::scope` with `workers = threads.min(tiles)` (excess
//! spawns measured 2-9x slower on the ts1 box — gemm.rs comment),
//! an `AtomicUsize` tile dispenser (work stealing absorbs heterogeneous
//! /throttled cores), and cache-padded per-tile result slots so workers
//! never false-share a line.  A tile's content is a pure function of the
//! tile, never of which thread computed it or in what order — the lock
//! on each slot guards HANDOFF only (taken exactly once per tile).

use std::ops::Range;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

use crate::reduce::{pairwise_fold, Reduce};
use crate::TILE_I;

/// Hand-rolled cache-line padding for the per-tile slots (64 B: x86-64
/// line; fs-alloc uses a 128 B superset but is not in this crate's
/// dependency closure — padding is bit-neutral, only contention-relevant).
#[repr(align(64))]
pub struct CachePadded<T>(pub T);

/// Number of i-slab tiles for an extent of `m` planes: `ceil(m / TILE_I)`.
/// A pure function of `m` — the property the bit contract rests on.
#[must_use]
pub fn tile_count(m: usize) -> usize {
    m.div_ceil(TILE_I)
}

/// Plane range of tile `t` (ascending, disjoint, covering `0..m`).
#[must_use]
pub fn tile_range(t: usize, m: usize) -> Range<usize> {
    let lo = t * TILE_I;
    lo..m.min(lo + TILE_I)
}

/// All tile ranges of extent `m`, in tile order.
#[must_use]
pub fn tile_ranges(m: usize) -> Vec<Range<usize>> {
    (0..tile_count(m)).map(|t| tile_range(t, m)).collect()
}

/// Run one job per tile and return the per-tile results IN TILE ORDER.
///
/// `threads <= 1` (or a single job) runs serially in tile order; more
/// threads run the dispenser/slot machinery described in the module docs.
/// Both paths call the same `work` closure on the same `(t, job)` pairs,
/// so per-tile results are identical; only the completion ORDER differs,
/// and the returned Vec re-establishes tile order — this is the
/// bit-identity-by-construction argument, in one place.
pub fn tiled_run<J, P, F>(jobs: Vec<J>, threads: usize, work: F) -> Vec<P>
where
    J: Send,
    P: Send,
    F: Fn(usize, J) -> P + Sync,
{
    let nt = jobs.len();
    let workers = threads.max(1).min(nt);
    if workers <= 1 {
        return jobs.into_iter().enumerate().map(|(t, j)| work(t, j)).collect();
    }
    let job_cells: Vec<Mutex<Option<J>>> = jobs.into_iter().map(|j| Mutex::new(Some(j))).collect();
    let slots: Vec<CachePadded<Mutex<Option<P>>>> =
        (0..nt).map(|_| CachePadded(Mutex::new(None))).collect();
    let next = AtomicUsize::new(0);
    std::thread::scope(|s| {
        for _ in 0..workers {
            s.spawn(|| loop {
                let t = next.fetch_add(1, Ordering::Relaxed);
                if t >= nt {
                    break;
                }
                let job = job_cells[t]
                    .lock()
                    .expect("job cell")
                    .take()
                    .expect("fetch_add dispenses each tile exactly once");
                let p = work(t, job);
                *slots[t].0.lock().expect("slot") = Some(p);
            });
        }
    });
    slots
        .into_iter()
        .map(|c| c.0.into_inner().expect("slot").expect("every tile completed"))
        .collect()
}

/// Tiled reduction over an extent of `m` planes: per-tile partials from
/// `per_tile` (which must accumulate in ascending (i, j, k) order within
/// its plane range — the caller-side half of the bit contract), combined
/// up the fixed pairwise tree in tile order.
pub fn sweep_reduce<P, F>(m: usize, threads: usize, per_tile: F) -> P
where
    P: Reduce + Send,
    F: Fn(Range<usize>) -> P + Sync,
{
    pairwise_fold(tiled_run(tile_ranges(m), threads, |_t, r| per_tile(r)))
}

/// Split a component-major buffer (`comps` equal components of
/// `comp_len = data.len() / comps` elements, each an array of planes of
/// `unit` elements) into per-tile jobs: `jobs[t][a]` is the mutable
/// sub-slice of component `a` covering planes `bounds[t]` (bounds are in
/// plane units, ascending and disjoint; planes outside every bound —
/// e.g. ghost planes — are simply not handed out).  Pure safe splitting:
/// the disjointness the threaded map needs is established here once, by
/// construction, with no unsafe anywhere.
#[must_use]
pub fn split_tiles<'a>(
    mut data: &'a mut [f64],
    comps: usize,
    unit: usize,
    bounds: &[Range<usize>],
) -> Vec<Vec<&'a mut [f64]>> {
    let comp_len = data.len() / comps;
    let mut jobs: Vec<Vec<&'a mut [f64]>> =
        (0..bounds.len()).map(|_| Vec::with_capacity(comps)).collect();
    for _a in 0..comps {
        let (comp, rest) = data.split_at_mut(comp_len);
        data = rest;
        let mut comp = comp;
        let mut pos = 0usize;
        for (t, b) in bounds.iter().enumerate() {
            let (start, end) = (b.start * unit, b.end * unit);
            debug_assert!(start >= pos && end <= pos + comp.len());
            let (_gap, tail) = comp.split_at_mut(start - pos);
            let (chunk, rest2) = tail.split_at_mut(end - start);
            comp = rest2;
            pos = end;
            jobs[t].push(chunk);
        }
    }
    jobs
}

/// f64 max partial for tiled max-reductions (drift-max).  Max is
/// order-insensitive over a fixed value set, but it rides the same tile
/// protocol as everything else (roadmap G1: "keep the tile protocol
/// anyway").
#[derive(Clone, Copy, Debug)]
pub struct MaxP(pub f64);

impl Reduce for MaxP {
    fn identity() -> Self {
        MaxP(0.0)
    }
    fn merge(self, other: Self) -> Self {
        MaxP(self.0.max(other.0))
    }
}

/// Plain f64 sum partial (per-tile plain accumulation — current physics
/// semantics; see `reduce::Compensated` for the compensated alternative).
#[derive(Clone, Copy, Debug)]
pub struct SumP(pub f64);

impl Reduce for SumP {
    fn identity() -> Self {
        SumP(0.0)
    }
    fn merge(self, other: Self) -> Self {
        SumP(self.0 + other.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tiles_cover_the_extent_disjointly() {
        for m in 1..40usize {
            let rs = tile_ranges(m);
            assert_eq!(rs.len(), m.div_ceil(TILE_I));
            let mut next = 0usize;
            for r in &rs {
                assert_eq!(r.start, next);
                assert!(r.end > r.start && r.end - r.start <= TILE_I);
                next = r.end;
            }
            assert_eq!(next, m);
        }
    }

    #[test]
    fn tiled_run_is_thread_count_invariant_and_ordered() {
        let jobs: Vec<usize> = (0..13).collect();
        let serial = tiled_run(jobs.clone(), 1, |t, j| (t, j * j));
        for threads in 2..=4 {
            let par = tiled_run(jobs.clone(), threads, |t, j| (t, j * j));
            assert_eq!(serial, par, "threads={threads}");
        }
        assert!(serial.iter().enumerate().all(|(i, &(t, _))| i == t));
    }

    #[test]
    fn sweep_reduce_matches_serial_tile_fold_bitwise() {
        // Ill-conditioned per-plane values: any order change would move bits.
        let val = |i: usize| if i % 3 == 0 { 1.0e16 } else { -0.3 * i as f64 };
        let per_tile = |r: Range<usize>| {
            let mut s = 0.0f64;
            for i in r {
                s += val(i);
            }
            SumP(s)
        };
        for m in [1usize, 4, 5, 17, 33, 97] {
            let s1 = sweep_reduce(m, 1, per_tile).0;
            for threads in 2..=4 {
                let st = sweep_reduce(m, threads, per_tile).0;
                assert_eq!(s1.to_bits(), st.to_bits(), "m={m} threads={threads}");
            }
        }
    }

    #[test]
    fn split_tiles_hands_out_disjoint_component_chunks() {
        let n = 6usize; // 6 planes of 4 elements, 2 components, tiles of 4
        let mut data = vec![0.0f64; 2 * n * 4];
        let bounds = tile_ranges(n);
        let jobs = split_tiles(&mut data, 2, 4, &bounds);
        assert_eq!(jobs.len(), 2);
        assert_eq!(jobs[0].len(), 2);
        assert_eq!(jobs[0][0].len(), 16);
        assert_eq!(jobs[1][0].len(), 8);
        // Write through every chunk, then check full coverage.
        let mut jobs = jobs;
        for (t, job) in jobs.iter_mut().enumerate() {
            for (a, chunk) in job.iter_mut().enumerate() {
                for v in chunk.iter_mut() {
                    *v = (t * 10 + a) as f64 + 1.0;
                }
            }
        }
        assert!(data.iter().all(|&v| v != 0.0));
    }
}
