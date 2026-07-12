#!/bin/sh
# Phase G4 task farm: 23 independent relaxations, 4 OS processes at a time
# (the survey's measured 3.96x multiprocess scaling route; per-run
# bit-identity intact — parallelism is process scheduling only).
# Order: single-knot reference + attractive channel FIRST (they carry the
# bond equation), then aligned, then repulsive.
# Usage: sh farm.sh [outdir] [maxit]
set -e
cd "$(dirname "$0")"
OUT="${1:-runs}"
MAXIT="${2:-400}"
BIN=./target/release/twoknot_run
mkdir -p "$OUT"

run4() {
  for spec in "$@"; do
    ch="${spec%%:*}"
    m="${spec##*:}"
    if [ "$ch" = single ]; then log="$OUT/log_single.txt"; else log="$OUT/log_${ch}_m${m}.txt"; fi
    "$BIN" "$ch" "$m" "$MAXIT" "$OUT" > "$log" 2>&1 &
  done
  wait
}

echo "batch 1/6 (single + attract 14,16,18)"; run4 single:0 attract:14 attract:16 attract:18
echo "batch 2/6 (attract 20,23,26,29)";       run4 attract:20 attract:23 attract:26 attract:29
echo "batch 3/6 (attract 38 + align 14,16,18)"; run4 attract:38 align:14 align:16 align:18
echo "batch 4/6 (align 20,23,26,29)";         run4 align:20 align:23 align:26 align:29
echo "batch 5/6 (repulse 14,16,18,20)";       run4 repulse:14 repulse:16 repulse:18 repulse:20
echo "batch 6/6 (repulse 23,26,29)";          run4 repulse:23 repulse:26 repulse:29
echo "farm complete"
