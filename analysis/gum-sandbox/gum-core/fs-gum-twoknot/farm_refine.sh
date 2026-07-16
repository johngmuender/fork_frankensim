#!/bin/sh
# Phase G5b refine farm: ATTRACTIVE channel only, well window + far anchor +
# single-knot reference, at a 3x iteration cap (default 1200; the 1600 plan
# was cut to 1200 after a container restart killed the first farm and G5a's
# 4-thread n192_fr5 took the box — max 3 concurrent processes here).
# G4's runs/ stay untouched; output goes to runs1200/ by default.
# Banking order: the anchor pair (single + far) and the well point (m=18)
# run FIRST so the key observable survives another restart.
# Usage: sh farm_refine.sh [outdir] [maxit]
set -e
cd "$(dirname "$0")"
OUT="${1:-runs1200}"
MAXIT="${2:-1200}"
BIN=./target/release/twoknot_run
mkdir -p "$OUT"

runN() {
  for spec in "$@"; do
    ch="${spec%%:*}"
    m="${spec##*:}"
    if [ "$ch" = single ]; then log="$OUT/log_single.txt"; else log="$OUT/log_${ch}_m${m}.txt"; fi
    "$BIN" "$ch" "$m" "$MAXIT" "$OUT" > "$log" 2>&1 &
  done
  wait
}

echo "batch 1/2 (single + attract 38,18 — bank the anchor pair + well first)"
runN single:0 attract:38 attract:18
echo "batch 2/2 (attract 16,20,23)"
runN attract:16 attract:20 attract:23
echo "refine farm complete"
