#!/bin/sh
# Phase G5b refine farm: ATTRACTIVE channel only, well window + far anchor +
# single-knot reference, at a 4x iteration cap (default 1600).  G4's runs/
# stay untouched; output goes to runs1600/ by default.
# Usage: sh farm_refine.sh [outdir] [maxit]
set -e
cd "$(dirname "$0")"
OUT="${1:-runs1600}"
MAXIT="${2:-1600}"
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

echo "batch 1/2 (single + attract 16,18,20)"; run4 single:0 attract:16 attract:18 attract:20
echo "batch 2/2 (attract 23,38)";             run4 attract:23 attract:38
echo "refine farm complete"
