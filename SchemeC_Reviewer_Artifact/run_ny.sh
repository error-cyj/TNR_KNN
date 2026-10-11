#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PHASE="${1:-density}"
case "$PHASE" in all|density|k) ;; *) echo "Usage: bash run_ny.sh [density|k|all]" >&2; exit 2;; esac
for f in "$ROOT/data/NY/Ein.txt" "$ROOT/Object data/NY/data_0.001.txt" "$ROOT/Queary data/NY/_k20_Query.txt"; do
    if [[ ! -s "$f" ]]; then echo "Missing required NY experimental input: $f" >&2; exit 2; fi
done
if [[ "$PHASE" == "density" || "$PHASE" == "all" ]]; then
    for d in 0.1 0.01 0.001 0.0001; do
        [[ -s "$ROOT/Object data/NY/data_${d}.txt" ]] || { echo "Missing object file for density $d" >&2; exit 2; }
    done
fi
if [[ "$PHASE" == "k" || "$PHASE" == "all" ]]; then
    for k in 10 20 30 40; do
        [[ -s "$ROOT/Queary data/NY/_k${k}_Query.txt" ]] || { echo "Missing query file for k=$k" >&2; exit 2; }
    done
fi
if [[ ! -x "$ROOT/build/v89_final_benchmark" ]]; then
    echo "Build first: cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j2" >&2
    exit 2
fi
exec "$ROOT/build/v89_final_benchmark" --base "$ROOT" --output-dir "$ROOT/results/ny_${PHASE}" --phase "$PHASE" --datasets NY --warmup 100
