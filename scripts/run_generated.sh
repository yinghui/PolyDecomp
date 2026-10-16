#!/usr/bin/env bash
# SPDX-License-Identifier: BSD-3-Clause
# Generates COUNT inputs per family (seeds 1..COUNT), runs the exact pipeline on
# each with the exact certificate, and summarises the per-input statistics.
#   usage: scripts/run_generated.sh [COUNT] [BUILD_DIR] [STAGE]
#   (defaults: 1000 per family, build/gcc-release, both; STAGE holes runs hole elimination
#   only; output in BUILD_DIR/generated/, pieces as WKT and as exact rationals)
set -euo pipefail
cd "$(dirname "$0")/.."
count="${1:-1000}"
build="${2:-build/gcc-release}"
stage="${3:-both}"
out="${build}/generated"
mkdir -p "${out}"

families=(row star c cluster grid)
for f in "${families[@]}"; do
  "${build}/pd_gen" --family "${f}" --count "${count}" --seed 1 --out "${out}/${f}.wkt"
done

# One process per family, in parallel.
pids=()
for f in "${families[@]}"; do
  "${build}/pd_decompose" --kernel exact --verify --quiet --stage "${stage}" \
    --wkt "${out}/${f}.wkt" --out "${out}/${f}.pieces.wkt" --rational "${out}/${f}.pieces.rat" \
    --stats "${out}/${f}.jsonl" &
  pids+=($!)
done
status=0
for p in "${pids[@]}"; do wait "${p}" || status=1; done

python3 scripts/summarize_stats.py "${out}"/{row,star,c,cluster,grid}.jsonl | tee "${out}/SUMMARY.txt"
exit "${status}"
