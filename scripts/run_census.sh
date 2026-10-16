#!/usr/bin/env bash
# SPDX-License-Identifier: BSD-3-Clause
# Floating-point robustness census: for every generator family, runs the exact
# pipeline and a floating-point one (the robust double kernel, and the naive
# unprotected baseline) on the same inputs and classifies each run.
#   usage: scripts/run_census.sh [COUNT_REGULAR] [COUNT_ADVERSARIAL] [BUILD_DIR] [STAGE]
#   (defaults: 1000, 300, build/gcc-release, both; STAGE holes compares hole elimination only;
#   output in BUILD_DIR/census/)
set -euo pipefail
cd "$(dirname "$0")/.."
nreg="${1:-1000}"
nadv="${2:-300}"
build="${3:-build/gcc-release}"
stage="${4:-both}"
out="${build}/census"
mkdir -p "${out}"

regular=(row star c cluster grid)
adversarial=(nearcollinear neartouch snapped)
for f in "${regular[@]}"; do "${build}/pd_gen" --family "${f}" --count "${nreg}" --seed 1 --out "${out}/${f}.wkt"; done
for f in "${adversarial[@]}"; do "${build}/pd_gen" --family "${f}" --count "${nadv}" --seed 1 --out "${out}/${f}.wkt"; done

pids=()
for f in "${regular[@]}" "${adversarial[@]}"; do
  for k in double naive; do
    "${build}/pd_decompose" --quiet --stage "${stage}" --wkt "${out}/${f}.wkt" --compare "${k}" --stats "${out}/${f}.${k}.jsonl" &
    pids+=($!)
  done
done
status=0
for p in "${pids[@]}"; do wait "${p}" || status=1; done

{
  echo "## robust double kernel vs exact"
  python3 scripts/summarize_census.py "${out}"/{row,star,c,cluster,grid,nearcollinear,neartouch,snapped}.double.jsonl
  echo
  echo "## naive double (no protections) vs exact"
  python3 scripts/summarize_census.py "${out}"/{row,star,c,cluster,grid,nearcollinear,neartouch,snapped}.naive.jsonl
} | tee "${out}/CENSUS.txt"
exit "${status}"
