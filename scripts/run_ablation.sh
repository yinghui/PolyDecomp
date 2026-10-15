#!/usr/bin/env bash
# SPDX-License-Identifier: BSD-3-Clause
# Compares the cut rules of hole elimination (pd::CutRule, pd_decompose --cut-rule) with the
# exact pipeline and the certificate (--verify) on the generated families written by
# run_generated.sh and on random_holes with 64, 256 and 1024 holes (seeds 1 to 5).
# Stats go to BUILD_DIR/ablation/; the summary goes to stdout. Fails if any result is not
# certified.
#   usage: scripts/run_ablation.sh [BUILD_DIR] [CPU]
set -euo pipefail
b="${1:-build/gcc-release}"
cpu="${2:-2}"
d="${b}/ablation"
G="${b}/generated"
mkdir -p "${d}"
for h in 64 256 1024; do
  for s in 1 2 3 4 5; do "${b}/pd_gen" --random-holes "${h}" --seed "${s}"; done > "${d}/holes-${h}.wkt"
done
for rule in farthest closest boundary horizontal; do
  for f in row star c cluster grid; do
    taskset -c "${cpu}" "${b}/pd_decompose" --quiet --verify --cut-rule "${rule}" --wkt "${G}/${f}.wkt" \
      --out /dev/null --stats "${d}/${rule}.${f}.jsonl"
  done
  for h in 64 256 1024; do
    taskset -c "${cpu}" "${b}/pd_decompose" --quiet --verify --cut-rule "${rule}" --wkt "${d}/holes-${h}.wkt" \
      --out /dev/null --stats "${d}/${rule}.holes-${h}.jsonl"
  done
done
python3 "$(dirname "$0")/summarize_ablation.py" "${d}"
