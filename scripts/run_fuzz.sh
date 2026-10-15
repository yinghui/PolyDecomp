#!/usr/bin/env bash
# SPDX-License-Identifier: BSD-3-Clause
# The fuzzing campaign reported in the paper: both libFuzzer targets (clang-fuzz
# preset, ASan + UBSan), each for SECONDS seconds with 8 forked workers, run
# concurrently. Logs: build/clang-fuzz/fuzz_{wkt,grid}.log; the last status line
# of each gives the executions and "oom/timeout/crash" counts.
#   usage: scripts/run_fuzz.sh [SECONDS]   (default 600)
set -euo pipefail
cd "$(dirname "$0")/.."
secs="${1:-600}"
cmake --preset clang-fuzz > /dev/null && cmake --build --preset clang-fuzz
cd build/clang-fuzz
mkdir -p corpus_wkt corpus_grid artifacts
./fuzz_wkt corpus_wkt -fork=8 -max_total_time="${secs}" -max_len=4096 \
  -artifact_prefix=artifacts/wkt- > fuzz_wkt.log 2>&1 &
p1=$!
./fuzz_grid corpus_grid -fork=8 -max_total_time="${secs}" -max_len=512 \
  -artifact_prefix=artifacts/grid- > fuzz_grid.log 2>&1 &
p2=$!
status=0
wait "${p1}" || status=1
wait "${p2}" || status=1
grep -E '^#[0-9]+: cov:' fuzz_wkt.log | tail -1
grep -E '^#[0-9]+: cov:' fuzz_grid.log | tail -1
# In fork mode a crash stops the campaign and prints a report; list any.
if grep -E 'FUZZ FINDING|ERROR: libFuzzer|ERROR: AddressSanitizer|runtime error:' fuzz_wkt.log fuzz_grid.log; then
  status=1
fi
exit "${status}"
