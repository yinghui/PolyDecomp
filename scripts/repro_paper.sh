#!/usr/bin/env bash
# SPDX-License-Identifier: BSD-3-Clause
# Regenerates every result reported in the paper, in build/paper-results/:
#   1. unit tests on all presets              (scripts/check.sh)
#   2. legacy input validation                (VALIDATION.txt)
#   3. 5 x 1000 generated inputs, exact pipeline + exact certificate (GENERATED.txt)
#   4. the same outputs through the independent GEOS oracle          (ORACLE.txt)
#   5. floating-point census: double and naive vs exact              (CENSUS.txt)
#   6. double-kernel outputs through the GEOS oracle                 (ORACLE_DOUBLE.txt)
#   7. double kernel with OnFailure::RetryExact on the census inputs  (RETRY.txt)
#   8. unit-test case count                                          (TESTS.txt)
#      and the hand-made degenerate set, exact + certificate         (DEGENERATE.txt)
#   9. benchmarks against CGAL                                       (BENCH.csv, BENCH.md)
# Needs: GMP, CGAL (libcgal-dev), and the oracle venv (scripts/oracle/setup_venv.sh).
# Timings depend on the machine; the script records its CPU and frequency policy.
set -euo pipefail
cd "$(dirname "$0")/.."
out=build/paper-results
mkdir -p "${out}"
py="${PD_ORACLE_VENV:-$HOME/.venvs/pd-oracle}/bin/python"
[ -x "${py}" ] || scripts/oracle/setup_venv.sh

{
  echo "date: $(date -Iseconds)"
  echo "host: $(uname -srm)"
  lscpu | grep -E 'Model name|Socket|Core|Thread|MHz' || true
  echo "governor: $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo unknown)"
  g++ --version | head -1
} > "${out}/MACHINE.txt"

scripts/check.sh
cmake --preset gcc-bench > /dev/null && cmake --build --preset gcc-bench

scripts/validate_legacy.sh build/gcc-release
cp data/legacy/VALIDATION.txt "${out}/"

scripts/run_generated.sh 1000 build/gcc-release > "${out}/GENERATED.txt"
G=build/gcc-release/generated
for f in row star c cluster grid; do
  "${py}" scripts/oracle/oracle.py "${G}/${f}.wkt" "${G}/${f}.pieces.wkt" --convex | sed "s/^/${f}: /"
done > "${out}/ORACLE.txt"

scripts/run_census.sh 1000 300 build/gcc-release > "${out}/CENSUS.txt"
C=build/gcc-release/census
for f in row star c cluster grid nearcollinear neartouch snapped; do
  build/gcc-release/pd_decompose --quiet --kernel double --wkt "${C}/${f}.wkt" --out "${C}/${f}.double.pieces.wkt" || true
  "${py}" scripts/oracle/oracle.py "${C}/${f}.wkt" "${C}/${f}.double.pieces.wkt" --convex | sed "s/^/${f}: /"
done > "${out}/ORACLE_DOUBLE.txt"

retry=()
for f in row star c cluster grid nearcollinear neartouch snapped; do
  build/gcc-release/pd_decompose --quiet --kernel double --retry-exact --wkt "${C}/${f}.wkt" --out "${C}/${f}.retry.pieces.wkt" --stats "${C}/${f}.retry.jsonl" || true
  retry+=("${f}=${C}/${f}.retry.jsonl")
done
python3 scripts/summarize_retry.py "${retry[@]}" > "${out}/RETRY.txt"
build/gcc-release/tests/pd_tests 2>&1 | grep 'test cases' > "${out}/TESTS.txt"
build/gcc-release/pd_gen --degenerate --out build/gcc-release/degenerate.wkt
build/gcc-release/pd_decompose --verify --wkt build/gcc-release/degenerate.wkt \
  --out build/gcc-release/degenerate.pieces.wkt --stats "${out}/DEGENERATE.jsonl" > "${out}/DEGENERATE.txt" 2>&1

taskset -c 2 build/gcc-bench/pd_bench --repeats 5 --csv "${out}/BENCH.csv"
python3 scripts/summarize_bench.py "${out}/BENCH.csv" > "${out}/BENCH.md"
echo "results in ${out}/"
