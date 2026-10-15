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
#   9. benchmarks against CGAL and PolyPartition                     (BENCH.csv, BENCH.md)
#  10. scaling: stars and combs up to 2^16 vertices, up to 2^14 holes  (SCALE.csv, SCALE.md)
#  11. cut-rule ablation, exact + certificate                         (ABLATION.txt)
# Needs: GMP, CGAL (libcgal-dev), and the oracle venv (scripts/oracle/setup_venv.sh).
# Timings depend on the machine. The benchmarks run on CPU ${PD_BENCH_CPU:-2}; the
# script records its frequency policy in MACHINE.txt and warns unless its governor
# is "performance" (a fixed frequency needs root, for example
#   echo performance | sudo tee /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor
# with turbo disabled).
set -euo pipefail
cd "$(dirname "$0")/.."
out=build/paper-results
mkdir -p "${out}"
py="${PD_ORACLE_VENV:-$HOME/.venvs/pd-oracle}/bin/python"
cpu="${PD_BENCH_CPU:-2}"
freq=/sys/devices/system/cpu/cpu${cpu}/cpufreq
[ -x "${py}" ] || scripts/oracle/setup_venv.sh

{
  echo "date: $(date -Iseconds)"
  echo "host: $(uname -srm)"
  lscpu | grep -E 'Model name|Socket|Core|Thread|MHz' || true
  echo "governor: $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo unknown)"
  echo "bench cpu: ${cpu}"
  echo "bench governor: $(cat ${freq}/scaling_governor 2>/dev/null || echo unknown)"
  echo "bench scaling min/max kHz: $(cat ${freq}/scaling_min_freq 2>/dev/null || echo ?) $(cat ${freq}/scaling_max_freq 2>/dev/null || echo ?)"
  echo "no_turbo: $(cat /sys/devices/system/cpu/intel_pstate/no_turbo 2>/dev/null || echo unknown)"
  g++ --version | head -1
} > "${out}/MACHINE.txt"
if [ "$(cat ${freq}/scaling_governor 2>/dev/null)" != performance ]; then
  echo "warning: CPU ${cpu} governor is not 'performance'; timings will vary with frequency scaling" >&2
fi

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

taskset -c "${cpu}" build/gcc-bench/pd_bench --repeats 5 --csv "${out}/BENCH.csv"
echo "bench avg kHz after BENCH: $(cat ${freq}/cpuinfo_avg_freq 2>/dev/null || echo ?)" >> "${out}/MACHINE.txt"
python3 scripts/summarize_bench.py "${out}/BENCH.csv" > "${out}/BENCH.md"
taskset -c "${cpu}" build/gcc-bench/pd_bench --suite scale --repeats 3 --csv "${out}/SCALE.csv"
echo "bench avg kHz after SCALE: $(cat ${freq}/cpuinfo_avg_freq 2>/dev/null || echo ?)" >> "${out}/MACHINE.txt"
python3 scripts/summarize_bench.py "${out}/SCALE.csv" > "${out}/SCALE.md"
scripts/run_ablation.sh build/gcc-release "${cpu}" > "${out}/ABLATION.txt"
echo "results in ${out}/"
