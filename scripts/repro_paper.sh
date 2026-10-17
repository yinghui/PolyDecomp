#!/usr/bin/env bash
# SPDX-License-Identifier: BSD-3-Clause
# Regenerates every result reported in the paper, in build/paper-results/:
#   1. unit tests on all presets              (scripts/check.sh)
#   2. legacy input validation                (VALIDATION.txt)
#   3. 5 x 1000 generated inputs, exact hole elimination + exact certificate (GENERATED.txt)
#   4. the same hole-free pieces through the independent GEOS oracle  (ORACLE.txt)
#   5. floating-point census of hole elimination: double and naive vs exact (CENSUS.txt)
#   6. double-kernel hole-free pieces through the GEOS oracle         (ORACLE_DOUBLE.txt)
#   7. double kernel with OnFailure::RetryExact on the census inputs  (RETRY.txt)
#   8. unit-test case count                                           (TESTS.txt)
#      and the hand-made degenerate set, hole elimination + certificate (DEGENERATE.txt)
#   9. scaling of hole elimination: random holes up to 2^14, cheese up to 2^12 holes, against
#      PolyPartition's and CGAL's bridging and the triangulations of Earcut and CGAL
#                                                                     (HOLE_SCALE.csv, .md)
#  10. cut-rule ablation, hole elimination + certificate               (ABLATION.txt)
#  11. cutting along lines: 200 cuts per sequence, certificate after every cut, with the
#      Clipper2, CGAL Boolean and CGAL arrangement baselines          (SLICE.csv, SLICE.md)
#      and Clipper2 at five integer scales                            (CLIPPER_SCALE.csv)
#  12. strokes on moving pieces: 300 strokes on cheese with 64 and 256 holes (STROKES.csv)
#  13. state sizes after the last cut of every sequence of step 11, ours against the CGAL
#      arrangement (deterministic, no timing)                          (SLICE_SIZE.csv)
#  14. triangulations of the census inputs (Earcut, CGAL's constrained triangulation) under the
#      exact certificate (no timing)                                  (TRIANGULATE.csv)
# Needs: GMP, CGAL (libcgal-dev), and the oracle venv (scripts/oracle/setup_venv.sh).
# Timings depend on the machine. The benchmarks run on CPU ${PD_BENCH_CPU:-2}; the
# script records its frequency policy in MACHINE.txt, samples its average frequency
# every half second while each benchmark runs (mean and minimum in MACHINE.txt),
# and warns unless its governor is "performance" (a fixed frequency needs root, for example
#   echo performance | sudo tee /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor
# with turbo disabled).
set -euo pipefail
cd "$(dirname "$0")/.."
out=build/paper-results
mkdir -p "${out}"
# Results of an earlier run (also of older versions of this script) must not survive into this one.
rm -f "${out}"/*.txt "${out}"/*.csv "${out}"/*.md "${out}"/*.jsonl
py="${PD_ORACLE_VENV:-$HOME/.venvs/pd-oracle}/bin/python"
cpu="${PD_BENCH_CPU:-2}"
freq=/sys/devices/system/cpu/cpu${cpu}/cpufreq
[ -x "${py}" ] || scripts/oracle/setup_venv.sh

# bench NAME COMMAND...: runs a benchmark while another CPU samples the benchmark CPU's average
# frequency (cpuinfo_avg_freq, kHz) every half second; appends the mean and the minimum to MACHINE.txt.
bench() {
  local name=$1 status=0 samples sampler other=0
  shift
  [ "${cpu}" = 0 ] && other=1
  samples=$(mktemp)
  taskset -c "${other}" sh -c 'while :; do cat "$1" 2>/dev/null || echo 0; sleep 0.5; done' _ \
    "${freq}/cpuinfo_avg_freq" > "${samples}" &
  sampler=$!
  "$@" || status=$?
  kill "${sampler}" 2>/dev/null || true
  wait "${sampler}" 2>/dev/null || true
  awk -v n="${name}" '{ s += $1; c++; if (c == 1 || $1 < m) m = $1 }
    END { printf "bench kHz during %s: mean %d min %d samples %d\n", n, (c ? s / c : 0), m, c }' \
    "${samples}" >> "${out}/MACHINE.txt"
  rm -f "${samples}"
  return "${status}"
}

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

scripts/run_generated.sh 1000 build/gcc-release holes > "${out}/GENERATED.txt"
G=build/gcc-release/generated
for f in row star c cluster grid; do
  "${py}" scripts/oracle/oracle.py "${G}/${f}.wkt" "${G}/${f}.pieces.wkt" --rational "${G}/${f}.pieces.rat" | sed "s/^/${f}: /"
done > "${out}/ORACLE.txt"

scripts/run_census.sh 1000 300 build/gcc-release holes > "${out}/CENSUS.txt"
C=build/gcc-release/census
for f in row star c cluster grid nearcollinear neartouch snapped; do
  build/gcc-release/pd_decompose --quiet --stage holes --kernel double --wkt "${C}/${f}.wkt" --out "${C}/${f}.double.pieces.wkt" || true
  "${py}" scripts/oracle/oracle.py "${C}/${f}.wkt" "${C}/${f}.double.pieces.wkt" | sed "s/^/${f}: /"
done > "${out}/ORACLE_DOUBLE.txt"

retry=()
for f in row star c cluster grid nearcollinear neartouch snapped; do
  build/gcc-release/pd_decompose --quiet --stage holes --kernel double --retry-exact --wkt "${C}/${f}.wkt" --out "${C}/${f}.retry.pieces.wkt" --stats "${C}/${f}.retry.jsonl" || true
  retry+=("${f}=${C}/${f}.retry.jsonl")
done
python3 scripts/summarize_retry.py "${retry[@]}" > "${out}/RETRY.txt"
build/gcc-release/tests/pd_tests 2>&1 | grep 'test cases' > "${out}/TESTS.txt"
build/gcc-release/pd_gen --degenerate --out build/gcc-release/degenerate.wkt
build/gcc-release/pd_decompose --verify --stage holes --wkt build/gcc-release/degenerate.wkt \
  --out build/gcc-release/degenerate.pieces.wkt --stats "${out}/DEGENERATE.jsonl" > "${out}/DEGENERATE.txt" 2>&1

bench HOLESCALE taskset -c "${cpu}" build/gcc-bench/pd_bench --suite hole-scale --repeats 3 --csv "${out}/HOLE_SCALE.csv"
python3 scripts/summarize_bench.py "${out}/HOLE_SCALE.csv" > "${out}/HOLE_SCALE.md"
scripts/run_ablation.sh build/gcc-release "${cpu}" holes > "${out}/ABLATION.txt"
bench SLICE taskset -c "${cpu}" build/gcc-bench/pd_bench --suite slice --cuts 200 --csv "${out}/SLICE.csv"
python3 scripts/summarize_slice.py "${out}/SLICE.csv" > "${out}/SLICE.md"
taskset -c "${cpu}" build/gcc-bench/pd_bench --suite clipper-scale --csv "${out}/CLIPPER_SCALE.csv"
bench STROKES taskset -c "${cpu}" build/gcc-bench/pd_bench --suite strokes --cuts 300 --csv "${out}/STROKES.csv"
build/gcc-bench/pd_bench --suite slice-size --cuts 200 --csv "${out}/SLICE_SIZE.csv"
tri=()
for f in row star c cluster grid nearcollinear neartouch snapped; do tri+=(--wkt "${C}/${f}.wkt"); done
build/gcc-bench/pd_bench --suite triangulate "${tri[@]}" --csv "${out}/TRIANGULATE.csv"
echo "results in ${out}/"
