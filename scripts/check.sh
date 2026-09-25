#!/usr/bin/env bash
# SPDX-License-Identifier: BSD-3-Clause
# Configure, build and test polydecomp on every preset (or the ones given as
# arguments). Stops at the first failure.
set -euo pipefail
cd "$(dirname "$0")/.."

presets=("$@")
if [ ${#presets[@]} -eq 0 ]; then
  presets=(gcc-debug gcc-release clang-debug clang-release gcc-asan clang-asan gcc-nogmp)
fi

# Some systems ship a newer GCC than its libstdc++ development files, so clang
# cannot link. If so, point clang at the newest GCC toolchain that works
# (override by setting PD_CLANG_EXTRA_FLAGS yourself).
if [ -z "${PD_CLANG_EXTRA_FLAGS+x}" ] && command -v clang++ > /dev/null; then
  tmp=$(mktemp -d)
  printf 'int main(){}\n' > "${tmp}/t.cpp"
  if ! clang++ "${tmp}/t.cpp" -o "${tmp}/t" 2> /dev/null; then
    for d in $(ls -d /usr/lib/gcc/*/* 2> /dev/null | sort -V -r); do
      if clang++ --gcc-install-dir="${d}" "${tmp}/t.cpp" -o "${tmp}/t" 2> /dev/null; then
        export PD_CLANG_EXTRA_FLAGS="--gcc-install-dir=${d}"
        echo "note: clang cannot link with its default GCC toolchain; using ${PD_CLANG_EXTRA_FLAGS}"
        break
      fi
    done
  fi
  rm -rf "${tmp}"
fi

for p in "${presets[@]}"; do
  echo "=== ${p}"
  cmake --preset "${p}" > /dev/null
  cmake --build --preset "${p}"
  ctest --preset "${p}"
done
if [[ " ${presets[*]} " == *" gcc-release "* ]]; then
  echo "=== install"
  scripts/test_install.sh build/gcc-release
fi
echo "=== all presets passed: ${presets[*]}"
