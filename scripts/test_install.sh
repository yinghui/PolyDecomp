#!/usr/bin/env bash
# SPDX-License-Identifier: BSD-3-Clause
# Installs a build into a temporary prefix, then configures, builds and runs an
# external project (examples/consumer) against it with find_package(polydecomp).
#   usage: scripts/test_install.sh [BUILD_DIR]   (default build/gcc-release)
set -euo pipefail
cd "$(dirname "$0")/.."
build="${1:-build/gcc-release}"
tmp=$(mktemp -d)
trap 'rm -rf "${tmp}"' EXIT
cmake --install "${build}" --prefix "${tmp}/prefix" > /dev/null
cmake -S examples/consumer -B "${tmp}/consumer" -G Ninja -DCMAKE_PREFIX_PATH="${tmp}/prefix" > /dev/null
cmake --build "${tmp}/consumer" > /dev/null
"${tmp}/consumer/consumer"
echo "install test passed"
