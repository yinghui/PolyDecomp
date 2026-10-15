# polydecomp

A C++17 library that decomposes a polygon with holes into convex pieces, in two stages:

1. **Hole elimination (line cut).** Repeatedly cut along a line derived from the hole geometry
   until every piece is a simple, hole-free polygon.
2. **Convex decomposition (chord split).** Split each simple polygon at a reflex vertex along the
   extension of its incoming edge until every piece is convex.

Both stages run on an exact rational kernel (GMP) or a fast double-precision kernel.

Version 1.2.0. Released under the BSD 3-Clause License (see `LICENSE`); cite it as described in
`CITATION.cff`. Repository: <https://github.com/yinghui/PolyDecomp>.

**Documentation:** the user guide (`docs/USER_GUIDE.md`) covers installation, the input class,
the kernels and their guarantees, the API, the tools and how to reproduce the paper's results.
The API reference is built with `cmake --build --preset gcc-release --target docs` (Doxygen).
**Examples:** `examples/` (each is also run as a test).

## Build

Requirements: CMake ≥ 3.20, Ninja, a C++17 compiler (GCC or Clang), and GMP with its C++
bindings (`libgmp-dev`) for the exact kernel. Without GMP (`-DPD_WITH_GMP=OFF`, preset
`gcc-nogmp`) only the double-precision kernels are built.

```sh
cmake --preset gcc-release
cmake --build --preset gcc-release
ctest --preset gcc-release
```

`scripts/check.sh` builds and tests every preset (GCC and Clang in Debug, Release and
AddressSanitizer/UndefinedBehaviorSanitizer builds, and a build without GMP), then installs the
release build and compiles an external project against it (`scripts/test_install.sh`).

Install with `cmake --install build/gcc-release --prefix <prefix>`; other CMake projects then use
`find_package(polydecomp 1.0 REQUIRED)` and link `pd::polydecomp`.

If Clang fails to link (`cannot find -lstdc++`), it has picked a GCC toolchain whose `libstdc++`
development files are missing. Either install them (e.g. `libstdc++-14-dev`) or set
`PD_CLANG_EXTRA_FLAGS=--gcc-install-dir=/usr/lib/gcc/x86_64-linux-gnu/13` before configuring;
`scripts/check.sh` does the latter automatically.

## Verification and benchmarks

- `scripts/oracle/setup_venv.sh` creates a Python environment with shapely/GEOS;
  `scripts/oracle/oracle.py INPUTS.wkt PIECES.wkt [--convex]` checks results independently of the
  library (valid pieces, inside the input, pairwise disjoint, total area, convexity), using GEOS
  snap-rounded overlay.
- Fuzzing (clang): `cmake --preset clang-fuzz && cmake --build --preset clang-fuzz`, then run
  `build/clang-fuzz/fuzz_wkt` (arbitrary WKT text) or `build/clang-fuzz/fuzz_grid` (polygons on
  a small integer grid). On every valid input the exact result must pass the exact certificate
  and the double pipeline (with exact retry) must succeed.
- Benchmarks: `cmake --preset gcc-bench && cmake --build --preset gcc-bench` (needs CGAL),
  `build/gcc-bench/pd_bench --csv bench.csv`, `scripts/summarize_bench.py bench.csv`.
- `scripts/repro_paper.sh` regenerates every result of the accompanying paper in
  `build/paper-results/`.

## Third-party code

- `tests/third_party/doctest.h`: doctest 2.4.9 (MIT license), used only by the test suite.
- `third_party/polypartition/`: PolyPartition by Ivan Fratrić, commit `b000a4a` (MIT license),
  compiled only into `pd_bench` as a baseline. Not part of the library and not installed.
