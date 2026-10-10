# polydecomp user guide

polydecomp decomposes a **polygon with holes** into **convex pieces** in two stages:

1. **Hole elimination** (`pd::eliminate_holes`) cuts the polygon along lines derived from the
   hole geometry until every piece is a simple polygon without holes.
2. **Convex decomposition** (`pd::convex_decompose`) splits each simple polygon at a reflex
   vertex along the extension of the vertex's incoming edge, until every piece is convex.

`pd::decompose` runs both stages. Every function works with an **exact** kernel (GMP rationals)
or a **double-precision** kernel, and every exact result can be checked with an independent
**exact certificate** (`pd::verify_partition`).

## Contents

1. Building and installing
2. Input
3. Quick start
4. Kernels and what they guarantee
5. The algorithms and their guarantees
6. Results
7. Errors
8. Input and output formats
9. Command-line tools
10. Performance
11. Limitations
12. Reproducing the paper's results

## 1. Building and installing

Requirements: a C++17 compiler (GCC or Clang), CMake 3.20 or newer, and GMP with its C++
bindings for the exact kernel (`libgmp-dev`; optional, see below). Ninja is used by the presets.

```sh
cmake --preset gcc-release
cmake --build --preset gcc-release
ctest --preset gcc-release
cmake --install build/gcc-release --prefix /usr/local
```

Other projects then use the library through CMake:

```cmake
find_package(polydecomp 1.0 REQUIRED)
target_link_libraries(my_target PRIVATE pd::polydecomp)
```

CMake options: `PD_WITH_GMP` (exact kernel; default ON), `PD_BUILD_TESTS`, `PD_BUILD_EXAMPLES`,
`PD_WITH_CGAL` (CGAL baselines in `pd_bench`), `PD_FUZZ` (libFuzzer targets, clang). Without GMP
the library provides the double-precision kernels only (preset `gcc-nogmp`); there the double
kernel's ordering of two chord crossings along a ray has no exact fallback and decides near ties
in floating point (orientation tests stay exact).

Presets: `gcc-debug`, `gcc-release`, `clang-debug`, `clang-release`, `gcc-asan`, `clang-asan`
(AddressSanitizer and UndefinedBehaviorSanitizer), `gcc-nogmp`, `gcc-bench` (with CGAL) and
`clang-fuzz`. `scripts/check.sh` builds and tests the first six. The API documentation is built
with `cmake --build --preset gcc-release --target docs` (needs Doxygen).

## 2. Input

A polygon with holes (`pd::PolygonWithHoles<FT>`) is an outer ring and a list of hole rings. A
ring (`pd::Ring<FT>`) lists its vertices once; the closing edge back to the first vertex is
implicit. The supported **input class** is:

- every ring has at least three vertices and no two equal consecutive vertices;
- every ring is simple (no two edges touch, except neighbouring edges at their shared vertex);
- every hole lies in the interior of the outer ring;
- the holes are pairwise disjoint as closed sets (no two hole rings touch) and not nested.

Ring orientations are arbitrary; the library normalises them. Holes may interleave freely
("entangled" holes, for example a hole inside the pocket of a C-shaped hole). `pd::validate`
checks the input class exactly and reports every violation (`pd::Issue`); the algorithms call
it first unless told not to (`validate_input = false`).

With the double-precision kernel every coordinate must also be 0 or have magnitude between
2^-400 and 2^400 (about 1e-120 to 1e120), so that the exact-sign predicates cannot overflow
or underflow. The exact kernel has no such limit.

## 3. Quick start

```cpp
#include "pd/pd.hpp"

pd::PolygonWithHoles<double> input{
    {{0, 0}, {10, 0}, {10, 10}, {0, 10}},   // outer ring
    {{{3, 3}, {6, 3}, {6, 6}, {3, 6}}}};    // one hole

// Exact: convert the doubles exactly to rationals, decompose, certify.
using K = pd::ExactKernel;
auto exact = pd::convert_polygon<K>(input);
auto result = pd::decompose<K>(exact);
if (result.ok()) {
  auto pieces = result.leaf_rings();                        // convex, counter-clockwise
  auto cert = pd::verify_partition<K>(exact, pieces, true);  // cert.ok == true
}

// Fast: double precision, with an exact rerun if a numerical problem is detected.
pd::PipelineOptions opt;
opt.on_failure = pd::OnFailure::RetryExact;
auto fast = pd::decompose<pd::DoubleKernel>(input, opt);
```

`examples/` holds three complete programs: `holes_to_convex.cpp` (exact pipeline with the
certificate), `double_with_retry.cpp` (double precision on a nearly degenerate input) and
`convex_and_tree.cpp` (convex decomposition and its tree).

## 4. Kernels and what they guarantee

A kernel supplies the number type `FT` and the orientation predicate. Every algorithm is a
template on the kernel.

| Kernel | Coordinates | Decisions | Constructions |
|---|---|---|---|
| `pd::ExactKernel` | GMP rationals (`mpq_class`) | exact | exact |
| `pd::SymbolicKernel` | rationals with an infinitesimal ε (`pd::EpsQ`) | exact | exact |
| `pd::DoubleKernel` | `double` | exact on stored points | rounded once |
| `pd::NaiveDoubleKernel` | `double` | plain floating point | rounded |

**ExactKernel.** Every predicate and every constructed point is exact, so the theorems below
hold without any general-position assumption. Input doubles are converted exactly
(`pd::convert_polygon`); a double such as 0.1 becomes its exact binary value, not 1/10. Every
constructed point is the intersection of two lines through input data, so the size of its
exact representation is bounded independently of the recursion depth.

**SymbolicKernel** is used automatically, never directly: when an exact run meets a dangling
contact (see Section 5), `eliminate_holes` and `decompose` redo the input with this kernel and
evaluate the result at ε = 0 (`stats.symbolic` is set). Its values are rational functions of ε,
compared as ε → 0+.

**DoubleKernel.** Coordinates are doubles, but every decision about stored points (the side of
a line a point lies on, and the orderings derived from it) is exact: a floating-point filter
with a proven error bound, falling back to an exact expansion computation. Each constructed
point is rounded once. To keep rounded constructions consistent the kernel additionally
- keeps every cut line at least 1e-9 times the coordinate scale away from every vertex that
  does not lie exactly on it, and replaces a line with a dangling contact by the fallback line
  (no symbolic split in floating point);
- does not split vertices whose turn is below 1e-12 radians (pieces are convex up to that
  angle), because the needle-shaped pieces this would create cannot be represented in double
  precision;
- checks each result before returning it: every piece simple and counter-clockwise, convex up
  to 1e-12 radians, and the piece areas summing to the input area within a relative 1e-9.
A detected problem is reported as `pd::Status::NumericalFailure` rather than returned as a
wrong result; with `PipelineOptions::on_failure = OnFailure::RetryExact` the input is re-run
with the exact kernel instead (requires GMP; `stats.retried_exact` is set).

**NaiveDoubleKernel** evaluates the predicate in plain floating point and has none of these
protections. It exists only to measure what they are worth.

## 5. The algorithms and their guarantees

### Hole elimination

A node is an outer ring with the input holes it contains. Each cut is a line:

- one hole: through the anchor `a = (v_{m-1} + v_m + v_{m+1}) / 3`, where `v_m` is the outer
  vertex nearest to the hole's vertex-average centroid `c`, and through `c`;
- two or more holes: through the vertex-average centroids of the two holes whose centroids are
  farthest apart.

The line may pass through vertices and run along edges (*contacts*). The pieces on each side
S are traced as if the line were shifted infinitesimally into S: an edge contributes a
crossing to S's list if the line crosses its interior, or if one endpoint lies on the line and
the other strictly on S, in which case the crossing point is that endpoint. Hence

- a vertex whose two edges lie on opposite sides is a crossing on both sides;
- a vertex whose two edges both lie on S is listed twice on S (it splits the material on S when
  the angle between the edges is not material) and not on the other side;
- a run of edges on the line contributes its end vertices to the sides of its outer edges.

Each side's list is walked along the line and the rings; pierced holes become notches of the
pieces, the other holes go to the piece that contains them. `stats.vertex_contacts` and
`stats.edge_contacts` count the contacts.

A hole that meets the line in a single vertex v and lies otherwise on one side (a *dangling
contact*) would give a piece that touches itself at v. The exact path then splits that hole
with the point c = v + ε d on the line (ε an infinitesimal), so that the contact becomes an edge
on the line, runs both stages on the perturbed polygon, and evaluates the result at ε = 0
(`stats.dangles`, `stats.symbolic`). The double kernel uses the fallback line instead.

The line is replaced by a fallback only if it has no direction (two centroids coincide, or
`a == c`): a line through the centroid of one of those holes (or, if that centroid is a vertex,
the midpoint of the hole's longest edge) in a direction, taken in the widest gap between the
directions to the vertices, that avoids every vertex. Every cut crosses at least one hole.

With the exact kernel: every cut removes at least one hole, so there are at most `h` cuts for
`h` holes; the pieces are hole-free and tile the input (their interiors are pairwise disjoint
and their closures cover it). They are simple, except after a dangling contact, where a piece
touches itself at that vertex (it is weakly simple).

### Convex decomposition

While a piece has a reflex vertex, the first one `v_i` is removed by a chord along the
extension of the incoming edge `(v_{i-1}, v_i)`, up to the first boundary contact. The contact
is inside an edge (a new vertex) or exactly at a vertex (a diagonal).

With the exact kernel: a polygon with `n_c` reflex vertices gives at most `n_c + 1` convex
pieces. Equality holds unless a chord ends at a reflex vertex that becomes convex on both
sides; `stats.double_resolutions` counts those. By Chazelle's lower bound `ceil(n_c / 2) + 1`
on the number of convex pieces of any partition, this is always within a factor 2 of the
optimum. The pieces tile the polygon. Chords end at points inside edges (Steiner points), so
the pieces need not meet edge to edge: a vertex of one piece can lie inside an edge of another
(a T-junction).

### The pipeline

Stage 1 produces simple counter-clockwise pieces, which is exactly Stage 2's input class, so
the pipeline's output is a set of convex pieces that tile the input polygon with holes.
Hole-free pieces that are already convex pass through unchanged. After a dangling contact both
stages run on the perturbed (simple) pieces; at ε = 0 no convex piece contains both v and c,
so the convex pieces are simple in every case.

### The certificate

`pd::verify_partition(input, pieces, require_convex, weakly_simple = false)` accepts the pieces
if and only if every piece is simple (weakly simple, if so requested: it may pass through a
point twice without crossing itself) and counter-clockwise (and convex, if required) and the oriented boundaries of
the pieces add up exactly, line by line, to the boundary of the input. By the additivity of
winding numbers this proves that the pieces tile the input. It shares no code with the
algorithms, allows T-junctions, and runs in O(E log E) for E edges in total. It is exact with
`ExactKernel`; for double results use `pd::check_pieces_approx` (necessary conditions only) or
the GEOS oracle in `scripts/oracle/`.

## 6. Results

All three entry points return a result with

- `pool`: every vertex, input vertices first, then constructed points;
- `leaves`: the output pieces as counter-clockwise lists of pool indices (`leaf_ring(i)` and
  `leaf_rings()` give the coordinates);
- `stats`: counts (cuts, fallbacks, splits, vertex hits, depths, ...);
- `diag`: the `pd::Status` and, on failure, a message and the validation issues.

`HoleElimination::tree` records every cut (line, pierced holes, children).
`ConvexDecomposition::tree` records every chord; `pd::reconstruct_ring` glues the pieces back
together into any node's ring, the root included. `PipelineResult` additionally keeps the
hole-free pieces (`hole_free`), which piece each convex piece came from (`leaf_source`) and the
Stage 1 tree.

## 7. Errors

The library does not throw from the algorithms; it reports through `diag.status`:

- `Status::InvalidInput`: the input is outside the input class (`diag.issues` lists why);
- `Status::NumericalFailure`: a floating-point kernel detected an inconsistency (retry with
  `ExactKernel`, or use `OnFailure::RetryExact`);
- `Status::InternalError`: an invariant failed in exact arithmetic. This indicates a bug;
  please report it with the input.

The readers throw `pd::ParseError` on malformed files; `pd::convert_polygon` throws
`std::domain_error` on a non-finite coordinate.

## 8. Input and output formats

- **WKT** (`pd::parse_wkt`, `pd::read_wkt_file`, `pd::write_wkt`,
  `pd::write_wkt_multipolygon`): POLYGON and MULTIPOLYGON. Numbers are read correctly rounded
  and written in the shortest form that reads back to the same double, so a round trip is
  bit-exact.
- **Exact rationals** (`pd::write_rational`, `pd::read_rational`,
  `pd::write_rational_rings`): one `p/q` pair per vertex; loses nothing.
- **Legacy OBJ** (`pd::read_obj_ring`, `pd::read_obj_rings`): rings exported from Blender, one
  ring per file or several as edge cycles; plane X-Z by default.

## 9. Command-line tools

- `pd_decompose`: decomposes every polygon in a WKT file, or an OBJ input. Options: `--stage
  both|holes|convex`, `--kernel exact|double`, `--verify` (exact certificate), `--retry-exact`,
  `--out` (WKT), `--rational`, `--stats` (JSON lines), `--compare double|naive` (classify a
  floating-point run against the exact one).
- `pd_validate`: checks the input class with both kernels.
- `pd_gen`: seeded test inputs (families `row`, `star`, `c`, `cluster`, `grid`, and the
  numerically adversarial `nearcollinear`, `neartouch`, `snapped`; combs; degenerate cases).
  The same family and seed give the same polygon on every platform.
- `pd_bench`: piece counts and timings, with CGAL's decompositions as baselines when built with
  `PD_WITH_CGAL`.

## 10. Performance

Convex decomposition performs `O(n (n_c + 1))` predicate evaluations for `n` vertices and
`n_c` reflex vertices, so it is quadratic in the worst case (for example on combs). Hole
elimination makes at most `h` cuts, each linear in the size of the current piece apart from
sorting the crossings and assigning the remaining holes. In the benchmarks of version 1.1.0 the
exact kernel is about 10 to 90 times slower than the double kernel (median about 40). Indicative
figures on one core of a Xeon E5-2667 v3 with frequency scaling active: a 2048-vertex
star-shaped polygon takes about 10 ms in double precision and 0.9 s exactly; a polygon with 256
holes and about 1000 vertices takes about 4 ms and 0.22 s.

## 11. Limitations

- The input class excludes self-intersecting rings, holes touching each other or the outer
  ring, and nested holes (islands inside holes).
- Pieces are not conforming: a vertex of one piece can lie inside an edge of a neighbour.
- The convex decomposition is not optimal. It is within a factor 2 of the optimum; on random
  star-shaped polygons Greene's optimal partition without Steiner points, as provided by CGAL,
  is about 10% smaller.
- Exact results rounded to doubles can contain pieces too small to survive rounding, when
  input vertices lie within a few units in the last place of each other.
- The double kernel's outputs are convex up to 1e-12 radians and can differ structurally from
  the exact outputs on degenerate data (they remain valid tilings).
- After a dangling contact, the hole-free pieces of `eliminate_holes` (exact) include a weakly
  simple piece that touches itself at the contact vertex; the pipeline's convex pieces are
  simple. The double kernel uses a fallback line for dangling contacts.

## 12. Reproducing the paper's results

`scripts/repro_paper.sh` rebuilds everything, runs the tests on all presets, validates the
legacy inputs, decomposes and certifies 5000 generated polygons, checks them with an
independent GEOS oracle, runs the floating-point census and the benchmarks against CGAL, and
writes all tables to `build/paper-results/`. It needs GMP, CGAL and the oracle's Python
environment (`scripts/oracle/setup_venv.sh`). Timings depend on the machine; the script
records the processor and its frequency policy. It also records the double kernel with exact
retry on the census inputs (`RETRY.txt`), the test count (`TESTS.txt`) and the degenerate set
(`DEGENERATE.txt`). The fuzzing campaign is run separately with `scripts/run_fuzz.sh [SECONDS]`
(clang, preset `clang-fuzz`; logs in `build/clang-fuzz/`).
