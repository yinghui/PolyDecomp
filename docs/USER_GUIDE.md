# polydecomp user guide

polydecomp decomposes a **polygon with holes** into **convex pieces** in two stages:

1. **Hole elimination** (`pd::eliminate_holes`) cuts the polygon along lines derived from the
   hole geometry until every piece is a simple polygon without holes.
2. **Convex decomposition** (`pd::convex_decompose`) splits each simple polygon at a reflex
   vertex along the extension of the vertex's incoming edge, until every piece is convex.

`pd::decompose` runs both stages. Every function works with an **exact** kernel (GMP rationals)
or a **double-precision** kernel, and every exact result can be checked with an independent
**exact certificate** (`pd::verify_partition`).

`pd::Slicer` cuts a polygon with holes along **arbitrary lines**, as often as wanted, in exact
arithmetic, for interactive cutting (for example a player slicing a cheese in a game); each piece
can then be decomposed into convex pieces.

## Contents

1. Building and installing
2. Input
3. Quick start
4. Kernels and what they guarantee
5. The algorithms and their guarantees
6. Cutting along any line
7. Results
8. Errors
9. Input and output formats
10. Command-line tools
11. Performance
12. Limitations
13. Reproducing the paper's results

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

`examples/` holds four complete programs: `holes_to_convex.cpp` (exact pipeline with the
certificate), `double_with_retry.cpp` (double precision on a nearly degenerate input),
`convex_and_tree.cpp` (convex decomposition and its tree) and `slice_cheese.cpp` (cutting along
arbitrary lines, Section 6).

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
evaluate the result at ε = 0 (`stats.symbolic` is set); `pd::Slicer` switches to it at its first
dangling contact (Section 6). Its values are rational functions of ε in lowest terms, compared as
ε → 0+. A value that does not depend on ε is stored as a plain rational, so the kernel runs at
about the speed of `ExactKernel` except near the split vertices.

**DoubleKernel.** Coordinates are doubles, but every decision about stored points (the side of
a line a point lies on, and the orderings derived from it) is exact: a floating-point filter
with a proven error bound, falling back to an exact expansion computation. Each constructed
point is rounded once. To keep rounded constructions consistent the kernel additionally
- keeps every cut line at least 1e-9 times the coordinate scale away from every vertex that
  does not lie exactly on it, and replaces a line with a dangling contact by the fallback line
  (no symbolic split in floating point);
- in the convex decomposition, does not split vertices whose turn is below 1e-12 radians
  (pieces are convex up to that angle), because the needle-shaped pieces this would create
  cannot be represented in double precision;
- checks each result before returning it: every piece simple and counter-clockwise, the piece
  areas summing to the input area within a relative 1e-9, and, for `decompose`, every piece
  convex up to 1e-12 radians. `eliminate_holes` does this check too (since 1.4.0,
  `HoleElimOptions::self_check`, on by default).
A detected problem is reported as `pd::Status::NumericalFailure` rather than returned as a
wrong result; with `on_failure = OnFailure::RetryExact` (in `PipelineOptions` or
`HoleElimOptions`) the input is re-run with the exact kernel instead and the result rounded to
doubles (requires GMP; `stats.retried_exact` is set). In the census of the paper (5900
inputs), the double kernel returned no wrong hole-free result, where plain floating point
returned 557, and its check caught 2 results whose areas were off (holes 2^-30 apart).

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
(`stats.dangles`, `stats.symbolic`). After such a split, a single hole anchor that would use a
split vertex averages the vertex values at ε = 0, so that every cut line is a fixed rational line.
A trace whose exact area is zero at ε = 0 encloses no material and is not reported
(`stats.collapsed`). The double kernel uses the fallback line instead.

The line is replaced by a fallback only if it has no direction (two centroids coincide, or
`a == c`): a line through the centroid of one of those holes (or, if that centroid is a vertex,
the midpoint of the hole's longest edge) in a direction, taken in the widest gap between the
directions to the vertices, that avoids every vertex. Every cut crosses at least one hole.

With the exact kernel: every cut removes at least one hole, so there are at most `h` cuts for
`h` holes; the pieces are hole-free and tile the input (their interiors are pairwise disjoint
and their closures cover it). They are simple, except after a dangling contact, where a piece
may touch itself at that vertex (it is weakly simple).

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

## 6. Cutting along any line

```cpp
#include "pd/pd.hpp"

using P = pd::Point<double>;
pd::Slicer s(pd::gen::cheese(64, 1));             // validates the shape (s.ok(), s.diag())
pd::SliceResult r = s.cut(P{-9, 0.5}, P{9, -1});  // every piece the line crosses
r = s.cut(0, P{0, -9}, P{0.2, 9});                // only piece 0
for (std::size_t i = 0; i < s.size(); ++i) {
  pd::PolygonWithHoles<double> shown = s.piece(i);    // rounded once, for display
  pd::PipelineResult<mpq_class> parts = s.convex(i);  // convex pieces, e.g. for physics
}
```

`examples/slice_cheese.cpp` is a complete program (cuts, a cut through a vertex made by an earlier
cut, the certificate and the convex pieces).

A `Slicer` keeps one exact state for all its pieces: the points, the supporting line of every
edge (input edges and cut lines) and the holes. `cut(a, b)` splits every piece that the line
through `a` and `b` crosses, `cut(i, a, b)` only piece `i`, and `cut_exact(a, b[, i])` takes exact
rational points (for example a vertex made by an earlier cut). A cut uses the contact rules of hole elimination
(Section 5) with the player's line instead of a chosen one, so

- lines through vertices, along edges, through points made by earlier cuts, or repeating an
  earlier cut are all handled exactly; a piece that the line misses or only touches is kept;
- the pieces always tile the input, holes that the line misses stay holes of the piece that
  contains them, and crossed holes become notches;
- a hole that the line touches in a single vertex is split symbolically, and from then on the
  state is held in `SymbolicKernel` (`symbolic()`); pieces are reported at ε = 0, where a piece
  that runs round that hole touches itself at one point (weakly simple);
- every new point is the intersection of two supporting lines, so the size of every coordinate
  is bounded by the sizes of the input and of the cut lines, however many cuts follow
  (`max_bits()`).

`SliceResult` lists the new pieces (`pieces`, indices into the new list), how many pieces were cut,
the contacts and the dangling contacts; a failure is reported in `diag` and leaves the state
unchanged. `piece_exact(i)` gives exact coordinates, and `pd::verify_slices(input, pieces,
weakly_simple)` certifies pieces that still contain holes: every hole of a piece is an input hole,
used by one piece only, that lies inside that piece (a vertex strictly inside its outer ring and no
edge meeting that ring), and the outer rings tile the input with those holes filled. `convex(i)` runs the pipeline on piece `i`.

**Moving pieces.** In a game the pieces move after a cut. Keep a rigid motion (rotation and
translation) per piece in the application, and keep the slicer's pieces in the frame of the
input. To cut with a line drawn in the world, map its two points into the frame of each piece
that it crosses with the inverse motion and call `cut(i, a, b)` for that piece, from the highest
index down so that the indices of the other pieces stay valid; the new pieces inherit the motion
of their parent. Each cut is exact for the mapped line; rounding the mapped points only moves the
line by a few units in the last place relative to the piece, and the pieces still tile the input
in its frame. `pd_bench --suite strokes` implements this workload (finite strokes, pieces that
drift apart, the certificate after every stroke).

## 7. Results

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

## 8. Errors

The library does not throw from the algorithms; it reports through `diag.status`:

- `Status::InvalidInput`: the input is outside the input class (`diag.issues` lists why);
- `Status::NumericalFailure`: a floating-point kernel detected an inconsistency (retry with
  `ExactKernel`, or use `OnFailure::RetryExact`);
- `Status::InternalError`: an invariant failed in exact arithmetic. This indicates a bug;
  please report it with the input.

The readers throw `pd::ParseError` on malformed files; `pd::convert_polygon` throws
`std::domain_error` on a non-finite coordinate.

## 9. Input and output formats

- **WKT** (`pd::parse_wkt`, `pd::read_wkt_file`, `pd::write_wkt`,
  `pd::write_wkt_multipolygon`): POLYGON and MULTIPOLYGON. Numbers are read correctly rounded
  and written in the shortest form that reads back to the same double, so a round trip is
  bit-exact.
- **Exact rationals** (`pd::write_rational`, `pd::read_rational`,
  `pd::write_rational_rings`): one `p/q` pair per vertex; loses nothing.
- **Legacy OBJ** (`pd::read_obj_ring`, `pd::read_obj_rings`): rings exported from Blender, one
  ring per file or several as edge cycles; plane X-Z by default.

## 10. Command-line tools

- `pd_decompose`: decomposes every polygon in a WKT file, or an OBJ input. Options: `--stage
  both|holes|convex`, `--kernel exact|double`, `--verify` (exact certificate), `--retry-exact`,
  `--out` (WKT), `--rational`, `--stats` (JSON lines), `--compare double|naive` (classify a
  floating-point run against the exact one, for the pipeline or, with `--stage holes`, for hole
  elimination alone; the records include `other_fallbacks`, the fallback lines of the
  floating-point run).
- `pd_validate`: checks the input class with both kernels.
- `pd_slice`: cuts the first polygon of a WKT file along lines (`--cut X0 Y0 X1 Y1`, `--graze H J`
  for the exact tangent at vertex J of hole H, or `--through H1 J1 H2 J2` for the line through two
  input vertices, ring -1 being the outer ring), writes the pieces as WKT with their remaining holes
  (`--out`) and, with `--convex`, the convex pieces of every piece; `--verify` runs the
  certificates.
- `pd_gen`: seeded test inputs (families `row`, `star`, `c`, `cluster`, `grid`, and the
  numerically adversarial `nearcollinear`, `neartouch`, `snapped`; combs; random star-shaped
  polygons `--star N`; many random holes `--random-holes H`; a cheese with H round holes
  `--cheese H`; degenerate cases).
  The same family and seed give the same polygon on every platform.
- `pd_decompose --cut-rule farthest|closest|boundary|horizontal` chooses the hole-elimination
  cut (`pd::CutRule`): through the vertex averages of the two holes farthest apart (default) or
  closest together, the one-hole rule for the hole nearest the outer ring, or a horizontal line
  through a hole's average. The alternatives exist for comparison; only the default carries the
  bounded-coordinates guarantee.
- `pd_bench`: piece counts and timings, with PolyPartition's Hertel and Mehlhorn partition
  (vendored) and, when built with `PD_WITH_CGAL`, CGAL's decompositions and a Hertel and Mehlhorn
  merge of CGAL's constrained triangulation as baselines. Results whose pieces use only input
  vertices are checked by the exact certificate (column `certified`). `--suite scale` measures
  how the running times grow, up to 2^16 vertices and 2^14 holes, with a time budget per method
  (`--budget-ms`), and also times input validation. `--suite slice [--cuts N]` cuts cheese, random
  holes and the sample polygon N times (default 200) by random and by degenerate lines and
  records, per cut, the time, the pieces, the contacts, the time of the convex decomposition of
  the new pieces, the largest coordinate size and the exact certificate, for our slicer and for
  three baselines on the same lines (CGAL's arrangement of segments and Boolean set operations,
  both exact, and Clipper2 on integers). `--suite hole-scale` times hole elimination with both
  kernels, PolyPartition's bridging, Earcut's triangulation and, with CGAL, CGAL's `connect_holes` and
  constrained triangulation on random holes and cheeses with up to 2^14 holes. `--suite triangulate
  --wkt FILE ...` checks Earcut's and CGAL's triangulations of every polygon in the files with the
  exact certificate.
  `--suite strokes [--cuts N] [--world FILE]` draws N finite strokes on cheeses whose pieces move
  apart (Section 6, *Moving pieces*), with the certificate after every stroke. `--suite slice-size`
  compares the size of the slicer's state after the last cut of every slice sequence with that of the
  CGAL arrangement (vertices, edges, and faces in the material, in holes and outside the shape).

## 11. Performance

Figures from the paper's run of version 1.4.0, on one core of a Xeon E5-2667 v3 at a fixed 3.2 GHz.

**Cutting.** A slicer cut compares the line with each piece's bounding box and walks only the pieces
it may cross, in time linear in their size apart from sorting the crossings. Over 200 full lines
across cheeses with 16 to 256 holes and other shapes, a cut took 1.3 ms in the median, including
cuts on more than ten thousand pieces at the end of a sequence (200 random lines across a cheese
give about ten thousand pieces), and coordinates never needed more than 350 bits. A CGAL arrangement
of the same lines took 0.6 to 1.8 times as long, CGAL's Boolean set operations 8 to 17 times as
long. A finite stroke on pieces that move apart (`pd_bench --suite strokes`) took 0.10 ms in the
median and at most 5.4 ms, since it cuts only the pieces it crosses.

**Hole elimination** makes at most `h` cuts, each linear in the size of the current piece apart
from sorting the crossings and assigning the remaining holes. With the double kernel it removes
16384 random holes (about 106000 vertices) in about 2.9 s and the 4096 holes of a cheese (about
49000 vertices) in about 0.13 s; its time grows with slopes 1.5 to 2.1 on logarithmic axes. The
exact kernel is 22 to 233 times slower from 64 holes on (2048 random holes in 11.5 s). Input
validation takes about 0.5 s on the largest random-holes input.

**Convex decomposition** performs `O(n (n_c + 1))` predicate evaluations for `n` vertices and
`n_c` reflex vertices, so it is quadratic in the worst case (for example on combs). In the
benchmarks of version 1.2.0 (`pd_bench --suite all|scale`), a 2048-vertex star-shaped polygon took
about 10 ms in double precision and 0.9 s exactly, and the double kernel decomposed a
65536-vertex star-shaped polygon in about 7.5 s.

## 12. Limitations

- The input class excludes self-intersecting rings, holes touching each other or the outer
  ring, and nested holes (islands inside holes).
- Pieces are not conforming: a vertex of one piece can lie inside an edge of a neighbour.
- Input validation (on by default in `decompose` and `eliminate_holes`) and the simplicity test of
  `verify_partition` sweep edge bounding boxes in O(E log E + K) time, K being the number of edge
  pairs whose x ranges overlap. K is O(E) for convex rings but can grow quadratically for large
  non-convex rings (for example random star-shaped polygons); `validate_input = false` skips the
  check for inputs known to be valid.
- The convex decomposition is not optimal. It is within a factor 2 of the optimum; on random
  star-shaped polygons Greene's optimal partition without Steiner points, as provided by CGAL,
  is about 10% smaller.
- Exact results rounded to doubles can contain pieces too small to survive rounding, when
  input vertices lie within a few units in the last place of each other, and a piece can touch
  itself after rounding when a cut passes a vertex closer than a unit in the last place (its exact
  ring is simple; `oracle.py --rational` recognises this).
- The double kernel's outputs can differ structurally from the exact outputs on degenerate data
  (they remain valid tilings, and its convex pieces are convex up to 1e-12 radians).
- After a dangling contact, the hole-free pieces of `eliminate_holes` (exact) include a weakly
  simple piece that touches itself at the contact vertex; the pipeline's convex pieces are
  simple. The double kernel uses a fallback line for dangling contacts.
- `pd::Slicer` works in exact arithmetic only (the double kernel would move a line that comes too
  close to a vertex). A cut is a full line across the pieces it is applied to; to cut along a
  finite stroke, apply the line to the pieces the stroke crosses. Convex pieces can have more
  vertices than some physics engines accept per polygon (Box2D's default is 8).

## 13. Reproducing the paper's results

`scripts/repro_paper.sh` rebuilds everything, runs the tests on all presets, validates the
legacy inputs, runs exact hole elimination on 5000 generated polygons with the certificate,
checks the hole-free pieces with the independent GEOS oracle (with the exact recheck of pieces
that are invalid only after rounding), runs the floating-point census of hole elimination (double
kernel and plain floating point against the exact result, the double results through the
oracle, and the exact retry), the degenerate set, hole elimination at scale (`HOLE_SCALE.csv`),
the cut-rule ablation (`ABLATION.txt`), the cutting experiment with its three baselines
(`SLICE.csv`, `SLICE.md`, and Clipper2 at several integer scales, `CLIPPER_SCALE.csv`) and the
strokes on moving pieces (`STROKES.csv`) and the sizes of the slicer's and the arrangement's states
(`SLICE_SIZE.csv`) and the exact check of two triangulations of the census inputs
(`TRIANGULATE.csv`), and writes everything to `build/paper-results/`. It
needs GMP, CGAL and the oracle's Python environment (`scripts/oracle/setup_venv.sh`). Timings
depend on the machine; the script records the processor and the frequency policy of the
benchmark CPU (`PD_BENCH_CPU`, default 2), samples its frequency during each benchmark, and warns
unless its governor is `performance`. For stable timings fix the frequency first (root), for
example `echo performance | sudo tee /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor` with
turbo disabled. The convex benchmarks of earlier versions (`pd_bench --suite all|scale`) are no
longer part of the script. The fuzzing campaign (three targets, the last one for the slicer) is
run separately with `scripts/run_fuzz.sh [SECONDS]` (clang, preset `clang-fuzz`; logs in
`build/clang-fuzz/`).
