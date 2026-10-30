# Changelog

## 1.4.3 (2026-10-30)

- Automatic cut lines no longer depend on the infinitesimal of a symbolic run. After a dangling
  contact, the anchor of the single hole rule (and of `CutRule::NearestToBoundary`) averages the
  values of its vertices at ε = 0. This fixes three cases found in proof reviews: a collapsed
  two vertex piece reported as success with the default rule, a piece pinched along an outer edge
  that the weakly simple certificate rejected, and an internal error of `NearestToBoundary` after
  a later dangling split.
- A symbolic result omits every trace whose exact area is zero at ε = 0 (`stats.collapsed`
  counts them, and the tree's leaf indices are renumbered). `pd_decompose` reports `collapsed`.
- `Slicer`: the floating point box filter is no longer used for a line whose coefficient rounds
  to zero or to a subnormal double. Such a line could skip a piece that it crosses.
- Five regression tests (`tests/test_stage1.cpp`, `tests/test_slice.cpp`).

## 1.4.2 (2026-10-25)

- Fix retained-hole matching in `verify_slices`: compare canonical cyclic boundaries,
  allowing rotation and reversal, rather than unordered vertex sets. Different simple
  holes and self-crossing permutations with the same vertices are now rejected.
  Regression tests also cover valid rotations and reversals.

## 1.4.1 (2026-10-17)

- `pd_bench --suite slice-size [--cuts N]`: the sizes of the states after the last cut of every
  sequence of the slice suite, deterministic and untimed: the pieces of `pd::Slicer`, their distinct
  points and edges, against the CGAL arrangement's vertices, edges and faces, its faces split into
  material faces, faces inside holes and faces outside the shape. `scripts/repro_paper.sh` runs it
  as step 13 (`SLICE_SIZE.csv`). The arrangement baseline marks faces inside holes with 2 (outside
  the shape stays 0); timings and pieces are unchanged.
- Stronger baselines for removing holes. `third_party/earcut/`: earcut.hpp 2.2.4 (Mapbox, ISC
  license), used only by `pd_bench`. `pd_bench --suite hole-scale` now also runs Earcut's
  triangulation and, with CGAL, CGAL's constrained triangulation (`cgal-cdt`) and `connect_holes`
  (bridging by vertical segments, `cgal-connect-holes`); the triangles of both triangulations are
  checked by the exact certificate (zero-area triangles dropped).
- `pd_bench --suite triangulate --wkt FILE [--wkt FILE ...]`: Earcut and CGAL's constrained
  triangulation on every polygon of the files, each result under the exact certificate, with the
  exact relative area error (`TRIANGULATE.csv`, repro step 14, on the census inputs). On the 5900
  census inputs CGAL's triangles always pass; Earcut's fail on 2 near-collinear inputs, where its
  floating-point collinearity test drops vertices a few units in the last place off an edge and
  its triangles miss slivers (relative area below 1e-17).

## 1.4.0 (2026-10-17)

The paper now covers exact cutting and hole elimination; the convex decomposition stays in the
library, unchanged.


- `eliminate_holes` with a floating-point kernel now checks its result (`HoleElimOptions::self_check`,
  on by default: simple counter-clockwise pieces whose areas sum to the input's within a relative
  1e-9) and can answer a numerical failure with an exact rerun
  (`HoleElimOptions::on_failure = OnFailure::RetryExact`, `stats.retried_exact`), as `decompose`
  does. Before, a rounded construction could make the double kernel return hole-free pieces whose
  areas miss the input area without a report (2 of 5900 census inputs, both near-touching holes).
  `OnFailure` moved to `pd/hole_elim.hpp` (still available through `pd/pipeline.hpp`).
- `pd_decompose --stage holes` accepts `--retry-exact` and `--compare double|naive` (the census of
  hole elimination alone: agree, valid deviation, wrong, reported failure).
- `pd_bench --suite hole-scale`: hole elimination with the exact and the double kernel and
  PolyPartition's `RemoveHoles` (bridging), on random holes with 2^4 to 2^14 holes and cheese with
  2^4 to 2^12 holes, with the time budget of the scale suite, and input validation.
- `pd_bench --suite slice`: a fourth method, `cgal-arrangement`, cuts with a CGAL arrangement of
  segments on Epeck (each line inserted once, face flags kept by an observer); its material faces
  are certified after the last cut.
- `pd_bench --suite strokes`: a game workload. Pieces carry rigid motions and drift apart after
  every cut; a finite stroke drawn in the world cuts the pieces it crosses, along its line mapped
  into each piece's frame (`Slicer::cut(i, a, b)`); the pieces pass the exact certificate in their
  own frames after every stroke.
- `pd_bench --suite strokes --world FILE` writes the pieces of the first shape after the last
  stroke, placed in the world by their motions (WKT, one polygon per line), for figures.
- `pd_decompose --compare` records `other_fallbacks`, the cuts of the floating-point run that used a
  fallback line, so that deviations of the double kernel can be traced to its clearance.
- The GEOS oracle (`scripts/oracle/oracle.py`) takes `--rational FILE` (the exact pieces of
  `pd_decompose --rational`): a piece that GEOS finds invalid is checked again in exact rational
  arithmetic and counted as "invalid only after rounding" if its exact ring is simple. Rounding an
  exact result to doubles can merge two vertices closer than a unit in the last place (for example
  the two crossings of a cut that passes next to a hole corner), so that the rounded ring touches
  itself; 6 pieces of 5 generated grid inputs do this, while their exact rings are simple.
  `run_generated.sh` writes the exact pieces (`*.pieces.rat`) for this check.
- The fuzz targets also run hole elimination with the double kernel and `RetryExact` on every
  valid input and require it to succeed.
- `scripts/run_generated.sh`, `run_census.sh` and `run_ablation.sh` take an optional stage
  (`holes` runs hole elimination only); `scripts/repro_paper.sh` regenerates the results of the
  paper on exact cutting and hole elimination (steps 1 to 12; the convex benchmarks, `BENCH.csv`
  and `SCALE.csv`, are no longer part of it and remain available through `pd_bench`). It now clears
  the results of an earlier run first.

## 1.3.0 (2026-10-16)

- New: `pd::Slicer` (`pd/slice.hpp`, included by `pd/pd.hpp`; needs GMP) cuts a polygon with holes
  along arbitrary lines, as often as wanted, in exact arithmetic, for interactive cutting (for
  example a player slicing a cheese with many holes in a game). `cut(a, b)` splits every piece that
  the line through a and b crosses, `cut(i, a, b)` only piece i, and `cut_exact` takes rational
  points. Cuts follow the contact rules of hole elimination: lines through vertices, along edges,
  repeating an earlier cut or passing through a point made by one are all handled, holes that the
  line misses stay holes of their piece, and a hole that the line touches in a single vertex is
  split symbolically (the state switches to `SymbolicKernel` at the first such contact). Every new
  point is the intersection of two supporting lines (input edge lines and cut lines), so
  coordinates stay bounded however many cuts follow (`max_bits()`). Results: `piece(i)` (doubles,
  for display), `piece_exact(i)` and `convex(i)` (the pipeline on one piece). The double kernel is
  not offered for slicing, since its fallback would move the player's line.
- New tool `pd_slice`: cuts a WKT polygon along given lines (`--cut`, `--graze H J` for the exact
  tangent at a hole vertex, `--through` for the line through two input vertices), writes the
  pieces and optionally their convex pieces, and
  certifies them (`--verify`).
- New: `verify_slices(input, pieces, weakly_simple)` certifies pieces that still contain holes:
  every hole of a piece is an input hole, used by one piece only, that lies inside that piece (a
  vertex strictly inside its outer ring and no edge meeting that ring), and the outer rings tile the
  input with those holes filled (`verify_partition`).
- Hole elimination: `HoleEliminator::run` is split into `init` and the cut loop, and `cut_node`
  into `choose_cut` and `cut_along(work, cut, children, require_pierce)`, which the slicer calls
  with an arbitrary line. Automatic elimination is unchanged: the generated inputs, the degenerate
  set, the benchmarks, the scaling suite and the ablation give the same pieces as 1.2.0.
- `SymbolicKernel` is faster. Every value is now kept in lowest terms (the polynomial gcd of its
  numerator and denominator cancelled), so that, as with `mpq_class`, its representation does not
  depend on how it was computed. Before, values computed through several constructions carried
  ever larger polynomials: after a dangling contact, the convex decomposition of slices built
  points of degree 10 with coefficients of 4500 bits whose true value has degree 1 or none, and
  decomposing the new pieces of one cut of a cheese with 256 holes took up to 23 s (now at most
  0.43 s, as with the exact kernel). A value that does not depend on ε is stored as a plain rational,
  and the orientation test of three such points uses the exact kernel's code, which reuses its
  temporaries; before, a symbolic run took about 1.7 times as long as an exact one even where no ε
  occurred, now about as long. The evaluation of a symbolic pipeline result at ε = 0 is factored
  out (`pipeline_limit`).
- Generator `pd::gen::cheese(h, seed)` (`pd_gen --cheese H`): a roughly round outer ring with h
  round holes (8 to 16 vertices) of random radii at random positions, at least 0.15 apart.
- `pd_bench --suite slice [--cuts N]`: cheese with 16, 64 and 256 holes, 64 random holes and the
  bundled sample polygon, each cut N times (default 200) by random lines, axis-parallel lines
  through vertices, lines through two vertices, lines that repeat earlier cuts and (on cheese)
  lines tangent to a hole at a vertex. Per cut: time, pieces, contacts, dangling contacts, the
  time to decompose the new pieces into convex pieces, the largest coordinate size, and the exact
  certificate of all pieces. The same lines cut the shape with two baselines: Clipper2 (vendored
  unmodified in `third_party/clipper2/`, Boost Software License, 64-bit integer coordinates) and,
  with CGAL, CGAL's Boolean set operations on Epeck with segment traits (with the default polyline
  traits, CGAL's own output after a cut through a hole vertex violates the input precondition of
  the next operation). `--suite clipper-scale` runs Clipper2 at five integer scales. Step 12 of
  `scripts/repro_paper.sh` (`SLICE.csv`, `CLIPPER_SCALE.csv`, summarized by
  `scripts/summarize_slice.py` in `SLICE.md`).
- `scripts/repro_paper.sh` samples the average frequency of the benchmark CPU every half second while
  each benchmark runs and records the mean and the minimum in `MACHINE.txt`. Before, it read the
  frequency once after each benchmark, which could hit a momentary drop below the nominal frequency
  (seen once, with the `performance` governor and turbo off).
- Fuzzing: a third libFuzzer target, `fuzz_slice`, cuts a grid polygon (the decoder of `fuzz_grid`,
  now in `fuzz/fuzz_common.hpp`) up to eight times, mostly through vertices, along edges and
  tangent to holes, and certifies the pieces after every cut and the convex pieces of every piece
  at the end. `scripts/run_fuzz.sh` runs all three targets.

## 1.2.0 (2026-10-15)

- Fix (found by `fuzz_grid`, present since 1.1.0): after a dangling contact, a hole-free piece that
  runs along the edge from the dangling vertex to its split vertex repeated that vertex at ε = 0, so
  `eliminate_holes` and `PipelineResult::hole_free` could return a ring that the certificate rejects.
  Coincident consecutive vertices are now removed when the symbolic result is evaluated. The convex
  pieces were not affected. Regression tests at both levels.
- `scripts/run_fuzz.sh` reports crash reports from the logs (in fork mode a crash stops the campaign
  while the status counters still read 0).
- Hole elimination: option `cut_rule` (`HoleElimOptions`, `PipelineOptions`; `pd_decompose
  --cut-rule farthest|closest|boundary|horizontal`) to compare the choice of the cut line. The
  default `FarthestPair` is the rule of earlier versions, and results are unchanged. Every rule
  draws the line through a hole's vertex average, so every cut crosses a hole; only the default is
  covered by the bounded-coordinates guarantee. `scripts/run_ablation.sh` compares the rules
  (step 11 of `repro_paper.sh`, `ABLATION.txt`).
- Generators: `random_holes(h, seed)`, h random star-shaped holes (3 to 10 vertices) in a grid of
  unit cells, valid for every h; `pd_gen --random-holes H` and `pd_gen --star N` (with `--seed`).
- `validate`: the hole-in-hole test first compares the point with the hole's bounding box (exact),
  so the quadratic number of point-in-ring tests is avoided; 16384 random holes now validate in
  0.7 s instead of 21 s. The edge sweep is unchanged: it takes O(E log E + K) time, K being the
  number of edge pairs whose x ranges overlap (O(E) for convex rings, up to quadratic otherwise).
- `pd_bench --suite scale [--budget-ms T]`: random star polygons and combs with 2^10 to 2^16
  vertices and `random_holes` with 2^4 to 2^14 holes; within each series a method is run on the
  next size only while its median time stays at most T ms (default 5000). The suite also times
  input validation (`validate-double`, `validate-exact`).
- `pd_bench`: the untimed first run is now the one checked by the exact certificate (before, an
  extra run after timing), so every method runs once less.
- `pd_bench`: `cgal-cdt-hm` keeps its pieces as linked half-edges, so each merge takes constant
  time. Before, every merge copied a piece, which was quadratic on combs, where one piece collects
  the whole base. Same merge order and tests, so the pieces are unchanged.
- `scripts/repro_paper.sh`: records the frequency policy of the benchmark CPU (governor, min/max,
  turbo, average frequency after each benchmark), warns unless the governor is `performance`, and
  adds step 10, the scaling suite (`SCALE.csv`, `SCALE.md`).

## 1.1.2 (2026-10-14)

Benchmark harness only; the library is unchanged.

- `pd_bench`: two Hertel and Mehlhorn baselines that accept holes. `polypartition-hm` is
  PolyPartition's `ConvexPartition_HM` (holes bridged, ear clipping, merge; double precision),
  vendored unmodified in `third_party/polypartition/` (MIT). `cgal-cdt-hm` merges the triangles of
  CGAL's constrained triangulation, removing every diagonal whose endpoints stay convex (CGAL's
  exact orientation test).
- `pd_bench`: new CSV column `certified`. After timing, one untimed run of every method whose
  pieces use only input vertices (`ours-exact`, `polypartition-hm`, `cgal-hertel-mehlhorn`,
  `cgal-greene-approx`, `cgal-optimal`, `cgal-cdt`, `cgal-cdt-hm`) is checked by the exact
  certificate with convexity required.

## 1.1.1 (2026-10-11)

- Degenerate set: new case `ring-edge-contacts` (edges of the outer ring on the cut line, a wing
  and a wedge), with a test.

## 1.1.0 (2026-10-10)

- Hole elimination: cut lines may now pass through vertices and run along edges. Each side's
  crossings are those of the line shifted infinitesimally into that side, so every contact is
  classified exactly and the line through the hole centroids is kept; the fallback line is used
  only when that line has no direction (and, in floating point, for near-contacts within the
  clearance and for dangling contacts). New statistics: `vertex_contacts`, `edge_contacts`.
- Dangling contacts (a hole meeting the line in a single vertex): the exact path splits the hole
  symbolically with c = v + ε d, reruns with the new `SymbolicKernel` (`EpsQ`: rationals with an
  infinitesimal ε) and evaluates the result at ε = 0. Statistics `dangles`, `symbolic`. The
  hole-free pieces are then weakly simple; the convex pieces are simple.
  A convex piece whose area vanishes at ε = 0 is dropped (`stats.vanished`; none occurred in
  the tests or the paper's runs).
- `verify_partition` gained a `weakly_simple` mode (for hole-free pieces after a dangling
  contact); `pd_decompose --stage holes --verify` uses it when needed.
- Degenerate set: new cases `vertex-contacts`, `edge-contacts`, `dangling-contacts`.
- Tests for every contact type, the symbolic kernel, and symbolic-vs-exact agreement.

## 1.0.1 (2026-09-25)

- Fix: `pd_gen --family nearcollinear` failed for some seeds. Nudging a coordinate 0 by one ulp
  gave a subnormal, which the double kernel's coordinate-range check (added in 1.0.0) rejects.
  The family now lies in [1, 11]^2. A test generates every census seed of all adversarial
  families.
- `scripts/repro_paper.sh` also records the double kernel with exact retry on the census inputs
  (`RETRY.txt`, via the new `scripts/summarize_retry.py`), the unit-test count (`TESTS.txt`) and
  the degenerate set with certificate (`DEGENERATE.txt`, `DEGENERATE.jsonl`).
- New `scripts/run_fuzz.sh`: the fuzzing campaign reported in the paper.
- User guide: reproduction steps updated; the limitation of the build without GMP (the ordering of
  two chord crossings has no exact fallback near ties) is documented.

## 1.0.0 (2026-09-25)

First public version.

- Hole elimination (`eliminate_holes`): line cuts that are generic by construction; at most h
  cuts for h holes, with no general-position assumption.
- Convex decomposition (`convex_decompose`): chord splitting with exact vertex-hit handling; at
  most n_c + 1 convex pieces for n_c reflex vertices (a factor two of the optimum).
- The pipeline (`decompose`) from a polygon with holes to convex pieces.
- Exact kernel (GMP rationals) and a robust double kernel (exact predicates on stored points,
  self-check, optional exact retry); a naive kernel for measurements.
- Exact partition certificate (`verify_partition`), input validation (`validate`).
- I/O: WKT, legacy OBJ, exact rational text.
- Tools: `pd_decompose`, `pd_validate`, `pd_gen`, `pd_bench`; fuzz targets; GEOS oracle and
  scripts that regenerate the paper's results.
