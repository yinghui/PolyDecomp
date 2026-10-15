# Changelog

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
