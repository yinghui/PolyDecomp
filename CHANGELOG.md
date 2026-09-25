# Changelog

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
