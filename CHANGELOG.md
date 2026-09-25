# Changelog

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
