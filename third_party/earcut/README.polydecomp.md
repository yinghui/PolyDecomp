# earcut.hpp (vendored, benchmark baseline only)

Unmodified copy of `include/mapbox/earcut.hpp` and `LICENSE` from
https://github.com/mapbox/earcut.hpp at tag `v2.2.4`, commit
`4811a2b69b91f6127a75e780de6e2113609ddabb` (2022-09-06), ISC License (see `LICENSE`).

It is used only by `pd_bench`, where Earcut's triangulation of a polygon with holes (holes bridged
into the outer ring, then ear clipping in double precision) serves as a baseline for removing holes
(`--suite hole-scale`) and its triangles are checked by the exact certificate (`--suite
triangulate`). It is not part of the polydecomp library and is not installed.
