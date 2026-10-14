# PolyPartition (vendored, benchmark baseline only)

Unmodified copy of `src/polypartition.h`, `src/polypartition.cpp` and `LICENSE.md` from
https://github.com/ivanfratric/polypartition at commit
`b000a4a2a72b46e1305fb6e95b080448d7c12049` (2024-11-25), MIT License (see `LICENSE.md`).

It is compiled only into `pd_bench`, where `TPPLPartition::ConvexPartition_HM` (holes removed by
bridging, ear clipping, then Hertel and Mehlhorn's merge) serves as a baseline. It is not part of
the polydecomp library and is not installed.
