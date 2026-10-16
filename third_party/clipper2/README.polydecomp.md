# Clipper2 (vendored, benchmark baseline only)

Unmodified copy of `CPP/Clipper2Lib/include/clipper2/*.h`, `CPP/Clipper2Lib/src/*.cpp` and `LICENSE`
from https://github.com/AngusJohnson/Clipper2 at tag `Clipper2_2.0.1`, commit
`21ebba05db8894f0c7217ad35ea518080f324946` (2025-12-20), Boost Software License 1.0 (see `LICENSE`).

It is compiled only into `pd_bench`, where `Clipper64` (polygon clipping on 64-bit integer
coordinates) serves as a baseline for cutting along lines. It is not part of the polydecomp library
and is not installed.
