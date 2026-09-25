// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#pragma once
/// \file
/// Basic geometric value types, parameterised by the number type FT.

#include <vector>

namespace pd {

template <class FT>
struct Point {
  FT x, y;
};

template <class FT>
bool operator==(const Point<FT>& a, const Point<FT>& b)
{
  return a.x == b.x && a.y == b.y;
}
template <class FT>
bool operator!=(const Point<FT>& a, const Point<FT>& b)
{
  return !(a == b);
}

/// A closed polygonal ring given by its vertices; the closing edge from the
/// last vertex back to the first is implicit (the first vertex is not repeated).
template <class FT>
using Ring = std::vector<Point<FT>>;

/// An outer ring with zero or more hole rings (a "complex polygon").
template <class FT>
struct PolygonWithHoles {
  Ring<FT> outer;
  std::vector<Ring<FT>> holes;
};

} // namespace pd
