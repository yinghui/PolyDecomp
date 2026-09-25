// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#pragma once
/// \file
/// Geometric routines templated on a kernel K, which provides the number type
/// K::FT and the orientation predicate K::orient. All decisions go through
/// K::orient and exact comparisons of stored coordinates, so they are exact in
/// ExactKernel and DoubleKernel alike.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>

#include "pd/geometry_types.hpp"

namespace pd {

// ---------------------------------------------------------------- basics

inline std::size_t next_index(std::size_t i, std::size_t n) { return i + 1 == n ? 0 : i + 1; }
inline std::size_t prev_index(std::size_t i, std::size_t n) { return i == 0 ? n - 1 : i - 1; }

/// Lexicographic order (x, then y).
template <class FT>
bool lex_less(const Point<FT>& a, const Point<FT>& b)
{
  return a.x < b.x || (a.x == b.x && a.y < b.y);
}

/// Twice the signed area (shoelace); positive for counter-clockwise rings.
/// Exact in ExactKernel; rounded in the double kernels.
template <class FT>
FT signed_area2(const Ring<FT>& r)
{
  FT sum = 0;
  const std::size_t n = r.size();
  for (std::size_t i = 0; i < n; ++i) {
    const Point<FT>& p = r[i];
    const Point<FT>& q = r[next_index(i, n)];
    FT term = p.x * q.y;
    term -= p.y * q.x;
    sum += term;
  }
  return sum;
}

/// Orientation of a simple ring: +1 counter-clockwise, -1 clockwise, 0 if the
/// ring is degenerate. Decided exactly by the turn at the lexicographically
/// smallest vertex, which is a strictly convex vertex of any simple ring.
template <class K>
int ring_orientation(const Ring<typename K::FT>& r)
{
  const std::size_t n = r.size();
  if (n < 3) return 0;
  std::size_t m = 0;
  for (std::size_t i = 1; i < n; ++i)
    if (lex_less(r[i], r[m])) m = i;
  return K::orient(r[prev_index(m, n)], r[m], r[next_index(m, n)]);
}

/// Reverses r in place if needed so that it is counter-clockwise (or clockwise).
template <class K>
void make_ccw(Ring<typename K::FT>& r)
{
  if (ring_orientation<K>(r) < 0) std::reverse(r.begin(), r.end());
}
template <class K>
void make_cw(Ring<typename K::FT>& r)
{
  if (ring_orientation<K>(r) > 0) std::reverse(r.begin(), r.end());
}

// ---------------------------------------------------------------- segments

/// True if p lies in the axis-aligned bounding box of segment ab. Combined
/// with orient(a, b, p) == 0 this means p lies on the closed segment ab.
template <class FT>
bool in_box(const Point<FT>& a, const Point<FT>& b, const Point<FT>& p)
{
  const bool x_ok = (a.x <= b.x) ? (a.x <= p.x && p.x <= b.x) : (b.x <= p.x && p.x <= a.x);
  const bool y_ok = (a.y <= b.y) ? (a.y <= p.y && p.y <= b.y) : (b.y <= p.y && p.y <= a.y);
  return x_ok && y_ok;
}

/// True if p lies on the closed segment ab.
template <class K>
bool on_segment(const Point<typename K::FT>& a, const Point<typename K::FT>& b,
                const Point<typename K::FT>& p)
{
  return K::orient(a, b, p) == 0 && in_box(a, b, p);
}

/// True if the closed segments ab and cd share at least one point.
template <class K>
bool segments_touch(const Point<typename K::FT>& a, const Point<typename K::FT>& b,
                    const Point<typename K::FT>& c, const Point<typename K::FT>& d)
{
  const int o1 = K::orient(a, b, c);
  const int o2 = K::orient(a, b, d);
  const int o3 = K::orient(c, d, a);
  const int o4 = K::orient(c, d, b);
  if (o1 * o2 < 0 && o3 * o4 < 0) return true;  // proper crossing
  if (o1 == 0 && in_box(a, b, c)) return true;
  if (o2 == 0 && in_box(a, b, d)) return true;
  if (o3 == 0 && in_box(c, d, a)) return true;
  if (o4 == 0 && in_box(c, d, b)) return true;
  return false;
}

// ---------------------------------------------------------------- point location

enum class Location { Outside, Inside, OnBoundary };

/// Exact location of q relative to the region bounded by the ring r (either
/// orientation). Uses the crossing-number rule with half-open edges in y, so
/// rays through vertices and horizontal edges are handled without any
/// general-position assumption.
template <class K>
Location locate(const Point<typename K::FT>& q, const Ring<typename K::FT>& r)
{
  const std::size_t n = r.size();
  bool inside = false;
  for (std::size_t i = 0; i < n; ++i) {
    const Point<typename K::FT>& a = r[i];
    const Point<typename K::FT>& b = r[next_index(i, n)];
    const int o = K::orient(a, b, q);
    if (o == 0 && in_box(a, b, q)) return Location::OnBoundary;
    // Edge straddles the horizontal line through q (half-open in y).
    if ((a.y > q.y) != (b.y > q.y)) {
      // The crossing lies to the right of q iff q is left of an upward edge
      // or right of a downward edge.
      const bool upward = b.y > a.y;
      if ((upward && o > 0) || (!upward && o < 0)) inside = !inside;
    }
  }
  return inside ? Location::Inside : Location::Outside;
}

// ---------------------------------------------------------------- conversion

/// Converts a double ring to kernel KTo's number type (exactly, for the exact
/// kernel). Throws std::domain_error on a non-finite coordinate.
template <class KTo>
Ring<typename KTo::FT> convert_ring(const Ring<double>& r)
{
  Ring<typename KTo::FT> out;
  out.reserve(r.size());
  for (const Point<double>& p : r) {
    const double x = p.x, y = p.y;
    if (!std::isfinite(x) || !std::isfinite(y))
      throw std::domain_error("polydecomp: non-finite coordinate");
    out.push_back({KTo::from_double(x), KTo::from_double(y)});
  }
  return out;
}

template <class KTo>
PolygonWithHoles<typename KTo::FT> convert_polygon(const PolygonWithHoles<double>& p)
{
  PolygonWithHoles<typename KTo::FT> out;
  out.outer = convert_ring<KTo>(p.outer);
  for (const Ring<double>& h : p.holes) out.holes.push_back(convert_ring<KTo>(h));
  return out;
}

/// Rounds a ring of any kernel to doubles (nearest).
template <class K>
Ring<double> to_double_ring(const Ring<typename K::FT>& r)
{
  Ring<double> out;
  out.reserve(r.size());
  for (const auto& p : r) out.push_back({K::to_double(p.x), K::to_double(p.y)});
  return out;
}

} // namespace pd
