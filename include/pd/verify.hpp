// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#pragma once
/// \file
/// Exact certificate that a set of rings partitions a polygon with holes.
///
/// verify_partition accepts the pieces iff
///   (1) every piece is a simple, counter-clockwise ring (convex, if required),
///   (2) the pieces' boundaries sum, as oriented 1-chains, to the boundary of
///       the input (outer ring counter-clockwise, holes clockwise).
/// (2) is checked line by line: every edge lies on a line, and along each line
/// the signed coverage by piece edges must equal that by input edges on every
/// sub-interval. Since winding numbers add over sums of chains, (2) gives, for
/// every point x off the edges, sum_j wind(piece_j, x) = wind(input, x), which
/// is 1 inside the input region and 0 elsewhere. By (1) each term is 0 or 1, so
/// every such point lies in exactly one piece if it is in the region and in
/// none otherwise: the pieces tile the region with disjoint interiors. T-
/// junctions (a vertex of one piece inside an edge of another) are allowed.
///
/// The certificate is exact with ExactKernel. With a double kernel the chain
/// test needs exact collinearity of rounded points and is only indicative.

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "pd/geometry.hpp"
#include "pd/validate.hpp"

namespace pd {

struct VerifyReport {
  bool ok = true;
  std::string failure;  ///< first failed condition when !ok
};

namespace detail {

template <class FT>
struct LineKey {
  FT a, b, c;  ///< normalised line a x + b y = c, with a == 1, or a == 0 and b == 1
  bool operator<(const LineKey& o) const
  {
    if (a != o.a) return a < o.a;
    if (b != o.b) return b < o.b;
    return c < o.c;
  }
};

} // namespace detail

template <class K>
VerifyReport verify_partition(const PolygonWithHoles<typename K::FT>& input,
                              const std::vector<Ring<typename K::FT>>& pieces, bool require_convex)
{
  using FT = typename K::FT;
  VerifyReport rep;
  auto fail = [&](std::string msg) {
    rep.ok = false;
    rep.failure = std::move(msg);
    return rep;
  };

  // (1) each piece simple and counter-clockwise (and convex if required)
  for (std::size_t k = 0; k < pieces.size(); ++k) {
    const Ring<FT>& r = pieces[k];
    const std::string name = "piece " + std::to_string(k);
    if (r.size() < 3) return fail(name + " has fewer than 3 vertices");
    for (std::size_t i = 0; i < r.size(); ++i)
      if (r[i] == r[next_index(i, r.size())]) return fail(name + " repeats a vertex");
    if (!conflicting_edge_pairs<K>(PolygonWithHoles<FT>{r, {}}, true, 1).empty())
      return fail(name + " is not simple");
    if (ring_orientation<K>(r) <= 0) return fail(name + " is not counter-clockwise");
    if (require_convex)
      for (std::size_t i = 0; i < r.size(); ++i)
        if (K::orient(r[prev_index(i, r.size())], r[i], r[next_index(i, r.size())]) < 0)
          return fail(name + " is not convex at vertex " + std::to_string(i));
  }

  // (2) boundary chains agree, line by line
  std::map<detail::LineKey<FT>, std::vector<std::pair<FT, int>>> events;
  auto add_edge = [&](const Point<FT>& p, const Point<FT>& q, int mult) {
    FT a = q.y - p.y, b = p.x - q.x;
    FT c = a * p.x + b * p.y;
    if (a != 0) {
      b /= a;
      c /= a;
      a = 1;
    } else {
      c /= b;
      b = 1;
    }
    const bool vertical = (b == 0);  // then parametrise by y, otherwise by x
    const FT& s0 = vertical ? p.y : p.x;
    const FT& s1 = vertical ? q.y : q.x;
    std::vector<std::pair<FT, int>>& ev = events[detail::LineKey<FT>{a, b, c}];
    if (s0 < s1) {  // runs along the line's increasing direction
      ev.emplace_back(s0, mult);
      ev.emplace_back(s1, -mult);
    } else {
      ev.emplace_back(s1, -mult);
      ev.emplace_back(s0, mult);
    }
  };
  auto add_ring = [&](const Ring<FT>& r, int mult) {
    for (std::size_t i = 0; i < r.size(); ++i) add_edge(r[i], r[next_index(i, r.size())], mult);
  };
  for (const Ring<FT>& r : pieces) add_ring(r, +1);
  Ring<FT> outer = input.outer;
  make_ccw<K>(outer);
  add_ring(outer, -1);
  for (const Ring<FT>& h : input.holes) {
    Ring<FT> hole = h;
    make_cw<K>(hole);
    add_ring(hole, -1);
  }
  for (auto& entry : events) {
    std::vector<std::pair<FT, int>>& ev = entry.second;
    std::sort(ev.begin(), ev.end(),
              [](const std::pair<FT, int>& x, const std::pair<FT, int>& y) { return x.first < y.first; });
    int running = 0;
    for (std::size_t i = 0; i < ev.size();) {
      std::size_t j = i;
      while (j < ev.size() && ev[j].first == ev[i].first) running += ev[j++].second;
      if (running != 0)
        return fail("piece boundaries do not match the input boundary along a line (near " +
                    std::to_string(K::to_double(ev[i].first)) + ")");
      i = j;
    }
  }
  return rep;
}

/// Necessary conditions for a floating-point result, checked in its own
/// arithmetic: every piece is simple and counter-clockwise (exact tests on the
/// stored coordinates), convex up to angle_tol (radians) if required, and the
/// piece areas sum to the input's area within a relative area_tol. Unlike
/// verify_partition this does not prove the pieces tile the input; it catches
/// the failures rounding typically causes (self-touching or folded pieces,
/// lost or doubled area).
template <class K>
VerifyReport check_pieces_approx(const PolygonWithHoles<typename K::FT>& input,
                                 const std::vector<Ring<typename K::FT>>& pieces, bool require_convex,
                                 double angle_tol, double area_tol)
{
  using FT = typename K::FT;
  VerifyReport rep;
  auto fail = [&](std::string msg) {
    rep.ok = false;
    rep.failure = std::move(msg);
    return rep;
  };
  double area = 0;
  for (std::size_t k = 0; k < pieces.size(); ++k) {
    const Ring<FT>& r = pieces[k];
    const std::string name = "piece " + std::to_string(k);
    if (!validate<K>(PolygonWithHoles<FT>{r, {}}).ok()) return fail(name + " is not simple");
    if (ring_orientation<K>(r) <= 0) return fail(name + " is not counter-clockwise");
    if (require_convex)
      for (std::size_t i = 0; i < r.size(); ++i) {
        const Point<FT>& a = r[prev_index(i, r.size())];
        const Point<FT>& b = r[i];
        const Point<FT>& c = r[next_index(i, r.size())];
        if (K::orient(a, b, c) >= 0) continue;
        const double ux = K::to_double(b.x) - K::to_double(a.x), uy = K::to_double(b.y) - K::to_double(a.y);
        const double vx = K::to_double(c.x) - K::to_double(b.x), vy = K::to_double(c.y) - K::to_double(b.y);
        if (!(ux * vx + uy * vy > 0 &&
              std::abs(ux * vy - uy * vx) <= angle_tol * std::hypot(ux, uy) * std::hypot(vx, vy)))
          return fail(name + " is not convex at vertex " + std::to_string(i));
      }
    area += K::to_double(signed_area2(r));
  }
  double expected = std::abs(K::to_double(signed_area2(input.outer)));
  for (const Ring<FT>& h : input.holes) expected -= std::abs(K::to_double(signed_area2(h)));
  if (std::abs(area - expected) > area_tol * std::abs(expected))
    return fail("piece areas do not sum to the input area");
  return rep;
}

} // namespace pd
