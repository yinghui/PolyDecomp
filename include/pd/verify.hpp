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
///
/// With weakly_simple = true, (1) is relaxed to weakly simple rings: a ring may
/// pass through one point more than once (as a hole-free piece does after a
/// dangling contact), provided its edges never cross or overlap and its
/// interior angles at that point are disjoint. Winding numbers are still 1
/// inside and 0 outside such a ring, so the conclusion is unchanged.

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

/// True if the direction x - v lies strictly inside the counter-clockwise arc
/// of directions from u - v to w - v.
template <class K>
bool in_open_arc(const Point<typename K::FT>& v, const Point<typename K::FT>& u,
                 const Point<typename K::FT>& w, const Point<typename K::FT>& x)
{
  const int uw = K::orient(v, u, w);
  if (uw > 0) return K::orient(v, u, x) > 0 && K::orient(v, x, w) > 0;
  if (uw < 0) return !(K::orient(v, w, x) >= 0 && K::orient(v, x, u) >= 0);
  return K::orient(v, u, x) > 0;  // u, w opposite: a half-plane
}

/// Weak simplicity of a counter-clockwise ring without repeated consecutive
/// vertices: non-adjacent edges may only meet at a common endpoint, without
/// overlapping, and at a point visited twice the interior angles are disjoint.
template <class K>
bool weakly_simple(const Ring<typename K::FT>& r)
{
  using FT = typename K::FT;
  const std::size_t n = r.size();
  for (const EdgePair& pr : conflicting_edge_pairs<K>(PolygonWithHoles<FT>{r, {}}, true, 0)) {
    const std::size_t e = static_cast<std::size_t>(pr.first.edge), f = static_cast<std::size_t>(pr.second.edge);
    const Point<FT>& a = r[e];
    const Point<FT>& b = r[next_index(e, n)];
    const Point<FT>& c = r[f];
    const Point<FT>& d = r[next_index(f, n)];
    const Point<FT>* shared = nullptr;
    const Point<FT>*u = nullptr, *w = nullptr;
    if (a == c) shared = &a, u = &b, w = &d;
    else if (a == d) shared = &a, u = &b, w = &c;
    else if (b == c) shared = &b, u = &a, w = &d;
    else if (b == d) shared = &b, u = &a, w = &c;
    if (!shared) return false;  // they cross or touch elsewhere
    if (K::orient(*u, *shared, *w) == 0 && (in_box(*shared, *u, *w) || in_box(*shared, *w, *u)))
      return false;  // collinear overlap
  }
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = i + 1; j < n; ++j) {
      if (!(r[i] == r[j])) continue;
      const Point<FT>& v = r[i];
      const Point<FT>& ai = r[prev_index(i, n)];
      const Point<FT>& bi = r[next_index(i, n)];
      const Point<FT>& aj = r[prev_index(j, n)];
      const Point<FT>& bj = r[next_index(j, n)];
      // Interior at an occurrence: the CCW arc from its next to its previous vertex.
      if (in_open_arc<K>(v, bi, ai, aj) || in_open_arc<K>(v, bi, ai, bj) ||
          in_open_arc<K>(v, bj, aj, ai) || in_open_arc<K>(v, bj, aj, bi))
        return false;
    }
  return true;
}

} // namespace detail

template <class K>
VerifyReport verify_partition(const PolygonWithHoles<typename K::FT>& input,
                              const std::vector<Ring<typename K::FT>>& pieces, bool require_convex,
                              bool weakly_simple = false)
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
    if (weakly_simple ? !detail::weakly_simple<K>(r)
                      : !conflicting_edge_pairs<K>(PolygonWithHoles<FT>{r, {}}, true, 1).empty())
      return fail(name + (weakly_simple ? " is not weakly simple" : " is not simple"));
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

/// Certificate for pieces that may still contain holes, such as those of a Slicer. It accepts iff
///   (1) every hole of a piece is a hole of the input (the same vertices), kept by exactly one piece,
///       has a vertex strictly inside that piece's outer ring, and has no edge that meets an edge of
///       that ring, and
///   (2) the outer rings of the pieces tile the input with the kept holes filled in
///       (verify_partition, weakly simple rings allowed if asked).
/// Then the pieces, each its outer ring minus its holes, tile the input: by (2) the outer rings have
/// pairwise disjoint interiors that cover the region with the kept holes filled; by (1) the boundary
/// of a kept hole avoids its piece's outer ring and has a point inside it, so the whole hole lies
/// inside that ring (and in no other piece), and removing it from that piece leaves a tiling of the
/// input. Holes are matched in O(h log h) comparisons of vertex lists; the test of (1) compares the
/// edges of each kept hole with the edges of its piece's outer ring whose boxes meet the hole's box.
template <class K>
VerifyReport verify_slices(const PolygonWithHoles<typename K::FT>& input,
                           const std::vector<PolygonWithHoles<typename K::FT>>& pieces,
                           bool weakly_simple = false)
{
  using FT = typename K::FT;
  using Key = std::vector<Point<FT>>;
  VerifyReport rep;
  auto fail = [&](std::string msg) {
    rep.ok = false;
    rep.failure = std::move(msg);
    return rep;
  };
  auto less_point = [](const Point<FT>& a, const Point<FT>& b) { return lex_less(a, b); };
  auto key = [&](const Ring<FT>& r) {
    Key v(r.begin(), r.end());
    std::sort(v.begin(), v.end(), less_point);
    return v;
  };
  auto less_key = [&](const Key& a, const Key& b) {
    if (a.size() != b.size()) return a.size() < b.size();
    return std::lexicographical_compare(a.begin(), a.end(), b.begin(), b.end(), less_point);
  };
  std::map<Key, std::size_t, decltype(less_key)> input_keys(less_key);
  for (std::size_t j = 0; j < input.holes.size(); ++j) input_keys.emplace(key(input.holes[j]), j);
  std::vector<bool> kept(input.holes.size(), false);
  std::vector<Ring<FT>> outers;
  for (std::size_t k = 0; k < pieces.size(); ++k) {
    const PolygonWithHoles<FT>& p = pieces[k];
    for (const Ring<FT>& h : p.holes) {
      const auto it = h.empty() ? input_keys.end() : input_keys.find(key(h));
      if (it == input_keys.end()) return fail("piece " + std::to_string(k) + " has a hole that is not an input hole");
      if (kept[it->second]) return fail("input hole " + std::to_string(it->second) + " is a hole of two pieces");
      kept[it->second] = true;
      if (locate<K>(h.front(), p.outer) != Location::Inside)
        return fail("a hole of piece " + std::to_string(k) + " is not inside it");
      // No edge of the hole meets an edge of the outer ring (closed segments, exact).
      FT x0 = h[0].x, x1 = x0, y0 = h[0].y, y1 = y0;
      for (const auto& v : h) {
        if (v.x < x0) x0 = v.x;
        if (x1 < v.x) x1 = v.x;
        if (v.y < y0) y0 = v.y;
        if (y1 < v.y) y1 = v.y;
      }
      const std::size_t m = p.outer.size();
      for (std::size_t i = 0; i < m; ++i) {
        const Point<FT>& a = p.outer[i];
        const Point<FT>& b = p.outer[(i + 1) % m];
        if ((a.x < x0 && b.x < x0) || (x1 < a.x && x1 < b.x) || (a.y < y0 && b.y < y0) || (y1 < a.y && y1 < b.y))
          continue;  // the edge's box misses the hole's box
        for (std::size_t j = 0; j < h.size(); ++j)
          if (segments_touch<K>(a, b, h[j], h[(j + 1) % h.size()]))
            return fail("a hole of piece " + std::to_string(k) + " meets its outer ring");
      }
    }
    outers.push_back(p.outer);
  }
  PolygonWithHoles<FT> filled{input.outer, {}};
  for (std::size_t j = 0; j < input.holes.size(); ++j)
    if (!kept[j]) filled.holes.push_back(input.holes[j]);
  return verify_partition<K>(filled, outers, false, weakly_simple);
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
