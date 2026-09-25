// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#pragma once
/// \file
/// Exact validation of the input class: a simple outer ring and simple holes
/// that lie inside it, pairwise disjoint (no two rings touch) and not nested.

#include <algorithm>
#include <cstddef>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "pd/geometry.hpp"

namespace pd {

enum class IssueKind {
  TooFewVertices,      ///< ring_a has fewer than 3 vertices
  RepeatedVertex,      ///< ring_a has two equal consecutive vertices at edge_a
  EdgesIntersect,      ///< edge_a of ring_a touches edge_b of ring_b
  HoleNotInsideOuter,  ///< hole ring_a is not inside the outer ring
  HoleInsideHole,      ///< hole ring_a lies inside hole ring_b
  CoordinateOutOfRange,  ///< vertex edge_a of ring_a is outside the kernel's safe range
};

/// Ring index 0 is the outer ring; ring k >= 1 is hole k - 1. Edge i of a ring
/// runs from vertex i to vertex i + 1 (cyclically).
struct Issue {
  IssueKind kind;
  int ring_a = -1, edge_a = -1;
  int ring_b = -1, edge_b = -1;
};

std::string describe(const Issue& issue);

struct ValidationReport {
  std::vector<Issue> issues;
  bool ok() const { return issues.empty(); }
};

struct EdgeRef {
  int ring;
  int edge;
};
inline bool operator<(const EdgeRef& a, const EdgeRef& b)
{
  return std::tie(a.ring, a.edge) < std::tie(b.ring, b.edge);
}
inline bool operator==(const EdgeRef& a, const EdgeRef& b)
{
  return a.ring == b.ring && a.edge == b.edge;
}
using EdgePair = std::pair<EdgeRef, EdgeRef>;  // first < second

namespace detail {

template <class FT>
const Ring<FT>& ring_at(const PolygonWithHoles<FT>& p, int k)
{
  return k == 0 ? p.outer : p.holes[static_cast<std::size_t>(k - 1)];
}

/// True if two distinct edges conflict: edges of different rings, or
/// non-adjacent edges of one ring, must not touch; adjacent edges of one ring
/// may share only their common vertex (no fold-back along a common line).
/// Assumes no ring has repeated consecutive vertices.
template <class K>
bool edges_conflict(const PolygonWithHoles<typename K::FT>& p, EdgeRef e, EdgeRef f)
{
  using FT = typename K::FT;
  const Ring<FT>& re = ring_at(p, e.ring);
  const Ring<FT>& rf = ring_at(p, f.ring);
  const std::size_t ne = re.size(), nf = rf.size();
  const auto ei = static_cast<std::size_t>(e.edge), fi = static_cast<std::size_t>(f.edge);
  const Point<FT>& a = re[ei];
  const Point<FT>& b = re[next_index(ei, ne)];
  const Point<FT>& c = rf[fi];
  const Point<FT>& d = rf[next_index(fi, nf)];
  if (e.ring == f.ring) {
    // Edge ei ends where edge fi starts (or vice versa): they share one vertex.
    const bool e_then_f = next_index(ei, ne) == fi;
    const bool f_then_e = next_index(fi, nf) == ei;
    if (e_then_f || f_then_e) {
      // Shared vertex v with outer endpoints u (on e or f) and w (on the other).
      const Point<FT>& u = e_then_f ? a : c;
      const Point<FT>& v = e_then_f ? b : d;
      const Point<FT>& w = e_then_f ? d : b;
      if (e_then_f && f_then_e) return true;  // two-vertex ring; excluded earlier
      // They overlap beyond v iff collinear and u, w on the same side of v.
      return K::orient(u, v, w) == 0 && (in_box(v, u, w) || in_box(v, w, u));
    }
  }
  return segments_touch<K>(a, b, c, d);
}

template <class FT>
struct EdgeBox {
  FT xmin, xmax, ymin, ymax;
  EdgeRef ref;
};

} // namespace detail

/// All conflicting edge pairs (see detail::edges_conflict), sorted. The sweep
/// sorts edges by their minimum x and tests only pairs whose bounding boxes
/// overlap; with use_sweep = false every pair is tested (for cross-checking).
/// Stops after max_pairs pairs if max_pairs > 0.
template <class K>
std::vector<EdgePair> conflicting_edge_pairs(const PolygonWithHoles<typename K::FT>& p,
                                             bool use_sweep = true, std::size_t max_pairs = 0)
{
  using FT = typename K::FT;
  std::vector<detail::EdgeBox<FT>> boxes;
  const int nrings = 1 + static_cast<int>(p.holes.size());
  for (int k = 0; k < nrings; ++k) {
    const Ring<FT>& r = detail::ring_at(p, k);
    const std::size_t n = r.size();
    for (std::size_t i = 0; i < n; ++i) {
      const Point<FT>& a = r[i];
      const Point<FT>& b = r[next_index(i, n)];
      detail::EdgeBox<FT> box{a.x < b.x ? a.x : b.x, a.x < b.x ? b.x : a.x,
                              a.y < b.y ? a.y : b.y, a.y < b.y ? b.y : a.y,
                              EdgeRef{k, static_cast<int>(i)}};
      boxes.push_back(std::move(box));
    }
  }

  std::vector<EdgePair> pairs;
  auto record = [&](EdgeRef e, EdgeRef f) {
    if (f < e) std::swap(e, f);
    pairs.emplace_back(e, f);
    return max_pairs > 0 && pairs.size() >= max_pairs;
  };

  // Returns early once max_pairs pairs have been recorded.
  auto scan = [&]() {
    if (!use_sweep) {
      for (std::size_t i = 0; i < boxes.size(); ++i)
        for (std::size_t j = i + 1; j < boxes.size(); ++j)
          if (detail::edges_conflict<K>(p, boxes[i].ref, boxes[j].ref) &&
              record(boxes[i].ref, boxes[j].ref))
            return;
      return;
    }
    std::vector<std::size_t> order(boxes.size());
    for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(),
              [&](std::size_t i, std::size_t j) { return boxes[i].xmin < boxes[j].xmin; });
    std::vector<std::size_t> active;
    for (std::size_t idx : order) {
      const detail::EdgeBox<FT>& cur = boxes[idx];
      active.erase(std::remove_if(active.begin(), active.end(),
                                  [&](std::size_t k) { return boxes[k].xmax < cur.xmin; }),
                   active.end());
      for (std::size_t k : active) {
        const detail::EdgeBox<FT>& other = boxes[k];
        if (other.ymin <= cur.ymax && cur.ymin <= other.ymax &&
            detail::edges_conflict<K>(p, other.ref, cur.ref) && record(other.ref, cur.ref))
          return;
      }
      active.push_back(idx);
    }
  };
  scan();
  std::sort(pairs.begin(), pairs.end());
  return pairs;
}

struct ValidateOptions {
  bool use_sweep = true;
  std::size_t max_issues = 100;  ///< stop collecting after this many
};

/// Checks that p is a valid complex polygon: every ring has at least three
/// vertices and no repeated consecutive vertex, no two edges conflict (so every
/// ring is simple and the rings are pairwise disjoint), every hole lies inside
/// the outer ring, and no hole lies inside another. Ring orientation is not
/// checked; the algorithms normalise it.
template <class K>
ValidationReport validate(const PolygonWithHoles<typename K::FT>& p, const ValidateOptions& opt = {})
{
  using FT = typename K::FT;
  ValidationReport rep;
  const int nrings = 1 + static_cast<int>(p.holes.size());

  // 0. Coordinates the kernel can handle (floating-point kernels have a range).
  for (int k = 0; k < nrings; ++k) {
    const Ring<FT>& r = detail::ring_at(p, k);
    for (std::size_t i = 0; i < r.size(); ++i)
      if (!K::coordinate_ok(r[i].x) || !K::coordinate_ok(r[i].y))
        rep.issues.push_back({IssueKind::CoordinateOutOfRange, k, static_cast<int>(i)});
  }
  if (!rep.ok()) return rep;

  // 1. Per-ring structure. Later checks assume it, so stop here on failure.
  for (int k = 0; k < nrings; ++k) {
    const Ring<FT>& r = detail::ring_at(p, k);
    if (r.size() < 3) {
      rep.issues.push_back({IssueKind::TooFewVertices, k});
      continue;
    }
    for (std::size_t i = 0; i < r.size(); ++i)
      if (r[i] == r[next_index(i, r.size())])
        rep.issues.push_back({IssueKind::RepeatedVertex, k, static_cast<int>(i)});
  }
  if (!rep.ok()) return rep;

  // 2. No conflicting edges: rings simple and pairwise disjoint.
  for (const EdgePair& e : conflicting_edge_pairs<K>(p, opt.use_sweep, opt.max_issues))
    rep.issues.push_back(
        {IssueKind::EdgesIntersect, e.first.ring, e.first.edge, e.second.ring, e.second.edge});
  if (!rep.ok()) return rep;

  // 3. Containment. Boundaries are disjoint, so one vertex decides for a whole ring.
  for (int k = 1; k < nrings; ++k)
    if (locate<K>(detail::ring_at(p, k).front(), p.outer) != Location::Inside)
      rep.issues.push_back({IssueKind::HoleNotInsideOuter, k});
  for (int a = 1; a < nrings; ++a)
    for (int b = 1; b < nrings; ++b)
      if (a != b && locate<K>(detail::ring_at(p, a).front(), detail::ring_at(p, b)) ==
                        Location::Inside)
        rep.issues.push_back({IssueKind::HoleInsideHole, a, -1, b});
  return rep;
}

} // namespace pd
