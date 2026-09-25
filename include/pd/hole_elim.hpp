// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#pragma once
/// \file
/// Stage 1: hole elimination by line cuts.
///
/// A node is an outer ring (counter-clockwise) with the input holes it contains.
/// Each cut is a line that pierces at least one of the node's holes:
///  - one hole: the line through the anchor a = (v_{m-1} + v_m + v_{m+1}) / 3,
///    v_m being the outer vertex nearest to the hole's vertex-average centroid
///    c, and c itself;
///  - two or more holes: the line through the centroids of the two holes whose
///    centroids are farthest apart.
/// If that preferred line passes exactly through a vertex of the node, or has
/// no direction (coincident centroids, a == c), it is replaced by a line through
/// a point q of one of those holes: its centroid, or, if the centroid is itself
/// a vertex, the midpoint of the hole's longest edge. The direction is taken in
/// the widest angular gap between the directions from q to the node's vertices.
/// A cut line therefore never passes through a vertex. It crosses the ring of
/// the hole through q at least twice (q is interior to the hole's convex hull,
/// or a transverse crossing point), so every cut pierces a hole and at most h
/// cuts occur, with no general-position assumption.
///
/// The line meets the rings at crossings x_1 < ... < x_m (m even) along its
/// direction; the region occupies the segments (x_{2i-1}, x_{2i}). The pieces on
/// the left of the line are traced by walking a segment x_{2i-1} -> x_{2i} and
/// then the ring arc that leaves x_{2i} to the left, up to the next crossing on
/// that ring; each closed walk is one simple piece. The pieces on the right are
/// traced by the same walk on the reversed line. Pierced holes become notches of
/// the pieces; the other holes go to the piece that contains them.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "pd/detail/lines.hpp"
#include "pd/geometry.hpp"
#include "pd/status.hpp"
#include "pd/validate.hpp"

namespace pd {

enum class CutFallback {
  None,          ///< the preferred line was used
  Rotated,       ///< a line through the hole centroid in a vertex-free direction
  EdgeMidpoint,  ///< the centroid is a vertex: a line through an edge midpoint
};

struct HoleElimOptions {
  bool validate_input = true;  ///< check the input class first
  bool check = false;          ///< extra invariant checks (piece orientation, unique hole owner)
  bool keep_tree = true;       ///< record the decomposition tree
};

struct HoleElimStats {
  std::size_t holes = 0;
  std::size_t cuts = 0;
  std::size_t max_depth = 0;
  std::size_t max_crossings = 0;          ///< most crossings of a single cut
  std::size_t fallback_rotated = 0;       ///< cuts that used CutFallback::Rotated
  std::size_t fallback_midpoint = 0;      ///< cuts that used CutFallback::EdgeMidpoint
  std::size_t degenerate_directions = 0;  ///< preferred line had no direction
};

template <class FT>
struct HoleNode {
  int parent = -1;
  std::vector<int> children;
  int leaf = -1;              ///< index into leaves, for a node without holes
  std::vector<int> holes;     ///< input hole indices contained in this node
  bool has_cut = false;
  Point<FT> cut_point{}, cut_dir{};  ///< the cut line: cut_point + t cut_dir
  CutFallback fallback = CutFallback::None;
  std::vector<int> pierced;   ///< input hole indices dissolved by the cut
};

template <class FT>
struct HoleElimination {
  std::vector<Point<FT>> pool;                     ///< input vertices (outer, then holes), then crossings
  std::vector<std::vector<std::uint32_t>> leaves;  ///< hole-free pieces, CCW, as pool ids
  std::vector<HoleNode<FT>> tree;                  ///< tree[0] is the root (if kept)
  HoleElimStats stats;
  Diagnostics diag;

  bool ok() const { return diag.status == Status::Ok; }
  Ring<FT> leaf_ring(std::size_t i) const
  {
    Ring<FT> r;
    r.reserve(leaves[i].size());
    for (std::uint32_t id : leaves[i]) r.push_back(pool[id]);
    return r;
  }
  std::vector<Ring<FT>> leaf_rings() const
  {
    std::vector<Ring<FT>> out;
    out.reserve(leaves.size());
    for (std::size_t i = 0; i < leaves.size(); ++i) out.push_back(leaf_ring(i));
    return out;
  }
};

namespace detail {

struct Corner1 {
  std::uint32_t v;     ///< pool id
  std::uint32_t line;  ///< line containing the edge from this vertex to the next
};

template <class K>
class HoleEliminator {
public:
  using FT = typename K::FT;

  HoleEliminator(HoleElimination<FT>& out, const HoleElimOptions& opt)
      : out_(out), opt_(opt),
        broken_(is_exact_kernel<K>::value ? Status::InternalError : Status::NumericalFailure)
  {}

  void run(PolygonWithHoles<FT> in);

private:
  struct Cut {
    Point<FT> q, d;
    CutFallback fallback = CutFallback::None;
    bool degenerate_direction = false;
  };
  struct Crossing {
    FT t;               // position along the cut direction
    std::uint32_t v;    // pool id of the crossing point
    int ring;           // -1: the node's outer ring; otherwise an input hole index
    std::uint32_t edge; // edge of that ring containing the crossing
  };
  struct Work {
    std::vector<Corner1> outer;
    std::vector<int> holes;
    int node;
    std::size_t depth;
  };

  HoleElimination<FT>& out_;
  const HoleElimOptions& opt_;
  const Status broken_;
  std::vector<Line<FT>> lines_;
  std::vector<std::uint32_t> hole_off_, hole_len_;
  double scale_ = 1;  ///< largest absolute input coordinate (at least 1)
  double tol_ = 0;    ///< K::clearance * scale_ (0: exact incidence only)

  /// Distance from w to the line through q and q2, in double precision.
  double distance_to_line(const Point<FT>& q, const Point<FT>& q2, const Point<FT>& w) const
  {
    const double qx = K::to_double(q.x), qy = K::to_double(q.y);
    const double dx = K::to_double(q2.x) - qx, dy = K::to_double(q2.y) - qy;
    return std::abs(dx * (K::to_double(w.y) - qy) - dy * (K::to_double(w.x) - qx)) / std::hypot(dx, dy);
  }

  bool fail(Status s, std::string msg)
  {
    out_.diag = {s, std::move(msg), {}};
    return false;
  }
  const Point<FT>& pt(std::uint32_t id) const { return out_.pool[id]; }

  template <class F>
  void for_each_vertex(const Work& w, F f) const
  {
    for (const Corner1& c : w.outer) f(c.v);
    for (int h : w.holes)
      for (std::uint32_t k = 0; k < hole_len_[h]; ++k) f(hole_off_[h] + k);
  }

  Point<FT> centroid(int h) const
  {
    FT sx = 0, sy = 0;
    for (std::uint32_t k = 0; k < hole_len_[h]; ++k) {
      sx += pt(hole_off_[h] + k).x;
      sy += pt(hole_off_[h] + k).y;
    }
    const FT n = static_cast<long>(hole_len_[h]);
    return {sx / n, sy / n};
  }

  /// True if no vertex of the node lies on the line through q and q2 (or, in a
  /// floating-point kernel, within tol_ of it).
  bool line_avoids_vertices(const Work& w, const Point<FT>& q, const Point<FT>& q2) const
  {
    bool ok = true;
    for_each_vertex(w, [&](std::uint32_t id) {
      ok = ok && K::orient(q, q2, pt(id)) != 0 && (tol_ == 0 || distance_to_line(q, q2, pt(id)) > tol_);
    });
    return ok;
  }

  /// True if q is a vertex of the node (or, in a floating-point kernel, within tol_ of one).
  bool is_vertex(const Work& w, const Point<FT>& q) const
  {
    bool hit = false;
    for_each_vertex(w, [&](std::uint32_t id) {
      const Point<FT>& v = pt(id);
      hit = hit || v == q ||
            (tol_ > 0 && std::hypot(K::to_double(v.x) - K::to_double(q.x),
                                    K::to_double(v.y) - K::to_double(q.y)) <= tol_);
    });
    return hit;
  }

  Cut choose_cut(const Work& w) const;
  Point<FT> vertex_free_direction(const Work& w, const Point<FT>& q, Point<FT> base) const;
  bool cut_node(Work& w, std::vector<Work>& stack);
  bool trace_side(int sign, std::uint32_t line_id, const Point<FT>& q, const Point<FT>& q2,
                  const std::vector<Crossing>& xs, const std::vector<std::size_t>& order,
                  const std::vector<std::size_t>& ring_next, const std::vector<std::size_t>& ring_prev,
                  const Work& w, std::vector<std::vector<Corner1>>& pieces);
};

template <class K>
typename HoleEliminator<K>::Cut HoleEliminator<K>::choose_cut(const Work& w) const
{
  Point<FT> q0, d0;
  int hole;
  if (w.holes.size() == 1) {
    hole = w.holes[0];
    const Point<FT> c = centroid(hole);
    // Outer vertex nearest to c (the last one in ring order on ties).
    const std::size_t m = w.outer.size();
    std::size_t best = 0;
    FT best_d2 = 0;
    for (std::size_t i = 0; i < m; ++i) {
      const Point<FT>& v = pt(w.outer[i].v);
      const FT dx = v.x - c.x, dy = v.y - c.y;
      const FT d2 = dx * dx + dy * dy;
      if (i == 0 || d2 <= best_d2) {
        best = i;
        best_d2 = d2;
      }
    }
    const Point<FT>& a0 = pt(w.outer[prev_index(best, m)].v);
    const Point<FT>& a1 = pt(w.outer[best].v);
    const Point<FT>& a2 = pt(w.outer[next_index(best, m)].v);
    q0 = {(a0.x + a1.x + a2.x) / 3, (a0.y + a1.y + a2.y) / 3};
    d0 = {c.x - q0.x, c.y - q0.y};
  } else {
    // Farthest pair of centroids (the last one in pair order on ties).
    std::vector<Point<FT>> cs;
    for (int h : w.holes) cs.push_back(centroid(h));
    std::size_t r = 0, s = 1;
    FT best_d2 = -1;
    for (std::size_t i = 0; i + 1 < cs.size(); ++i)
      for (std::size_t j = i + 1; j < cs.size(); ++j) {
        const FT dx = cs[j].x - cs[i].x, dy = cs[j].y - cs[i].y;
        const FT d2 = dx * dx + dy * dy;
        if (d2 >= best_d2) {
          r = i;
          s = j;
          best_d2 = d2;
        }
      }
    hole = w.holes[r];
    q0 = cs[r];
    d0 = {cs[s].x - cs[r].x, cs[s].y - cs[r].y};
  }

  Cut cut;
  cut.degenerate_direction = (d0.x == 0 && d0.y == 0);
  if (!cut.degenerate_direction &&
      line_avoids_vertices(w, q0, Point<FT>{q0.x + d0.x, q0.y + d0.y})) {
    cut.q = q0;
    cut.d = d0;
    return cut;
  }

  // Fallback: a vertex-free line through a point of `hole`.
  cut.q = centroid(hole);
  cut.fallback = CutFallback::Rotated;
  if (is_vertex(w, cut.q)) {
    // Midpoint of the hole's longest edge (the first on ties); never a vertex.
    std::uint32_t best = 0;
    FT best_len = -1;
    for (std::uint32_t k = 0; k < hole_len_[hole]; ++k) {
      const Point<FT>& a = pt(hole_off_[hole] + k);
      const Point<FT>& b = pt(hole_off_[hole] + (k + 1) % hole_len_[hole]);
      const FT dx = b.x - a.x, dy = b.y - a.y;
      const FT len = dx * dx + dy * dy;
      if (len > best_len) {
        best = k;
        best_len = len;
      }
    }
    const Point<FT>& a = pt(hole_off_[hole] + best);
    const Point<FT>& b = pt(hole_off_[hole] + (best + 1) % hole_len_[hole]);
    cut.q = {(a.x + b.x) / 2, (a.y + b.y) / 2};
    cut.fallback = CutFallback::EdgeMidpoint;
  }
  cut.d = vertex_free_direction(w, cut.q, cut.degenerate_direction ? Point<FT>{FT(1), FT(0)} : d0);
  return cut;
}

/// A direction base + s * perp(base), s in (-1, 1), such that the line through
/// q in that direction avoids every vertex of the node. The finitely many
/// forbidden values of s split (-1, 1) into gaps; the widest gap is tried first,
/// at its simplest dyadic point. In exact arithmetic that line always avoids the
/// vertices; with rounding, candidates are re-checked (including the clearance
/// tol_) and further points and gaps are tried. q must not be a vertex.
template <class K>
Point<typename K::FT> HoleEliminator<K>::vertex_free_direction(const Work& w, const Point<FT>& q,
                                                               Point<FT> base) const
{
  const Point<FT> nrm{-base.y, base.x};
  std::vector<FT> forbidden;
  for_each_vertex(w, [&](std::uint32_t id) {
    const FT ux = pt(id).x - q.x, uy = pt(id).y - q.y;
    const FT den = nrm.x * uy - nrm.y * ux;
    if (den == 0) return;  // no s makes the line pass through this vertex
    const FT s = -(base.x * uy - base.y * ux) / den;
    if (s > -1 && s < 1) forbidden.push_back(s);
  });
  std::sort(forbidden.begin(), forbidden.end());
  forbidden.erase(std::unique(forbidden.begin(), forbidden.end()), forbidden.end());
  std::vector<std::pair<FT, FT>> gaps;  // (lo, hi), widest first (stable: earlier first on ties)
  FT prev = -1;
  for (std::size_t i = 0; i <= forbidden.size(); ++i) {
    const FT next = i < forbidden.size() ? forbidden[i] : FT(1);
    gaps.emplace_back(prev, next);
    prev = next;
  }
  std::stable_sort(gaps.begin(), gaps.end(), [](const std::pair<FT, FT>& x, const std::pair<FT, FT>& y) {
    return x.second - x.first > y.second - y.first;
  });
  auto direction = [&](const FT& s) { return Point<FT>{base.x + s * nrm.x, base.y + s * nrm.y}; };
  for (const auto& g : gaps) {
    const FT s0 = dyadic_between(g.first, g.second);
    for (const FT& s : {s0, FT((g.first + s0) / 2), FT((s0 + g.second) / 2)}) {
      const Point<FT> d = direction(s);
      if (line_avoids_vertices(w, q, Point<FT>{q.x + d.x, q.y + d.y})) return d;
    }
  }
  return direction(dyadic_between(gaps.front().first, gaps.front().second));  // best effort
}

template <class K>
bool HoleEliminator<K>::trace_side(int sign, std::uint32_t line_id, const Point<FT>& q,
                                   const Point<FT>& q2, const std::vector<Crossing>& xs,
                                   const std::vector<std::size_t>& order,
                                   const std::vector<std::size_t>& ring_next,
                                   const std::vector<std::size_t>& ring_prev, const Work& w,
                                   std::vector<std::vector<Corner1>>& pieces)
{
  const std::size_t m = xs.size();
  std::vector<std::size_t> seq(order);
  if (sign < 0) std::reverse(seq.begin(), seq.end());
  std::vector<std::size_t> pos(m);
  for (std::size_t p = 0; p < m; ++p) pos[seq[p]] = p;
  std::vector<bool> used(m, false);

  auto ring_len = [&](int ring) -> std::uint32_t {
    return ring < 0 ? static_cast<std::uint32_t>(w.outer.size()) : hole_len_[ring];
  };
  auto vid = [&](int ring, std::uint32_t k) -> std::uint32_t {
    return ring < 0 ? w.outer[k].v : hole_off_[ring] + k;
  };
  auto lid = [&](int ring, std::uint32_t k) -> std::uint32_t {
    return ring < 0 ? w.outer[k].line : hole_off_[ring] + k;  // input hole edge k has line id = its start vertex id
  };

  for (std::size_t start = 0; start < m; start += 2) {
    if (used[seq[start]]) continue;
    std::vector<Corner1> piece;
    std::size_t p = start;
    for (std::size_t steps = 0;; ++steps) {
      if (steps > m) return fail(broken_, "a piece walk does not close");
      const std::size_t xa = seq[p], xb = seq[p + 1];
      if (used[xa] || used[xb]) return fail(broken_, "a piece walk reuses a crossing");
      used[xa] = used[xb] = true;
      piece.push_back({xs[xa].v, line_id});  // material segment xa -> xb on the cut line

      // The ring arc leaving xb on this side, up to the next crossing on that ring.
      const Crossing& cb = xs[xb];
      const std::uint32_t len = ring_len(cb.ring);
      const std::uint32_t e = cb.edge, en = (e + 1) % len;
      const bool forward = sign * K::orient(q, q2, pt(vid(cb.ring, en))) > 0;
      piece.push_back({cb.v, lid(cb.ring, e)});
      std::size_t xc;
      if (forward) {
        xc = ring_next[xb];
        const std::uint32_t f = xs[xc].edge;
        for (std::uint32_t k = en;; k = (k + 1) % len) {
          piece.push_back({vid(cb.ring, k), lid(cb.ring, k)});
          if (k == f) break;
        }
      } else {
        xc = ring_prev[xb];
        const std::uint32_t stop = (xs[xc].edge + 1) % len;
        for (std::uint32_t k = e;; k = (k + len - 1) % len) {
          piece.push_back({vid(cb.ring, k), lid(cb.ring, (k + len - 1) % len)});
          if (k == stop) break;
        }
      }
      p = pos[xc];
      if (p % 2 != 0) return fail(broken_, "a ring arc ends at the far end of a material segment");
      if (xc == seq[start]) break;  // back at the start: the piece is closed
    }
    pieces.push_back(std::move(piece));
  }
  return true;
}

template <class K>
bool HoleEliminator<K>::cut_node(Work& w, std::vector<Work>& stack)
{
  const Cut cut = choose_cut(w);
  const Point<FT> q2{cut.q.x + cut.d.x, cut.q.y + cut.d.y};
  const auto line_id = static_cast<std::uint32_t>(lines_.size());
  lines_.push_back({cut.q, cut.d});

  // ---- crossings, ring by ring in edge order
  std::vector<Crossing> xs;
  std::vector<std::pair<std::size_t, std::size_t>> groups;  // [begin, end) of each crossed ring
  auto scan = [&](int ring, std::uint32_t len, auto vid, auto lid) {
    const std::size_t begin = xs.size();
    for (std::uint32_t e = 0; e < len; ++e) {
      const int oa = K::orient(cut.q, q2, pt(vid(e)));
      const int ob = K::orient(cut.q, q2, pt(vid((e + 1) % len)));
      if (oa == 0 || ob == 0) return fail(broken_, "the cut line passes through a vertex");
      if (oa == ob) continue;
      Point<FT> x = intersect_lines(lines_[line_id], lines_[lid(e)]);
      if constexpr (!is_exact_kernel<K>::value) {
        const Point<FT>& ea = pt(vid(e));
        const Point<FT>& eb = pt(vid((e + 1) % len));
        if (!point_finite(x) || !near_box(ea, eb, x, 1e-9 * scale_))
          x = intersect_through(cut.q, q2, ea, eb);
      }
      if (!point_finite(x)) return fail(broken_, "a constructed point is not finite");
      const FT t = (x.x - cut.q.x) * cut.d.x + (x.y - cut.q.y) * cut.d.y;
      xs.push_back({t, static_cast<std::uint32_t>(out_.pool.size()), ring, e});
      out_.pool.push_back(x);
    }
    if (xs.size() > begin) groups.emplace_back(begin, xs.size());
    return true;
  };
  if (!scan(-1, static_cast<std::uint32_t>(w.outer.size()),
            [&](std::uint32_t k) { return w.outer[k].v; },
            [&](std::uint32_t k) { return w.outer[k].line; }))
    return false;
  std::vector<int> pierced;
  for (int h : w.holes) {
    const std::size_t before = xs.size();
    if (!scan(h, hole_len_[h], [&](std::uint32_t k) { return hole_off_[h] + k; },
              [&](std::uint32_t k) { return hole_off_[h] + k; }))
      return false;
    if (xs.size() > before) pierced.push_back(h);
  }

  const std::size_t m = xs.size();
  if (pierced.empty()) return fail(broken_, "a cut pierces no hole");
  if (m % 2 != 0) return fail(broken_, "odd number of crossings");
  for (const auto& g : groups)
    if ((g.second - g.first) % 2 != 0) return fail(broken_, "a ring is crossed an odd number of times");

  std::vector<std::size_t> order(m);
  for (std::size_t i = 0; i < m; ++i) order[i] = i;
  std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return xs[a].t < xs[b].t; });
  for (std::size_t i = 0; i + 1 < m; ++i)
    if (!(xs[order[i]].t < xs[order[i + 1]].t)) return fail(broken_, "two crossings coincide");
  if (xs[order.front()].ring != -1 || xs[order.back()].ring != -1)
    return fail(broken_, "the outermost crossings are not on the outer ring");

  std::vector<std::size_t> ring_next(m), ring_prev(m);
  for (const auto& g : groups)
    for (std::size_t i = g.first; i < g.second; ++i) {
      ring_next[i] = i + 1 < g.second ? i + 1 : g.first;
      ring_prev[i] = i > g.first ? i - 1 : g.second - 1;
    }

  // ---- pieces: left of the line, then right
  std::vector<std::vector<Corner1>> pieces;
  if (!trace_side(+1, line_id, cut.q, q2, xs, order, ring_next, ring_prev, w, pieces)) return false;
  const std::size_t n_left = pieces.size();
  if (!trace_side(-1, line_id, cut.q, q2, xs, order, ring_next, ring_prev, w, pieces)) return false;

  auto piece_ring = [&](std::size_t c) {
    Ring<FT> r;
    r.reserve(pieces[c].size());
    for (const Corner1& k : pieces[c]) r.push_back(pt(k.v));
    return r;
  };
  if (opt_.check)
    for (std::size_t c = 0; c < pieces.size(); ++c)
      if (ring_orientation<K>(piece_ring(c)) <= 0)
        return fail(broken_, "a piece is not counter-clockwise");

  // ---- un-pierced holes go to the piece containing them
  std::vector<std::vector<int>> child_holes(pieces.size());
  for (int h : w.holes) {
    if (std::find(pierced.begin(), pierced.end(), h) != pierced.end()) continue;
    const Point<FT>& u = pt(hole_off_[h]);
    const bool left = K::orient(cut.q, q2, u) > 0;
    const std::size_t lo = left ? 0 : n_left, hi = left ? n_left : pieces.size();
    int owner = -1;
    if (hi - lo == 1 && !opt_.check) {
      owner = static_cast<int>(lo);
    } else {
      for (std::size_t c = lo; c < hi; ++c) {
        if (locate<K>(u, piece_ring(c)) != Location::Inside) continue;
        if (owner >= 0) return fail(broken_, "a hole lies in two pieces");
        owner = static_cast<int>(c);
        if (!opt_.check) break;
      }
    }
    if (owner < 0) return fail(broken_, "a hole lies in no piece");
    child_holes[static_cast<std::size_t>(owner)].push_back(h);
  }

  // ---- record, then continue depth-first with the pieces in order
  ++out_.stats.cuts;
  out_.stats.max_crossings = std::max(out_.stats.max_crossings, m);
  out_.stats.fallback_rotated += (cut.fallback == CutFallback::Rotated);
  out_.stats.fallback_midpoint += (cut.fallback == CutFallback::EdgeMidpoint);
  out_.stats.degenerate_directions += cut.degenerate_direction;
  std::vector<int> child_nodes(pieces.size(), -1);
  if (w.node >= 0) {
    for (std::size_t c = 0; c < pieces.size(); ++c) {
      child_nodes[c] = static_cast<int>(out_.tree.size());
      HoleNode<FT> child;
      child.parent = w.node;
      child.holes = child_holes[c];
      out_.tree.push_back(std::move(child));
    }
    HoleNode<FT>& nd = out_.tree[static_cast<std::size_t>(w.node)];
    nd.children = child_nodes;
    nd.has_cut = true;
    nd.cut_point = cut.q;
    nd.cut_dir = cut.d;
    nd.fallback = cut.fallback;
    nd.pierced = pierced;
  }
  for (std::size_t c = pieces.size(); c-- > 0;)
    stack.push_back({std::move(pieces[c]), std::move(child_holes[c]), child_nodes[c], w.depth + 1});
  return true;
}

template <class K>
void HoleEliminator<K>::run(PolygonWithHoles<FT> in)
{
  if (opt_.validate_input) {
    ValidationReport rep = validate<K>(in);
    if (!rep.ok()) {
      out_.diag = {Status::InvalidInput, "input is not in the supported class", std::move(rep.issues)};
      return;
    }
  } else if (in.outer.size() < 3) {
    out_.diag = {Status::InvalidInput, "outer ring has fewer than 3 vertices", {}};
    return;
  }
  make_ccw<K>(in.outer);
  for (const Point<FT>& v : in.outer)
    scale_ = std::max({scale_, std::abs(K::to_double(v.x)), std::abs(K::to_double(v.y))});
  tol_ = K::clearance * scale_;

  // Pool and lines: the input edge starting at vertex id v has line id v.
  auto add_ring = [&](const Ring<FT>& r) {
    const auto off = static_cast<std::uint32_t>(out_.pool.size());
    for (std::size_t k = 0; k < r.size(); ++k) {
      const Point<FT>& a = r[k];
      const Point<FT>& b = r[next_index(k, r.size())];
      out_.pool.push_back(a);
      lines_.push_back({a, {b.x - a.x, b.y - a.y}});
    }
    return off;
  };
  const std::uint32_t outer_off = add_ring(in.outer);
  for (const Ring<FT>& h : in.holes) {
    hole_off_.push_back(add_ring(h));
    hole_len_.push_back(static_cast<std::uint32_t>(h.size()));
  }
  out_.stats.holes = in.holes.size();

  Work root;
  for (std::uint32_t k = 0; k < in.outer.size(); ++k) root.outer.push_back({outer_off + k, outer_off + k});
  for (std::size_t h = 0; h < in.holes.size(); ++h) root.holes.push_back(static_cast<int>(h));
  root.node = opt_.keep_tree ? 0 : -1;
  root.depth = 0;
  if (opt_.keep_tree) {
    out_.tree.emplace_back();
    out_.tree[0].holes = root.holes;
  }

  std::vector<Work> stack;
  stack.push_back(std::move(root));
  while (!stack.empty()) {
    Work w = std::move(stack.back());
    stack.pop_back();
    out_.stats.max_depth = std::max(out_.stats.max_depth, w.depth);
    if (w.holes.empty()) {
      std::vector<std::uint32_t> ids;
      ids.reserve(w.outer.size());
      for (const Corner1& c : w.outer) ids.push_back(c.v);
      if (w.node >= 0) out_.tree[static_cast<std::size_t>(w.node)].leaf = static_cast<int>(out_.leaves.size());
      out_.leaves.push_back(std::move(ids));
      continue;
    }
    if (out_.stats.cuts >= out_.stats.holes) {
      fail(broken_, "more cuts than holes");
      return;
    }
    if (!cut_node(w, stack)) return;
  }
}

} // namespace detail

/// Eliminates the holes of a polygon: returns simple, hole-free pieces that
/// tile it. The input must be a valid complex polygon (see validate); ring
/// orientations are arbitrary.
template <class K>
HoleElimination<typename K::FT> eliminate_holes(PolygonWithHoles<typename K::FT> input,
                                                const HoleElimOptions& opt = {})
{
  HoleElimination<typename K::FT> out;
  detail::HoleEliminator<K>(out, opt).run(std::move(input));
  return out;
}

} // namespace pd
