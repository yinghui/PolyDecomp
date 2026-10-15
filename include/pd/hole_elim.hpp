// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#pragma once
/// \file
/// Stage 1: hole elimination by line cuts.
///
/// A node is an outer ring (counter-clockwise) with the input holes it contains.
/// Each cut is a line that crosses at least one of the node's holes:
///  - one hole: the line through the anchor a = (v_{m-1} + v_m + v_{m+1}) / 3,
///    v_m being the outer vertex nearest to the hole's vertex-average centroid
///    c, and c itself;
///  - two or more holes: the line through the centroids of the two holes whose
///    centroids are farthest apart.
/// The line may pass through vertices and run along edges ("contacts"). Only
/// if it has no direction (coincident centroids, a == c) is it replaced by a
/// line through a point q of a hole in a vertex-free direction (the fallback);
/// floating-point kernels also fall back when a vertex lies within their
/// clearance of the line without lying on it, and for a dangling contact.
///
/// Pieces on each side S of the line are traced as if the line were shifted
/// infinitesimally into S. An edge then contributes a crossing to the list of
/// side S if the line crosses its interior, or if one endpoint lies on the line
/// and the other strictly on S; in the second case the crossing point is that
/// endpoint. So a vertex whose edges lie on both sides is listed once on each
/// side; a vertex whose edges both lie on S is listed twice on S and not on the
/// other side; and a run of edges on the line contributes its end vertices to
/// the sides its two outer edges lie on. Two copies of one vertex are ordered
/// by the direction of their edges. Each list is walked as before: from a
/// crossing x_{2i-1} along the line to x_{2i}, then along the ring towards S
/// up to the next crossing on that ring; each closed walk is one piece, and
/// repeated consecutive vertices (zero-length segments) are removed.
///
/// A hole that meets the line in a single vertex v, lying otherwise on one side
/// (a dangling contact), would make a piece that touches itself at v. The exact
/// path splits such a hole with the point c = v + ε d (ε an infinitesimal, see
/// kernel/symbolic.hpp), turning the contact into an edge on the line; the
/// result is evaluated at ε = 0, where that piece is weakly simple (it touches
/// itself at v). Convex decomposition of the perturbed pieces then yields
/// convex pieces that tile the input (see decompose).

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "pd/detail/lines.hpp"
#include "pd/geometry.hpp"
#include "pd/status.hpp"
#include "pd/validate.hpp"
#ifdef PD_HAVE_GMP
#include <gmpxx.h>

#include "pd/kernel/exact.hpp"
#include "pd/kernel/symbolic.hpp"
#endif

namespace pd {

enum class CutFallback {
  None,          ///< the preferred line was used
  Rotated,       ///< a line through the hole centroid in a vertex-free direction
  EdgeMidpoint,  ///< the centroid is a vertex: a line through an edge midpoint
};

/// How the preferred cut line of a node is chosen. Every rule gives a line
/// through the vertex average of one of the node's holes, so every cut crosses a
/// hole and hole elimination terminates after at most h cuts. Only FarthestPair
/// (the default) is covered by the bounded-coordinates guarantee: the other rules
/// exist to compare the choice of the cut (the paper's ablation).
enum class CutRule {
  FarthestPair,       ///< one hole: anchor rule; more: averages of the two holes farthest apart
  ClosestPair,        ///< one hole: anchor rule; more: averages of the two holes closest together
  NearestToBoundary,  ///< anchor rule for the hole whose average is nearest to an outer vertex
  Horizontal,         ///< horizontal line through the average of the node's first hole
};

struct HoleElimOptions {
  bool validate_input = true;  ///< check the input class first
  bool check = false;          ///< extra invariant checks (piece orientation, unique hole owner)
  bool keep_tree = true;       ///< record the decomposition tree
  CutRule cut_rule = CutRule::FarthestPair;  ///< choice of the preferred cut line
};

struct HoleElimStats {
  std::size_t holes = 0;
  std::size_t cuts = 0;
  std::size_t max_depth = 0;
  std::size_t max_crossings = 0;          ///< most crossings on one side of a single cut
  std::size_t vertex_contacts = 0;        ///< vertices on a cut line (not on an edge run)
  std::size_t edge_contacts = 0;          ///< runs of edges lying on a cut line
  std::size_t dangles = 0;                ///< holes split symbolically (dangling contacts)
  std::size_t fallback_rotated = 0;       ///< cuts that used CutFallback::Rotated
  std::size_t fallback_midpoint = 0;      ///< cuts that used CutFallback::EdgeMidpoint
  std::size_t degenerate_directions = 0;  ///< preferred line had no direction
  bool symbolic = false;                  ///< computed with the symbolic kernel (then ε = 0)
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
  std::vector<Point<FT>> pool;                     ///< input vertices (outer, then holes), then constructed points
  std::vector<std::vector<std::uint32_t>> leaves;  ///< hole-free pieces, CCW, as pool ids
  std::vector<HoleNode<FT>> tree;                  ///< tree[0] is the root (if kept)
  HoleElimStats stats;
  Diagnostics diag;
  /// Internal: an exact run met a dangling contact and must be redone with the
  /// symbolic kernel (eliminate_holes and decompose do this automatically).
  bool needs_symbolic = false;

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

/// Removes consecutive repeated vertices (zero-length edges) from a closed
/// walk, keeping the entry whose outgoing edge has positive length.
inline void drop_repeats(std::vector<Corner1>& piece)
{
  std::vector<Corner1> out;
  out.reserve(piece.size());
  for (std::size_t i = 0; i < piece.size(); ++i) {
    if (i + 1 < piece.size() && piece[i + 1].v == piece[i].v) continue;
    out.push_back(piece[i]);
  }
  while (out.size() > 1 && out.back().v == out.front().v) out.pop_back();
  piece = std::move(out);
}

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
    FT t;                 // position along the cut direction
    std::uint32_t v;      // pool id of the crossing point (the vertex itself for a contact)
    int ring;             // -1: the node's outer ring; otherwise an input hole index
    std::uint32_t edge;   // edge of that ring containing the crossing
    bool contact = false; // the crossing point is an endpoint of the edge, on the line
    std::uint32_t far = 0;  // contact: pool id of the edge's other endpoint (strictly on the side)
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
  std::vector<std::vector<std::uint32_t>> hole_line_;  ///< line id of each hole edge
  double scale_ = 1;  ///< largest absolute input coordinate (at least 1)
  double tol_ = 0;    ///< K::clearance * scale_ (0: exact incidence only)

  static constexpr bool floating = !is_exact_kernel<K>::value;

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

  // Ring access: ring -1 is the node's outer ring, otherwise an input hole index.
  std::uint32_t ring_len(const Work& w, int ring) const
  {
    return ring < 0 ? static_cast<std::uint32_t>(w.outer.size()) : hole_len_[ring];
  }
  std::uint32_t vid(const Work& w, int ring, std::uint32_t k) const
  {
    return ring < 0 ? w.outer[k].v : hole_off_[ring] + k;
  }
  std::uint32_t lid(const Work& w, int ring, std::uint32_t k) const
  {
    return ring < 0 ? w.outer[k].line : hole_line_[ring][k];
  }

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
  /// floating-point kernel, within tol_ of it). Used for fallback lines.
  bool line_avoids_vertices(const Work& w, const Point<FT>& q, const Point<FT>& q2) const
  {
    bool ok = true;
    for_each_vertex(w, [&](std::uint32_t id) {
      ok = ok && K::orient(q, q2, pt(id)) != 0 && (tol_ == 0 || distance_to_line(q, q2, pt(id)) > tol_);
    });
    return ok;
  }

  /// Index (in ring order) of the single vertex of hole h on the line, if the
  /// hole meets the line only there and lies otherwise on one side (a dangling
  /// contact); -1 otherwise. `sign` receives that side.
  long dangling_vertex(int h, const Point<FT>& q, const Point<FT>& q2, int& sign) const
  {
    long on = -1;
    int side = 0;
    for (std::uint32_t k = 0; k < hole_len_[h]; ++k) {
      const int o = K::orient(q, q2, pt(hole_off_[h] + k));
      if (o == 0) {
        if (on >= 0) return -1;  // two or more vertices on the line
        on = static_cast<long>(k);
      } else if (side == 0) {
        side = o;
      } else if (o != side) {
        return -1;  // the line crosses the hole
      }
    }
    sign = side;
    return on;
  }

  /// Whether the preferred line can be used as it is. Exact kernels handle every
  /// contact; floating-point kernels need every vertex either exactly on the line
  /// or beyond their clearance, and no dangling contact (no symbolic split).
  bool line_usable(const Work& w, const Point<FT>& q, const Point<FT>& q2) const
  {
    if constexpr (!floating) {
      return true;
    } else {
      bool ok = true;
      for_each_vertex(w, [&](std::uint32_t id) {
        ok = ok && (tol_ == 0 || K::orient(q, q2, pt(id)) == 0 || distance_to_line(q, q2, pt(id)) > tol_);
      });
      int side;
      for (int h : w.holes) ok = ok && dangling_vertex(h, q, q2, side) < 0;
      return ok;
    }
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

  /// For two contact crossings at the same vertex v on the side `sign`, with far
  /// endpoints a and b: true if the copy on the edge to a comes first along d
  /// when the line is shifted infinitesimally into that side. The shifted line
  /// meets the edge to w at the offset ((w - v) . d) / |d x (w - v)| along d.
  bool shifted_before(int sign, const Point<FT>& v, const Point<FT>& a, const Point<FT>& b,
                      const Point<FT>& d) const
  {
    auto decide = [&](auto X, auto conv) {
      using T = decltype(X);
      const T dx = conv(d.x), dy = conv(d.y), vx = conv(v.x), vy = conv(v.y);
      const T ax = conv(a.x) - vx, ay = conv(a.y) - vy, bx = conv(b.x) - vx, by = conv(b.y) - vy;
      const T da = ax * dx + ay * dy, db = bx * dx + by * dy;  // along d
      const T ca = dx * ay - dy * ax, cb = dx * by - dy * bx;  // across (sign = side)
      const T det = da * cb - db * ca;
      return (det > 0) - (det < 0);
    };
    int s;
    if constexpr (!floating) {
      s = decide(FT(0), [](const FT& x) { return x; });
    } else {
#ifdef PD_HAVE_GMP
      s = decide(mpq_class(0), [](const FT& x) { return mpq_class(double(x)); });
#else
      s = decide(0.0L, [](const FT& x) { return static_cast<long double>(x); });
#endif
    }
    return sign * s < 0;
  }

  /// Splits the dangling hole h at its vertex k (on the line, the rest of the
  /// hole on side `sign`) with c = v + ε d, inserted next to the neighbour whose
  /// direction is closer to d, so that the quad stays simple. Symbolic kernel only.
  void split_dangle(int h, std::uint32_t k, int sign, const Point<FT>& d, std::uint32_t cut_line)
  {
    if constexpr (is_symbolic_kernel<K>::value) {
      const std::uint32_t len = hole_len_[h], off = hole_off_[h];
      const std::uint32_t kp = (k + len - 1) % len, kn = (k + 1) % len;
      const Point<FT> v = pt(off + k);
      const Point<FT> c{v.x + FT::eps() * d.x, v.y + FT::eps() * d.y};
      // The neighbour closer in angle to d borders the sliver triangle (v, c, w).
      const bool after = sign * K::orient(v, pt(off + kn), pt(off + kp)) > 0;
      std::vector<Point<FT>> ring;
      std::vector<std::uint32_t> lines;
      for (std::uint32_t j = 0; j < len; ++j) {
        const std::uint32_t old_line = hole_line_[h][j];
        if (j == k && after) {  // ..., v, c, kn, ...: v -> c on the cut line, c -> kn new
          ring.push_back(v);
          lines.push_back(cut_line);
          ring.push_back(c);
          const Point<FT>& wn = pt(off + kn);
          lines.push_back(static_cast<std::uint32_t>(lines_.size()));
          lines_.push_back({c, {wn.x - c.x, wn.y - c.y}});
        } else if (j == kp && !after) {  // ..., kp, c, v, ...: kp -> c new, c -> v on the cut line
          const Point<FT>& wp = pt(off + kp);
          ring.push_back(wp);
          lines.push_back(static_cast<std::uint32_t>(lines_.size()));
          lines_.push_back({wp, {c.x - wp.x, c.y - wp.y}});
          ring.push_back(c);
          lines.push_back(cut_line);
        } else {
          ring.push_back(pt(off + j));
          lines.push_back(old_line);
        }
      }
      hole_off_[h] = static_cast<std::uint32_t>(out_.pool.size());
      hole_len_[h] = static_cast<std::uint32_t>(ring.size());
      for (Point<FT>& p : ring) out_.pool.push_back(std::move(p));
      hole_line_[h] = std::move(lines);
      ++out_.stats.dangles;
    } else {
      (void)h, (void)k, (void)sign, (void)d, (void)cut_line;
    }
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
  const std::size_t m = w.outer.size();
  // Outer vertex nearest to c (the last one in ring order on ties), and its squared distance.
  const auto nearest_outer = [&](const Point<FT>& c, FT& best_d2) {
    std::size_t best = 0;
    for (std::size_t i = 0; i < m; ++i) {
      const Point<FT>& v = pt(w.outer[i].v);
      const FT dx = v.x - c.x, dy = v.y - c.y;
      const FT d2 = dx * dx + dy * dy;
      if (i == 0 || d2 <= best_d2) {
        best = i;
        best_d2 = d2;
      }
    }
    return best;
  };
  if (opt_.cut_rule == CutRule::Horizontal) {
    hole = w.holes[0];
    q0 = centroid(hole);
    d0 = {FT(1), FT(0)};
  } else if (w.holes.size() == 1 || opt_.cut_rule == CutRule::NearestToBoundary) {
    // Anchor rule, for the node's only hole or (NearestToBoundary) for the hole whose
    // average is nearest to an outer vertex (the first one on ties).
    hole = w.holes[0];
    FT best_d2 = 0;
    std::size_t best = nearest_outer(centroid(hole), best_d2);
    for (std::size_t k = 1; k < w.holes.size(); ++k) {
      FT d2 = 0;
      const std::size_t i = nearest_outer(centroid(w.holes[k]), d2);
      if (d2 < best_d2) {
        hole = w.holes[k];
        best = i;
        best_d2 = d2;
      }
    }
    const Point<FT> c = centroid(hole);
    const Point<FT>& a0 = pt(w.outer[prev_index(best, m)].v);
    const Point<FT>& a1 = pt(w.outer[best].v);
    const Point<FT>& a2 = pt(w.outer[next_index(best, m)].v);
    q0 = {(a0.x + a1.x + a2.x) / 3, (a0.y + a1.y + a2.y) / 3};
    d0 = {c.x - q0.x, c.y - q0.y};
  } else {
    // Farthest (or, for ClosestPair, closest) pair of centroids, the last one in pair
    // order on ties.
    const bool closest = opt_.cut_rule == CutRule::ClosestPair;
    std::vector<Point<FT>> cs;
    for (int h : w.holes) cs.push_back(centroid(h));
    std::size_t r = 0, s = 1;
    FT best_d2 = -1;
    for (std::size_t i = 0; i + 1 < cs.size(); ++i)
      for (std::size_t j = i + 1; j < cs.size(); ++j) {
        const FT dx = cs[j].x - cs[i].x, dy = cs[j].y - cs[i].y;
        const FT d2 = dx * dx + dy * dy;
        if (closest ? (best_d2 < 0 || d2 <= best_d2) : d2 >= best_d2) {
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
  if (!cut.degenerate_direction && line_usable(w, q0, Point<FT>{q0.x + d0.x, q0.y + d0.y})) {
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

/// Traces the pieces on one side (sign +1: left of the line, -1: right) from
/// that side's crossing list (`order`: ascending along the cut direction).
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

  // Side of a vertex relative to the line shifted into this side: a vertex on
  // the line lies on the far side.
  auto on_this_side = [&](std::uint32_t id) {
    const int o = K::orient(q, q2, pt(id));
    return (o != 0 ? o : -sign) == sign;
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
      piece.push_back({xs[xa].v, line_id});  // segment xa -> xb on the cut line

      // The ring arc leaving xb on this side, up to the next crossing on that ring.
      const Crossing& cb = xs[xb];
      const std::uint32_t len = ring_len(w, cb.ring);
      const std::uint32_t e = cb.edge, en = (e + 1) % len;
      const bool forward = on_this_side(vid(w, cb.ring, en));
      piece.push_back({cb.v, lid(w, cb.ring, e)});
      std::size_t xc;
      if (forward) {
        xc = ring_next[xb];
        const std::uint32_t f = xs[xc].edge;
        for (std::uint32_t k = en;; k = (k + 1) % len) {
          piece.push_back({vid(w, cb.ring, k), lid(w, cb.ring, k)});
          if (k == f) break;
        }
      } else {
        xc = ring_prev[xb];
        const std::uint32_t stop = (xs[xc].edge + 1) % len;
        for (std::uint32_t k = e;; k = (k + len - 1) % len) {
          piece.push_back({vid(w, cb.ring, k), lid(w, cb.ring, (k + len - 1) % len)});
          if (k == stop) break;
        }
      }
      p = pos[xc];
      if (p % 2 != 0) return fail(broken_, "a ring arc ends at the far end of a segment on the line");
      if (xc == seq[start]) break;  // back at the start: the piece is closed
    }
    drop_repeats(piece);
    if (piece.size() < 3) return fail(broken_, "a piece has fewer than 3 vertices");
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

  // ---- dangling contacts: split symbolically (exact kernels), never in floating ones
  for (int h : w.holes) {
    int side = 0;
    const long k = dangling_vertex(h, cut.q, q2, side);
    if (k < 0) continue;
    if constexpr (floating) {
      return fail(broken_, "a cut line touches a hole at a single vertex");
    } else if constexpr (!is_symbolic_kernel<K>::value) {
      out_.needs_symbolic = true;
      return fail(Status::InternalError, "a dangling contact needs the symbolic kernel");
    } else {
      split_dangle(h, static_cast<std::uint32_t>(k), side, cut.d, line_id);
    }
  }

  // ---- side of every vertex, contacts
  std::vector<int> rings{-1};
  for (int h : w.holes) rings.push_back(h);
  std::vector<std::vector<int>> side_of(rings.size());
  for (std::size_t r = 0; r < rings.size(); ++r) {
    const std::uint32_t len = ring_len(w, rings[r]);
    side_of[r].resize(len);
    for (std::uint32_t k = 0; k < len; ++k) side_of[r][k] = K::orient(cut.q, q2, pt(vid(w, rings[r], k)));
    // Runs of consecutive vertices on the line: single vertices and edge runs.
    for (std::uint32_t k = 0; k < len; ++k) {
      if (side_of[r][k] != 0 || side_of[r][(k + len - 1) % len] == 0) continue;  // not the start of a run
      std::uint32_t run = 0;
      while (run < len && side_of[r][(k + run) % len] == 0) ++run;
      if (run == 1) ++out_.stats.vertex_contacts;
      else ++out_.stats.edge_contacts;
    }
  }
  auto t_of = [&](const Point<FT>& x) -> FT { return (x.x - cut.q.x) * cut.d.x + (x.y - cut.q.y) * cut.d.y; };

  // ---- proper crossings (shared by both sides), constructed once
  std::vector<std::vector<long>> proper(rings.size());
  std::vector<std::vector<FT>> proper_t(rings.size());
  for (std::size_t r = 0; r < rings.size(); ++r) {
    const int ring = rings[r];
    const std::uint32_t len = ring_len(w, ring);
    proper[r].assign(len, -1);
    proper_t[r].resize(len);
    for (std::uint32_t e = 0; e < len; ++e) {
      const int oa = side_of[r][e], ob = side_of[r][(e + 1) % len];
      if (oa * ob >= 0) continue;
      Point<FT> x = intersect_lines(lines_[line_id], lines_[lid(w, ring, e)]);
      if constexpr (floating) {
        const Point<FT>& ea = pt(vid(w, ring, e));
        const Point<FT>& eb = pt(vid(w, ring, (e + 1) % len));
        if (!point_finite(x) || !near_box(ea, eb, x, 1e-9 * scale_)) x = intersect_through(cut.q, q2, ea, eb);
      }
      if (!point_finite(x)) return fail(broken_, "a constructed point is not finite");
      proper[r][e] = static_cast<long>(out_.pool.size());
      proper_t[r][e] = t_of(x);
      out_.pool.push_back(x);
    }
  }

  // ---- each side: crossing list of the line shifted into that side, then the walk
  std::vector<std::vector<Corner1>> pieces;
  std::size_t n_left = 0;
  std::vector<int> pierced;
  for (int sign : {+1, -1}) {
    std::vector<Crossing> xs;
    std::vector<std::pair<std::size_t, std::size_t>> groups;  // [begin, end) per crossed ring
    for (std::size_t r = 0; r < rings.size(); ++r) {
      const int ring = rings[r];
      const std::uint32_t len = ring_len(w, ring);
      const std::size_t begin = xs.size();
      for (std::uint32_t e = 0; e < len; ++e) {
        const std::uint32_t en = (e + 1) % len;
        const int oa = side_of[r][e], ob = side_of[r][en];
        if (proper[r][e] >= 0) {
          xs.push_back({proper_t[r][e], static_cast<std::uint32_t>(proper[r][e]), ring, e});
        } else if (oa == 0 && ob == sign) {
          const std::uint32_t a = vid(w, ring, e);
          xs.push_back({t_of(pt(a)), a, ring, e, true, vid(w, ring, en)});
        } else if (ob == 0 && oa == sign) {
          const std::uint32_t b = vid(w, ring, en);
          xs.push_back({t_of(pt(b)), b, ring, e, true, vid(w, ring, e)});
        }
      }
      if (xs.size() > begin) {
        groups.emplace_back(begin, xs.size());
        if (ring >= 0 && std::find(pierced.begin(), pierced.end(), ring) == pierced.end()) pierced.push_back(ring);
      }
    }
    const std::size_t m = xs.size();
    out_.stats.max_crossings = std::max(out_.stats.max_crossings, m);
    if (m % 2 != 0) return fail(broken_, "odd number of crossings");
    for (const auto& g : groups)
      if ((g.second - g.first) % 2 != 0) return fail(broken_, "a ring is crossed an odd number of times");

    std::vector<std::size_t> order(m);
    for (std::size_t i = 0; i < m; ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
      if (xs[a].t != xs[b].t) return xs[a].t < xs[b].t;
      if (xs[a].v != xs[b].v || !xs[a].contact || !xs[b].contact) return a < b;  // rejected below
      return shifted_before(sign, pt(xs[a].v), pt(xs[a].far), pt(xs[b].far), cut.d);
    });
    for (std::size_t i = 0; i + 1 < m; ++i) {
      const Crossing& x0 = xs[order[i]];
      const Crossing& x1 = xs[order[i + 1]];
      if (x0.t < x1.t) continue;
      if (!(x0.v == x1.v && x0.contact && x1.contact)) return fail(broken_, "two crossings coincide");
    }
    if (m > 0 && (xs[order.front()].ring != -1 || xs[order.back()].ring != -1))
      return fail(broken_, "the outermost crossings are not on the outer ring");

    std::vector<std::size_t> ring_next(m), ring_prev(m);
    for (const auto& g : groups)
      for (std::size_t i = g.first; i < g.second; ++i) {
        ring_next[i] = i + 1 < g.second ? i + 1 : g.first;
        ring_prev[i] = i > g.first ? i - 1 : g.second - 1;
      }
    if (!trace_side(sign, line_id, cut.q, q2, xs, order, ring_next, ring_prev, w, pieces)) return false;
    if (sign > 0) n_left = pieces.size();
  }
  if (pierced.empty()) return fail(broken_, "a cut pierces no hole");

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

  // ---- un-pierced holes (strictly on one side) go to the piece containing them
  std::vector<std::vector<int>> child_holes(pieces.size());
  for (int h : w.holes) {
    if (std::find(pierced.begin(), pierced.end(), h) != pierced.end()) continue;
    const Point<FT>& u = pt(hole_off_[h]);
    const int o = K::orient(cut.q, q2, u);
    if (o == 0) return fail(broken_, "an un-pierced hole touches the cut line");
    const bool left = o > 0;
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
    const std::uint32_t off = add_ring(h);
    hole_off_.push_back(off);
    hole_len_.push_back(static_cast<std::uint32_t>(h.size()));
    std::vector<std::uint32_t> ids(h.size());
    for (std::uint32_t k = 0; k < h.size(); ++k) ids[k] = off + k;
    hole_line_.push_back(std::move(ids));
  }
  out_.stats.holes = in.holes.size();
  out_.stats.symbolic = is_symbolic_kernel<K>::value;

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

/// Runs the eliminator, without the symbolic rerun.
template <class K>
HoleElimination<typename K::FT> eliminate_holes_once(PolygonWithHoles<typename K::FT> input,
                                                     const HoleElimOptions& opt)
{
  HoleElimination<typename K::FT> out;
  HoleEliminator<K>(out, opt).run(std::move(input));
  return out;
}

#ifdef PD_HAVE_GMP
inline PolygonWithHoles<EpsQ> to_symbolic(const PolygonWithHoles<mpq_class>& p)
{
  auto ring = [](const Ring<mpq_class>& r) {
    Ring<EpsQ> out;
    out.reserve(r.size());
    for (const auto& v : r) out.push_back({EpsQ(v.x), EpsQ(v.y)});
    return out;
  };
  PolygonWithHoles<EpsQ> out{ring(p.outer), {}};
  for (const auto& h : p.holes) out.holes.push_back(ring(h));
  return out;
}

/// Removes consecutive pool ids whose points coincide. At ε = 0 the edge from a
/// dangling vertex v1 to its split vertex v = v1 + ε d has length zero; a piece
/// that runs along it would otherwise repeat the vertex.
inline void drop_coincident(std::vector<std::uint32_t>& g, const std::vector<Point<mpq_class>>& pool)
{
  std::vector<std::uint32_t> out;
  out.reserve(g.size());
  for (std::uint32_t id : g)
    if (out.empty() || !(pool[out.back()] == pool[id])) out.push_back(id);
  while (out.size() > 1 && pool[out.back()] == pool[out.front()]) out.pop_back();
  g = std::move(out);
}

/// The result of a symbolic run, evaluated at ε = 0.
inline HoleElimination<mpq_class> symbolic_limit(HoleElimination<EpsQ>&& s)
{
  HoleElimination<mpq_class> r;
  r.pool.reserve(s.pool.size());
  for (const auto& p : s.pool) r.pool.push_back(limit_point(p));
  r.leaves = std::move(s.leaves);
  for (auto& g : r.leaves) drop_coincident(g, r.pool);
  for (const auto& n : s.tree) {
    HoleNode<mpq_class> m;
    m.parent = n.parent;
    m.children = n.children;
    m.leaf = n.leaf;
    m.holes = n.holes;
    m.has_cut = n.has_cut;
    m.cut_point = limit_point(n.cut_point);
    m.cut_dir = limit_point(n.cut_dir);
    m.fallback = n.fallback;
    m.pierced = n.pierced;
    r.tree.push_back(std::move(m));
  }
  r.stats = s.stats;
  r.diag = std::move(s.diag);
  return r;
}
#endif

} // namespace detail

/// Eliminates the holes of a polygon: returns hole-free pieces that tile it.
/// The input must be a valid complex polygon (see validate); ring orientations
/// are arbitrary. The pieces are simple and counter-clockwise, except after a
/// dangling contact (stats.dangles > 0, exact kernel), where a piece may touch
/// itself at that vertex (it is weakly simple); decompose's convex pieces are
/// simple in every case.
template <class K>
HoleElimination<typename K::FT> eliminate_holes(PolygonWithHoles<typename K::FT> input,
                                                const HoleElimOptions& opt = {})
{
#ifdef PD_HAVE_GMP
  if constexpr (std::is_same_v<K, ExactKernel>) {
    HoleElimination<mpq_class> out = detail::eliminate_holes_once<K>(input, opt);
    if (!out.needs_symbolic) return out;
    HoleElimOptions s = opt;
    s.validate_input = false;  // already validated
    return detail::symbolic_limit(detail::eliminate_holes_once<SymbolicKernel>(detail::to_symbolic(input), s));
  }
#endif
  return detail::eliminate_holes_once<K>(std::move(input), opt);
}

} // namespace pd
