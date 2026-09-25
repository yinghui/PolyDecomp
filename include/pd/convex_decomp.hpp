// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#pragma once
/// \file
/// Stage 2: convex decomposition of a simple polygon by chord splitting.
///
/// The polygon is taken counter-clockwise. While a piece has a reflex vertex,
/// the first one, v_i, is removed by a chord: the incoming edge (v_{i-1}, v_i)
/// is extended beyond v_i to its first contact with the piece's boundary.
///  - Proper hit, inside edge (v_k, v_{k+1}) at a new point p:
///      P1 = (v_{k+1}, ..., v_{i-1}, p),   P2 = (v_i, ..., v_k, p).
///  - Vertex hit, exactly at a vertex w:
///      P1 = (w_{+1}, ..., v_{i-1}, w),    P2 = (v_i, ..., w).
/// v_i is convex in P2 and straight (hence dropped) in P1. Reflex flags are
/// inherited by all other vertices; only a hit vertex w has its two new angles
/// recomputed. Every split lowers the total reflex count by at least one, so a
/// polygon with n_c reflex vertices yields at most n_c + 1 convex pieces
/// (exactly n_c + 1 when no chord ends at a vertex that becomes convex on both
/// sides).
///
/// Every edge carries the id of the input edge whose supporting line contains
/// it; chords lie on such lines too. A new point p is computed as the
/// intersection of two input lines, so exact coordinates never grow with the
/// depth of the recursion, and in double precision each point is rounded once.

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

struct ConvexOptions {
  bool validate_input = true;  ///< check that the input is a simple polygon first
  bool check = false;          ///< after each split, recompute all reflex flags and
                               ///< compare them with the inherited ones
  bool keep_tree = true;       ///< record the decomposition tree
};

struct ConvexStats {
  std::size_t input_vertices = 0;
  std::size_t reflex_vertices = 0;     ///< n_c of the input
  std::size_t splits = 0;
  std::size_t vertex_hits = 0;         ///< chords ending exactly at a vertex
  std::size_t double_resolutions = 0;  ///< vertex hits that also resolved a reflex hit vertex
  std::size_t max_depth = 0;
};

/// A node of the decomposition tree. Internal nodes record their split; leaves
/// point into ConvexDecomposition::leaves.
struct ConvexNode {
  int parent = -1;
  int child[2] = {-1, -1};  ///< [0]: P1 (without the reflex vertex), [1]: P2
  int leaf = -1;
  std::uint32_t reflex_pos = 0;  ///< position of the split reflex vertex in this node's ring
  bool vertex_hit = false;
  std::uint32_t chord_from = 0;  ///< pool id of the reflex vertex
  std::uint32_t chord_to = 0;    ///< pool id of the hit point or hit vertex
};

template <class FT>
struct ConvexDecomposition {
  std::vector<Point<FT>> pool;                  ///< input vertices, then constructed points
  std::vector<std::vector<std::uint32_t>> leaves;  ///< convex pieces, CCW, as pool ids
  std::vector<ConvexNode> tree;                 ///< tree[0] is the root (if kept)
  ConvexStats stats;
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

/// Rebuilds the ring (as pool ids, in its original vertex order) of tree node
/// `node` by gluing its descendants' leaves back together along their chords.
inline std::vector<std::uint32_t> reconstruct_ring(const std::vector<ConvexNode>& tree,
                                                   const std::vector<std::vector<std::uint32_t>>& leaves,
                                                   int node)
{
  const ConvexNode& nd = tree[static_cast<std::size_t>(node)];
  if (nd.leaf >= 0) return leaves[static_cast<std::size_t>(nd.leaf)];
  const std::vector<std::uint32_t> p1 = reconstruct_ring(tree, leaves, nd.child[0]);
  const std::vector<std::uint32_t> p2 = reconstruct_ring(tree, leaves, nd.child[1]);
  // Both children end with the chord's far end; P2 starts at the reflex vertex.
  std::vector<std::uint32_t> ring(p2.begin(), nd.vertex_hit ? p2.end() : p2.end() - 1);
  ring.insert(ring.end(), p1.begin(), p1.end() - 1);
  // The joined ring starts at the reflex vertex, which sat at reflex_pos.
  std::rotate(ring.begin(), ring.end() - nd.reflex_pos, ring.end());
  return ring;
}

namespace detail {

/// A candidate first contact of the ray o -> v_i: the vertex at position pos,
/// or a proper crossing of edge (pos, pos + 1) whose start lies on side sigma
/// (= orient(o, v_i, start)) of the ray's line. Along the line the orientation
/// of a point p(s) = o + s (v_i - o) relative to such an edge is
/// sign(sigma * (s - s_X)), s_X being the crossing's parameter, which turns
/// every comparison into exact orientation tests, except crossing vs crossing.
struct Hit {
  bool vertex;
  std::size_t pos;
  int sigma;
};

/// True if the crossing of the ray o -> v with edge (a1, b1) comes strictly
/// before its crossing with edge (a2, b2); sigma_j = orient(o, v, a_j).
/// Compares the parameters s_j = n_j / d_j, where d_j has sign -sigma_j. Exact
/// in the exact kernel; for robust floating-point kernels a filter with a
/// conservative error bound, falling back to exact rationals near ties.
template <class K, class FT>
bool crossing_before(const Point<FT>& o, const Point<FT>& v, const Point<FT>& a1, const Point<FT>& b1,
                     int sigma1, const Point<FT>& a2, const Point<FT>& b2, int sigma2)
{
  auto sign_of_det = [&](auto X, auto conv) {
    using T = decltype(X);
    const T ox = conv(o.x), oy = conv(o.y), dx = conv(v.x) - ox, dy = conv(v.y) - oy;
    const T e1x = conv(b1.x) - conv(a1.x), e1y = conv(b1.y) - conv(a1.y);
    const T e2x = conv(b2.x) - conv(a2.x), e2y = conv(b2.y) - conv(a2.y);
    const T n1 = (conv(a1.x) - ox) * e1y - (conv(a1.y) - oy) * e1x, d1 = dx * e1y - dy * e1x;
    const T n2 = (conv(a2.x) - ox) * e2y - (conv(a2.y) - oy) * e2x, d2 = dx * e2y - dy * e2x;
    const T det = n1 * d2 - n2 * d1;
    return (det > 0) - (det < 0);
  };
  int s;  // sign of n1 d2 - n2 d1
  if constexpr (is_exact_kernel<K>::value || !K::robust) {
    s = sign_of_det(FT(0), [](const FT& x) { return x; });
  } else {
    const double ox = double(o.x), oy = double(o.y), dx = double(v.x) - ox, dy = double(v.y) - oy;
    const double e1x = double(b1.x) - double(a1.x), e1y = double(b1.y) - double(a1.y);
    const double e2x = double(b2.x) - double(a2.x), e2y = double(b2.y) - double(a2.y);
    const double p1 = (double(a1.x) - ox) * e1y, q1 = (double(a1.y) - oy) * e1x;
    const double p2 = (double(a2.x) - ox) * e2y, q2 = (double(a2.y) - oy) * e2x;
    const double r1 = dx * e1y, t1 = dy * e1x, r2 = dx * e2y, t2 = dy * e2x;
    const double n1 = p1 - q1, d1 = r1 - t1, n2 = p2 - q2, d2 = r2 - t2;
    const double det = n1 * d2 - n2 * d1;
    // Each of n_j, d_j carries a relative error of a few ulps of |p| + |q|;
    // 32 * 2^-53 bounds the error of det with room to spare.
    const double bound = 32 * 1.1102230246251565e-16 *
                         ((std::abs(p1) + std::abs(q1)) * (std::abs(r2) + std::abs(t2)) +
                          (std::abs(p2) + std::abs(q2)) * (std::abs(r1) + std::abs(t1)));
    if (det > bound || -det > bound) {
      s = (det > 0) - (det < 0);
    } else {
#ifdef PD_HAVE_GMP
      s = sign_of_det(mpq_class(0), [](const FT& x) { return mpq_class(double(x)); });
#else
      s = (det > 0) - (det < 0);  // no exact fallback without GMP
#endif
    }
  }
  // s1 < s2  <=>  (n1 d2 - n2 d1) / (d1 d2) < 0, and sign(d1 d2) = sigma1 * sigma2.
  return s * sigma1 * sigma2 < 0;
}

struct Corner {
  std::uint32_t v;     ///< pool id
  std::uint32_t line;  ///< input line containing the edge from this vertex to the next
  bool reflex;
};

} // namespace detail

template <class K>
ConvexDecomposition<typename K::FT> convex_decompose(Ring<typename K::FT> ring,
                                                     const ConvexOptions& opt = {})
{
  using FT = typename K::FT;
  using detail::Corner;
  ConvexDecomposition<FT> out;
  const Status broken = is_exact_kernel<K>::value ? Status::InternalError : Status::NumericalFailure;

  // ---- input
  if (opt.validate_input) {
    ValidationReport rep = validate<K>(PolygonWithHoles<FT>{ring, {}});
    if (!rep.ok()) {
      out.diag = {Status::InvalidInput, "input is not a simple polygon", std::move(rep.issues)};
      return out;
    }
  } else if (ring.size() < 3) {
    out.diag = {Status::InvalidInput, "input has fewer than 3 vertices", {}};
    return out;
  }
  const int orientation = ring_orientation<K>(ring);
  if (orientation == 0) {
    out.diag = {Status::InvalidInput, "input ring is degenerate", {}};
    return out;
  }
  if (orientation < 0) std::reverse(ring.begin(), ring.end());

  const std::size_t n = ring.size();
  double tol = 0;  // floating-point kernels: tolerance for constructed points
  if constexpr (!is_exact_kernel<K>::value) {
    double scale = 1;
    for (const Point<FT>& v : ring) scale = std::max({scale, std::abs(double(v.x)), std::abs(double(v.y))});
    tol = 1e-9 * scale;
  }
  out.pool = std::move(ring);
  const std::vector<Point<FT>>& pool = out.pool;  // grows; never hold references across a push
  std::vector<detail::Line<FT>> lines(n);
  for (std::size_t j = 0; j < n; ++j) {
    const Point<FT>& a = out.pool[j];
    const Point<FT>& b = out.pool[next_index(j, n)];
    lines[j] = {a, {b.x - a.x, b.y - a.y}};
  }

  auto is_reflex = [&](std::uint32_t a, std::uint32_t b, std::uint32_t c) {
    const Point<FT>& pa = out.pool[a];
    const Point<FT>& pb = out.pool[b];
    const Point<FT>& pc = out.pool[c];
    if (K::orient(pa, pb, pc) >= 0) return false;
    if constexpr (K::angle_tolerance > 0) {
      // Nearly straight: |sin(turn)| <= tolerance.
      const double ux = double(pb.x) - double(pa.x), uy = double(pb.y) - double(pa.y);
      const double vx = double(pc.x) - double(pb.x), vy = double(pc.y) - double(pb.y);
      const double dot = ux * vx + uy * vy;
      if (dot > 0 && std::abs(ux * vy - uy * vx) <= K::angle_tolerance * std::hypot(ux, uy) * std::hypot(vx, vy))
        return false;
    }
    return true;
  };

  std::vector<Corner> root(n);
  for (std::size_t j = 0; j < n; ++j) {
    const auto id = static_cast<std::uint32_t>(j);
    root[j] = {id, id,
               is_reflex(static_cast<std::uint32_t>(prev_index(j, n)), id,
                         static_cast<std::uint32_t>(next_index(j, n)))};
    out.stats.reflex_vertices += root[j].reflex;
  }
  out.stats.input_vertices = n;

  struct Work {
    std::vector<Corner> ring;
    int node;
    std::size_t depth;
  };
  std::vector<Work> stack;
  if (opt.keep_tree) out.tree.push_back(ConvexNode{});
  stack.push_back({std::move(root), opt.keep_tree ? 0 : -1, 0});

  while (!stack.empty()) {
    Work w = std::move(stack.back());
    stack.pop_back();
    const std::vector<Corner>& q = w.ring;
    const std::size_t m = q.size();
    out.stats.max_depth = std::max(out.stats.max_depth, w.depth);

    // ---- first reflex vertex, or a convex leaf
    std::size_t i = 0;
    while (i < m && !q[i].reflex) ++i;
    if (i == m) {
      std::vector<std::uint32_t> ids(m);
      for (std::size_t j = 0; j < m; ++j) ids[j] = q[j].v;
      if (w.node >= 0) out.tree[static_cast<std::size_t>(w.node)].leaf = static_cast<int>(out.leaves.size());
      out.leaves.push_back(std::move(ids));
      continue;
    }

    // ---- first boundary contact of the ray from v_{i-1} through v_i, beyond v_i.
    // All decisions use exact predicates on stored points (see detail::Hit).
    const std::size_t ip = prev_index(i, m);
    const Point<FT> o = pool[q[ip].v];
    const Point<FT> vi = pool[q[i].v];
    const Point<FT> dir{vi.x - o.x, vi.y - o.y};
    const std::size_t ex0 = prev_index(ip, m);  // edges incident to v_{i-1} or v_i
    // p1 strictly before p2 along the ray, for two points on the ray's line:
    // along a non-vertical line x is monotone, along a vertical one y is.
    auto along_less = [&](const Point<FT>& p1, const Point<FT>& p2) {
      if (dir.x != 0) return dir.x > 0 ? p1.x < p2.x : p2.x < p1.x;
      return dir.y > 0 ? p1.y < p2.y : p2.y < p1.y;
    };
    auto edge_a = [&](const detail::Hit& h) -> const Point<FT>& { return pool[q[h.pos].v]; };
    auto edge_b = [&](const detail::Hit& h) -> const Point<FT>& { return pool[q[next_index(h.pos, m)].v]; };
    auto before = [&](const detail::Hit& h1, const detail::Hit& h2) {
      if (h1.vertex && h2.vertex) return along_less(pool[q[h1.pos].v], pool[q[h2.pos].v]);
      if (!h1.vertex && h2.vertex) return K::orient(edge_a(h1), edge_b(h1), pool[q[h2.pos].v]) == h1.sigma;
      if (h1.vertex && !h2.vertex) return K::orient(edge_a(h2), edge_b(h2), pool[q[h1.pos].v]) == -h2.sigma;
      return detail::crossing_before<K>(o, vi, edge_a(h1), edge_b(h1), h1.sigma, edge_a(h2), edge_b(h2),
                                        h2.sigma);
    };
    bool found = false;
    detail::Hit best{};
    auto offer = [&](const detail::Hit& h) {
      if (!found || before(h, best)) {
        found = true;
        best = h;
      }
    };
    for (std::size_t j = 0; j < m; ++j) {
      if (j == ex0 || j == ip || j == i) continue;
      const std::size_t jn = next_index(j, m);
      const Point<FT>& a = pool[q[j].v];
      const Point<FT>& b = pool[q[jn].v];
      const int oa = K::orient(o, vi, a), ob = K::orient(o, vi, b);
      if (oa == 0 && along_less(vi, a)) offer({true, j, 0});
      if (ob == 0 && along_less(vi, b)) offer({true, jn, 0});
      // A proper crossing of the ray's line lies beyond v_i iff orient(a, b, v_i) == -oa.
      if (oa * ob < 0 && K::orient(a, b, vi) == -oa) offer({false, j, oa});
    }
    if (!found) {
      out.diag = {broken, "no boundary contact beyond a reflex vertex", {}};
      return out;
    }
    const bool hit_is_vertex = best.vertex;
    const std::size_t hit_pos = best.pos;

    // ---- split. In exact arithmetic every split lowers the total reflex count,
    // so there are at most n_c; the budget only stops floating-point runaways.
    if (out.stats.splits >= 2 * n + 16) {
      out.diag = {broken, "split budget exceeded", {}};
      return out;
    }
    std::vector<Corner> p1, p2;
    const Corner cprev = q[ip];
    std::uint32_t chord_to;
    if (!hit_is_vertex) {
      const std::size_t k = hit_pos;
      chord_to = static_cast<std::uint32_t>(out.pool.size());
      // From the two input lines (bounded size, one rounding); if rounding made
      // them parallel, from the stored points, which lie strictly on opposite
      // sides of the ray.
      Point<FT> p = detail::intersect_lines(lines[cprev.line], lines[q[k].line]);
      if constexpr (!is_exact_kernel<K>::value) {
        const Point<FT>& ka = pool[q[k].v];
        const Point<FT>& kb = pool[q[next_index(k, m)].v];
        if (!detail::point_finite(p) || !detail::near_box(ka, kb, p, tol))
          p = detail::intersect_through(o, vi, ka, kb);
      }
      if (!detail::point_finite(p)) {
        out.diag = {broken, "a constructed point is not finite", {}};
        return out;
      }
      out.pool.push_back(std::move(p));
      for (std::size_t pos = next_index(k, m); pos != i; pos = next_index(pos, m)) p1.push_back(q[pos]);
      p1.push_back({chord_to, q[k].line, false});
      for (std::size_t pos = i;; pos = next_index(pos, m)) {
        p2.push_back(q[pos]);
        if (pos == k) break;
      }
      p2.push_back({chord_to, cprev.line, false});
    } else {
      const std::size_t wp = hit_pos;
      chord_to = q[wp].v;
      for (std::size_t pos = next_index(wp, m); pos != i; pos = next_index(pos, m)) p1.push_back(q[pos]);
      p1.push_back(q[wp]);
      for (std::size_t pos = i;; pos = next_index(pos, m)) {
        p2.push_back(q[pos]);
        if (pos == wp) break;
      }
      p2.back().line = cprev.line;  // w's outgoing edge is now the chord (w, v_i)
      ++out.stats.vertex_hits;
    }
    if (p1.size() < 3 || p2.size() < 3) {
      out.diag = {broken, "a split produced a piece with fewer than 3 vertices", {}};
      return out;
    }
    // Recompute the flags of the vertices next to the chord (its ends and their
    // neighbours in each piece). In exact arithmetic this gives what the theory
    // predicts (v_i and p convex, the hit vertex's two new angles, the rest
    // unchanged); with rounded constructions it keeps the flags consistent with
    // the stored points.
    auto refresh = [&](std::vector<Corner>& c, std::size_t j) {
      c[j].reflex = is_reflex(c[prev_index(j, c.size())].v, c[j].v, c[next_index(j, c.size())].v);
    };
    for (std::vector<Corner>* piece : {&p1, &p2}) {
      refresh(*piece, 0);
      refresh(*piece, piece->size() - 2);
      refresh(*piece, piece->size() - 1);
    }
    if (hit_is_vertex)
      out.stats.double_resolutions += (q[hit_pos].reflex && !p1.back().reflex && !p2.back().reflex);

    if (opt.check) {
      for (const std::vector<Corner>* piece : {&p1, &p2}) {
        const std::vector<Corner>& c = *piece;
        for (std::size_t j = 0; j < c.size(); ++j)
          if (c[j].reflex !=
              is_reflex(c[prev_index(j, c.size())].v, c[j].v, c[next_index(j, c.size())].v)) {
            out.diag = {broken, "an inherited reflex flag differs from the recomputed one", {}};
            return out;
          }
      }
    }

    // ---- record and recurse (P1 first, as in a depth-first recursion)
    ++out.stats.splits;
    int n1 = -1, n2 = -1;
    if (w.node >= 0) {
      n1 = static_cast<int>(out.tree.size());
      n2 = n1 + 1;
      ConvexNode child;
      child.parent = w.node;
      out.tree.push_back(child);
      out.tree.push_back(child);
      ConvexNode& nd = out.tree[static_cast<std::size_t>(w.node)];
      nd.child[0] = n1;
      nd.child[1] = n2;
      nd.reflex_pos = static_cast<std::uint32_t>(i);
      nd.vertex_hit = hit_is_vertex;
      nd.chord_from = q[i].v;
      nd.chord_to = chord_to;
    }
    const std::size_t depth = w.depth + 1;
    stack.push_back({std::move(p2), n2, depth});
    stack.push_back({std::move(p1), n1, depth});
  }
  return out;
}

} // namespace pd
