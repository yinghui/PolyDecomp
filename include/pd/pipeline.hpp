// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#pragma once
/// \file
/// The pipeline: a polygon with holes -> simple hole-free pieces (Stage 1,
/// eliminate_holes) -> convex pieces (Stage 2, convex_decompose on each piece).
/// Stage 1 emits simple counter-clockwise pieces, exactly the input class of
/// Stage 2, so the stages compose; pieces that are already convex pass through.

#include <algorithm>
#include <cstdint>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "pd/convex_decomp.hpp"
#include "pd/hole_elim.hpp"
#include "pd/status.hpp"
#include "pd/verify.hpp"
#ifdef PD_HAVE_GMP
#include "pd/kernel/exact.hpp"
#endif

namespace pd {

/// What a floating-point run does when it detects a numerical failure.
enum class OnFailure {
  Report,      ///< return the NumericalFailure diagnostics
  RetryExact,  ///< rerun the whole input with ExactKernel and return that result, rounded
};

struct PipelineOptions {
  bool validate_input = true;  ///< check the input class first
  bool check = false;          ///< extra invariant checks in both stages
  bool keep_tree = true;       ///< keep the Stage 1 decomposition tree
  /// Floating-point kernels only: check the output's necessary conditions
  /// (check_pieces_approx) and report a NumericalFailure if they fail.
  bool self_check = true;
  OnFailure on_failure = OnFailure::Report;  ///< floating-point kernels only
};

struct PipelineStats {
  HoleElimStats holes;                ///< Stage 1
  std::size_t hole_free_pieces = 0;   ///< Stage 1 output
  std::size_t already_convex = 0;     ///< hole-free pieces with no reflex vertex
  std::size_t reflex_vertices = 0;    ///< total over the hole-free pieces
  std::size_t splits = 0;             ///< Stage 2, totals over the hole-free pieces
  std::size_t vertex_hits = 0;
  std::size_t double_resolutions = 0;
  std::size_t convex_pieces = 0;
  bool retried_exact = false;         ///< the result comes from an exact rerun
};

template <class FT>
struct PipelineResult {
  std::vector<Point<FT>> pool;                        ///< all vertices, shared by both stages
  std::vector<std::vector<std::uint32_t>> hole_free;  ///< Stage 1 pieces (CCW pool ids)
  std::vector<std::vector<std::uint32_t>> leaves;     ///< convex pieces (CCW pool ids)
  std::vector<std::uint32_t> leaf_source;             ///< hole_free index each convex piece came from
  std::vector<HoleNode<FT>> hole_tree;                ///< Stage 1 tree (if kept)
  PipelineStats stats;
  Diagnostics diag;

  bool ok() const { return diag.status == Status::Ok; }
  Ring<FT> ring_of(const std::vector<std::uint32_t>& ids) const
  {
    Ring<FT> r;
    r.reserve(ids.size());
    for (std::uint32_t id : ids) r.push_back(pool[id]);
    return r;
  }
  Ring<FT> leaf_ring(std::size_t i) const { return ring_of(leaves[i]); }
  std::vector<Ring<FT>> leaf_rings() const
  {
    std::vector<Ring<FT>> out;
    for (const auto& l : leaves) out.push_back(ring_of(l));
    return out;
  }
  std::vector<Ring<FT>> hole_free_rings() const
  {
    std::vector<Ring<FT>> out;
    for (const auto& l : hole_free) out.push_back(ring_of(l));
    return out;
  }
};

namespace detail {

template <class K>
PipelineResult<typename K::FT> decompose_once(PolygonWithHoles<typename K::FT> input,
                                              const PipelineOptions& opt)
{
  using FT = typename K::FT;
  PipelineResult<FT> out;

  HoleElimOptions s1;
  s1.validate_input = opt.validate_input;
  s1.check = opt.check;
  s1.keep_tree = opt.keep_tree;
  const PolygonWithHoles<FT> original = input;
  HoleElimination<FT> h = eliminate_holes<K>(std::move(input), s1);
  out.stats.holes = h.stats;
  if (!h.ok()) {
    out.diag = std::move(h.diag);
    return out;
  }
  out.pool = std::move(h.pool);
  out.hole_free = std::move(h.leaves);
  out.hole_tree = std::move(h.tree);
  out.stats.hole_free_pieces = out.hole_free.size();

  ConvexOptions s2;
  s2.validate_input = opt.check;  // Stage 1 output is simple by construction
  s2.check = opt.check;
  s2.keep_tree = false;
  for (std::size_t i = 0; i < out.hole_free.size(); ++i) {
    const std::vector<std::uint32_t>& ids = out.hole_free[i];
    Ring<FT> ring = out.ring_of(ids);
    if (ring_orientation<K>(ring) <= 0) {
      out.diag = {is_exact_kernel<K>::value ? Status::InternalError : Status::NumericalFailure,
                  "a hole-free piece is not counter-clockwise", {}};
      return out;
    }
    ConvexDecomposition<FT> c = convex_decompose<K>(std::move(ring), s2);
    if (!c.ok()) {
      out.diag = std::move(c.diag);
      out.diag.message = "hole-free piece " + std::to_string(i) + ": " + out.diag.message;
      return out;
    }
    // Local ids below ids.size() are the piece's own vertices (in order, as the
    // piece is already counter-clockwise); the rest are new points.
    std::vector<std::uint32_t> map(c.pool.size());
    for (std::size_t k = 0; k < c.pool.size(); ++k) {
      if (k < ids.size()) {
        map[k] = ids[k];
      } else {
        map[k] = static_cast<std::uint32_t>(out.pool.size());
        out.pool.push_back(std::move(c.pool[k]));
      }
    }
    for (const std::vector<std::uint32_t>& leaf : c.leaves) {
      std::vector<std::uint32_t> g;
      g.reserve(leaf.size());
      for (std::uint32_t k : leaf) g.push_back(map[k]);
      out.leaves.push_back(std::move(g));
      out.leaf_source.push_back(static_cast<std::uint32_t>(i));
    }
    out.stats.already_convex += (c.stats.reflex_vertices == 0);
    out.stats.reflex_vertices += c.stats.reflex_vertices;
    out.stats.splits += c.stats.splits;
    out.stats.vertex_hits += c.stats.vertex_hits;
    out.stats.double_resolutions += c.stats.double_resolutions;
  }
  out.stats.convex_pieces = out.leaves.size();
  if constexpr (!is_exact_kernel<K>::value) {
    if (opt.self_check) {
      const VerifyReport rep = check_pieces_approx<K>(original, out.leaf_rings(), true,
                                                      std::max(K::angle_tolerance, 1e-12), 1e-9);
      if (!rep.ok) out.diag = {Status::NumericalFailure, "self-check: " + rep.failure, {}};
    }
  }
  return out;
}

} // namespace detail

/// Decomposes a polygon with holes into convex pieces. With a floating-point
/// kernel and opt.on_failure == RetryExact, a numerical failure (reported by an
/// invariant check or the self-check) is answered by an exact rerun whose
/// result is rounded to doubles (requires GMP; stats.retried_exact is set).
template <class K>
PipelineResult<typename K::FT> decompose(PolygonWithHoles<typename K::FT> input,
                                         const PipelineOptions& opt = {})
{
  using FT = typename K::FT;
  if constexpr (!is_exact_kernel<K>::value && std::is_same_v<FT, double>) {
#ifdef PD_HAVE_GMP
    if (opt.on_failure == OnFailure::RetryExact) {
      PipelineResult<FT> out = detail::decompose_once<K>(input, opt);
      if (out.diag.status != Status::NumericalFailure) return out;
      PipelineResult<mpq_class> ex = detail::decompose_once<ExactKernel>(convert_polygon<ExactKernel>(input), opt);
      PipelineResult<FT> r;
      r.pool.reserve(ex.pool.size());
      for (const auto& p : ex.pool) r.pool.push_back({p.x.get_d(), p.y.get_d()});
      r.hole_free = std::move(ex.hole_free);
      r.leaves = std::move(ex.leaves);
      r.leaf_source = std::move(ex.leaf_source);
      r.stats = ex.stats;
      r.stats.retried_exact = true;
      r.diag = std::move(ex.diag);
      for (const auto& n : ex.hole_tree) {  // rounded copy of the Stage 1 tree
        HoleNode<FT> m;
        m.parent = n.parent;
        m.children = n.children;
        m.leaf = n.leaf;
        m.holes = n.holes;
        m.has_cut = n.has_cut;
        m.cut_point = {n.cut_point.x.get_d(), n.cut_point.y.get_d()};
        m.cut_dir = {n.cut_dir.x.get_d(), n.cut_dir.y.get_d()};
        m.fallback = n.fallback;
        m.pierced = n.pierced;
        r.hole_tree.push_back(std::move(m));
      }
      return r;
    }
#endif
  }
  return detail::decompose_once<K>(std::move(input), opt);
}

} // namespace pd
