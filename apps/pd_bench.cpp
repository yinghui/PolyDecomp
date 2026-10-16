// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
// pd_bench: piece counts and timings of polydecomp, of PolyPartition's
// Hertel-Mehlhorn partition, and of CGAL's convex decompositions when built
// with PD_WITH_CGAL.
//
//   pd_bench [--suite simple|holes|all|scale|slice|clipper-scale] [--repeats R] [--csv FILE] [--quick]
//            [--budget-ms T] [--cuts N] [--baseline-budget-s S]
//
// Every (input, method) pair gives one CSV row:
//   suite,input,n,holes,reflex,lower_bound,method,status,pieces,median_ms,certified
// n counts all vertices, reflex the reflex vertices of the outer ring and holes
// (for simple polygons lower_bound is Chazelle's ceil(r/2) + 1, otherwise 0).
// Times are medians of R runs after one untimed run. That untimed run is checked
// by the exact certificate (verify_partition with convexity required) where this
// is possible: certified is 1 if it passes, 0 if it fails, and empty for methods
// whose pieces contain rounded constructed points (ours-double, cgal-bisector,
// cgal-vertical) or are not convex (ours-holes-exact). Baseline pieces are made
// counterclockwise first.
//
// The scale suite (not part of "all") runs random star polygons and combs with
// 2^10 to 2^16 vertices and random_holes with 2^4 to 2^14 holes, in increasing
// size. Within each of these three series a method is run on the next size only
// while its median time stays at most T ms (default 5000); larger sizes of a
// method that exceeded it get no row. The scale suite also times input
// validation, which decompose() performs by default: validate-double and
// validate-exact (pieces 0).
//
// The slice suite (not part of "all") cuts shapes with pd::Slicer, N cuts per sequence (default
// 200): cheese with 16, 64 and 256 holes, random_holes with 64 holes and the bundled sample polygon,
// each with the sequences random (lines through two random points of the bounding box), axis
// (alternately vertical and horizontal lines through random input vertices), vertex (lines through
// two random input vertices), repeat (every other cut repeats a random earlier line) and, for the
// cheese, graze (the tangent at a random vertex of a random hole, through that vertex and parallel
// to the chord of its neighbours). One row per cut:
//   shape,n,holes,method,sequence,cut,pieces,cut_pieces,vertex_contacts,edge_contacts,dangles,
//   symbolic,cut_ms,convex_ms,convex_pieces,max_bits,certified,tiny,area_err,split_pieces
// cut_ms times the cut alone, convex_ms the convex decomposition of the pieces it created, and
// certified is the exact certificate (verify_slices) of all pieces after the cut.
// The same lines then cut the shape with two baselines (column method): clipper2, Clipper2's
// Clipper64 on the shape scaled by a power of two to 64-bit integers below 2^40 (vendored), and,
// with CGAL, cgal-boolean, CGAL's Boolean set operations on Epeck with segment traits. Both cut a
// piece by computing piece ∩ H and piece − H for the half-plane H left of the line (a large
// rectangle clipped by the line) and skip pieces whose bounding box lies on one side. For clipper2,
// certified checks the pieces against Clipper2's own rounded input after the first and the last
// cut, tiny counts
// pieces with less than 1e-9 of the input's area, and area_err is the relative difference between
// the pieces' total area and that input's (exact, in integers). cgal-boolean stops once its cut
// times add up to the budget (--baseline-budget-s, default 120 s); its lazy numbers are evaluated
// only after the last cut (certified, max_bits, tiny, and split_pieces, the pieces counted as the
// slicer reports them: a hole touching the outer ring merged into it, two regions touching at a
// point split; see CgalSlicer::exact_pieces).
//   ...,certified,tiny,area_err,split_pieces
// The clipper-scale suite cuts cheese with 16 and 64 holes 60 times per sequence with Clipper2 at
// the scales 2^20, 2^26, 2^32, 2^40 and 2^50, against the slicer's piece count:
//   shape,sequence,cuts,bits,exact_pieces,pieces,tiny,area_err
//
// Methods: ours-exact / ours-double (the full pipeline, exact or double
// arithmetic; the double run retries exactly on a numerical failure),
// ours-holes-exact (Stage 1 only: hole-free pieces), polypartition-hm
// (PolyPartition's ConvexPartition_HM: holes bridged, ear clipping, Hertel and
// Mehlhorn's merge, in double precision). With CGAL:
// cgal-hertel-mehlhorn (approx_convex_partition_2), cgal-greene-approx,
// cgal-optimal (Greene's O(n^4) optimum, n <= 64), cgal-bisector
// (Small_side_angle_bisector_decomposition_2), cgal-cdt (constrained Delaunay
// triangles), cgal-cdt-hm (Hertel and Mehlhorn's merge applied to those
// triangles, polygons with or without holes) and, for polygons with holes,
// cgal-vertical (vertical decomposition).

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iostream>
#include <list>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "pd/generate.hpp"
#include "pd/hole_elim.hpp"
#include "pd/io/obj.hpp"
#include "pd/kernel/double.hpp"
#include "pd/kernel/exact.hpp"
#include "pd/pipeline.hpp"
#include "pd/slice.hpp"
#include "pd/verify.hpp"

#include "clipper2/clipper.engine.h"
#include "polypartition.h"

#ifdef PD_HAVE_CGAL
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_triangulation_decomposition_2.h>
#include <CGAL/Polygon_vertical_decomposition_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Small_side_angle_bisector_decomposition_2.h>
#include <CGAL/partition_2.h>
#endif

namespace {

struct Input {
  std::string suite, name;
  pd::PolygonWithHoles<double> poly;
  std::string series = {};  ///< scale suite: star, comb or holes (sizes increase within a series)
};

struct Result {
  bool ok = false;
  std::size_t pieces = 0;
  int certified = -1;  ///< 1 / 0: passes / fails the exact certificate, -1: not checked
};

/// A method runs on an input. With certify = true (the untimed run) it also
/// sets Result::certified where its pieces can be checked exactly.
using Method = std::pair<std::string, std::function<Result(const pd::PolygonWithHoles<double>&, bool)>>;

/// Median of `repeats` timed runs (the caller makes the untimed first run).
double median_ms(int repeats, const std::function<void()>& f)
{
  std::vector<double> t;
  for (int r = 0; r < repeats; ++r) {
    const auto t0 = std::chrono::steady_clock::now();
    f();
    t.push_back(std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count());
  }
  std::sort(t.begin(), t.end());
  return t[t.size() / 2];
}

std::size_t reflex_count(const pd::PolygonWithHoles<double>& p)
{
  using DK = pd::DoubleKernel;
  auto count = [](pd::Ring<double> r, bool outer) {
    // Outer ring: reflex = right turn when CCW. Hole: the region is outside,
    // so a reflex vertex of the region is a convex vertex of the hole.
    pd::make_ccw<DK>(r);
    std::size_t c = 0;
    for (std::size_t i = 0; i < r.size(); ++i) {
      const int o = DK::orient(r[pd::prev_index(i, r.size())], r[i], r[pd::next_index(i, r.size())]);
      c += outer ? (o < 0) : (o > 0);
    }
    return c;
  };
  std::size_t r = count(p.outer, true);
  for (const auto& h : p.holes) r += count(h, false);
  return r;
}

/// The exact certificate for pieces whose vertices are input points, given in
/// double precision (converted exactly). Pieces are made counterclockwise first.
int certify_rings(const pd::PolygonWithHoles<double>& input, std::vector<pd::Ring<double>> rings)
{
  using EK = pd::ExactKernel;
  std::vector<pd::Ring<mpq_class>> pieces;
  pieces.reserve(rings.size());
  for (pd::Ring<double>& r : rings) {
    pd::make_ccw<pd::DoubleKernel>(r);
    pieces.push_back(pd::convert_ring<EK>(r));
  }
  return pd::verify_partition<EK>(pd::convert_polygon<EK>(input), pieces, true).ok ? 1 : 0;
}

std::vector<Method> our_methods()
{
  using EK = pd::ExactKernel;
  std::vector<Method> m;
  m.push_back({"ours-exact", [](const pd::PolygonWithHoles<double>& p, bool certify) {
                 pd::PipelineOptions o;
                 o.validate_input = false;
                 o.keep_tree = false;
                 const auto in = pd::convert_polygon<EK>(p);
                 const auto r = pd::decompose<EK>(in, o);
                 Result res{r.ok(), r.leaves.size()};
                 if (certify && r.ok()) res.certified = pd::verify_partition<EK>(in, r.leaf_rings(), true).ok ? 1 : 0;
                 return res;
               }});
  m.push_back({"ours-double", [](const pd::PolygonWithHoles<double>& p, bool) {
                 pd::PipelineOptions o;
                 o.validate_input = false;
                 o.keep_tree = false;
                 o.on_failure = pd::OnFailure::RetryExact;
                 const auto r = pd::decompose<pd::DoubleKernel>(p, o);
                 return Result{r.ok(), r.leaves.size()};
               }});
  m.push_back({"polypartition-hm", [](const pd::PolygonWithHoles<double>& p, bool certify) {
                 // Outer ring counterclockwise, holes clockwise and marked, as PolyPartition requires.
                 TPPLPolyList in, parts;
                 auto add = [&in](pd::Ring<double> r, bool hole) {
                   if (hole) pd::make_cw<pd::DoubleKernel>(r);
                   else pd::make_ccw<pd::DoubleKernel>(r);
                   TPPLPoly t;
                   t.Init(static_cast<long>(r.size()));
                   for (std::size_t i = 0; i < r.size(); ++i) {
                     t[static_cast<int>(i)].x = r[i].x;
                     t[static_cast<int>(i)].y = r[i].y;
                   }
                   t.SetHole(hole);
                   in.push_back(t);
                 };
                 add(p.outer, false);
                 for (const auto& h : p.holes) add(h, true);
                 TPPLPartition pp;
                 Result res{pp.ConvexPartition_HM(&in, &parts) == 1, parts.size()};
                 if (certify && res.ok) {
                   std::vector<pd::Ring<double>> rings;
                   for (TPPLPoly& t : parts) {
                     pd::Ring<double> r;
                     for (long i = 0; i < t.GetNumPoints(); ++i) r.push_back({t[static_cast<int>(i)].x, t[static_cast<int>(i)].y});
                     rings.push_back(std::move(r));
                   }
                   res.certified = certify_rings(p, std::move(rings));
                 }
                 return res;
               }});
  return m;
}

#ifdef PD_HAVE_CGAL
using Epick = CGAL::Exact_predicates_inexact_constructions_kernel;
using Epeck = CGAL::Exact_predicates_exact_constructions_kernel;

template <class K>
CGAL::Polygon_2<K> cgal_ring(pd::Ring<double> r, bool ccw)
{
  if (ccw) pd::make_ccw<pd::DoubleKernel>(r);
  else pd::make_cw<pd::DoubleKernel>(r);
  CGAL::Polygon_2<K> out;
  for (const auto& p : r) out.push_back(typename K::Point_2(p.x, p.y));
  return out;
}

template <class K>
CGAL::Polygon_with_holes_2<K> cgal_polygon(const pd::PolygonWithHoles<double>& p)
{
  CGAL::Polygon_with_holes_2<K> out(cgal_ring<K>(p.outer, true));
  for (const auto& h : p.holes) out.add_hole(cgal_ring<K>(h, false));
  return out;
}

template <class Range>
std::vector<pd::Ring<double>> rings_of(const Range& polygons)
{
  std::vector<pd::Ring<double>> rings;
  for (const auto& poly : polygons) {
    pd::Ring<double> r;
    for (auto v = poly.vertices_begin(); v != poly.vertices_end(); ++v) r.push_back({v->x(), v->y()});
    rings.push_back(std::move(r));
  }
  return rings;
}

template <class F>
Method partition_method(const char* name, F partition)
{
  return {name, [partition](const pd::PolygonWithHoles<double>& p, bool certify) {
            using Traits = CGAL::Partition_traits_2<Epick>;
            const auto poly = cgal_ring<Epick>(p.outer, true);
            std::list<Traits::Polygon_2> parts;
            partition(poly.vertices_begin(), poly.vertices_end(), std::back_inserter(parts));
            Result res{true, parts.size()};
            if (certify) res.certified = certify_rings(p, rings_of(parts));
            return res;
          }};
}

/// Hertel and Mehlhorn's merge on a triangulation whose vertices are input
/// points. Every diagonal (an edge shared by two pieces) is removed, in the order
/// in which it is found (by triangle, then by edge), if both of its endpoints stay
/// convex (a straight angle counts as convex); the test is CGAL's exact orientation
/// on the stored doubles. The pieces are kept as linked half-edges, so each removal
/// takes constant time and the merge is linear after the triangulation.
std::vector<std::vector<std::size_t>> hertel_mehlhorn(const std::vector<Epick::Point_2>& pts,
                                                      const std::vector<std::vector<std::size_t>>& tris)
{
  std::vector<std::size_t> origin, next, prev, piece;  // per half-edge
  std::unordered_map<std::uint64_t, std::size_t> edge;  // directed edge (u, v) -> half-edge
  const auto key = [&pts](std::size_t u, std::size_t v) { return static_cast<std::uint64_t>(u) * pts.size() + v; };
  for (std::size_t k = 0; k < tris.size(); ++k) {
    const std::size_t base = origin.size(), m = tris[k].size();
    for (std::size_t i = 0; i < m; ++i) {
      origin.push_back(tris[k][i]);
      next.push_back(base + (i + 1) % m);
      prev.push_back(base + (i + m - 1) % m);
      piece.push_back(k);
      edge[key(tris[k][i], tris[k][(i + 1) % m])] = base + i;
    }
  }
  const auto head = [&](std::size_t e) { return origin[next[e]]; };
  std::vector<std::pair<std::size_t, std::size_t>> diagonals;  // (u -> v, v -> u) with u < v
  for (std::size_t e = 0; e < origin.size(); ++e) {
    const auto twin = edge.find(key(head(e), origin[e]));
    if (origin[e] < head(e) && twin != edge.end()) diagonals.push_back({e, twin->second});
  }
  std::vector<std::size_t> parent(tris.size());
  for (std::size_t k = 0; k < parent.size(); ++k) parent[k] = k;
  const auto find = [&parent](std::size_t k) {
    while (parent[k] != k) k = parent[k] = parent[parent[k]];
    return k;
  };
  std::vector<char> removed(origin.size(), 0);
  for (const auto& [e, t] : diagonals) {  // e = u -> v in piece A, t = v -> u in piece B
    const std::size_t a = find(piece[e]), b = find(piece[t]);
    if (a == b) continue;
    const std::size_t u = origin[e], v = origin[t];
    // Merged angles: at u between the vertex before u in A and the one after u in B,
    // at v between the vertex before v in B and the one after v in A.
    if (CGAL::orientation(pts[origin[prev[e]]], pts[u], pts[head(next[t])]) == CGAL::RIGHT_TURN ||
        CGAL::orientation(pts[origin[prev[t]]], pts[v], pts[head(next[e])]) == CGAL::RIGHT_TURN)
      continue;
    const std::size_t pe = prev[e], ne = next[e], pt = prev[t], nt = next[t];
    next[pe] = nt;
    prev[nt] = pe;
    next[pt] = ne;
    prev[ne] = pt;
    removed[e] = removed[t] = 1;
    parent[b] = a;
  }
  std::vector<std::vector<std::size_t>> out;
  std::vector<char> seen(origin.size(), 0);
  for (std::size_t e = 0; e < origin.size(); ++e) {
    if (removed[e] || seen[e]) continue;
    std::vector<std::size_t> ring;
    for (std::size_t f = e; !seen[f]; f = next[f]) {
      seen[f] = 1;
      ring.push_back(origin[f]);
    }
    out.push_back(std::move(ring));
  }
  return out;
}

std::vector<Method> cgal_methods(bool simple, std::size_t n)
{
  std::vector<Method> m;
  if (simple) {
    m.push_back(partition_method("cgal-hertel-mehlhorn", [](auto b, auto e, auto out) {
      CGAL::approx_convex_partition_2(b, e, out);
    }));
    m.push_back(partition_method("cgal-greene-approx", [](auto b, auto e, auto out) {
      CGAL::greene_approx_convex_partition_2(b, e, out);
    }));
    if (n <= 64)
      m.push_back(partition_method("cgal-optimal", [](auto b, auto e, auto out) {
        CGAL::optimal_convex_partition_2(b, e, out);
      }));
    m.push_back({"cgal-bisector", [](const pd::PolygonWithHoles<double>& p, bool) {
                   std::vector<CGAL::Polygon_2<Epeck>> parts;
                   CGAL::Small_side_angle_bisector_decomposition_2<Epeck> d;
                   d(cgal_ring<Epeck>(p.outer, true), std::back_inserter(parts));
                   return Result{true, parts.size()};
                 }});
  } else {
    m.push_back({"cgal-vertical", [](const pd::PolygonWithHoles<double>& p, bool) {
                   std::vector<CGAL::Polygon_2<Epeck>> parts;
                   CGAL::Polygon_vertical_decomposition_2<Epeck> d;
                   d(cgal_polygon<Epeck>(p), std::back_inserter(parts));
                   return Result{true, parts.size()};
                 }});
  }
  m.push_back({"cgal-cdt", [](const pd::PolygonWithHoles<double>& p, bool certify) {
                 std::vector<CGAL::Polygon_2<Epick>> parts;
                 CGAL::Polygon_triangulation_decomposition_2<Epick> d;
                 d(cgal_polygon<Epick>(p), std::back_inserter(parts));
                 Result res{true, parts.size()};
                 if (certify) res.certified = certify_rings(p, rings_of(parts));
                 return res;
               }});
  m.push_back({"cgal-cdt-hm", [](const pd::PolygonWithHoles<double>& p, bool certify) {
                 std::vector<CGAL::Polygon_2<Epick>> tris;
                 CGAL::Polygon_triangulation_decomposition_2<Epick> d;
                 d(cgal_polygon<Epick>(p), std::back_inserter(tris));
                 // Number the vertices (all are input points) and orient the triangles counterclockwise.
                 std::map<std::pair<double, double>, std::size_t> id;
                 std::vector<Epick::Point_2> pts;
                 std::vector<std::vector<std::size_t>> pieces;
                 for (auto& t : tris) {
                   if (t.orientation() == CGAL::CLOCKWISE) t.reverse_orientation();
                   std::vector<std::size_t> piece;
                   for (auto v = t.vertices_begin(); v != t.vertices_end(); ++v) {
                     const auto [it, fresh] = id.emplace(std::make_pair(v->x(), v->y()), pts.size());
                     if (fresh) pts.push_back(*v);
                     piece.push_back(it->second);
                   }
                   pieces.push_back(std::move(piece));
                 }
                 pieces = hertel_mehlhorn(pts, pieces);
                 Result res{true, pieces.size()};
                 if (certify) {
                   std::vector<pd::Ring<double>> rings;
                   for (const auto& piece : pieces) {
                     pd::Ring<double> r;
                     for (std::size_t k : piece) r.push_back({pts[k].x(), pts[k].y()});
                     rings.push_back(std::move(r));
                   }
                   res.certified = certify_rings(p, std::move(rings));
                 }
                 return res;
               }});
  return m;
}
#endif

std::vector<Input> simple_inputs(bool quick)
{
  std::vector<Input> in;
  const std::string data = PD_DATA_DIR;
  for (const char* name : {"cpolygon0", "polygon", "test5", "test1", "pentagon"})
    in.push_back({"simple", std::string("legacy-") + name,
                  {pd::read_obj_ring(data + "/legacy/simp/" + name + ".obj"), {}}});
  const std::vector<int> sizes = quick ? std::vector<int>{32, 128, 512} : std::vector<int>{32, 64, 128, 256, 512, 1024, 2048};
  for (int n : sizes)
    for (std::uint64_t s = 1; s <= (quick ? 1u : 3u); ++s)
      in.push_back({"simple", "star-" + std::to_string(n) + "-" + std::to_string(s), {pd::gen::star_polygon(n, s), {}}});
  const std::vector<int> teeth = quick ? std::vector<int>{16, 64} : std::vector<int>{16, 64, 256, 1024};
  for (int t : teeth) in.push_back({"simple", "comb-" + std::to_string(4 * t), {pd::gen::comb(t), {}}});
  return in;
}

std::vector<Input> hole_inputs(bool quick)
{
  std::vector<Input> in;
  const std::string data = PD_DATA_DIR;
  {
    pd::PolygonWithHoles<double> p{pd::read_obj_ring(data + "/legacy/compgon/cpolygon0.obj"), {}};
    for (int i = 0; i <= 6; ++i)
      p.holes.push_back(pd::read_obj_ring(data + "/legacy/compgon/hole0-" + std::to_string(i) + ".obj"));
    in.push_back({"holes", "legacy-bundled", p});
  }
  for (pd::gen::Family f : pd::gen::all_families())
    for (std::uint64_t s = 1; s <= (quick ? 2u : 10u); ++s)
      in.push_back({"holes", std::string(pd::gen::name(f)) + "-" + std::to_string(s), pd::gen::generate(f, s)});
  const std::vector<int> hs = quick ? std::vector<int>{16, 64} : std::vector<int>{8, 16, 32, 64, 128, 256};
  for (int h : hs) in.push_back({"holes", "spiral-" + std::to_string(h), pd::gen::spiral_holes(h, 1)});
  return in;
}

std::vector<Input> scale_inputs()
{
  std::vector<Input> in;
  for (int e = 10; e <= 16; ++e)
    in.push_back({"scale", "star-" + std::to_string(1 << e), {pd::gen::star_polygon(1 << e, 1), {}}, "star"});
  for (int e = 10; e <= 16; ++e)
    in.push_back({"scale", "comb-" + std::to_string(1 << e), {pd::gen::comb((1 << e) / 4), {}}, "comb"});
  for (int e = 4; e <= 14; ++e)
    in.push_back({"scale", "holes-" + std::to_string(1 << e), pd::gen::random_holes(1 << e, 1), "holes"});
  return in;
}

int usage()
{
  std::cerr << "usage: pd_bench [--suite simple|holes|all|scale|slice|clipper-scale] [--repeats R] [--csv FILE] [--quick]\n"
               "                [--budget-ms T] [--cuts N] [--baseline-budget-s S]\n";
  return 2;
}

/// Uniform in [0, 1) from the raw 64-bit output, the same on every platform.
double unit(std::mt19937_64& g) { return static_cast<double>(g() >> 11) * 0x1.0p-53; }

using Q = pd::Point<mpq_class>;
struct SliceLine {
  Q a, b;  ///< the cut is the line through a and b
};

/// The seeded cut sequence `seq` for a shape (see the file comment).
std::vector<SliceLine> slice_lines(const std::string& name, const std::string& seq,
                                   const pd::PolygonWithHoles<double>& shape, int cuts)
{
  auto exact = [](const pd::Point<double>& p) { return Q{mpq_class(p.x), mpq_class(p.y)}; };
  std::vector<pd::Point<double>> verts(shape.outer.begin(), shape.outer.end());
  for (const auto& h : shape.holes) verts.insert(verts.end(), h.begin(), h.end());
  double x0 = verts[0].x, x1 = x0, y0 = verts[0].y, y1 = y0;
  for (const auto& v : verts) {
    x0 = std::min(x0, v.x), x1 = std::max(x1, v.x), y0 = std::min(y0, v.y), y1 = std::max(y1, v.y);
  }
  std::mt19937_64 g(std::hash<std::string>{}(name + "/" + seq) & 0xffffffffu);
  std::vector<SliceLine> done;
  for (int k = 0; k < cuts; ++k) {
    Q a, b;
    auto random_point = [&]() {
      return pd::Point<double>{x0 + (x1 - x0) * (unit(g) * 1.1 - 0.05), y0 + (y1 - y0) * (unit(g) * 1.1 - 0.05)};
    };
    auto vertex = [&]() { return verts[static_cast<std::size_t>(unit(g) * verts.size())]; };
    if (seq == "random" || (seq == "repeat" && (k % 2 == 0 || done.empty()))) {
      a = exact(random_point());
      do b = exact(random_point()); while (b == a);
    } else if (seq == "repeat") {
      const SliceLine& l = done[static_cast<std::size_t>(unit(g) * done.size())];
      a = l.a, b = l.b;
    } else if (seq == "axis") {
      const auto v = vertex();
      a = exact(v);
      b = (k % 2 == 0) ? exact({v.x, v.y + 1}) : exact({v.x + 1, v.y});
    } else if (seq == "vertex") {
      a = exact(vertex());
      do b = exact(vertex()); while (b == a);
    } else {  // graze: tangent at a vertex of a (convex, regular) cheese hole
      const auto& h = shape.holes[static_cast<std::size_t>(unit(g) * shape.holes.size())];
      const std::size_t m = h.size(), j = static_cast<std::size_t>(unit(g) * m);
      const auto &v = h[j], &vp = h[(j + m - 1) % m], &vn = h[(j + 1) % m];
      a = exact(v);
      b = Q{a.x + mpq_class(vn.x) - mpq_class(vp.x), a.y + mpq_class(vn.y) - mpq_class(vp.y)};
    }
    done.push_back({a, b});
  }
  return done;
}

/// One CSV row of the slice suite; empty optional fields are written as empty cells.
struct SliceRow {
  std::string shape, method, seq;
  std::size_t n = 0, holes = 0;
  int cut = 0;
  std::size_t pieces = 0, cut_pieces = 0;
  std::string vertex_contacts, edge_contacts, dangles, symbolic, convex_ms, convex_pieces, max_bits,
      certified, tiny, area_err, split_pieces;
  double cut_ms = 0;
};
void write_row(std::ostream& csv, const SliceRow& r)
{
  csv << r.shape << ',' << r.n << ',' << r.holes << ',' << r.method << ',' << r.seq << ',' << r.cut << ','
      << r.pieces << ',' << r.cut_pieces << ',' << r.vertex_contacts << ',' << r.edge_contacts << ','
      << r.dangles << ',' << r.symbolic << ',' << r.cut_ms << ',' << r.convex_ms << ',' << r.convex_pieces << ','
      << r.max_bits << ',' << r.certified << ',' << r.tiny << ',' << r.area_err << ',' << r.split_pieces << '\n';
}

__extension__ typedef __int128 i128;  // exact areas of 64-bit integer polygons

/// Pieces whose area is below 1e-9 times the shape's area (slivers).
constexpr double kTiny = 1e-9;

double area_of(const pd::PolygonWithHoles<double>& p)
{
  double a = std::abs(pd::signed_area2(p.outer));
  for (const auto& h : p.holes) a -= std::abs(pd::signed_area2(h));
  return a / 2;
}

/// Cutting with Clipper2 (vendored): the shape is scaled by a power of two and rounded to 64-bit
/// integers once, and every cut computes piece ∩ H and piece − H for the half-plane H left of the line,
/// a large rectangle clipped by the line with its two corners on the line rounded once.
struct ClipperSlicer {
  using P64 = Clipper2Lib::Point64;
  using Piece = Clipper2Lib::Paths64;  ///< outer ring first, then the holes
  double scale = 1;
  std::vector<Piece> pieces;
  std::vector<Clipper2Lib::Rect64> boxes;
  double bx0 = 0, by0 = 0, bx1 = 0, by1 = 0;  // the rectangle that H is cut from

  P64 pt(double x, double y) const { return P64(std::llround(x * scale), std::llround(y * scale)); }
  Clipper2Lib::Path64 path(const pd::Ring<double>& r, bool positive) const
  {
    Clipper2Lib::Path64 p;
    for (const auto& v : r) p.push_back(pt(v.x, v.y));
    if ((area2(p) > 0) != positive) std::reverse(p.begin(), p.end());
    return p;
  }
  static i128 area2(const Clipper2Lib::Path64& p)
  {
    i128 a = 0;
    for (std::size_t i = 0, j = p.size() - 1; i < p.size(); j = i++)
      a += static_cast<i128>(p[j].x) * p[i].y - static_cast<i128>(p[i].x) * p[j].y;
    return a;
  }
  static i128 area2(const Piece& p)
  {
    i128 a = 0;
    for (const auto& r : p) a += area2(r);  // holes are negative
    return a;
  }
  explicit ClipperSlicer(const pd::PolygonWithHoles<double>& s, int bits = 40)
  {
    double m = 0;
    bx0 = bx1 = s.outer[0].x, by0 = by1 = s.outer[0].y;
    for (const auto& v : s.outer) {
      m = std::max({m, std::abs(v.x), std::abs(v.y)});
      bx0 = std::min(bx0, v.x), bx1 = std::max(bx1, v.x), by0 = std::min(by0, v.y), by1 = std::max(by1, v.y);
    }
    const double w = std::max(bx1 - bx0, by1 - by0);
    bx0 -= w, by0 -= w, bx1 += w, by1 += w;
    m = std::max({m, std::abs(bx0), std::abs(bx1), std::abs(by0), std::abs(by1)});
    scale = std::ldexp(1.0, bits - std::ilogb(m) - 1);  // |coordinates| < 2^bits after scaling
    Piece p{path(s.outer, true)};
    for (const auto& h : s.holes) p.push_back(path(h, false));
    pieces = {p};
    boxes = {Clipper2Lib::GetBounds(p)};
  }
  /// The input as the exact rationals that Clipper2 works on.
  pd::PolygonWithHoles<mpq_class> exact_piece(const Piece& p) const
  {
    const mpq_class inv = 1 / mpq_class(scale);
    auto ring = [&](const Clipper2Lib::Path64& r) {
      pd::Ring<mpq_class> out;
      for (const auto& v : r) out.push_back({mpq_class(static_cast<long>(v.x)) * inv, mpq_class(static_cast<long>(v.y)) * inv});
      return out;
    };
    pd::PolygonWithHoles<mpq_class> q{ring(p[0]), {}};
    for (std::size_t i = 1; i < p.size(); ++i) q.holes.push_back(ring(p[i]));
    return q;
  }
  static void collect(const Clipper2Lib::PolyPath64& node, std::vector<Piece>& out)
  {
    for (std::size_t i = 0; i < node.Count(); ++i) {
      const Clipper2Lib::PolyPath64* outer = node.Child(i);
      Piece piece{outer->Polygon()};
      for (std::size_t j = 0; j < outer->Count(); ++j) {
        piece.push_back(outer->Child(j)->Polygon());
        collect(*outer->Child(j), out);  // islands inside holes
      }
      out.push_back(std::move(piece));
    }
  }
  /// Cuts every piece that the line through a and b crosses; returns the number of pieces cut.
  std::size_t cut(const pd::Point<double>& a, const pd::Point<double>& b)
  {
    // H: the rectangle's part left of the line (Sutherland and Hodgman), its two corners on the line
    // rounded once.
    const double dx = b.x - a.x, dy = b.y - a.y;
    auto side = [&](double x, double y) { return dx * (y - a.y) - dy * (x - a.x); };
    const double cx[4] = {bx0, bx1, bx1, bx0}, cy[4] = {by0, by0, by1, by1};
    Clipper2Lib::Path64 H;
    std::vector<P64> on;
    for (int i = 0; i < 4; ++i) {
      const int j = (i + 1) % 4;
      const double si = side(cx[i], cy[i]), sj = side(cx[j], cy[j]);
      if (si >= 0) H.push_back(pt(cx[i], cy[i]));
      if ((si > 0 && sj < 0) || (si < 0 && sj > 0)) {
        const double t = si / (si - sj);
        H.push_back(pt(cx[i] + t * (cx[j] - cx[i]), cy[i] + t * (cy[j] - cy[i])));
        on.push_back(H.back());
      }
    }
    if (on.size() != 2) return 0;  // the line misses the rectangle (cannot happen for these shapes)
    auto orient = [&](const P64& q) {
      const i128 c = static_cast<i128>(on[1].x - on[0].x) * (q.y - on[0].y) -
                         static_cast<i128>(on[1].y - on[0].y) * (q.x - on[0].x);
      return (c > 0) - (c < 0);
    };
    std::vector<Piece> next;
    std::vector<Clipper2Lib::Rect64> next_boxes;
    std::size_t cut_pieces = 0;
    for (std::size_t i = 0; i < pieces.size(); ++i) {
      const auto& r = boxes[i];
      const int s0 = orient(P64(r.left, r.top));
      if (s0 != 0 && orient(P64(r.right, r.top)) == s0 && orient(P64(r.right, r.bottom)) == s0 &&
          orient(P64(r.left, r.bottom)) == s0) {
        next.push_back(std::move(pieces[i]));
        next_boxes.push_back(r);
        continue;
      }
      std::vector<Piece> parts[2];
      for (int op = 0; op < 2; ++op) {
        Clipper2Lib::Clipper64 c;
        c.AddSubject(pieces[i]);
        c.AddClip({H});
        Clipper2Lib::PolyTree64 tree;
        c.Execute(op == 0 ? Clipper2Lib::ClipType::Intersection : Clipper2Lib::ClipType::Difference,
                  Clipper2Lib::FillRule::NonZero, tree);
        collect(tree, parts[op]);
      }
      if (parts[0].empty() || parts[1].empty()) {  // the line misses the piece or only touches it
        next.push_back(std::move(pieces[i]));
        next_boxes.push_back(r);
        continue;
      }
      ++cut_pieces;
      for (auto& part : parts)
        for (auto& p : part) {
          next_boxes.push_back(Clipper2Lib::GetBounds(p));
          next.push_back(std::move(p));
        }
    }
    pieces = std::move(next);
    boxes = std::move(next_boxes);
    return cut_pieces;
  }
};

#ifdef PD_HAVE_CGAL
/// Cutting with CGAL's Boolean set operations on the exact kernel (Epeck): every cut computes
/// piece ∩ H and piece − H for the half-plane H left of the line, a large rectangle clipped exactly.
struct CgalSlicer {
  using Pt = Epeck::Point_2;
  using Poly = CGAL::Polygon_2<Epeck>;
  using Pwh = CGAL::Polygon_with_holes_2<Epeck>;
  std::vector<Pwh> pieces;
  std::vector<CGAL::Bbox_2> boxes;
  Pt corner[4];

  using ET = std::decay_t<decltype(CGAL::exact(std::declval<Epeck::FT>()))>;  // mpq_class or Boost's mpq
  template <class E = ET>
  static Epeck::FT ft(const mpq_class& v)
  {
    if constexpr (std::is_same_v<E, mpq_class>) return Epeck::FT(v);
    else return Epeck::FT(E(v.get_mpq_t()));
  }
  template <class E = ET>
  static mpq_class q(const Epeck::FT& x)
  {
    const E e = CGAL::exact(x);
    if constexpr (std::is_same_v<E, mpq_class>) return e;
    else return mpq_class(e.backend().data());
  }
  explicit CgalSlicer(const pd::PolygonWithHoles<double>& s)
  {
    auto poly = [](const pd::Ring<double>& r, bool ccw) {
      Poly p;
      for (const auto& v : r) p.push_back(Pt(v.x, v.y));
      if (p.is_counterclockwise_oriented() != ccw) p.reverse_orientation();
      return p;
    };
    Pwh p(poly(s.outer, true));
    for (const auto& h : s.holes) p.add_hole(poly(h, false));
    pieces = {p};
    boxes = {p.outer_boundary().bbox()};
    const CGAL::Bbox_2 b = boxes[0];
    const double w = std::max(b.xmax() - b.xmin(), b.ymax() - b.ymin());
    const double x0 = b.xmin() - w, x1 = b.xmax() + w, y0 = b.ymin() - w, y1 = b.ymax() + w;
    corner[0] = Pt(x0, y0), corner[1] = Pt(x1, y0), corner[2] = Pt(x1, y1), corner[3] = Pt(x0, y1);
  }
  std::size_t cut(const Q& a, const Q& b)
  {
    const Epeck::Line_2 l(Pt(ft(a.x), ft(a.y)), Pt(ft(b.x), ft(b.y)));
    auto value = [&](const Pt& p) { return l.a() * p.x() + l.b() * p.y() + l.c(); };
    Poly H;
    for (int i = 0; i < 4; ++i) {
      const Pt &p = corner[i], &r = corner[(i + 1) % 4];
      const CGAL::Oriented_side si = l.oriented_side(p), sj = l.oriented_side(r);
      if (si != CGAL::ON_NEGATIVE_SIDE) H.push_back(p);
      if (si != CGAL::ON_ORIENTED_BOUNDARY && sj != CGAL::ON_ORIENTED_BOUNDARY && si != sj) {
        const Epeck::FT vp = value(p), t = vp / (vp - value(r));  // exact crossing of the edge p r
        H.push_back(Pt(p.x() + t * (r.x() - p.x()), p.y() + t * (r.y() - p.y())));
      }
    }
    std::vector<Pwh> next;
    std::vector<CGAL::Bbox_2> next_boxes;
    std::size_t cut_pieces = 0;
    for (std::size_t i = 0; i < pieces.size(); ++i) {
      const CGAL::Bbox_2& r = boxes[i];
      const CGAL::Oriented_side s0 = l.oriented_side(Pt(r.xmin(), r.ymin()));
      if (s0 != CGAL::ON_ORIENTED_BOUNDARY && l.oriented_side(Pt(r.xmax(), r.ymin())) == s0 &&
          l.oriented_side(Pt(r.xmax(), r.ymax())) == s0 && l.oriented_side(Pt(r.xmin(), r.ymax())) == s0) {
        next.push_back(std::move(pieces[i]));
        next_boxes.push_back(r);
        continue;
      }
      std::vector<Pwh> in, out;
      // Segment traits (Tag_false): with the default polyline traits, a polygon whose hole touches
      // its outer ring at a vertex, as these operations return after a cut through a hole vertex,
      // violates the input precondition of the next operation.
      CGAL::intersection(pieces[i], H, std::back_inserter(in), CGAL::Tag_false());
      CGAL::difference(pieces[i], H, std::back_inserter(out), CGAL::Tag_false());
      if (in.empty() || out.empty()) {
        next.push_back(std::move(pieces[i]));
        next_boxes.push_back(r);
        continue;
      }
      ++cut_pieces;
      for (auto* part : {&in, &out})
        for (auto& p : *part) {
          next_boxes.push_back(p.outer_boundary().bbox());
          next.push_back(std::move(p));
        }
    }
    pieces = std::move(next);
    boxes = std::move(next_boxes);
    return cut_pieces;
  }
  /// The pieces in exact rationals (evaluates CGAL's lazy numbers), in the form the slicer reports
  /// them: a hole that touches the outer ring at a vertex is merged into that ring (which then
  /// passes through the vertex twice), and an outer ring that passes through a point twice with
  /// material on both lobes (two regions that touch at a point) is split into two pieces. The point
  /// sets are unchanged.
  std::vector<pd::PolygonWithHoles<mpq_class>> exact_pieces() const
  {
    using R = pd::Ring<mpq_class>;
    auto ring = [](const Poly& p) {
      R r;
      for (const auto& v : p.vertices()) r.push_back({q(v.x()), q(v.y())});
      return r;
    };
    std::vector<pd::PolygonWithHoles<mpq_class>> out;
    for (const auto& p : pieces) {
      R outer = ring(p.outer_boundary());
      std::vector<R> holes;
      for (const auto& h : p.holes()) {
        R hr = ring(h);
        std::size_t i = 0, j = 0;
        bool touch = false;
        for (i = 0; i < outer.size() && !touch; ++i)
          for (j = 0; j < hr.size() && !touch; ++j) touch = outer[i] == hr[j];
        if (!touch) {
          holes.push_back(std::move(hr));
          continue;
        }
        --i, --j;  // outer[i] == hr[j]: walk round the hole from there (both keep the material on the left)
        R merged(outer.begin(), outer.begin() + static_cast<std::ptrdiff_t>(i) + 1);
        for (std::size_t k = 1; k <= hr.size(); ++k) merged.push_back(hr[(j + k) % hr.size()]);
        merged.insert(merged.end(), outer.begin() + static_cast<std::ptrdiff_t>(i) + 1, outer.end());
        outer = std::move(merged);
      }
      // Split at points visited twice when both lobes are counterclockwise (material on both).
      std::vector<R> todo{outer}, outers;
      while (!todo.empty()) {
        R r = std::move(todo.back());
        todo.pop_back();
        bool split = false;
        for (std::size_t a = 0; a < r.size() && !split; ++a)
          for (std::size_t b = a + 1; b < r.size() && !split; ++b) {
            if (!(r[a] == r[b])) continue;
            R one(r.begin() + static_cast<std::ptrdiff_t>(a), r.begin() + static_cast<std::ptrdiff_t>(b));
            R two(r.begin() + static_cast<std::ptrdiff_t>(b), r.end());
            two.insert(two.end(), r.begin(), r.begin() + static_cast<std::ptrdiff_t>(a));
            if (sgn(pd::signed_area2(one)) > 0 && sgn(pd::signed_area2(two)) > 0) {
              todo.push_back(std::move(one)), todo.push_back(std::move(two));
              split = true;
            }
          }
        if (!split) outers.push_back(std::move(r));
      }
      for (std::size_t k = 0; k < outers.size(); ++k) {
        pd::PolygonWithHoles<mpq_class> e{outers[k], {}};
        for (const R& h : holes)  // a hole goes to the lobe that contains it
          if (outers.size() == 1 || pd::locate<pd::ExactKernel>(h.front(), outers[k]) == pd::Location::Inside)
            e.holes.push_back(h);
        out.push_back(std::move(e));
      }
    }
    return out;
  }
};
#endif

[[maybe_unused]] std::size_t max_bits_of(const std::vector<pd::PolygonWithHoles<mpq_class>>& ps)
{
  std::size_t b = 0;
  for (const auto& p : ps) {
    for (const auto& v : p.outer) b = std::max({b, pd::detail::bit_size(v.x), pd::detail::bit_size(v.y)});
    for (const auto& h : p.holes)
      for (const auto& v : h) b = std::max({b, pd::detail::bit_size(v.x), pd::detail::bit_size(v.y)});
  }
  return b;
}

/// Clipper2 at several integer scales (2^bits bounds the coordinates) on cheese with 16 and 64
/// holes, 60 cuts of every sequence; compared with the exact pieces of the slicer.
int run_clipper_scale_suite(std::ostream& csv)
{
  csv << "shape,sequence,cuts,bits,exact_pieces,pieces,tiny,area_err\n";
  for (int h : {16, 64}) {
    const std::string name = "cheese-" + std::to_string(h);
    const auto shape = pd::gen::cheese(h, 1);
    for (const std::string seq : {"random", "axis", "vertex", "repeat", "graze"}) {
      const auto lines = slice_lines(name, seq, shape, 60);
      pd::Slicer ours(shape);
      for (const auto& l : lines) ours.cut_exact(l.a, l.b);
      for (int bits : {20, 26, 32, 40, 50}) {
        ClipperSlicer c(shape, bits);
        const i128 in2 = ClipperSlicer::area2(c.pieces[0]);
        for (const auto& l : lines) c.cut({l.a.x.get_d(), l.a.y.get_d()}, {l.b.x.get_d(), l.b.y.get_d()});
        i128 total = 0;
        std::size_t tiny = 0;
        for (const auto& p : c.pieces) {
          const i128 a2 = ClipperSlicer::area2(p);
          total += a2;
          tiny += static_cast<double>(a2) < kTiny * static_cast<double>(in2);
        }
        const i128 d = total - in2;
        csv << name << ',' << seq << ',' << lines.size() << ',' << bits << ',' << ours.size() << ',' << c.pieces.size()
            << ',' << tiny << ',' << static_cast<double>(d < 0 ? -d : d) / static_cast<double>(in2) << '\n';
      }
    }
  }
  return 0;
}

int run_slice_suite(std::ostream& csv, int cuts, bool quick, double baseline_budget_s)
{
  (void)baseline_budget_s;  // used by the CGAL baseline
  std::vector<std::pair<std::string, pd::PolygonWithHoles<double>>> shapes;
  for (int h : quick ? std::vector<int>{16} : std::vector<int>{16, 64, 256})
    shapes.push_back({"cheese-" + std::to_string(h), pd::gen::cheese(h, 1)});
  if (!quick) {
    shapes.push_back({"holes-64", pd::gen::random_holes(64, 1)});
    for (const Input& in : hole_inputs(true))
      if (in.name == "legacy-bundled") shapes.push_back({"bundled", in.poly});
  }
  csv << "shape,n,holes,method,sequence,cut,pieces,cut_pieces,vertex_contacts,edge_contacts,dangles,symbolic,"
         "cut_ms,convex_ms,convex_pieces,max_bits,certified,tiny,area_err,split_pieces\n";
  int failures = 0;
  using clk = std::chrono::steady_clock;
  auto ms = [](clk::time_point t0, clk::time_point t1) { return std::chrono::duration<double, std::milli>(t1 - t0).count(); };
  for (const auto& [name, shape] : shapes) {
    std::size_t n = shape.outer.size();
    for (const auto& h : shape.holes) n += h.size();
    const bool cheese = name.rfind("cheese", 0) == 0;
    std::vector<std::string> seqs{"random", "axis", "vertex", "repeat"};
    if (cheese) seqs.push_back("graze");
    const auto exact_in = pd::convert_polygon<pd::ExactKernel>(shape);
    const double shape_area = area_of(shape);
    for (const std::string& seq : seqs) {
      const std::vector<SliceLine> lines = slice_lines(name, seq, shape, cuts);
      SliceRow row;
      row.shape = name, row.seq = seq, row.n = n, row.holes = shape.holes.size();

      // ours: pd::Slicer, certificate after every cut
      row.method = "ours";
      pd::Slicer s(shape);
      if (!s.ok()) {
        std::cerr << name << ": invalid shape\n";
        return 1;
      }
      for (int k = 0; k < cuts; ++k) {
        const auto t0 = clk::now();
        const pd::SliceResult r = s.cut_exact(lines[k].a, lines[k].b);
        const auto t1 = clk::now();
        if (!r.ok()) {
          std::cerr << name << "/" << seq << " cut " << k << ": " << r.diag.message << '\n';
          ++failures;
          break;
        }
        std::size_t convex_pieces = 0;
        for (std::size_t i : r.pieces) convex_pieces += s.convex(i).leaves.size();
        const auto t2 = clk::now();
        const bool cert = pd::verify_slices<pd::ExactKernel>(exact_in, s.pieces_exact(), s.symbolic()).ok;
        if (!cert) ++failures;
        row.cut = k, row.pieces = s.size(), row.cut_pieces = r.cut_pieces;
        row.vertex_contacts = std::to_string(r.vertex_contacts), row.edge_contacts = std::to_string(r.edge_contacts);
        row.dangles = std::to_string(r.dangles), row.symbolic = s.symbolic() ? "1" : "0";
        row.cut_ms = ms(t0, t1), row.convex_ms = std::to_string(ms(t1, t2));
        row.convex_pieces = std::to_string(convex_pieces), row.max_bits = std::to_string(s.max_bits());
        row.certified = cert ? "1" : "0";
        if (k == cuts - 1) {  // slivers, at the end
          std::size_t tiny = 0;
          for (std::size_t i = 0; i < s.size(); ++i) tiny += area_of(s.piece(i)) < kTiny * shape_area;
          row.tiny = std::to_string(tiny);
        }
        write_row(csv, row);
      }
      row = SliceRow();
      row.shape = name, row.seq = seq, row.n = n, row.holes = shape.holes.size();

      // Clipper2 on 64-bit integers; certificate against its own rounded input after every cut
      {
        row.method = "clipper2";
        ClipperSlicer c(shape);
        pd::PolygonWithHoles<mpq_class> rounded_in = c.exact_piece(c.pieces[0]);
        const i128 in_area2 = ClipperSlicer::area2(c.pieces[0]);
        for (int k = 0; k < cuts; ++k) {
          const pd::Point<double> a{lines[k].a.x.get_d(), lines[k].a.y.get_d()};
          const pd::Point<double> b{lines[k].b.x.get_d(), lines[k].b.y.get_d()};
          const auto t0 = clk::now();
          const std::size_t cut_pieces = c.cut(a, b);
          const auto t1 = clk::now();
          i128 total = 0;
          std::size_t tiny = 0;
          for (const auto& p : c.pieces) {
            const i128 a2 = ClipperSlicer::area2(p);
            total += a2;
            tiny += static_cast<double>(a2) < kTiny * static_cast<double>(in_area2);
          }
          row.cut = k, row.pieces = c.pieces.size(), row.cut_pieces = cut_pieces, row.cut_ms = ms(t0, t1);
          row.tiny = std::to_string(tiny), row.certified.clear();
          if (k == 0 || k == cuts - 1) {  // the certificate after the first and the last cut
            std::vector<pd::PolygonWithHoles<mpq_class>> ex;
            for (const auto& p : c.pieces) ex.push_back(c.exact_piece(p));
            row.certified = pd::verify_slices<pd::ExactKernel>(rounded_in, ex, true).ok ? "1" : "0";
          }
          const i128 d = total - in_area2;
          std::ostringstream err;
          err << static_cast<double>(d < 0 ? -d : d) / static_cast<double>(in_area2);
          row.area_err = err.str();
          write_row(csv, row);
        }
      }
#ifdef PD_HAVE_CGAL
      // CGAL Boolean set operations (Epeck); stops when its cut times exceed the budget; certificate
      // and coordinate sizes after the last cut (evaluating the lazy numbers earlier would change
      // its timings)
      {
        row = SliceRow();
        row.shape = name, row.method = "cgal-boolean", row.seq = seq, row.n = n, row.holes = shape.holes.size();
        CgalSlicer c(shape);
        double total_ms = 0;
        for (int k = 0; k < cuts; ++k) {
          const auto t0 = clk::now();
          const std::size_t cut_pieces = c.cut(lines[k].a, lines[k].b);
          const auto t1 = clk::now();
          total_ms += ms(t0, t1);
          row.cut = k, row.pieces = c.pieces.size(), row.cut_pieces = cut_pieces, row.cut_ms = ms(t0, t1);
          const bool last = k == cuts - 1 || total_ms > 1000 * baseline_budget_s;
          if (last) {
            const auto ex = c.exact_pieces();
            row.split_pieces = std::to_string(ex.size());
            const auto rep = pd::verify_slices<pd::ExactKernel>(exact_in, ex, true);
            const bool cert = rep.ok;
            if (!cert) std::cerr << '\n' << name << '/' << seq << " cgal-boolean cut " << k << ": " << rep.failure << '\n';
            row.certified = cert ? "1" : "0", row.max_bits = std::to_string(max_bits_of(ex));
            std::size_t tiny = 0;
            for (const auto& p : c.pieces) {
              double a2 = CGAL::to_double(p.outer_boundary().area());
              for (const auto& h : p.holes()) a2 += CGAL::to_double(h.area());  // holes are clockwise
              tiny += a2 < kTiny * shape_area;
            }
            row.tiny = std::to_string(tiny);
          }
          write_row(csv, row);
          if (last) break;
        }
      }
#endif
      csv.flush();
      std::cerr << '.' << std::flush;
    }
  }
  std::cerr << '\n';
  return failures == 0 ? 0 : 1;
}

} // namespace

int main(int argc, char** argv)
{
  std::string suite = "all", csv_path;
  int repeats = 5, cuts = 200;
  double budget_ms = 5000, baseline_budget_s = 120;
  bool quick = false;
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    if (a == "--suite" && i + 1 < argc) suite = argv[++i];
    else if (a == "--repeats" && i + 1 < argc) repeats = std::stoi(argv[++i]);
    else if (a == "--csv" && i + 1 < argc) csv_path = argv[++i];
    else if (a == "--quick") quick = true;
    else if (a == "--budget-ms" && i + 1 < argc) budget_ms = std::stod(argv[++i]);
    else if (a == "--cuts" && i + 1 < argc) cuts = std::stoi(argv[++i]);
    else if (a == "--baseline-budget-s" && i + 1 < argc) baseline_budget_s = std::stod(argv[++i]);
    else return usage();
  }
  if (suite != "simple" && suite != "holes" && suite != "all" && suite != "scale" && suite != "slice" &&
      suite != "clipper-scale")
    return usage();
  if (suite == "clipper-scale") {
    std::ofstream file;
    if (!csv_path.empty()) file.open(csv_path);
    return run_clipper_scale_suite(csv_path.empty() ? std::cout : file);
  }
  if (suite == "slice") {
    std::ofstream file;
    if (!csv_path.empty()) file.open(csv_path);
    return run_slice_suite(csv_path.empty() ? std::cout : file, cuts, quick, baseline_budget_s);
  }

  std::vector<Input> inputs;
  if (suite == "scale") inputs = scale_inputs();
  if (suite == "simple" || suite == "all") inputs = simple_inputs(quick);
  if (suite == "holes" || suite == "all") {
    auto h = hole_inputs(quick);
    inputs.insert(inputs.end(), h.begin(), h.end());
  }
  std::set<std::pair<std::string, std::string>> over_budget;  // (series, method)

  std::ofstream file;
  if (!csv_path.empty()) file.open(csv_path);
  std::ostream& csv = csv_path.empty() ? std::cout : file;
  csv << "suite,input,n,holes,reflex,lower_bound,method,status,pieces,median_ms,certified\n";

  for (const Input& in : inputs) {
    std::size_t n = in.poly.outer.size();
    for (const auto& h : in.poly.holes) n += h.size();
    const std::size_t reflex = reflex_count(in.poly);
    const bool simple = in.poly.holes.empty();
    const std::size_t lower = simple ? (reflex + 1) / 2 + 1 : 0;
    std::vector<Method> methods = our_methods();
    if (!simple && in.series.empty())
      methods.push_back({"ours-holes-exact", [](const pd::PolygonWithHoles<double>& p, bool) {
                           pd::HoleElimOptions o;
                           o.validate_input = false;
                           o.keep_tree = false;
                           const auto r = pd::eliminate_holes<pd::ExactKernel>(
                               pd::convert_polygon<pd::ExactKernel>(p), o);
                           return Result{r.ok(), r.leaves.size()};
                         }});
#ifdef PD_HAVE_CGAL
    for (auto& m : cgal_methods(simple, n)) methods.push_back(std::move(m));
#endif
    if (!in.series.empty()) {
      methods.push_back({"validate-double", [](const pd::PolygonWithHoles<double>& p, bool) {
                           return Result{pd::validate<pd::DoubleKernel>(p).ok(), 0};
                         }});
      methods.push_back({"validate-exact", [](const pd::PolygonWithHoles<double>& p, bool) {
                           return Result{pd::validate<pd::ExactKernel>(pd::convert_polygon<pd::ExactKernel>(p)).ok(), 0};
                         }});
    }
    for (const Method& m : methods) {
      if (over_budget.count({in.series, m.first})) continue;
      Result res = m.second(in.poly, true);  // untimed, checked by the certificate
      const int certified = res.certified;
      const double ms = median_ms(repeats, [&]() { res = m.second(in.poly, false); });
      if (!in.series.empty() && ms > budget_ms) over_budget.insert({in.series, m.first});
      csv << in.suite << ',' << in.name << ',' << n << ',' << in.poly.holes.size() << ',' << reflex << ','
          << lower << ',' << m.first << ',' << (res.ok ? "ok" : "failed") << ',' << res.pieces << ',' << ms
          << ',' << (certified < 0 ? "" : certified ? "1" : "0") << '\n';
      csv.flush();
    }
    if (!csv_path.empty()) std::cerr << "." << std::flush;
  }
  if (!csv_path.empty()) std::cerr << '\n';
  return 0;
}
