// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
// pd_bench: piece counts and timings of polydecomp, of PolyPartition's
// Hertel-Mehlhorn partition, and of CGAL's convex decompositions when built
// with PD_WITH_CGAL.
//
//   pd_bench [--suite simple|holes|all|scale] [--repeats R] [--csv FILE] [--quick]
//            [--budget-ms T]
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
#include <set>
#include <string>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "pd/generate.hpp"
#include "pd/hole_elim.hpp"
#include "pd/io/obj.hpp"
#include "pd/kernel/double.hpp"
#include "pd/kernel/exact.hpp"
#include "pd/pipeline.hpp"
#include "pd/verify.hpp"

#include "polypartition.h"

#ifdef PD_HAVE_CGAL
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
  std::cerr << "usage: pd_bench [--suite simple|holes|all|scale] [--repeats R] [--csv FILE] [--quick]\n"
               "                [--budget-ms T]\n";
  return 2;
}

} // namespace

int main(int argc, char** argv)
{
  std::string suite = "all", csv_path;
  int repeats = 5;
  double budget_ms = 5000;
  bool quick = false;
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    if (a == "--suite" && i + 1 < argc) suite = argv[++i];
    else if (a == "--repeats" && i + 1 < argc) repeats = std::stoi(argv[++i]);
    else if (a == "--csv" && i + 1 < argc) csv_path = argv[++i];
    else if (a == "--quick") quick = true;
    else if (a == "--budget-ms" && i + 1 < argc) budget_ms = std::stod(argv[++i]);
    else return usage();
  }
  if (suite != "simple" && suite != "holes" && suite != "all" && suite != "scale") return usage();

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
