// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#include "doctest.h"

#include <cmath>
#include <random>
#include <string>

#include "pd/generate.hpp"
#include "pd/hole_elim.hpp"
#include "pd/io/obj.hpp"
#include "pd/kernel/double.hpp"
#ifdef PD_HAVE_GMP
#include "pd/kernel/exact.hpp"
#include "pd/verify.hpp"
#endif

namespace {

const std::string kData = PD_DATA_DIR;
std::string legacy(const std::string& rel) { return kData + "/legacy/" + rel; }

pd::Ring<double> R(std::initializer_list<std::pair<double, double>> pts)
{
  pd::Ring<double> r;
  for (const auto& p : pts) r.push_back({p.first, p.second});
  return r;
}

[[maybe_unused]] pd::Ring<double> square(double x0, double y0, double x1, double y1)
{
  return R({{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}});
}

// The bundled legacy case and the test3..test10 directories.
pd::PolygonWithHoles<double> bundled()
{
  pd::PolygonWithHoles<double> p;
  p.outer = pd::read_obj_ring(legacy("compgon/cpolygon0.obj"));
  for (int i = 0; i <= 6; ++i)
    p.holes.push_back(pd::read_obj_ring(legacy("compgon/hole0-" + std::to_string(i) + ".obj")));
  return p;
}

pd::PolygonWithHoles<double> legacy_dir(const std::string& dir, const std::string& outer,
                                        std::initializer_list<const char*> holes)
{
  pd::PolygonWithHoles<double> p;
  p.outer = pd::read_obj_ring(legacy("compgon/" + dir + "/" + outer));
  for (const char* h : holes) p.holes.push_back(pd::read_obj_ring(legacy("compgon/" + dir + "/" + h)));
  return p;
}

std::vector<std::pair<std::string, pd::PolygonWithHoles<double>>> legacy_cases()
{
  std::vector<std::pair<std::string, pd::PolygonWithHoles<double>>> v;
  v.emplace_back("bundled", bundled());
  v.emplace_back("test3", legacy_dir("test3", "p3.obj", {"h3.obj"}));
  v.emplace_back("test4", legacy_dir("test4", "p.obj", {"h1.obj", "h2.obj"}));
  v.emplace_back("test5", legacy_dir("test5", "p5.obj", {"h5.obj"}));
  v.emplace_back("test6", legacy_dir("test6", "p.obj", {"h1.obj", "h2.obj"}));
  v.emplace_back("test7", legacy_dir("test7", "p.obj", {"h1.obj", "h2.obj"}));
  v.emplace_back("test8", legacy_dir("test8", "p.obj", {"h.obj"}));
  v.emplace_back("test9", legacy_dir("test9", "p.obj", {"h.obj"}));
  v.emplace_back("test10", legacy_dir("test10", "p.obj",
                                      {"h1.obj", "h2.obj", "h3.obj", "h4.obj", "h5.obj", "h6.obj", "h7.obj"}));
  // ctest1.obj: three rings in one file; the largest is the outer ring.
  auto rings = pd::read_obj_rings(legacy("simp/ctest1.obj"));
  std::size_t big = 0;
  for (std::size_t i = 1; i < rings.size(); ++i)
    if (std::abs(pd::signed_area2(rings[i])) > std::abs(pd::signed_area2(rings[big]))) big = i;
  pd::PolygonWithHoles<double> c{rings[big], {}};
  for (std::size_t i = 0; i < rings.size(); ++i)
    if (i != big) c.holes.push_back(rings[i]);
  v.emplace_back("ctest1", c);
  return v;
}

} // namespace

#ifdef PD_HAVE_GMP

namespace {

using EK = pd::ExactKernel;

pd::HoleElimination<EK::FT> run_exact(const pd::PolygonWithHoles<double>& in)
{
  pd::HoleElimOptions opt;
  opt.check = true;
  return pd::eliminate_holes<EK>(pd::convert_polygon<EK>(in), opt);
}

// Runs the exact elimination and certifies the result (weakly simple pieces are
// allowed only after a dangling contact).
pd::HoleElimination<EK::FT> run_certified(const pd::PolygonWithHoles<double>& in)
{
  auto r = run_exact(in);
  REQUIRE_MESSAGE(r.ok(), r.diag.message);
  const auto rep =
      pd::verify_partition<EK>(pd::convert_polygon<EK>(in), r.leaf_rings(), false, r.stats.dangles > 0);
  CHECK_MESSAGE(rep.ok, rep.failure);
  CHECK(r.stats.cuts <= in.holes.size());
  return r;
}

} // namespace

TEST_CASE("no holes: the outer ring is the only piece")
{
  const auto r = run_certified({square(0, 0, 10, 10), {}});
  CHECK(r.leaves.size() == 1);
  CHECK(r.stats.cuts == 0);
}

TEST_CASE("one square hole")
{
  const auto r = run_certified({square(0, 0, 10, 10), {square(3, 4, 5, 6)}});
  CHECK(r.stats.cuts == 1);
  CHECK(r.leaves.size() >= 2);
}

TEST_CASE("concentric square-in-square (segfaulted the old code)")
{
  // The preferred line runs along a diagonal through four corners, each with its
  // edges on both sides: it is used as it is.
  const auto r = run_certified({square(-2, -2, 2, 2), {square(-0.5, -0.5, 0.5, 0.5)}});
  CHECK(r.stats.cuts == 1);
  CHECK(r.stats.fallback_rotated == 0);
  CHECK(r.stats.vertex_contacts == 4);
  CHECK(r.leaves.size() == 2);
}

TEST_CASE("a hole whose centroid is one of its own vertices")
{
  // Arrowhead with vertex average (0, -1) = its reflex vertex.
  const auto dart = R({{0, 3}, {3, -3}, {0, -1}, {-3, -3}});
  const auto r = run_certified({square(-10, -10, 10, 10), {dart}});
  CHECK(r.stats.fallback_midpoint == 0);
  CHECK(r.stats.vertex_contacts >= 1);
}

TEST_CASE("a hole centroid on another hole's vertex (entangled holes)")
{
  // C-shaped hole with centroid (5.5, 5) in its pocket; a triangle hole has a
  // vertex exactly there, and the preferred line runs through it.
  const auto c_hole = R({{1, 1}, {9, 1}, {9, 3}, {3, 3}, {3, 7}, {9, 7}, {9, 9}, {1, 9}});
  const auto tri = R({{5.5, 5}, {8, 4}, {8, 6}});
  const auto r = run_certified({square(0, 0, 10, 10), {c_hole, tri}});
  CHECK(r.stats.fallback_midpoint == 0);
  CHECK(r.stats.vertex_contacts >= 1);
}

pd::PolygonWithHoles<double> degenerate(const std::string& name)
{
  for (const auto& d : pd::gen::degenerate_set())
    if (d.first == name) return d.second;
  FAIL("no degenerate case " << name);
  return {};
}

TEST_CASE("a cut line through vertices of every kind")
{
  const auto r = run_certified(degenerate("vertex-contacts"));
  CHECK(r.stats.cuts == 1);  // the line y = 0 crosses all four holes
  CHECK(r.stats.vertex_contacts == 7);
  CHECK(r.stats.edge_contacts == 0);
  CHECK(r.stats.fallback_rotated + r.stats.fallback_midpoint == 0);
  CHECK(r.leaves.size() == 4);
}

TEST_CASE("a cut line along edges")
{
  const auto r = run_certified(degenerate("edge-contacts"));
  CHECK(r.stats.cuts == 1);
  CHECK(r.stats.edge_contacts == 2);
  CHECK(r.stats.fallback_rotated + r.stats.fallback_midpoint == 0);
  CHECK(r.leaves.size() == 2);
}

TEST_CASE("a cut line along edges of the outer ring")
{
  // A wing (a step of the outer ring) and a wedge (the bottom of a slot) on the line y = 0.
  const auto r = run_certified(degenerate("ring-edge-contacts"));
  CHECK(r.stats.cuts == 1);
  CHECK(r.stats.edge_contacts == 2);
  CHECK(r.stats.vertex_contacts == 4);
  CHECK(r.stats.fallback_rotated + r.stats.fallback_midpoint == 0);
  CHECK(r.leaves.size() == 3);  // above: left and right of the slot; below: one piece
}

TEST_CASE("after a dangling contact no hole-free piece repeats a vertex (fuzz finding)")
{
  // Found by fuzz_grid in 1.2.0 (present since 1.1.0): a later piece runs along the edge from the
  // dangling vertex to its split vertex, whose length is zero at the limit.
  const auto in = pd::PolygonWithHoles<double>{
      R({{120, 72}, {80, 126}, {62, 75}, {8, 42}, {94, 16}}),
      {R({{41, 41}, {40, 41}, {47, 42}}), R({{51, 48}, {50, 55}, {47, 48}, {40, 49}, {47, 45}}),
       R({{43, 42}, {42, 42}, {45, 44}}), R({{48, 40}, {47, 41}, {46, 40}, {46, 39}, {47, 39}, {48, 39}})}};
  const auto r = run_certified(in);
  CHECK(r.stats.dangles >= 1);
  for (std::size_t i = 0; i < r.leaves.size(); ++i) {
    const auto ring = r.leaf_ring(i);
    for (std::size_t k = 0; k < ring.size(); ++k) CHECK(!(ring[k] == ring[(k + 1) % ring.size()]));
  }
}

TEST_CASE("dangling contacts are split symbolically")
{
  const auto in = degenerate("dangling-contacts");
  const auto r = run_certified(in);  // weakly simple
  CHECK(r.stats.symbolic);
  CHECK(r.stats.dangles == 2);
  CHECK(r.stats.cuts == 1);
  CHECK(r.stats.fallback_rotated + r.stats.fallback_midpoint == 0);
  // The pieces touch themselves at the dangling vertices: not strictly simple.
  CHECK_FALSE(pd::verify_partition<EK>(pd::convert_polygon<EK>(in), r.leaf_rings(), false).ok);
}

TEST_CASE("the double kernel falls back for dangling contacts")
{
  const auto r = pd::eliminate_holes<pd::DoubleKernel>(degenerate("dangling-contacts"));
  REQUIRE_MESSAGE(r.ok(), r.diag.message);
  CHECK(r.stats.dangles == 0);
  CHECK(r.stats.fallback_rotated + r.stats.fallback_midpoint >= 1);
}

TEST_CASE("two holes with the same centroid")
{
  const auto c_hole = R({{1, 1}, {9, 1}, {9, 3}, {3, 3}, {3, 7}, {9, 7}, {9, 9}, {1, 9}});
  const auto inner = square(5, 4.5, 6, 5.5);  // centred at (5.5, 5), the C's centroid
  const auto r = run_certified({square(0, 0, 10, 10), {c_hole, inner}});
  CHECK(r.stats.degenerate_directions >= 1);
}

TEST_CASE("one hole whose centroid equals the anchor")
{
  // Triangle outer: the anchor is the triangle's centroid (4, 4), also the hole's.
  const auto r = run_certified({R({{0, 0}, {12, 0}, {0, 12}}), {square(3, 3, 5, 5)}});
  CHECK(r.stats.degenerate_directions == 1);
}

TEST_CASE("legacy inputs are eliminated and certified")
{
  for (const auto& c : legacy_cases()) {
    CAPTURE(c.first);
    const auto r = run_certified(c.second);
    MESSAGE(c.first << ": holes " << c.second.holes.size() << ", cuts " << r.stats.cuts
                    << ", pieces " << r.leaves.size() << ", depth " << r.stats.max_depth
                    << ", max crossings " << r.stats.max_crossings << ", fallbacks "
                    << r.stats.fallback_rotated + r.stats.fallback_midpoint);
  }
}

TEST_CASE("bundled input: same tree shape as the legacy implementation")
{
  // Legacy run: 3 cuts, 16 hole-free pieces, depth 2.
  const auto r = run_certified(bundled());
  CHECK(r.stats.cuts == 3);
  CHECK(r.leaves.size() == 16);
  CHECK(r.stats.max_depth == 2);
}

TEST_CASE("random polygons with many holes are certified")
{
  // Star-shaped outer ring; square and triangular holes on a jittered grid,
  // kept only if the input is valid.
  std::mt19937_64 rng(2026);
  std::uniform_real_distribution<double> jit(-0.3, 0.3), size(0.2, 0.45);
  int tested = 0;
  for (int it = 0; it < 60; ++it) {
    pd::PolygonWithHoles<double> p;
    const int nv = 24;
    for (int k = 0; k < nv; ++k) {
      const double a = 6.283185307179586 * k / nv;
      const double r = 8 + 2 * std::sin(3 * a) + jit(rng);
      p.outer.push_back({r * std::cos(a), r * std::sin(a)});
    }
    for (int gx = -4; gx <= 4; gx += 2)
      for (int gy = -4; gy <= 4; gy += 2) {
        const double cx = gx + jit(rng), cy = gy + jit(rng), s = size(rng);
        if ((gx + gy) % 4 == 0) p.holes.push_back(square(cx - s, cy - s, cx + s, cy + s));
        else p.holes.push_back(R({{cx - s, cy - s}, {cx + s, cy - s * 0.5}, {cx, cy + s}}));
      }
    if (!pd::validate<EK>(pd::convert_polygon<EK>(p)).ok()) continue;
    ++tested;
    CAPTURE(it);
    run_certified(p);
  }
  CHECK(tested > 40);
}

#endif // PD_HAVE_GMP

TEST_CASE("double kernel runs the legacy inputs")
{
  for (const auto& c : legacy_cases()) {
    CAPTURE(c.first);
    const auto r = pd::eliminate_holes<pd::DoubleKernel>(c.second);
    CHECK_MESSAGE(r.ok(), r.diag.message);
#ifdef PD_HAVE_GMP
    const auto e = pd::eliminate_holes<pd::ExactKernel>(pd::convert_polygon<pd::ExactKernel>(c.second));
    CHECK(r.leaves.size() == e.leaves.size());
    CHECK(r.stats.cuts == e.stats.cuts);
#endif
  }
}

#ifdef PD_HAVE_GMP
TEST_CASE("eliminate_holes with the double kernel: self-check and exact retry")
{
  // Near-touching holes, seed 49: a rounded construction makes the piece areas miss the input area,
  // which the self-check reports (found by the floating-point census of hole elimination).
  const auto in = pd::gen::generate(pd::gen::Family::NearTouch, 49);
  pd::HoleElimOptions opt;
  const auto checked = pd::eliminate_holes<pd::DoubleKernel>(in, opt);
  CHECK(checked.diag.status == pd::Status::NumericalFailure);
  CHECK(checked.diag.message.find("self-check") != std::string::npos);

  opt.self_check = false;  // unchecked, the failure goes unnoticed
  CHECK(pd::eliminate_holes<pd::DoubleKernel>(in, opt).ok());

  opt.self_check = true;
  opt.on_failure = pd::OnFailure::RetryExact;
  const auto retried = pd::eliminate_holes<pd::DoubleKernel>(in, opt);
  REQUIRE(retried.ok());
  CHECK(retried.stats.retried_exact);
  const auto exact = pd::eliminate_holes<pd::ExactKernel>(pd::convert_polygon<pd::ExactKernel>(in));
  REQUIRE(exact.ok());
  CHECK(retried.leaves == exact.leaves);  // the exact result, rounded once

  // An ordinary input passes the self-check without a retry.
  const auto fine = pd::eliminate_holes<pd::DoubleKernel>(pd::gen::cheese(16, 1), opt);
  CHECK(fine.ok());
  CHECK_FALSE(fine.stats.retried_exact);
}
#endif

#ifdef PD_HAVE_GMP
// Regression tests for 1.4.3. Before it, a single hole anchor could average a split vertex
// v1 + eps d, so a later automatic cut line depended on eps. These inputs come from proof reviews.
namespace {

using Q = mpq_class;
using QRing = pd::Ring<Q>;
using QPoly = pd::PolygonWithHoles<Q>;

QRing qbox(long x0, long y0, long x1, long y1) { return {{Q(x0), Q(y0)}, {Q(x1), Q(y0)}, {Q(x1), Q(y1)}, {Q(x0), Q(y1)}}; }

// Every reported piece has positive area, and the pieces pass the (weakly simple) certificate.
void certify_symbolic(const QPoly& in, const pd::HoleElimination<Q>& r)
{
  REQUIRE_MESSAGE(r.ok(), r.diag.message);
  for (std::size_t i = 0; i < r.leaves.size(); ++i) {
    CAPTURE(i);
    CHECK(pd::signed_area2(r.leaf_ring(i)) > 0);
  }
  const auto rep = pd::verify_partition<pd::ExactKernel>(in, r.leaf_rings(), false, r.stats.symbolic);
  CHECK_MESSAGE(rep.ok, rep.failure);
}

} // namespace

TEST_CASE("an anchor next to a split vertex gives no collapsed piece")
{
  // The first cut, y = 0, splits the triangle hole H2 at (0, 0). The anchor of H3 is then next to
  // the split vertex, and an eps dependent line would cut off a sliver along the edge (40,16)(46,18)
  // that collapses to a two vertex piece at eps = 0.
  QPoly in;
  for (auto [x, y] : std::vector<std::pair<long, long>>{
           {-220, -20}, {220, -20}, {220, 12}, {-44, 12}, {-44, 14}, {40, 16}, {46, 18}, {-44, 20}, {-220, 20}})
    in.outer.push_back({Q(x), Q(y)});
  in.holes = {qbox(-110, -2, -66, 2), qbox(66, -2, 110, 2), {{Q(0), Q(0)}, {Q(22), Q(4)}, {Q(-22), Q(4)}},
              qbox(11, 9, 33, 11)};
  REQUIRE(pd::validate<pd::ExactKernel>(in).ok());
  const auto r = pd::eliminate_holes<pd::ExactKernel>(in);
  CHECK(r.stats.dangles == 1);
  CHECK(r.stats.symbolic);
  CHECK(r.stats.collapsed == 0);
  certify_symbolic(in, r);
}

TEST_CASE("an anchor next to a split vertex gives no corridor pinched along an edge")
{
  QPoly in;
  for (auto [x, y] : std::vector<std::pair<long, long>>{
           {-60, -80}, {60, -80}, {60, 80}, {30, 80}, {30, 0}, {10, 0}, {10, 80}, {-60, 80}})
    in.outer.push_back({Q(x), Q(y)});
  in.holes = {{{Q(-21), Q(-60)}, {Q(-19), Q(-61)}, {Q(-20), Q(-59)}},
              {{Q(-21), Q(60)}, {Q(-19), Q(61)}, {Q(-20), Q(59)}},
              {{Q(-20), Q(-30)}, {Q(-19), Q(10)}, {Q(-18), Q(20)}},
              {{Q(-11), Q(0)}, {Q(-9), Q(1)}, {Q(-10), Q(-1)}}};
  REQUIRE(pd::validate<pd::ExactKernel>(in).ok());
  const auto r = pd::eliminate_holes<pd::ExactKernel>(in);
  CHECK(r.stats.symbolic);
  certify_symbolic(in, r);
}

TEST_CASE("peeling: a later dangling split does not cross an earlier anchor cut")
{
  // With an eps dependent anchor cut, the retained hole H3 came within O(eps) of that cut, and a
  // later split of H3 crossed it ("the outermost crossings are not on the outer ring").
  const auto q = [](long n, long d = 1) { return Q(n, d); };
  QPoly in;
  in.outer = {{q(-100), q(-100)}, {q(100), q(-100)}, {q(100), q(100)}, {q(-78), q(100)},
              {q(-79), q(98)},    {q(-80), q(100)},  {q(-100), q(100)}};
  in.holes.push_back({{q(-9993, 100), q(-999, 10)}, {q(-9987, 100), q(-9993, 100)}, {q(-999, 10), q(-9987, 100)}});
  in.holes.push_back({{q(0), q(0)}, {q(-2), q(1)}, {q(-1), q(2)}});
  const pd::Point<Q> c2{q(-1), q(1, 5)};
  in.holes.push_back({{q(-101, 100), q(1, 5)}, {q(-99, 100), q(21, 100)}, {q(-1), q(19, 100)}});
  const pd::Point<Q> a2{q(-10189, 300), q(-9889, 300)};
  const pd::Point<Q> p3{c2.x + (c2.x - a2.x) / 2, c2.y + (c2.y - a2.y) / 2};
  in.holes.push_back({p3, {p3.x - q(1, 2), p3.y + q(1)}, {p3.x - q(1), p3.y + q(3)}});
  const pd::Point<Q> a3{q(-79), q(298, 3)};
  const pd::Point<Q> c4{a3.x + (p3.x - a3.x) / 20, a3.y + (p3.y - a3.y) / 20};
  in.holes.push_back({{c4.x - q(1, 10), c4.y}, {c4.x + q(1, 10), c4.y + q(1, 10)}, {c4.x, c4.y - q(1, 10)}});
  REQUIRE(pd::validate<pd::ExactKernel>(in).ok());
  pd::HoleElimOptions opt;
  opt.cut_rule = pd::CutRule::NearestToBoundary;
  const auto r = pd::eliminate_holes<pd::ExactKernel>(in, opt);
  CHECK(r.stats.dangles >= 1);
  certify_symbolic(in, r);
}

TEST_CASE("the symbolic limit omits zero area traces and renumbers the tree")
{
  using pd::EpsQ;
  const EpsQ e = EpsQ::eps();
  pd::HoleElimination<EpsQ> s;
  // 0..2: a triangle. 3..5: collinear at eps = 0 but not for eps > 0. 6, 7: a two vertex trace.
  s.pool = {{EpsQ(0), EpsQ(0)}, {EpsQ(4), EpsQ(0)}, {EpsQ(0), EpsQ(4)}, {EpsQ(1), EpsQ(1)},
            {EpsQ(3), EpsQ(3) + e},  {EpsQ(2), EpsQ(2)}, {EpsQ(5), EpsQ(5)}, {EpsQ(6), EpsQ(6)}};
  s.leaves = {{0, 1, 2}, {3, 4, 5}, {6, 7}, {1, 2, 0}};
  s.tree.resize(5);
  s.tree[0].children = {1, 2, 3, 4};
  for (int i = 1; i <= 4; ++i) s.tree[static_cast<std::size_t>(i)].parent = 0, s.tree[static_cast<std::size_t>(i)].leaf = i - 1;
  const auto r = pd::detail::symbolic_limit(std::move(s));
  REQUIRE(r.leaves.size() == 2);
  CHECK(r.stats.collapsed == 2);
  CHECK(r.tree[1].leaf == 0);
  CHECK(r.tree[2].leaf == -1);
  CHECK(r.tree[3].leaf == -1);
  CHECK(r.tree[4].leaf == 1);
  CHECK(r.tree[0].leaf == -1);
}
#endif
