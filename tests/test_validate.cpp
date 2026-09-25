// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#include "doctest.h"

#include <algorithm>
#include <cmath>
#include <random>

#include "pd/kernel/double.hpp"
#include "pd/validate.hpp"
#ifdef PD_HAVE_GMP
#include "pd/kernel/exact.hpp"
#define PD_TEST_KERNELS pd::DoubleKernel, pd::ExactKernel
#else
#define PD_TEST_KERNELS pd::DoubleKernel
#endif

using pd::IssueKind;

namespace {

pd::Ring<double> R(std::initializer_list<std::pair<double, double>> pts)
{
  pd::Ring<double> r;
  for (const auto& p : pts) r.push_back({p.first, p.second});
  return r;
}

template <class K>
pd::ValidationReport check(const pd::PolygonWithHoles<double>& p)
{
  return pd::validate<K>(pd::convert_polygon<K>(p));
}

bool has(const pd::ValidationReport& rep, IssueKind k)
{
  return std::any_of(rep.issues.begin(), rep.issues.end(),
                     [&](const pd::Issue& i) { return i.kind == k; });
}

const pd::Ring<double> kSquare = R({{0, 0}, {10, 0}, {10, 10}, {0, 10}});

} // namespace

TEST_CASE_TEMPLATE("valid inputs", K, PD_TEST_KERNELS)
{
  CHECK(check<K>({kSquare, {}}).ok());
  CHECK(check<K>({kSquare, {R({{2, 2}, {4, 2}, {4, 4}, {2, 4}})}}).ok());
  // Entangled holes: a C-shaped hole with a second hole in its pocket.
  const auto c_hole = R({{1, 1}, {9, 1}, {9, 3}, {3, 3}, {3, 7}, {9, 7}, {9, 9}, {1, 9}});
  const auto pocket = R({{5, 4}, {8, 4}, {8, 6}, {5, 6}});
  CHECK(check<K>({kSquare, {c_hole, pocket}}).ok());
  // Orientation is not part of validity.
  pd::Ring<double> cw(kSquare.rbegin(), kSquare.rend());
  CHECK(check<K>({cw, {}}).ok());
}

TEST_CASE_TEMPLATE("invalid inputs are reported", K, PD_TEST_KERNELS)
{
  CHECK(has(check<K>({R({{0, 0}, {1, 0}}), {}}), IssueKind::TooFewVertices));
  CHECK(has(check<K>({R({{0, 0}, {10, 0}, {10, 0}, {0, 10}}), {}}), IssueKind::RepeatedVertex));
  // Bow tie.
  CHECK(has(check<K>({R({{0, 0}, {10, 10}, {10, 0}, {0, 10}}), {}}), IssueKind::EdgesIntersect));
  // Ring touching itself at a vertex (figure eight through (5, 5)).
  CHECK(has(check<K>({R({{0, 0}, {5, 5}, {10, 0}, {10, 10}, {5, 5}, {0, 10}}), {}}),
            IssueKind::EdgesIntersect));
  // Adjacent edges folding back along a common line (a spike).
  const auto spike = check<K>({R({{0, 0}, {4, 0}, {4, 4}, {4, 5}, {4, 4.5}, {0, 4}}), {}});
  CHECK(std::any_of(spike.issues.begin(), spike.issues.end(), [](const pd::Issue& i) {
    return i.kind == IssueKind::EdgesIntersect && i.ring_a == 0 && i.edge_a == 2 &&
           i.ring_b == 0 && i.edge_b == 3;
  }));
  // Hole touching the outer ring at one point.
  CHECK(has(check<K>({kSquare, {R({{0, 5}, {3, 4}, {3, 6}})}}), IssueKind::EdgesIntersect));
  // Two holes touching at a vertex (excluded: holes must be disjoint).
  CHECK(has(check<K>({kSquare, {R({{2, 2}, {5, 2}, {5, 5}}), R({{5, 5}, {8, 5}, {8, 8}})}}),
            IssueKind::EdgesIntersect));
  // Hole outside the outer ring.
  CHECK(has(check<K>({kSquare, {R({{12, 2}, {14, 2}, {14, 4}})}}), IssueKind::HoleNotInsideOuter));
  // Hole inside another hole.
  const auto rep = check<K>({kSquare, {R({{1, 1}, {9, 1}, {9, 9}, {1, 9}}),
                                       R({{4, 4}, {6, 4}, {6, 6}, {4, 6}})}});
  REQUIRE(has(rep, IssueKind::HoleInsideHole));
  CHECK(pd::describe(rep.issues.front()) == "hole 1 lies inside hole 0");
}

namespace {

// Random rings on a small integer grid: many collinear and touching edges.
pd::PolygonWithHoles<double> random_grid_polygon(std::mt19937_64& rng)
{
  std::uniform_int_distribution<int> coord(0, 12), nverts(3, 12), nholes(0, 3);
  auto ring = [&]() {
    pd::Ring<double> r(static_cast<std::size_t>(nverts(rng)));
    for (auto& p : r) p = {double(coord(rng)), double(coord(rng))};
    return r;
  };
  pd::PolygonWithHoles<double> p{ring(), {}};
  for (int h = nholes(rng); h > 0; --h) p.holes.push_back(ring());
  return p;
}

// Random star-shaped ring (valid) with non-integer coordinates.
pd::Ring<double> random_star(std::mt19937_64& rng, double cx, double cy, double rmax, int n)
{
  const double two_pi = 6.283185307179586;
  std::uniform_real_distribution<double> ang(0, two_pi), rad(0.3 * rmax, rmax);
  std::vector<double> a(static_cast<std::size_t>(n));
  for (double& v : a) v = ang(rng);
  std::sort(a.begin(), a.end());
  pd::Ring<double> r;
  for (double v : a) {
    const double rr = rad(rng);
    r.push_back({cx + rr * std::cos(v), cy + rr * std::sin(v)});
  }
  return r;
}

} // namespace

TEST_CASE_TEMPLATE("sweep and naive edge-pair search agree", K, PD_TEST_KERNELS)
{
  std::mt19937_64 rng(7);
  int invalid = 0;
  for (int i = 0; i < 300; ++i) {
    const auto p = pd::convert_polygon<K>(random_grid_polygon(rng));
    const auto sweep = pd::conflicting_edge_pairs<K>(p, true);
    const auto naive = pd::conflicting_edge_pairs<K>(p, false);
    REQUIRE(sweep.size() == naive.size());
    CHECK(sweep == naive);
    invalid += !sweep.empty();
  }
  for (int i = 0; i < 100; ++i) {
    pd::PolygonWithHoles<double> d{random_star(rng, 0, 0, 10, 40), {}};
    d.holes.push_back(random_star(rng, 0, 0, 2, 12));
    const auto p = pd::convert_polygon<K>(d);
    CHECK(pd::conflicting_edge_pairs<K>(p, true) == pd::conflicting_edge_pairs<K>(p, false));
  }
  MESSAGE("random grid polygons with conflicts: " << invalid << " of 300");
  CHECK(invalid > 0);
}

#ifdef PD_HAVE_GMP
TEST_CASE("double and exact kernels agree on grid inputs")
{
  // Integer coordinates make every double predicate exact, so the two kernels
  // must report identical conflicts.
  std::mt19937_64 rng(11);
  for (int i = 0; i < 300; ++i) {
    const auto d = random_grid_polygon(rng);
    CHECK(pd::conflicting_edge_pairs<pd::DoubleKernel>(d) ==
          pd::conflicting_edge_pairs<pd::ExactKernel>(pd::convert_polygon<pd::ExactKernel>(d)));
  }
}
#endif

TEST_CASE("coordinate range of the double kernel")
{
  const pd::Ring<double> big{{0, 0}, {1e130, 0}, {0, 1e130}};
  const pd::Ring<double> ok{{0, 0}, {1e100, 0}, {0, 1e100}};
  const pd::Ring<double> tiny{{0, 0}, {1e-130, 0}, {0, 1e-130}};
  const auto r = pd::validate<pd::DoubleKernel>({big, {}});
  REQUIRE_FALSE(r.ok());
  CHECK(r.issues.front().kind == IssueKind::CoordinateOutOfRange);
  CHECK(pd::validate<pd::DoubleKernel>({ok, {}}).ok());
  CHECK_FALSE(pd::validate<pd::DoubleKernel>({tiny, {}}).ok());
  CHECK_FALSE(pd::validate<pd::DoubleKernel>({{{0, 0}, {INFINITY, 0}, {0, 1}}, {}}).ok());
#ifdef PD_HAVE_GMP
  CHECK(pd::validate<pd::ExactKernel>(pd::convert_polygon<pd::ExactKernel>({big, {}})).ok());
#endif
}
