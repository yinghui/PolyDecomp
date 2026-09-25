// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#include "doctest.h"

#include <cmath>
#include <random>

#include "pd/convex_decomp.hpp"
#include "pd/generate.hpp"
#include "pd/kernel/double.hpp"
#include "pd/pipeline.hpp"
#include "pd/verify.hpp"
#ifdef PD_HAVE_GMP
#include "pd/kernel/exact.hpp"
#endif

#ifdef PD_HAVE_GMP

namespace {

using EK = pd::ExactKernel;
using DK = pd::DoubleKernel;

pd::Point<EK::FT> ex(const pd::Point<double>& p) { return {EK::from_double(p.x), EK::from_double(p.y)}; }

} // namespace

TEST_CASE("the filtered crossing order agrees with exact arithmetic near ties")
{
  // Two edges crossing the ray o -> v at almost the same point: X on the ray,
  // the second edge through X shifted by a few ulps along the ray.
  std::mt19937_64 rng(99);
  std::uniform_real_distribution<double> u(-1, 1);
  const int n = 200000;
  int wrong = 0, tested = 0;
  for (int i = 0; i < n; ++i) {
    const pd::Point<double> o{u(rng), u(rng)}, v{o.x + u(rng), o.y + u(rng)};
    const double s = 1.5 + std::abs(u(rng));
    const pd::Point<double> x{o.x + s * (v.x - o.x), o.y + s * (v.y - o.y)};
    pd::Point<double> x2 = x;
    x2.x = std::nextafter(x2.x, (i % 2) ? INFINITY : -INFINITY);
    auto edge_through = [&](const pd::Point<double>& c, pd::Point<double>& a, pd::Point<double>& b) {
      const double dx = u(rng), dy = u(rng);
      a = {c.x - dx, c.y - dy};
      b = {c.x + dx, c.y + dy};
    };
    pd::Point<double> a1, b1, a2, b2;
    edge_through(x, a1, b1);
    edge_through(x2, a2, b2);
    const int s1 = DK::orient(o, v, a1), s2 = DK::orient(o, v, a2);
    if (s1 * DK::orient(o, v, b1) >= 0 || s2 * DK::orient(o, v, b2) >= 0) continue;  // not proper crossings
    ++tested;
    const bool fast = pd::detail::crossing_before<DK>(o, v, a1, b1, s1, a2, b2, s2);
    const bool exact = pd::detail::crossing_before<EK>(ex(o), ex(v), ex(a1), ex(b1), s1, ex(a2), ex(b2), s2);
    wrong += (fast != exact);
  }
  MESSAGE("near-tie crossing comparisons: " << tested);
  CHECK(tested > n / 2);
  CHECK(wrong == 0);
}

TEST_CASE("double pipeline equals exact on generic families")
{
  for (pd::gen::Family f : pd::gen::all_families()) {
    if (f == pd::gen::Family::Grid) continue;  // exact incidences: see the next test
    for (std::uint64_t seed = 1; seed <= 15; ++seed) {
      CAPTURE(pd::gen::name(f));
      CAPTURE(seed);
      const auto in = pd::gen::generate(f, seed);
      const auto e = pd::decompose<EK>(pd::convert_polygon<EK>(in));
      const auto d = pd::decompose<DK>(in);
      REQUIRE(e.ok());
      REQUIRE(d.ok());
      CHECK(d.leaves == e.leaves);
      REQUIRE(d.pool.size() == e.pool.size());
      double worst = 0;
      for (std::size_t i = 0; i < e.pool.size(); ++i)
        worst = std::max({worst, std::abs(d.pool[i].x - e.pool[i].x.get_d()),
                          std::abs(d.pool[i].y - e.pool[i].y.get_d())});
      CHECK(worst < 1e-8);  // coordinates are at most 10 in magnitude
    }
  }
}

TEST_CASE("degenerate and adversarial inputs: double reports or succeeds, never returns garbage")
{
  std::vector<pd::gen::Family> fams = pd::gen::adversarial_families();
  fams.push_back(pd::gen::Family::Grid);
  for (pd::gen::Family f : fams)
    for (std::uint64_t seed = 1; seed <= 25; ++seed) {
      CAPTURE(pd::gen::name(f));
      CAPTURE(seed);
      const auto in = pd::gen::generate(f, seed);
      const auto d = pd::decompose<DK>(in);
      CHECK((d.diag.status == pd::Status::Ok || d.diag.status == pd::Status::NumericalFailure));
      if (d.ok()) {
        const auto rep = pd::check_pieces_approx<DK>(in, d.leaf_rings(), true, 1e-12, 1e-9);
        CHECK_MESSAGE(rep.ok, rep.failure);
      }
      // With RetryExact every input is decomposed.
      pd::PipelineOptions opt;
      opt.on_failure = pd::OnFailure::RetryExact;
      const auto r = pd::decompose<DK>(in, opt);
      CHECK(r.ok());
      CHECK(r.stats.retried_exact == !d.ok());
      // The naive kernel may fail or be wrong, but must not crash or hang.
      pd::PipelineOptions naive;
      naive.self_check = false;
      (void)pd::decompose<pd::NaiveDoubleKernel>(in, naive);
    }
}

TEST_CASE("check_pieces_approx rejects the typical rounding failures")
{
  using R = pd::Ring<double>;
  const pd::PolygonWithHoles<double> sq{{{0, 0}, {4, 0}, {4, 4}, {0, 4}}, {}};
  const R left{{0, 0}, {2, 0}, {2, 4}, {0, 4}}, right{{2, 0}, {4, 0}, {4, 4}, {2, 4}};
  CHECK(pd::check_pieces_approx<DK>(sq, {left, right}, true, 1e-12, 1e-9).ok);
  CHECK_FALSE(pd::check_pieces_approx<DK>(sq, {left}, true, 1e-12, 1e-9).ok);            // lost area
  const R bowtie{{2, 0}, {4, 4}, {4, 0}, {2, 4}};
  CHECK_FALSE(pd::check_pieces_approx<DK>(sq, {left, bowtie}, true, 1e-12, 1e-9).ok);    // not simple
  const R notch{{2, 0}, {4, 0}, {4, 4}, {3, 2}, {2, 4}};
  CHECK_FALSE(pd::check_pieces_approx<DK>(sq, {left, notch}, true, 1e-12, 1e-9).ok);     // not convex
  // Not convex is fine if convexity is not required (the notch loses 2 of 16 in area).
  CHECK(pd::check_pieces_approx<DK>(sq, {left, notch}, false, 1e-12, 0.2).ok);
}

#endif // PD_HAVE_GMP
