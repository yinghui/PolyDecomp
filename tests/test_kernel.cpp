// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#include "doctest.h"

#include <cmath>
#include <random>

#include "pd/kernel/double.hpp"
#ifdef PD_HAVE_GMP
#include "pd/kernel/exact.hpp"
#endif

using pd::Point;

#ifdef PD_HAVE_GMP
#define PD_ALL_KERNELS pd::DoubleKernel, pd::NaiveDoubleKernel, pd::ExactKernel
#else
#define PD_ALL_KERNELS pd::DoubleKernel, pd::NaiveDoubleKernel
#endif

TEST_CASE_TEMPLATE("orient on clear-cut inputs", K, PD_ALL_KERNELS)
{
  using FT = typename K::FT;
  const Point<FT> a{FT(0), FT(0)}, b{FT(4), FT(0)};
  CHECK(K::orient(a, b, {FT(1), FT(3)}) == 1);
  CHECK(K::orient(a, b, {FT(1), FT(-3)}) == -1);
  CHECK(K::orient(a, b, {FT(9), FT(0)}) == 0);
  CHECK(K::orient(a, a, {FT(9), FT(5)}) == 0);
}

#ifdef PD_HAVE_GMP

namespace {

int exact_orient(const Point<double>& a, const Point<double>& b, const Point<double>& c)
{
  using E = pd::ExactKernel;
  return E::orient({E::from_double(a.x), E::from_double(a.y)},
                   {E::from_double(b.x), E::from_double(b.y)},
                   {E::from_double(c.x), E::from_double(c.y)});
}

double nudge(double v, int ulps)
{
  for (; ulps > 0; --ulps) v = std::nextafter(v, INFINITY);
  for (; ulps < 0; ++ulps) v = std::nextafter(v, -INFINITY);
  return v;
}

} // namespace

TEST_CASE("filtered double orient equals the exact sign on near-degenerate triples")
{
  // c is the rounded point a + t (b - a), nudged by up to two ulps, so the
  // triples are collinear or within a few ulps of it.
  std::mt19937_64 rng(20260925);
  std::uniform_real_distribution<double> u(-1.0, 1.0), t01(0.0, 1.0);
  const double scales[] = {1.0, 1e-3, 1e3, 1e8};

  const int n = 1000000;
  int filtered_wrong = 0, exact_path_wrong = 0, naive_wrong = 0, collinear = 0;
  for (int i = 0; i < n; ++i) {
    const double s = scales[i % 4];
    const Point<double> a{u(rng) * s, u(rng) * s}, b{u(rng) * s, u(rng) * s};
    const double t = t01(rng);
    Point<double> c{a.x + t * (b.x - a.x), a.y + t * (b.y - a.y)};
    c.x = nudge(c.x, i % 5 - 2);

    const int want = exact_orient(a, b, c);
    collinear += (want == 0);
    filtered_wrong += (pd::DoubleKernel::orient(a, b, c) != want);
    exact_path_wrong += (pd::detail::orient2d_exact(a.x, a.y, b.x, b.y, c.x, c.y) != want);
    naive_wrong += (pd::NaiveDoubleKernel::orient(a, b, c) != want);
  }
  MESSAGE("near-degenerate triples: " << n << ", exactly collinear: " << collinear
                                      << ", naive-double wrong signs: " << naive_wrong);
  CHECK(filtered_wrong == 0);
  CHECK(exact_path_wrong == 0);
}

TEST_CASE("filtered double orient is exact on collinear lattice triples and their nudges")
{
  // a, a + d and a + t d on an integer lattice, scaled by a power of two, so
  // the triple is exactly collinear; then nudged by one ulp in x or y.
  std::mt19937_64 rng(42);
  std::uniform_int_distribution<long> coord(-(1L << 20), 1L << 20);
  std::uniform_int_distribution<int> tdist(-4, 4), expo(-60, 60);
  const int n = 200000;
  int filtered_wrong = 0, zero_seen = 0, naive_wrong = 0;
  for (int i = 0; i < n; ++i) {
    const double sc = std::ldexp(1.0, expo(rng));
    const double ax = double(coord(rng)), ay = double(coord(rng));
    const double dx = double(coord(rng)), dy = double(coord(rng));
    const int t = tdist(rng);
    const Point<double> a{ax * sc, ay * sc}, b{(ax + dx) * sc, (ay + dy) * sc};
    Point<double> c{(ax + t * dx) * sc, (ay + t * dy) * sc};
    if (i % 3 == 1) c.x = nudge(c.x, (i % 2) ? 1 : -1);
    if (i % 3 == 2) c.y = nudge(c.y, (i % 2) ? 1 : -1);
    const int want = exact_orient(a, b, c);
    zero_seen += (want == 0);
    filtered_wrong += (pd::DoubleKernel::orient(a, b, c) != want);
    naive_wrong += (pd::NaiveDoubleKernel::orient(a, b, c) != want);
  }
  MESSAGE("lattice triples: " << n << ", exactly collinear: " << zero_seen
                              << ", naive-double wrong signs: " << naive_wrong);
  CHECK(zero_seen > n / 4);
  CHECK(filtered_wrong == 0);
}

TEST_CASE("classroom example: 256 x 256 grid near the line y = x")
{
  // Kettner et al., "Classroom examples of robustness problems in geometric
  // computations" (2008): p = (0.5 + i u, 0.5 + j u), q = (12, 12), r = (24, 24).
  const double ulp = std::ldexp(1.0, -53);
  const Point<double> q{12.0, 12.0}, r{24.0, 24.0};
  int filtered_wrong = 0, naive_wrong = 0;
  for (int i = 0; i < 256; ++i)
    for (int j = 0; j < 256; ++j) {
      const Point<double> p{0.5 + i * ulp, 0.5 + j * ulp};
      const int want = exact_orient(p, q, r);
      filtered_wrong += (pd::DoubleKernel::orient(p, q, r) != want);
      naive_wrong += (pd::NaiveDoubleKernel::orient(p, q, r) != want);
    }
  MESSAGE("classroom grid: naive-double wrong signs: " << naive_wrong << " of 65536");
  CHECK(filtered_wrong == 0);
  CHECK(naive_wrong > 0);  // the example is known to defeat plain floating point
}

#endif // PD_HAVE_GMP
