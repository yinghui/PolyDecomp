// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#include "doctest.h"

#include "pd/geometry.hpp"
#include "pd/kernel/double.hpp"
#ifdef PD_HAVE_GMP
#include "pd/kernel/exact.hpp"
#define PD_TEST_KERNELS pd::DoubleKernel, pd::ExactKernel
#else
#define PD_TEST_KERNELS pd::DoubleKernel
#endif

namespace {

template <class K>
pd::Ring<typename K::FT> ring(std::initializer_list<std::pair<double, double>> pts)
{
  pd::Ring<double> r;
  for (const auto& p : pts) r.push_back({p.first, p.second});
  return pd::convert_ring<K>(r);
}

template <class K>
pd::Point<typename K::FT> pt(double x, double y)
{
  return {K::from_double(x), K::from_double(y)};
}

} // namespace

TEST_CASE_TEMPLATE("ring orientation and area", K, PD_TEST_KERNELS)
{
  using FT = typename K::FT;
  auto sq = ring<K>({{0, 0}, {4, 0}, {4, 4}, {0, 4}});
  CHECK(pd::ring_orientation<K>(sq) == 1);
  CHECK(pd::signed_area2(sq) == FT(32));
  pd::Ring<FT> rev(sq.rbegin(), sq.rend());
  CHECK(pd::ring_orientation<K>(rev) == -1);
  pd::make_ccw<K>(rev);
  CHECK(pd::ring_orientation<K>(rev) == 1);
  pd::make_cw<K>(rev);
  CHECK(pd::ring_orientation<K>(rev) == -1);

  // Concave ring whose lexicographically smallest vertex is still convex.
  auto c = ring<K>({{0, 0}, {4, 0}, {4, 4}, {2, 1}, {0, 4}});
  CHECK(pd::ring_orientation<K>(c) == 1);
}

TEST_CASE_TEMPLATE("point location, including boundary and degenerate rays", K, PD_TEST_KERNELS)
{
  using pd::Location;
  auto sq = ring<K>({{0, 0}, {4, 0}, {4, 4}, {0, 4}});
  CHECK(pd::locate<K>(pt<K>(2, 2), sq) == Location::Inside);
  CHECK(pd::locate<K>(pt<K>(5, 2), sq) == Location::Outside);
  CHECK(pd::locate<K>(pt<K>(-1, 2), sq) == Location::Outside);
  CHECK(pd::locate<K>(pt<K>(0, 0), sq) == Location::OnBoundary);   // vertex
  CHECK(pd::locate<K>(pt<K>(2, 0), sq) == Location::OnBoundary);   // horizontal edge
  CHECK(pd::locate<K>(pt<K>(4, 2), sq) == Location::OnBoundary);   // vertical edge

  // Notched ring: the ray from (1, 2) passes exactly through the reflex vertex (2, 2).
  auto notch = ring<K>({{0, 0}, {4, 0}, {4, 4}, {2, 2}, {0, 4}});
  CHECK(pd::locate<K>(pt<K>(1, 2), notch) == Location::Inside);
  CHECK(pd::locate<K>(pt<K>(3, 3.5), notch) == Location::Outside);
  CHECK(pd::locate<K>(pt<K>(3, 3), notch) == Location::OnBoundary);

  // L-shape: the ray from (1, 2) runs along the horizontal edge (4,2)-(2,2).
  auto ell = ring<K>({{0, 0}, {4, 0}, {4, 2}, {2, 2}, {2, 4}, {0, 4}});
  CHECK(pd::locate<K>(pt<K>(1, 2), ell) == Location::Inside);
  CHECK(pd::locate<K>(pt<K>(3, 2), ell) == Location::OnBoundary);
  CHECK(pd::locate<K>(pt<K>(3, 3), ell) == Location::Outside);
  CHECK(pd::locate<K>(pt<K>(-1, 2), ell) == Location::Outside);

  // Orientation does not matter.
  pd::Ring<typename K::FT> rev(ell.rbegin(), ell.rend());
  CHECK(pd::locate<K>(pt<K>(1, 2), rev) == Location::Inside);
  CHECK(pd::locate<K>(pt<K>(3, 3), rev) == Location::Outside);
}

TEST_CASE_TEMPLATE("segment contact", K, PD_TEST_KERNELS)
{
  auto p = [](double x, double y) { return pt<K>(x, y); };
  CHECK(pd::segments_touch<K>(p(0, 0), p(4, 4), p(0, 4), p(4, 0)));      // proper crossing
  CHECK(pd::segments_touch<K>(p(0, 0), p(4, 0), p(4, 0), p(4, 4)));      // shared endpoint
  CHECK(pd::segments_touch<K>(p(0, 0), p(4, 0), p(2, 0), p(2, 3)));      // T-junction
  CHECK(pd::segments_touch<K>(p(0, 0), p(4, 0), p(3, 0), p(6, 0)));      // collinear overlap
  CHECK_FALSE(pd::segments_touch<K>(p(0, 0), p(4, 0), p(5, 0), p(6, 0))); // collinear, apart
  CHECK_FALSE(pd::segments_touch<K>(p(0, 0), p(4, 0), p(0, 1), p(4, 1))); // parallel
  CHECK_FALSE(pd::segments_touch<K>(p(0, 0), p(4, 0), p(2, 1), p(3, 5))); // apart
  CHECK(pd::on_segment<K>(p(0, 0), p(4, 4), p(1, 1)));
  CHECK_FALSE(pd::on_segment<K>(p(0, 0), p(4, 4), p(5, 5)));
}

#ifdef PD_HAVE_GMP
TEST_CASE("exact conversion keeps the binary value of a double")
{
  const pd::Ring<double> r{{0.1, 0.2}, {1e-300, 3.5}, {-7.25, 1e300}};
  const auto e = pd::convert_ring<pd::ExactKernel>(r);
  const auto back = pd::to_double_ring<pd::ExactKernel>(e);
  for (std::size_t i = 0; i < r.size(); ++i) {
    CHECK(back[i].x == r[i].x);
    CHECK(back[i].y == r[i].y);
  }
  CHECK(e[0].x != mpq_class(1, 10));
  const pd::Ring<double> bad{{0.0, 0.0}, {1.0, std::numeric_limits<double>::quiet_NaN()}};
  CHECK_THROWS_AS(pd::convert_ring<pd::ExactKernel>(bad), std::domain_error);
}
#endif
