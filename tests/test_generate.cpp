// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#include "doctest.h"

#include <sstream>

#include "pd/generate.hpp"
#include "pd/geometry.hpp"
#include "pd/io/wkt.hpp"
#include "pd/kernel/double.hpp"
#include "pd/validate.hpp"

namespace {

std::string wkt(const pd::PolygonWithHoles<double>& p)
{
  std::ostringstream os;
  pd::write_wkt(os, p);
  return os.str();
}

} // namespace

TEST_CASE("generators are deterministic and produce valid inputs")
{
  for (pd::gen::Family f : pd::gen::all_families()) {
    CAPTURE(pd::gen::name(f));
    CHECK(pd::gen::family_from_name(pd::gen::name(f)) == f);
    for (std::uint64_t seed = 1; seed <= 100; ++seed) {
      const auto p = pd::gen::generate(f, seed);
      CHECK(wkt(p) == wkt(pd::gen::generate(f, seed)));
      CHECK_FALSE(p.holes.empty());
      CHECK(pd::validate<pd::DoubleKernel>(p).ok());
    }
  }
  CHECK(wkt(pd::gen::generate(pd::gen::Family::StarHoles, 7)) !=
        wkt(pd::gen::generate(pd::gen::Family::CShapes, 7)));
  CHECK_THROWS_AS(pd::gen::family_from_name("nope"), std::invalid_argument);
}

TEST_CASE("C-shaped holes have their vertex-average centroid outside the ring")
{
  for (std::uint64_t seed = 1; seed <= 50; ++seed)
    for (const auto& h : pd::gen::generate(pd::gen::Family::CShapes, seed).holes) {
      double sx = 0, sy = 0;
      for (const auto& p : h) {
        sx += p.x;
        sy += p.y;
      }
      const pd::Point<double> c{sx / double(h.size()), sy / double(h.size())};
      CHECK(pd::locate<pd::DoubleKernel>(c, h) == pd::Location::Outside);
    }
}

TEST_CASE("comb and degenerate inputs")
{
  CHECK(pd::gen::comb(5).size() == 20);
  CHECK(pd::validate<pd::DoubleKernel>({pd::gen::comb(50), {}}).ok());
  for (const auto& d : pd::gen::degenerate_set()) {
    CAPTURE(d.first);
    CHECK(pd::validate<pd::DoubleKernel>(d.second).ok());
  }
}

TEST_CASE("star polygons and spiral holes are valid and deterministic")
{
  for (int n : {3, 10, 200})
    for (std::uint64_t s = 1; s <= 5; ++s) {
      const auto r = pd::gen::star_polygon(n, s);
      CHECK(r.size() <= static_cast<std::size_t>(n));
      CHECK(pd::validate<pd::DoubleKernel>({r, {}}).ok());
      CHECK(wkt({r, {}}) == wkt({pd::gen::star_polygon(n, s), {}}));
    }
  const auto sp = pd::gen::spiral_holes(40, 3);
  CHECK(sp.holes.size() == 40);
  CHECK(pd::validate<pd::DoubleKernel>(sp).ok());
  CHECK(wkt(sp) == wkt(pd::gen::spiral_holes(40, 3)));
}
