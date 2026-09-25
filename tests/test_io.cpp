// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#include "doctest.h"

#include <random>
#include <sstream>
#include <string>

#include "pd/io/leaves.hpp"
#include "pd/io/obj.hpp"
#include "pd/io/wkt.hpp"
#ifdef PD_HAVE_GMP
#include "pd/io/rational.hpp"
#endif

namespace {

const std::string kData = PD_DATA_DIR;
std::string legacy(const std::string& rel) { return kData + "/legacy/" + rel; }

} // namespace

TEST_CASE("WKT parsing")
{
  const auto polys = pd::parse_wkt(
      "POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0), (2 2, 4 2, 4 4, 2 2))\n"
      "multipolygon (((0 0 5, 1 0 5, 0 1 5, 0 0 5)), ((5 5, 6 5, 5 6)))\n"
      "POLYGON EMPTY");
  REQUIRE(polys.size() == 4);
  CHECK(polys[0].outer.size() == 4);  // closing vertex dropped
  REQUIRE(polys[0].holes.size() == 1);
  CHECK(polys[0].holes[0].size() == 3);
  CHECK(polys[1].outer.size() == 3);  // Z ignored
  CHECK(polys[2].outer.size() == 3);  // ring left open in the text
  CHECK(polys[3].outer.empty());
  CHECK(polys[0].outer[2].x == 10.0);

  CHECK_THROWS_AS(pd::parse_wkt("LINESTRING (0 0, 1 1)"), pd::ParseError);
  CHECK_THROWS_AS(pd::parse_wkt("POLYGON ((0 0, 1 x))"), pd::ParseError);
  CHECK_THROWS_AS(pd::parse_wkt("POLYGON ((0 0, 1 1)"), pd::ParseError);
}

TEST_CASE("WKT round trip is bit-exact")
{
  std::mt19937_64 rng(3);
  std::uniform_real_distribution<double> u(-1e6, 1e6);
  pd::PolygonWithHoles<double> p;
  for (int i = 0; i < 50; ++i) p.outer.push_back({u(rng), u(rng)});
  p.outer.push_back({0.1, -0.30000000000000004});
  p.outer.push_back({1e-300, 5e-324});
  p.holes.push_back({{1, 2}, {3, 4}, {5, 7}});
  std::ostringstream os;
  pd::write_wkt(os, p);
  const auto back = pd::parse_wkt(os.str());
  REQUIRE(back.size() == 1);
  REQUIRE(back[0].outer.size() == p.outer.size());
  for (std::size_t i = 0; i < p.outer.size(); ++i) {
    CHECK(back[0].outer[i].x == p.outer[i].x);
    CHECK(back[0].outer[i].y == p.outer[i].y);
  }
  REQUIRE(back[0].holes.size() == 1);
  CHECK(back[0].holes[0].size() == 3);

  std::ostringstream ms;
  pd::write_wkt_multipolygon(ms, {p.holes[0], p.holes[0]});
  CHECK(pd::parse_wkt(ms.str()).size() == 2);
}

TEST_CASE("legacy hole-elimination inputs load")
{
  const auto outer = pd::read_obj_ring(legacy("compgon/cpolygon0.obj"));
  CHECK(outer.size() == 150);
  for (int i = 0; i <= 6; ++i)
    CHECK(pd::read_obj_ring(legacy("compgon/hole0-" + std::to_string(i) + ".obj")).size() >= 3);
  CHECK(pd::read_obj_ring(legacy("compgon/test10/p.obj")).size() >= 3);
  CHECK(pd::read_obj_ring(legacy("compgon/test3/h3.obj")).size() >= 3);
}

TEST_CASE("legacy convex-decomposition inputs load")
{
  CHECK(pd::read_obj_ring(legacy("simp/cpolygon0.obj")).size() == 150);
  CHECK(pd::read_obj_ring(legacy("simp/pentagon.obj")).size() == 5);
  CHECK(pd::read_obj_ring(legacy("simp/polygon.obj")).size() == 43);
  CHECK(pd::read_obj_ring(legacy("simp/test5.obj")).size() == 150);
  // ctest1 holds three rings in one file.
  CHECK(pd::read_obj_rings(legacy("simp/ctest1.obj")).size() == 3);
  CHECK_THROWS_AS(pd::read_obj_ring(legacy("simp/ctest1.obj")), pd::ParseError);
  CHECK_THROWS_AS(pd::read_obj_ring(legacy("simp/no-such-file.obj")), pd::ParseError);
}

TEST_CASE("legacy regression leaves load")
{
  CHECK(pd::read_leaves(legacy("simp/cpolygon0.obj.leaves")).size() == 71);
  CHECK(pd::read_leaves(legacy("simp/pentagon.obj.leaves")).size() == 2);
  CHECK(pd::read_leaves(legacy("simp/polygon.obj.leaves")).size() == 17);
  CHECK(pd::read_leaves(legacy("simp/test5.obj.leaves")).size() == 74);
  CHECK(pd::read_leaves(legacy("simp/test1.obj.leaves")).size() == 14);
}

#ifdef PD_HAVE_GMP
TEST_CASE("rational format round trip")
{
  pd::PolygonWithHoles<mpq_class> p;
  p.outer = {{mpq_class(0), mpq_class(0)}, {mpq_class(1, 3), mpq_class(0)},
             {mpq_class(1, 3), mpq_class(-7, 2)}};
  p.holes.push_back({{mpq_class(1, 10), mpq_class(1, 10)},
                     {mpq_class(2, 10), mpq_class(1, 10)},
                     {mpq_class(2, 10), mpq_class(3, 10)}});
  std::stringstream ss;
  pd::write_rational(ss, p);
  const auto back = pd::read_rational(ss);
  REQUIRE(back.outer.size() == 3);
  CHECK(back.outer[2].y == mpq_class(-7, 2));
  REQUIRE(back.holes.size() == 1);
  CHECK(back.holes[0][0].x == mpq_class(1, 10));

  std::istringstream bad1("polydecomp-rational 1 polygon\n1\n3\n0 0\n1/0 0\n0 1\n");
  CHECK_THROWS_AS(pd::read_rational(bad1), pd::ParseError);
  std::istringstream bad2("polydecomp-rational 1 polygon\n1\n3\n0 0\nabc 0\n0 1\n");
  CHECK_THROWS_AS(pd::read_rational(bad2), pd::ParseError);
}
#endif
