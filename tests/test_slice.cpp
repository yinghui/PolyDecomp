// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
// Tests for pd::Slicer: exact cuts along arbitrary lines, any number of times.
#include "doctest.h"

#include <cmath>
#include <cstdint>
#include <random>
#include <vector>

#include "pd/generate.hpp"
#include "pd/geometry.hpp"
#include "pd/kernel/double.hpp"
#include "pd/validate.hpp"
#include "pd/verify.hpp"

#ifdef PD_HAVE_GMP
#include "pd/kernel/exact.hpp"
#include "pd/slice.hpp"

namespace {

using EK = pd::ExactKernel;
using P = pd::Point<double>;

pd::Ring<double> square(double x0, double y0, double x1, double y1)
{
  return {{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}};
}

// The pieces tile the input, each a polygon with the holes it still contains.
void certify(const pd::Slicer& s, const pd::PolygonWithHoles<double>& in)
{
  const auto rep = pd::verify_slices<EK>(pd::convert_polygon<EK>(in), s.pieces_exact(), s.symbolic());
  CHECK_MESSAGE(rep.ok, rep.failure);
}

// The convex decomposition of every piece is a certified convex partition of that piece.
void certify_convex(const pd::Slicer& s)
{
  for (std::size_t i = 0; i < s.size(); ++i) {
    CAPTURE(i);
    const auto r = s.convex(i);
    REQUIRE_MESSAGE(r.ok(), r.diag.message);
    const auto rep = pd::verify_partition<EK>(s.piece_exact(i), r.leaf_rings(), true);
    CHECK_MESSAGE(rep.ok, rep.failure);
  }
}

std::size_t holes_left(const pd::Slicer& s)
{
  std::size_t n = 0;
  for (const auto& p : s.pieces_exact()) n += p.holes.size();
  return n;
}

std::size_t input_bits(const pd::PolygonWithHoles<double>& in)
{
  std::size_t b = 0;
  auto ring = [&](const pd::Ring<double>& r) {
    for (const auto& v : r)
      b = std::max({b, pd::detail::bit_size(mpq_class(v.x)), pd::detail::bit_size(mpq_class(v.y))});
  };
  ring(in.outer);
  for (const auto& h : in.holes) ring(h);
  return b;
}

const pd::PolygonWithHoles<double> kSquareHole{square(0, 0, 10, 10), {square(4, 4, 6, 6)}};

} // namespace

TEST_CASE("slicer: a cut through a hole dissolves it")
{
  pd::Slicer s(kSquareHole);
  REQUIRE(s.ok());
  const auto r = s.cut(P{-1, 5.3}, P{11, 5.3});
  REQUIRE_MESSAGE(r.ok(), r.diag.message);
  CHECK(r.cut_pieces == 1);
  CHECK(r.pieces.size() == 2);
  CHECK(s.size() == 2);
  CHECK(holes_left(s) == 0);
  certify(s, kSquareHole);
  certify_convex(s);
}

TEST_CASE("slicer: a hole that the line misses stays with its piece")
{
  pd::Slicer s(kSquareHole);
  const auto r = s.cut(P{-1, 2}, P{11, 2});
  REQUIRE(r.ok());
  CHECK(s.size() == 2);
  CHECK(holes_left(s) == 1);
  certify(s, kSquareHole);
  certify_convex(s);
}

TEST_CASE("slicer: lines that miss the shape or only touch its boundary change nothing")
{
  pd::Slicer s(kSquareHole);
  for (const auto& [a, b] : std::vector<std::pair<P, P>>{{{-1, 20}, {11, 20}},     // misses
                                                         {{-1, 0}, {11, 0}},       // along an outer edge
                                                         {{-1, 1}, {1, -1}},       // through a corner from outside
                                                         {{10, -5}, {10, 15}}}) {  // along another outer edge
    const auto r = s.cut(a, b);
    REQUIRE_MESSAGE(r.ok(), r.diag.message);
    CHECK(r.cut_pieces == 0);
    CHECK(r.pieces.empty());
    CHECK(s.size() == 1);
  }
  certify(s, kSquareHole);
}

TEST_CASE("slicer: contacts with vertices and edges")
{
  SUBCASE("through vertices of the hole and of the outer ring")
  {
    pd::Slicer s(kSquareHole);
    const auto r = s.cut(P{0, 0}, P{10, 10});
    REQUIRE(r.ok());
    CHECK(r.vertex_contacts == 4);
    CHECK(s.size() == 2);
    certify(s, kSquareHole);
    certify_convex(s);
  }
  SUBCASE("along an edge of the hole")
  {
    pd::Slicer s(kSquareHole);
    const auto r = s.cut(P{-1, 4}, P{11, 4});
    REQUIRE(r.ok());
    CHECK(r.edge_contacts == 1);
    CHECK(s.size() == 2);
    CHECK(holes_left(s) == 0);
    certify(s, kSquareHole);
    certify_convex(s);
  }
}

TEST_CASE("slicer: repeating a cut, and cutting through a point made by an earlier cut")
{
  pd::Slicer s(kSquareHole);
  REQUIRE(s.cut(P{-1, 5.3}, P{11, 5.3}).ok());
  const auto again = s.cut(P{-1, 5.3}, P{11, 5.3});  // runs along the new edges only
  REQUIRE(again.ok());
  CHECK(again.cut_pieces == 0);
  CHECK(s.size() == 2);
  // The line through the outer corner (0, 0) and the crossing (4, 5.3) of the first cut with the hole
  // passes through vertices of both pieces.
  const auto r = s.cut(P{0, 0}, P{4, 5.3});
  REQUIRE(r.ok());
  CHECK(r.cut_pieces == 2);
  CHECK(r.vertex_contacts >= 2);
  // x = 4 runs along the hole's left edge, which both pieces now share with the line.
  const auto along = s.cut(P{4, -1}, P{4, 11});
  REQUIRE(along.ok());
  CHECK(along.edge_contacts >= 2);
  certify(s, kSquareHole);
  certify_convex(s);
}

TEST_CASE("slicer: a line grazing a hole at one vertex switches to the symbolic kernel")
{
  const pd::PolygonWithHoles<double> in{square(0, 0, 10, 10), {{{5, 7}, {4, 6}, {5, 5}, {6, 6}}}};
  pd::Slicer s(in);
  CHECK_FALSE(s.symbolic());
  const auto r = s.cut(P{-1, 7}, P{11, 7});
  REQUIRE_MESSAGE(r.ok(), r.diag.message);
  CHECK(r.dangles == 1);
  CHECK(s.symbolic());
  CHECK(s.size() == 2);
  CHECK(holes_left(s) == 0);
  certify(s, in);  // weakly simple pieces allowed
  certify_convex(s);
  // Later cuts work on the perturbed state.
  REQUIRE(s.cut(P{5.5, -1}, P{5.5, 11}).ok());
  REQUIRE(s.cut(P{-1, 3}, P{11, 8}).ok());
  certify(s, in);
  certify_convex(s);
}

TEST_CASE("slicer: one piece at a time, and exact lines")
{
  pd::Slicer s(kSquareHole);
  REQUIRE(s.cut(P{-1, 5.3}, P{11, 5.3}).ok());
  const std::size_t before = s.size();
  const auto r = s.cut(0, P{3, -1}, P{3, 11});  // only piece 0
  REQUIRE(r.ok());
  CHECK(r.cut_pieces == 1);
  CHECK(s.size() == before + 1);
  certify(s, kSquareHole);
  const auto bad = s.cut(99, P{3, -1}, P{3, 11});
  CHECK(bad.diag.status == pd::Status::InvalidInput);
  const auto same = s.cut(P{1, 1}, P{1, 1});
  CHECK(same.diag.status == pd::Status::InvalidInput);
  // An exact line through two points that are not doubles.
  const pd::Point<mpq_class> a{mpq_class(1, 3), mpq_class(0)}, b{mpq_class(1, 3), mpq_class(1)};
  REQUIRE(s.cut_exact(a, b).ok());
  // Points given as braced lists are doubles.
  REQUIRE(s.cut({-1, 7.5}, {11, 7.5}).ok());
  certify(s, kSquareHole);
}

TEST_CASE("slicer: lines that pass extremely close to a piece's bounding box")
{
  // The box test runs in floating point with an exact fallback; lines near box corners must
  // still cut, or keep, exactly the right pieces.
  const double t = std::ldexp(1.0, -40);
  pd::Slicer s(kSquareHole);
  const auto r = s.cut(P{10 - t, 10}, P{10, 10 - t});  // cuts a corner of size 2^-40 off the square
  REQUIRE(r.ok());
  CHECK(r.cut_pieces == 1);
  CHECK(s.size() == 2);
  const auto touch = s.cut(P{10 + t, 10 - t}, P{10 - t, 10 + t});  // through the corner (10, 10)
  REQUIRE(touch.ok());
  CHECK(touch.cut_pieces == 0);
  const pd::Point<mpq_class> a{mpq_class(10), mpq_class(0)}, b{mpq_class(10) + mpq_class(1, 1) / (mpq_class(1) << 200),
                                                               mpq_class(1)};
  const auto steep = s.cut_exact(a, b);  // through the corner (10, 0), leaving the square at once
  REQUIRE(steep.ok());
  CHECK(steep.cut_pieces == 0);
  certify(s, kSquareHole);
  certify_convex(s);
}

TEST_CASE("verify_slices rejects a hole that leaves its piece")
{
  using R = pd::Ring<mpq_class>;
  auto ring = [](std::initializer_list<std::pair<int, int>> pts) {
    R r;
    for (const auto& [x, y] : pts) r.push_back({mpq_class(x), mpq_class(y)});
    return r;
  };
  const pd::PolygonWithHoles<mpq_class> in{ring({{0, 0}, {10, 0}, {10, 10}, {0, 10}}),
                                           {ring({{4, 4}, {6, 4}, {6, 6}, {4, 6}})}};
  // The halves x <= 5 and x >= 5 tile the square, but the hole straddles x = 5, so giving it to the
  // left half alone does not tile the region, although a vertex of the hole lies inside that half.
  const pd::PolygonWithHoles<mpq_class> left{ring({{0, 0}, {5, 0}, {5, 10}, {0, 10}}),
                                             {ring({{4, 4}, {6, 4}, {6, 6}, {4, 6}})}};
  const pd::PolygonWithHoles<mpq_class> right{ring({{5, 0}, {10, 0}, {10, 10}, {5, 10}}), {}};
  const auto rep = pd::verify_slices<EK>(in, {left, right});
  CHECK_FALSE(rep.ok);
  // The same hole inside a piece that contains it is accepted.
  const pd::PolygonWithHoles<mpq_class> lower{ring({{0, 0}, {10, 0}, {10, 7}, {0, 7}}),
                                              {ring({{4, 4}, {6, 4}, {6, 6}, {4, 6}})}};
  const pd::PolygonWithHoles<mpq_class> upper{ring({{0, 7}, {10, 7}, {10, 10}, {0, 10}}), {}};
  CHECK(pd::verify_slices<EK>(in, {lower, upper}).ok);
  // Two pieces that both claim the hole are rejected.
  CHECK_FALSE(pd::verify_slices<EK>(in, {lower, pd::PolygonWithHoles<mpq_class>{upper.outer, lower.holes}}).ok);
}

TEST_CASE("slicer: an invalid shape is reported")
{
  const pd::PolygonWithHoles<double> touching{square(0, 0, 10, 10), {square(0, 4, 6, 6)}};
  pd::Slicer s(touching);
  CHECK_FALSE(s.ok());
  CHECK(s.diag().status == pd::Status::InvalidInput);
  CHECK(s.cut(P{-1, 5}, P{11, 5}).diag.status == pd::Status::InvalidInput);
}

TEST_CASE("slicer: many random cuts on cheese keep the tiling and bounded coordinates")
{
  const auto in = pd::gen::cheese(16, 3);
  pd::Slicer s(in);
  REQUIRE(s.ok());
  std::mt19937_64 rng(7);
  std::uniform_real_distribution<double> u(-8, 8);
  std::size_t line_bits = 0;
  for (int k = 0; k < 200; ++k) {
    const P a{u(rng), u(rng)}, b{u(rng), u(rng)};
    line_bits = std::max({line_bits, pd::detail::bit_size(mpq_class(a.x)), pd::detail::bit_size(mpq_class(a.y)),
                          pd::detail::bit_size(mpq_class(b.x)), pd::detail::bit_size(mpq_class(b.y))});
    const auto r = s.cut(a, b);
    REQUIRE_MESSAGE(r.ok(), r.diag.message);
    if (k % 20 == 19) certify(s, in);
  }
  certify(s, in);
  CHECK(s.size() > 100);
  // Every vertex is an input vertex or the intersection of two lines through input or cut points.
  const std::size_t L = std::max(input_bits(in), line_bits);
  CHECK(s.max_bits() <= 8 * L);
}

TEST_CASE("slicer: degenerate cut sequences on a grid of holes and on cheese")
{
  SUBCASE("axis-aligned lines through vertex coordinates")
  {
    const auto in = pd::gen::random_holes(16, 2);
    pd::Slicer s(in);
    REQUIRE(s.ok());
    std::vector<P> vs(in.outer.begin(), in.outer.end());
    for (const auto& h : in.holes) vs.insert(vs.end(), h.begin(), h.end());
    for (std::size_t k = 0; k < 40; ++k) {
      const P& v = vs[(k * 37) % vs.size()];
      const auto r = (k % 2 == 0) ? s.cut(P{v.x, -100}, P{v.x, 100}) : s.cut(P{-100, v.y}, P{100, v.y});
      REQUIRE_MESSAGE(r.ok(), r.diag.message);
    }
    certify(s, in);
    certify_convex(s);
  }
  SUBCASE("tangent lines at hole vertices (dangling contacts)")
  {
    const auto in = pd::gen::cheese(8, 5);
    pd::Slicer s(in);
    for (std::size_t j = 0; j < in.holes.size(); j += 2) {
      const auto& h = in.holes[j];
      const std::size_t m = h.size();
      // Through vertex 0, parallel to the chord of its neighbours: the convex hole lies on one side.
      const pd::Point<mpq_class> v{mpq_class(h[0].x), mpq_class(h[0].y)};
      const pd::Point<mpq_class> dir{mpq_class(h[1].x) - mpq_class(h[m - 1].x), mpq_class(h[1].y) - mpq_class(h[m - 1].y)};
      const auto r = s.cut_exact(v, {v.x + dir.x, v.y + dir.y});
      REQUIRE_MESSAGE(r.ok(), r.diag.message);
    }
    CHECK(s.symbolic());
    certify(s, in);
    certify_convex(s);
  }
}
TEST_CASE("slicer: a line whose coefficients underflow in doubles still cuts the piece")
{
  // The box filter rounds the exact line to doubles. Here a coefficient rounds to zero, so the
  // filter has no error bound and must not skip the piece that the line crosses.
  const auto pow2 = [](int k) {
    mpz_class z = 1;
    if (k >= 0) return mpq_class(z << k);
    return mpq_class(mpz_class(1), mpz_class(z << -k));
  };
  const mpq_class x0 = pow2(100), x1 = pow2(101), y0(3, 2), y1(5, 2);
  const pd::PolygonWithHoles<mpq_class> in{{{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}}, {}};
  REQUIRE(pd::validate<EK>(in).ok());
  const pd::Point<mpq_class> q{0, 1}, d{-3 * pow2(-999), -pow2(-1098)};
  const pd::Point<mpq_class> q2{q.x + d.x, q.y + d.y};
  REQUIRE(EK::orient(q, q2, {x1, y0}) * EK::orient(q, q2, {x0, y1}) < 0);  // the line crosses the box
  pd::Slicer s(in);
  const auto res = s.cut_exact(q, q2);
  REQUIRE_MESSAGE(res.ok(), res.diag.message);
  CHECK(s.size() == 2);
  const auto rep = pd::verify_slices<EK>(in, s.pieces_exact(), s.symbolic());
  CHECK_MESSAGE(rep.ok, rep.failure);
}
#endif
