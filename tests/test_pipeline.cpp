// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#include "doctest.h"

#include <string>

#include "pd/generate.hpp"
#include "pd/io/obj.hpp"
#include "pd/kernel/double.hpp"
#include "pd/pipeline.hpp"
#ifdef PD_HAVE_GMP
#include "pd/kernel/exact.hpp"
#include "pd/verify.hpp"
#endif

#ifdef PD_HAVE_GMP

namespace {

using EK = pd::ExactKernel;
const std::string kData = PD_DATA_DIR;

// Runs the exact pipeline and certifies a convex partition of the input.
pd::PipelineResult<EK::FT> certified(const pd::PolygonWithHoles<double>& in)
{
  pd::PipelineOptions opt;
  opt.check = true;
  const auto e = pd::convert_polygon<EK>(in);
  auto r = pd::decompose<EK>(e, opt);
  REQUIRE_MESSAGE(r.ok(), r.diag.message);
  const auto rep = pd::verify_partition<EK>(e, r.leaf_rings(), /*convex=*/true);
  CHECK_MESSAGE(rep.ok, rep.failure);
  // Stage 1 output alone is a hole-free partition as well.
  const auto rep1 = pd::verify_partition<EK>(e, r.hole_free_rings(), false);
  CHECK_MESSAGE(rep1.ok, rep1.failure);
  CHECK(r.stats.convex_pieces == r.leaves.size());
  CHECK(r.stats.convex_pieces <= r.stats.reflex_vertices + r.stats.hole_free_pieces);
  return r;
}

} // namespace

TEST_CASE("pipeline on the bundled legacy input")
{
  pd::PolygonWithHoles<double> p;
  p.outer = pd::read_obj_ring(kData + "/legacy/compgon/cpolygon0.obj");
  for (int i = 0; i <= 6; ++i)
    p.holes.push_back(pd::read_obj_ring(kData + "/legacy/compgon/hole0-" + std::to_string(i) + ".obj"));
  const auto r = certified(p);
  CHECK(r.stats.hole_free_pieces == 16);
  MESSAGE("bundled: " << r.stats.hole_free_pieces << " hole-free pieces, " << r.stats.convex_pieces
                      << " convex pieces, " << r.stats.reflex_vertices << " reflex vertices");
}

TEST_CASE("pipeline on the degenerate set")
{
  for (const auto& d : pd::gen::degenerate_set()) {
    CAPTURE(d.first);
    certified(d.second);
  }
}

TEST_CASE("pipeline on generated inputs of every family")
{
#ifdef NDEBUG
  const std::uint64_t per_family = 60;
#else
  const std::uint64_t per_family = 8;
#endif
  for (pd::gen::Family f : pd::gen::all_families())
    for (std::uint64_t seed = 1; seed <= per_family; ++seed) {
      CAPTURE(pd::gen::name(f));
      CAPTURE(seed);
      certified(pd::gen::generate(f, seed));
    }
}

TEST_CASE("without holes the pipeline is the convex decomposition")
{
  const auto comb = pd::gen::comb(20);
  const auto r = certified({comb, {}});
  const auto c = pd::convex_decompose<EK>(pd::convert_ring<EK>(comb));
  CHECK(r.stats.hole_free_pieces == 1);
  CHECK(r.leaves == c.leaves);
}

TEST_CASE("invalid input is reported by the pipeline")
{
  const pd::PolygonWithHoles<double> p{{{0, 0}, {10, 0}, {10, 10}, {0, 10}},
                                       {{{2, 2}, {5, 2}, {5, 5}}, {{5, 5}, {8, 5}, {8, 8}}}};
  const auto r = pd::decompose<EK>(pd::convert_polygon<EK>(p));
  CHECK(r.diag.status == pd::Status::InvalidInput);
  CHECK(r.leaves.empty());
}

#endif // PD_HAVE_GMP

TEST_CASE("double-kernel pipeline on generated inputs")
{
  std::size_t ok = 0, total = 0;
  for (pd::gen::Family f : pd::gen::all_families())
    for (std::uint64_t seed = 1; seed <= 8; ++seed) {
      ++total;
      ok += pd::decompose<pd::DoubleKernel>(pd::gen::generate(f, seed)).ok();
    }
  MESSAGE("double kernel: " << ok << " of " << total << " generated inputs decomposed");
  CHECK(ok == total);
}
