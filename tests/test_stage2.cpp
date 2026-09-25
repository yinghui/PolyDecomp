// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#include "doctest.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>

#include "pd/convex_decomp.hpp"
#include "pd/io/leaves.hpp"
#include "pd/io/obj.hpp"
#include "pd/kernel/double.hpp"
#ifdef PD_HAVE_GMP
#include "pd/kernel/exact.hpp"
#include "pd/verify.hpp"
#endif

namespace {

const std::string kData = PD_DATA_DIR;
[[maybe_unused]] std::string legacy(const std::string& rel) { return kData + "/legacy/" + rel; }

pd::Ring<double> R(std::initializer_list<std::pair<double, double>> pts)
{
  pd::Ring<double> r;
  for (const auto& p : pts) r.push_back({p.first, p.second});
  return r;
}

// Axis-aligned comb with `teeth` teeth of height h, counter-clockwise, n = 4 teeth.
// The ray from each left gap corner runs along y = 1 and ends exactly at the
// next gap's corner, so the decomposition takes the vertex-hit path.
[[maybe_unused]] pd::Ring<double> comb(int teeth, double h = 3)
{
  pd::Ring<double> r{{0, 0}, {2.0 * teeth - 1, 0}};
  for (int j = teeth - 1; j >= 0; --j) {
    r.push_back({2.0 * j + 1, h});
    r.push_back({2.0 * j, h});
    if (j > 0) {
      r.push_back({2.0 * j, 1});
      r.push_back({2.0 * j - 1, 1});
    }
  }
  return r;
}

const char* const kLegacyRings[] = {
    "cpolygon0", "hole0-0", "hole0-1", "hole0-2", "hole0-3", "hole0-4", "hole0-5",
    "hole0-6",   "hole1-0", "hole1-1", "pentagon", "polygon", "polygon2", "test1",
    "test2",     "test3",   "test4",   "test5",   "test6",   "test7"};

// Rings equal up to rotation, within tol per coordinate.
[[maybe_unused]] bool same_up_to_rotation(const pd::Ring<double>& a, const pd::Ring<double>& b, double tol)
{
  if (a.size() != b.size()) return false;
  const std::size_t n = a.size();
  for (std::size_t s = 0; s < n; ++s) {
    bool all = true;
    for (std::size_t i = 0; i < n && all; ++i)
      all = std::abs(a[i].x - b[(i + s) % n].x) <= tol && std::abs(a[i].y - b[(i + s) % n].y) <= tol;
    if (all) return true;
  }
  return false;
}

// Drops consecutive vertices closer than eps (the old code could emit them).
[[maybe_unused]] pd::Ring<double> dedup(const pd::Ring<double>& r, double eps = 1e-6)
{
  pd::Ring<double> out;
  for (const auto& p : r)
    if (out.empty() || std::hypot(p.x - out.back().x, p.y - out.back().y) > eps) out.push_back(p);
  while (out.size() > 1 && std::hypot(out.front().x - out.back().x, out.front().y - out.back().y) <= eps)
    out.pop_back();
  return out;
}

template <class K>
pd::ConvexDecomposition<typename K::FT> run(const pd::Ring<double>& in, bool check = true)
{
  pd::ConvexOptions opt;
  opt.check = check;
  return pd::convex_decompose<K>(pd::convert_ring<K>(in), opt);
}

} // namespace

#ifdef PD_HAVE_GMP
#define PD_TEST_KERNELS pd::DoubleKernel, pd::ExactKernel
#else
#define PD_TEST_KERNELS pd::DoubleKernel
#endif

TEST_CASE_TEMPLATE("convex inputs are returned whole", K, PD_TEST_KERNELS)
{
  const auto d = run<K>(R({{0, 0}, {4, 0}, {4, 4}, {0, 4}}));
  REQUIRE(d.ok());
  CHECK(d.leaves.size() == 1);
  CHECK(d.stats.splits == 0);
  REQUIRE(d.tree.size() == 1);
  CHECK(d.tree[0].leaf == 0);
  // Clockwise input is reoriented; still one piece.
  CHECK(run<K>(R({{0, 0}, {0, 4}, {4, 4}, {4, 0}})).leaves.size() == 1);
}

TEST_CASE_TEMPLATE("invalid inputs are rejected, not decomposed", K, PD_TEST_KERNELS)
{
  const auto bowtie = run<K>(R({{0, 0}, {4, 4}, {4, 0}, {0, 4}}));
  CHECK(bowtie.diag.status == pd::Status::InvalidInput);
  CHECK_FALSE(bowtie.diag.issues.empty());
  CHECK(bowtie.leaves.empty());
  CHECK(run<K>(R({{0, 0}, {4, 4}})).diag.status == pd::Status::InvalidInput);
}

TEST_CASE_TEMPLATE("the leaf count follows the reflex count", K, PD_TEST_KERNELS)
{
  // A vertex hit that resolves the hit vertex too: in this H shape the ray from
  // (4, 2) ends at the reflex vertex (2, 2), which becomes convex on both sides.
  const auto h = R({{0, 0}, {6, 0}, {6, 2}, {4, 2}, {4, 4}, {6, 4}, {6, 6}, {0, 6},
                    {0, 4}, {2, 4}, {2, 2}, {0, 2}});
  const auto d = run<K>(h);
  REQUIRE(d.ok());
  CHECK(d.stats.reflex_vertices == 4);
  CHECK(d.stats.double_resolutions >= 1);
  CHECK(d.leaves.size() == d.stats.reflex_vertices + 1 - d.stats.double_resolutions);

  // Straight vertices hit exactly by chords.
  const auto u = R({{0, 0}, {4, 0}, {6, 0}, {6, 4}, {4, 4}, {4, 2}, {2, 2}, {2, 4}, {0, 4}, {0, 2}});
  const auto e = run<K>(u);
  REQUIRE(e.ok());
  CHECK(e.stats.vertex_hits >= 1);
  CHECK(e.leaves.size() == e.stats.reflex_vertices + 1 - e.stats.double_resolutions);
}

#ifdef PD_HAVE_GMP

TEST_CASE("exact decomposition of every legacy ring is a certified convex partition")
{
  using K = pd::ExactKernel;
  for (const char* name : kLegacyRings) {
    CAPTURE(name);
    const auto in = pd::read_obj_ring(legacy(std::string("simp/") + name + ".obj"));
    const auto d = run<K>(in);
    REQUIRE(d.ok());
    CHECK(d.leaves.size() == d.stats.reflex_vertices + 1 - d.stats.double_resolutions);
    CHECK(d.leaves.size() <= d.stats.reflex_vertices + 1);
    const auto rep =
        pd::verify_partition<K>({pd::convert_ring<K>(in), {}}, d.leaf_rings(), /*convex=*/true);
    CHECK_MESSAGE(rep.ok, rep.failure);
    // Gluing the leaves back along the chords gives the input ring.
    std::vector<std::uint32_t> ids(in.size());
    std::iota(ids.begin(), ids.end(), 0u);
    CHECK(pd::reconstruct_ring(d.tree, d.leaves, 0) == ids);
  }
}

TEST_CASE("regression against the legacy implementation's leaves")
{
  using K = pd::ExactKernel;
  const std::pair<const char*, std::size_t> cases[] = {
      {"cpolygon0", 71}, {"pentagon", 2}, {"polygon", 17}, {"test5", 74}, {"test1", 14}};
  for (const auto& c : cases) {
    CAPTURE(c.first);
    const auto in = pd::read_obj_ring(legacy(std::string("simp/") + c.first + ".obj"));
    const auto old = pd::read_leaves(legacy(std::string("simp/") + c.first + ".obj.leaves"));
    const auto d = run<K>(in);
    REQUIRE(d.ok());
    CHECK(d.leaves.size() == c.second);
    REQUIRE(old.size() == c.second);
    // Every old leaf is one of the new leaves (float legacy output: tolerance 1e-4).
    std::vector<bool> used(d.leaves.size(), false);
    std::size_t matched = 0;
    for (const auto& o : old) {
      const auto oc = dedup(o);
      for (std::size_t j = 0; j < d.leaves.size(); ++j) {
        if (used[j]) continue;
        if (same_up_to_rotation(oc, pd::to_double_ring<K>(d.leaf_ring(j)), 1e-4)) {
          used[j] = true;
          ++matched;
          break;
        }
      }
    }
    CHECK(matched == c.second);
  }
}

TEST_CASE("double and exact kernels agree on the legacy rings")
{
  for (const char* name : kLegacyRings) {
    CAPTURE(name);
    const auto in = pd::read_obj_ring(legacy(std::string("simp/") + name + ".obj"));
    const auto e = run<pd::ExactKernel>(in);
    const auto d = run<pd::DoubleKernel>(in, false);
    REQUIRE(e.ok());
    REQUIRE(d.ok());
    REQUIRE(d.leaves.size() == e.leaves.size());
    CHECK(d.leaves == e.leaves);  // same pieces, as vertex ids
    REQUIRE(d.pool.size() == e.pool.size());
    double worst = 0;  // largest coordinate difference over all (input and constructed) points
    for (std::size_t j = 0; j < e.pool.size(); ++j)
      worst = std::max({worst, std::abs(d.pool[j].x - e.pool[j].x.get_d()),
                        std::abs(d.pool[j].y - e.pool[j].y.get_d())});
    CHECK(worst < 1e-9);
  }
}

TEST_CASE("axis-aligned combs: vertex hits at scale")
{
  using K = pd::ExactKernel;
#ifdef NDEBUG
  const int sizes[] = {2, 10, 100, 2500};  // n up to 10 000
#else
  const int sizes[] = {2, 10, 100, 250};
#endif
  for (int teeth : sizes) {
    CAPTURE(teeth);
    const auto in = comb(teeth);
    const auto d = run<K>(in, teeth <= 100);
    REQUIRE(d.ok());
    CHECK(d.stats.input_vertices == static_cast<std::size_t>(4 * teeth));
    CHECK(d.stats.reflex_vertices == static_cast<std::size_t>(2 * (teeth - 1)));
    if (teeth >= 3) CHECK(d.stats.vertex_hits > 0);  // needs a next gap corner to hit
    CHECK(d.leaves.size() == d.stats.reflex_vertices + 1 - d.stats.double_resolutions);
    const auto rep = pd::verify_partition<K>({pd::convert_ring<K>(in), {}}, d.leaf_rings(), true);
    CHECK_MESSAGE(rep.ok, rep.failure);
    // Integer data: the double kernel makes the same decisions.
    const auto dd = run<pd::DoubleKernel>(in, false);
    REQUIRE(dd.ok());
    CHECK(dd.leaves == d.leaves);
  }
}

#endif // PD_HAVE_GMP
