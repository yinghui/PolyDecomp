// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
//
// Example: convex decomposition of a simple polygon, its statistics, and the
// decomposition tree, from which the input is rebuilt by gluing the pieces.
#include <iostream>
#include <numeric>

#include "pd/pd.hpp"

int main()
{
  // A comb with 5 teeth: 20 vertices, 8 reflex vertices.
  const pd::Ring<double> comb = pd::gen::comb(5);

  const auto d = pd::convex_decompose<pd::DoubleKernel>(comb);
  if (!d.ok()) return 1;
  std::cout << "reflex vertices: " << d.stats.reflex_vertices << ", splits: " << d.stats.splits
            << ", chords ending at a vertex: " << d.stats.vertex_hits
            << ", convex pieces: " << d.leaves.size() << '\n';

  // The root of the tree glues back to the input ring (pool ids 0 .. n-1).
  std::vector<std::uint32_t> ids(comb.size());
  std::iota(ids.begin(), ids.end(), 0u);
  const bool rebuilt = pd::reconstruct_ring(d.tree, d.leaves, 0) == ids;
  std::cout << "input rebuilt from the tree: " << (rebuilt ? "yes" : "no") << '\n';
  return rebuilt ? 0 : 1;
}
