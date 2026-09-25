// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
//
// Example: decompose a polygon with holes into convex pieces with exact
// arithmetic, certify the result, and print the pieces as WKT.
#include <iostream>

#include "pd/pd.hpp"

int main()
{
  // A square with two holes; the second sits in the pocket of the first
  // (entangled holes). Coordinates are doubles; ring orientation is free.
  const pd::PolygonWithHoles<double> input{
      {{0, 0}, {10, 0}, {10, 10}, {0, 10}},
      {{{1, 1}, {9, 1}, {9, 3}, {3, 3}, {3, 7}, {9, 7}, {9, 9}, {1, 9}},
       {{5, 4}, {8, 4}, {8, 6}, {5, 6}}}};

  using K = pd::ExactKernel;
  const auto exact_input = pd::convert_polygon<K>(input);  // exact: no rounding
  const auto result = pd::decompose<K>(exact_input);
  if (!result.ok()) {
    std::cerr << "failed: " << pd::to_string(result.diag.status) << ": " << result.diag.message << '\n';
    return 1;
  }

  // Independent exact certificate: the pieces are convex and tile the input.
  const auto cert = pd::verify_partition<K>(exact_input, result.leaf_rings(), /*require_convex=*/true);
  std::cout << result.stats.hole_free_pieces << " hole-free pieces, " << result.leaves.size()
            << " convex pieces, certified: " << (cert.ok ? "yes" : "no") << '\n';

  std::vector<pd::Ring<double>> pieces;
  for (const auto& r : result.leaf_rings()) pieces.push_back(pd::to_double_ring<K>(r));
  pd::write_wkt_multipolygon(std::cout, pieces);
  std::cout << '\n';
  return cert.ok ? 0 : 1;
}
