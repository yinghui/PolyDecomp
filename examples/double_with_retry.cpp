// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
//
// Example: the fast double-precision pipeline on a hard input (edges bent by a
// few units in the last place), with an exact rerun if it detects a problem.
#include <cmath>
#include <iostream>

#include "pd/pd.hpp"

int main()
{
  // Outer ring whose bottom edge has a vertex one ulp above the line y = 1,
  // and a hole whose top edge is split by a vertex one ulp below y = 6.
  // (Coordinates must be 0 or between 2^-400 and 2^400 for the double kernel.)
  pd::PolygonWithHoles<double> input{{{0, 1}, {5, std::nextafter(1.0, 2.0)}, {10, 1}, {10, 10}, {0, 10}},
                                     {{{2, 2}, {8, 2}, {8, 6}, {5, std::nextafter(6.0, 0.0)}, {2, 6}}}};

  pd::PipelineOptions opt;
#ifdef PD_HAVE_GMP
  opt.on_failure = pd::OnFailure::RetryExact;  // needs the exact kernel
#endif
  const auto result = pd::decompose<pd::DoubleKernel>(input, opt);
  std::cout << "status: " << pd::to_string(result.diag.status) << ", convex pieces: " << result.leaves.size()
            << (result.stats.retried_exact ? " (exact rerun)" : "") << '\n';
  return result.ok() ? 0 : 1;
}
