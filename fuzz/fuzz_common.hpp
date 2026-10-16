// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#pragma once
// Shared oracle for the fuzz targets: on every valid input the exact pipeline
// must succeed and pass the exact certificate (convex and hole-free pieces), and
// the double pipeline and double hole elimination, each with RetryExact, must
// succeed. Any violation aborts (a finding for libFuzzer).

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>

#include "pd/kernel/double.hpp"
#include "pd/kernel/exact.hpp"
#include "pd/pipeline.hpp"
#include "pd/validate.hpp"
#include "pd/verify.hpp"

namespace pd_fuzz {

[[noreturn]] inline void bug(const char* what, const std::string& detail)
{
  std::fprintf(stderr, "FUZZ FINDING: %s: %s\n", what, detail.c_str());
  std::abort();
}

/// Reads the fuzzer's bytes in order; past the end it yields zeros.
struct Reader {
  const std::uint8_t* d;
  std::size_t n, i = 0;
  bool more() const { return i < n; }
  std::uint8_t next() { return i < n ? d[i++] : 0; }
};

/// A star-shaped ring on the integer grid (format in fuzz_grid.cpp): a centre byte
/// (x = low 4 bits, y = high 4 bits, on a 4x grid), a vertex-count byte (3..14), and one
/// byte per vertex (angle step = low 4 bits + 1 in units of 1/64 turn, radius = high 4
/// bits + 1).
inline pd::Ring<double> grid_ring(Reader& in, double scale, double cx0, double cy0)
{
  const std::uint8_t c = in.next();
  const double cx = cx0 + (c & 15) * scale / 4, cy = cy0 + (c >> 4) * scale / 4;
  const int k = 3 + in.next() % 12;
  pd::Ring<double> r;
  int angle = 0;
  for (int j = 0; j < k && angle < 64; ++j) {
    const std::uint8_t b = in.next();
    const double rad = ((b >> 4) + 1) * scale / 4;
    const double a = 6.283185307179586 * angle / 64;
    // Round to the grid: exact integer coordinates.
    r.push_back({std::round(cx + rad * std::cos(a)), std::round(cy + rad * std::sin(a))});
    angle += (b & 15) + 1;
  }
  return r;
}

/// A polygon with up to four holes from the bytes (byte 0: number of holes).
inline pd::PolygonWithHoles<double> grid_polygon(Reader& in)
{
  const int holes = in.next() % 5;
  pd::PolygonWithHoles<double> p;
  p.outer = grid_ring(in, 16, 64, 64);
  for (int h = 0; h < holes && in.more(); ++h) p.holes.push_back(grid_ring(in, 2, 40, 40));
  return p;
}

inline void check_polygon(const pd::PolygonWithHoles<double>& p)
{
  using EK = pd::ExactKernel;
  using DK = pd::DoubleKernel;
  if (!pd::validate<DK>(p).ok()) return;  // also rejects non-finite and out-of-range coordinates
  const auto e = pd::convert_polygon<EK>(p);
  if (!pd::validate<EK>(e).ok()) bug("kernels disagree on validity", "");

  const auto r = pd::decompose<EK>(e);
  if (!r.ok()) bug("exact pipeline failed on a valid input", r.diag.message);
  const auto cert = pd::verify_partition<EK>(e, r.leaf_rings(), true);
  if (!cert.ok) bug("exact result not certified", cert.failure);
  const auto cert1 = pd::verify_partition<EK>(e, r.hole_free_rings(), false, r.stats.holes.dangles > 0);
  if (!cert1.ok) bug("hole-free pieces not certified", cert1.failure);

  pd::PipelineOptions opt;
  opt.on_failure = pd::OnFailure::RetryExact;
  const auto d = pd::decompose<DK>(p, opt);
  if (!d.ok()) bug("double pipeline with RetryExact failed", d.diag.message);

  pd::HoleElimOptions hopt;  // hole elimination alone: self-check, then the exact rerun
  hopt.on_failure = pd::OnFailure::RetryExact;
  const auto h = pd::eliminate_holes<DK>(p, hopt);
  if (!h.ok()) bug("double hole elimination with RetryExact failed", h.diag.message);
}

} // namespace pd_fuzz
