// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
// libFuzzer target: bytes -> a polygon with holes on the integer grid (as in
// fuzz_grid.cpp) -> up to eight cuts with pd::Slicer. After every cut the pieces
// must pass the exact certificate (verify_slices), and at the end the convex
// decomposition of every piece must be a certified convex partition of it.
//
// After the polygon, each cut reads a mode byte and four bytes. The modes make
// degenerate cuts frequent: 0, a line through two grid points near the holes;
// 1, a line through an input vertex in a small integer direction; 2, the line
// through two input vertices (often along an edge or through several vertices);
// 3, the line through a hole vertex parallel to the chord of its neighbours
// (tangent at a convex vertex: a dangling contact).
#include <cstddef>
#include <cstdint>
#include <vector>

#include "fuzz_common.hpp"
#include "pd/slice.hpp"

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size)
{
  using EK = pd::ExactKernel;
  using P = pd::Point<double>;
  if (size < 4 || size > 512) return 0;
  pd_fuzz::Reader in{data, size};
  const pd::PolygonWithHoles<double> p = pd_fuzz::grid_polygon(in);
  if (!pd::validate<pd::DoubleKernel>(p).ok()) return 0;
  pd::Slicer s(p);
  if (!s.ok()) pd_fuzz::bug("the slicer rejects a valid polygon", s.diag().message);
  const auto e = pd::convert_polygon<EK>(p);
  std::vector<P> vs(p.outer.begin(), p.outer.end());
  for (const auto& h : p.holes) vs.insert(vs.end(), h.begin(), h.end());

  for (int k = 0; k < 8 && in.i + 5 <= in.n; ++k) {
    const int mode = in.next() % 4;
    const std::uint8_t b1 = in.next(), b2 = in.next(), b3 = in.next(), b4 = in.next();
    P a, b;
    if (mode == 0) {
      a = {24.0 + b1 % 64, 24.0 + b2 % 64};
      b = {24.0 + b3 % 64, 24.0 + b4 % 64};
    } else if (mode == 1) {
      a = vs[b1 % vs.size()];
      b = {a.x + (b2 % 16) - 8, a.y + (b3 % 16) - 8};
    } else if (mode == 2) {
      a = vs[b1 % vs.size()];
      b = vs[(b1 + 1 + b2) % vs.size()];
    } else {
      if (p.holes.empty()) continue;
      const auto& h = p.holes[b1 % p.holes.size()];
      const std::size_t m = h.size(), j = b2 % m;
      a = h[j];
      const P& prev = h[(j + m - 1) % m];
      const P& next = h[(j + 1) % m];
      b = {a.x + next.x - prev.x, a.y + next.y - prev.y};
    }
    if (a == b) continue;
    const pd::SliceResult r = s.cut(a, b);
    if (!r.ok()) pd_fuzz::bug("a cut failed", r.diag.message);
    const auto cert = pd::verify_slices<EK>(e, s.pieces_exact(), s.symbolic());
    if (!cert.ok) pd_fuzz::bug("pieces not certified after a cut", cert.failure);
  }
  for (std::size_t i = 0; i < s.size(); ++i) {
    const auto c = s.convex(i);
    if (!c.ok()) pd_fuzz::bug("the convex decomposition of a piece failed", c.diag.message);
    const auto cert = pd::verify_partition<EK>(s.piece_exact(i), c.leaf_rings(), true);
    if (!cert.ok) pd_fuzz::bug("the convex pieces of a piece are not certified", cert.failure);
  }
  return 0;
}
