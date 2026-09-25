// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
// libFuzzer target: bytes -> star-shaped rings on a small integer grid (many
// collinear and nearly touching features) -> validation -> both pipelines.
//
// Byte 0: number of holes (0..4). Each ring then reads: a centre byte
// (x = low 4 bits, y = high 4 bits, on a 4x grid), a vertex-count byte
// (3..14), and one byte per vertex (angle step = low 4 bits + 1, radius = high
// 4 bits + 1). Angles accumulate in units of 1/64 turn; the outer ring is
// scaled up so that holes can fit inside it.
#include <cmath>
#include <cstddef>
#include <cstdint>

#include "fuzz_common.hpp"

namespace {

struct Reader {
  const std::uint8_t* d;
  std::size_t n, i = 0;
  bool more() const { return i < n; }
  std::uint8_t next() { return i < n ? d[i++] : 0; }
};

pd::Ring<double> ring(Reader& in, double scale, double cx0, double cy0)
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

} // namespace

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size)
{
  if (size < 4 || size > 512) return 0;
  Reader in{data, size};
  const int holes = in.next() % 5;
  pd::PolygonWithHoles<double> p;
  p.outer = ring(in, 16, 64, 64);
  for (int h = 0; h < holes && in.more(); ++h) p.holes.push_back(ring(in, 2, 40, 40));
  pd_fuzz::check_polygon(p);
  return 0;
}
