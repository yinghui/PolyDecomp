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
#include <cstddef>
#include <cstdint>

#include "fuzz_common.hpp"

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size)
{
  if (size < 4 || size > 512) return 0;
  pd_fuzz::Reader in{data, size};
  pd_fuzz::check_polygon(pd_fuzz::grid_polygon(in));
  return 0;
}
