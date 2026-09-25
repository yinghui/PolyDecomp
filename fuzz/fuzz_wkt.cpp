// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
// libFuzzer target: arbitrary text -> WKT parser -> validation -> both pipelines.
#include <cstddef>
#include <cstdint>
#include <string>

#include "fuzz_common.hpp"
#include "pd/io/wkt.hpp"

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size)
{
  if (size > 4096) return 0;
  std::vector<pd::PolygonWithHoles<double>> polys;
  try {
    polys = pd::parse_wkt(std::string(reinterpret_cast<const char*>(data), size));
  } catch (const pd::ParseError&) {
    return 0;
  }
  for (const auto& p : polys) {
    std::size_t n = p.outer.size();
    for (const auto& h : p.holes) n += h.size();
    if (n <= 80 && p.holes.size() <= 6) pd_fuzz::check_polygon(p);
  }
  return 0;
}
