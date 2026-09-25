// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#include "pd/io/leaves.hpp"

#include <fstream>

namespace pd {

std::vector<Ring<double>> read_leaves(const std::string& path)
{
  std::ifstream in(path);
  if (!in) throw ParseError(path + ": cannot open file");
  std::size_t count = 0;
  if (!(in >> count)) throw ParseError(path + ": missing leaf count");
  std::vector<Ring<double>> leaves(count);
  for (std::size_t i = 0; i < count; ++i) {
    std::size_t m = 0;
    if (!(in >> m) || m < 3)
      throw ParseError(path + ": bad vertex count for leaf " + std::to_string(i));
    leaves[i].resize(m);
    for (Point<double>& p : leaves[i])
      if (!(in >> p.x >> p.y))
        throw ParseError(path + ": truncated leaf " + std::to_string(i));
  }
  return leaves;
}

} // namespace pd
