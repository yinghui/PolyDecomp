// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#pragma once
/// \file
/// Reader for the legacy `.leaves` regression files written by the old
/// convex-decomposition harness: the leaf count, then one line per leaf,
/// "m x1 y1 ... xm ym".

#include <string>
#include <vector>

#include "pd/geometry_types.hpp"
#include "pd/io/io_error.hpp"

namespace pd {

std::vector<Ring<double>> read_leaves(const std::string& path);

} // namespace pd
