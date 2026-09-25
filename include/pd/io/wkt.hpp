// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#pragma once
/// \file
/// Well-known text (WKT) for polygons with holes: POLYGON and MULTIPOLYGON.
/// Coordinates are read correctly rounded and written in the shortest form that
/// reads back to the same double, so a write/read round trip is bit-exact.

#include <iosfwd>
#include <string>
#include <vector>

#include "pd/geometry_types.hpp"
#include "pd/io/io_error.hpp"

namespace pd {

/// Parses a sequence of POLYGON / MULTIPOLYGON geometries (separated by
/// whitespace); each MULTIPOLYGON contributes one entry per member polygon.
/// Keywords are case-insensitive; a third (Z) coordinate is ignored; the
/// repeated closing vertex of each ring is dropped. Throws ParseError.
std::vector<PolygonWithHoles<double>> parse_wkt(const std::string& text);
std::vector<PolygonWithHoles<double>> read_wkt_file(const std::string& path);

void write_wkt(std::ostream& os, const PolygonWithHoles<double>& p);
/// Writes the rings as the members of one MULTIPOLYGON (each without holes).
void write_wkt_multipolygon(std::ostream& os, const std::vector<Ring<double>>& rings);

/// Shortest decimal representation of v that reads back to exactly v.
std::string format_double(double v);

} // namespace pd
