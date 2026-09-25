// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#pragma once
/// \file
/// Seeded generators for test and benchmark inputs. The same (family, seed)
/// gives the same polygon on every platform: the random stream is the standard
/// mt19937_64 mapped to [0, 1) by hand, and coordinates are rounded to 1e-5
/// (so last-bit differences between math libraries do not change the result).
/// Every generated polygon with holes is validated before it is returned.

#include <cstdint>
#include <string>
#include <vector>

#include "pd/geometry_types.hpp"

namespace pd::gen {

enum class Family {
  ConvexRow,  ///< convex holes along a line through the centre
  StarHoles,  ///< concave star-shaped holes at random positions
  CShapes,    ///< C-shaped holes whose vertex-average centroid lies outside the ring
  Cluster,    ///< a tight cluster of mixed convex, star and C holes
  Grid,       ///< axis-aligned square holes on a grid (many exact collinearities)
  /// Numerically adversarial families (valid inputs, hard for floating point),
  /// not rounded to 1e-5. NearCollinear and NearTouch use only grid coordinates
  /// and one-ulp nudges (no libm), so they are bit-for-bit reproducible.
  NearCollinear,  ///< edges subdivided by vertices a few ulps off the edge line
  NearTouch,      ///< rectangular holes 2^-30 apart from each other and the boundary
  Snapped,        ///< a mixed cluster snapped to a 1/8 grid (exact collinearities)
};

/// The regular families (ConvexRow .. Grid), in order.
const std::vector<Family>& all_families();
/// The adversarial families (NearCollinear .. Snapped), in order.
const std::vector<Family>& adversarial_families();
const char* name(Family f);
/// Parses a family name as printed by name(); throws std::invalid_argument.
Family family_from_name(const std::string& s);

/// A valid polygon with holes of the given family, inside the square
/// [0, 10] x [0, 10]. Deterministic in (family, seed).
PolygonWithHoles<double> generate(Family f, std::uint64_t seed);

/// Axis-aligned comb with `teeth` teeth of height h (no holes), counter-clockwise,
/// 4 * teeth vertices. Chords from its gap corners hit vertices exactly.
Ring<double> comb(int teeth, double h = 3);

/// A random star-shaped simple polygon with n vertices (angles sorted, radii in
/// [3, 10] around (10, 10)), coordinates rounded to 1e-5. Deterministic in (n, seed).
Ring<double> star_polygon(int n, std::uint64_t seed);

/// h small quadrilateral holes on a jittered spiral inside [-100, 100]^2; the
/// farthest-pair cuts leave most holes on one side, so the recursion is deep.
PolygonWithHoles<double> spiral_holes(int h, std::uint64_t seed);

/// Hand-made degenerate inputs: {name, polygon}. Each forces a fallback or
/// special case of the hole-elimination stage.
std::vector<std::pair<std::string, PolygonWithHoles<double>>> degenerate_set();

} // namespace pd::gen
