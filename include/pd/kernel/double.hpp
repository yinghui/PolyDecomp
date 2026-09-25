// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#pragma once
/// \file
/// Double-precision kernels.
///
/// DoubleKernel stores coordinates as doubles but evaluates the orientation
/// predicate exactly (filtered, see orient_filter.hpp), so every decision made
/// on stored points is consistent. Constructed points (intersections) are
/// rounded to double.
///
/// NaiveDoubleKernel evaluates the predicate in plain floating point. It is kept
/// only to demonstrate and measure the failures that robust predicates avoid.

#include "pd/geometry_types.hpp"
#include "pd/kernel/orient_filter.hpp"

namespace pd {

struct DoubleKernel {
  using FT = double;
  static constexpr const char* name = "double";
  /// A cut line must pass at least clearance * (coordinate scale) from every
  /// vertex; otherwise rounded crossings near a vertex can merge with it.
  static constexpr double clearance = 1e-9;
  static constexpr bool robust = true;  ///< derived predicates are exact on stored points
  /// Convex decomposition treats a vertex as reflex only if it turns by more than
  /// this angle (radians); nearly straight vertices would otherwise produce
  /// needle pieces too thin for double precision. Output pieces are convex up to it.
  static constexpr double angle_tolerance = 1e-12;

  static int orient(const Point<FT>& a, const Point<FT>& b, const Point<FT>& c)
  {
    return detail::orient2d_filtered(a.x, a.y, b.x, b.y, c.x, c.y);
  }
  static FT from_double(double v) { return v; }
  static double to_double(FT v) { return v; }
  /// The exact-sign predicates need products of coordinate differences to
  /// neither overflow nor underflow: every coordinate must be 0 or have
  /// magnitude in [2^-400, 2^400] (about 1e-120 to 1e120).
  static bool coordinate_ok(FT v)
  {
    const double a = v < 0 ? -v : v;
    return a == 0 || (a >= 0x1p-400 && a <= 0x1p400);
  }
};

struct NaiveDoubleKernel {
  using FT = double;
  static constexpr const char* name = "naive-double";
  static constexpr double clearance = 0.0;  ///< unhardened baseline
  static constexpr bool robust = false;     ///< derived predicates in plain floating point
  static constexpr double angle_tolerance = 0.0;
  static bool coordinate_ok(FT v) { return v - v == 0; }  // finite

  static int orient(const Point<FT>& a, const Point<FT>& b, const Point<FT>& c)
  {
    const double det = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    return (det > 0.0) - (det < 0.0);
  }
  static FT from_double(double v) { return v; }
  static double to_double(FT v) { return v; }
};

} // namespace pd
