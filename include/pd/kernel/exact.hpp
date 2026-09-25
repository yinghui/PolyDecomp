// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#pragma once
/// \file
/// Exact rational kernel (GMP). Every predicate and construction is exact, so
/// the algorithms' theorems hold with no general-position caveat.

#ifndef PD_HAVE_GMP
#error "pd/kernel/exact.hpp requires GMP (configure with PD_WITH_GMP=ON)"
#endif

#include <gmpxx.h>

#include "pd/geometry_types.hpp"
#include "pd/status.hpp"

namespace pd {

struct ExactKernel {
  using FT = mpq_class;
  static constexpr const char* name = "exact";
  /// Minimum distance (relative to the input's coordinate scale) a cut line must
  /// keep from every vertex. Exact arithmetic needs none: only exact incidence matters.
  static constexpr double clearance = 0.0;
  static constexpr bool robust = true;
  static constexpr double angle_tolerance = 0.0;  ///< every reflex vertex is resolved exactly

  static int orient(const Point<FT>& a, const Point<FT>& b, const Point<FT>& c)
  {
    // Reused temporaries avoid an allocation per call.
    thread_local mpq_class l, r, t;
    l = b.x - a.x;
    t = c.y - a.y;
    l *= t;
    r = b.y - a.y;
    t = c.x - a.x;
    r *= t;
    const int s = cmp(l, r);
    return (s > 0) - (s < 0);
  }
  /// Every finite double is a dyadic rational, so this conversion is exact.
  static FT from_double(double v) { return FT(v); }
  static double to_double(const FT& v) { return v.get_d(); }
  static bool coordinate_ok(const FT&) { return true; }  // no range limit
};

template <>
struct is_exact_kernel<ExactKernel> {
  static constexpr bool value = true;
};

} // namespace pd
