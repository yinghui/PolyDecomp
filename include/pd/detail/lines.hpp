// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#pragma once
/// \file
/// Supporting lines and number-type helpers shared by the two stages.

#include <algorithm>
#include <cmath>

#include "pd/geometry_types.hpp"
#ifdef PD_HAVE_GMP
#include <gmpxx.h>
#endif

namespace pd::detail {

template <class FT>
struct Line {
  Point<FT> p, d;  ///< the points p + t d
};

/// Intersection of two non-parallel lines.
template <class FT>
Point<FT> intersect_lines(const Line<FT>& a, const Line<FT>& b)
{
  const FT den = a.d.x * b.d.y - a.d.y * b.d.x;
  const FT wx = b.p.x - a.p.x, wy = b.p.y - a.p.y;
  const FT t = (wx * b.d.y - wy * b.d.x) / den;
  return {a.p.x + t * a.d.x, a.p.y + t * a.d.y};
}

/// Intersection of the line through a0, a1 with the line through b0, b1 (not parallel).
template <class FT>
Point<FT> intersect_through(const Point<FT>& a0, const Point<FT>& a1, const Point<FT>& b0,
                            const Point<FT>& b1)
{
  return intersect_lines(Line<FT>{a0, {a1.x - a0.x, a1.y - a0.y}}, Line<FT>{b0, {b1.x - b0.x, b1.y - b0.y}});
}

/// False only for a floating-point point with a NaN or infinite coordinate.
inline bool ft_finite(double v) { return std::isfinite(v); }
#ifdef PD_HAVE_GMP
inline bool ft_finite(const mpq_class&) { return true; }
#endif
template <class FT>
bool point_finite(const Point<FT>& p)
{
  return ft_finite(p.x) && ft_finite(p.y);
}

/// True if x lies in the bounding box of segment ab enlarged by tol (all in
/// double precision). Used to catch floating-point constructions that land far
/// from the edge they should lie on.
template <class FT>
bool near_box(const Point<FT>& a, const Point<FT>& b, const Point<FT>& x, double tol)
{
  auto d = [](const FT& v) { return static_cast<double>(v); };
  return d(x.x) >= std::min(d(a.x), d(b.x)) - tol && d(x.x) <= std::max(d(a.x), d(b.x)) + tol &&
         d(x.y) >= std::min(d(a.y), d(b.y)) - tol && d(x.y) <= std::max(d(a.y), d(b.y)) + tol;
}

/// Largest integer not above v.
inline double ft_floor(double v) { return std::floor(v); }
#ifdef PD_HAVE_GMP
inline mpq_class ft_floor(const mpq_class& v)
{
  mpz_class z;
  mpz_fdiv_q(z.get_mpz_t(), v.get_num_mpz_t(), v.get_den_mpz_t());
  return mpq_class(z);
}
#endif

/// A dyadic rational j / 2^k strictly between lo < hi, with k as small as
/// possible (so its exact representation stays short). Falls back to the
/// midpoint if the gap is too narrow for the number type to resolve.
template <class FT>
FT dyadic_between(const FT& lo, const FT& hi)
{
  FT scale = 1;
  for (int k = 0; k < 1100; ++k) {
    const FT j = ft_floor(lo * scale) + FT(1);
    const FT s = j / scale;
    if (lo < s && s < hi) return s;
    scale *= 2;
  }
  return (lo + hi) / 2;
}

} // namespace pd::detail
