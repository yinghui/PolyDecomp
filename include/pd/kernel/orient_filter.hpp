// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#pragma once
/// \file
/// Exact sign of the 2D orientation determinant for double coordinates, using a
/// floating-point filter with an exact fallback (Shewchuk, "Adaptive Precision
/// Floating-Point Arithmetic and Fast Robust Geometric Predicates", 1997).
///
/// The fallback expands the determinant into a sum of products, makes each
/// product error-free with a fused multiply-add, and sums the terms exactly as
/// a floating-point expansion. It is exact for finite inputs whose products do
/// not overflow or underflow, which holds for coordinates of magnitude between
/// about 1e-150 and 1e150 (or exactly zero).

#include <cmath>

namespace pd::detail {

/// s + e == a + b exactly, with s = fl(a + b) (Knuth's TwoSum).
inline void two_sum(double a, double b, double& s, double& e)
{
  s = a + b;
  const double bv = s - a;
  e = (a - (s - bv)) + (b - bv);
}

/// p + e == a * b exactly, with p = fl(a * b).
inline void two_product(double a, double b, double& p, double& e)
{
  p = a * b;
  e = std::fma(a, b, -p);
}

/// Adds b to the nonoverlapping expansion e[0..n) (ordered by increasing
/// magnitude) and writes the zero-free result to h, returning its length.
/// Shewchuk's GROW-EXPANSION with zero elimination; h may not alias e.
inline int grow_expansion_zeroelim(int n, const double* e, double b, double* h)
{
  int len = 0;
  double q = b;
  for (int i = 0; i < n; ++i) {
    double s, err;
    two_sum(q, e[i], s, err);
    q = s;
    if (err != 0.0) h[len++] = err;
  }
  if (q != 0.0 || len == 0) h[len++] = q;
  return len;
}

/// Exact sign of the sum of count terms (count <= 16).
inline int exact_sign_of_sum(const double* terms, int count)
{
  double a[17], b[17];
  double* cur = a;
  double* next = b;
  int len = 1;
  cur[0] = terms[0];
  for (int i = 1; i < count; ++i) {
    len = grow_expansion_zeroelim(len, cur, terms[i], next);
    double* t = cur; cur = next; next = t;
  }
  // The largest-magnitude component of a nonoverlapping expansion carries its sign.
  const double top = cur[len - 1];
  return (top > 0.0) - (top < 0.0);
}

/// Exact sign of (b - a) x (c - a).
inline int orient2d_exact(double ax, double ay, double bx, double by, double cx, double cy)
{
  // (bx-ax)(cy-ay) - (by-ay)(cx-ax) expanded; the ax*ay terms cancel.
  const double f[6][2] = {{bx, cy}, {-bx, ay}, {-ax, cy}, {-by, cx}, {by, ax}, {ay, cx}};
  double terms[12];
  for (int i = 0; i < 6; ++i) two_product(f[i][0], f[i][1], terms[2 * i], terms[2 * i + 1]);
  return exact_sign_of_sum(terms, 12);
}

/// Sign of (b - a) x (c - a): +1 left turn (counter-clockwise), -1 right turn,
/// 0 collinear. Exact under the magnitude condition stated above.
inline int orient2d_filtered(double ax, double ay, double bx, double by, double cx, double cy)
{
  constexpr double eps = 1.1102230246251565e-16;              // 2^-53
  constexpr double errbound_a = (3.0 + 16.0 * eps) * eps;     // Shewchuk's ccwerrboundA

  const double detleft = (ax - cx) * (by - cy);
  const double detright = (ay - cy) * (bx - cx);
  const double det = detleft - detright;
  double detsum;
  if (detleft > 0.0) {
    if (detright <= 0.0) return (det > 0.0) - (det < 0.0);
    detsum = detleft + detright;
  } else if (detleft < 0.0) {
    if (detright >= 0.0) return (det > 0.0) - (det < 0.0);
    detsum = -detleft - detright;
  } else {
    // detleft is exactly zero, so det = -detright, whose sign is exact:
    // rounding preserves the signs of the differences and of their product.
    return (det > 0.0) - (det < 0.0);
  }
  const double errbound = errbound_a * detsum;
  if (det >= errbound || -det >= errbound) return (det > 0.0) - (det < 0.0);
  return orient2d_exact(ax, ay, bx, by, cx, cy);
}

} // namespace pd::detail
