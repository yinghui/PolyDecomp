// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#pragma once
/// \file
/// Exact symbolic kernel: rationals extended by one positive infinitesimal ε.
///
/// A value is a rational function N(ε)/D(ε) with rational coefficients, and
/// comparisons are decided in the limit ε -> 0+, i.e. by the sign of the
/// lowest-order coefficients. Hole elimination uses it for one degenerate case
/// only: a hole that meets a cut line in a single vertex v is split into a
/// polygon with the extra vertex c = v + ε d on the line (see hole_elim.hpp).
/// Every decision is then the limit of decisions on simple polygons, and the
/// result is evaluated at ε = 0 (EpsQ::limit) at the end.

#ifndef PD_HAVE_GMP
#error "pd/kernel/symbolic.hpp requires GMP (configure with PD_WITH_GMP=ON)"
#endif

#include <gmpxx.h>

#include <stdexcept>
#include <utility>
#include <vector>

#include "pd/geometry_types.hpp"
#include "pd/kernel/exact.hpp"
#include "pd/status.hpp"

namespace pd {

/// A rational function of ε with rational coefficients, ordered as ε -> 0+.
class EpsQ {
public:
  using Poly = std::vector<mpq_class>;  ///< coefficient k multiplies ε^k

  EpsQ() : n_{mpq_class(0)}, d_{mpq_class(1)} {}
  EpsQ(int v) : n_{mpq_class(v)}, d_{mpq_class(1)} {}
  EpsQ(long v) : n_{mpq_class(v)}, d_{mpq_class(1)} {}
  EpsQ(double v) : n_{mpq_class(v)}, d_{mpq_class(1)} {}  ///< exact: doubles are dyadic
  EpsQ(const mpq_class& v) : n_{v}, d_{mpq_class(1)} {}

  /// The infinitesimal ε itself.
  static EpsQ eps()
  {
    EpsQ e;
    e.n_ = {mpq_class(0), mpq_class(1)};
    return e;
  }

  /// True if the value does not depend on ε.
  bool is_constant() const { return n_.size() == 1 && d_.size() == 1; }

  /// Value at ε = 0 (the limit ε -> 0+); throws if the limit is infinite.
  mpq_class limit() const
  {
    if (sgn(d_[0]) == 0) throw std::domain_error("EpsQ: infinite limit");
    return n_[0] / d_[0];
  }

  /// Sign of the value for all sufficiently small ε > 0.
  int sign() const { return low_sign(n_) * low_sign(d_); }

  friend EpsQ operator+(const EpsQ& a, const EpsQ& b) { return combine(a, b, +1); }
  friend EpsQ operator-(const EpsQ& a, const EpsQ& b) { return combine(a, b, -1); }
  friend EpsQ operator*(const EpsQ& a, const EpsQ& b)
  {
    EpsQ r;
    r.n_ = mul(a.n_, b.n_);
    r.d_ = mul(a.d_, b.d_);
    r.normalize();
    return r;
  }
  friend EpsQ operator/(const EpsQ& a, const EpsQ& b)
  {
    if (b.sign() == 0) throw std::domain_error("EpsQ: division by zero");
    EpsQ r;
    r.n_ = mul(a.n_, b.d_);
    r.d_ = mul(a.d_, b.n_);
    r.normalize();
    return r;
  }
  EpsQ operator-() const
  {
    EpsQ r = *this;
    for (mpq_class& c : r.n_) c = -c;
    return r;
  }
  EpsQ& operator+=(const EpsQ& b) { return *this = *this + b; }
  EpsQ& operator-=(const EpsQ& b) { return *this = *this - b; }
  EpsQ& operator*=(const EpsQ& b) { return *this = *this * b; }
  EpsQ& operator/=(const EpsQ& b) { return *this = *this / b; }

  friend bool operator==(const EpsQ& a, const EpsQ& b) { return (a - b).sign() == 0; }
  friend bool operator!=(const EpsQ& a, const EpsQ& b) { return !(a == b); }
  friend bool operator<(const EpsQ& a, const EpsQ& b) { return (a - b).sign() < 0; }
  friend bool operator>(const EpsQ& a, const EpsQ& b) { return b < a; }
  friend bool operator<=(const EpsQ& a, const EpsQ& b) { return !(b < a); }
  friend bool operator>=(const EpsQ& a, const EpsQ& b) { return !(a < b); }

private:
  Poly n_, d_;  // d_ is never the zero polynomial; both have at least one coefficient

  static int low_sign(const Poly& p)
  {
    for (const mpq_class& c : p)
      if (sgn(c) != 0) return sgn(c);
    return 0;
  }
  static Poly mul(const Poly& a, const Poly& b)
  {
    Poly r(a.size() + b.size() - 1, mpq_class(0));
    for (std::size_t i = 0; i < a.size(); ++i) {
      if (sgn(a[i]) == 0) continue;
      for (std::size_t j = 0; j < b.size(); ++j) r[i + j] += a[i] * b[j];
    }
    return r;
  }
  static Poly add(const Poly& a, const Poly& b, int s)
  {
    Poly r(std::max(a.size(), b.size()), mpq_class(0));
    for (std::size_t i = 0; i < a.size(); ++i) r[i] += a[i];
    for (std::size_t i = 0; i < b.size(); ++i) r[i] += s > 0 ? b[i] : mpq_class(-b[i]);
    return r;
  }
  static EpsQ combine(const EpsQ& a, const EpsQ& b, int s)
  {
    EpsQ r;
    if (a.d_.size() == 1 && b.d_.size() == 1 && a.d_[0] == b.d_[0]) {  // common denominator
      r.n_ = add(a.n_, b.n_, s);
      r.d_ = a.d_;
    } else {
      r.n_ = add(mul(a.n_, b.d_), mul(b.n_, a.d_), s);
      r.d_ = mul(a.d_, b.d_);
    }
    r.normalize();
    return r;
  }
  /// Strips high-order zeros, cancels a common power of ε, and makes a constant
  /// denominator equal to 1, so that constants stay constants.
  void normalize()
  {
    auto strip = [](Poly& p) {
      while (p.size() > 1 && sgn(p.back()) == 0) p.pop_back();
    };
    strip(n_);
    strip(d_);
    if (low_sign(n_) == 0) {  // zero
      n_ = {mpq_class(0)};
      d_ = {mpq_class(1)};
      return;
    }
    std::size_t k = 0;
    while (k < n_.size() && k < d_.size() && sgn(n_[k]) == 0 && sgn(d_[k]) == 0) ++k;
    if (k > 0) {
      n_.erase(n_.begin(), n_.begin() + static_cast<std::ptrdiff_t>(k));
      d_.erase(d_.begin(), d_.begin() + static_cast<std::ptrdiff_t>(k));
    }
    if (d_.size() == 1 && d_[0] != 1) {
      const mpq_class inv = 1 / d_[0];
      for (mpq_class& c : n_) c *= inv;
      d_ = {mpq_class(1)};
    }
  }
};

/// |v| in the limit order (used by generic code paths; exact).
inline EpsQ abs(const EpsQ& v) { return v.sign() < 0 ? -v : v; }
/// Helpers found by argument-dependent lookup from pd::detail (see detail/lines.hpp):
/// every value is finite, and the floor of a value is that of its limit (only used
/// to pick simple dyadic directions, which are then checked exactly).
inline bool ft_finite(const EpsQ&) { return true; }
inline EpsQ ft_floor(const EpsQ& v)
{
  const mpq_class l = v.limit();
  mpz_class z;
  mpz_fdiv_q(z.get_mpz_t(), l.get_num_mpz_t(), l.get_den_mpz_t());
  return EpsQ(mpq_class(z));
}

/// The exact kernel over EpsQ. Only used, automatically, when an exact run
/// meets a dangling contact (see eliminate_holes and decompose).
struct SymbolicKernel {
  using FT = EpsQ;
  static constexpr const char* name = "symbolic";
  static constexpr double clearance = 0.0;
  static constexpr bool robust = true;
  static constexpr double angle_tolerance = 0.0;

  static int orient(const Point<FT>& a, const Point<FT>& b, const Point<FT>& c)
  {
    return ((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)).sign();
  }
  static FT from_double(double v) { return FT(v); }
  static double to_double(const FT& v) { return v.limit().get_d(); }
  static bool coordinate_ok(const FT&) { return true; }
};

template <>
struct is_exact_kernel<SymbolicKernel> {
  static constexpr bool value = true;
};
template <>
struct is_symbolic_kernel<SymbolicKernel> {
  static constexpr bool value = true;
};

/// Exact point at ε = 0.
inline Point<mpq_class> limit_point(const Point<EpsQ>& p) { return {p.x.limit(), p.y.limit()}; }

} // namespace pd
