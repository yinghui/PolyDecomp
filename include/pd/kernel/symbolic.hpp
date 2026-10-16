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

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include "pd/geometry_types.hpp"
#include "pd/kernel/exact.hpp"
#include "pd/status.hpp"

namespace pd {

/// A rational function of ε with rational coefficients, ordered as ε -> 0+.
///
/// Almost every value of a computation does not depend on ε (only a split vertex and what is built
/// from it does). A constant is therefore stored as a plain rational, and only a value that depends
/// on ε carries its polynomials N and D, shared between copies since values are never modified.
class EpsQ {
public:
  using Poly = std::vector<mpq_class>;  ///< coefficient k multiplies ε^k

  EpsQ() = default;
  EpsQ(int v) : c_(v) {}
  EpsQ(long v) : c_(v) {}
  EpsQ(double v) : c_(v) {}  ///< exact: doubles are dyadic
  EpsQ(const mpq_class& v) : c_(v) {}
  EpsQ(mpq_class&& v) : c_(std::move(v)) {}

  /// The infinitesimal ε itself.
  static EpsQ eps() { return EpsQ(Poly{mpq_class(0), mpq_class(1)}, Poly{mpq_class(1)}); }

  /// True if the value does not depend on ε.
  bool is_constant() const { return !f_; }
  /// The value of a constant (is_constant()).
  const mpq_class& constant() const { return c_; }

  /// Value at ε = 0 (the limit ε -> 0+); throws if the limit is infinite.
  mpq_class limit() const
  {
    if (!f_) return c_;
    if (sgn(f_->d[0]) == 0) throw std::domain_error("EpsQ: infinite limit");
    return f_->n[0] / f_->d[0];
  }

  /// Sign of the value for all sufficiently small ε > 0.
  int sign() const { return f_ ? low_sign(f_->n) * low_sign(f_->d) : sgn(c_); }

  // Constants (the common case) take fast paths on their rational values.
  friend EpsQ operator+(const EpsQ& a, const EpsQ& b)
  {
    if (!a.f_ && !b.f_) return EpsQ(mpq_class(a.c_ + b.c_));
    return combine(a, b, +1);
  }
  friend EpsQ operator-(const EpsQ& a, const EpsQ& b)
  {
    if (!a.f_ && !b.f_) return EpsQ(mpq_class(a.c_ - b.c_));
    return combine(a, b, -1);
  }
  friend EpsQ operator*(const EpsQ& a, const EpsQ& b)
  {
    if (!a.f_ && !b.f_) return EpsQ(mpq_class(a.c_ * b.c_));
    return EpsQ(mul(a.num(), b.num()), mul(a.den(), b.den()));
  }
  friend EpsQ operator/(const EpsQ& a, const EpsQ& b)
  {
    if (b.sign() == 0) throw std::domain_error("EpsQ: division by zero");
    if (!a.f_ && !b.f_) return EpsQ(mpq_class(a.c_ / b.c_));
    return EpsQ(mul(a.num(), b.den()), mul(a.den(), b.num()));
  }
  EpsQ operator-() const
  {
    if (!f_) return EpsQ(mpq_class(-c_));
    Poly n = f_->n;
    for (mpq_class& c : n) c = -c;
    return EpsQ(std::move(n), f_->d);
  }
  EpsQ& operator+=(const EpsQ& b) { return *this = *this + b; }
  EpsQ& operator-=(const EpsQ& b) { return *this = *this - b; }
  EpsQ& operator*=(const EpsQ& b) { return *this = *this * b; }
  EpsQ& operator/=(const EpsQ& b) { return *this = *this / b; }

  friend bool operator==(const EpsQ& a, const EpsQ& b)
  {
    if (!a.f_ && !b.f_) return a.c_ == b.c_;
    return (a - b).sign() == 0;
  }
  friend bool operator!=(const EpsQ& a, const EpsQ& b) { return !(a == b); }
  friend bool operator<(const EpsQ& a, const EpsQ& b)
  {
    if (!a.f_ && !b.f_) return a.c_ < b.c_;
    return (a - b).sign() < 0;
  }
  friend bool operator>(const EpsQ& a, const EpsQ& b) { return b < a; }
  friend bool operator<=(const EpsQ& a, const EpsQ& b) { return !(b < a); }
  friend bool operator>=(const EpsQ& a, const EpsQ& b) { return !(a < b); }

private:
  struct Fn {
    Poly n, d;  // N and D: D is never the zero polynomial, and N/D is never a constant
  };
  mpq_class c_;                  // the value of a constant
  std::shared_ptr<const Fn> f_;  // N/D of a value that depends on ε, empty for a constant

  /// The value N/D, normalized; a constant is stored as a constant.
  EpsQ(Poly n, Poly d) { normalize(std::move(n), std::move(d)); }
  Poly num() const { return f_ ? f_->n : Poly{c_}; }
  Poly den() const { return f_ ? f_->d : Poly{mpq_class(1)}; }

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
    Poly ad = a.den(), bd = b.den();
    if (ad.size() == 1 && bd.size() == 1 && ad[0] == bd[0])  // common denominator
      return EpsQ(add(a.num(), b.num(), s), std::move(ad));
    return EpsQ(add(mul(a.num(), bd), mul(b.num(), ad), s), mul(ad, bd));
  }
  /// Removes high-order zero coefficients (the zero polynomial becomes empty).
  static void trim(Poly& p)
  {
    while (!p.empty() && sgn(p.back()) == 0) p.pop_back();
  }
  /// Remainder of a divided by b (b trimmed and not zero); trimmed.
  static Poly remainder(Poly a, const Poly& b)
  {
    trim(a);
    while (a.size() >= b.size()) {
      const mpq_class f = a.back() / b.back();
      const std::size_t s = a.size() - b.size();
      for (std::size_t j = 0; j < b.size(); ++j) a[s + j] -= f * b[j];
      a.pop_back();  // the leading coefficient is now zero
      trim(a);
    }
    return a;
  }
  /// The quotient a / b of a polynomial a that b divides (b trimmed and not zero).
  static Poly quotient(Poly a, const Poly& b)
  {
    Poly q(a.size() - b.size() + 1, mpq_class(0));
    for (std::size_t i = q.size(); i-- > 0;) {
      q[i] = a[i + b.size() - 1] / b.back();
      for (std::size_t j = 0; j < b.size(); ++j) a[i + j] -= q[i] * b[j];
    }
    return q;
  }
  /// A greatest common divisor of two nonzero trimmed polynomials (Euclid).
  static Poly gcd(Poly a, Poly b)
  {
    while (!b.empty()) {
      Poly r = remainder(std::move(a), b);
      a = std::move(b);
      b = std::move(r);
    }
    return a;
  }
  /// Brings N/D to lowest terms (common powers of ε and the polynomial gcd cancelled, the lowest
  /// nonzero coefficient of D equal to 1), so that, as with mpq_class, a value has one representation
  /// whichever way it was computed; a value that does not depend on ε becomes a constant.
  void normalize(Poly n, Poly d)
  {
    trim(n);
    trim(d);
    if (n.empty()) return;  // zero: c_ is 0
    std::size_t k = 0;
    while (sgn(n[k]) == 0 && sgn(d[k]) == 0) ++k;
    if (k > 0) {
      n.erase(n.begin(), n.begin() + static_cast<std::ptrdiff_t>(k));
      d.erase(d.begin(), d.begin() + static_cast<std::ptrdiff_t>(k));
    }
    if (n.size() > 1 && d.size() > 1) {
      const Poly g = gcd(n, d);
      if (g.size() > 1) {
        n = quotient(std::move(n), g);
        d = quotient(std::move(d), g);
      }
    }
    std::size_t low = 0;
    while (sgn(d[low]) == 0) ++low;
    if (d[low] != 1) {
      const mpq_class inv = 1 / d[low];
      for (mpq_class& c : n) c *= inv;
      for (mpq_class& c : d) c *= inv;
    }
    if (n.size() == 1 && d.size() == 1) {
      c_ = std::move(n[0]);
      return;
    }
    f_ = std::make_shared<const Fn>(Fn{std::move(n), std::move(d)});
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
    if (a.x.is_constant() && a.y.is_constant() && b.x.is_constant() && b.y.is_constant() &&
        c.x.is_constant() && c.y.is_constant())
      return ExactKernel::orient(a.x.constant(), a.y.constant(), b.x.constant(), b.y.constant(),
                                 c.x.constant(), c.y.constant());
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
