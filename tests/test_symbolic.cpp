// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#include "doctest.h"

#ifdef PD_HAVE_GMP
#include "pd/kernel/symbolic.hpp"

using pd::EpsQ;

TEST_CASE("EpsQ: arithmetic and the limit order")
{
  const EpsQ e = EpsQ::eps();
  const EpsQ one(1), two(2);
  CHECK(e > 0);
  CHECK(e < mpq_class(1, 1000000000));  // smaller than every positive rational
  CHECK(-e < 0);
  CHECK(one + e > one);
  CHECK(one - e < one);
  CHECK((one + e) * (one - e) < one);  // 1 - e^2
  CHECK((one + e) * (one - e) > one - e);
  CHECK(e * e < e);
  CHECK(e * e > 0);
  CHECK((two * e) / e == two);  // common powers of e cancel
  CHECK((e / (e * e)) > EpsQ(1000000));  // 1/e is infinitely large
  CHECK((one + e).limit() == 1);
  CHECK(((two * e + e * e) / e).limit() == 2);
  CHECK((one / (one + e)).limit() == 1);
  CHECK(EpsQ(0.1).limit() == mpq_class(0.1));  // doubles convert exactly
  CHECK((e - e).sign() == 0);
  CHECK(e != 0);
  CHECK(EpsQ(3) == 3);
  CHECK_THROWS((one / e).limit());
}

TEST_CASE("EpsQ: values are kept in lowest terms")
{
  const EpsQ e = EpsQ::eps();
  const EpsQ one(1), two(2);
  const EpsQ q = (one + e) / (one + e);
  CHECK(q.is_constant());
  CHECK(q.constant() == 1);
  // (e + 1)(e + 2) / (e + 1) = e + 2: the common factor cancels, so the difference is a constant.
  const EpsQ r = (e * e + 3 * e + two) / (e + one);
  CHECK(r == e + two);
  CHECK((r - e).is_constant());
  CHECK((r - e).constant() == 2);
  // A point on the line through c = (1 + e, 0) and (0, 1), at x = 1/2, computed along two paths.
  const EpsQ cx = one + e, t = (cx - mpq_class(1, 2)) / cx;  // y = t on that line
  const EpsQ y1 = t, y2 = one - mpq_class(1, 2) / cx;
  CHECK(y1 == y2);
  CHECK((y1 - y2).is_constant());
  CHECK((y1 * cx - (cx - mpq_class(1, 2))).is_constant());
}

TEST_CASE("symbolic kernel: orientation with a perturbed point")
{
  using K = pd::SymbolicKernel;
  const EpsQ e = EpsQ::eps();
  const pd::Point<EpsQ> a{0, 0}, b{2, 0};
  const pd::Point<EpsQ> on{1, 0}, up{1, e}, down{1, -e};
  CHECK(K::orient(a, b, on) == 0);
  CHECK(K::orient(a, b, up) > 0);
  CHECK(K::orient(a, b, down) < 0);
  // c = v + e d stays exactly on the line through v in direction d.
  const pd::Point<EpsQ> v{3, 0}, c{3 + 2 * e, 0};
  CHECK(K::orient(a, b, c) == 0);
  CHECK(pd::is_exact_kernel<K>::value);
  CHECK(pd::is_symbolic_kernel<K>::value);
}
#endif
