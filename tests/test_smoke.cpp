// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#include "doctest.h"
#include "pd/version.hpp"

#include <string>

#ifdef PD_HAVE_GMP
#include <gmpxx.h>
#endif

TEST_CASE("library reports a version")
{
  CHECK(std::string(pd::version()) == PD_VERSION_STRING);
  CHECK(std::string(pd::version()).find('.') != std::string::npos);
}

#ifdef PD_HAVE_GMP
TEST_CASE("GMP rationals are exact")
{
  CHECK(pd::has_exact_kernel());
  mpq_class a(1, 3), b(1, 6);
  CHECK(a + b == mpq_class(1, 2));

  // 0.1 is not representable in binary; converting the double to mpq keeps
  // its exact binary value, which is what the exact kernel relies on.
  mpq_class tenth(0.1);
  CHECK(tenth != mpq_class(1, 10));
  CHECK(tenth.get_d() == 0.1);
}
#endif
