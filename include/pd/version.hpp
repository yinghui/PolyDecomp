// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#pragma once

namespace pd {

/// Library version, e.g. "1.0.0".
const char* version();

/// True when the library was built with the exact (GMP) kernel.
constexpr bool has_exact_kernel()
{
#ifdef PD_HAVE_GMP
  return true;
#else
  return false;
#endif
}

} // namespace pd
