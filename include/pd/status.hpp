// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#pragma once
/// \file
/// Outcome reporting shared by the algorithms.

#include <string>
#include <vector>

#include "pd/validate.hpp"

namespace pd {

enum class Status {
  Ok,
  InvalidInput,      ///< the input is outside the supported class (see issues)
  NumericalFailure,  ///< a floating-point kernel broke an invariant; retry with ExactKernel
  InternalError,     ///< an invariant failed in exact arithmetic (a bug; please report)
};

inline const char* to_string(Status s)
{
  switch (s) {
  case Status::Ok: return "ok";
  case Status::InvalidInput: return "invalid input";
  case Status::NumericalFailure: return "numerical failure";
  case Status::InternalError: return "internal error";
  }
  return "unknown";
}

struct Diagnostics {
  Status status = Status::Ok;
  std::string message;        ///< human-readable detail when status != Ok
  std::vector<Issue> issues;  ///< validation issues when status == InvalidInput
};

/// True for kernels whose arithmetic is exact, so that a broken invariant is a
/// bug rather than a rounding effect. Specialised by the exact kernel.
template <class K>
struct is_exact_kernel {
  static constexpr bool value = false;
};

} // namespace pd
