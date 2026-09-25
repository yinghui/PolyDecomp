// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#pragma once

#include <stdexcept>
#include <string>

namespace pd {

/// Thrown by the readers on malformed or unsupported input.
class ParseError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

} // namespace pd
