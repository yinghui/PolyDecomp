// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#pragma once
// Shared oracle for the fuzz targets: on every valid input the exact pipeline
// must succeed and pass the exact certificate, and the double pipeline with
// RetryExact must succeed. Any violation aborts (a finding for libFuzzer).

#include <cstdio>
#include <cstdlib>

#include "pd/kernel/double.hpp"
#include "pd/kernel/exact.hpp"
#include "pd/pipeline.hpp"
#include "pd/validate.hpp"
#include "pd/verify.hpp"

namespace pd_fuzz {

[[noreturn]] inline void bug(const char* what, const std::string& detail)
{
  std::fprintf(stderr, "FUZZ FINDING: %s: %s\n", what, detail.c_str());
  std::abort();
}

inline void check_polygon(const pd::PolygonWithHoles<double>& p)
{
  using EK = pd::ExactKernel;
  using DK = pd::DoubleKernel;
  if (!pd::validate<DK>(p).ok()) return;  // also rejects non-finite and out-of-range coordinates
  const auto e = pd::convert_polygon<EK>(p);
  if (!pd::validate<EK>(e).ok()) bug("kernels disagree on validity", "");

  const auto r = pd::decompose<EK>(e);
  if (!r.ok()) bug("exact pipeline failed on a valid input", r.diag.message);
  const auto cert = pd::verify_partition<EK>(e, r.leaf_rings(), true);
  if (!cert.ok) bug("exact result not certified", cert.failure);
  const auto cert1 = pd::verify_partition<EK>(e, r.hole_free_rings(), false);
  if (!cert1.ok) bug("hole-free pieces not certified", cert1.failure);

  pd::PipelineOptions opt;
  opt.on_failure = pd::OnFailure::RetryExact;
  const auto d = pd::decompose<DK>(p, opt);
  if (!d.ok()) bug("double pipeline with RetryExact failed", d.diag.message);
}

} // namespace pd_fuzz
