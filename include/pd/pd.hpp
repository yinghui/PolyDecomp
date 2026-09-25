// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#pragma once
/// \file
/// Umbrella header for polydecomp.

#include "pd/version.hpp"
#include "pd/geometry_types.hpp"
#include "pd/geometry.hpp"
#include "pd/validate.hpp"
#include "pd/status.hpp"
#include "pd/convex_decomp.hpp"
#include "pd/hole_elim.hpp"
#include "pd/pipeline.hpp"
#include "pd/generate.hpp"
#include "pd/verify.hpp"
#include "pd/kernel/double.hpp"
#include "pd/io/wkt.hpp"
#include "pd/io/obj.hpp"
#include "pd/io/leaves.hpp"
#ifdef PD_HAVE_GMP
#include "pd/kernel/exact.hpp"
#include "pd/io/rational.hpp"
#endif
