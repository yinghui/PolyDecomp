// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#pragma once
/// \file
/// Reader for the legacy Wavefront OBJ ring files (Blender exports: one ring per
/// file, vertices in a plane of the 3D file, edges as `l` lines or one `f` face).

#include <string>
#include <vector>

#include "pd/geometry_types.hpp"
#include "pd/io/io_error.hpp"

namespace pd {

/// Which two file coordinates form the plane. The legacy data lies in the
/// file's X-Z plane (file Y is 0): XZ maps (X, Y, Z) to (X, Z).
enum class ObjPlane { XZ, XY };

/// Reads every ring in the file. With `l` or `f` elements, the rings are the
/// cycles of that edge graph (every vertex must have degree 2); without them,
/// all vertices in file order form one ring. Each ring starts at its
/// lowest-numbered vertex and follows file order where it can.
/// Throws ParseError.
std::vector<Ring<double>> read_obj_rings(const std::string& path, ObjPlane plane = ObjPlane::XZ);

/// Reads a file that must contain exactly one ring.
Ring<double> read_obj_ring(const std::string& path, ObjPlane plane = ObjPlane::XZ);

} // namespace pd
