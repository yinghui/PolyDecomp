// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#include "pd/validate.hpp"

namespace pd {

namespace {

std::string ring_name(int k)
{
  return k == 0 ? std::string("outer ring") : "hole " + std::to_string(k - 1);
}

} // namespace

std::string describe(const Issue& issue)
{
  switch (issue.kind) {
  case IssueKind::TooFewVertices:
    return ring_name(issue.ring_a) + " has fewer than 3 vertices";
  case IssueKind::RepeatedVertex:
    return ring_name(issue.ring_a) + " repeats vertex " + std::to_string(issue.edge_a) +
           " (zero-length edge)";
  case IssueKind::EdgesIntersect:
    return "edge " + std::to_string(issue.edge_a) + " of " + ring_name(issue.ring_a) +
           " touches edge " + std::to_string(issue.edge_b) + " of " + ring_name(issue.ring_b);
  case IssueKind::HoleNotInsideOuter:
    return ring_name(issue.ring_a) + " is not inside the outer ring";
  case IssueKind::HoleInsideHole:
    return ring_name(issue.ring_a) + " lies inside " + ring_name(issue.ring_b);
  case IssueKind::CoordinateOutOfRange:
    return "vertex " + std::to_string(issue.edge_a) + " of " + ring_name(issue.ring_a) +
           " is outside the kernel's coordinate range";
  }
  return "unknown issue";
}

} // namespace pd
