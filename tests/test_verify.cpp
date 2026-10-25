// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#include "doctest.h"

#ifdef PD_HAVE_GMP

#include "pd/kernel/exact.hpp"
#include "pd/verify.hpp"

namespace {

using K = pd::ExactKernel;
using FT = K::FT;

pd::Ring<FT> R(std::initializer_list<std::pair<double, double>> pts)
{
  pd::Ring<double> r;
  for (const auto& p : pts) r.push_back({p.first, p.second});
  return pd::convert_ring<K>(r);
}

pd::VerifyReport verify(const pd::Ring<FT>& outer, std::vector<pd::Ring<FT>> holes,
                        std::vector<pd::Ring<FT>> pieces, bool convex = false)
{
  return pd::verify_partition<K>({outer, std::move(holes)}, pieces, convex);
}

const auto kSquare = R({{0, 0}, {10, 0}, {10, 10}, {0, 10}});

} // namespace

TEST_CASE("verify accepts valid partitions")
{
  // Two halves.
  CHECK(verify(kSquare, {}, {R({{0, 0}, {5, 0}, {5, 10}, {0, 10}}),
                             R({{5, 0}, {10, 0}, {10, 10}, {5, 10}})}).ok);
  // T-junction: (5, 5) is a vertex of the right pieces but not of the left one.
  CHECK(verify(kSquare, {}, {R({{0, 0}, {5, 0}, {5, 10}, {0, 10}}),
                             R({{5, 0}, {10, 0}, {10, 5}, {5, 5}}),
                             R({{5, 5}, {10, 5}, {10, 10}, {5, 10}})}, true).ok);
  // Square with a square hole, cut into four pieces around it.
  const auto hole = R({{4, 4}, {6, 4}, {6, 6}, {4, 6}});
  CHECK(verify(kSquare, {hole},
               {R({{0, 0}, {10, 0}, {10, 4}, {0, 4}}), R({{0, 6}, {10, 6}, {10, 10}, {0, 10}}),
                R({{0, 4}, {4, 4}, {4, 6}, {0, 6}}), R({{6, 4}, {10, 4}, {10, 6}, {6, 6}})},
               true)
            .ok);
  // Input orientation does not matter.
  pd::Ring<FT> cw(kSquare.rbegin(), kSquare.rend());
  CHECK(verify(cw, {}, {R({{0, 0}, {10, 0}, {10, 10}, {0, 10}})}).ok);
}

TEST_CASE("verify rejects overlaps, gaps, spill-over and bad pieces")
{
  const auto left = R({{0, 0}, {5, 0}, {5, 10}, {0, 10}});
  const auto right = R({{5, 0}, {10, 0}, {10, 10}, {5, 10}});
  CHECK_FALSE(verify(kSquare, {}, {left, left, right}).ok);  // overlap
  CHECK_FALSE(verify(kSquare, {}, {left}).ok);               // gap
  // Same total area as the L-shaped input, but one square lies in the notch.
  const auto ell = R({{0, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 2}, {0, 2}});
  CHECK_FALSE(verify(ell, {}, {R({{0, 0}, {1, 0}, {1, 1}, {0, 1}}),
                               R({{1, 0}, {2, 0}, {2, 1}, {1, 1}}),
                               R({{1, 1}, {2, 1}, {2, 2}, {1, 2}})})
                   .ok);
  CHECK(verify(ell, {}, {R({{0, 0}, {1, 0}, {1, 1}, {0, 1}}), R({{1, 0}, {2, 0}, {2, 1}, {1, 1}}),
                         R({{0, 1}, {1, 1}, {1, 2}, {0, 2}})})
            .ok);
  // A clockwise piece.
  pd::Ring<FT> left_cw(left.rbegin(), left.rend());
  CHECK_FALSE(verify(kSquare, {}, {left_cw, right}).ok);
  // A non-convex piece: fine unless convexity is required.
  CHECK(verify(ell, {}, {ell}).ok);
  const auto rep = verify(ell, {}, {ell}, true);
  CHECK_FALSE(rep.ok);
  CHECK(rep.failure.find("not convex") != std::string::npos);
  // A piece that ignores the hole.
  const auto hole = R({{4, 4}, {6, 4}, {6, 6}, {4, 6}});
  CHECK_FALSE(verify(kSquare, {hole}, {kSquare}).ok);
}


TEST_CASE("verify_slices matches retained holes up to cyclic rotation and reversal")
{
  const auto hole = R({{3, 3}, {7, 3}, {7, 7}, {5, 5}, {3, 7}});
  const pd::PolygonWithHoles<FT> input{kSquare, {hole}};
  for (std::size_t shift = 0; shift < hole.size(); ++shift) {
    auto moved = hole;
    std::rotate(moved.begin(), moved.begin() + static_cast<std::ptrdiff_t>(shift), moved.end());
    CHECK(pd::verify_slices<K>(input, {{kSquare, {moved}}}).ok);
    std::reverse(moved.begin(), moved.end());
    CHECK(pd::verify_slices<K>(input, {{kSquare, {moved}}}).ok);
  }
}

TEST_CASE("verify_slices rejects a self-crossing permutation of a retained hole")
{
  const auto hole = R({{3, 3}, {6, 3}, {6, 6}, {3, 6}});
  const auto bow_tie = R({{3, 3}, {6, 6}, {6, 3}, {3, 6}});
  const pd::PolygonWithHoles<FT> input{kSquare, {hole}};
  CHECK_FALSE(pd::verify_slices<K>(input, {{kSquare, {bow_tie}}}).ok);
}

TEST_CASE("verify_slices rejects a different simple hole with the same vertices and area")
{
  const auto hole = R({{3, 3}, {7, 3}, {7, 7}, {5, 5}, {3, 7}});
  const auto changed = R({{3, 3}, {5, 5}, {7, 3}, {7, 7}, {3, 7}});
  REQUIRE(pd::conflicting_edge_pairs<K>({changed, {}}, true, 1).empty());
  REQUIRE(pd::signed_area2(hole) == pd::signed_area2(changed));
  const pd::PolygonWithHoles<FT> input{kSquare, {hole}};
  CHECK_FALSE(pd::verify_slices<K>(input, {{kSquare, {changed}}}).ok);
}

#endif // PD_HAVE_GMP
