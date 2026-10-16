// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
//
// Example: slice a cheese with 64 holes along a few lines, as a player would in a
// game, including a line through a vertex made by an earlier cut. Every cut is
// exact; the pieces are certified and each is decomposed into convex pieces.
#include <iostream>

#include "pd/pd.hpp"

int main()
{
  using P = pd::Point<double>;
  const pd::PolygonWithHoles<double> cheese = pd::gen::cheese(64, 1);
  pd::Slicer s(cheese);
  if (!s.ok()) {
    std::cerr << "invalid shape: " << s.diag().message << '\n';
    return 1;
  }
  for (const auto& [a, b] : {std::pair{P{-9, 0.5}, P{9, -1}}, std::pair{P{0, -9}, P{0.2, 9}},
                             std::pair{P{-9, -9}, P{9, 9}}}) {
    const pd::SliceResult r = s.cut(a, b);  // every piece the line crosses
    if (!r.ok()) {
      std::cerr << "cut failed: " << r.diag.message << '\n';
      return 1;
    }
    std::cout << "cut " << r.cut_pieces << " pieces into " << r.pieces.size() << "; " << s.size()
              << " pieces now\n";
  }
  // A cut of piece 0 only, through one of its vertices (exact coordinates).
  const pd::PolygonWithHoles<mpq_class> p0 = s.piece_exact(0);
  const pd::Point<mpq_class> v = p0.outer[0];
  const pd::SliceResult r = s.cut_exact(v, {v.x + 1, v.y + 3}, 0);
  std::cout << "through a vertex of piece 0: " << r.vertex_contacts << " vertex contacts; " << s.size()
            << " pieces now\n";

  // The pieces still contain the holes no line crossed; certify them exactly.
  using K = pd::ExactKernel;
  const auto cert = pd::verify_slices<K>(pd::convert_polygon<K>(cheese), s.pieces_exact(), s.symbolic());
  std::size_t convex = 0;
  for (std::size_t i = 0; i < s.size(); ++i) {
    const auto c = s.convex(i);  // convex pieces of piece i, for example for a physics engine
    if (!c.ok()) return 1;
    convex += c.leaves.size();
  }
  std::cout << s.size() << " pieces, certified: " << (cert.ok ? "yes" : "no") << "; " << convex
            << " convex pieces; largest coordinate " << s.max_bits() << " bits\n";
  return cert.ok ? 0 : 1;
}
