// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#pragma once
/// \file
/// Exact text format for rational coordinates (requires GMP).
///
/// \verbatim
///   polydecomp-rational 1 <polygon|rings>
///   <number of rings>
///   <n>                     then n lines "x y", each a GMP rational "p/q" or "p"
///   ...
/// \endverbatim
///
/// Kind "polygon" is a polygon with holes (first ring outer, the rest holes);
/// kind "rings" is a list of independent rings (e.g. decomposition leaves).

#ifndef PD_HAVE_GMP
#error "pd/io/rational.hpp requires GMP (configure with PD_WITH_GMP=ON)"
#endif

#include <gmpxx.h>

#include <istream>
#include <ostream>
#include <string>
#include <vector>

#include "pd/geometry_types.hpp"
#include "pd/io/io_error.hpp"

namespace pd {

inline void write_rational_rings(std::ostream& os, const std::vector<Ring<mpq_class>>& rings,
                                 const char* kind = "rings")
{
  os << "polydecomp-rational 1 " << kind << '\n' << rings.size() << '\n';
  for (const Ring<mpq_class>& r : rings) {
    os << r.size() << '\n';
    for (const Point<mpq_class>& p : r) os << p.x.get_str() << ' ' << p.y.get_str() << '\n';
  }
}

inline void write_rational(std::ostream& os, const PolygonWithHoles<mpq_class>& p)
{
  std::vector<Ring<mpq_class>> rings{p.outer};
  rings.insert(rings.end(), p.holes.begin(), p.holes.end());
  write_rational_rings(os, rings, "polygon");
}

/// Reads either kind; returns the rings and sets kind to "polygon" or "rings".
inline std::vector<Ring<mpq_class>> read_rational_rings(std::istream& is, std::string& kind)
{
  std::string magic;
  int version = 0;
  if (!(is >> magic >> version >> kind) || magic != "polydecomp-rational" || version != 1 ||
      (kind != "polygon" && kind != "rings"))
    throw ParseError("rational: bad header");
  std::size_t nrings = 0;
  if (!(is >> nrings)) throw ParseError("rational: missing ring count");
  std::vector<Ring<mpq_class>> rings(nrings);
  for (Ring<mpq_class>& r : rings) {
    std::size_t n = 0;
    if (!(is >> n)) throw ParseError("rational: missing vertex count");
    r.resize(n);
    for (Point<mpq_class>& p : r) {
      std::string sx, sy;
      if (!(is >> sx >> sy)) throw ParseError("rational: truncated ring");
      // set_str signals failure by its return value; a zero denominator
      // must be rejected before canonicalize() would divide by it.
      if (p.x.set_str(sx, 10) != 0 || p.y.set_str(sy, 10) != 0 ||
          sgn(p.x.get_den()) == 0 || sgn(p.y.get_den()) == 0)
        throw ParseError("rational: bad number '" + sx + "' or '" + sy + "'");
      p.x.canonicalize();
      p.y.canonicalize();
    }
  }
  return rings;
}

inline PolygonWithHoles<mpq_class> read_rational(std::istream& is)
{
  std::string kind;
  std::vector<Ring<mpq_class>> rings = read_rational_rings(is, kind);
  if (kind != "polygon" || rings.empty()) throw ParseError("rational: not a polygon");
  PolygonWithHoles<mpq_class> p;
  p.outer = std::move(rings.front());
  p.holes.assign(std::make_move_iterator(rings.begin() + 1), std::make_move_iterator(rings.end()));
  return p;
}

} // namespace pd
