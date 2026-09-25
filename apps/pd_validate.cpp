// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
// pd_validate: checks whether inputs belong to the supported class (simple
// outer ring; simple, pairwise-disjoint, non-nested holes inside it).
//
//   pd_validate [--plane xz|xy] --outer OUTER.obj [--hole HOLE.obj ...]
//   pd_validate [--plane xz|xy] --multi RINGS.obj   (largest ring = outer, rest = holes)
//   pd_validate --wkt FILE.wkt
//
// Exit status: 0 all inputs valid, 1 some input invalid, 2 usage or read error.

#include <cmath>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "pd/geometry.hpp"
#include "pd/io/obj.hpp"
#include "pd/io/wkt.hpp"
#include "pd/kernel/double.hpp"
#include "pd/validate.hpp"
#ifdef PD_HAVE_GMP
#include "pd/kernel/exact.hpp"
#endif

namespace {

int usage()
{
  std::cerr << "usage: pd_validate [--plane xz|xy] --outer OUTER.obj [--hole HOLE.obj ...]\n"
               "       pd_validate [--plane xz|xy] --multi RINGS.obj\n"
               "       pd_validate --wkt FILE.wkt\n";
  return 2;
}

const char* orientation_name(int o) { return o > 0 ? "CCW" : o < 0 ? "CW" : "degenerate"; }

template <class K>
bool report(const char* label, const pd::PolygonWithHoles<typename K::FT>& p)
{
  const pd::ValidationReport rep = pd::validate<K>(p);
  if (rep.ok()) {
    std::cout << "  " << label << ": IN CLASS\n";
    return true;
  }
  std::cout << "  " << label << ": OUT OF CLASS (" << rep.issues.size() << " issue"
            << (rep.issues.size() == 1 ? "" : "s") << ")\n";
  for (const pd::Issue& is : rep.issues) std::cout << "    - " << pd::describe(is) << '\n';
  return false;
}

bool check(const std::string& name, const pd::PolygonWithHoles<double>& p)
{
  std::cout << name << ": outer " << p.outer.size() << " vertices ("
            << orientation_name(pd::ring_orientation<pd::DoubleKernel>(p.outer)) << ")";
  if (!p.holes.empty()) {
    std::cout << ", " << p.holes.size() << " hole" << (p.holes.size() == 1 ? "" : "s") << " [";
    for (std::size_t i = 0; i < p.holes.size(); ++i)
      std::cout << (i ? " " : "") << p.holes[i].size() << ':'
                << orientation_name(pd::ring_orientation<pd::DoubleKernel>(p.holes[i]));
    std::cout << ']';
  }
  std::cout << '\n';
  bool ok = report<pd::DoubleKernel>("double", p);
#ifdef PD_HAVE_GMP
  const bool exact_ok = report<pd::ExactKernel>("exact ", pd::convert_polygon<pd::ExactKernel>(p));
  if (exact_ok != ok) std::cout << "  WARNING: kernels disagree\n";
  ok = exact_ok;
#endif
  return ok;
}

} // namespace

int main(int argc, char** argv)
{
  pd::ObjPlane plane = pd::ObjPlane::XZ;
  std::string outer, wkt, multi;
  std::vector<std::string> holes;
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    const bool has_value = i + 1 < argc;
    if (a == "--plane" && has_value) {
      const std::string v = argv[++i];
      if (v == "xz") plane = pd::ObjPlane::XZ;
      else if (v == "xy") plane = pd::ObjPlane::XY;
      else return usage();
    } else if (a == "--outer" && has_value) {
      outer = argv[++i];
    } else if (a == "--hole" && has_value) {
      holes.push_back(argv[++i]);
    } else if (a == "--wkt" && has_value) {
      wkt = argv[++i];
    } else if (a == "--multi" && has_value) {
      multi = argv[++i];
    } else {
      return usage();
    }
  }
  if (int(!outer.empty()) + int(!wkt.empty()) + int(!multi.empty()) != 1) return usage();

  try {
    bool all_ok = true;
    if (!multi.empty()) {
      // The ring of largest absolute area is the outer ring.
      std::vector<pd::Ring<double>> rings = pd::read_obj_rings(multi, plane);
      std::size_t best = 0;
      for (std::size_t i = 1; i < rings.size(); ++i)
        if (std::abs(pd::signed_area2(rings[i])) > std::abs(pd::signed_area2(rings[best])))
          best = i;
      pd::PolygonWithHoles<double> p;
      p.outer = rings[best];
      for (std::size_t i = 0; i < rings.size(); ++i)
        if (i != best) p.holes.push_back(rings[i]);
      all_ok = check(multi, p);
    } else if (!wkt.empty()) {
      const std::vector<pd::PolygonWithHoles<double>> polys = pd::read_wkt_file(wkt);
      for (std::size_t i = 0; i < polys.size(); ++i)
        all_ok &= check(wkt + "#" + std::to_string(i), polys[i]);
    } else {
      pd::PolygonWithHoles<double> p;
      p.outer = pd::read_obj_ring(outer, plane);
      for (const std::string& h : holes) p.holes.push_back(pd::read_obj_ring(h, plane));
      all_ok = check(outer, p);
    }
    return all_ok ? 0 : 1;
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << '\n';
    return 2;
  }
}
