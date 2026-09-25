// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
// pd_gen: writes seeded test inputs as WKT, one POLYGON per line.
//
//   pd_gen --family NAME|all|adversarial [--count N] [--seed S] [--out FILE]
//     NAME: row, star, c, cluster, grid (regular) or nearcollinear, neartouch,
//     snapped (numerically adversarial); all = the regular families
//   pd_gen --comb TEETH [--out FILE]
//   pd_gen --degenerate [--out FILE]
//
// Case i of a family uses seed S + i, so any case can be regenerated alone.

#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "pd/generate.hpp"
#include "pd/io/wkt.hpp"

namespace {

int usage()
{
  std::cerr << "usage: pd_gen --family NAME|all|adversarial [--count N] [--seed S] [--out FILE]\n"
               "       pd_gen --comb TEETH [--out FILE]\n"
               "       pd_gen --degenerate [--out FILE]\n";
  return 2;
}

} // namespace

int main(int argc, char** argv)
{
  std::string family, out_path;
  long count = 1, teeth = 0;
  std::uint64_t seed = 1;
  bool degenerate = false;
  try {
    for (int i = 1; i < argc; ++i) {
      const std::string a = argv[i];
      const bool v = i + 1 < argc;
      if (a == "--family" && v) family = argv[++i];
      else if (a == "--count" && v) count = std::stol(argv[++i]);
      else if (a == "--seed" && v) seed = std::stoull(argv[++i]);
      else if (a == "--out" && v) out_path = argv[++i];
      else if (a == "--comb" && v) teeth = std::stol(argv[++i]);
      else if (a == "--degenerate") degenerate = true;
      else return usage();
    }
    if (int(!family.empty()) + int(teeth > 0) + int(degenerate) != 1 || count < 1) return usage();

    std::ofstream file;
    if (!out_path.empty()) {
      file.open(out_path);
      if (!file) throw std::runtime_error("cannot write " + out_path);
    }
    std::ostream& os = out_path.empty() ? std::cout : file;

    if (teeth > 0) {
      pd::write_wkt(os, {pd::gen::comb(static_cast<int>(teeth)), {}});
      os << '\n';
    } else if (degenerate) {
      for (const auto& d : pd::gen::degenerate_set()) {
        pd::write_wkt(os, d.second);
        os << '\n';
      }
    } else {
      std::vector<pd::gen::Family> fams;
      if (family == "all") fams = pd::gen::all_families();
      else if (family == "adversarial") fams = pd::gen::adversarial_families();
      else fams.push_back(pd::gen::family_from_name(family));
      for (pd::gen::Family f : fams)
        for (long i = 0; i < count; ++i) {
          pd::write_wkt(os, pd::gen::generate(f, seed + static_cast<std::uint64_t>(i)));
          os << '\n';
        }
    }
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << '\n';
    return 2;
  }
}
