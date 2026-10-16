// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
// pd_slice: cuts a polygon with holes (the first polygon of a WKT file) along lines with
// pd::Slicer, in exact arithmetic, and writes the pieces as WKT, one POLYGON per line (holes
// that no line crossed stay holes of their piece). Optionally certifies the pieces and writes the
// convex pieces of every piece.
//
//   pd_slice --wkt IN.wkt (--cut X0 Y0 X1 Y1 | --graze H J | --through H1 J1 H2 J2)...
//            [--out PIECES.wkt] [--convex CONVEX.wkt] [--verify]
//
// Each --cut is the line through (X0, Y0) and (X1, Y1); --graze H J is the line through vertex J
// of hole H of the input parallel to the chord of its two neighbours, computed exactly (for a
// convex vertex, a line that touches the hole only there); --through H1 J1 H2 J2 is the line
// through vertex J1 of ring H1 and vertex J2 of ring H2 (ring -1 is the outer ring, ring k >= 0
// hole k). Cuts are applied in the order given.
// Without --out the pieces go to stdout.
// The exit status is 0 if every cut succeeded (and, with --verify, the certificate passed).
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "pd/pd.hpp"

namespace {

int usage()
{
  std::cerr << "usage: pd_slice --wkt IN.wkt (--cut X0 Y0 X1 Y1 | --graze H J | --through H1 J1 H2 J2)...\n"
               "                [--out PIECES.wkt] [--convex CONVEX.wkt] [--verify]\n";
  return 2;
}

} // namespace

int main(int argc, char** argv)
{
  std::string in_path, out_path, convex_path;
  struct CutArg {
    int kind;                // 0 --cut, 1 --graze, 2 --through
    pd::Point<double> a, b;  // --cut
    long r1, j1, r2, j2;     // rings and vertex indices (--graze uses r1, j1)
  };
  std::vector<CutArg> cuts;
  bool verify = false;
  try {
    for (int i = 1; i < argc; ++i) {
      const std::string a = argv[i];
      if (a == "--wkt" && i + 1 < argc) in_path = argv[++i];
      else if (a == "--out" && i + 1 < argc) out_path = argv[++i];
      else if (a == "--convex" && i + 1 < argc) convex_path = argv[++i];
      else if (a == "--verify") verify = true;
      else if (a == "--cut" && i + 4 < argc) {
        const double x0 = std::stod(argv[i + 1]), y0 = std::stod(argv[i + 2]);
        const double x1 = std::stod(argv[i + 3]), y1 = std::stod(argv[i + 4]);
        cuts.push_back({0, {x0, y0}, {x1, y1}, 0, 0, 0, 0});
        i += 4;
      } else if (a == "--graze" && i + 2 < argc) {
        cuts.push_back({1, {}, {}, std::stol(argv[i + 1]), std::stol(argv[i + 2]), 0, 0});
        i += 2;
      } else if (a == "--through" && i + 4 < argc) {
        cuts.push_back({2, {}, {}, std::stol(argv[i + 1]), std::stol(argv[i + 2]), std::stol(argv[i + 3]),
                        std::stol(argv[i + 4])});
        i += 4;
      } else return usage();
    }
    if (in_path.empty()) return usage();
    const auto shapes = pd::read_wkt_file(in_path);
    if (shapes.empty()) {
      std::cerr << in_path << ": no polygon\n";
      return 1;
    }
    pd::Slicer s(shapes[0]);
    if (!s.ok()) {
      std::cerr << "invalid input: " << s.diag().message << '\n';
      for (const auto& issue : s.diag().issues) std::cerr << "  " << pd::describe(issue) << '\n';
      return 1;
    }
    for (std::size_t k = 0; k < cuts.size(); ++k) {
      pd::SliceResult r;
      const auto& shape = shapes[0];
      auto ring = [&](long h) -> const pd::Ring<double>& {
        if (h < -1 || h >= static_cast<long>(shape.holes.size())) throw std::invalid_argument("no ring " + std::to_string(h));
        return h < 0 ? shape.outer : shape.holes[static_cast<std::size_t>(h)];
      };
      auto vertex = [&](long h, long j) {
        const auto& rg = ring(h);
        const auto& v = rg[static_cast<std::size_t>(j) % rg.size()];
        return pd::Point<mpq_class>{mpq_class(v.x), mpq_class(v.y)};
      };
      if (cuts[k].kind == 0) {
        r = s.cut(cuts[k].a, cuts[k].b);
      } else if (cuts[k].kind == 2) {
        r = s.cut_exact(vertex(cuts[k].r1, cuts[k].j1), vertex(cuts[k].r2, cuts[k].j2));
      } else {
        const auto& h = ring(cuts[k].r1);
        const std::size_t m = h.size(), j = static_cast<std::size_t>(cuts[k].j1) % m;
        const pd::Point<mpq_class> v{mpq_class(h[j].x), mpq_class(h[j].y)};
        const pd::Point<mpq_class> w{v.x + mpq_class(h[(j + 1) % m].x) - mpq_class(h[(j + m - 1) % m].x),
                                     v.y + mpq_class(h[(j + 1) % m].y) - mpq_class(h[(j + m - 1) % m].y)};
        r = s.cut_exact(v, w);
      }
      if (!r.ok()) {
        std::cerr << "cut " << k << " failed: " << r.diag.message << '\n';
        return 1;
      }
      std::cerr << "cut " << k << ": " << r.cut_pieces << " pieces cut, " << s.size() << " pieces, "
                << r.vertex_contacts << " vertex contacts, " << r.edge_contacts << " edge contacts, "
                << r.dangles << " dangling contacts\n";
    }
    int status = 0;
    if (verify) {
      const auto rep = pd::verify_slices<pd::ExactKernel>(pd::convert_polygon<pd::ExactKernel>(shapes[0]),
                                                          s.pieces_exact(), s.symbolic());
      std::cerr << "certificate: " << (rep.ok ? "ok" : "FAILED: " + rep.failure) << '\n';
      if (!rep.ok) status = 1;
    }
    std::ofstream file;
    if (!out_path.empty()) file.open(out_path);
    std::ostream& os = out_path.empty() ? std::cout : file;
    for (std::size_t i = 0; i < s.size(); ++i) {
      pd::write_wkt(os, s.piece(i));
      os << '\n';
    }
    if (!convex_path.empty()) {
      std::ofstream cf(convex_path);
      for (std::size_t i = 0; i < s.size(); ++i) {
        const auto c = s.convex(i);
        if (!c.ok()) {
          std::cerr << "convex decomposition of piece " << i << " failed: " << c.diag.message << '\n';
          return 1;
        }
        if (verify) {
          const auto rep = pd::verify_partition<pd::ExactKernel>(s.piece_exact(i), c.leaf_rings(), true);
          if (!rep.ok) {
            std::cerr << "convex pieces of piece " << i << " not certified: " << rep.failure << '\n';
            status = 1;
          }
        }
        std::vector<pd::Ring<double>> rings;
        for (const auto& r : c.leaf_rings()) rings.push_back(pd::to_double_ring<pd::ExactKernel>(r));
        pd::write_wkt_multipolygon(cf, rings);
        cf << '\n';
      }
    }
    return status;
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << '\n';
    return 1;
  }
}
