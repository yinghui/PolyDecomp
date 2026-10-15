// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
// pd_decompose: decomposes polygons with holes.
//
//   pd_decompose [options] --wkt FILE                 (every POLYGON in the file)
//   pd_decompose [options] --obj OUTER.obj [--hole H.obj ...] [--plane xz|xy]
//   pd_decompose [options] --multi RINGS.obj [--plane xz|xy]
//
// Options:
//   --stage both|holes|convex   pipeline (default), hole elimination only, or
//                               convex decomposition only (input without holes)
//   --kernel exact|double       arithmetic (default: exact; double if built without GMP)
//   --check                     extra invariant checks inside the algorithms
//   --verify                    certify each result with the exact partition test
//   --out FILE                  pieces as WKT, one MULTIPOLYGON line per input (default stdout)
//   --rational FILE             pieces with exact rational coordinates (exact kernel)
//   --stats FILE                one JSON object per input with counts and timing
//   --retry-exact               with --kernel double: answer a numerical failure by an
//                               exact rerun (pipeline only)
//   --cut-rule farthest|closest|boundary|horizontal
//                               choice of the hole-elimination cut (default farthest;
//                               the others are for comparison, see pd::CutRule)
//   --quiet                     no summary on stderr
//   --compare double|naive      run the exact pipeline and a floating-point one on
//                               each input and classify the latter as agree /
//                               deviate / reported_failure (stats in --stats)
//
// Exit status: 0 every input decomposed (and certified, with --verify),
// 1 some input invalid, failed or not certified, 2 usage or read error.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "pd/convex_decomp.hpp"
#include "pd/hole_elim.hpp"
#include "pd/io/obj.hpp"
#include "pd/io/wkt.hpp"
#include "pd/kernel/double.hpp"
#include "pd/pipeline.hpp"
#include "pd/validate.hpp"
#include "pd/version.hpp"
#ifdef PD_HAVE_GMP
#include "pd/io/rational.hpp"
#include "pd/kernel/exact.hpp"
#include "pd/verify.hpp"
#endif

namespace {

enum class Stage { Both, Holes, Convex };

enum class Compare { None, Double, Naive };

struct Config {
  Stage stage = Stage::Both;
  bool exact = pd::has_exact_kernel();  // default: exact if built with GMP, else double
  bool check = false;
  bool verify = false;
  Compare compare = Compare::None;
  bool retry_exact = false;
  pd::CutRule cut_rule = pd::CutRule::FarthestPair;
  std::ostream* rational = nullptr;
};

// What one input produced, in kernel-independent form.
struct Outcome {
  pd::Status status = pd::Status::Ok;
  std::string message;
  std::vector<pd::Ring<double>> pieces;              // rounded to double
  std::vector<std::pair<const char*, std::size_t>> counts;
  double ms = 0;
  int verified = -1;  // -1 not run, 0 failed, 1 certified
  std::string verify_failure;
  std::size_t max_bits = 0;  // exact kernel: largest coordinate size (numerator + denominator bits)
};

int usage()
{
  std::cerr << "usage: pd_decompose [--stage both|holes|convex] [--kernel exact|double] [--check]\n"
               "                    [--verify] [--out FILE] [--rational FILE] [--stats FILE] [--quiet]\n"
               "                    (--wkt FILE | --obj OUTER.obj [--hole H.obj ...] | --multi RINGS.obj)\n"
               "                    [--plane xz|xy] [--compare double|naive]\n"
               "                    [--cut-rule farthest|closest|boundary|horizontal]\n";
  return 2;
}

#ifdef PD_HAVE_GMP
std::size_t bits(const mpq_class& q)
{
  return mpz_sizeinbase(q.get_num_mpz_t(), 2) + mpz_sizeinbase(q.get_den_mpz_t(), 2);
}
#endif

template <class K>
Outcome run(const pd::PolygonWithHoles<double>& in, const Config& cfg)
{
  using FT = typename K::FT;
  Outcome o;
  const pd::PolygonWithHoles<FT> input = pd::convert_polygon<K>(in);
  std::vector<pd::Ring<FT>> pieces;
  const std::vector<pd::Point<FT>>* pool = nullptr;

  const auto t0 = std::chrono::steady_clock::now();
  pd::Diagnostics diag;
  pd::PipelineResult<FT> pr;
  pd::HoleElimination<FT> hr;
  pd::ConvexDecomposition<FT> cr;
  if (cfg.stage == Stage::Both) {
    pd::PipelineOptions opt;
    opt.check = cfg.check;
    opt.cut_rule = cfg.cut_rule;
    if (cfg.retry_exact) opt.on_failure = pd::OnFailure::RetryExact;
    pr = pd::decompose<K>(input, opt);
    diag = pr.diag;
    pieces = pr.leaf_rings();
    pool = &pr.pool;
    const auto& s = pr.stats;
    o.counts = {{"cuts", s.holes.cuts},
                {"hole_free_pieces", s.hole_free_pieces},
                {"already_convex", s.already_convex},
                {"reflex_vertices", s.reflex_vertices},
                {"splits", s.splits},
                {"vertex_hits", s.vertex_hits},
                {"double_resolutions", s.double_resolutions},
                {"convex_pieces", s.convex_pieces},
                {"vertex_contacts", s.holes.vertex_contacts},
                {"edge_contacts", s.holes.edge_contacts},
                {"dangles", s.holes.dangles},
                {"fallback_rotated", s.holes.fallback_rotated},
                {"fallback_midpoint", s.holes.fallback_midpoint},
                {"max_depth_holes", s.holes.max_depth},
                {"retried_exact", s.retried_exact ? 1u : 0u},
                {"symbolic", s.symbolic ? 1u : 0u},
                {"vanished", s.vanished}};
  } else if (cfg.stage == Stage::Holes) {
    pd::HoleElimOptions opt;
    opt.check = cfg.check;
    opt.cut_rule = cfg.cut_rule;
    hr = pd::eliminate_holes<K>(input, opt);
    diag = hr.diag;
    pieces = hr.leaf_rings();
    pool = &hr.pool;
    const auto& s = hr.stats;
    o.counts = {{"cuts", s.cuts},
                {"hole_free_pieces", hr.leaves.size()},
                {"max_crossings", s.max_crossings},
                {"vertex_contacts", s.vertex_contacts},
                {"edge_contacts", s.edge_contacts},
                {"dangles", s.dangles},
                {"fallback_rotated", s.fallback_rotated},
                {"fallback_midpoint", s.fallback_midpoint},
                {"degenerate_directions", s.degenerate_directions},
                {"max_depth", s.max_depth},
                {"symbolic", s.symbolic ? 1u : 0u}};
  } else {
    if (!in.holes.empty()) {
      o.status = pd::Status::InvalidInput;
      o.message = "--stage convex needs an input without holes";
      return o;
    }
    pd::ConvexOptions opt;
    opt.check = cfg.check;
    cr = pd::convex_decompose<K>(input.outer, opt);
    diag = cr.diag;
    pieces = cr.leaf_rings();
    pool = &cr.pool;
    const auto& s = cr.stats;
    o.counts = {{"reflex_vertices", s.reflex_vertices},
                {"splits", s.splits},
                {"vertex_hits", s.vertex_hits},
                {"double_resolutions", s.double_resolutions},
                {"convex_pieces", cr.leaves.size()},
                {"max_depth", s.max_depth}};
  }
  o.ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
  o.status = diag.status;
  o.message = diag.message;
  for (const pd::Issue& is : diag.issues) o.message += "; " + pd::describe(is);
  if (o.status != pd::Status::Ok) return o;

  for (const auto& r : pieces) o.pieces.push_back(pd::to_double_ring<K>(r));

#ifdef PD_HAVE_GMP
  if constexpr (pd::is_exact_kernel<K>::value) {
    for (const auto& p : *pool) o.max_bits = std::max({o.max_bits, bits(p.x), bits(p.y)});
    if (cfg.rational) pd::write_rational_rings(*cfg.rational, pieces);
  }
  if (cfg.verify) {
    // The certificate is exact; double results are converted exactly first.
    using EK = pd::ExactKernel;
    std::vector<pd::Ring<mpq_class>> exact_pieces;
    if constexpr (pd::is_exact_kernel<K>::value) exact_pieces = pieces;
    else
      for (const auto& r : o.pieces) exact_pieces.push_back(pd::convert_ring<EK>(r));
    // Hole-free pieces may be weakly simple after a dangling contact.
    const pd::VerifyReport rep = pd::verify_partition<EK>(pd::convert_polygon<EK>(in), exact_pieces,
                                                          cfg.stage != Stage::Holes,
                                                          cfg.stage == Stage::Holes && hr.stats.dangles > 0);
    o.verified = rep.ok ? 1 : 0;
    o.verify_failure = rep.failure;
  }
#else
  (void)pool;
#endif
  return o;
}

#ifdef PD_HAVE_GMP
// Result of running the exact pipeline and a floating-point one on one input.
struct Comparison {
  std::string category;  // agree, deviate_valid, deviate_suspect, reported_failure or exact_failed
  pd::Status exact_status = pd::Status::Ok, other_status = pd::Status::Ok;
  std::string other_message;
  bool same_structure = false;
  double max_rel_diff = -1;  // largest coordinate difference / input scale (same structure only)
  double area_rel_err = -1;  // |sum of piece areas - input area| / input area (other kernel)
  bool pieces_simple = false;  // every piece simple (exact test on its stored doubles)
  bool pieces_convex = false;  // every piece convex up to the kernel's angle tolerance
  double exact_ms = 0, other_ms = 0;
  std::size_t exact_pieces = 0, other_pieces = 0;
};

template <class KO>
Comparison compare(const pd::PolygonWithHoles<double>& in, bool check)
{
  using EK = pd::ExactKernel;
  Comparison c;
  pd::PipelineOptions opt;
  opt.check = check;
  opt.keep_tree = false;
  auto t0 = std::chrono::steady_clock::now();
  const auto e = pd::decompose<EK>(pd::convert_polygon<EK>(in), opt);
  auto t1 = std::chrono::steady_clock::now();
  opt.validate_input = false;  // validated (exactly) by the run above
  opt.self_check = KO::robust;  // the naive baseline runs without any protection
  const auto d = e.ok() ? pd::decompose<KO>(in, opt) : pd::PipelineResult<double>{};
  auto t2 = std::chrono::steady_clock::now();
  c.exact_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
  c.other_ms = std::chrono::duration<double, std::milli>(t2 - t1).count();
  c.exact_status = e.diag.status;
  c.other_status = d.diag.status;
  c.other_message = d.diag.message;
  c.exact_pieces = e.leaves.size();
  c.other_pieces = d.leaves.size();
  if (!e.ok()) {
    c.category = "exact_failed";
    return c;
  }
  if (!d.ok()) {
    c.category = "reported_failure";
    return c;
  }
  double scale = 1, area_in = std::abs(pd::signed_area2(in.outer));
  for (const auto& p : in.outer) scale = std::max({scale, std::abs(p.x), std::abs(p.y)});
  for (const auto& h : in.holes) area_in -= std::abs(pd::signed_area2(h));
  double area_out = 0;
  for (const auto& l : d.leaves) area_out += pd::signed_area2(d.ring_of(l));
  c.area_rel_err = std::abs(area_out - area_in) / area_in;
  c.same_structure = d.pool.size() == e.pool.size() && d.leaves == e.leaves && d.hole_free == e.hole_free;
  if (c.same_structure) {
    c.max_rel_diff = 0;
    for (std::size_t i = 0; i < e.pool.size(); ++i)
      c.max_rel_diff = std::max({c.max_rel_diff, std::abs(d.pool[i].x - e.pool[i].x.get_d()) / scale,
                                 std::abs(d.pool[i].y - e.pool[i].y.get_d()) / scale});
  }
  c.pieces_simple = c.pieces_convex = true;
  const double tol = std::max(KO::angle_tolerance, 1e-12);
  for (const auto& l : d.leaves) {
    const pd::Ring<double> r = d.ring_of(l);
    c.pieces_simple = c.pieces_simple && pd::validate<pd::DoubleKernel>({r, {}}).ok() &&
                      pd::ring_orientation<pd::DoubleKernel>(r) > 0;
    for (std::size_t k = 0; k < r.size() && c.pieces_convex; ++k) {
      const auto& a = r[pd::prev_index(k, r.size())];
      const auto& b = r[k];
      const auto& n3 = r[pd::next_index(k, r.size())];
      if (pd::DoubleKernel::orient(a, b, n3) >= 0) continue;
      const double ux = b.x - a.x, uy = b.y - a.y, vx = n3.x - b.x, vy = n3.y - b.y;
      c.pieces_convex = ux * vx + uy * vy > 0 &&
                        std::abs(ux * vy - uy * vx) <= tol * std::hypot(ux, uy) * std::hypot(vx, vy);
    }
  }
  if (c.same_structure && c.max_rel_diff <= 1e-9) c.category = "agree";
  else if (c.pieces_simple && c.pieces_convex && c.area_rel_err <= 1e-9) c.category = "deviate_valid";
  else c.category = "deviate_suspect";
  return c;
}
#endif

// A JSON number, or null if it is not finite.
[[maybe_unused]] std::string json_num(double v)
{
  if (!std::isfinite(v)) return "null";
  std::ostringstream os;
  os << v;
  return os.str();
}

std::string json_escape(const std::string& s)
{
  std::string out;
  for (char c : s) {
    if (c == '"' || c == '\\') out += '\\';
    out += c;
  }
  return out;
}

void write_stats(std::ostream& os, std::size_t index, const pd::PolygonWithHoles<double>& in,
                 const Outcome& o)
{
  std::size_t n = in.outer.size();
  for (const auto& h : in.holes) n += h.size();
  os << "{\"index\":" << index << ",\"status\":\"" << pd::to_string(o.status) << "\""
     << ",\"vertices\":" << n << ",\"holes\":" << in.holes.size() << ",\"pieces\":" << o.pieces.size();
  for (const auto& c : o.counts) os << ",\"" << c.first << "\":" << c.second;
  if (o.max_bits) os << ",\"max_bits\":" << o.max_bits;
  os << ",\"time_ms\":" << o.ms;
  if (o.verified >= 0) os << ",\"verified\":" << (o.verified ? "true" : "false");
  if (!o.message.empty()) os << ",\"message\":\"" << json_escape(o.message) << "\"";
  if (!o.verify_failure.empty()) os << ",\"verify_failure\":\"" << json_escape(o.verify_failure) << "\"";
  os << "}\n";
}

} // namespace

int main(int argc, char** argv)
{
  Config cfg;
  std::string wkt, obj_outer, multi, out_path, rational_path, stats_path;
  std::vector<std::string> obj_holes;
  pd::ObjPlane plane = pd::ObjPlane::XZ;
  bool quiet = false;
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    const bool v = i + 1 < argc;
    if (a == "--stage" && v) {
      const std::string s = argv[++i];
      if (s == "both") cfg.stage = Stage::Both;
      else if (s == "holes") cfg.stage = Stage::Holes;
      else if (s == "convex") cfg.stage = Stage::Convex;
      else return usage();
    } else if (a == "--kernel" && v) {
      const std::string k = argv[++i];
      if (k == "exact") cfg.exact = true;
      else if (k == "double") cfg.exact = false;
      else return usage();
    } else if (a == "--plane" && v) {
      const std::string p = argv[++i];
      if (p == "xz") plane = pd::ObjPlane::XZ;
      else if (p == "xy") plane = pd::ObjPlane::XY;
      else return usage();
    } else if (a == "--compare" && v) {
      const std::string k = argv[++i];
      if (k == "double") cfg.compare = Compare::Double;
      else if (k == "naive") cfg.compare = Compare::Naive;
      else return usage();
    } else if (a == "--cut-rule" && v) {
      const std::string r = argv[++i];
      if (r == "farthest") cfg.cut_rule = pd::CutRule::FarthestPair;
      else if (r == "closest") cfg.cut_rule = pd::CutRule::ClosestPair;
      else if (r == "boundary") cfg.cut_rule = pd::CutRule::NearestToBoundary;
      else if (r == "horizontal") cfg.cut_rule = pd::CutRule::Horizontal;
      else return usage();
    } else if (a == "--check") cfg.check = true;
    else if (a == "--retry-exact") cfg.retry_exact = true;
    else if (a == "--verify") cfg.verify = true;
    else if (a == "--quiet") quiet = true;
    else if (a == "--out" && v) out_path = argv[++i];
    else if (a == "--rational" && v) rational_path = argv[++i];
    else if (a == "--stats" && v) stats_path = argv[++i];
    else if (a == "--wkt" && v) wkt = argv[++i];
    else if (a == "--obj" && v) obj_outer = argv[++i];
    else if (a == "--hole" && v) obj_holes.push_back(argv[++i]);
    else if (a == "--multi" && v) multi = argv[++i];
    else return usage();
  }
  if (int(!wkt.empty()) + int(!obj_outer.empty()) + int(!multi.empty()) != 1) return usage();
#ifndef PD_HAVE_GMP
  if (cfg.exact || cfg.verify || !rational_path.empty()) {
    std::cerr << "error: built without GMP; use --kernel double and no --verify/--rational\n";
    return 2;
  }
#endif
#ifndef PD_HAVE_GMP
  if (cfg.compare != Compare::None) {
    std::cerr << "error: --compare needs GMP\n";
    return 2;
  }
#endif
  if (!rational_path.empty() && !cfg.exact) {
    std::cerr << "error: --rational needs --kernel exact\n";
    return 2;
  }

  try {
    std::vector<pd::PolygonWithHoles<double>> inputs;
    if (!wkt.empty()) {
      inputs = pd::read_wkt_file(wkt);
    } else if (!obj_outer.empty()) {
      pd::PolygonWithHoles<double> p;
      p.outer = pd::read_obj_ring(obj_outer, plane);
      for (const auto& h : obj_holes) p.holes.push_back(pd::read_obj_ring(h, plane));
      inputs.push_back(std::move(p));
    } else {
      std::vector<pd::Ring<double>> rings = pd::read_obj_rings(multi, plane);
      std::size_t big = 0;  // the largest ring is the outer one
      for (std::size_t i = 1; i < rings.size(); ++i)
        if (std::abs(pd::signed_area2(rings[i])) > std::abs(pd::signed_area2(rings[big]))) big = i;
      pd::PolygonWithHoles<double> p{rings[big], {}};
      for (std::size_t i = 0; i < rings.size(); ++i)
        if (i != big) p.holes.push_back(rings[i]);
      inputs.push_back(std::move(p));
    }

    std::ofstream out_file, rational_file, stats_file;
    auto open = [](std::ofstream& f, const std::string& path) {
      f.open(path);
      if (!f) throw std::runtime_error("cannot write " + path);
    };
    if (!out_path.empty()) open(out_file, out_path);
    if (!rational_path.empty()) {
      open(rational_file, rational_path);
      cfg.rational = &rational_file;
    }
    if (!stats_path.empty()) open(stats_file, stats_path);
    std::ostream& out = out_path.empty() ? std::cout : out_file;

#ifdef PD_HAVE_GMP
    if (cfg.compare != Compare::None) {
      std::map<std::string, std::size_t> tally;
      for (std::size_t i = 0; i < inputs.size(); ++i) {
        const Comparison c = cfg.compare == Compare::Double ? compare<pd::DoubleKernel>(inputs[i], cfg.check)
                                                            : compare<pd::NaiveDoubleKernel>(inputs[i], cfg.check);
        ++tally[c.category];
        if (stats_file.is_open()) {
          stats_file << "{\"index\":" << i << ",\"category\":\"" << c.category << "\""
                     << ",\"exact_status\":\"" << pd::to_string(c.exact_status) << "\""
                     << ",\"other_status\":\"" << pd::to_string(c.other_status) << "\""
                     << ",\"same_structure\":" << (c.same_structure ? "true" : "false")
                     << ",\"max_rel_diff\":" << json_num(c.max_rel_diff)
                     << ",\"area_rel_err\":" << json_num(c.area_rel_err)
                     << ",\"pieces_simple\":" << (c.pieces_simple ? "true" : "false")
                     << ",\"pieces_convex\":" << (c.pieces_convex ? "true" : "false")
                     << ",\"exact_pieces\":" << c.exact_pieces << ",\"other_pieces\":" << c.other_pieces
                     << ",\"exact_ms\":" << c.exact_ms << ",\"other_ms\":" << c.other_ms;
          if (!c.other_message.empty())
            stats_file << ",\"other_message\":\"" << json_escape(c.other_message) << "\"";
          stats_file << "}\n";
        }
      }
      if (!quiet) {
        std::cerr << "inputs " << inputs.size();
        for (const auto& t : tally) std::cerr << ", " << t.first << " " << t.second;
        std::cerr << '\n';
      }
      return tally.count("exact_failed") ? 1 : 0;
    }
#endif
    std::size_t ok = 0, invalid = 0, failed = 0, certified = 0, not_certified = 0, pieces = 0;
    double total_ms = 0;
    for (std::size_t i = 0; i < inputs.size(); ++i) {
      Outcome o;
#ifdef PD_HAVE_GMP
      o = cfg.exact ? run<pd::ExactKernel>(inputs[i], cfg) : run<pd::DoubleKernel>(inputs[i], cfg);
#else
      o = run<pd::DoubleKernel>(inputs[i], cfg);
#endif
      total_ms += o.ms;
      if (o.status == pd::Status::Ok) {
        ++ok;
        pieces += o.pieces.size();
        pd::write_wkt_multipolygon(out, o.pieces);
        out << '\n';
      } else {
        (o.status == pd::Status::InvalidInput ? invalid : failed)++;
        out << "MULTIPOLYGON EMPTY\n";
      }
      if (o.verified == 1) ++certified;
      if (o.verified == 0) ++not_certified;
      if (stats_file.is_open()) write_stats(stats_file, i, inputs[i], o);
      if (!quiet && o.status != pd::Status::Ok)
        std::cerr << "input " << i << ": " << pd::to_string(o.status) << ": " << o.message << '\n';
      if (!quiet && o.verified == 0)
        std::cerr << "input " << i << ": not certified: " << o.verify_failure << '\n';
    }
    if (!quiet) {
      std::cerr << "inputs " << inputs.size() << ", ok " << ok << ", invalid " << invalid
                << ", failed " << failed;
      if (cfg.verify) std::cerr << ", certified " << certified << ", not certified " << not_certified;
      std::cerr << ", pieces " << pieces << ", time " << total_ms << " ms\n";
    }
    return (ok == inputs.size() && not_certified == 0) ? 0 : 1;
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << '\n';
    return 2;
  }
}
