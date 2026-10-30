// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#pragma once
/// \file
/// Exact cutting of a polygon with holes along arbitrary lines, as often as wanted (for example,
/// slicing a shape in a game).
///
/// A Slicer keeps one state for all its pieces: the points, the supporting line of every edge (input
/// edges and cut lines) and the holes. A cut along a line splits every piece it crosses with the rules
/// of hole elimination: contacts (the line through vertices or along edges) are classified as if the
/// line were shifted infinitesimally into each side, and a hole that the line touches in a single
/// vertex is split symbolically. Every new point is the intersection of two supporting lines, so the
/// pieces always tile the input, holes that the line misses stay holes of their piece, and the
/// coordinates stay bounded by the sizes of the input and of the cut lines, however many cuts follow.
///
/// All arithmetic is exact. The state is rational and switches to the symbolic kernel (rationals
/// with an infinitesimal ε) at the first dangling contact; pieces are reported at ε = 0, where a piece
/// that runs round such a hole touches itself at one point (it is weakly simple).

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "pd/geometry.hpp"
#include "pd/hole_elim.hpp"
#include "pd/pipeline.hpp"
#include "pd/status.hpp"

#ifdef PD_HAVE_GMP
#include "pd/kernel/exact.hpp"
#include "pd/kernel/symbolic.hpp"

namespace pd {

/// Outcome of one cut.
struct SliceResult {
  Diagnostics diag;                 ///< Ok, InvalidInput (the line or the piece index) or InternalError
  std::vector<std::size_t> pieces;  ///< indices, in the new piece list, of the pieces the cut created
  std::size_t cut_pieces = 0;       ///< pieces that the line split (a piece it misses or only touches is kept)
  std::size_t vertex_contacts = 0;  ///< vertices of the pieces on the line
  std::size_t edge_contacts = 0;    ///< runs of edges on the line
  std::size_t dangles = 0;          ///< holes touched in a single vertex (split symbolically)
  bool ok() const { return diag.status == Status::Ok; }
};

namespace detail {

/// The persistent state of a Slicer in one kernel (ExactKernel or SymbolicKernel).
template <class K>
struct SliceState {
  using FT = typename K::FT;
  using Elim = HoleEliminator<K>;
  using Work = typename Elim::Work;
  using Cut = typename Elim::Cut;

  HoleElimOptions opt;
  HoleElimination<FT> out;
  Elim elim{out, opt};
  std::vector<Work> pieces;
  /// A box that contains a piece in its interior, in doubles rounded outwards.
  struct Box {
    double x0, y0, x1, y1;
  };
  std::vector<Box> boxes;  ///< one per piece

  /// The line {p : a p.x + b p.y + c = 0} with its coefficients rounded to doubles, for filtering.
  struct LineD {
    double a, b, c;
  };

  Box box_of(const Work& w) const
  {
    double x0 = HUGE_VAL, y0 = HUGE_VAL, x1 = -HUGE_VAL, y1 = -HUGE_VAL;
    for (const Corner1& c : w.outer) {
      const Point<FT>& p = out.pool[c.v];
      const double x = K::to_double(p.x), y = K::to_double(p.y);
      x0 = std::min(x0, x), x1 = std::max(x1, x), y0 = std::min(y0, y), y1 = std::max(y1, y);
    }
    // The conversion truncates towards zero, so each coordinate (and, in the symbolic kernel, each
    // value near its limit) lies strictly between the neighbours of its double.
    return {std::nextafter(x0, -HUGE_VAL), std::nextafter(y0, -HUGE_VAL), std::nextafter(x1, HUGE_VAL),
            std::nextafter(y1, HUGE_VAL)};
  }
  /// The line through q in direction d: a x + b y + c is positive on its left, computed exactly and
  /// then rounded (truncated) to doubles, each with relative error below 2^-52.
  static LineD line_of(const Point<FT>& q, const Point<FT>& d)
  {
    const FT a = -d.y, b = d.x, c = d.y * q.x - d.x * q.y;
    const LineD l{K::to_double(a), K::to_double(b), K::to_double(c)};
    // A nonzero coefficient that rounds to zero or to a subnormal double has no relative error
    // bound, so the filter could wrongly skip a piece that the line crosses. Such a line is never
    // used to skip a piece: NaN coefficients make misses() return false.
    const auto lost = [](const FT& exact, double x) {
      return !(exact == FT(0)) && !(std::abs(x) >= std::numeric_limits<double>::min());
    };
    if (lost(a, l.a) || lost(b, l.b) || lost(c, l.c)) {
      const double nan = std::numeric_limits<double>::quiet_NaN();
      return {nan, nan, nan};
    }
    return l;
  }
  /// True if the whole box lies strictly on one side of the line, certified by a floating point
  /// filter: with exact coefficients the value at a corner differs from the computed one by less
  /// than err (the rounding of the coefficients and of the three operations, plus a margin for
  /// underflow). An uncertain or non-finite case returns false, and the piece is then cut normally.
  static bool misses(const Box& r, const LineD& l)
  {
    const double xs[4] = {r.x0, r.x1, r.x1, r.x0}, ys[4] = {r.y0, r.y0, r.y1, r.y1};
    int side = 0;
    for (int k = 0; k < 4; ++k) {
      const double ax = l.a * xs[k], by = l.b * ys[k];
      const double f = ax + by + l.c;
      const double err = 0x1p-49 * (std::abs(ax) + std::abs(by) + std::abs(l.c)) + 0x1p-1000;
      if (!(std::abs(f) > err)) return false;  // also false for NaN and infinities
      const int t = f > 0 ? 1 : -1;
      if (side != 0 && t != side) return false;
      side = t;
    }
    return true;
  }
  void rebuild_boxes()
  {
    boxes.clear();
    for (const Work& w : pieces) boxes.push_back(box_of(w));
  }

  SliceState() { opt.keep_tree = false; }
  SliceState(const SliceState&) = delete;  // elim refers to out and opt
  SliceState& operator=(const SliceState&) = delete;

  bool init(const PolygonWithHoles<FT>& shape)
  {
    Work root;
    if (!elim.init(shape, root)) return false;
    pieces.assign(1, std::move(root));
    rebuild_boxes();
    return true;
  }

  struct Snapshot {
    std::size_t pool, lines;
    std::vector<std::uint32_t> hole_off, hole_len;
    std::vector<std::vector<std::uint32_t>> hole_line;
    HoleElimStats stats;
  };
  Snapshot snapshot() const
  {
    return {out.pool.size(), elim.lines_.size(), elim.hole_off_, elim.hole_len_, elim.hole_line_, out.stats};
  }
  void restore(const Snapshot& s)
  {
    out.pool.resize(s.pool);
    elim.lines_.resize(s.lines);
    elim.hole_off_ = s.hole_off;
    elim.hole_len_ = s.hole_len;
    elim.hole_line_ = s.hole_line;
    out.stats = s.stats;
    out.diag = {};
    out.needs_symbolic = false;
  }

  /// Cuts piece `which` (every piece if which == npos) along the line q + t d. Returns 1 on success,
  /// 0 on failure (res.diag set, state unchanged), -1 if a dangling contact needs the symbolic kernel
  /// (state unchanged).
  int cut(const Point<FT>& q, const Point<FT>& d, std::size_t which, SliceResult& res)
  {
    constexpr std::size_t npos = static_cast<std::size_t>(-1);
    const Snapshot snap = snapshot();
    const LineD line = line_of(q, d);
    // First cut every piece that the line may cross, keeping the old pieces intact, so that a
    // failure leaves the state unchanged; then replace each cut piece by its parts in place.
    std::vector<std::pair<std::size_t, std::vector<Work>>> parts;
    for (std::size_t i = 0; i < pieces.size(); ++i) {
      if ((which != npos && i != which) || misses(boxes[i], line)) continue;
      Work w = pieces[i];
      std::vector<Work> stack;
      Cut c;
      c.q = q;
      c.d = d;
      if (!elim.cut_along(w, c, stack, false)) {
        const bool symbolic = out.needs_symbolic;
        const Diagnostics diag = out.diag;
        restore(snap);
        if (symbolic) return -1;
        res.diag = diag;
        return 0;
      }
      if (stack.size() < 2) continue;  // the line misses the piece or only touches its boundary
      std::reverse(stack.begin(), stack.end());  // cut_along pushes the pieces in reverse order
      parts.push_back({i, std::move(stack)});
    }
    std::vector<std::size_t> created;
    if (!parts.empty()) {
      std::vector<Work> next;
      std::vector<Box> next_boxes;
      next.reserve(pieces.size() + parts.size());
      next_boxes.reserve(pieces.size() + parts.size());
      std::size_t j = 0;
      for (std::size_t i = 0; i < pieces.size(); ++i) {
        if (j < parts.size() && parts[j].first == i) {
          for (Work& p : parts[j].second) {
            p.node = -1;
            created.push_back(next.size());
            next_boxes.push_back(box_of(p));
            next.push_back(std::move(p));
          }
          ++j;
        } else {
          next.push_back(std::move(pieces[i]));
          next_boxes.push_back(boxes[i]);
        }
      }
      pieces = std::move(next);
      boxes = std::move(next_boxes);
    }
    res.pieces = std::move(created);
    res.cut_pieces = parts.size();
    res.vertex_contacts = out.stats.vertex_contacts - snap.stats.vertex_contacts;
    res.edge_contacts = out.stats.edge_contacts - snap.stats.edge_contacts;
    res.dangles = out.stats.dangles - snap.stats.dangles;
    return 1;
  }

  /// Copies an exact state into this symbolic one (at the first dangling contact). Every
  /// SliceState is a friend of every HoleEliminator, so it may read the other kernel's tables.
  template <class E>
  void copy_from(const SliceState<E>& e)
  {
    auto pt = [](const auto& p) { return Point<FT>{FT(p.x), FT(p.y)}; };
    out.pool.clear();
    out.pool.reserve(e.out.pool.size());
    for (const auto& p : e.out.pool) out.pool.push_back(pt(p));
    elim.lines_.clear();
    for (const auto& l : e.elim.lines_) elim.lines_.push_back({pt(l.p), pt(l.d)});
    elim.hole_off_ = e.elim.hole_off_;
    elim.hole_len_ = e.elim.hole_len_;
    elim.hole_line_ = e.elim.hole_line_;
    elim.scale_ = e.elim.scale_;
    elim.tol_ = e.elim.tol_;
    out.stats = e.out.stats;
    out.stats.symbolic = is_symbolic_kernel<K>::value;
    pieces.clear();
    for (const auto& w : e.pieces) pieces.push_back({w.outer, w.holes, w.node, w.depth});
    rebuild_boxes();
  }

  /// Piece i in the state's number type: its outer ring and the holes it still contains.
  PolygonWithHoles<FT> piece(std::size_t i) const
  {
    PolygonWithHoles<FT> p;
    for (const Corner1& c : pieces[i].outer) p.outer.push_back(out.pool[c.v]);
    for (int h : pieces[i].holes) {
      Ring<FT> r;
      for (std::uint32_t k = 0; k < elim.hole_len_[h]; ++k) r.push_back(out.pool[elim.hole_off_[h] + k]);
      p.holes.push_back(std::move(r));
    }
    return p;
  }
};

inline std::size_t bit_size(const mpq_class& q)
{
  return mpz_sizeinbase(q.get_num_mpz_t(), 2) + mpz_sizeinbase(q.get_den_mpz_t(), 2);
}

} // namespace detail

/// Cuts a polygon with holes along arbitrary lines, exactly, as often as wanted. See the file
/// comment. The input is checked against the input class (see validate) when the Slicer is built.
class Slicer {
public:
  static constexpr std::size_t npos = static_cast<std::size_t>(-1);

  explicit Slicer(const PolygonWithHoles<double>& shape) : Slicer(convert_polygon<ExactKernel>(shape)) {}

  explicit Slicer(const PolygonWithHoles<mpq_class>& shape)
      : exact_(std::make_unique<detail::SliceState<ExactKernel>>())
  {
    if (!exact_->init(shape)) diag_ = exact_->out.diag;
  }

  bool ok() const { return diag_.status == Status::Ok; }
  const Diagnostics& diag() const { return diag_; }  ///< InvalidInput if the shape is not valid
  std::size_t size() const { return symbolic_ ? symbolic_->pieces.size() : exact_->pieces.size(); }
  bool symbolic() const { return symbolic_ != nullptr; }  ///< a dangling contact has occurred

  /// Cuts every piece that the line through a and b crosses.
  SliceResult cut(const Point<double>& a, const Point<double>& b) { return cut_exact(exact(a), exact(b), npos); }
  /// Cuts one piece along the line through a and b.
  SliceResult cut(std::size_t piece, const Point<double>& a, const Point<double>& b)
  {
    return cut_exact(exact(a), exact(b), piece);
  }
  /// The same with an exact line (for example through a vertex created by an earlier cut): every
  /// piece it crosses, or only `piece`.
  SliceResult cut_exact(const Point<mpq_class>& a, const Point<mpq_class>& b, std::size_t piece = npos)
  {
    SliceResult res;
    if (!ok()) {
      res.diag = diag_;
      return res;
    }
    if (a == b) {
      res.diag = {Status::InvalidInput, "the two points of the cut line coincide", {}};
      return res;
    }
    if (piece != npos && piece >= size()) {
      res.diag = {Status::InvalidInput, "no such piece", {}};
      return res;
    }
    const Point<mpq_class> d{b.x - a.x, b.y - a.y};
    if (!symbolic_) {
      const int r = exact_->cut(a, d, piece, res);
      if (r >= 0) return res;
      symbolic_ = std::make_unique<detail::SliceState<SymbolicKernel>>();  // the first dangling contact
      symbolic_->copy_from(*exact_);
      exact_.reset();
    }
    const int r = symbolic_->cut({EpsQ(a.x), EpsQ(a.y)}, {EpsQ(d.x), EpsQ(d.y)}, piece, res);
    if (r < 0) res.diag = {Status::InternalError, "a dangling contact in the symbolic kernel", {}};
    return res;
  }

  /// Piece i with exact coordinates (at ε = 0 after a dangling contact).
  PolygonWithHoles<mpq_class> piece_exact(std::size_t i) const
  {
    if (!symbolic_) return exact_->piece(i);
    const PolygonWithHoles<EpsQ> s = symbolic_->piece(i);
    PolygonWithHoles<mpq_class> p;
    for (const auto& v : s.outer) {
      const Point<mpq_class> l = limit_point(v);
      if (p.outer.empty() || !(p.outer.back() == l)) p.outer.push_back(l);
    }
    while (p.outer.size() > 1 && p.outer.back() == p.outer.front()) p.outer.pop_back();
    for (const auto& h : s.holes) {
      Ring<mpq_class> r;
      for (const auto& v : h) r.push_back(limit_point(v));
      p.holes.push_back(std::move(r));
    }
    return p;
  }

  /// Piece i rounded to doubles (for display; the state itself stays exact).
  PolygonWithHoles<double> piece(std::size_t i) const
  {
    const PolygonWithHoles<mpq_class> e = piece_exact(i);
    auto ring = [](const Ring<mpq_class>& r) {
      Ring<double> out;
      for (const auto& v : r) out.push_back({v.x.get_d(), v.y.get_d()});
      return out;
    };
    PolygonWithHoles<double> p{ring(e.outer), {}};
    for (const auto& h : e.holes) p.holes.push_back(ring(h));
    return p;
  }

  std::vector<PolygonWithHoles<mpq_class>> pieces_exact() const
  {
    std::vector<PolygonWithHoles<mpq_class>> out;
    for (std::size_t i = 0; i < size(); ++i) out.push_back(piece_exact(i));
    return out;
  }

  /// Convex pieces of piece i, from the pipeline (decompose) with exact arithmetic. Pieces are valid
  /// by construction, so the options default to no input validation.
  PipelineResult<mpq_class> convex(std::size_t i, PipelineOptions opt = no_validation()) const
  {
    if (!symbolic_) return decompose<ExactKernel>(exact_->piece(i), opt);
    return detail::pipeline_limit(detail::decompose_once<SymbolicKernel>(symbolic_->piece(i), opt));
  }

  /// Largest bit size (numerator plus denominator) of a coordinate of a piece vertex.
  std::size_t max_bits() const
  {
    std::size_t b = 0;
    for (std::size_t i = 0; i < size(); ++i) {
      const PolygonWithHoles<mpq_class> p = piece_exact(i);
      for (const auto& v : p.outer) b = std::max({b, detail::bit_size(v.x), detail::bit_size(v.y)});
      for (const auto& h : p.holes)
        for (const auto& v : h) b = std::max({b, detail::bit_size(v.x), detail::bit_size(v.y)});
    }
    return b;
  }

private:
  static Point<mpq_class> exact(const Point<double>& p) { return {mpq_class(p.x), mpq_class(p.y)}; }
  static PipelineOptions no_validation()
  {
    PipelineOptions o;
    o.validate_input = false;
    o.keep_tree = false;
    return o;
  }

  std::unique_ptr<detail::SliceState<ExactKernel>> exact_;
  std::unique_ptr<detail::SliceState<SymbolicKernel>> symbolic_;
  Diagnostics diag_;
};

} // namespace pd
#endif
