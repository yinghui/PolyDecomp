// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#include "pd/generate.hpp"

#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>

#include "pd/geometry.hpp"
#include "pd/kernel/double.hpp"
#include "pd/validate.hpp"

namespace pd::gen {

namespace {

constexpr double kSize = 10.0;     // the outer square is [0, kSize]^2
constexpr double kMargin = 0.05;   // minimum distance between any two rings
constexpr double kTwoPi = 6.283185307179586;

class Rng {
public:
  explicit Rng(std::uint64_t seed) : g_(seed) {}
  double uniform() { return static_cast<double>(g_() >> 11) * 0x1.0p-53; }  // [0, 1)
  double uniform(double a, double b) { return a + (b - a) * uniform(); }
  int integer(int lo, int hi)  // inclusive
  {
    return lo + static_cast<int>(g_() % static_cast<std::uint64_t>(hi - lo + 1));
  }

private:
  std::mt19937_64 g_;
};

double round5(double v) { return std::round(v * 1e5) / 1e5; }

Ring<double> rounded(const Ring<double>& r)
{
  Ring<double> out;
  for (const auto& p : r) out.push_back({round5(p.x), round5(p.y)});
  return out;
}

// ---- hole shapes centred at (cx, cy) with size r and rotation a

Ring<double> convex_hole(Rng& g, double cx, double cy, double r)
{
  const int n = g.integer(3, 8);
  const double a0 = g.uniform(0, kTwoPi);
  Ring<double> h;
  for (int k = 0; k < n; ++k) {  // points on a circle in angular order: convex
    const double a = a0 + kTwoPi * (k + g.uniform(-0.3, 0.3)) / n;
    h.push_back({cx + r * std::cos(a), cy + r * std::sin(a)});
  }
  return h;
}

Ring<double> star_hole(Rng& g, double cx, double cy, double r)
{
  const int spikes = g.integer(3, 8);
  const double a0 = g.uniform(0, kTwoPi);
  Ring<double> h;
  for (int k = 0; k < 2 * spikes; ++k) {
    const double a = a0 + kTwoPi * k / (2 * spikes);
    const double rr = (k % 2 == 0) ? r * g.uniform(0.8, 1.0) : r * g.uniform(0.3, 0.6);
    h.push_back({cx + rr * std::cos(a), cy + rr * std::sin(a)});
  }
  return h;
}

// A C: a band between radii r and r - t over an arc of 250..310 degrees. Its
// vertex average lies in the empty middle, outside the ring.
Ring<double> c_hole(Rng& g, double cx, double cy, double r)
{
  const double t = r * g.uniform(0.25, 0.4);
  const double span = g.uniform(250, 310) * kTwoPi / 360;
  const double a0 = g.uniform(0, kTwoPi);
  const int s = g.integer(4, 8);
  Ring<double> h;
  for (int k = 0; k <= s; ++k) {
    const double a = a0 + span * k / s;
    h.push_back({cx + r * std::cos(a), cy + r * std::sin(a)});
  }
  for (int k = s; k >= 0; --k) {
    const double a = a0 + span * k / s;
    h.push_back({cx + (r - t) * std::cos(a), cy + (r - t) * std::sin(a)});
  }
  return h;
}

// ---- placement tests (double precision; the final polygon is validated exactly)

double point_segment_distance(const Point<double>& p, const Point<double>& a, const Point<double>& b)
{
  const double dx = b.x - a.x, dy = b.y - a.y;
  const double len2 = dx * dx + dy * dy;
  double t = len2 > 0 ? ((p.x - a.x) * dx + (p.y - a.y) * dy) / len2 : 0;
  t = std::clamp(t, 0.0, 1.0);
  return std::hypot(p.x - (a.x + t * dx), p.y - (a.y + t * dy));
}

double ring_distance(const Ring<double>& r, const Ring<double>& s)
{
  double d = INFINITY;
  for (std::size_t i = 0; i < r.size(); ++i)
    for (std::size_t j = 0; j < s.size(); ++j) {
      const auto& a = r[i];
      const auto& b = r[next_index(i, r.size())];
      const auto& c = s[j];
      const auto& e = s[next_index(j, s.size())];
      if (segments_touch<DoubleKernel>(a, b, c, e)) return 0;
      d = std::min({d, point_segment_distance(a, c, e), point_segment_distance(b, c, e),
                    point_segment_distance(c, a, b), point_segment_distance(e, a, b)});
    }
  return d;
}

bool fits(const Ring<double>& h, const std::vector<Ring<double>>& placed)
{
  for (const auto& p : h)
    if (p.x < kMargin || p.x > kSize - kMargin || p.y < kMargin || p.y > kSize - kMargin) return false;
  for (const auto& other : placed) {
    if (ring_distance(h, other) < kMargin) return false;
    // No nesting either way (boundaries are apart, so one vertex decides).
    if (locate<DoubleKernel>(h[0], other) != Location::Outside) return false;
    if (locate<DoubleKernel>(other[0], h) != Location::Outside) return false;
  }
  return true;
}

Ring<double> outer_square() { return {{0, 0}, {kSize, 0}, {kSize, kSize}, {0, kSize}}; }

Ring<double> shape(Rng& g, int kind, double cx, double cy, double r)
{
  if (kind == 0) return convex_hole(g, cx, cy, r);
  if (kind == 1) return star_hole(g, cx, cy, r);
  return c_hole(g, cx, cy, r);
}

// Places up to `count` holes by rejection sampling at random centres in a disc
// (centre c, radius spread), with sizes in [rmin, rmax] and shapes from `kind`
// (-1: random mix).
std::vector<Ring<double>> scatter(Rng& g, int count, int kind, double rmin, double rmax,
                                  double spread)
{
  std::vector<Ring<double>> holes;
  for (int attempt = 0; attempt < 200 * count && static_cast<int>(holes.size()) < count; ++attempt) {
    const double a = g.uniform(0, kTwoPi), d = spread * std::sqrt(g.uniform());
    const double cx = kSize / 2 + d * std::cos(a), cy = kSize / 2 + d * std::sin(a);
    const int k = kind >= 0 ? kind : g.integer(0, 2);
    Ring<double> h = rounded(shape(g, k, cx, cy, g.uniform(rmin, rmax)));
    if (k == 2 && locate<DoubleKernel>(Point<double>{cx, cy}, h) != Location::Outside) continue;
    if (fits(h, holes)) holes.push_back(std::move(h));
  }
  return holes;
}

PolygonWithHoles<double> attempt(Family f, Rng& g)
{
  PolygonWithHoles<double> p{outer_square(), {}};
  switch (f) {
  case Family::ConvexRow: {
    const int k = g.integer(2, 8);
    const double a = g.uniform(0, kTwoPi / 2);
    const double half = 4.0, step = 2 * half / k;
    const double r = std::min(0.4 * step, g.uniform(0.3, 1.0));
    for (int i = 0; i < k; ++i) {
      const double s = -half + step * (i + 0.5);
      const double cx = kSize / 2 + s * std::cos(a), cy = kSize / 2 + s * std::sin(a);
      Ring<double> h = rounded(convex_hole(g, cx, cy, r));
      if (fits(h, p.holes)) p.holes.push_back(std::move(h));
    }
    break;
  }
  case Family::StarHoles: p.holes = scatter(g, g.integer(2, 6), 1, 0.6, 1.6, 3.2); break;
  case Family::CShapes: p.holes = scatter(g, g.integer(2, 5), 2, 0.8, 1.8, 3.0); break;
  case Family::Cluster: p.holes = scatter(g, g.integer(3, 8), -1, 0.3, 0.8, 1.8); break;
  case Family::Grid: {
    // m x m aligned squares of half-size s on a regular grid: squares in the
    // same row or column share coordinates exactly, so many vertices and
    // centroids are collinear.
    const int m = g.integer(2, 5);
    const double step = kSize / (m + 1), s = step * (g.integer(1, 3) / 8.0);
    for (int i = 1; i <= m; ++i)
      for (int j = 1; j <= m; ++j) {
        if (g.uniform() < 0.25) continue;  // leave some cells empty
        const double cx = step * i, cy = step * j;
        p.holes.push_back({{cx - s, cy - s}, {cx + s, cy - s}, {cx + s, cy + s}, {cx - s, cy + s}});
      }
    if (p.holes.empty()) p.holes.push_back({{4, 4}, {6, 4}, {6, 6}, {4, 6}});
    break;
  }
  default: break;
  }
  return p;
}

// Nudges v by k ulps (k may be negative).
double nudge(double v, int k)
{
  for (; k > 0; --k) v = std::nextafter(v, INFINITY);
  for (; k < 0; ++k) v = std::nextafter(v, -INFINITY);
  return v;
}

// Replaces every edge of r by `parts` edges whose inner vertices lie a few ulps
// off the edge line (on either side), so their angles are within ~1e-15 of 180.
Ring<double> subdivide_nudged(Rng& g, const Ring<double>& r, int parts)
{
  Ring<double> out;
  for (std::size_t i = 0; i < r.size(); ++i) {
    const auto& a = r[i];
    const auto& b = r[next_index(i, r.size())];
    out.push_back(a);
    for (int k = 1; k < parts; ++k) {
      const double t = double(k) / parts;
      Point<double> m{a.x + t * (b.x - a.x), a.y + t * (b.y - a.y)};
      const int u = g.integer(-3, 3);
      if (a.x == b.x) m.x = nudge(m.x, u);  // vertical edge: move across it
      else m.y = nudge(m.y, u);
      out.push_back(m);
    }
  }
  return out;
}

PolygonWithHoles<double> adversarial_attempt(Family f, Rng& g)
{
  PolygonWithHoles<double> p{outer_square(), {}};
  switch (f) {
  case Family::NearCollinear: {
    // Shifted to [1, 11]^2: a one-ulp nudge of a coordinate 0 would be a
    // subnormal, outside the double kernel's coordinate range.
    p.outer = subdivide_nudged(g, {{1, 1}, {1 + kSize, 1}, {1 + kSize, 1 + kSize}, {1, 1 + kSize}},
                               g.integer(2, 6));
    const int m = g.integer(2, 4);
    const double step = kSize / (m + 1), s = step / 4;
    for (int i = 1; i <= m; ++i)
      for (int j = 1; j <= m; ++j) {
        if (g.uniform() < 0.3) continue;
        const double cx = 1 + step * i, cy = 1 + step * j;
        p.holes.push_back(subdivide_nudged(
            g, {{cx - s, cy - s}, {cx + s, cy - s}, {cx + s, cy + s}, {cx - s, cy + s}}, g.integer(2, 5)));
      }
    break;
  }
  case Family::NearTouch: {
    // Columns of rectangles separated by gaps of 2^-30.
    const double gap = std::ldexp(1.0, -30);
    const int cols = g.integer(2, 6);
    const double w = (kSize - (cols + 1) * gap) / cols;
    for (int c = 0; c < cols; ++c) {
      const double x0 = gap + c * (w + gap), x1 = x0 + w;
      const double y0 = g.uniform() < 0.5 ? gap : 1.0 + g.integer(0, 3);
      const double y1 = g.uniform() < 0.5 ? kSize - gap : kSize - 1.0 - g.integer(0, 3);
      p.holes.push_back({{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}});
    }
    break;
  }
  case Family::Snapped: {
    std::vector<Ring<double>> holes = scatter(g, g.integer(3, 8), -1, 0.4, 1.0, 2.5);
    for (Ring<double>& h : holes) {
      Ring<double> snapped;
      for (const auto& q : h) {
        const Point<double> s{std::round(q.x * 8) / 8, std::round(q.y * 8) / 8};
        if (snapped.empty() || snapped.back() != s) snapped.push_back(s);
      }
      while (snapped.size() > 1 && snapped.front() == snapped.back()) snapped.pop_back();
      if (snapped.size() >= 3) p.holes.push_back(std::move(snapped));
    }
    break;
  }
  default: break;
  }
  return p;
}

} // namespace

const std::vector<Family>& adversarial_families()
{
  static const std::vector<Family> v{Family::NearCollinear, Family::NearTouch, Family::Snapped};
  return v;
}

const std::vector<Family>& all_families()
{
  static const std::vector<Family> v{Family::ConvexRow, Family::StarHoles, Family::CShapes,
                                     Family::Cluster, Family::Grid};
  return v;
}

const char* name(Family f)
{
  switch (f) {
  case Family::ConvexRow: return "row";
  case Family::StarHoles: return "star";
  case Family::CShapes: return "c";
  case Family::Cluster: return "cluster";
  case Family::Grid: return "grid";
  case Family::NearCollinear: return "nearcollinear";
  case Family::NearTouch: return "neartouch";
  case Family::Snapped: return "snapped";
  }
  return "?";
}

Family family_from_name(const std::string& s)
{
  for (Family f : all_families())
    if (s == name(f)) return f;
  for (Family f : adversarial_families())
    if (s == name(f)) return f;
  throw std::invalid_argument("unknown family '" + s + "'");
}

PolygonWithHoles<double> generate(Family f, std::uint64_t seed)
{
  // Mix the family into the seed so families with equal seeds differ.
  Rng g(seed * 0x9E3779B97F4A7C15ULL + static_cast<std::uint64_t>(f) + 1);
  const bool adversarial = f == Family::NearCollinear || f == Family::NearTouch || f == Family::Snapped;
  for (int tries = 0; tries < 1000; ++tries) {
    PolygonWithHoles<double> p = adversarial ? adversarial_attempt(f, g) : attempt(f, g);
    if (!p.holes.empty() && validate<DoubleKernel>(p).ok()) return p;
  }
  throw std::runtime_error("generator failed to produce a valid polygon");
}

Ring<double> comb(int teeth, double h)
{
  Ring<double> r{{0, 0}, {2.0 * teeth - 1, 0}};
  for (int j = teeth - 1; j >= 0; --j) {
    r.push_back({2.0 * j + 1, h});
    r.push_back({2.0 * j, h});
    if (j > 0) {
      r.push_back({2.0 * j, 1});
      r.push_back({2.0 * j - 1, 1});
    }
  }
  return r;
}

Ring<double> star_polygon(int n, std::uint64_t seed)
{
  Rng g(seed * 0x9E3779B97F4A7C15ULL + 0x5a17);
  for (int tries = 0; tries < 100; ++tries) {
    std::vector<double> a(static_cast<std::size_t>(n));
    for (double& v : a) v = g.uniform(0, kTwoPi);
    std::sort(a.begin(), a.end());
    Ring<double> r;
    for (double v : a) {
      const double rad = g.uniform(3, 10);
      const Point<double> p{round5(10 + rad * std::cos(v)), round5(10 + rad * std::sin(v))};
      if (r.empty() || r.back() != p) r.push_back(p);
    }
    if (r.size() >= 3 && validate<DoubleKernel>({r, {}}).ok()) return r;
  }
  throw std::runtime_error("star_polygon: no valid polygon");
}

PolygonWithHoles<double> spiral_holes(int h, std::uint64_t seed)
{
  Rng g(seed * 0x9E3779B97F4A7C15ULL + 0x5b1a);
  for (int tries = 0; tries < 100; ++tries) {
    PolygonWithHoles<double> p{{{-100, -100}, {100, -100}, {100, 100}, {-100, 100}}, {}};
    for (int k = 0; k < h; ++k) {
      const double a = 0.7 * k, r = 3 + 1.2 * k * (64.0 / h);
      const double cx = r * std::cos(a) + g.uniform(-0.05, 0.05), cy = r * std::sin(a) + g.uniform(-0.05, 0.05);
      p.holes.push_back(rounded({{cx - 0.4, cy - 0.3 + g.uniform(-0.05, 0.05)},
                                 {cx + 0.4 + g.uniform(-0.05, 0.05), cy - 0.4},
                                 {cx + 0.3, cy + 0.4 + g.uniform(-0.05, 0.05)},
                                 {cx - 0.4 + g.uniform(-0.05, 0.05), cy + 0.3}}));
    }
    if (validate<DoubleKernel>(p).ok()) return p;
  }
  throw std::runtime_error("spiral_holes: no valid polygon");
}

PolygonWithHoles<double> random_holes(int h, std::uint64_t seed)
{
  if (h < 1) throw std::invalid_argument("random_holes: h must be at least 1");
  Rng g(seed * 0x9E3779B97F4A7C15ULL + 0x7a3d);
  const int k = static_cast<int>(std::ceil(std::sqrt(static_cast<double>(h))));
  const double side = k;
  PolygonWithHoles<double> p{{{0, 0}, {side, 0}, {side, side}, {0, side}}, {}};
  p.holes.reserve(static_cast<std::size_t>(h));
  for (int c = 0; c < h; ++c) {
    // Centre within 0.05 of the cell centre, radii in [0.15, 0.4]: the hole stays
    // inside [0.05, 0.95] of its cell. Angles increase strictly (jitter below half
    // a step), so the ring is star-shaped about its centre and simple.
    const double cx = c % k + 0.5 + g.uniform(-0.05, 0.05), cy = c / k + 0.5 + g.uniform(-0.05, 0.05);
    const int m = g.integer(3, 10);
    const double a0 = g.uniform(0, kTwoPi), step = kTwoPi / m;
    Ring<double> r;
    for (int i = 0; i < m; ++i) {
      const double a = a0 + (i + g.uniform(-0.3, 0.3)) * step, rad = g.uniform(0.15, 0.4);
      r.push_back({cx + rad * std::cos(a), cy + rad * std::sin(a)});
    }
    p.holes.push_back(rounded(r));
  }
  if (!validate<DoubleKernel>(p).ok()) throw std::runtime_error("random_holes: invalid polygon");
  return p;
}

std::vector<std::pair<std::string, PolygonWithHoles<double>>> degenerate_set()
{
  auto sq = [](double x0, double y0, double x1, double y1) {
    return Ring<double>{{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}};
  };
  const Ring<double> c_hole{{1, 1}, {9, 1}, {9, 3}, {3, 3}, {3, 7}, {9, 7}, {9, 9}, {1, 9}};
  const Ring<double> diamond_l{{-9, 0}, {-8, -1}, {-7, 0}, {-8, 1}};
  const Ring<double> diamond_r{{7, 0}, {8, -1}, {9, 0}, {8, 1}};
  return {
      {"concentric-squares", {sq(-2, -2, 2, 2), {sq(-0.5, -0.5, 0.5, 0.5)}}},
      {"centroid-is-own-vertex", {sq(-10, -10, 10, 10), {{{0, 3}, {3, -3}, {0, -1}, {-3, -3}}}}},
      {"centroid-on-other-hole-vertex", {sq(0, 0, 10, 10), {c_hole, {{5.5, 5}, {8, 4}, {8, 6}}}}},
      {"equal-centroids", {sq(0, 0, 10, 10), {c_hole, sq(5, 4.5, 6, 5.5)}}},
      {"anchor-equals-centroid", {{{0, 0}, {12, 0}, {0, 12}}, {sq(3, 3, 5, 5)}}},
      {"diagonal-squares", {sq(0, 0, 10, 10), {sq(1, 1, 2, 2), sq(4.5, 4.5, 5.5, 5.5), sq(8, 8, 9, 9)}}},
      // The cut line y = 0 through the diamonds' centroids (-8, 0) and (8, 0) passes
      // through vertices of every kind: diamond corners (edges on both sides), the
      // outer ring's notch (0, 0), a hole tip (3, 0) and a hole's notch (-4.5, 0).
      {"vertex-contacts",
       {{{-12, -6}, {12, -6}, {12, 6}, {1, 6}, {0, 0}, {-1, 6}, {-12, 6}},
        {diamond_l, diamond_r, {{2, -3}, {6, -3}, {6, 1}, {5, 1}, {4, -1}, {3, 0}, {2, -1}},
         {{-6, -2}, {-3, -2}, {-3, 2}, {-4, 2}, {-4.5, 0}, {-5, 2}, {-6, 2}}}}},
      // ... and along edges: one with its neighbours on both sides, one on one side.
      {"edge-contacts",
       {sq(-12, -6, 12, 6),
        {diamond_l, diamond_r, {{-5, -2}, {-3, 0}, {-1, 0}, {0, 2}, {-4, 2}}, {{2, 0}, {4, 0}, {5, 2}, {1, 2}}}}},
      // ... and along edges of the outer ring: a step (-12,0)-(-10,0), its neighbours on both
      // sides (a wing), and the bottom (1,0)-(3,0) of a slot, its neighbours on one side (a wedge).
      {"ring-edge-contacts",
       {{{-10, -6}, {12, -6}, {12, 6}, {3, 6}, {3, 0}, {1, 0}, {1, 6}, {-12, 6}, {-12, 0}, {-10, 0}},
        {diamond_l, diamond_r}}},
      // ... and through single vertices of two triangles that lie otherwise on one
      // side (dangling contacts: split symbolically by the exact kernel).
      {"dangling-contacts",
       {sq(-12, -6, 12, 6), {diamond_l, diamond_r, {{0, 0}, {1, 2}, {-1, 2}}, {{3, 0}, {4, -2}, {2, -2}}}}},
  };
}

} // namespace pd::gen
