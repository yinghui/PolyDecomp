// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#include "pd/io/wkt.hpp"

#include <cctype>
#include <charconv>
#include <fstream>
#include <ostream>
#include <sstream>

namespace pd {

namespace {

class WktParser {
public:
  explicit WktParser(const std::string& text) : s_(text) {}

  std::vector<PolygonWithHoles<double>> parse_all()
  {
    std::vector<PolygonWithHoles<double>> out;
    skip_ws();
    while (pos_ < s_.size()) {
      const std::string kw = keyword();
      if (kw == "POLYGON") {
        out.push_back(polygon_body());
      } else if (kw == "MULTIPOLYGON") {
        if (try_empty()) continue;
        expect('(');
        do out.push_back(polygon_body());
        while (accept(','));
        expect(')');
      } else {
        fail("expected POLYGON or MULTIPOLYGON, found '" + kw + "'");
      }
      skip_ws();
    }
    return out;
  }

private:
  const std::string& s_;
  std::size_t pos_ = 0;

  [[noreturn]] void fail(const std::string& msg) const
  {
    throw ParseError("WKT: " + msg + " at offset " + std::to_string(pos_));
  }
  void skip_ws()
  {
    while (pos_ < s_.size() && std::isspace(static_cast<unsigned char>(s_[pos_]))) ++pos_;
  }
  bool accept(char c)
  {
    skip_ws();
    if (pos_ < s_.size() && s_[pos_] == c) {
      ++pos_;
      return true;
    }
    return false;
  }
  void expect(char c)
  {
    if (!accept(c)) fail(std::string("expected '") + c + "'");
  }
  std::string keyword()
  {
    skip_ws();
    std::string kw;
    while (pos_ < s_.size() && std::isalpha(static_cast<unsigned char>(s_[pos_])))
      kw += static_cast<char>(std::toupper(static_cast<unsigned char>(s_[pos_++])));
    return kw;
  }
  // Consumes "EMPTY" if present.
  bool try_empty()
  {
    const std::size_t save = pos_;
    if (keyword() == "EMPTY") return true;
    pos_ = save;
    return false;
  }
  bool at_number()
  {
    skip_ws();
    if (pos_ >= s_.size()) return false;
    const char c = s_[pos_];
    return std::isdigit(static_cast<unsigned char>(c)) || c == '-' || c == '+' || c == '.';
  }
  double number()
  {
    skip_ws();
    const char* first = s_.data() + pos_;
    const char* last = s_.data() + s_.size();
    if (first < last && *first == '+') ++first;  // from_chars rejects a leading '+'
    double v = 0;
    const auto res = std::from_chars(first, last, v);
    if (res.ec != std::errc()) fail("invalid number");
    pos_ = static_cast<std::size_t>(res.ptr - s_.data());
    return v;
  }
  Ring<double> ring()
  {
    Ring<double> r;
    expect('(');
    do {
      const double x = number();
      const double y = number();
      if (at_number()) number();  // ignore Z
      r.push_back({x, y});
    } while (accept(','));
    expect(')');
    if (r.size() > 1 && r.front() == r.back()) r.pop_back();
    return r;
  }
  PolygonWithHoles<double> polygon_body()
  {
    PolygonWithHoles<double> p;
    if (try_empty()) return p;
    expect('(');
    p.outer = ring();
    while (accept(',')) p.holes.push_back(ring());
    expect(')');
    return p;
  }
};

void write_ring(std::ostream& os, const Ring<double>& r)
{
  os << '(';
  for (std::size_t i = 0; i <= r.size(); ++i) {  // repeat the first vertex to close
    const Point<double>& p = r[i % r.size()];
    if (i) os << ", ";
    os << format_double(p.x) << ' ' << format_double(p.y);
  }
  os << ')';
}

void write_polygon_body(std::ostream& os, const PolygonWithHoles<double>& p)
{
  os << '(';
  write_ring(os, p.outer);
  for (const Ring<double>& h : p.holes) {
    os << ", ";
    write_ring(os, h);
  }
  os << ')';
}

} // namespace

std::string format_double(double v)
{
  char buf[64];
  const auto res = std::to_chars(buf, buf + sizeof buf, v);
  return std::string(buf, res.ptr);
}

std::vector<PolygonWithHoles<double>> parse_wkt(const std::string& text)
{
  return WktParser(text).parse_all();
}

std::vector<PolygonWithHoles<double>> read_wkt_file(const std::string& path)
{
  std::ifstream in(path);
  if (!in) throw ParseError("cannot open " + path);
  std::ostringstream ss;
  ss << in.rdbuf();
  return parse_wkt(ss.str());
}

void write_wkt(std::ostream& os, const PolygonWithHoles<double>& p)
{
  if (p.outer.empty()) {
    os << "POLYGON EMPTY";
    return;
  }
  os << "POLYGON ";
  write_polygon_body(os, p);
}

void write_wkt_multipolygon(std::ostream& os, const std::vector<Ring<double>>& rings)
{
  if (rings.empty()) {
    os << "MULTIPOLYGON EMPTY";
    return;
  }
  os << "MULTIPOLYGON (";
  for (std::size_t i = 0; i < rings.size(); ++i) {
    if (i) os << ", ";
    os << '(';
    write_ring(os, rings[i]);
    os << ')';
  }
  os << ')';
}

} // namespace pd
