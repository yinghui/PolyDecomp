// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Yinghui Zhang
#include "pd/io/obj.hpp"

#include <algorithm>
#include <fstream>
#include <sstream>

namespace pd {

namespace {

[[noreturn]] void fail(const std::string& path, int line, const std::string& msg)
{
  throw ParseError(path + (line > 0 ? ":" + std::to_string(line) : std::string()) + ": " + msg);
}

// Converts an OBJ vertex reference ("7", "7/1/2", "-1") to a 0-based index.
std::size_t vertex_index(const std::string& tok, std::size_t nverts, const std::string& path,
                         int line)
{
  const std::string head = tok.substr(0, tok.find('/'));
  long v = 0;
  try {
    v = std::stol(head);
  } catch (const std::exception&) {
    fail(path, line, "bad vertex reference '" + tok + "'");
  }
  const long idx = v > 0 ? v - 1 : static_cast<long>(nverts) + v;  // negative = relative
  if (v == 0 || idx < 0 || idx >= static_cast<long>(nverts))
    fail(path, line, "vertex reference out of range '" + tok + "'");
  return static_cast<std::size_t>(idx);
}

} // namespace

std::vector<Ring<double>> read_obj_rings(const std::string& path, ObjPlane plane)
{
  std::ifstream in(path);
  if (!in) fail(path, 0, "cannot open file");

  std::vector<Point<double>> verts;
  std::vector<std::pair<std::size_t, std::size_t>> edges;
  std::string text;
  int line = 0;
  while (std::getline(in, text)) {
    ++line;
    std::istringstream ls(text);
    std::string tag;
    if (!(ls >> tag)) continue;
    if (tag == "v") {
      double x, y, z;
      if (!(ls >> x >> y >> z)) fail(path, line, "expected three vertex coordinates");
      verts.push_back(plane == ObjPlane::XZ ? Point<double>{x, z} : Point<double>{x, y});
    } else if (tag == "l" || tag == "f") {
      std::vector<std::size_t> idx;
      std::string tok;
      while (ls >> tok) idx.push_back(vertex_index(tok, verts.size(), path, line));
      if (idx.size() < 2) fail(path, line, "element with fewer than two vertices");
      for (std::size_t i = 0; i + 1 < idx.size(); ++i) edges.emplace_back(idx[i], idx[i + 1]);
      if (tag == "f") edges.emplace_back(idx.back(), idx.front());  // faces are closed
    }
  }
  if (verts.empty()) fail(path, 0, "no vertices");
  if (edges.empty()) return {Ring<double>(verts.begin(), verts.end())};

  // Adjacency: every vertex must lie on exactly two edges.
  const std::size_t n = verts.size();
  std::vector<std::vector<std::size_t>> adj(n);
  for (const auto& e : edges) {
    if (e.first == e.second) fail(path, 0, "degenerate edge");
    adj[e.first].push_back(e.second);
    adj[e.second].push_back(e.first);
  }
  for (std::size_t v = 0; v < n; ++v)
    if (adj[v].size() != 2)
      fail(path, 0,
           "vertex " + std::to_string(v + 1) + " has degree " + std::to_string(adj[v].size()) +
               " (the edges do not form disjoint cycles)");

  std::vector<Ring<double>> rings;
  std::vector<bool> seen(n, false);
  for (std::size_t s = 0; s < n; ++s) {
    if (seen[s]) continue;
    // Follow file order when the successor s+1 is a neighbour.
    std::size_t prev = s;
    std::size_t cur = adj[s][1] == s + 1 ? adj[s][1] : adj[s][0];
    Ring<double> r{verts[s]};
    seen[s] = true;
    while (cur != s) {
      if (seen[cur]) fail(path, 0, "malformed edge cycle");
      seen[cur] = true;
      r.push_back(verts[cur]);
      const std::size_t nxt = adj[cur][0] == prev ? adj[cur][1] : adj[cur][0];
      prev = cur;
      cur = nxt;
    }
    rings.push_back(std::move(r));
  }
  return rings;
}

Ring<double> read_obj_ring(const std::string& path, ObjPlane plane)
{
  std::vector<Ring<double>> rings = read_obj_rings(path, plane);
  if (rings.size() != 1)
    fail(path, 0, "expected one ring, found " + std::to_string(rings.size()));
  return std::move(rings.front());
}

} // namespace pd
