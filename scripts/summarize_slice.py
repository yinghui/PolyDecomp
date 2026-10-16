#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""Summarizes the slicing experiment (pd_bench --suite slice, SLICE.csv) as Markdown tables: one row
per shape, cut sequence and method (ours, clipper2, cgal-boolean).

usage: summarize_slice.py SLICE.csv
"""
import csv
import math
import statistics
import sys


def pct(values, q):
    """Nearest-rank percentile."""
    s = sorted(values)
    return s[max(0, min(len(s) - 1, math.ceil(q / 100 * len(s)) - 1))]


def main():
    rows = list(csv.DictReader(open(sys.argv[1])))
    complete = [r for r in rows if None not in r.values()]
    if len(complete) < len(rows):
        print(f"warning: skipped {len(rows) - len(complete)} incomplete rows", file=sys.stderr)
    groups = {}
    for r in complete:
        groups.setdefault((r["shape"], r["sequence"], r["method"]), []).append(r)
    for (shape, seq, method), g in groups.items():
        g.sort(key=lambda r: int(r["cut"]))

    print("## Ours (pd::Slicer)\n")
    print("| shape | n | h | sequence | cuts | pieces | cut pieces | vertex contacts | edge contacts "
          "| dangles | symbolic from cut | cut ms median | cut ms p99 | convex ms median | convex ms p99 "
          "| max bits (cut 10 / end) | certified |")
    print("|" + "---|" * 17)
    for (shape, seq, method), g in groups.items():
        if method != "ours":
            continue
        cut = [float(r["cut_ms"]) for r in g]
        conv = [float(r["convex_ms"]) for r in g]
        sym = [int(r["cut"]) for r in g if r["symbolic"] == "1"]
        total = lambda k: sum(int(r[k]) for r in g)
        print(f"| {shape} | {g[0]['n']} | {g[0]['holes']} | {seq} | {len(g)} | {g[-1]['pieces']} "
              f"| {total('cut_pieces')} | {total('vertex_contacts')} | {total('edge_contacts')} "
              f"| {total('dangles')} | {sym[0] if sym else '-'} "
              f"| {statistics.median(cut):.2f} | {pct(cut, 99):.2f} | {statistics.median(conv):.2f} "
              f"| {pct(conv, 99):.2f} | {g[min(9, len(g) - 1)]['max_bits']} / {g[-1]['max_bits']} "
              f"| {sum(r['certified'] == '1' for r in g)}/{len(g)} |")

    print("\n## Baselines on the same lines\n")
    print("| shape | sequence | method | cuts | pieces (ours) | pieces | pieces as the slicer counts them "
          "| tiny | area error | certified | cut ms median | cut ms p99 | total s (ours) | total s |")
    print("|" + "---|" * 14)
    for (shape, seq, method), g in groups.items():
        if method == "ours":
            continue
        ours = groups[(shape, seq, "ours")]
        k = len(g)
        cut = [float(r["cut_ms"]) for r in g]
        ours_total = sum(float(r["cut_ms"]) for r in ours[:k]) / 1000
        checked = [r["certified"] for r in g if r["certified"]]
        cert = f"{checked.count('1')}/{len(checked)}" if checked else "-"
        print(f"| {shape} | {seq} | {method} | {k} | {ours[k - 1]['pieces']} | {g[-1]['pieces']} "
              f"| {g[-1]['split_pieces'] or '-'} | {g[-1]['tiny'] or '-'} | {g[-1]['area_err'] or '-'} | {cert} "
              f"| {statistics.median(cut):.2f} | {pct(cut, 99):.2f} | {ours_total:.2f} | {sum(cut) / 1000:.2f} |")


if __name__ == "__main__":
    main()
