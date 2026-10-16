#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""Tables from a pd_bench CSV: piece counts and median times per input group and method.

usage: summarize_bench.py BENCH.csv

Inputs are grouped by name (star-512-2 -> star-512, row-7 -> row, comb-1024,
spiral-64, legacy-*). For simple polygons the piece count is also given relative
to Chazelle's lower bound ceil(r/2) + 1 (summed over the group).
"""

import csv
import re
import statistics
import sys
from collections import OrderedDict, defaultdict


def group_of(suite, name):
    """Simple suite: star-<n>-<seed> -> star-<n>, comb-<n> kept. Holes suite:
    spiral-<h> kept, <family>-<seed> -> <family>. Legacy inputs and every input of the
    scale suite (star-<n>, comb-<n>, holes-<h>) are their own group."""
    if name.startswith("legacy-") or suite in ("scale", "hole-scale"):
        return name
    if suite == "simple":
        m = re.match(r"^(star|comb)-(\d+)", name)
        return f"{m.group(1)}-{m.group(2)}" if m else name
    if name.startswith("spiral-"):
        return name
    return re.sub(r"-\d+$", "", name)


def main(path):
    rows = list(csv.DictReader(open(path)))
    methods = list(OrderedDict.fromkeys(r["method"] for r in rows))
    groups = OrderedDict()
    for r in rows:
        g = (r["suite"], group_of(r["suite"], r["input"]))
        groups.setdefault(g, []).append(r)

    failed = [r for r in rows if r["status"] != "ok"]
    print(f"rows {len(rows)}, failed runs {len(failed)}")
    for suite in ("simple", "holes", "scale", "hole-scale"):
        keys = [g for g in groups if g[0] == suite]
        if not keys:
            continue
        print(f"\n## {suite}: total pieces per group" + (" (/ lower bound)" if suite == "simple" else ""))
        cols = [m for m in methods if any(r["method"] == m for g in keys for r in groups[g])]
        print("| group | inputs | n | " + ("lower bound | " if suite == "simple" else "holes | ") + " | ".join(cols) + " |")
        print("|" + "---|" * (4 + len(cols)))
        for g in keys:
            rs = groups[g]
            inputs = sorted({r["input"] for r in rs})
            per = defaultdict(int)
            present = set()
            for r in rs:
                per[r["method"]] += int(r["pieces"])
                present.add(r["method"])
            one = [r for r in rs if r["method"] == rs[0]["method"]]
            n = round(statistics.mean(int(r["n"]) for r in one))
            if suite == "simple":
                lb = sum(int(r["lower_bound"]) for r in one)
                cells = [f"{per[m]} ({per[m] / lb:.2f})" if m in present else "-" for m in cols]
                extra = str(lb)
            else:
                cells = [str(per[m]) if m in present else "-" for m in cols]
                extra = str(round(statistics.mean(int(r["holes"]) for r in one)))
            print(f"| {g[1]} | {len(inputs)} | {n} | {extra} | " + " | ".join(cells) + " |")

        print(f"\n## {suite}: median time per input (ms)")
        print("| group | " + " | ".join(cols) + " |")
        print("|" + "---|" * (1 + len(cols)))
        for g in keys:
            rs = groups[g]
            cells = []
            for m in cols:
                t = [float(r["median_ms"]) for r in rs if r["method"] == m]
                cells.append(f"{statistics.median(t):.3g}" if t else "-")
            print(f"| {g[1]} | " + " | ".join(cells) + " |")


if __name__ == "__main__":
    main(sys.argv[1])
