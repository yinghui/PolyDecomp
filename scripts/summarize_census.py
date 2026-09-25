#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""Summarises pd_decompose --compare stats files: outcome counts per file.

Categories: agree (same structure as the exact run, coordinates within 1e-9 of
the scale), deviate_valid (different, but every piece simple, convex up to the
angle tolerance, and areas sum correctly), deviate_suspect (a necessary
condition fails: a silently wrong result), reported_failure (the kernel
reported a numerical failure), exact_failed (the exact run failed).
"""

import json
import sys

CATEGORIES = ["agree", "deviate_valid", "deviate_suspect", "reported_failure", "exact_failed"]


def main(paths):
    header = ["family", "inputs"] + CATEGORIES
    rows = []
    for path in paths:
        recs = [json.loads(line) for line in open(path) if line.strip()]
        family = path.rsplit("/", 1)[-1].split(".")[0]
        counts = [sum(1 for r in recs if r["category"] == c) for c in CATEGORIES]
        rows.append([family, len(recs)] + counts)
    totals = ["TOTAL", sum(r[1] for r in rows)] + [sum(r[i] for r in rows) for i in range(2, len(header))]
    widths = [max(len(str(x)) for x in col) for col in zip(header, *rows, totals)]
    for r in [header] + rows + [totals]:
        print("  ".join(str(x).rjust(w) for x, w in zip(r, widths)))


if __name__ == "__main__":
    main(sys.argv[1:])
