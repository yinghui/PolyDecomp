#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""Summary of run_ablation.sh: per cut rule and input set, totals over the inputs.

usage: summarize_ablation.py DIR   (files RULE.SET.jsonl written by pd_decompose --stats)

Sets: `generated` (the five generated families together, also listed one by one) and
`holes-64`, `holes-256`, `holes-1024` (random_holes, five seeds each).
"""

import json
import sys
from pathlib import Path

RULES = ["farthest", "closest", "boundary", "horizontal"]
FAMILIES = ["row", "star", "c", "cluster", "grid"]
SETS = [("generated", FAMILIES)] + [(f, [f]) for f in FAMILIES] + \
       [(f"holes-{h}", [f"holes-{h}"]) for h in (64, 256, 1024)]


def main(d):
    d = Path(d)
    cols = ["rule", "set", "inputs", "ok", "verified", "holes", "cuts", "max_depth", "hole_free",
            "convex", "time_ms"]
    print(" ".join(f"{c:>10}" for c in cols))
    for rule in RULES:
        for name, parts in SETS:
            rows = [json.loads(line) for p in parts for line in open(d / f"{rule}.{p}.jsonl") if line.strip()]
            vals = [rule, name, len(rows), sum(r["status"] == "ok" for r in rows),
                    sum(bool(r.get("verified")) for r in rows), sum(r["holes"] for r in rows),
                    sum(r["cuts"] for r in rows), max(r["max_depth_holes"] for r in rows),
                    sum(r["hole_free_pieces"] for r in rows), sum(r["convex_pieces"] for r in rows),
                    f"{sum(r['time_ms'] for r in rows):.0f}"]
            print(" ".join(f"{v:>10}" for v in vals))


if __name__ == "__main__":
    main(sys.argv[1])
