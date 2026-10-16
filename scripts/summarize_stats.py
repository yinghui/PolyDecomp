#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""Summarises pd_decompose --stats files (JSON lines), one table row per file."""

import json
import statistics
import sys


def summarize(path):
    rows = [json.loads(line) for line in open(path) if line.strip()]
    ok = [r for r in rows if r["status"] == "ok"]
    certified = sum(1 for r in rows if r.get("verified") is True)
    not_certified = sum(1 for r in rows if r.get("verified") is False)

    def total(key):
        return sum(r.get(key, 0) for r in ok)

    def most(key):
        return max((r.get(key, 0) for r in ok), default=0)

    times = [r["time_ms"] for r in ok]
    return {
        "file": path.rsplit("/", 1)[-1],
        "inputs": len(rows),
        "ok": len(ok),
        "certified": certified,
        "not_certified": not_certified,
        "holes": total("holes"),
        "cuts": total("cuts"),
        "hole_free": total("hole_free_pieces"),
        "convex": total("convex_pieces"),
        "reflex": total("reflex_vertices"),
        "vertex_hits": total("vertex_hits"),
        "double_res": total("double_resolutions"),
        "v_contacts": total("vertex_contacts"),
        "e_contacts": total("edge_contacts"),
        "dangles": total("dangles"),
        "symbolic": total("symbolic"),
        "fb_rotated": total("fallback_rotated"),
        "fb_midpoint": total("fallback_midpoint"),
        "max_depth": max(most("max_depth_holes"), most("max_depth")),
        "max_bits": most("max_bits"),
        "median_ms": statistics.median(times) if times else 0.0,
        "max_ms": max(times, default=0.0),
    }


def main(paths):
    rows = [summarize(p) for p in paths]
    keys = list(rows[0].keys())
    widths = {k: max(len(k), *(len(f"{r[k]:.2f}" if isinstance(r[k], float) else str(r[k])) for r in rows))
              for k in keys}
    fmt = lambda v: f"{v:.2f}" if isinstance(v, float) else str(v)
    print("  ".join(k.rjust(widths[k]) for k in keys))
    for r in rows:
        print("  ".join(fmt(r[k]).rjust(widths[k]) for k in keys))
    all_ok = sum(r["ok"] for r in rows)
    all_in = sum(r["inputs"] for r in rows)
    all_cert = sum(r["certified"] for r in rows)
    print(f"TOTAL inputs {all_in}, ok {all_ok}, certified {all_cert}, "
          f"not certified {sum(r['not_certified'] for r in rows)}")


if __name__ == "__main__":
    main(sys.argv[1:])
