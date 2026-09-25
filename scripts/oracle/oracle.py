#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""Independent GEOS oracle for polydecomp outputs (shares no code with the library).

usage: oracle.py INPUTS.wkt PIECES.wkt [--convex] [--tol T] [--jsonl OUT] [--workers N]

INPUTS.wkt holds one POLYGON per line; PIECES.wkt the corresponding line of
pd_decompose --out (one MULTIPOLYGON per line, or MULTIPOLYGON EMPTY for an
input that was not decomposed, which is reported and skipped).

A case passes if, with A the input area and tolerance tol (default 1e-9):
  valid    every non-negligible piece (area > tol * A) is a valid polygon
           without holes; negligible pieces are counted, not failed (exact
           results can contain pieces too small to survive rounding to doubles)
  inside   sum over pieces of area(piece - input) / A        <= tol
  overlap  sum over piece pairs of area(intersection) / A    <= tol
  area     |sum of piece areas - A| / A                       <= tol
  convex   (with --convex) every non-negligible piece: (hull area - area) / area <= tol
Pieces inside the input, pairwise almost disjoint and of the right total area
cover the input up to the tolerance (inclusion-exclusion), so no global union
is needed.

All overlay operations use GEOS snap-rounding (grid 1e-12 times the
coordinate scale), GEOS's robust fixed-precision mode. Its default
floating-point overlay is not reliable on these outputs, where many pieces meet
along nearly collinear edges with T-junctions: it has returned a union with a
spurious hole and an intersection equal to one operand while the difference
was that operand too. For information, the oracle still evaluates the
floating-point union coverage and pairwise overlap and counts the cases in
which they exceed the tolerance ("floating-overlay anomalies").

Pieces are rounded to doubles, so tiny residuals are expected.
"""

import argparse
import json
import sys
from multiprocessing import Pool

import shapely
from shapely import wkt
from shapely.ops import unary_union
from shapely.strtree import STRtree

ARGS = None


def check(item):
    idx, in_text, out_text = item
    rec = {"index": idx}
    poly = wkt.loads(in_text)
    out = wkt.loads(out_text)
    if out.is_empty:
        rec["skipped"] = True
        return rec
    pieces = list(out.geoms)
    area = poly.area
    tol = ARGS.tol
    xs, ys = poly.exterior.xy
    grid = 1e-12 * max(1.0, max(map(abs, xs)), max(map(abs, ys)))
    real = [p for p in pieces if p.area > tol * area]
    rec["pieces"] = len(pieces)
    rec["negligible"] = len(pieces) - len(real)
    rec["valid"] = all(p.is_valid and len(p.interiors) == 0 for p in real)
    rec["inside"] = sum(shapely.difference(p, poly, grid_size=grid).area for p in pieces) / area
    tree = STRtree(pieces)
    overlap = float_overlap = 0.0
    for i, p in enumerate(pieces):
        for j in tree.query(p):
            if j > i:
                overlap += shapely.intersection(p, pieces[j], grid_size=grid).area
                float_overlap += p.intersection(pieces[j]).area
    rec["overlap"] = overlap / area
    rec["area"] = abs(sum(p.area for p in pieces) - area) / area
    if ARGS.convex:
        rec["convex_excess"] = max(((p.convex_hull.area - p.area) / p.area for p in real), default=0.0)
    rec["float_overlap"] = float_overlap / area
    rec["float_union_coverage"] = unary_union(pieces).symmetric_difference(poly).area / area
    rec["ok"] = (rec["valid"] and rec["inside"] <= tol and rec["overlap"] <= tol and rec["area"] <= tol
                 and (not ARGS.convex or rec["convex_excess"] <= tol))
    return rec


def main():
    global ARGS
    ap = argparse.ArgumentParser()
    ap.add_argument("inputs")
    ap.add_argument("pieces")
    ap.add_argument("--convex", action="store_true")
    ap.add_argument("--tol", type=float, default=1e-9)
    ap.add_argument("--jsonl")
    ap.add_argument("--workers", type=int, default=16)
    ARGS = ap.parse_args()

    ins = [l for l in open(ARGS.inputs) if l.strip()]
    outs = [l for l in open(ARGS.pieces) if l.strip()]
    if len(ins) != len(outs):
        sys.exit(f"{len(ins)} inputs but {len(outs)} outputs")
    with Pool(ARGS.workers) as pool:
        recs = pool.map(check, [(i, a, b) for i, (a, b) in enumerate(zip(ins, outs))], chunksize=16)
    if ARGS.jsonl:
        with open(ARGS.jsonl, "w") as f:
            for r in recs:
                f.write(json.dumps(r) + "\n")

    done = [r for r in recs if not r.get("skipped")]
    bad = [r for r in done if not r["ok"]]

    def worst(key):
        vals = [r.get(key, 0.0) for r in done]
        return max(vals, default=0.0)

    anomalies = sum(1 for r in done if r["float_union_coverage"] > ARGS.tol or r["float_overlap"] > ARGS.tol)
    print(f"GEOS {shapely.geos_version_string}: cases {len(recs)}, checked {len(done)}, "
          f"skipped {len(recs) - len(done)}, pass {len(done) - len(bad)}, fail {len(bad)} "
          f"(invalid {sum(1 for r in done if not r['valid'])}); "
          f"negligible pieces {sum(r['negligible'] for r in done)}; max inside {worst('inside'):.1e}, "
          f"max overlap {worst('overlap'):.1e}, max area {worst('area'):.1e}"
          + (f", max convex excess {worst('convex_excess'):.1e}" if ARGS.convex else "")
          + f"; GEOS floating-overlay anomalies: {anomalies}")
    for r in bad[:5]:
        print("  fail:", json.dumps(r))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
