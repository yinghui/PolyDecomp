#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""Independent GEOS oracle for polydecomp outputs (shares no code with the library).

usage: oracle.py INPUTS.wkt PIECES.wkt [--rational PIECES.rat] [--convex] [--tol T]
                 [--jsonl OUT] [--workers N]

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

Pieces are rounded to doubles, so tiny residuals are expected. Rounding can
also make a piece invalid: two exact vertices closer than a unit in the last
place, such as the two crossings of a cut that passes next to a hole corner,
round to one point, and the ring touches itself. With --rational (the exact
pieces from pd_decompose --rational, one block per decomposed input, in the
order of PIECES.wkt), a piece that GEOS finds invalid is checked again in exact
rational arithmetic (Python fractions); if its exact ring is simple, it counts
as "invalid only after rounding" instead of as a failure.
"""

import argparse
import json
import sys
from fractions import Fraction
from multiprocessing import Pool

import shapely
from shapely import wkt
from shapely.ops import unary_union
from shapely.strtree import STRtree

ARGS = None


def orient(a, b, c):
    d = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    return (d > 0) - (d < 0)


def on_segment(a, b, p):
    return min(a[0], b[0]) <= p[0] <= max(a[0], b[0]) and min(a[1], b[1]) <= p[1] <= max(a[1], b[1])


def segments_meet(a, b, c, d):
    o1, o2, o3, o4 = orient(a, b, c), orient(a, b, d), orient(c, d, a), orient(c, d, b)
    if o1 != o2 and o3 != o4 and o1 * o2 <= 0 and o3 * o4 <= 0:
        return True
    return ((o1 == 0 and on_segment(a, b, c)) or (o2 == 0 and on_segment(a, b, d))
            or (o3 == 0 and on_segment(c, d, a)) or (o4 == 0 and on_segment(c, d, b)))


def exact_simple(ring):
    """True if the ring (list of Fraction pairs, not closed) is a simple polygon."""
    n = len(ring)
    if n < 3 or len(set(ring)) != n:
        return False
    edges = [(ring[i], ring[(i + 1) % n]) for i in range(n)]
    for i in range(n):
        a, b = edges[i]
        c = edges[(i + 1) % n][1]
        if orient(a, b, c) == 0 and not on_segment(a, c, b):
            return False            # the next edge folds back over this one
        for j in range(i + 2, n):
            if i == 0 and j == n - 1:
                continue            # adjacent through the closing vertex
            if segments_meet(a, b, *edges[j]):
                return False
    return True


def read_rational(path):
    """Blocks of exact rings, one per decomposed input, as lists of Fraction pairs."""
    tok = open(path).read().split()
    blocks, k = [], 0
    while k < len(tok):
        if tok[k] != "polydecomp-rational" or tok[k + 1] != "1" or tok[k + 2] != "rings":
            sys.exit(f"{path}: bad header at token {k}")
        nrings = int(tok[k + 3])
        k += 4
        rings = []
        for _ in range(nrings):
            n = int(tok[k])
            k += 1
            rings.append([(Fraction(tok[k + 2 * v]), Fraction(tok[k + 2 * v + 1])) for v in range(n)])
            k += 2 * n
        blocks.append(rings)
    return blocks


def check(item):
    idx, in_text, out_text, exact = item
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
    if exact is not None and len(exact) != len(pieces):
        sys.exit(f"case {idx}: {len(pieces)} pieces but {len(exact)} exact rings")
    rec["valid"] = True
    rec["rounding_invalid"] = 0
    for k, p in enumerate(pieces):
        if p.area <= tol * area or (p.is_valid and len(p.interiors) == 0):
            continue
        if exact is not None and len(p.interiors) == 0 and exact_simple(exact[k]):
            rec["rounding_invalid"] += 1
        else:
            rec["valid"] = False
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
    ap.add_argument("--rational")
    ap.add_argument("--convex", action="store_true")
    ap.add_argument("--tol", type=float, default=1e-9)
    ap.add_argument("--jsonl")
    ap.add_argument("--workers", type=int, default=16)
    ARGS = ap.parse_args()

    ins = [l for l in open(ARGS.inputs) if l.strip()]
    outs = [l for l in open(ARGS.pieces) if l.strip()]
    if len(ins) != len(outs):
        sys.exit(f"{len(ins)} inputs but {len(outs)} outputs")
    exact = [None] * len(outs)
    if ARGS.rational:
        blocks = iter(read_rational(ARGS.rational))
        for i, b in enumerate(outs):
            if not wkt.loads(b).is_empty:
                exact[i] = next(blocks, None)
                if exact[i] is None:
                    sys.exit(f"{ARGS.rational}: fewer blocks than decomposed inputs")
        if next(blocks, None) is not None:
            sys.exit(f"{ARGS.rational}: more blocks than decomposed inputs")
    with Pool(ARGS.workers) as pool:
        recs = pool.map(check, [(i, a, b, e) for i, (a, b, e) in enumerate(zip(ins, outs, exact))],
                        chunksize=16)
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
          + (f"pieces invalid only after rounding {sum(r['rounding_invalid'] for r in done)} "
             f"in {sum(1 for r in done if r['rounding_invalid'])} cases; " if ARGS.rational else "")
          + f"negligible pieces {sum(r['negligible'] for r in done)}; max inside {worst('inside'):.1e}, "
          f"max overlap {worst('overlap'):.1e}, max area {worst('area'):.1e}"
          + (f", max convex excess {worst('convex_excess'):.1e}" if ARGS.convex else "")
          + f"; GEOS floating-overlay anomalies: {anomalies}")
    for r in bad[:5]:
        print("  fail:", json.dumps(r))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
