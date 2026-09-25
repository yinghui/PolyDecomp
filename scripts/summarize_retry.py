#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
# Summarizes pd_decompose --stats files of double-kernel runs with --retry-exact.
#   usage: scripts/summarize_retry.py FAMILY=FILE.jsonl ...
import json
import sys

for arg in sys.argv[1:]:
    family, path = arg.split("=", 1)
    runs = [json.loads(line) for line in open(path) if line.strip()]
    ok = sum(1 for r in runs if r["status"] == "ok")
    retried = sum(r.get("retried_exact", 0) for r in runs)
    print(f"{family}: inputs {len(runs)}, ok {ok}, retried_exact {retried}")
