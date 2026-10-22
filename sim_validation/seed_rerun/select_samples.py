#!/usr/bin/env python3
"""
select_samples.py — pick samples for the per-sample-seed Go rerun and write its task table.

The test5 Go runs used -r 101/102/103 for every sample, and Go seeds chain c with
seed + 1000*c, so every sample replayed the same 15 RNG streams. This rerun gives
each sample its own seed, 1_000_000 + 100_000*k, so no two samples (and none of
the 101-103 batches) share a chain stream.

Picks N samples evenly spaced by SSM count among those with all 45 test5 Go chains,
and looks up each sample's inputs in the test5 class table of its seed-101 run.

Writes <out>/tasks.tsv: i, k, cls, n_ssm, seed, ssm_path, cnv_path
"""

import argparse
import json
import re
from pathlib import Path

SEED_BASE, SEED_STEP, CHAINS = 1_000_000, 100_000, 15


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--manifest", required=True, help="posterior_cmp manifest_all.json")
    ap.add_argument("--test5", required=True, help="test5 dir with class_<cls>.tsv")
    ap.add_argument("--n", type=int, default=20)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    assert CHAINS * 1000 < SEED_STEP, "chain streams of neighbouring samples would overlap"

    complete = sorted((e for e in json.load(open(args.manifest)) if len(e["go"]) == 45),
                      key=lambda e: e["n_ssm"])
    if len(complete) < args.n:
        raise ValueError(f"only {len(complete)} complete samples < --n {args.n}")
    picks = [complete[round(i * (len(complete) - 1) / (args.n - 1))] for i in range(args.n)]
    if len({e["k"] for e in picks}) != args.n:
        raise ValueError("duplicate picks")

    rows = []
    for i, e in enumerate(picks):
        # seed-101 run dir is <cls>_<idx> (seed-102/103 dirs are <cls>102_<idx>)
        run = next(m.group(1, 2) for p in e["go"]
                   if (m := re.search(r"/runs/([A-D])_(\d{3})/", p)))
        table = Path(args.test5) / f"class_{run[0]}.tsv"
        hit = [ln.rstrip("\n").split("\t") for ln in open(table)
               if ln.split("\t", 1)[0] == str(int(run[1]))]
        if len(hit) != 1:
            raise ValueError(f"k={e['k']}: {len(hit)} rows for idx {run[1]} in {table.name}")
        _, ssm, cnv, _ = hit[0]
        rows.append([i, e["k"], e["cls"], e["n_ssm"], SEED_BASE + SEED_STEP * e["k"], ssm, cnv])

    Path(args.out).mkdir(parents=True, exist_ok=True)
    with open(Path(args.out) / "tasks.tsv", "w") as f:
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")
    print("picked k:", [r[1] for r in rows])
    print("classes:", {c: sum(r[2] == c for r in rows) for c in "ABCD"})


if __name__ == "__main__":
    main()
