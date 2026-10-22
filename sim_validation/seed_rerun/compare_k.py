#!/usr/bin/env python3
"""
compare_k.py — per-sample mean K: per-sample-seed Go rerun vs Python vs the test5 Go batches.

For each rerun sample, mean K (populations with >= 1 SSM, raw trees) is averaged
per chain, then over the 15 chains. Python and the test5 Go seed batches (101/102/103)
come from the per-chain table kstats/<k>.json; batches are labelled from the
manifest's chain paths (<cls>_<idx> = 101, <cls>102_<idx>, <cls>103_<idx>).

Prediction if the K lean was the shared-seed artefact: rerun Go vs Python shows no
lean (Go > Py in about half of samples, mean diff near 0).

Note: the SE for the test5 batches is not valid (their chain streams are shared
across samples); it is valid for the rerun, whose streams are all distinct.
"""

import argparse
import json
import re
import statistics as st
from math import comb
from pathlib import Path

N_CHAINS, N_TREES = 15, 2500


def chain_mean_k(path):
    ks = [sum(p["num_ssms"] > 0 for p in json.loads(ln)["populations"].values())
          for ln in open(path) if ln.strip()]
    if len(ks) != N_TREES:
        raise ValueError(f"{path}: {len(ks)} trees, expected {N_TREES}")
    return sum(ks) / len(ks)


def sign_p(wins, n):
    """Two-sided exact sign test."""
    lo = min(wins, n - wins)
    return min(1.0, 2 * sum(comb(n, i) for i in range(lo + 1)) / 2 ** n)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--rerun", required=True, help="rerun dir (tasks.tsv, runs/k<k>/chains)")
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--kstats", required=True, help="dir of per-chain kstats/<k>.json")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    man = {e["k"]: e for e in json.load(open(args.manifest))}
    rows = []
    for ln in open(Path(args.rerun) / "tasks.tsv"):
        _, k, cls, n_ssm, seed = ln.split("\t")[:5]
        k = int(k)
        chains = sorted((Path(args.rerun) / f"runs/k{k}/chains").glob("chain_*_trees.ndjson"))
        if len(chains) != N_CHAINS:
            raise ValueError(f"k={k}: {len(chains)} rerun chains, expected {N_CHAINS}")
        ks = json.load(open(Path(args.kstats) / f"{k}.json"))
        batch = [re.search(r"/runs/[A-D](\d{3})?_\d{3}/", p).group(1) or "101"
                 for p in man[k]["go"]]
        if len(batch) != len(ks["go"]) or len(ks["py"]) != N_CHAINS:
            raise ValueError(f"k={k}: kstats/manifest chain counts disagree")
        row = {"k": k, "cls": cls, "n_ssm": int(n_ssm),
               "rerun": st.mean(chain_mean_k(c) for c in chains),
               "py": st.mean(c["meanK"] for c in ks["py"])}
        for b in ("101", "102", "103"):
            v = [c["meanK"] for c, bb in zip(ks["go"], batch) if bb == b]
            if len(v) != N_CHAINS:
                raise ValueError(f"k={k}: batch {b} has {len(v)} chains")
            row[b] = st.mean(v)
        rows.append(row)

    print(f"{'k':>4} {'cls':>3} {'n_ssm':>6} {'rerun':>7} {'py':>7} {'101':>7} {'102':>7} {'103':>7}")
    for r in rows:
        print(f"{r['k']:>4} {r['cls']:>3} {r['n_ssm']:>6} " + " ".join(
            f"{r[c]:>7.2f}" for c in ("rerun", "py", "101", "102", "103")))

    n = len(rows)
    summary = {}
    print(f"\nGo set vs Python over {n} samples (ties dropped from sign test):")
    for g in ("rerun", "101", "102", "103"):
        d = [r[g] - r["py"] for r in rows]
        wins, losses = sum(x > 0 for x in d), sum(x < 0 for x in d)
        s = {"go_gt_py": wins, "go_lt_py": losses, "sign_p": sign_p(wins, wins + losses),
             "median_diff": st.median(d), "mean_diff": st.mean(d),
             "se": st.stdev(d) / n ** 0.5}
        summary[g] = s
        print(f"  {g:>5}: Go>Py {wins}/{wins + losses}  p={s['sign_p']:.3g}  "
              f"median {s['median_diff']:+.3f}  mean {s['mean_diff']:+.3f} ± {s['se']:.3f}"
              + ("" if g == "rerun" else "  (SE invalid: shared streams)"))
    json.dump({"rows": rows, "summary": summary}, open(args.out, "w"), indent=1)


if __name__ == "__main__":
    main()
