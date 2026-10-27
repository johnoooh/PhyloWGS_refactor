#!/usr/bin/env python3
"""
summarize.py — summarize NeoQual pool_compare results (out/<k>.json).

Per metric: in how many samples go-py is the largest of py-py / go-go / go-py
(about 1/3 expected if Go behaves like one more Python pool), and the median of
go-py / mean(py-py, go-go). Mean K Go - Python: sign test across samples, valid
here because every sample has its own Go seeds. Breakdown by size class.
"""

import argparse
import collections
import glob
import json
import statistics as st
from math import comb
from pathlib import Path

M = ["coclust", "ancestry", "phi", "K_tv", "Kcnv_tv", "Kbig_tv", "branching"]


def sign_p(w, n):
    lo = min(w, n - w)
    return min(1.0, 2 * sum(comb(n, i) for i in range(lo + 1)) / 2 ** n) if n else 1.0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("work", help="work dir with manifest.json and out/")
    args = ap.parse_args()
    W = Path(args.work)
    man = json.load(open(W / "manifest.json"))
    rows = [(man[int(Path(f).stem)], json.load(open(f))) for f in glob.glob(str(W / "out/*.json"))]
    print(f"samples compared: {len(rows)} of {len(man)} in manifest")
    if not rows:
        raise SystemExit("no out/*.json results")
    print("Go chains used:", dict(collections.Counter(d["n_go_chains"] for _, d in rows)))

    top, n, ratio, flags = collections.Counter(), collections.Counter(), collections.defaultdict(list), 0
    for e, d in rows:
        g = d["groups"]
        for x in M:
            pp, gg, gp = (g[q][x]["mean"] for q in ("py-py", "go-go", "go-py"))
            if max(pp, gg, gp) == 0:
                continue
            n[x] += 1
            if (pp + gg) > 0:
                ratio[x].append(gp / ((pp + gg) / 2))
            if gp > max(pp, gg):
                top[x] += 1
                sd = max(g["py-py"][x]["sd"], g["go-go"][x]["sd"]) or 1e-9
                flags += (gp - max(pp, gg)) / sd > 2
    print("\nmetric     go-py largest   median go-py/within (IQR)")
    for x in M:
        r = sorted(ratio[x])
        print(f"{x:<10} {top[x]:>5}/{n[x]:<6}   {st.median(r):.2f} ({r[len(r)//4]:.2f}-{r[3*len(r)//4]:.2f})")
    print(f"sample-metrics with go-py > 2 SD beyond both baselines: {flags}")

    print("\nmean K, Go - Python (sign test across samples):")
    for kk in ("K", "K_cnv", "K_big"):
        for cls in ["all"] + sorted({e["cls"] for e, _ in rows}):
            v = [d["mean_K"][kk]["go"] - d["mean_K"][kk]["py"] for e, d in rows if cls in ("all", e["cls"])]
            w, l = sum(x > 0 for x in v), sum(x < 0 for x in v)
            print(f"  {kk:<6} {cls:<4} n={len(v):<5} Go>Py {w}/{w + l}  p={sign_p(w, w + l):.3g}  "
                  f"median {st.median(v):+.3f}  mean {st.mean(v):+.3f} ± {st.stdev(v) / len(v) ** .5 if len(v) > 1 else 0:.3f}")


if __name__ == "__main__":
    main()
