#!/usr/bin/env python3
"""
summarize.py — summarize NeoQual pool_compare results (out/<k>.json).

Per metric: in how many samples go-py is the largest of py-py / go-go / go-py
(about 1/3 expected if Go behaves like one more Python pool), and the median of
go-py / mean(py-py, go-go). Mean K Go - Python: sign test across samples, valid
here because every sample has its own Go seeds. Breakdown by size class. Same for mean tree depth (deepest SSM-carrying
node). Runtime (out/runtime_<k>.json from runtime.py): Go/Python wall and CPU
ratios by size class. Writes <work>/summary_table.tsv, one row per sample keyed by k.
"""

import argparse
import collections
import glob
import json
import statistics as st
from math import comb
from pathlib import Path

M = ["coclust", "ancestry", "phi", "K_tv", "Kcnv_tv", "Kbig_tv", "branching", "depth_tv"]


def sign_p(w, n):
    lo = min(w, n - w)
    return min(1.0, 2 * sum(comb(n, i) for i in range(lo + 1)) / 2 ** n) if n else 1.0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("work", help="work dir with manifest.json and out/")
    args = ap.parse_args()
    W = Path(args.work)
    man = json.load(open(W / "manifest.json"))
    rows = [(man[int(Path(f).stem)], json.load(open(f))) for f in glob.glob(str(W / "out/[0-9]*.json"))]
    print(f"samples compared: {len(rows)} of {len(man)} in manifest")
    if not rows:
        raise SystemExit("no out/*.json results")
    print("Go chains used:", dict(collections.Counter(d["n_go_chains"] for _, d in rows)))

    top, n, ratio, flags = collections.Counter(), collections.Counter(), collections.defaultdict(list), 0
    for e, d in rows:
        g = d["groups"]
        for x in M:
            if x not in g["go-py"]:
                continue
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
        if not r:
            print(f"{x:<10} n/a (no results carry it)")
            continue
        print(f"{x:<10} {top[x]:>5}/{n[x]:<6}   {st.median(r):.2f} ({r[len(r)//4]:.2f}-{r[3*len(r)//4]:.2f})")
    print(f"sample-metrics with go-py > 2 SD beyond both baselines: {flags}")

    print("\nmean K, Go - Python (sign test across samples):")
    for kk in ("K", "K_cnv", "K_big", "depth"):
        for cls in ["all"] + sorted({e["cls"] for e, _ in rows}):
            v = [d["mean_K"][kk]["go"] - d["mean_K"][kk]["py"] for e, d in rows
                 if cls in ("all", e["cls"]) and kk in d["mean_K"]]
            if not v:
                continue
            w, l = sum(x > 0 for x in v), sum(x < 0 for x in v)
            print(f"  {kk:<6} {cls:<4} n={len(v):<5} Go>Py {w}/{w + l}  p={sign_p(w, w + l):.3g}  "
                  f"median {st.median(v):+.3f}  mean {st.mean(v):+.3f} ± {st.stdev(v) / len(v) ** .5 if len(v) > 1 else 0:.3f}")

    runtime(W, man, {e["k"]: d for e, d in rows})


def runtime(W, man, cmp):
    rt = {int(Path(f).stem.split("_")[1]): json.load(open(f)) for f in glob.glob(str(W / "out/runtime_*.json"))}
    print(f"\nruntime records: {len(rt)};  go_status:", dict(collections.Counter(r["go_status"] for r in rt.values())))
    print("by go_src:", dict(collections.Counter(r["go_src"] for r in rt.values())))
    print("cls  n_ok  median Py wall h  median Go wall h  Go/Py wall (IQR)       Go/Py CPU (IQR)")
    for cls in ["all", "S", "M", "L", "X"]:
        ok = [r for r in rt.values() if r["go_status"] == "ok" and cls in ("all", r["cls"])]
        if not ok:
            continue
        q = lambda v: (st.median(v), sorted(v)[len(v) // 4], sorted(v)[3 * len(v) // 4])
        w, c = q([r["wall_ratio_go_py"] for r in ok]), q([r["cpu_ratio_go_py"] for r in ok])
        print(f"{cls:<4} {len(ok):<5} {st.median([r['py_wall_s'] for r in ok]) / 3600:>10.2f} "
              f"{st.median([r['go_wall_s'] for r in ok]) / 3600:>17.2f}   {w[0]:.3f} ({w[1]:.3f}-{w[2]:.3f})"
              f"   {c[0]:.3f} ({c[1]:.3f}-{c[2]:.3f})")
    cols = ["k", "cls", "n_ssm", "go_src", "go_status", "py_wall_s", "py_cpu_s", "go_wall_s", "go_cpu_s",
            "go_cores", "wall_ratio_go_py", "cpu_ratio_go_py", "n_go_chains", "py_mean_depth", "go_mean_depth",
            "py_max_depth", "go_max_depth", "py_mean_K", "go_mean_K"] + [f"{x}_gp_over_within" for x in M]
    n = 0
    with open(W / "summary_table.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for e in man:
            r, d = rt.get(e["k"], {}), cmp.get(e["k"])
            if not r and d is None:
                continue
            row = {**{c: e.get(c) for c in ("k", "cls", "n_ssm", "go_src")}, **r}
            if d is not None:
                row["n_go_chains"] = d["n_go_chains"]
                for eng in ("py", "go"):
                    row[f"{eng}_mean_K"] = d["mean_K"]["K"][eng]
                    row[f"{eng}_mean_depth"] = d["mean_K"].get("depth", {}).get(eng)
                    row[f"{eng}_max_depth"] = d.get("max_depth", {}).get(eng)
                g = d["groups"]
                for x in M:
                    if x in g["go-py"]:
                        within = (g["py-py"][x]["mean"] + g["go-go"][x]["mean"]) / 2
                        row[f"{x}_gp_over_within"] = g["go-py"][x]["mean"] / within if within else None
            f.write("\t".join("" if row.get(c) is None else
                               (f"{row[c]:.4g}" if isinstance(row[c], float) else str(row[c])) for c in cols) + "\n")
            n += 1
    print(f"wrote {n} rows to {W / 'summary_table.tsv'}")


if __name__ == "__main__":
    main()
