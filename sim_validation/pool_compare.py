#!/usr/bin/env python3
"""
pool_compare.py — Go vs Python posterior comparison on real data, by chain pools.

On real data each chain freezes into one mode early, so a single chain is
close to a point mass and the posterior is the mixture over chains. This
compares pooled posteriors of equal-size, disjoint, random chain subsets:

  py-py  two disjoint pools of Python chains   (Python sampling noise)
  go-go  two disjoint pools of Go chains        (Go sampling noise)
  go-py  a Go pool vs a Python pool             (implementation difference)

Metrics are those of compare_posteriors.py (coclust, ancestry, phi, K_tv,
branching). Each chain is summarized once; a pool is the tree-count-weighted
mean of its chains (chains are thinned to the same tree count).

Inputs per chain:
  Go:     chain_<k>_trees.ndjson (one snapshot per line, post-burn-in)
  Python: a directory holding write_results.py output for that chain's
          trees.zip (tree_summaries.json.gz + mutass.zip), made with
          --include-multiprimary --keep-superclones --min-ssms 0

Usage:
    python pool_compare.py --go run/chains/chain_*_trees.ndjson \
        --py pyconv/chain_* --pool 7 --draws 20 --out sample.json
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_posteriors import _iter_py_trees, distances, summarize_run  # noqa: E402

METRICS = ["coclust", "ancestry", "phi", "K_tv", "Kcnv_tv", "Kbig_tv", "branching"]


def load_go_chain(path):
    with open(path) as f:
        trees = [(t["structure"], t["populations"], t["mut_assignments"])
                 for t in map(json.loads, f) if t]
    if not trees:
        raise ValueError(f"{path}: no trees")
    return trees


def load_py_chain(path):
    trees = list(_iter_py_trees(Path(path)))
    if not trees:
        raise ValueError(f"{path}: no trees")
    return trees


def ssm_set(trees):
    return {s for _, _, ma in trees for m in ma.values() for s in m.get("ssms", [])}


def pool(summaries):
    return {
        "C": np.mean([s["C"] for s in summaries], axis=0),
        "A": np.mean([s["A"] for s in summaries], axis=0),
        "phi": np.mean([s["phi"] for s in summaries], axis=0),
        "K": np.concatenate([s["K"] for s in summaries]),
        "K_cnv": np.concatenate([s["K_cnv"] for s in summaries]),
        "K_big": np.concatenate([s["K_big"] for s in summaries]),
        "branching": float(np.mean([s["branching"] for s in summaries])),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--go", nargs="+", required=True, help="Go chain ndjson files")
    ap.add_argument("--py", nargs="+", required=True, help="Python per-chain result dirs")
    ap.add_argument("--pool", type=int, default=7, help="Chains per pool")
    ap.add_argument("--draws", type=int, default=20, help="Random pool draws per group")
    ap.add_argument("--trees-per-chain", type=int, default=250,
                    help="Thin every chain to this many trees")
    ap.add_argument("--max-ssms", type=int, default=1500,
                    help="Score a fixed random subset of this many SSMs (bounds M x M memory;"
                         " pair-averaged metrics stay unbiased)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    for impl, n in (("go", len(args.go)), ("py", len(args.py))):
        if n < 2 * args.pool:
            raise ValueError(f"{n} {impl} chains < 2 x pool ({args.pool})")

    chains = [("go", p, load_go_chain(p)) for p in args.go] + \
             [("py", p, load_py_chain(p)) for p in args.py]
    for impl, p, trees in chains:
        if len(trees) < args.trees_per_chain:
            raise ValueError(f"{p}: {len(trees)} trees < --trees-per-chain")
    sets = [ssm_set(t) for _, _, t in chains]
    if any(s != sets[0] for s in sets):
        raise ValueError("Chains do not share the same SSM set")
    ssm_ids = sorted(sets[0], key=lambda s: (len(s), s))
    n_ssm_total = len(ssm_ids)
    if n_ssm_total > args.max_ssms:
        keep = np.random.default_rng(args.seed).choice(n_ssm_total, args.max_ssms, replace=False)
        ssm_ids = [ssm_ids[i] for i in sorted(keep)]

    summ = {"go": [], "py": []}
    for impl, p, trees in chains:
        summ[impl].append(summarize_run(trees, ssm_ids, args.trees_per_chain))
    del chains

    rng = np.random.default_rng(args.seed)
    P = args.pool
    rows = []
    for _ in range(args.draws):
        g = rng.permutation(len(summ["go"]))
        p = rng.permutation(len(summ["py"]))
        go_a, go_b = pool([summ["go"][i] for i in g[:P]]), pool([summ["go"][i] for i in g[P:2 * P]])
        py_a, py_b = pool([summ["py"][i] for i in p[:P]]), pool([summ["py"][i] for i in p[P:2 * P]])
        rows.append(("py-py", distances(py_a, py_b)))
        rows.append(("go-go", distances(go_a, go_b)))
        rows.append(("go-py", distances(go_a, py_a)))

    out = {"M": n_ssm_total, "M_scored": len(ssm_ids), "n_go_chains": len(summ["go"]), "n_py_chains": len(summ["py"]),
           "pool": P, "draws": args.draws,
           "mean_K": {kk: {k: float(np.mean([s[kk].mean() for s in v])) for k, v in summ.items()}
                      for kk in ("K", "K_cnv", "K_big")},
           "groups": {}}
    for grp in ["py-py", "go-go", "go-py"]:
        out["groups"][grp] = {
            m: {"mean": float(np.mean([d[m] for gg, d in rows if gg == grp])),
                "sd": float(np.std([d[m] for gg, d in rows if gg == grp]))}
            for m in METRICS}
    with open(args.out, "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out["groups"]), file=sys.stderr)


if __name__ == "__main__":
    main()
