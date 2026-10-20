#!/usr/bin/env python3
"""
compare_posteriors.py — Compare Go vs Python PhyloWGS posterior tree distributions.

Compares the posterior distributions, not single best trees. Two runs of the
same sampler with different seeds will not pick the same best tree. Every
metric is label-invariant, averaged over posterior trees, and computed per SSM:

  coclust   P(i and j in the same node)
  ancestry  P(node(i) is a strict ancestor of node(j))
  phi       posterior mean cellular prevalence of node(i), per sample
  K         distribution of the number of non-empty populations
  branching fraction of trees in which some node has >= 2 non-empty children

For each metric, the distance between every pair of runs is grouped into
py-py, go-go (seed noise) and go-py (implementation difference). Python is the
reference: each run's mean distance to the Python runs is listed, and a Go run
is flagged when it sits well outside the Python-to-Python spread. Needs >= 2
runs per implementation (different seeds); 3+ Python runs make a usable baseline.

Run inputs:
  Go:     <dir>/trees.zip (JSON snapshots, written by phylowgs-go)
  Python: <dir>/tree_summaries.json.gz + <dir>/mutass.zip (from write_results.py;
          the raw pickled trees.zip from multievolve.py is not read)

CNVs are ignored; only SSM placement is compared.

Usage:
    python compare_posteriors.py --go runs/go_s1 runs/go_s2 runs/go_s3 \
        --py runs/py_s1 runs/py_s2 runs/py_s3 --out comparison.tsv
"""

import argparse
import gzip
import itertools
import json
import sys
import zipfile
from pathlib import Path

import numpy as np

# Arbitrary thresholds for the per-run '*' flag; read the numbers, not just the flag.
FLAG_FACTOR = 2.0
FLAG_FLOOR = 0.01


def _iter_go_trees(run_dir):
    with zipfile.ZipFile(run_dir / "trees.zip") as zf:
        for name in zf.namelist():
            if name.startswith("tree_"):
                t = json.loads(zf.read(name))
                yield t["structure"], t["populations"], t["mut_assignments"]


def _iter_py_trees(run_dir):
    with gzip.open(run_dir / "tree_summaries.json.gz", "rt") as f:
        trees = json.load(f)["trees"]
    with zipfile.ZipFile(run_dir / "mutass.zip") as zf:
        for idx, t in trees.items():
            mutass = json.loads(zf.read(f"{idx}.json"))["mut_assignments"]
            yield t["structure"], t["populations"], mutass


def load_trees(run_dir, impl):
    run_dir = Path(run_dir)
    it = _iter_go_trees(run_dir) if impl == "go" else _iter_py_trees(run_dir)
    trees = list(it)
    if not trees:
        raise ValueError(f"{run_dir}: no posterior trees found ({impl})")
    return trees


def summarize_run(trees, ssm_ids, max_trees):
    """Average label-invariant per-SSM features over (a thinned subset of) trees."""
    if len(trees) > max_trees:
        keep = np.linspace(0, len(trees) - 1, max_trees).round().astype(int)
        trees = [trees[i] for i in keep]
    M = len(ssm_ids)
    pos = {s: i for i, s in enumerate(ssm_ids)}
    C = np.zeros((M, M))
    A = np.zeros((M, M))
    phi = None
    Ks, branching = [], []
    for structure, pops, mutass in trees:
        nodes = sorted(pops, key=int)
        nidx = {n: k for k, n in enumerate(nodes)}
        lab = np.full(M, -1)
        for node, muts in mutass.items():
            for s in muts.get("ssms", []):
                lab[pos[s]] = nidx[node]
        if (lab < 0).any():
            raise ValueError(f"{(lab < 0).sum()} SSMs unassigned in a tree")
        # anc[k, l] = 1 if node k is a strict ancestor of node l
        parent = {str(c): p for p, cs in structure.items() for c in cs}
        anc = np.zeros((len(nodes), len(nodes)))
        for n in nodes:
            p = parent.get(n)
            while p is not None:
                if p in nidx:
                    anc[nidx[p], nidx[n]] = 1
                p = parent.get(str(p))
        C += lab[:, None] == lab[None, :]
        A += anc[lab[:, None], lab[None, :]]
        cp = np.array([pops[n]["cellular_prevalence"] for n in nodes])
        phi = cp[lab] if phi is None else phi + cp[lab]
        nonempty = {n for n in nodes if pops[n]["num_ssms"] > 0}
        Ks.append(len(nonempty))
        branching.append(any(
            sum(str(c) in nonempty for c in cs) >= 2 for cs in structure.values()))
    T = len(trees)
    return {"C": C / T, "A": A / T, "phi": phi / T,
            "K": np.array(Ks), "branching": float(np.mean(branching)), "T": T}


def distances(a, b):
    M = a["C"].shape[0]
    offdiag = ~np.eye(M, dtype=bool)
    kmax = max(a["K"].max(), b["K"].max()) + 1
    pa = np.bincount(a["K"], minlength=kmax) / len(a["K"])
    pb = np.bincount(b["K"], minlength=kmax) / len(b["K"])
    return {
        "coclust": np.abs(a["C"] - b["C"])[offdiag].mean(),
        "ancestry": np.abs(a["A"] - b["A"])[offdiag].mean(),
        "phi": np.abs(a["phi"] - b["phi"]).mean(),
        "K_tv": 0.5 * np.abs(pa - pb).sum(),
        "branching": abs(a["branching"] - b["branching"]),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--go", nargs="+", required=True, help="Go run directories")
    ap.add_argument("--py", nargs="+", required=True, help="Python run directories")
    ap.add_argument("--max-trees", type=int, default=1000,
                    help="Evenly thin each run to at most this many trees")
    ap.add_argument("--out", help="Write per-pair distances as TSV")
    args = ap.parse_args()

    runs = [("go", d) for d in args.go] + [("py", d) for d in args.py]
    loaded = [(impl, d, load_trees(d, impl)) for impl, d in runs]

    ssm_sets = [{s for _, _, ma in trees for m in ma.values() for s in m.get("ssms", [])}
                for _, _, trees in loaded]
    if any(s != ssm_sets[0] for s in ssm_sets):
        raise ValueError("Runs do not share the same SSM set; were they run on the same input?")
    ssm_ids = sorted(ssm_sets[0], key=lambda s: (len(s), s))

    summaries = []
    for impl, d, trees in loaded:
        s = summarize_run(trees, ssm_ids, args.max_trees)
        print(f"{impl}  {d}: {len(trees)} trees (used {s['T']}), "
              f"mean K={s['K'].mean():.2f}, branching={s['branching']:.2f}", file=sys.stderr)
        summaries.append((impl, d, s))

    rows = []
    for (ia, da, sa), (ib, db, sb) in itertools.combinations(summaries, 2):
        group = "-".join(sorted([ia, ib]))
        rows.append({"group": group, "run_a": da, "run_b": db, **distances(sa, sb)})

    metrics = ["coclust", "ancestry", "phi", "K_tv", "branching"]
    print(f"\nM={len(ssm_ids)} SSMs. Mean distance per pair group (max in parentheses):")
    print(f"{'metric':<10} {'py-py':>16} {'go-go':>16} {'go-py':>16}")
    for m in metrics:
        cells = []
        for g in ["py-py", "go-go", "go-py"]:
            v = [r[m] for r in rows if r["group"] == g]
            cells.append(f"{np.mean(v):.4f} ({max(v):.4f})" if v else "n/a")
        print(f"{m:<10} {cells[0]:>16} {cells[1]:>16} {cells[2]:>16}")

    # Per-run view with Python as the reference: each run's mean distance to
    # the (other) Python runs. A pooled go-go baseline lets one bad Go run
    # inflate the "noise" and hide itself, so flag Go runs individually:
    # beyond FLAG_FACTOR x the largest Python value and above FLAG_FLOOR.
    def to_py(d):
        pairs = [r for r in rows if d in (r["run_a"], r["run_b"])
                 and (r["run_b"] if r["run_a"] == d else r["run_a"]) in args.py]
        return {m: np.mean([r[m] for r in pairs]) for m in metrics}
    per_run = [(impl, d, to_py(d)) for impl, d, _ in summaries]
    py_max = {m: max(v[m] for impl, _, v in per_run if impl == "py") for m in metrics}
    print(f"\nPer run: mean distance to Python runs ('*' = Go run > {FLAG_FACTOR:g}x "
          f"every Python run and > {FLAG_FLOOR:g})")
    print(f"{'run':<24}" + "".join(f"{m:>12}" for m in metrics))
    for impl, d, v in per_run:
        print(f"{impl + ' ' + Path(d).name:<24}" + "".join(
            f"{v[m]:>11.4f}" + ("*" if impl == "go" and v[m] > max(FLAG_FACTOR * py_max[m], FLAG_FLOOR)
               else " ")
            for m in metrics))
    if len(args.go) < 2 or len(args.py) < 2:
        print("\nWARNING: <2 runs for an implementation; no seed-noise baseline for it.",
              file=sys.stderr)

    if args.out:
        with open(args.out, "w") as f:
            f.write("\t".join(["group", "run_a", "run_b"] + metrics) + "\n")
            for r in rows:
                f.write("\t".join([r["group"], r["run_a"], r["run_b"]]
                                  + [f"{r[m]:.6g}" for m in metrics]) + "\n")


if __name__ == "__main__":
    main()
