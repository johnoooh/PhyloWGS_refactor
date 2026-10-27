#!/usr/bin/env python3
"""
mk_manifest.py — manifest for the NeoQual Go-vs-Python posterior comparison.

Takes every sample under <nq>/<cohort>/phylo_checkpoints/<sample>/ that has all
15 Python chains complete (trees.zip present, mcmc_samples.txt with 1000 burn-in
+ 2500 samples + header = 3501 lines) and its own ssm_data.txt/cnv_data.txt,
which are Python's exact inputs. Go runs on those same files.

Each sample k gets its own Go seed, 1_000_000 + 100_000*k (chain c uses
seed + 1000*c, so no two samples share an RNG stream). Python chains already
converted by the earlier 74-sample run are reused (--reuse manifest_all.json).

Writes <out>/manifest.json (list of dicts, index = k) and
<out>/tasks_<cls>.txt (one k per line) for size classes by SSM count:
S < 300, M < 1000, L < 3000, X >= 3000.
"""

import argparse
import collections
import json
import os
from pathlib import Path

SEED_BASE, SEED_STEP, N_CHAINS, MCMC_LINES = 1_000_000, 100_000, 15, 3501
CLASSES = (("S", 300), ("M", 1000), ("L", 3000), ("X", float("inf")))


def n_lines(p):
    with open(p, "rb") as f:
        return sum(1 for _ in f)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--nq", required=True, help="NeoQual_test root")
    ap.add_argument("--out", required=True, help="work dir for this run")
    ap.add_argument("--reuse", help="earlier manifest_all.json whose py_dir conversions to reuse")
    args = ap.parse_args()
    out = Path(args.out).resolve()

    reuse = {}
    if args.reuse:
        for e in json.load(open(args.reuse)):
            reuse[os.path.abspath(Path(e["py_src"][0]).parent.parent.parent)] = e["py_dir"]

    man, why = [], collections.Counter()
    for ck in sorted(Path(os.path.abspath(args.nq)).glob("*/phylo_checkpoints/*")):
        if not ck.is_dir():
            continue
        chains = [ck / f"chains/chain_{c}" for c in range(N_CHAINS)]
        if not all((c / "trees.zip").exists() for c in chains):
            why["missing_py_chain"] += 1
            continue
        if any(n_lines(c / "mcmc_samples.txt") != MCMC_LINES for c in chains):
            why["py_chain_incomplete"] += 1
            continue
        ssm, cnv = ck / "ssm_data.txt", ck / "cnv_data.txt"
        if not (ssm.exists() and cnv.exists()):
            why["missing_inputs"] += 1
            continue
        n_ssm = n_lines(ssm) - 1
        if n_ssm < 2:
            why["lt_2_ssms"] += 1
            continue
        k = len(man)
        cls = next(c for c, hi in CLASSES if n_ssm < hi)
        man.append({
            "k": k, "cls": cls, "n_ssm": n_ssm, "cohort": ck.parent.parent.name,
            "seed": SEED_BASE + SEED_STEP * k, "ssm": str(ssm), "cnv": str(cnv),
            "py_src": [str(c / "trees.zip") for c in chains],
            "py_dir": reuse.get(str(ck), str(out / f"py/{k}")),  # ck is absolute, not symlink-resolved
            "go_dir": str(out / f"go/{k}"),
        })
        why["included"] += 1

    if not man:
        raise SystemExit(f"no samples found under {args.nq}: {dict(why)}")
    out.mkdir(parents=True, exist_ok=True)
    json.dump(man, open(out / "manifest.json", "w"), indent=1)
    for c, _ in CLASSES:
        ks = [str(e["k"]) for e in man if e["cls"] == c]
        (out / f"tasks_{c}.txt").write_text("\n".join(ks) + ("\n" if ks else ""))
    print(dict(why))
    print("by class:", dict(collections.Counter(e["cls"] for e in man)))
    print("by cohort:", dict(collections.Counter(e["cohort"] for e in man)))
    print("python conversions reused:", sum(not e["py_dir"].startswith(str(out)) for e in man))


if __name__ == "__main__":
    main()
