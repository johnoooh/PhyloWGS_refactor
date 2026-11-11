#!/usr/bin/env python3
"""
mk_manifest.py — manifest for the NeoQual Go-vs-Python posterior comparison.

Takes every sample under <nq>/<cohort>/phylo_checkpoints/<sample>/ that has all
15 Python chains complete (trees.zip present, mcmc_samples.txt with 1000 burn-in
+ 2500 samples + header = 3501 lines) and its own ssm_data.txt/cnv_data.txt,
which are Python's exact inputs. Go runs on those same files.

Each new Go run for sample k gets its own seed, 100_000_000 + 100_000*k (chain c
uses seed + 1000*c), a range disjoint from the per-sample seeds of the earlier
seed_rerun (1_000_000 + 100_000*k_old), so no two samples share an RNG stream.

Reuse (no rerun of what was already run correctly):
  --reuse     earlier manifest_all.json: its Python conversions (py_dir)
  --reuse-go  an earlier per-sample-seed Go run dir (seed_rerun: tasks.tsv +
              runs/k<k_old>/chains, k_old indexing --reuse). A sample's Go chains
              are reused only if all 15 chain_<c>_samples.txt have 2501 lines, all
              15 trees.ndjson and summary.json exist, and that run's ssm/cnv inputs
              are byte-identical to the checkpoint's. go_job is the latest COMPLETED
              sacct task of that run (read by runtime.py).
Every entry records go_src = new | reused:<run dir name>.

Writes <out>/manifest.json (list of dicts, index = k) and chunked task lists
<out>/tasks_<cls><i>.txt (one k per line). Chunk 0 of each class is the pilot
(--pilot-n per class: up to half reused-Go samples, the rest new); later chunks
hold <= --chunk samples. Size classes by SSM count: S < 300, M < 1000,
L < 3000, X >= 3000.
"""

import argparse
import collections
import filecmp
import json
import os
import subprocess
from pathlib import Path

SEED_BASE, SEED_STEP, N_CHAINS, MCMC_LINES, GO_LINES = 100_000_000, 100_000, 15, 3501, 2501
CLASSES = (("S", 300), ("M", 1000), ("L", 3000), ("X", float("inf")))


def n_lines(p):
    with open(p, "rb") as f:
        return sum(1 for _ in f)


def ckpt_of(py_src0):
    """Checkpoint sample dir from a chain's trees.zip path."""
    return os.path.abspath(Path(py_src0).parent.parent.parent)


def go_incomplete(chains):
    """Why a Go chains dir is not reusable, or None if it is."""
    if not (chains / "summary.json").exists():
        return "no_summary"
    for c in range(N_CHAINS):
        s, t = chains / f"chain_{c}_samples.txt", chains / f"chain_{c}_trees.ndjson"
        if not (s.exists() and t.exists()):
            return "missing_chain"
        if n_lines(s) != GO_LINES:
            return "short_chain"
    return None


def sacct_tasks(job_name):
    """{array index: latest COMPLETED 'jobid_index'} for a job name."""
    out = subprocess.run(["sacct", "-n", "-P", "-X", "-S", "2026-01-01", "--name", job_name,
                          "-o", "JobID,State,Start"], capture_output=True, text=True, check=True).stdout
    best = {}
    for ln in out.splitlines():
        jid, state, start = ln.split("|")
        if state != "COMPLETED" or "_" not in jid:
            continue
        i = int(jid.split("_")[1])
        if i not in best or start > best[i][1]:
            best[i] = (jid, start)
    if not best:
        raise SystemExit(f"sacct found no COMPLETED tasks named {job_name}")
    return {i: j for i, (j, _) in best.items()}


def load_reuse_go(run_dir, manifest_all, job_name):
    """{ckpt dir: (chains dir, seed, go_job)} for reusable samples of an earlier Go run."""
    run_dir = Path(run_dir).resolve()
    by_k = {e["k"]: ckpt_of(e["py_src"][0]) for e in manifest_all}
    jobs = sacct_tasks(job_name)
    ok, why = {}, collections.Counter()
    rows = [ln.split("\t") for ln in open(run_dir / "tasks.tsv").read().splitlines() if ln.strip()]
    if not rows:
        raise SystemExit(f"{run_dir}/tasks.tsv has no rows")
    for i, k_old, _cls, _n, seed, ssm, cnv in rows:
        ck = by_k[int(k_old)]
        chains = run_dir / f"runs/k{k_old}/chains"
        bad = go_incomplete(chains)
        if bad is None and not (filecmp.cmp(ssm, f"{ck}/ssm_data.txt", shallow=False)
                                and filecmp.cmp(cnv, f"{ck}/cnv_data.txt", shallow=False)):
            bad = "inputs_differ"
        if bad is None and int(i) not in jobs:
            bad = "no_completed_job"
        if bad:
            why[bad] += 1
            continue
        ok[ck] = (str(chains), int(seed), jobs[int(i)])
        why["reusable"] += 1
    print(f"reuse-go {run_dir.name}: {len(rows)} rows, {dict(why)}")
    return ok, run_dir.name


def chunks_for(entries, n_pilot, size):
    """Pilot chunk (up to half reused, rest new) then <= size chunks of the rest, by k."""
    reused = [e["k"] for e in entries if e["go_src"] != "new"]
    new = [e["k"] for e in entries if e["go_src"] == "new"]
    p = reused[:max(1, n_pilot // 2)] if n_pilot else []
    p += new[:n_pilot - len(p)]
    rest = sorted(set(reused + new) - set(p))
    return [p] + [rest[i:i + size] for i in range(0, len(rest), size)]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--nq", required=True, help="NeoQual_test root")
    ap.add_argument("--out", required=True, help="work dir for this run")
    ap.add_argument("--reuse", help="earlier manifest_all.json whose py_dir conversions to reuse")
    ap.add_argument("--reuse-go", help="earlier per-sample-seed Go run dir (needs --reuse)")
    ap.add_argument("--reuse-go-job", default="go_seedrerun", help="its Slurm job name")
    ap.add_argument("--chunk", type=int, default=300, help="max samples per array chunk")
    ap.add_argument("--pilot-n", default="S6,M4,L3,X1", help="chunk-0 (pilot) size per class")
    args = ap.parse_args()
    out = Path(args.out).resolve()

    reuse, reuse_go, go_run = {}, {}, None
    if args.reuse:
        man_all = json.load(open(args.reuse))
        reuse = {ckpt_of(e["py_src"][0]): e["py_dir"] for e in man_all}
        if args.reuse_go:
            reuse_go, go_run = load_reuse_go(args.reuse_go, man_all, args.reuse_go_job)
    elif args.reuse_go:
        raise SystemExit("--reuse-go needs --reuse (it maps the old run's k to samples)")

    man, why = [], collections.Counter()
    for ck in sorted(Path(os.path.abspath(args.nq)).glob("*/phylo_checkpoints/*")):
        if not ck.is_dir():
            continue
        why["checkpoint_dirs"] += 1
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
        e = {"k": k, "cls": cls, "n_ssm": n_ssm, "cohort": ck.parent.parent.name,
             "ssm": str(ssm), "cnv": str(cnv),
             "py_src": [str(c / "trees.zip") for c in chains],
             "py_dir": reuse.get(str(ck), str(out / f"py/{k}")),  # ck is absolute, not symlink-resolved
             "seed": SEED_BASE + SEED_STEP * k, "go_dir": str(out / f"go/{k}"),
             "go_src": "new", "go_job": None}
        if str(ck) in reuse_go:
            chains_dir, seed, job = reuse_go[str(ck)]
            e.update(seed=seed, go_dir=str(Path(chains_dir).parent), go_src=f"reused:{go_run}", go_job=job)
        man.append(e)
        why["included"] += 1

    if not man:
        raise SystemExit(f"no samples found under {args.nq}: {dict(why)}")
    n_go_reused = sum(e["go_src"] != "new" for e in man)
    if n_go_reused != len(reuse_go):
        raise SystemExit(f"{len(reuse_go)} reusable Go samples but {n_go_reused} matched a manifest entry")
    out.mkdir(parents=True, exist_ok=True)
    json.dump(man, open(out / "manifest.json", "w"), indent=1)

    pilot = {p[0]: int(p[1:]) for p in args.pilot_n.split(",")}
    for c, _ in CLASSES:
        chunks = chunks_for([e for e in man if e["cls"] == c], pilot.get(c, 0), args.chunk)
        for i, ch in enumerate(chunks):
            if ch:
                (out / f"tasks_{c}{i}.txt").write_text("\n".join(map(str, ch)) + "\n")
        print(f"class {c}: chunk sizes {[len(ch) for ch in chunks]}")
    print(dict(why))
    print("by cohort:", dict(collections.Counter(e["cohort"] for e in man)))
    print("python conversions reused:", sum(not e["py_dir"].startswith(str(out)) for e in man))
    print("go chains reused:", n_go_reused,
          "| new go by class:", dict(collections.Counter(e["cls"] for e in man if e["go_src"] == "new")))


if __name__ == "__main__":
    main()
