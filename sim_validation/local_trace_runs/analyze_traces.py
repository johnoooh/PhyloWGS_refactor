#!/usr/bin/env python3
"""
Analyze local Go trace runs against Python reference results.

For each fixture, computes:
  - Per-iter spawn dynamics (spawns_in_find_node, spawns_in_stick_orders,
    cap_hits_find_node) — burnin vs sample phase
  - K_before_cull and K_after_cull distributions per phase
  - Hyperparameter trajectories (dp_alpha, dp_gamma, alpha_decay)
  - LLH trajectory and best-llh
  - Final tree count distribution from chain_*_trees.ndjson
  - Side-by-side with Python reference K distribution from
    simulation_validation/results/original-python/<fixture>/tree_summaries.json.gz

Writes a markdown report to local_trace_runs/REPORT.md
"""

import gzip
import json
import statistics
from collections import Counter
from pathlib import Path

REPO_ROOT = Path("/sessions/dreamy-pensive-knuth/mnt/phylowgs_cowork")
RUNS_ROOT = REPO_ROOT / "PhyloWGS_refactor/sim_validation/local_trace_runs"
PY_REF_ROOT = REPO_ROOT / "simulation_validation/results/original-python"

FIXTURES = [
    "K3_S1_T200_M30_C2_rep0",
    "K3_S1_T200_M30_C0_rep0",
    "K5_S1_T200_M50_C2_rep0",
    "K5_S1_T200_M50_C0_rep0",
]

EXPERIMENTS = [
    # (label, description, fixtures, burnin_count)
    ("expA_baseline", "HPC-parity baseline (B=500, S=1000, 1 chain)", FIXTURES, 500),
    ("expB_long",     "Long convergence (B=2000, S=5000, 1 chain)",   FIXTURES, 2000),
    ("expC_multichain","Multi-chain variance K3_C2 (B=500, S=1000, 4 chains)",
        ["K3_S1_T200_M30_C2_rep0"], 500),
]


def load_trace(path: Path):
    """Robust against partial / mid-write trace files: skips malformed final lines."""
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                rows.append(json.loads(line))
            except json.JSONDecodeError:
                # Likely a partial flush from a still-running writer; stop reading.
                break
    return rows


def split_burnin(rows, burnin_count):
    """In main.go:2689 the trace iter is offset by +burnin, so trace iter is always
    in [0, burnin+samples). The first `burnin_count` records are burnin; the rest
    are sample iters."""
    burnin = rows[:burnin_count]
    sample = rows[burnin_count:]
    return burnin, sample


def stats(values):
    if not values:
        return {"n": 0}
    return {
        "n": len(values),
        "min": min(values),
        "max": max(values),
        "mean": round(statistics.mean(values), 3),
        "median": statistics.median(values),
        "stdev": round(statistics.stdev(values), 3) if len(values) > 1 else 0.0,
    }


def summarize_trace(rows, burnin_count):
    if not rows:
        return None
    burnin, sample = split_burnin(rows, burnin_count)
    out = {"total_iters": len(rows), "n_burnin": len(burnin), "n_sample": len(sample)}
    for phase, rs in (("burnin", burnin), ("sample", sample)):
        if not rs:
            continue
        out[phase] = {
            "K_before_cull": stats([r["K_before_cull"] for r in rs]),
            "K_after_cull":  stats([r["K_after_cull"] for r in rs]),
            "spawns_find":   stats([r["spawns_in_find_node"] for r in rs]),
            "spawns_stick":  stats([r["spawns_in_stick_orders"] for r in rs]),
            "cap_hits":      sum(r.get("cap_hits_find_node", 0) for r in rs),
            "kills_in_cull": stats([r["kills_in_cull"] for r in rs]),
            "dp_alpha":      stats([r["dp_alpha"] for r in rs]),
            "dp_gamma":      stats([r["dp_gamma"] for r in rs]),
            "alpha_decay":   stats([r["alpha_decay"] for r in rs]),
            "best_llh":      max(r["best_llh_this_iter"] for r in rs),
        }
    return out


def k_distribution_from_trees(ndjson_path: Path):
    pops = []
    llhs = []
    with open(ndjson_path) as f:
        for line in f:
            d = json.loads(line)
            p = d.get("num_populations") or len(d.get("populations", {}) or {})
            pops.append(p)
            llhs.append(d.get("llh", 0))
    return pops, llhs


def python_reference_distribution(fixture: str):
    p = PY_REF_ROOT / fixture / "tree_summaries.json.gz"
    if not p.exists():
        return None
    with gzip.open(p) as f:
        data = json.load(f)
    trees = data.get("trees", data)
    pops, llhs = [], []
    for k, v in trees.items():
        n = v.get("populations", v.get("num_populations"))
        if isinstance(n, dict):
            n = len(n)
        pops.append(n)
        llhs.append(v.get("llh", v.get("llh_complete", 0)))
    return pops, llhs


def fmt_dist(label, pops, llhs):
    c = Counter(pops)
    top = c.most_common(3)
    bi = max(range(len(llhs)), key=lambda i: llhs[i])
    return (
        f"{label:18s} N={len(pops):5d}  range={min(pops)}-{max(pops):2d}  "
        f"median={statistics.median(pops):>5}  mean={statistics.mean(pops):.2f}  "
        f"mode={top[0]}  best K={pops[bi]} llh={llhs[bi]:.2f}"
    )


def render():
    out = ["# Local trace baseline report", ""]
    out.append(f"Generated from `{RUNS_ROOT.relative_to(REPO_ROOT)}`")
    out.append("")
    out.append("Single-page calibration of Go-port spawn dynamics and K distributions")
    out.append("on a 4-fixture local slice, run BEFORE HPC results land. Used to")
    out.append("interpret tomorrow's HPC numbers.")
    out.append("")

    # Python reference table first
    out.append("## Python reference K distributions (same fixtures)")
    out.append("")
    out.append("```")
    for f in FIXTURES:
        ref = python_reference_distribution(f)
        if ref is None:
            out.append(f"{f}: [no python reference]")
            continue
        pops, llhs = ref
        out.append(fmt_dist(f, pops, llhs))
    out.append("```")
    out.append("")

    for label, desc, fixtures, burnin_count in EXPERIMENTS:
        out.append(f"## {label} — {desc}")
        out.append("")
        for fixture in fixtures:
            fdir = RUNS_ROOT / label / fixture
            if not fdir.exists():
                out.append(f"### {fixture}: [missing run dir]")
                continue
            out.append(f"### {fixture}")
            out.append("")

            # 1. Trace summary (single chain or chain 0)
            trace_files = sorted(fdir.glob("trace.ndjson*"))
            if not trace_files:
                out.append("_(no trace file written)_")
                out.append("")
                continue
            for tf in trace_files:
                rows = load_trace(tf)
                summary = summarize_trace(rows, burnin_count)
                if summary is None:
                    out.append(f"_(empty trace: {tf.name})_")
                    continue
                out.append(f"**trace: `{tf.name}`** total_iters={summary['total_iters']}  "
                           f"burnin={summary['n_burnin']}  sample={summary['n_sample']}")
                out.append("")
                out.append("```")
                for phase in ("burnin", "sample"):
                    if phase not in summary:
                        continue
                    s = summary[phase]
                    out.append(f"  [{phase}]")
                    out.append(f"    K_before_cull   {s['K_before_cull']}")
                    out.append(f"    K_after_cull    {s['K_after_cull']}")
                    out.append(f"    spawns_find     {s['spawns_find']}")
                    out.append(f"    spawns_stick    {s['spawns_stick']}")
                    out.append(f"    kills_in_cull   {s['kills_in_cull']}")
                    out.append(f"    cap_hits_find   {s['cap_hits']}")
                    out.append(f"    dp_alpha        {s['dp_alpha']}")
                    out.append(f"    dp_gamma        {s['dp_gamma']}")
                    out.append(f"    alpha_decay     {s['alpha_decay']}")
                    out.append(f"    best_llh        {s['best_llh']:.2f}")
                out.append("```")
                out.append("")

            # 2. Tree NDJSON distributions vs Python reference
            tree_files = sorted(fdir.glob("chain_*_trees.ndjson"))
            if tree_files:
                go_pops, go_llhs = [], []
                for tf in tree_files:
                    p, l = k_distribution_from_trees(tf)
                    go_pops += p
                    go_llhs += l
                out.append("**posterior K distribution (vs Python reference)**")
                out.append("")
                out.append("```")
                ref = python_reference_distribution(fixture)
                if ref is not None:
                    out.append(fmt_dist("Python (orig)", ref[0], ref[1]))
                out.append(fmt_dist("Go local", go_pops, go_llhs))
                out.append("```")
                out.append("")

    out_path = RUNS_ROOT / "REPORT.md"
    out_path.write_text("\n".join(out))
    print(f"Wrote {out_path}")
    return out_path


if __name__ == "__main__":
    render()
