#!/usr/bin/env python3
"""Per-iter trajectory exploration of K_before/after, spawn dynamics, and
hyperparameter correlations in the local trace runs. Reports trajectory
markers (early burnin, mid burnin, late burnin, early sample, late sample)
and Spearman-style rank correlations between K_after_cull and key hypers."""

import json
import statistics
from pathlib import Path

RUNS = Path("/sessions/dreamy-pensive-knuth/mnt/phylowgs_cowork/PhyloWGS_refactor/sim_validation/local_trace_runs")

FIXTURES = [
    "K3_S1_T200_M30_C2_rep0",
    "K3_S1_T200_M30_C0_rep0",
    "K5_S1_T200_M50_C2_rep0",
    "K5_S1_T200_M50_C0_rep0",
]


def load(path):
    rows = []
    with open(path) as f:
        for L in f:
            L = L.strip()
            if not L:
                continue
            try:
                rows.append(json.loads(L))
            except json.JSONDecodeError:
                break
    return rows


def windowed(rows, label, n_burnin):
    """Return slices for early/mid/late burnin and early/late sample."""
    burnin = rows[:n_burnin]
    sample = rows[n_burnin:]
    if not burnin or not sample:
        return None
    bn = len(burnin)
    sn = len(sample)
    return {
        "early_burnin": burnin[: max(1, bn // 5)],
        "late_burnin":  burnin[-max(1, bn // 5):],
        "early_sample": sample[: max(1, sn // 5)],
        "late_sample":  sample[-max(1, sn // 5):],
    }


def windowstat(window, key):
    if not window:
        return "—"
    vals = [r[key] for r in window]
    return f"{statistics.mean(vals):>6.1f}±{(statistics.stdev(vals) if len(vals)>1 else 0):>4.1f}"


def spearman(a, b):
    """Manual Spearman rank correlation, no scipy."""
    if len(a) != len(b) or len(a) < 3:
        return 0.0
    def rank(xs):
        # Simple average-rank
        sorted_idx = sorted(range(len(xs)), key=lambda i: xs[i])
        ranks = [0.0] * len(xs)
        for r, i in enumerate(sorted_idx):
            ranks[i] = r
        return ranks
    ra, rb = rank(a), rank(b)
    n = len(a)
    mean_a = sum(ra) / n
    mean_b = sum(rb) / n
    num = sum((ra[i] - mean_a) * (rb[i] - mean_b) for i in range(n))
    da = sum((ra[i] - mean_a) ** 2 for i in range(n)) ** 0.5
    db = sum((rb[i] - mean_b) ** 2 for i in range(n)) ** 0.5
    if da == 0 or db == 0:
        return 0.0
    return num / (da * db)


def analyze_trace(path: Path, n_burnin: int):
    rows = load(path)
    if not rows:
        return
    print(f"\n  ── {path.relative_to(RUNS)} ──  total={len(rows)}  burnin={n_burnin}")
    win = windowed(rows, "", n_burnin)
    if win is None:
        print("    (insufficient data)")
        return

    print("    K_before_cull   K_after_cull    spawns_find    spawns_stick    dp_alpha       alpha_decay     llh")
    for label, w in win.items():
        cells = [
            windowstat(w, "K_before_cull"),
            windowstat(w, "K_after_cull"),
            windowstat(w, "spawns_in_find_node"),
            windowstat(w, "spawns_in_stick_orders"),
            windowstat(w, "dp_alpha"),
            windowstat(w, "alpha_decay"),
            windowstat(w, "best_llh_this_iter"),
        ]
        print(f"    {label:13s}  " + "   ".join(cells))

    # Spearman correlations across the full sample phase
    sample = rows[n_burnin:]
    if len(sample) >= 20:
        K  = [r["K_after_cull"] for r in sample]
        Kb = [r["K_before_cull"] for r in sample]
        sf = [r["spawns_in_find_node"] for r in sample]
        ss = [r["spawns_in_stick_orders"] for r in sample]
        dp = [r["dp_alpha"] for r in sample]
        ad = [r["alpha_decay"] for r in sample]
        dg = [r["dp_gamma"] for r in sample]
        print("    SAMPLE-PHASE Spearman ρ vs K_after_cull:")
        print(f"        K_before_cull: {spearman(K, Kb):+.3f}")
        print(f"        spawns_find:   {spearman(K, sf):+.3f}")
        print(f"        spawns_stick:  {spearman(K, ss):+.3f}")
        print(f"        dp_alpha:      {spearman(K, dp):+.3f}")
        print(f"        alpha_decay:   {spearman(K, ad):+.3f}")
        print(f"        dp_gamma:      {spearman(K, dg):+.3f}")


def main():
    print("=" * 95)
    print("Exp A — HPC parity (B=500/S=1000)")
    print("=" * 95)
    for f in FIXTURES:
        p = RUNS / "expA_baseline" / f / "trace.ndjson"
        if p.exists():
            analyze_trace(p, 500)
    print()
    print("=" * 95)
    print("Exp B — long convergence (B=2000/S=5000) [if complete]")
    print("=" * 95)
    for f in FIXTURES:
        p = RUNS / "expB_long" / f / "trace.ndjson"
        if p.exists():
            analyze_trace(p, 2000)
    print()
    print("=" * 95)
    print("Exp C — multi-chain K3_C2")
    print("=" * 95)
    for tf in sorted((RUNS / "expC_multichain" / "K3_S1_T200_M30_C2_rep0").glob("trace.ndjson*")):
        analyze_trace(tf, 500)


if __name__ == "__main__":
    main()
