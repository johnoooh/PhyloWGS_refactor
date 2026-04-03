#!/usr/bin/env python3
"""
compare_implementations.py — Compare Go port vs original Python PhyloWGS results.

Reads scores.json from both implementations and produces:
  - comparison.tsv: side-by-side metrics per fixture
  - comparison_plots/: runtime speedup, K_error comparison, AUPRC comparison

Usage:
    python compare_implementations.py \
        --go-scores analysis/go-cpu/scores.json \
        --py-scores analysis/original-python/scores.json \
        --outdir analysis/comparison
"""

import argparse
import json
import os
import sys

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    print("ERROR: matplotlib required. Install with: pip install matplotlib")
    sys.exit(1)

plt.rcParams.update({
    "font.size": 11, "figure.dpi": 150, "savefig.dpi": 150,
    "savefig.bbox": "tight", "figure.facecolor": "white",
})


def load_scores(path):
    with open(path) as f:
        data = json.load(f)
    return {s["fixture"]: s for s in data["scores"]}


def main():
    parser = argparse.ArgumentParser(description="Compare Go vs Python PhyloWGS results")
    parser.add_argument("--go-scores", required=True)
    parser.add_argument("--py-scores", required=True)
    parser.add_argument("--outdir", default="analysis/comparison")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    go = load_scores(args.go_scores)
    py = load_scores(args.py_scores)

    common = sorted(set(go.keys()) & set(py.keys()))
    if not common:
        print("No common fixtures between Go and Python results.")
        return

    print(f"Comparing {len(common)} common fixtures")

    # ── Build comparison table ───────────────────────────────────────────
    rows = []
    for name in common:
        g, p = go[name], py[name]
        row = {
            "fixture": name,
            "K": g["K"], "S": g["S"], "T": g["T"], "M": g["M"],
            "go_time_s": g.get("total_time_s"),
            "py_time_s": p.get("total_time_s"),
            "go_K_error": g.get("K_error"),
            "py_K_error": p.get("K_error"),
            "go_inferred_K": g.get("inferred_K"),
            "py_inferred_K": p.get("inferred_K"),
            "go_best_llh": g.get("best_llh"),
            "py_best_llh": p.get("best_llh"),
            "go_auprc": g.get("cocluster_auprc"),
            "py_auprc": p.get("cocluster_auprc"),
        }
        if row["go_time_s"] and row["py_time_s"] and row["py_time_s"] > 0:
            row["speedup"] = row["py_time_s"] / row["go_time_s"]
        rows.append(row)

    # Write TSV
    tsv_path = os.path.join(args.outdir, "comparison.tsv")
    cols = list(rows[0].keys())
    with open(tsv_path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for row in rows:
            f.write("\t".join(str(row.get(c, "")) for c in cols) + "\n")
    print(f"  Wrote {tsv_path}")

    # ── Plot 1: Runtime speedup ──────────────────────────────────────────
    speedup_rows = [r for r in rows if r.get("speedup")]
    if speedup_rows:
        fig, ax = plt.subplots(figsize=(8, 5))
        Ms = [r["M"] for r in speedup_rows]
        speedups = [r["speedup"] for r in speedup_rows]
        Ks = [r["K"] for r in speedup_rows]
        unique_Ks = sorted(set(Ks))
        colors = plt.cm.Set2(np.linspace(0, 0.8, max(len(unique_Ks), 1)))
        K_color = {k: colors[i] for i, k in enumerate(unique_Ks)}

        for K in unique_Ks:
            subset_M = [m for m, k in zip(Ms, Ks) if k == K]
            subset_s = [s for s, k in zip(speedups, Ks) if k == K]
            ax.scatter(subset_M, subset_s, c=[K_color[K]], s=60, alpha=0.7,
                      label=f"K={K}", edgecolors="white")

        ax.axhline(y=1.0, color="red", linestyle="--", alpha=0.5, label="Parity")
        ax.set_xlabel("Number of SSMs (M)")
        ax.set_ylabel("Speedup (Python time / Go time)")
        ax.set_title("Go Port Speedup vs Original Python")
        ax.legend(title="Clones")
        mean_speedup = np.mean(speedups)
        ax.text(0.95, 0.05, f"Mean speedup: {mean_speedup:.1f}x",
               transform=ax.transAxes, ha="right", fontsize=11,
               bbox=dict(boxstyle="round", facecolor="#e8f5e9"))
        fig.savefig(os.path.join(args.outdir, "speedup.png"))
        plt.close(fig)
        print(f"  [saved] speedup.png (mean={mean_speedup:.1f}x)")

    # ── Plot 2: K error comparison ───────────────────────────────────────
    k_rows = [r for r in rows if r.get("go_K_error") is not None and r.get("py_K_error") is not None]
    if k_rows:
        fig, ax = plt.subplots(figsize=(6, 6))
        go_errs = [r["go_K_error"] for r in k_rows]
        py_errs = [r["py_K_error"] for r in k_rows]
        ax.scatter(py_errs, go_errs, c="#42a5f5", s=50, alpha=0.7, edgecolors="white")
        max_err = max(max(go_errs), max(py_errs), 1)
        ax.plot([0, max_err], [0, max_err], "r--", alpha=0.5, label="Equal")
        ax.set_xlabel("Python |K error|")
        ax.set_ylabel("Go |K error|")
        ax.set_title("Population Count Error: Go vs Python")
        ax.legend()
        fig.savefig(os.path.join(args.outdir, "k_error_comparison.png"))
        plt.close(fig)
        print("  [saved] k_error_comparison.png")

    # ── Plot 3: AUPRC comparison ─────────────────────────────────────────
    auprc_rows = [r for r in rows if r.get("go_auprc") is not None and r.get("py_auprc") is not None]
    if auprc_rows:
        fig, ax = plt.subplots(figsize=(6, 6))
        go_a = [r["go_auprc"] for r in auprc_rows]
        py_a = [r["py_auprc"] for r in auprc_rows]
        ax.scatter(py_a, go_a, c="#ab47bc", s=50, alpha=0.7, edgecolors="white")
        ax.plot([0, 1], [0, 1], "r--", alpha=0.5, label="Equal")
        ax.set_xlabel("Python AUPRC")
        ax.set_ylabel("Go AUPRC")
        ax.set_title("Co-clustering Accuracy: Go vs Python")
        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(-0.05, 1.05)
        ax.legend()
        fig.savefig(os.path.join(args.outdir, "auprc_comparison.png"))
        plt.close(fig)
        print("  [saved] auprc_comparison.png")

    # ── Plot 4: Summary dashboard ────────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Panel A: Runtime bars
    ax = axes[0]
    time_rows = [r for r in rows if r.get("go_time_s") and r.get("py_time_s")]
    if time_rows:
        go_times = [r["go_time_s"] for r in time_rows]
        py_times = [r["py_time_s"] for r in time_rows]
        ax.bar(["Go (mean)", "Python (mean)"],
               [np.mean(go_times), np.mean(py_times)],
               color=["#42a5f5", "#ef5350"], edgecolor="white")
        ax.set_ylabel("Wall Time (s)")
    ax.set_title("Mean Runtime")

    # Panel B: K error bars
    ax = axes[1]
    if k_rows:
        go_k = [r["go_K_error"] for r in k_rows]
        py_k = [r["py_K_error"] for r in k_rows]
        ax.bar(["Go (mean)", "Python (mean)"],
               [np.mean(go_k), np.mean(py_k)],
               color=["#42a5f5", "#ef5350"], edgecolor="white")
        ax.set_ylabel("Mean |K error|")
    ax.set_title("Population Count Error")

    # Panel C: AUPRC bars
    ax = axes[2]
    if auprc_rows:
        go_au = [r["go_auprc"] for r in auprc_rows]
        py_au = [r["py_auprc"] for r in auprc_rows]
        ax.bar(["Go (mean)", "Python (mean)"],
               [np.mean(go_au), np.mean(py_au)],
               color=["#42a5f5", "#ef5350"], edgecolor="white")
        ax.set_ylabel("Mean AUPRC")
        ax.set_ylim(0, 1.05)
    ax.set_title("Co-clustering Accuracy")

    fig.suptitle("Go Port vs Original Python — Simulation Validation", fontsize=14)
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "comparison_dashboard.png"))
    plt.close(fig)
    print("  [saved] comparison_dashboard.png")

    # ── Print summary ────────────────────────────────────────────────────
    print(f"\n{'='*50}")
    print("Go vs Python Comparison Summary")
    print(f"{'='*50}")
    print(f"Common fixtures: {len(common)}")
    if speedup_rows:
        print(f"Mean speedup:    {np.mean([r['speedup'] for r in speedup_rows]):.1f}x")
    if k_rows:
        print(f"Go mean K_err:   {np.mean([r['go_K_error'] for r in k_rows]):.2f}")
        print(f"Py mean K_err:   {np.mean([r['py_K_error'] for r in k_rows]):.2f}")
    if auprc_rows:
        print(f"Go mean AUPRC:   {np.mean([r['go_auprc'] for r in auprc_rows]):.3f}")
        print(f"Py mean AUPRC:   {np.mean([r['py_auprc'] for r in auprc_rows]):.3f}")


if __name__ == "__main__":
    main()
