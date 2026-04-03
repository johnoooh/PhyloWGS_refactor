#!/usr/bin/env python3
"""
plot_results.py — Visualize simulation validation results.

Reads scores.json from score_results.py and produces publication-quality figures:
  1. Population count accuracy heatmap (K × read depth)
  2. Co-clustering AUPRC by K and read depth
  3. Runtime scaling by fixture complexity
  4. Chain convergence (LLH std across chains)
  5. MCMC trace plots from chain sample files

Usage:
    python plot_results.py --analysis-dir /path/to/analysis \
                           --result-base /path/to/results \
                           --fixture-base /path/to/fixtures
"""

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm
    import matplotlib.gridspec as gridspec
except ImportError:
    print("ERROR: matplotlib required. Install with: pip install matplotlib")
    sys.exit(1)


# ── Style ────────────────────────────────────────────────────────────────────

plt.rcParams.update({
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.labelsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "figure.dpi": 150,
    "savefig.dpi": 150,
    "savefig.bbox": "tight",
    "figure.facecolor": "white",
})


# ── Data loading ─────────────────────────────────────────────────────────────

def load_scores(analysis_dir):
    """Load scores.json from analysis directory."""
    path = Path(analysis_dir) / "scores.json"
    with open(path) as f:
        data = json.load(f)
    return data["scores"], data.get("failures", [])


def load_chain_traces(result_dir):
    """Load MCMC LLH traces from chain files."""
    traces = {}
    result_path = Path(result_dir)
    for chain_file in sorted(result_path.glob("chain_*_samples.txt")):
        chain_id = chain_file.stem.split("_")[1]
        iters, llhs = [], []
        with open(chain_file) as f:
            f.readline()  # header
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    iters.append(int(parts[0]))
                    llhs.append(float(parts[1]))
        traces[chain_id] = (iters, llhs)
    return traces


# ── Plot 1: Population count accuracy heatmap ────────────────────────────────

def plot_k_accuracy(scores, outdir):
    """Heatmap of mean K error by true K and read depth T."""
    # Group by (K, T)
    groups = {}
    for s in scores:
        if "K_error" not in s:
            continue
        key = (s["K"], s["T"])
        groups.setdefault(key, []).append(s["K_error"])

    if not groups:
        print("  [skip] K accuracy — no K_error data")
        return

    Ks = sorted(set(k for k, _ in groups))
    Ts = sorted(set(t for _, t in groups))

    matrix = np.full((len(Ks), len(Ts)), np.nan)
    for i, K in enumerate(Ks):
        for j, T in enumerate(Ts):
            vals = groups.get((K, T), [])
            if vals:
                matrix[i, j] = np.mean(vals)

    fig, ax = plt.subplots(figsize=(max(5, len(Ts) * 1.5), max(4, len(Ks) * 1.2)))

    vmax = max(np.nanmax(matrix), 1)
    im = ax.imshow(matrix, cmap="RdYlGn_r", aspect="auto", vmin=0, vmax=vmax)

    ax.set_xticks(range(len(Ts)))
    ax.set_xticklabels([str(t) for t in Ts])
    ax.set_yticks(range(len(Ks)))
    ax.set_yticklabels([str(k) for k in Ks])
    ax.set_xlabel("Read Depth (T)")
    ax.set_ylabel("True Clones (K)")
    ax.set_title("Population Count Error (|inferred K - true K|)")

    # Annotate cells
    for i in range(len(Ks)):
        for j in range(len(Ts)):
            val = matrix[i, j]
            if not np.isnan(val):
                n = len(groups.get((Ks[i], Ts[j]), []))
                color = "white" if val > vmax * 0.6 else "black"
                ax.text(j, i, f"{val:.1f}\n(n={n})", ha="center", va="center",
                        fontsize=9, color=color)

    plt.colorbar(im, ax=ax, label="Mean |K error|")
    fig.savefig(os.path.join(outdir, "k_accuracy_heatmap.png"))
    plt.close(fig)
    print("  [saved] k_accuracy_heatmap.png")


# ── Plot 2: Co-clustering AUPRC ─────────────────────────────────────────────

def plot_auprc(scores, outdir):
    """Box/strip plot of AUPRC by K, colored by read depth."""
    filtered = [s for s in scores if "cocluster_auprc" in s]
    if not filtered:
        print("  [skip] AUPRC — no co-clustering data")
        return

    Ks = sorted(set(s["K"] for s in filtered))
    Ts = sorted(set(s["T"] for s in filtered))
    colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(Ts)))
    T_color = {t: colors[i] for i, t in enumerate(Ts)}

    fig, ax = plt.subplots(figsize=(max(6, len(Ks) * 2), 5))

    width = 0.7 / len(Ts)
    for ti, T in enumerate(Ts):
        positions = []
        data = []
        for ki, K in enumerate(Ks):
            vals = [s["cocluster_auprc"] for s in filtered if s["K"] == K and s["T"] == T]
            if vals:
                positions.append(ki + (ti - len(Ts) / 2 + 0.5) * width)
                data.append(vals)

        if data:
            bp = ax.boxplot(data, positions=positions, widths=width * 0.8,
                           patch_artist=True, showfliers=True,
                           medianprops={"color": "black", "linewidth": 1.5})
            for patch in bp["boxes"]:
                patch.set_facecolor(T_color[T])
                patch.set_alpha(0.7)

    ax.set_xticks(range(len(Ks)))
    ax.set_xticklabels([f"K={k}" for k in Ks])
    ax.set_ylabel("Co-clustering AUPRC")
    ax.set_title("Clustering Accuracy by Clone Count and Read Depth")
    ax.set_ylim(-0.05, 1.05)
    ax.axhline(y=1.0, color="gray", linestyle="--", alpha=0.3, label="_nolegend_")

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=T_color[t], alpha=0.7, label=f"T={t}") for t in Ts]
    ax.legend(handles=legend_elements, title="Read Depth", loc="lower left")

    fig.savefig(os.path.join(outdir, "auprc_by_k_and_depth.png"))
    plt.close(fig)
    print("  [saved] auprc_by_k_and_depth.png")


# ── Plot 3: Runtime scaling ──────────────────────────────────────────────────

def plot_runtime(scores, outdir):
    """Scatter plot of runtime vs. M (number of SSMs), colored by K."""
    filtered = [s for s in scores if s.get("total_time_s")]
    if not filtered:
        print("  [skip] Runtime — no timing data")
        return

    for s in filtered:
        try:
            s["_time"] = float(s["total_time_s"])
        except (ValueError, TypeError):
            s["_time"] = None
    filtered = [s for s in filtered if s["_time"] is not None]

    if not filtered:
        print("  [skip] Runtime — no parseable timing data")
        return

    Ks = sorted(set(s["K"] for s in filtered))
    colors = plt.cm.tab10(np.linspace(0, 0.8, max(len(Ks), 1)))
    K_color = {k: colors[i] for i, k in enumerate(Ks)}

    fig, ax = plt.subplots(figsize=(7, 5))

    for K in Ks:
        subset = [s for s in filtered if s["K"] == K]
        Ms = [s["M"] for s in subset]
        times = [s["_time"] for s in subset]
        ax.scatter(Ms, times, c=[K_color[K]], s=50, alpha=0.7, label=f"K={K}",
                  edgecolors="white", linewidths=0.5)

    ax.set_xlabel("Number of SSMs (M)")
    ax.set_ylabel("Wall Time (seconds)")
    ax.set_title("Runtime Scaling by Fixture Size")
    ax.legend(title="Clones (K)")

    if filtered:
        max_time = max(s["_time"] for s in filtered)
        if max_time > 100:
            ax.set_yscale("log")

    fig.savefig(os.path.join(outdir, "runtime_scaling.png"))
    plt.close(fig)
    print("  [saved] runtime_scaling.png")


# ── Plot 4: Chain convergence ────────────────────────────────────────────────

def plot_chain_convergence(scores, outdir):
    """Bar chart of inter-chain LLH std, grouped by K."""
    filtered = [s for s in scores if "llh_chain_std" in s]
    if not filtered:
        print("  [skip] Chain convergence — no multi-chain data")
        return

    # Sort by K then M
    filtered.sort(key=lambda s: (s["K"], s["M"], s["T"]))

    fig, ax = plt.subplots(figsize=(max(8, len(filtered) * 0.4), 5))

    labels = [s["fixture"].replace("_", "\n") for s in filtered]
    stds = [s["llh_chain_std"] for s in filtered]
    Ks = [s["K"] for s in filtered]

    unique_Ks = sorted(set(Ks))
    colors = plt.cm.Set2(np.linspace(0, 0.8, max(len(unique_Ks), 1)))
    K_color = {k: colors[i] for i, k in enumerate(unique_Ks)}
    bar_colors = [K_color[k] for k in Ks]

    bars = ax.bar(range(len(stds)), stds, color=bar_colors, edgecolor="white", linewidth=0.5)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=90, fontsize=7, ha="center")
    ax.set_ylabel("Std Dev of Final LLH Across Chains")
    ax.set_title("Chain Convergence (lower = better agreement)")

    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=K_color[k], label=f"K={k}") for k in unique_Ks]
    ax.legend(handles=legend_elements, title="Clones")

    fig.savefig(os.path.join(outdir, "chain_convergence.png"))
    plt.close(fig)
    print("  [saved] chain_convergence.png")


# ── Plot 5: MCMC trace plots ────────────────────────────────────────────────

def plot_traces(result_base, fixture_names, outdir, max_fixtures=6):
    """Plot MCMC LLH traces for a sample of fixtures."""
    result_path = Path(result_base)
    plotted = 0

    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    axes = axes.flatten()

    for name in fixture_names:
        if plotted >= max_fixtures:
            break

        rdir = result_path / name
        if not rdir.exists():
            continue

        traces = load_chain_traces(rdir)
        if not traces:
            continue

        ax = axes[plotted]
        for chain_id, (iters, llhs) in traces.items():
            ax.plot(iters, llhs, alpha=0.7, linewidth=0.8, label=f"Chain {chain_id}")

        ax.set_title(name.replace("_", " "), fontsize=9)
        ax.set_xlabel("Iteration")
        ax.set_ylabel("Log-likelihood")
        if plotted == 0:
            ax.legend(fontsize=7)
        plotted += 1

    # Hide unused axes
    for i in range(plotted, len(axes)):
        axes[i].set_visible(False)

    fig.suptitle("MCMC Log-Likelihood Traces (sample of fixtures)", fontsize=14)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "mcmc_traces.png"))
    plt.close(fig)
    print(f"  [saved] mcmc_traces.png ({plotted} fixtures)")


# ── Plot 6: Summary dashboard ────────────────────────────────────────────────

def plot_summary_dashboard(scores, failures, outdir):
    """Single-page overview with key metrics.

    Adapts panels based on what data is available — always shows useful
    content from chain files (K error, LLH, runtime, convergence) even
    when the Go port doesn't export assignments (no AUPRC).
    """
    fig = plt.figure(figsize=(14, 9))
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.4, wspace=0.35)

    n_scored = len(scores)
    n_failed = len(failures)

    k_errors = [s["K_error"] for s in scores if "K_error" in s]
    auprcs = [s["cocluster_auprc"] for s in scores if "cocluster_auprc" in s]
    times = [s["total_time_s"] for s in scores
             if isinstance(s.get("total_time_s"), (int, float))]
    best_llhs = [s["best_llh"] for s in scores if s.get("best_llh") is not None]
    chain_stds = [s["llh_chain_std"] for s in scores if "llh_chain_std" in s]

    # Panel A: pass/fail pie
    ax_pie = fig.add_subplot(gs[0, 0])
    if n_scored + n_failed > 0:
        sizes = [n_scored, n_failed]
        labels = [f"Scored ({n_scored})", f"Failed ({n_failed})"]
        colors_pie = ["#66bb6a", "#ef5350"] if n_failed > 0 else ["#66bb6a", "#e0e0e0"]
        ax_pie.pie(sizes, labels=labels, colors=colors_pie,
                   autopct="%1.0f%%", startangle=90)
    ax_pie.set_title(f"Fixture Results ({n_scored + n_failed} total)")

    # Panel B: K error distribution (from chain NumNodes)
    ax_kerr = fig.add_subplot(gs[0, 1])
    if k_errors:
        unique_errs = sorted(set(k_errors))
        counts = [k_errors.count(e) for e in unique_errs]
        bar_colors = ["#66bb6a" if e == 0 else "#ffa726" if e <= 1 else "#ef5350"
                      for e in unique_errs]
        ax_kerr.bar(unique_errs, counts, color=bar_colors, edgecolor="white")
        ax_kerr.set_xlabel("|K error|")
        ax_kerr.set_ylabel("Count")
        exact = sum(1 for e in k_errors if e == 0)
        ax_kerr.set_title(f"Pop Count Error ({exact}/{len(k_errors)} exact)")
    else:
        ax_kerr.text(0.5, 0.5, "No K data\n(chain files needed)",
                    ha="center", va="center", transform=ax_kerr.transAxes, color="gray")
        ax_kerr.set_title("Population Count Error")

    # Panel C: AUPRC or Best LLH distribution (adaptive)
    ax_c = fig.add_subplot(gs[0, 2])
    if auprcs:
        ax_c.hist(auprcs, bins=20, color="#42a5f5", edgecolor="white", alpha=0.8)
        ax_c.axvline(np.median(auprcs), color="red", linestyle="--",
                    label=f"median={np.median(auprcs):.3f}")
        ax_c.legend(fontsize=8)
        ax_c.set_xlabel("AUPRC")
        ax_c.set_ylabel("Count")
        ax_c.set_title("Co-clustering AUPRC")
    elif best_llhs:
        ax_c.hist(best_llhs, bins=20, color="#ab47bc", edgecolor="white", alpha=0.8)
        ax_c.axvline(np.median(best_llhs), color="red", linestyle="--",
                    label=f"median={np.median(best_llhs):.0f}")
        ax_c.legend(fontsize=8)
        ax_c.set_xlabel("Best Log-Likelihood")
        ax_c.set_ylabel("Count")
        ax_c.set_title("Best LLH Distribution")
    else:
        ax_c.text(0.5, 0.5, "No LLH data", ha="center", va="center",
                 transform=ax_c.transAxes, color="gray")
        ax_c.set_title("Log-Likelihood")

    # Panel D: K error by true K (or runtime by M)
    ax_d = fig.add_subplot(gs[1, 0])
    if k_errors:
        k_with_errors = [s for s in scores if "K_error" in s]
        Ks_plot = [s["K"] for s in k_with_errors]
        errs = [s["K_error"] for s in k_with_errors]
        ax_d.scatter(Ks_plot, errs, c="#ef5350", s=40, alpha=0.6, edgecolors="white")
        for K in sorted(set(Ks_plot)):
            vals = [e for k, e in zip(Ks_plot, errs) if k == K]
            ax_d.plot(K, np.mean(vals), "kv", markersize=10)
        ax_d.set_xlabel("True Clones (K)")
        ax_d.set_ylabel("|K error|")
        ax_d.set_title("K Error vs Clone Count")
    elif times:
        Ms = [s["M"] for s in scores if isinstance(s.get("total_time_s"), (int, float))]
        ax_d.scatter(Ms, times, c="#42a5f5", s=40, alpha=0.6, edgecolors="white")
        ax_d.set_xlabel("Number of SSMs (M)")
        ax_d.set_ylabel("Wall Time (s)")
        ax_d.set_title("Runtime vs Fixture Size")
    else:
        ax_d.text(0.5, 0.5, "No data", ha="center", va="center",
                 transform=ax_d.transAxes, color="gray")
        ax_d.set_title("K Error / Runtime")

    # Panel E: Chain convergence (LLH std)
    ax_e = fig.add_subplot(gs[1, 1])
    if chain_stds:
        ax_e.hist(chain_stds, bins=20, color="#ffa726", edgecolor="white", alpha=0.8)
        ax_e.axvline(np.median(chain_stds), color="red", linestyle="--",
                    label=f"median={np.median(chain_stds):.1f}")
        ax_e.legend(fontsize=8)
        ax_e.set_xlabel("Std Dev of Final LLH")
        ax_e.set_ylabel("Count")
        ax_e.set_title("Chain Convergence (lower=better)")
    else:
        ax_e.text(0.5, 0.5, "No multi-chain data", ha="center", va="center",
                 transform=ax_e.transAxes, color="gray")
        ax_e.set_title("Chain Convergence")

    # Panel F: text summary
    ax_text = fig.add_subplot(gs[1, 2])
    ax_text.axis("off")
    lines = ["Summary Statistics", "=" * 25]
    lines.append(f"Fixtures scored: {n_scored}")
    lines.append(f"Fixtures failed: {n_failed}")
    if k_errors:
        exact = sum(1 for e in k_errors if e == 0)
        lines.append(f"K exact match: {exact}/{len(k_errors)}")
        lines.append(f"Mean |K error|: {np.mean(k_errors):.2f}")
    if auprcs:
        lines.append(f"AUPRC mean: {np.mean(auprcs):.3f}")
        lines.append(f"AUPRC median: {np.median(auprcs):.3f}")
    else:
        lines.append("AUPRC: N/A (no assignments)")
    if best_llhs:
        lines.append(f"Best LLH median: {np.median(best_llhs):.0f}")
    if times:
        lines.append(f"Mean runtime: {np.mean(times):.1f}s")
        lines.append(f"Max runtime: {np.max(times):.1f}s")
    if chain_stds:
        lines.append(f"Chain std median: {np.median(chain_stds):.1f}")

    # CNV breakdown if present
    cnv_count = sum(1 for s in scores if s.get("has_cnvs"))
    ssm_count = sum(1 for s in scores if not s.get("has_cnvs"))
    if cnv_count > 0:
        lines.append(f"SSM-only: {ssm_count} | With CNV: {cnv_count}")

    ax_text.text(0.05, 0.95, "\n".join(lines), transform=ax_text.transAxes,
                fontsize=10, verticalalignment="top", fontfamily="monospace",
                bbox=dict(boxstyle="round,pad=0.5", facecolor="#f5f5f5", alpha=0.8))

    fig.suptitle("PhyloWGS Go Port — Simulation Validation Dashboard", fontsize=15, y=1.02)
    fig.savefig(os.path.join(outdir, "dashboard.png"))
    plt.close(fig)
    print("  [saved] dashboard.png")


# ── Plot 7: CNV vs SSM-only comparison ────────────────────────────────────────

def plot_cnv_comparison(scores, outdir):
    """Side-by-side box plots comparing AUPRC for SSM-only vs. SSM+CNV fixtures."""
    filtered = [s for s in scores if "cocluster_auprc" in s and "has_cnvs" in s]
    if not filtered:
        print("  [skip] CNV comparison — no data")
        return

    ssm_only = [s for s in filtered if not s["has_cnvs"]]
    with_cnv = [s for s in filtered if s["has_cnvs"]]

    if not ssm_only or not with_cnv:
        print("  [skip] CNV comparison — need both SSM-only and CNV fixtures")
        return

    Ks = sorted(set(s["K"] for s in filtered))

    fig, ax = plt.subplots(figsize=(max(6, len(Ks) * 2), 5))
    width = 0.3

    for ki, K in enumerate(Ks):
        ssm_vals = [s["cocluster_auprc"] for s in ssm_only if s["K"] == K]
        cnv_vals = [s["cocluster_auprc"] for s in with_cnv if s["K"] == K]

        if ssm_vals:
            bp1 = ax.boxplot([ssm_vals], positions=[ki - width/2], widths=width * 0.8,
                            patch_artist=True, showfliers=True,
                            medianprops={"color": "black", "linewidth": 1.5})
            bp1["boxes"][0].set_facecolor("#66bb6a")
            bp1["boxes"][0].set_alpha(0.7)

        if cnv_vals:
            bp2 = ax.boxplot([cnv_vals], positions=[ki + width/2], widths=width * 0.8,
                            patch_artist=True, showfliers=True,
                            medianprops={"color": "black", "linewidth": 1.5})
            bp2["boxes"][0].set_facecolor("#ef5350")
            bp2["boxes"][0].set_alpha(0.7)

    ax.set_xticks(range(len(Ks)))
    ax.set_xticklabels([f"K={k}" for k in Ks])
    ax.set_ylabel("Co-clustering AUPRC")
    ax.set_title("Clustering Accuracy: SSM-only vs. SSM+CNV")
    ax.set_ylim(-0.05, 1.05)

    from matplotlib.patches import Patch
    ax.legend(handles=[
        Patch(facecolor="#66bb6a", alpha=0.7, label="SSM-only"),
        Patch(facecolor="#ef5350", alpha=0.7, label="With CNVs"),
    ], loc="lower left")

    fig.savefig(os.path.join(outdir, "cnv_comparison.png"))
    plt.close(fig)
    print("  [saved] cnv_comparison.png")


def plot_cnv_impact_heatmap(scores, outdir):
    """Heatmap of mean AUPRC with axes (K, N_cnv)."""
    filtered = [s for s in scores if "cocluster_auprc" in s and "N_cnv" in s]
    if not filtered:
        print("  [skip] CNV impact heatmap — no data")
        return

    Ks = sorted(set(s["K"] for s in filtered))
    Ncnvs = sorted(set(s["N_cnv"] for s in filtered))

    if len(Ncnvs) < 2:
        print("  [skip] CNV impact heatmap — need multiple N_cnv values")
        return

    groups = {}
    for s in filtered:
        key = (s["K"], s["N_cnv"])
        groups.setdefault(key, []).append(s["cocluster_auprc"])

    matrix = np.full((len(Ks), len(Ncnvs)), np.nan)
    for i, K in enumerate(Ks):
        for j, nc in enumerate(Ncnvs):
            vals = groups.get((K, nc), [])
            if vals:
                matrix[i, j] = np.mean(vals)

    fig, ax = plt.subplots(figsize=(max(5, len(Ncnvs) * 1.5), max(4, len(Ks) * 1.2)))
    im = ax.imshow(matrix, cmap="RdYlGn", aspect="auto", vmin=0, vmax=1)

    ax.set_xticks(range(len(Ncnvs)))
    ax.set_xticklabels([str(nc) for nc in Ncnvs])
    ax.set_yticks(range(len(Ks)))
    ax.set_yticklabels([str(k) for k in Ks])
    ax.set_xlabel("Number of CNVs")
    ax.set_ylabel("True Clones (K)")
    ax.set_title("Mean AUPRC by Clone Count and CNV Count")

    for i in range(len(Ks)):
        for j in range(len(Ncnvs)):
            val = matrix[i, j]
            if not np.isnan(val):
                n = len(groups.get((Ks[i], Ncnvs[j]), []))
                color = "white" if val < 0.4 else "black"
                ax.text(j, i, f"{val:.2f}\n(n={n})", ha="center", va="center",
                        fontsize=9, color=color)

    plt.colorbar(im, ax=ax, label="Mean AUPRC")
    fig.savefig(os.path.join(outdir, "cnv_impact_heatmap.png"))
    plt.close(fig)
    print("  [saved] cnv_impact_heatmap.png")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Plot simulation validation results")
    parser.add_argument("--analysis-dir", required=True, help="Directory with scores.json")
    parser.add_argument("--result-base", help="Base dir with per-fixture results (for trace plots)")
    parser.add_argument("--fixture-base", help="Base dir with fixtures (for manifest)")
    parser.add_argument("--outdir", help="Output dir for plots (default: analysis-dir/plots)")

    args = parser.parse_args()
    outdir = args.outdir or os.path.join(args.analysis_dir, "plots")
    os.makedirs(outdir, exist_ok=True)

    print("Loading scores...")
    scores, failures = load_scores(args.analysis_dir)
    print(f"  {len(scores)} scored, {len(failures)} failed")

    if not scores:
        print("No scores to plot.")
        return

    print("\nGenerating plots:")
    plot_summary_dashboard(scores, failures, outdir)
    plot_k_accuracy(scores, outdir)
    plot_auprc(scores, outdir)
    plot_runtime(scores, outdir)
    plot_chain_convergence(scores, outdir)

    # CNV-specific plots (only if mixed SSM-only and CNV fixtures exist)
    has_cnv_data = any(s.get("has_cnvs") for s in scores)
    has_ssm_only = any(not s.get("has_cnvs") for s in scores)
    if has_cnv_data and has_ssm_only:
        plot_cnv_comparison(scores, outdir)
        plot_cnv_impact_heatmap(scores, outdir)

    # Trace plots if result-base provided
    if args.result_base:
        # Pick diverse fixtures: one per K value
        Ks_seen = set()
        sample_names = []
        for s in scores:
            if s["K"] not in Ks_seen:
                sample_names.append(s["fixture"])
                Ks_seen.add(s["K"])
        plot_traces(args.result_base, sample_names, outdir)

    print(f"\nAll plots saved to: {outdir}")


if __name__ == "__main__":
    main()
