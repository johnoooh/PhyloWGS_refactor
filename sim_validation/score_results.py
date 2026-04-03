#!/usr/bin/env python3
"""
score_results.py — Score PhyloWGS Go port results against simulation ground truth.

Metrics (extracted from actual Go port output):
  1. Population count error: |inferred_K - true_K| (from chain NumNodes column)
  2. Best log-likelihood and median LLH (from chain traces)
  3. Chain convergence: inter-chain LLH std
  4. Runtime (from summary.json total_time or timing.json)

Note: Co-clustering AUPRC requires the Go port to export SSM assignments,
which it does not currently do. That metric is computed when available.

Usage:
    python score_results.py --fixture-dir fixtures/K5_S3_T1000_M50_rep0 \
                            --result-dir results/K5_S3_T1000_M50_rep0
    python score_results.py --all --fixture-base fixtures/ --result-base results/ \
                            --outdir analysis/
"""

import argparse
import json
import os
import re
import sys
from pathlib import Path

import numpy as np


# ── Co-clustering matrix ─────────────────────────────────────────────────────

def build_cocluster_matrix(assignments, M):
    """Build binary co-clustering matrix from SSM assignments."""
    C = np.zeros((M, M), dtype=float)
    for i in range(M):
        for j in range(M):
            C[i, j] = 1.0 if assignments[i] == assignments[j] else 0.0
    return C


def compute_auprc(true_matrix, pred_matrix):
    """Compute area under the precision-recall curve for co-clustering."""
    M = true_matrix.shape[0]
    idx = np.triu_indices(M, k=1)
    y_true = true_matrix[idx]
    y_score = pred_matrix[idx]

    if y_true.sum() == 0 or y_true.sum() == len(y_true):
        return 1.0

    order = np.argsort(-y_score)
    y_true_sorted = y_true[order]

    tp = np.cumsum(y_true_sorted)
    fp = np.cumsum(1 - y_true_sorted)
    precision = tp / (tp + fp)
    recall = tp / y_true.sum()

    trapz = getattr(np, 'trapezoid', getattr(np, 'trapz', None))
    auprc = trapz(precision, recall)
    return float(auprc)


# ── Parse Go port output ─────────────────────────────────────────────────────

def load_truth(fixture_dir):
    """Load ground truth from fixture directory."""
    path = Path(fixture_dir) / "truth.json"
    with open(path) as f:
        return json.load(f)


def load_go_summary(result_dir):
    """Load Go port summary.json."""
    for name in ("summary.json", "results.json"):
        path = Path(result_dir) / name
        if path.exists():
            with open(path) as f:
                return json.load(f)
    return None


def parse_go_duration(s):
    """Parse Go duration string (e.g. '38.69s', '725.97ms', '1m30.5s') to seconds."""
    if not s:
        return None
    s = s.strip()
    # Try simple float ending in 's' (not 'ms')
    m = re.match(r'^([\d.]+)s$', s)
    if m:
        return float(m.group(1))
    # Milliseconds
    m = re.match(r'^([\d.]+)ms$', s)
    if m:
        return float(m.group(1)) / 1000.0
    # Minutes + seconds: "1m30.5s"
    m = re.match(r'^(?:(\d+)m)?([\d.]+)s$', s)
    if m:
        mins = int(m.group(1)) if m.group(1) else 0
        secs = float(m.group(2))
        return mins * 60 + secs
    # Hours
    m = re.match(r'^(?:(\d+)h)?(?:(\d+)m)?([\d.]+)s$', s)
    if m:
        hrs = int(m.group(1)) if m.group(1) else 0
        mins = int(m.group(2)) if m.group(2) else 0
        secs = float(m.group(3))
        return hrs * 3600 + mins * 60 + secs
    # Microseconds
    m = re.match(r'^([\d.]+)[µu]s$', s)
    if m:
        return float(m.group(1)) / 1e6
    return None


def load_chain_traces(result_dir):
    """Load MCMC traces from chain files.

    Supports both Go format (chain_*_samples.txt) and Python format (mcmc_samples.txt).
    Returns dict: chain_id -> {"llhs": [...], "num_nodes": [...]}
    """
    traces = {}
    result_path = Path(result_dir)

    # Go format: chain_N_samples.txt with Iteration/LLH/NumNodes
    for chain_file in sorted(result_path.glob("chain_*_samples.txt")):
        chain_id = chain_file.stem.split("_")[1]
        llhs = []
        num_nodes = []
        with open(chain_file) as f:
            header = f.readline().strip().split("\t")
            llh_col = 1
            nodes_col = 2 if len(header) > 2 else None

            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    llhs.append(float(parts[llh_col]))
                if nodes_col is not None and len(parts) > nodes_col:
                    try:
                        num_nodes.append(int(parts[nodes_col]))
                    except ValueError:
                        pass

        traces[chain_id] = {"llhs": llhs, "num_nodes": num_nodes}

    # Python format: mcmc_samples.txt with iteration/llh/time (single chain)
    if not traces:
        mcmc_path = result_path / "mcmc_samples.txt"
        if mcmc_path.exists():
            llhs = []
            with open(mcmc_path) as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) >= 2:
                        try:
                            llhs.append(float(parts[1]))
                        except ValueError:
                            pass
            if llhs:
                traces["0"] = {"llhs": llhs, "num_nodes": []}

    return traces


def load_best_tree(result_dir):
    """Load best tree output — supports Go (best_tree.json) and Python (tree_summaries.json.gz).

    Returns a normalized dict with: num_populations, populations, mut_assignments, llh
    """
    result_path = Path(result_dir)

    # Go port format
    go_path = result_path / "best_tree.json"
    if go_path.exists():
        with open(go_path) as f:
            return json.load(f)

    # Original Python format (from write_results.py)
    py_path = result_path / "tree_summaries.json.gz"
    if py_path.exists():
        import gzip
        with gzip.open(py_path, "rt") as f:
            summaries = json.load(f)
        # summaries["trees"] is a dict of tree_idx -> {llh, structure, populations}
        trees = summaries.get("trees", {})
        if not trees:
            return None
        # Pick tree with best LLH
        best_idx = max(trees.keys(), key=lambda k: trees[k].get("llh", float("-inf")))
        best = trees[best_idx]
        # Load mutation assignments from mutass.zip
        mutass_path = result_path / "mutass.zip"
        mut_assignments = {}
        if mutass_path.exists():
            import zipfile
            with zipfile.ZipFile(mutass_path) as zf:
                try:
                    with zf.open(f"{best_idx}.json") as mf:
                        mutdata = json.loads(mf.read())
                        mut_assignments = mutdata.get("mut_assignments", {})
                except KeyError:
                    pass

        num_pops = len([p for p in best.get("populations", {}).values()
                       if p.get("num_ssms", 0) > 0 or p.get("num_cnvs", 0) > 0])
        return {
            "num_populations": num_pops,
            "llh": best.get("llh"),
            "populations": best.get("populations", {}),
            "structure": best.get("structure", {}),
            "mut_assignments": mut_assignments,
        }

    return None


def load_timing(result_dir):
    """Load timing.json if present (written by SLURM job wrapper)."""
    path = Path(result_dir) / "timing.json"
    if path.exists():
        with open(path) as f:
            return json.load(f)
    return None


# ── Scoring ──────────────────────────────────────────────────────────────────

def score_fixture(fixture_dir, result_dir):
    """Score a single fixture's results against ground truth."""
    truth = load_truth(fixture_dir)
    summary = load_go_summary(result_dir)
    traces = load_chain_traces(result_dir)
    timing = load_timing(result_dir)

    params = truth["params"]
    true_K = params["K"]
    M = params["M"]

    scores = {
        "fixture": os.path.basename(fixture_dir),
        "K": true_K,
        "S": params["S"],
        "T": params["T"],
        "M": M,
        "has_cnvs": params.get("has_cnvs", False),
        "N_cnv": params.get("N_cnv", 0),
    }

    # ── Runtime ──────────────────────────────────────────────────────────
    # Try timing.json first (SLURM wrapper), then summary.json
    if timing and "wall_seconds" in timing:
        scores["total_time_s"] = float(timing["wall_seconds"])
    elif summary and "total_time" in summary:
        parsed = parse_go_duration(summary["total_time"])
        if parsed is not None:
            scores["total_time_s"] = parsed

    # ── Best log-likelihood ──────────────────────────────────────────────
    if summary:
        scores["best_llh"] = summary.get("best_llh", None)

    # ── Load best_tree.json (from Go port post-MCMC output) ────────────
    best_tree = load_best_tree(result_dir)

    # ── Population count ─────────────────────────────────────────────────
    if best_tree and "num_populations" in best_tree:
        inferred_K = best_tree["num_populations"]
        scores["inferred_K"] = inferred_K
        scores["K_error"] = abs(inferred_K - true_K)
    elif traces:
        # Fallback: use NumNodes from chain files (less accurate, includes
        # structural TSSB nodes — see INVESTIGATION_NumNodes.md)
        best_llh_val = float("-inf")
        best_nodes = None
        for trace in traces.values():
            llhs = trace["llhs"]
            nodes = trace["num_nodes"]
            if llhs and nodes and len(llhs) == len(nodes):
                max_idx = int(np.argmax(llhs))
                if llhs[max_idx] > best_llh_val:
                    best_llh_val = llhs[max_idx]
                    best_nodes = nodes[max_idx]
        if best_nodes is not None:
            scores["inferred_K_tssb"] = best_nodes
            scores["K_error_tssb"] = abs(best_nodes - 1 - true_K)

    # ── Median LLH (from post-burnin samples) ────────────────────────────
    if traces:
        all_llhs = []
        for trace in traces.values():
            llhs = trace["llhs"]
            if llhs:
                # Use last half as post-burnin
                n = len(llhs)
                post_burnin = llhs[n // 2:]
                all_llhs.extend(post_burnin)
        if all_llhs:
            scores["median_llh"] = float(np.median(all_llhs))

    # ── Chain convergence ────────────────────────────────────────────────
    if traces:
        final_llhs = []
        for trace in traces.values():
            llhs = trace["llhs"]
            if llhs:
                final_llhs.append(llhs[-1])
        if len(final_llhs) > 1:
            scores["llh_chain_std"] = float(np.std(final_llhs))
            scores["llh_chain_mean"] = float(np.mean(final_llhs))

    # ── Co-clustering AUPRC (from best_tree.json mut_assignments) ────────
    if best_tree and "mut_assignments" in best_tree:
        true_assignments = np.array(truth["assignments"])
        # Build inferred assignment: ssm_id -> population_id
        inferred_map = {}
        for pop_id, muts in best_tree["mut_assignments"].items():
            for ssm_id in muts.get("ssms", []):
                inferred_map[ssm_id] = int(pop_id)

        if len(inferred_map) == M:
            # Convert to array matching SSM order (s0, s1, s2, ...)
            inferred = np.array([inferred_map[f"s{m}"] for m in range(M)])
            true_C = build_cocluster_matrix(true_assignments, M)
            pred_C = build_cocluster_matrix(inferred, M)
            scores["cocluster_auprc"] = compute_auprc(true_C, pred_C)

    return scores


def score_all(fixture_base, result_base, outdir):
    """Score all fixtures and produce summary report."""
    fixture_base = Path(fixture_base)
    result_base = Path(result_base)
    os.makedirs(outdir, exist_ok=True)

    manifest_path = fixture_base / "manifest.json"
    if manifest_path.exists():
        with open(manifest_path) as f:
            manifest = json.load(f)
        fixture_names = [entry["name"] for entry in manifest]
    else:
        fixture_names = sorted(
            d.name for d in fixture_base.iterdir()
            if d.is_dir() and (d / "truth.json").exists()
        )

    all_scores = []
    failures = []

    for name in fixture_names:
        fixture_dir = fixture_base / name
        result_dir = result_base / name

        if not result_dir.exists():
            failures.append({"fixture": name, "error": "no result directory"})
            continue

        try:
            scores = score_fixture(str(fixture_dir), str(result_dir))
            all_scores.append(scores)
        except Exception as e:
            failures.append({"fixture": name, "error": str(e)})

    # Write detailed results
    results_path = os.path.join(outdir, "scores.json")
    with open(results_path, "w") as f:
        json.dump({"scores": all_scores, "failures": failures}, f, indent=2)

    # Write TSV summary
    tsv_path = os.path.join(outdir, "scores.tsv")
    if all_scores:
        # Collect all keys across all scores
        all_keys = []
        seen = set()
        for row in all_scores:
            for k in row:
                if k not in seen:
                    all_keys.append(k)
                    seen.add(k)
        with open(tsv_path, "w") as f:
            f.write("\t".join(all_keys) + "\n")
            for row in all_scores:
                f.write("\t".join(str(row.get(c, "")) for c in all_keys) + "\n")

    # Print summary
    print(f"\n{'='*60}")
    print(f"Simulation Validation Summary")
    print(f"{'='*60}")
    print(f"Scored:  {len(all_scores)}/{len(fixture_names)} fixtures")
    print(f"Failed:  {len(failures)}")

    if all_scores:
        k_errors = [s["K_error"] for s in all_scores if "K_error" in s]
        auprcs = [s["cocluster_auprc"] for s in all_scores if "cocluster_auprc" in s]
        times = [s["total_time_s"] for s in all_scores if "total_time_s" in s]
        best_llhs = [s["best_llh"] for s in all_scores if s.get("best_llh") is not None]

        if k_errors:
            print(f"\nPopulation count error (from chain NumNodes):")
            print(f"  Mean: {np.mean(k_errors):.2f}")
            print(f"  Exact: {sum(1 for e in k_errors if e == 0)}/{len(k_errors)}")

        if auprcs:
            print(f"\nCo-clustering AUPRC:")
            print(f"  Mean:   {np.mean(auprcs):.3f}")
            print(f"  Median: {np.median(auprcs):.3f}")
            print(f"  Min:    {np.min(auprcs):.3f}")

        if times:
            print(f"\nRuntime:")
            print(f"  Mean:  {np.mean(times):.1f}s")
            print(f"  Max:   {np.max(times):.1f}s")

        if best_llhs:
            print(f"\nBest LLH:")
            print(f"  Mean:   {np.mean(best_llhs):.1f}")
            print(f"  Median: {np.median(best_llhs):.1f}")

    if failures:
        print(f"\nFailures:")
        for fail in failures[:10]:
            print(f"  {fail['fixture']}: {fail['error']}")
        if len(failures) > 10:
            print(f"  ... and {len(failures) - 10} more")

    print(f"\nDetailed results: {results_path}")
    print(f"TSV summary:      {tsv_path}")

    return all_scores, failures


# ── CLI ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Score PhyloWGS results against simulation ground truth"
    )
    parser.add_argument("--all", action="store_true", help="Score all fixtures")
    parser.add_argument("--fixture-dir", help="Single fixture directory")
    parser.add_argument("--result-dir", help="Single result directory")
    parser.add_argument("--fixture-base", help="Base dir for all fixtures (with --all)")
    parser.add_argument("--result-base", help="Base dir for all results (with --all)")
    parser.add_argument("--outdir", default="analysis", help="Output directory for reports")

    args = parser.parse_args()

    if args.all:
        if not args.fixture_base or not args.result_base:
            parser.error("--all requires --fixture-base and --result-base")
        score_all(args.fixture_base, args.result_base, args.outdir)
    elif args.fixture_dir and args.result_dir:
        scores = score_fixture(args.fixture_dir, args.result_dir)
        print(json.dumps(scores, indent=2))
    else:
        parser.error("Provide --all or both --fixture-dir and --result-dir")


if __name__ == "__main__":
    main()
