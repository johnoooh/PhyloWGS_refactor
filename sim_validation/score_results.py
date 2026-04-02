#!/usr/bin/env python3
"""
score_results.py — Score PhyloWGS Go port results against simulation ground truth.

Metrics:
  1. Population count error: |inferred_K - true_K|
  2. Co-clustering AUPRC: do SSMs from the same clone cluster together?
  3. Phi MSE: mean squared error of inferred vs. true cellular prevalences
  4. Best log-likelihood: convergence quality

Reads:
  - truth.json                (from generate_fixtures.py)
  - summary.json              (from PhyloWGS Go port)
  - chain_*_samples.txt       (MCMC trace for convergence)

Usage:
    python score_results.py --fixture-dir fixtures/K5_S3_T1000_M50_rep0 \
                            --result-dir results/K5_S3_T1000_M50_rep0
    python score_results.py --all --fixture-base fixtures/ --result-base results/ \
                            --outdir analysis/
"""

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np


# ── Co-clustering matrix ─────────────────────────────────────────────────────

def build_cocluster_matrix(assignments, M):
    """Build binary co-clustering matrix from SSM assignments.

    C[i,j] = 1 if SSM i and SSM j are in the same clone.
    """
    C = np.zeros((M, M), dtype=float)
    for i in range(M):
        for j in range(M):
            C[i, j] = 1.0 if assignments[i] == assignments[j] else 0.0
    return C


def compute_auprc(true_matrix, pred_matrix):
    """Compute area under the precision-recall curve for co-clustering.

    Uses the upper triangle only (symmetric matrix, diagonal is trivial).
    """
    M = true_matrix.shape[0]
    idx = np.triu_indices(M, k=1)
    y_true = true_matrix[idx]
    y_score = pred_matrix[idx]

    if y_true.sum() == 0 or y_true.sum() == len(y_true):
        return 1.0  # trivial case

    # sort by score descending
    order = np.argsort(-y_score)
    y_true_sorted = y_true[order]

    tp = np.cumsum(y_true_sorted)
    fp = np.cumsum(1 - y_true_sorted)
    precision = tp / (tp + fp)
    recall = tp / y_true.sum()

    # AUPRC via trapezoidal rule
    auprc = np.trapz(precision, recall)
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


def load_chain_traces(result_dir):
    """Load MCMC log-likelihood traces from chain output files."""
    traces = {}
    result_path = Path(result_dir)
    for chain_file in sorted(result_path.glob("chain_*_samples.txt")):
        chain_id = chain_file.stem.split("_")[1]
        llhs = []
        with open(chain_file) as f:
            header = f.readline()
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    llhs.append(float(parts[1]))
        traces[chain_id] = llhs
    return traces


def infer_assignments_from_summary(summary, M):
    """Extract SSM-to-cluster assignments from Go port summary.

    The Go port outputs tree structure in summary.json. We parse the
    best tree's node assignments.
    """
    # The Go port stores per-chain results. Extract best tree.
    if "trees" in summary:
        # Find tree with best LLH
        best_tree = max(summary["trees"], key=lambda t: t.get("llh", float("-inf")))
        if "assignments" in best_tree:
            return np.array(best_tree["assignments"])

    # Fallback: try to extract from chain files
    return None


def count_populations_from_summary(summary):
    """Count the number of non-empty populations inferred."""
    if "trees" in summary:
        best_tree = max(summary["trees"], key=lambda t: t.get("llh", float("-inf")))
        return best_tree.get("num_nodes", None)
    if "num_chains" in summary:
        return summary.get("avg_nodes", None)
    return None


# ── Scoring ──────────────────────────────────────────────────────────────────

def score_fixture(fixture_dir, result_dir):
    """Score a single fixture's results against ground truth."""
    truth = load_truth(fixture_dir)
    summary = load_go_summary(result_dir)
    traces = load_chain_traces(result_dir)

    params = truth["params"]
    true_K = params["K"]
    M = params["M"]
    true_assignments = np.array(truth["assignments"])

    scores = {
        "fixture": os.path.basename(fixture_dir),
        "K": true_K,
        "S": params["S"],
        "T": params["T"],
        "M": M,
    }

    # Best log-likelihood
    if summary:
        scores["best_llh"] = summary.get("best_llh", None)
        scores["total_time_s"] = summary.get("total_time", "").rstrip("s")

    # Population count
    if summary:
        inferred_K = count_populations_from_summary(summary)
        if inferred_K is not None:
            # subtract 1 for root node
            scores["inferred_K"] = inferred_K - 1 if inferred_K > 0 else inferred_K
            scores["K_error"] = abs(scores["inferred_K"] - true_K)

    # Co-clustering AUPRC (if assignments available)
    if summary:
        inferred_assignments = infer_assignments_from_summary(summary, M)
        if inferred_assignments is not None:
            true_C = build_cocluster_matrix(true_assignments, M)
            pred_C = build_cocluster_matrix(inferred_assignments, M)
            scores["cocluster_auprc"] = compute_auprc(true_C, pred_C)

    # Chain convergence: variance of final LLH across chains
    if traces:
        final_llhs = [trace[-1] for trace in traces.values() if trace]
        if len(final_llhs) > 1:
            scores["llh_chain_std"] = float(np.std(final_llhs))
            scores["llh_chain_mean"] = float(np.mean(final_llhs))

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
        cols = list(all_scores[0].keys())
        with open(tsv_path, "w") as f:
            f.write("\t".join(cols) + "\n")
            for row in all_scores:
                f.write("\t".join(str(row.get(c, "")) for c in cols) + "\n")

    # Print summary
    print(f"\n{'='*60}")
    print(f"Simulation Validation Summary")
    print(f"{'='*60}")
    print(f"Scored:  {len(all_scores)}/{len(fixture_names)} fixtures")
    print(f"Failed:  {len(failures)}")

    if all_scores:
        k_errors = [s["K_error"] for s in all_scores if "K_error" in s]
        auprcs = [s["cocluster_auprc"] for s in all_scores if "cocluster_auprc" in s]

        if k_errors:
            print(f"\nPopulation count error:")
            print(f"  Mean: {np.mean(k_errors):.2f}")
            print(f"  Exact: {sum(1 for e in k_errors if e == 0)}/{len(k_errors)}")

        if auprcs:
            print(f"\nCo-clustering AUPRC:")
            print(f"  Mean:   {np.mean(auprcs):.3f}")
            print(f"  Median: {np.median(auprcs):.3f}")
            print(f"  Min:    {np.min(auprcs):.3f}")

    if failures:
        print(f"\nFailures:")
        for fail in failures:
            print(f"  {fail['fixture']}: {fail['error']}")

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
