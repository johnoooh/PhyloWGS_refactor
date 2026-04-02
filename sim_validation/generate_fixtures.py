#!/usr/bin/env python3
"""
generate_fixtures.py — Generate simulated PhyloWGS inputs with known ground truth.

Implements the pearsim generative model:
  1. Random tree topology (parents vector)
  2. Dirichlet-sampled subclone proportions (eta)
  3. Cellular prevalences via ancestral matrix (phi = Z @ eta)
  4. Binomial-sampled read counts

Outputs per fixture:
  - ssm_data.txt       (PhyloWGS SSM input format)
  - cnv_data.txt       (empty CNV file — SSM-only simulations)
  - truth.json         (ground truth: tree, phi, clusters, params)

Usage:
    python generate_fixtures.py --outdir fixtures/ --grid default
    python generate_fixtures.py --outdir fixtures/ -K 5 -S 3 -T 1000 -M 50 --seed 42
"""

import argparse
import json
import os
import sys

import numpy as np


# ── Tree generation ──────────────────────────────────────────────────────────

def make_parents(K, tree_type="random", rng=None):
    """Generate a random tree with K clones (excluding root=0).

    Returns parents vector of length K where parents[i] = parent of node i+1.
    Uses Morris lab convention: mu=0.75 probability of extending current branch.
    """
    if rng is None:
        rng = np.random.default_rng()

    parents = np.zeros(K, dtype=int)
    if K == 0:
        return parents

    parents[0] = 0  # node 1's parent is root

    current_branch_tip = 1
    for i in range(1, K):
        node_id = i + 1
        if rng.random() < 0.75:
            # extend current branch
            parents[i] = current_branch_tip
        else:
            # attach to random existing node (including root)
            parents[i] = rng.integers(0, node_id)
        current_branch_tip = node_id

    # enforce tree type constraints
    if tree_type == "monoprimary":
        # root has exactly one child — force all root children except first onto node 1
        root_children = [i for i in range(K) if parents[i] == 0]
        if len(root_children) > 1:
            for idx in root_children[1:]:
                parents[idx] = 1
    elif tree_type == "polyprimary":
        # root has multiple children — ensure at least 2
        root_children = [i for i in range(K) if parents[i] == 0]
        if len(root_children) < 2 and K >= 2:
            # pick a non-root-child and reattach to root
            non_root = [i for i in range(K) if parents[i] != 0]
            if non_root:
                parents[rng.choice(non_root)] = 0

    return parents


def parents_to_ancestral_matrix(parents):
    """Convert parents vector to ancestral matrix Z.

    Z[i,j] = 1 if node i is ancestor of node j (or i == j).
    Nodes are 0..K where 0 is root, 1..K are clones.
    """
    K = len(parents)
    n = K + 1  # include root
    Z = np.eye(n, dtype=float)

    for i in range(K):
        node = i + 1
        ancestor = parents[i]
        while ancestor >= 0:
            Z[ancestor, node] = 1.0
            if ancestor == 0:
                break
            ancestor = parents[ancestor - 1]

    return Z


# ── Population frequencies ───────────────────────────────────────────────────

def sample_eta(K, S, alpha=0.1, rng=None):
    """Sample subclone proportions from Dirichlet(alpha).

    Returns eta of shape (K+1, S) — includes root (normal) population.
    """
    if rng is None:
        rng = np.random.default_rng()

    n_pops = K + 1
    eta = np.zeros((n_pops, S))
    for s in range(S):
        raw = rng.dirichlet(np.full(n_pops, alpha))
        raw = np.maximum(raw, 1e-30)
        raw /= raw.sum()
        eta[:, s] = raw

    return eta


def compute_phi(Z, eta):
    """Compute cellular prevalences: phi = Z @ eta.

    phi[i, s] = fraction of cells in sample s that carry mutations from clone i.
    """
    return Z @ eta


# ── SSM generation ───────────────────────────────────────────────────────────

def assign_ssms(M, K, rng=None):
    """Assign M SSMs to K clones (1..K), guaranteeing each clone gets >= 1.

    Returns array of length M with values in [1, K].
    """
    if rng is None:
        rng = np.random.default_rng()

    if M < K:
        raise ValueError(f"Need at least K={K} mutations, got M={M}")

    # guarantee one per clone
    assignments = list(range(1, K + 1))
    # distribute remaining
    remaining = M - K
    if remaining > 0:
        weights = rng.dirichlet(np.ones(K))
        extra = rng.choice(range(1, K + 1), size=remaining, p=weights)
        assignments.extend(extra.tolist())

    rng.shuffle(assignments)
    return np.array(assignments)


def generate_reads(phi, assignments, T, omega=0.5, rng=None):
    """Generate binomial read counts for each SSM.

    Args:
        phi: (K+1, S) cellular prevalences
        assignments: (M,) clone assignment per SSM (1-indexed)
        T: total read depth per SSM per sample
        omega: variant read probability given variant present (0.5 for diploid)

    Returns:
        var_reads: (M, S) variant read counts
        total_reads: (M, S) total read counts (constant = T)
    """
    if rng is None:
        rng = np.random.default_rng()

    M = len(assignments)
    S = phi.shape[1]
    var_reads = np.zeros((M, S), dtype=int)
    total_reads = np.full((M, S), T, dtype=int)

    for m in range(M):
        clone = assignments[m]
        for s in range(S):
            p = omega * phi[clone, s]
            p = np.clip(p, 0.0, 1.0)
            var_reads[m, s] = rng.binomial(T, p)

    return var_reads, total_reads


# ── Output writers ───────────────────────────────────────────────────────────

def write_ssm_file(path, var_reads, total_reads, mu_r=0.999, mu_v=0.5):
    """Write PhyloWGS ssm_data.txt format."""
    M, S = var_reads.shape
    with open(path, "w") as f:
        f.write("id\tgene\ta\td\tmu_r\tmu_v\n")
        for m in range(M):
            sid = f"s{m}"
            gene = f"SIM_{m}"
            # 'a' = reference reads = total - variant (PhyloWGS convention)
            a_vals = ",".join(str(total_reads[m, s] - var_reads[m, s]) for s in range(S))
            d_vals = ",".join(str(total_reads[m, s]) for s in range(S))
            f.write(f"{sid}\t{gene}\t{a_vals}\t{d_vals}\t{mu_r}\t{mu_v}\n")


def write_cnv_file(path):
    """Write empty CNV file (SSM-only simulation)."""
    with open(path, "w") as f:
        f.write("cnv\ta\td\tssms\tphysical_cnvs\n")


def write_truth(path, parents, phi, eta, assignments, K, S, T, M, alpha, seed):
    """Write ground truth JSON."""
    # Build clusters: list of lists, clusters[k] = [ssm_ids in clone k+1]
    clusters = [[] for _ in range(K)]
    for m, clone in enumerate(assignments):
        clusters[clone - 1].append(f"s{m}")

    truth = {
        "parents": parents.tolist(),
        "phi": phi.tolist(),
        "eta": eta.tolist(),
        "assignments": assignments.tolist(),
        "clusters": clusters,
        "params": {
            "K": K,
            "S": S,
            "T": T,
            "M": M,
            "alpha": alpha,
            "seed": seed,
            "omega": 0.5,
        },
    }
    with open(path, "w") as f:
        json.dump(truth, f, indent=2)


# ── Main generation logic ────────────────────────────────────────────────────

def generate_fixture(outdir, K, S, T, M, alpha=0.1, tree_type="random", seed=None):
    """Generate one complete simulation fixture."""
    rng = np.random.default_rng(seed)
    os.makedirs(outdir, exist_ok=True)

    parents = make_parents(K, tree_type=tree_type, rng=rng)
    Z = parents_to_ancestral_matrix(parents)
    eta = sample_eta(K, S, alpha=alpha, rng=rng)
    phi = compute_phi(Z, eta)
    assignments = assign_ssms(M, K, rng=rng)
    var_reads, total_reads = generate_reads(phi, assignments, T, rng=rng)

    write_ssm_file(os.path.join(outdir, "ssm_data.txt"), var_reads, total_reads)
    write_cnv_file(os.path.join(outdir, "cnv_data.txt"))
    write_truth(
        os.path.join(outdir, "truth.json"),
        parents, phi, eta, assignments, K, S, T, M, alpha, seed,
    )

    return parents, phi, assignments


# ── Parameter grids ──────────────────────────────────────────────────────────

GRIDS = {
    "quick": {
        "K": [3, 5],
        "S": [1, 3],
        "T": [200, 1000],
        "M_per_K": [10],
        "replicates": 2,
    },
    "default": {
        "K": [3, 5, 10],
        "S": [1, 3],
        "T": [200, 1000],
        "M_per_K": [10, 50],
        "replicates": 3,
    },
    "thorough": {
        "K": [3, 5, 10, 20],
        "S": [1, 3, 10],
        "T": [50, 200, 1000],
        "M_per_K": [10, 50],
        "replicates": 4,
    },
}


def generate_grid(outdir, grid_name="default", base_seed=42):
    """Generate fixtures across a parameter grid."""
    grid = GRIDS[grid_name]
    rng = np.random.default_rng(base_seed)
    manifest = []

    for K in grid["K"]:
        for S in grid["S"]:
            for T in grid["T"]:
                for M_per_K in grid["M_per_K"]:
                    M = K * M_per_K
                    for rep in range(grid["replicates"]):
                        seed = int(rng.integers(0, 2**31))
                        name = f"K{K}_S{S}_T{T}_M{M}_rep{rep}"
                        fixture_dir = os.path.join(outdir, name)
                        print(f"  Generating {name} (seed={seed})")
                        generate_fixture(
                            fixture_dir, K=K, S=S, T=T, M=M, seed=seed,
                        )
                        manifest.append({
                            "name": name,
                            "K": K, "S": S, "T": T, "M": M,
                            "seed": seed, "rep": rep,
                        })

    manifest_path = os.path.join(outdir, "manifest.json")
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"\n  Wrote {len(manifest)} fixtures to {outdir}")
    print(f"  Manifest: {manifest_path}")
    return manifest


# ── CLI ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Generate simulated PhyloWGS inputs with known ground truth"
    )
    parser.add_argument("--outdir", required=True, help="Output directory for fixtures")

    # Grid mode
    parser.add_argument(
        "--grid", choices=list(GRIDS.keys()),
        help="Use a predefined parameter grid (quick/default/thorough)",
    )
    parser.add_argument("--base-seed", type=int, default=42, help="Base seed for grid")

    # Single fixture mode
    parser.add_argument("-K", type=int, help="Number of clones")
    parser.add_argument("-S", type=int, help="Number of samples/timepoints")
    parser.add_argument("-T", type=int, help="Read depth per SSM per sample")
    parser.add_argument("-M", type=int, help="Number of SSMs")
    parser.add_argument("--alpha", type=float, default=0.1, help="Dirichlet concentration")
    parser.add_argument("--tree-type", default="random", choices=["random", "monoprimary", "polyprimary"])
    parser.add_argument("--seed", type=int, default=42, help="Random seed")

    args = parser.parse_args()

    if args.grid:
        print(f"Generating '{args.grid}' grid:")
        generate_grid(args.outdir, args.grid, args.base_seed)
    elif args.K and args.S and args.T and args.M:
        print(f"Generating single fixture: K={args.K} S={args.S} T={args.T} M={args.M}")
        generate_fixture(
            args.outdir, K=args.K, S=args.S, T=args.T, M=args.M,
            alpha=args.alpha, tree_type=args.tree_type, seed=args.seed,
        )
    else:
        parser.error("Provide --grid or all of -K -S -T -M")


if __name__ == "__main__":
    main()
