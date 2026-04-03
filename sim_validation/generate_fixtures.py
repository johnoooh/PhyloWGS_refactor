#!/usr/bin/env python3
"""
generate_fixtures.py — Generate simulated PhyloWGS inputs with known ground truth.

Implements the pearsim generative model:
  1. Random tree topology (parents vector)
  2. Dirichlet-sampled subclone proportions (eta)
  3. Cellular prevalences via ancestral matrix (phi = Z @ eta)
  4. Binomial-sampled read counts
  5. Optional CNV events with timing-aware read generation

Outputs per fixture:
  - ssm_data.txt       (PhyloWGS SSM input format)
  - cnv_data.txt       (CNV file — empty for SSM-only, populated for CNV fixtures)
  - truth.json         (ground truth: tree, phi, clusters, CNVs, params)

Usage:
    python generate_fixtures.py --outdir fixtures/ --grid default
    python generate_fixtures.py --outdir fixtures/ -K 5 -S 3 -T 1000 -M 50 --seed 42
    python generate_fixtures.py --outdir fixtures/ -K 5 -S 3 -T 1000 -M 50 --N-cnv 2
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
        root_children = [i for i in range(K) if parents[i] == 0]
        if len(root_children) > 1:
            for idx in root_children[1:]:
                parents[idx] = 1
    elif tree_type == "polyprimary":
        root_children = [i for i in range(K) if parents[i] == 0]
        if len(root_children) < 2 and K >= 2:
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
    """Generate binomial read counts for SSMs without CNVs.

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


# ── CNV generation ───────────────────────────────────────────────────────────

# Common copy number states: (maternal_cn, paternal_cn), weight
CN_STATES = [
    ((1, 2), 3.0),  # single-copy gain
    ((2, 1), 3.0),  # single-copy gain (other allele)
    ((0, 1), 2.0),  # hemizygous deletion
    ((1, 0), 2.0),  # hemizygous deletion (other allele)
    ((0, 2), 1.5),  # CN-neutral LOH
    ((2, 0), 1.5),  # CN-neutral LOH (other allele)
    ((1, 3), 1.0),  # higher amplification
]


def generate_cnv_events(N_cnv, K, M, assignments, Z, cnv_overlap_frac=0.3, rng=None):
    """Generate CNV events with known tree relationships.

    Each CNV is assigned to a clone node, gets a copy number state, and
    affects a subset of SSMs. Timing relationships (ssm_before_cnv,
    cnv_before_ssm, same_node) are determined by the tree structure.

    Returns list of CNV dicts.
    """
    if rng is None:
        rng = np.random.default_rng()

    if N_cnv == 0:
        return []

    # Sample CN states
    states = [s for s, _ in CN_STATES]
    weights = np.array([w for _, w in CN_STATES])
    weights /= weights.sum()
    state_indices = rng.choice(len(states), size=N_cnv, p=weights)

    # Assign CNVs to clones (1..K)
    cnv_clones = rng.integers(1, K + 1, size=N_cnv)

    # Track which SSMs are already claimed by a CNV (max 1 CNV per SSM)
    ssm_claimed = set()

    # Generate fake genomic coordinates
    chroms = rng.integers(1, 23, size=N_cnv)

    cnvs = []
    for c in range(N_cnv):
        mat_cn, pat_cn = states[state_indices[c]]
        cnv_clone = int(cnv_clones[c])

        # Select SSMs to be affected — pick from unclaimed SSMs
        available = [m for m in range(M) if m not in ssm_claimed]
        n_affect = max(1, int(len(available) * cnv_overlap_frac))
        n_affect = min(n_affect, len(available))
        affected_indices = rng.choice(available, size=n_affect, replace=False).tolist()

        for idx in affected_indices:
            ssm_claimed.add(idx)

        # Classify timing for each affected SSM
        affected_ssms = []
        for m in affected_indices:
            ssm_clone = int(assignments[m])
            if ssm_clone == cnv_clone:
                timing = "same_node"
            elif Z[ssm_clone, cnv_clone] == 1.0:
                timing = "ssm_before_cnv"
            elif Z[cnv_clone, ssm_clone] == 1.0:
                timing = "cnv_before_ssm"
            else:
                # Neither is ancestor of the other — CNV doesn't affect
                # cells with this SSM in the Go model. The SSM is still
                # listed as affected in the file (CNV overlaps the locus),
                # but the timing means only some tree nodes see both.
                timing = "independent"

            affected_ssms.append({
                "ssm_idx": m,
                "ssm_id": f"s{m}",
                "ssm_clone": ssm_clone,
                "maternal_cn": mat_cn,
                "paternal_cn": pat_cn,
                "timing": timing,
            })

        start = int(rng.integers(1_000_000, 100_000_000))
        length = int(rng.integers(500_000, 5_000_000))

        cnvs.append({
            "id": f"c{c}",
            "clone": cnv_clone,
            "maternal_cn": mat_cn,
            "paternal_cn": pat_cn,
            "affected_ssms": affected_ssms,
            "chrom": int(chroms[c]),
            "start": start,
            "end": start + length,
        })

    return cnvs


def compute_cnv_aware_vaf(ssm_clone, cnv_clone, mat_cn, pat_cn, eta, Z, s):
    """Compute (nr, nv) pairs for a single SSM+CNV in one sample.

    Mirrors computeNGenomes (main.go:756-863) exactly.
    Uses eta[node, s] as the per-node mixture weight (pi).

    Returns list of (nr, nv) tuples — 4 if ssm_clone == cnv_clone, else 2.
    """
    n_nodes = eta.shape[0]  # K+1 including root
    nr1, nv1, nr2, nv2 = 0.0, 0.0, 0.0, 0.0
    nr3, nv3, nr4, nv4 = 0.0, 0.0, 0.0, 0.0

    for node in range(n_nodes):
        pi = eta[node, s]

        ssm_in_ancestors = Z[ssm_clone, node] == 1.0
        cnv_in_ancestors = Z[cnv_clone, node] == 1.0

        if not ssm_in_ancestors and not cnv_in_ancestors:
            # Case 1: no variant, no CNV — diploid
            nr1 += pi * 2
            nr2 += pi * 2
            nr3 += pi * 2
            nr4 += pi * 2
        elif ssm_in_ancestors and not cnv_in_ancestors:
            # Case 2: has variant, no CNV — heterozygous
            nr1 += pi; nv1 += pi
            nr2 += pi; nv2 += pi
            nr3 += pi; nv3 += pi
            nr4 += pi; nv4 += pi
        elif not ssm_in_ancestors and cnv_in_ancestors:
            # Case 3: no variant, has CNV — all copies reference
            total_cn = mat_cn + pat_cn
            nr1 += pi * total_cn
            nr2 += pi * total_cn
            nr3 += pi * total_cn
            nr4 += pi * total_cn
        else:
            # Case 4: has variant AND has CNV
            total_cn = mat_cn + pat_cn

            # Pairs 3,4: SSM occurred AFTER CNV (variant not amplified)
            nr3 += pi * max(0, total_cn - 1)
            nv3 += pi * min(1, total_cn)
            nr4 += pi * max(0, total_cn - 1)
            nv4 += pi * min(1, total_cn)

            # Pairs 1,2: check timing
            ssm_before_cnv = Z[ssm_clone, cnv_clone] == 1.0
            if ssm_before_cnv:
                # Variant gets amplified/deleted by CNV
                nr1 += pi * pat_cn; nv1 += pi * mat_cn
                nr2 += pi * mat_cn; nv2 += pi * pat_cn
            else:
                # CNV before or same node — variant not amplified
                nr1 += pi * max(0, total_cn - 1)
                nv1 += pi * min(1, total_cn)
                nr2 += pi * max(0, total_cn - 1)
                nv2 += pi * min(1, total_cn)

    if ssm_clone == cnv_clone:
        return [(nr1, nv1), (nr2, nv2), (nr3, nv3), (nr4, nv4)]
    return [(nr1, nv1), (nr2, nv2)]


def compute_expected_ref_frac(nr_nv_pairs, mu_r=0.999):
    """Compute expected reference read fraction from (nr, nv) pairs.

    Mirrors logLikelihoodWithCNVTree (main.go:865-914): averages over
    valid pairs (nv > 0) with uniform prior.
    """
    valid = [(nr, nv) for nr, nv in nr_nv_pairs if nv > 0]
    if not valid:
        return mu_r  # no variant signal — all reference

    mus = []
    for nr, nv in valid:
        total = nr + nv
        if total < 1e-15:
            total = 1e-15
        mu = (nr * mu_r + nv * (1 - mu_r)) / total
        mus.append(mu)

    return float(np.mean(mus))


def generate_reads_with_cnv(phi, eta, Z, assignments, cnvs, T, S, mu_r=0.999,
                            rng=None):
    """Generate read counts for all SSMs, accounting for CNV effects.

    SSMs without CNVs use simple Binomial(T, 0.5 * phi).
    SSMs with CNVs use the full computeNGenomes model.

    Returns:
        ref_reads: (M, S) reference read counts
        total_reads: (M, S) total read counts
    """
    if rng is None:
        rng = np.random.default_rng()

    M = len(assignments)
    ref_reads = np.zeros((M, S), dtype=int)
    total_reads = np.full((M, S), T, dtype=int)

    # Build SSM index → CNV mapping
    ssm_cnv_map = {}  # ssm_idx → (cnv_clone, mat_cn, pat_cn)
    for cnv in cnvs:
        for aff in cnv["affected_ssms"]:
            ssm_cnv_map[aff["ssm_idx"]] = (
                cnv["clone"],
                aff["maternal_cn"],
                aff["paternal_cn"],
            )

    for m in range(M):
        clone = int(assignments[m])
        for s in range(S):
            if m in ssm_cnv_map:
                cnv_clone, mat_cn, pat_cn = ssm_cnv_map[m]
                pairs = compute_cnv_aware_vaf(clone, cnv_clone, mat_cn, pat_cn,
                                             eta, Z, s)
                mu = compute_expected_ref_frac(pairs, mu_r)
                mu = np.clip(mu, 1e-6, 1 - 1e-6)
                ref_reads[m, s] = rng.binomial(T, mu)
            else:
                # Simple diploid model: ref_frac = 1 - 0.5 * phi
                p_ref = 1.0 - 0.5 * phi[clone, s]
                p_ref = np.clip(p_ref, 0.0, 1.0)
                ref_reads[m, s] = rng.binomial(T, p_ref)

    return ref_reads, total_reads


def generate_cnv_reads(cnvs, phi, T_cnv, S, rng=None):
    """Generate read counts for CNV regions.

    Models aggregate heterozygous SNP coverage across the CNV region.
    ref_frac = 1 - phi[cnv_clone, s] * (1 - minor_cn / total_cn)
    """
    if rng is None:
        rng = np.random.default_rng()

    for cnv in cnvs:
        clone = cnv["clone"]
        mat_cn = cnv["maternal_cn"]
        pat_cn = cnv["paternal_cn"]
        total_cn = mat_cn + pat_cn
        minor_cn = min(mat_cn, pat_cn)

        a = np.zeros(S, dtype=int)
        d = np.full(S, T_cnv, dtype=int)

        for s in range(S):
            if total_cn > 0:
                ref_frac = 1.0 - phi[clone, s] * (1.0 - minor_cn / total_cn)
            else:
                ref_frac = 1.0
            ref_frac = np.clip(ref_frac, 0.01, 0.99)
            a[s] = rng.binomial(T_cnv, ref_frac)

        cnv["a"] = a.tolist()
        cnv["d"] = d.tolist()


# ── Output writers ───────────────────────────────────────────────────────────

def write_ssm_file(path, ref_reads, total_reads, mu_r=0.999, mu_v=0.5):
    """Write PhyloWGS ssm_data.txt format."""
    M, S = ref_reads.shape
    with open(path, "w") as f:
        f.write("id\tgene\ta\td\tmu_r\tmu_v\n")
        for m in range(M):
            sid = f"s{m}"
            gene = f"SIM_{m}"
            a_vals = ",".join(str(ref_reads[m, s]) for s in range(S))
            d_vals = ",".join(str(total_reads[m, s]) for s in range(S))
            f.write(f"{sid}\t{gene}\t{a_vals}\t{d_vals}\t{mu_r}\t{mu_v}\n")


def write_cnv_file(path, cnvs=None):
    """Write cnv_data.txt — empty header for SSM-only, full data for CNV fixtures."""
    with open(path, "w") as f:
        f.write("cnv\ta\td\tssms\tphysical_cnvs\n")
        if not cnvs:
            return
        for cnv in cnvs:
            a_str = ",".join(str(x) for x in cnv["a"])
            d_str = ",".join(str(x) for x in cnv["d"])
            ssm_parts = []
            for aff in cnv["affected_ssms"]:
                ssm_parts.append(
                    f"{aff['ssm_id']},{aff['maternal_cn']},{aff['paternal_cn']}"
                )
            ssm_str = ";".join(ssm_parts)
            phys = (f"chrom={cnv['chrom']},start={cnv['start']},"
                    f"end={cnv['end']},"
                    f"major_cn={max(cnv['maternal_cn'], cnv['paternal_cn'])},"
                    f"minor_cn={min(cnv['maternal_cn'], cnv['paternal_cn'])},"
                    f"cell_prev=0.0000")
            f.write(f"{cnv['id']}\t{a_str}\t{d_str}\t{ssm_str}\t{phys}\n")


def write_truth(path, parents, phi, eta, assignments, K, S, T, M, alpha, seed,
                cnvs=None, N_cnv=0, T_cnv=0):
    """Write ground truth JSON."""
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
            "N_cnv": N_cnv,
            "T_cnv": T_cnv,
            "has_cnvs": N_cnv > 0,
        },
    }

    if cnvs:
        truth["cnvs"] = []
        for cnv in cnvs:
            truth["cnvs"].append({
                "id": cnv["id"],
                "clone": cnv["clone"],
                "maternal_cn": cnv["maternal_cn"],
                "paternal_cn": cnv["paternal_cn"],
                "chrom": cnv["chrom"],
                "start": cnv["start"],
                "end": cnv["end"],
                "affected_ssms": [
                    {
                        "ssm_id": a["ssm_id"],
                        "ssm_clone": a["ssm_clone"],
                        "maternal_cn": a["maternal_cn"],
                        "paternal_cn": a["paternal_cn"],
                        "timing": a["timing"],
                    }
                    for a in cnv["affected_ssms"]
                ],
            })

    with open(path, "w") as f:
        json.dump(truth, f, indent=2)


# ── Main generation logic ────────────────────────────────────────────────────

def generate_fixture(outdir, K, S, T, M, alpha=0.1, tree_type="random", seed=None,
                     N_cnv=0, cnv_overlap_frac=0.3, T_cnv=10000):
    """Generate one complete simulation fixture."""
    rng = np.random.default_rng(seed)
    os.makedirs(outdir, exist_ok=True)

    parents = make_parents(K, tree_type=tree_type, rng=rng)
    Z = parents_to_ancestral_matrix(parents)
    eta = sample_eta(K, S, alpha=alpha, rng=rng)
    phi = compute_phi(Z, eta)
    assignments = assign_ssms(M, K, rng=rng)

    cnvs = []
    if N_cnv > 0:
        cnvs = generate_cnv_events(N_cnv, K, M, assignments, Z,
                                   cnv_overlap_frac=cnv_overlap_frac, rng=rng)
        generate_cnv_reads(cnvs, phi, T_cnv, S, rng=rng)
        ref_reads, total_reads = generate_reads_with_cnv(
            phi, eta, Z, assignments, cnvs, T, S, rng=rng,
        )
    else:
        # SSM-only: generate variant reads, then convert to reference reads
        var_reads, total_reads = generate_reads(phi, assignments, T, rng=rng)
        ref_reads = total_reads - var_reads

    write_ssm_file(os.path.join(outdir, "ssm_data.txt"), ref_reads, total_reads)
    write_cnv_file(os.path.join(outdir, "cnv_data.txt"), cnvs if cnvs else None)
    write_truth(
        os.path.join(outdir, "truth.json"),
        parents, phi, eta, assignments, K, S, T, M, alpha, seed,
        cnvs=cnvs, N_cnv=N_cnv, T_cnv=T_cnv,
    )

    return parents, phi, assignments


# ── Parameter grids ──────────────────────────────────────────────────────────

GRIDS = {
    "quick": {
        "K": [3, 5],
        "S": [1, 3],
        "T": [200, 1000],
        "M_per_K": [10],
        "N_cnv": [0, 2],
        "replicates": 2,
    },
    "default": {
        "K": [3, 5, 10],
        "S": [1, 3],
        "T": [200, 1000],
        "M_per_K": [10, 50],
        "N_cnv": [0, 2, 5],
        "replicates": 3,
    },
    "thorough": {
        "K": [3, 5, 10, 20],
        "S": [1, 3, 10],
        "T": [50, 200, 1000],
        "M_per_K": [10, 50],
        "N_cnv": [0, 2, 5, 10],
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
                    for N_cnv in grid["N_cnv"]:
                        for rep in range(grid["replicates"]):
                            seed = int(rng.integers(0, 2**31))
                            name = f"K{K}_S{S}_T{T}_M{M}_C{N_cnv}_rep{rep}"
                            fixture_dir = os.path.join(outdir, name)
                            print(f"  Generating {name} (seed={seed})")
                            generate_fixture(
                                fixture_dir, K=K, S=S, T=T, M=M,
                                N_cnv=N_cnv, seed=seed,
                            )
                            manifest.append({
                                "name": name,
                                "K": K, "S": S, "T": T, "M": M,
                                "N_cnv": N_cnv,
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
    parser.add_argument("--tree-type", default="random",
                        choices=["random", "monoprimary", "polyprimary"])
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--N-cnv", type=int, default=0, help="Number of CNV events")
    parser.add_argument("--cnv-overlap-frac", type=float, default=0.3,
                        help="Fraction of SSMs affected per CNV")
    parser.add_argument("--T-cnv", type=int, default=10000,
                        help="Read depth for CNV regions")

    args = parser.parse_args()

    if args.grid:
        print(f"Generating '{args.grid}' grid:")
        generate_grid(args.outdir, args.grid, args.base_seed)
    elif args.K and args.S and args.T and args.M:
        print(f"Generating single fixture: K={args.K} S={args.S} T={args.T} M={args.M}"
              f" N_cnv={args.N_cnv}")
        generate_fixture(
            args.outdir, K=args.K, S=args.S, T=args.T, M=args.M,
            alpha=args.alpha, tree_type=args.tree_type, seed=args.seed,
            N_cnv=args.N_cnv, cnv_overlap_frac=args.cnv_overlap_frac,
            T_cnv=args.T_cnv,
        )
    else:
        parser.error("Provide --grid or all of -K -S -T -M")


if __name__ == "__main__":
    main()
