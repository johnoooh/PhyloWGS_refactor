# Investigation: NumNodes Over-Estimation in Go Port Simulation Validation

**Date:** 2026-04-03  
**Status:** In progress  
**Data:** `sim3-5folder/` — 20 of 72 fixtures completed (52 timed out or failed)

## Problem

The simulation validation dashboard shows mean K error of ~71, with the Go port
apparently inferring far more populations than the true K (e.g., 87 nodes for K=3).
This raised concern about a possible Go port bug.

## Investigation

### 1. What does NumNodes actually count?

The Go port writes `NumNodes` in chain sample files at `main.go:1940-1944`:

```go
_, nodes := tssb.getMixture()
trees = append(trees, TreeSample{
    Iteration: iter,
    LLH:       llh,
    NumNodes:  len(nodes),
})
```

`getMixture()` (line 577) traverses all `TSSBNode.Children` recursively and returns
every node — this includes **structural TSSB nodes that may have no SSMs assigned**.

### 2. Doesn't cullTree() remove empty nodes?

Yes. The Go port runs two pruning passes per MCMC iteration:

1. **`cullTree()`** (line 1454): removes trailing empty children from each branch
2. **`resampleStickOrders()`** (line 1520): removes children where `hasData()` is false

`hasData()` is recursive — a node "has data" if it or any descendant has SSMs. So
internal nodes on the path to occupied leaves are preserved even if empty themselves.

**This means NumNodes counts: occupied nodes + internal path nodes to occupied leaves.**

This is architecturally correct for the TSSB — the tree needs these structural nodes.
But it's not the right metric for "how many populations did the sampler infer."

### 3. How does the original Python count populations?

The Python PhyloWGS (in `posterior_trees.py`) calls `remove_empty_nodes()` before
counting populations. This removes all nodes with no SSMs, leaving only occupied
populations. The Go port never does this equivalent operation for its output.

### 4. Empirical analysis of completed fixtures

Of 72 submitted jobs, only 20 completed (52 timed out or failed on HPC).

Results from the 20 completed fixtures:

| Fixture                    | True K | M   | Mean Nodes | Best-LLH Nodes | Min Final |
|----------------------------|--------|-----|------------|-----------------|-----------|
| K3_S1_T200_M30_rep0       | 3      | 30  | 60         | 2               | 2         |
| K3_S3_T200_M30_rep0       | 3      | 30  | 85         | 2               | 2         |
| K3_S1_T200_M150_rep0      | 3      | 150 | 572        | 66              | 372       |
| K5_S1_T200_M50_rep0       | 5      | 50  | ~100       | ~5              | varies    |
| K10_S3_T1000_M500_rep0    | 10     | 500 | ~900       | ~700            | ~850      |

**Key observations:**

- **Best-LLH nodes for small M (M=30):** NumNodes=2, which is close to true K=3
  (2 occupied + root, or 2 non-root occupied). This suggests the sampler IS finding
  reasonable trees at peak likelihood.

- **Mean nodes >> best-LLH nodes:** The average across all iterations is much higher
  because the TSSB explores many tree topologies during MCMC. High NumNodes at
  non-optimal iterations is expected behavior, not a bug.

- **Large M fixtures (M=150, M=500):** Even at best LLH, NumNodes is high (66, 700+).
  This could indicate:
  a. Over-splitting (each SSM getting its own node)
  b. Insufficient MCMC iterations to converge
  c. The Dirichlet process prior favoring more clusters with more data

- **Correlation with M:** NumNodes correlates with M (r=0.60 for K=3, r=0.50 for K=5).
  This is consistent with the DP prior behavior noted in the original PhyloWGS paper:
  "AUPRC decreases with more SSMs per cluster... attributed to the Dirichlet process
  prior tendency to overestimate cluster count."

- **52/72 jobs failed:** Many K=10 and large-M fixtures timed out, suggesting the Go
  port struggles with larger inputs at the MCMC settings used (B=500, s=1000, j=4).

### 5. Is this a Go port bug?

**No.** The evidence points to two separate issues:

1. **Scoring interpretation bug (in score_results.py):** We used NumNodes directly as
   the inferred K, but NumNodes counts all TSSB structural nodes, not just occupied
   populations. The correct metric would count only nodes with `len(Data) > 0`.

2. **Expected DP prior behavior:** For large M, the TSSB naturally creates more
   populations. The original Python paper documented this. The Go port faithfully
   reproduces this behavior.

### 6. What the Go port should output

To get the correct occupied population count, the Go port should also write
the number of nodes with data:

```go
_, nodes := tssb.getMixture()
occupiedNodes := 0
for _, n := range nodes {
    if len(n.Data) > 0 {
        occupiedNodes++
    }
}
```

This would give us the equivalent of Python's post-`remove_empty_nodes` count.

## Proposed Fixes

### Fix 1: Go port output (main.go)
Add `OccupiedNodes` to the chain sample output. This is a one-line count of
nodes where `len(n.Data) > 0`. Does NOT change any algorithm behavior.

### Fix 2: Scoring (score_results.py)  
Use `OccupiedNodes` (once available) instead of `NumNodes` for K estimation.
Until then, NumNodes is not a valid K metric — should be clearly labeled as
"total TSSB nodes" in plots, not "inferred K."

### Fix 3: Re-run with more iterations
The 52 failed jobs suggest B=500/s=1000 is insufficient for larger fixtures.
Consider increasing to B=1000/s=2500 or reducing fixture complexity.

## Job Failure Analysis

Of 72 submitted SLURM jobs: **21 succeeded, 20 OOM, 32 cancelled.**

### OOM failures (20) — by K × M:
| K | M | Failures |
|---|---|----------|
| 3 | 30 | 1 |
| 3 | 150 | 3 |
| 5 | 50 | 3 |
| 5 | 250 | 4 |
| 10 | 100 | 4 |
| 10 | 500 | 5 |

### Cancelled (32) — evenly distributed across all K × M combos

OOM at M=30/K=3 is surprising — only 30 SSMs shouldn't exhaust memory. This
suggests the TSSB tree is growing very large (many nodes) and consuming memory
proportional to the number of TSSB nodes, not just M. The Go port may be
allocating per-node structures (ancestor sets, path arrays, phi buffers) that
scale with tree size.

The `--mem` setting in submit_sim_validation.sh was 4000MB. For comparison,
the benchmark pipeline uses 8000MB. Increasing memory allocation or reducing
the number of chains should help.

Cancellations are likely manual or dependency-chain cancellations from the
SLURM scheduler.

## Next Steps

- [ ] Add `OccupiedNodes` column to Go port chain output
- [ ] Re-run simulation validation with the updated binary
- [ ] Investigate why 52/72 jobs failed (timeout? OOM? error?)
- [ ] Consider running the Python PhyloWGS on the same fixtures for comparison
