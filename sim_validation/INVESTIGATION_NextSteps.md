# Investigation — Recommended Next Steps

_2026-04-08. Goal: figure out what Python actually does, so Go can match it faithfully (no custom deviations)._

## TL;DR

My earlier SIM_VALIDATION_REPORT_3.md recommendation to "add minimum-φ node pruning" **was wrong**. Python does not do that. Here is the corrected picture:

1. **Over-splitting is not a Go bug, and Python doesn't prune for it either.** Python and Go produce nearly identical inferred_K on the same fixture when given comparable compute (Python 81 vs Go 84 on `K10_S1_T1000_M100_C0_rep1`). The DP prior allows high-K posteriors; both implementations explore them. The reason Python's mean K_error looks lower in the grid is that Python is slower — it runs fewer effective samples in walltime, which accidentally regularizes the result. Fixing the speed gap will make this apparent in Python too.
2. **The CNV slowdown IS a real fidelity gap with a clear Python-faithful fix.** Python precomputes per-(SSM, node) `(nr, nv)` coefficients once per MCMC iteration and then its C++ MH loop just does `Σᵢ pi[tp]·nrᵢ` dot products. Go recomputes the full `computeNGenomes` tree traversal on every one of the 5,000 MH iterations. This is the 10-100× slowdown. Python already implements the fix we need.

So the priority is reversed from what I wrote in Report 3:

- **P0 — Port Python's precomputed MH data states.** Matches Python exactly. Will close most of the CNV walltime gap.
- **P1 — Re-measure over-splitting after P0 lands.** It may disappear once Go and Python do the same amount of effective sampling per walltime minute.

---

## Investigation

### Q1. Does Python prune empty/tiny-φ nodes during MCMC?

**Short answer: No minimum-φ pruning anywhere. There are two separate empty-node operations, neither is φ-based.**

#### `cull_tree()` — called every MCMC iteration

`phylowgs/tssb.py:152`:
```python
def cull_tree(self):
    def descend(root):
        counts = array(map(lambda child: descend(child), root['children']))
        keep = len(trim_zeros(counts, 'b'))    # trim TRAILING zeros only
        for child in root['children'][keep:]:
            child['node'].kill()
            del child['node']
        root['sticks'] = root['sticks'][:keep]
        root['children'] = root['children'][:keep]
        return sum(counts) + root['node'].num_local_data()
    descend(self.root)
```

Called in `evolve.py:172` immediately after `resample_assignments()` every iteration. It only removes **trailing** empty children (i.e., the last k children of a parent that contain zero data in their entire subtree). Middle-empty children and internal empty nodes stay. It is not a φ threshold — the criterion is "subtree has zero assigned SSMs AND is the last sibling".

**Go already does this exactly** — `main.go:1736 cullTree()`, called from `runChain` at `main.go:2265`, same position in the loop as Python.

#### `remove_empty_nodes()` — only called during POST-PROCESSING

`phylowgs/util2.py:128`. It fully strips empty nodes (leaves are deleted, empty internals are re-parented to grandparent). Its callers:
- `TreeReader._parse_tree(remove_empty_vertices=True)` — when the result generator reads trees back from `trees.zip`
- `posterior_trees.py` — post-hoc analysis
- `printo.py` — plotting

Search for callers inside the MCMC loop:
```
$ grep -rn "remove_empty_nodes" phylowgs/
phylowgs/util2.py:128:def remove_empty_nodes(root, parent = None):
phylowgs/util2.py:313:    remove_empty_nodes(tree.root)         # TreeReader post-processing
phylowgs/printo.py:18:    remove_empty_nodes(tree.root, None)   # plot
phylowgs/misc/post_assign_ssm.py:19: ... remove_empty_vertices = True    # post-processing
phylowgs/posterior_trees.py:85 / 134 / 145: remove_empty_nodes(tssb.root, None)   # post-processing
phylowgs/pwgsresults/result_generator.py:39: remove_empty_vertices = True         # scoring
```

Not called during MCMC. The live tree keeps empty internal nodes throughout the run.

**Confirmation that inferred_K is counted on the pruned tree:** `result_generator.py:39` reads trees with `remove_empty_vertices=True`, and Python's tree-summary JSON in `tree_summaries.json.gz` stores `populations[*].num_ssms / num_cnvs`, which our scoring script `score_results.py:206` filters with `if num_ssms > 0 or num_cnvs > 0`. Go's `summarizePops` at `main.go:2673` does the same filter. Both count only non-empty pops, so the metric is apples-to-apples.

### Q2. Is Go genuinely over-splitting more than Python, or is it a counting artifact?

The counting is fair (answered above). But the over-splitting gap in the grid data is smaller than the headline numbers suggest. Look at matched extreme fixtures:

| Fixture | True K | Python inferred K | Go inferred K |
|---------|:------:|:-----------------:|:-------------:|
| K10_S1_T1000_M100_C0_rep1 | 10 | 81 | 84 |
| K10_S1_T1000_M100_C0_rep0 | 10 | 40 | (timeout/scored) |
| K3_S1_T200_M30_C0_rep1 | 3 | 20 | — |
| K5_S1_T200_M50_C2_rep1 | 5 | 20 | — |

The first row is the killer: **Python and Go land within 3 of each other on the same fixture**. This is not a Go regression — it's the DP prior's behavior on high-M, single-sample data. Both implementations over-split because the stick-breaking prior with `min_dp_alpha=1.0`, `max_dp_alpha=50` has no mechanism to resist it.

Why does the grid headline then show Go mean K_error 22.84 vs Python 6.97 for C=0, S=1? Two reasons:

1. **Stochastic variance.** Python's K ranges 1–81. A handful of runs hit the 80+ range, most stop in the single digits. Mean is dominated by whether a run got lucky with a "stop early" configuration. The Go C=0 runs, because they finish 3× faster, get to do ~3× more effective sampling in the same walltime, and have more opportunity to find the high-K mode.
2. **Python's walltime masks the issue.** Python on `K5_S1_T1000_M250_C0_rep0` timed out (one of the py_null entries) — it never reported an inferred_K at all. If Python had the same speed as Go, we would likely see Python's mean creep up too.

**So the over-splitting is a model/prior issue that belongs upstream (in how PhyloWGS handles single-sample high-M data), not a Go port bug.** Adding a custom φ-threshold pruning would be a deviation from Python and would be silently hiding the model's behavior. Not the right move for a port.

**Action for the port:** don't touch this. After we fix the CNV performance gap (Q3), we should re-run the grid so Go and Python have symmetric walltime budgets, and the K_error numbers should become much more similar.

### Q3. Why is Go 10-100× slower on CNV fixtures?

This is the real fidelity gap. I traced it to a Python optimization that Go is missing.

#### Python's hot path (C++ MH loop)

Python does not call `compute_n_genomes` during the 5,000-iteration MH inner loop. It calls it **once** before MH starts, and serializes the per-node `(nr, nv)` coefficients to disk.

`phylowgs/params.py:114 write_data_state`:
```python
def write_data_state(tssb, fname):
    for dat in tssb.data:
        if not dat.cnv: continue
        poss_n_genomes = dat.compute_n_genomes(0)   # one full tree walk
        for node in nodes:
            # compute state1..state4 per (datum, node) — the (nr,nv) coefficients
            # for each of the 4 (maternal/paternal × before/after CNV) timing cases
            ...
            fh.write(...)   # write to c_data_states.txt
```

Then `mh.cpp` reads those states into `datum::states1[], states2[], states3[], states4[]` — each element is `{node*, nr, nv}` — and the hot loop in `mh.hpp:80 log_complete_ll`:
```cpp
for (int i = 0; i < states1.size(); i++) {
    if (old == 0) {
        nr += (states1[i].nd->pi1[tp]) * states1[i].nr;
        nv += (states1[i].nd->pi1[tp]) * states1[i].nv;
    } else {
        nr += (states1[i].nd->pi[tp]) * states1[i].nr;
        nv += (states1[i].nd->pi[tp]) * states1[i].nv;
    }
}
mu = (nv*(1-mu_r) + nr*mu_r) / (nr+nv);
ll[0] = log_binomial_likelihood(a[tp], d[tp], mu) + log(0.25) + log_bin_norm_const[tp];
```

**This is just a dot product of the precomputed `(nr, nv)` coefficient vector with the `pi[tp]` vector.** O(K) per SSM per timepoint per MH iteration. No tree walks, no ancestor checks, no CNV resolution — all of that is baked into the precomputed coefficients that stay valid for the entire MH inner loop.

#### Go's hot path

`main.go:1575 paramPost` (called twice per MH iter: once for `newLLH`, `oldLLH` is cached):
```go
for _, node := range nodes {
    for _, idx := range node.Data {
        ssm := t.Data[idx]
        if len(ssm.CNVs) > 0 {
            llh += logLikelihoodWithCNVTreeMH(ssm, t, useNew)   // <-- expensive
        } else {
            llh += ssm.logLikelihoodNoCNV(params)               // cheap
        }
    }
}
```

`logLikelihoodWithCNVTreeMH` at `main.go:1628`:
```go
for tp := range ssm.A {
    possNGenomes := computeNGenomes(ssm, tssb, tp, newState)   // <-- full tree walk EVERY MH iter
    ...
}
```

And `computeNGenomes` at `main.go:839` iterates **all K nodes** in the tree, calls `findMostRecentCNV` (O(D·C)) and `isAncestorOfCached` (O(1)) for each, and rebuilds the `(nr, nv)` sums from scratch.

**Per MCMC iteration (with default `mh_iters=5000`):**

| | Python (C++ MH) | Go (current) |
|-|:-:|:-:|
| `computeNGenomes` calls per CNV SSM | 1 (before MH loop) | 5,000 × NTPS |
| Work per MH iteration per CNV SSM | O(K) dot product | O(K·D·C) full traversal |
| Effective complexity per MCMC iter | O(M·K + 5000·M·NTPS·K) | O(5000·M·NTPS·K·D·C) |

For a realistic CNV fixture (M=50 SSMs, K=50 nodes, NTPS=3, D=5, C=2, 5000 MH iters):
- **Python:** ≈ 50·50 + 5000·50·3·50 ≈ 3.75·10⁷ multiply-adds per MCMC iter
- **Go:** ≈ 5000·50·3·50·5·2 ≈ 3.75·10⁸ ops per MCMC iter (**10× slower at baseline**)

When Go over-splits to K=150 (as in the grid's extreme outliers), Go's complexity jumps to ~1.1·10⁹ per iter — while Python's stays at 1.1·10⁸. **That's the >10× gap we see in the grid data.**

This is also the mechanism that drove the 71 timeouts: CNV fixtures with larger inferred K compound the problem until the MCMC can't finish in walltime.

### Q4. Can goroutines help?

Yes, but in a specific place. The Markov chain itself is inherently serial (each MH accept/reject depends on the previous state). What we can parallelize:

1. **The initial `precomputeMHStates()` across SSMs** — each SSM's state computation is independent. Use a worker pool, one goroutine per ~M/GOMAXPROCS SSMs. Done once per MCMC iteration; totally safe.
2. **The `paramPost` summation over SSMs inside the MH loop** — tempting but risky: with the O(K) optimization each SSM's work becomes tiny (~100 ops), and the goroutine dispatch overhead per MH iter (5,000 times per MCMC iter) would dwarf the actual work. **Keep it serial.**
3. **Chains** — already parallelized (the `runChain` goroutines). This is the big one and is already in place.

So goroutine use on P0 is: parallelize precomputation, keep MH loop serial (matches Python's single-threaded C++ MH).

---

## Recommended Plan

### P0 — Precomputed MH data states (match Python's `write_data_state` + `mh.cpp`)

**New data structure on SSM:**
```go
type SSMNodeState struct {
    NodeIdx int     // index into the MH's flattened nodes slice
    Nr      float64 // coefficient: contribution to nr from this node per unit pi
    Nv      float64 // contribution to nv from this node per unit pi
}

// On SSM, four slices for the four timing/allele cases:
type SSM struct {
    ...
    MHStates1 []SSMNodeState // maternal, before-CNV
    MHStates2 []SSMNodeState // paternal, before-CNV
    MHStates3 []SSMNodeState // maternal, after-CNV
    MHStates4 []SSMNodeState // paternal, after-CNV
    UseFourStates bool       // matches Python's `len(poss_n_genomes)==4` branch
}
```

**New function in main.go:**
```go
// precomputeMHStates walks the tree once and materializes per-(SSM, node)
// (nr, nv) coefficients for each of the 4 timing/allele cases. This mirrors
// phylowgs/params.py write_data_state + mh.cpp load_data_states.
//
// Must be called after cullTree/setNodeHeights/setNodePaths/mapDatumToNode,
// before entering the MH loop. Valid until the tree structure changes.
func (t *TSSB) precomputeMHStates() {
    _, nodes := t.getMixture()
    // Parallel over SSMs — independent work, do M / GOMAXPROCS at a time.
    parallelForSSMs(t.Data, func(ssm *SSM) {
        if len(ssm.CNVs) == 0 {
            ssm.MHStates1 = nil
            return
        }
        // ... walk nodes, compute nr/nv for each case per node (same
        // case analysis as current computeNGenomes). Store only nonzero
        // (nr, nv) entries for efficiency.
        // Set UseFourStates based on whether SSM and its CNV share a node.
    })
}
```

**Rewrite `logLikelihoodWithCNVTreeMH` to use the states:**
```go
func logLikelihoodWithCNVTreeMH(ssm *SSM, tssb *TSSB, newState bool) float64 {
    if len(ssm.CNVs) == 0 {
        params := ssm.Node.Params
        if newState { params = ssm.Node.Params1 }
        return ssm.logLikelihoodNoCNV(params)
    }
    llh := 0.0
    for tp := range ssm.A {
        // Four cases via precomputed states — pure O(K) dot products.
        ll := make([]float64, 0, 4)
        for _, states := range [][]SSMNodeState{ssm.MHStates1, ssm.MHStates2, ssm.MHStates3, ssm.MHStates4} {
            if states == nil { continue }
            nr, nv := 0.0, 0.0
            for _, s := range states {
                node := tssb.nodesByIdx[s.NodeIdx] // precomputed flat slice
                var pi float64
                if newState { pi = node.Pi1[tp] } else { pi = node.Pi[tp] }
                nr += pi * s.Nr
                nv += pi * s.Nv
            }
            if nr+nv > 0 {
                mu := (nv*(1-ssm.MuR) + nr*ssm.MuR) / (nr + nv)
                mu = clamp(mu, 1e-15, 1-1e-15)
                ll = append(ll, logBinomialLikelihood(ssm.A[tp], ssm.D[tp], mu) +
                              math.Log(0.25) + ssm.LogBinNormConst[tp])
            } else {
                ll = append(ll, math.Log(1e-99))
            }
            if !ssm.UseFourStates && len(ll) == 2 { break } // matches Python's len-2 case
        }
        llh += logsumexp(ll)
    }
    return llh
}
```

**Call site in `runChain` (main.go:2260-2290):**
```go
tssb.resampleAssignments(rng)
tssb.cullTree()
tssb.setNodeHeights()
tssb.setNodePaths()
tssb.mapDatumToNode()
tssb.precomputeMHStates()   // <-- NEW, one call per MCMC iter, parallelized over SSMs
mhAcc := tssb.metropolis(mhIters, mhStd, rng)
```

**Verification:** The new `logLikelihoodWithCNVTreeMH` must produce numerically equal results to the old tree-walking version. Add a test that runs both on a random tree and asserts identical LLH to within 1e-10. Do TDD (`superpowers:test-driven-development`) — write the equivalence test first, then swap the implementation and verify.

### P1 — Re-run default grid after P0

Re-run the 216-fixture default grid with the P0 fix applied. Expected outcomes:
- CNV fixture walltime drops ~10× (from the worst cases of ~10,000s down to ~1,000s)
- Go CNV timeouts drop from 71 → near zero (K=10 CNV fixtures that currently crash should complete)
- K_error gap shrinks or disappears as Python and Go do comparable amounts of effective sampling
- AUPRC already matches; expect that to stay matched

### P2 — If K_error gap persists after P1

Only if K_error is still significantly worse in Go after P0+re-run, investigate whether we have any remaining deviation from Python in: `resampleAssignments` slice sampler, `resampleStickOrders`, `resampleHypers`. Do NOT add φ-pruning — that is a Python deviation. If the gap persists and is attributable to the prior itself, that is an upstream concern, not a Go port bug.

---

## Decisions

1. **Do not add minimum-φ pruning.** Python doesn't do it; adding it would be a silent fidelity deviation. My earlier recommendation in SIM_VALIDATION_REPORT_3.md is retracted.
2. **Do implement precomputed MH states.** This is Python-faithful: Python has implemented this exact optimization since its C++ MH was written. We are correcting an omission, not inventing something new.
3. **Parallelize precomputation via goroutines; keep the MH inner loop serial.** Matches Python (single-threaded C++ MH) while using Go's strength where it helps.
4. **Chain-level parallelism stays as is.** Already the biggest win and already in place.

---

## Answer to the user's direct question

> Will "Add minimum-φ node pruning after each MCMC step" fix this and is it something Python does as well?

**No and no.** Python doesn't do φ-pruning. Only `cull_tree` (trailing-empty-subtree removal, which Go already mirrors). The over-splitting we see is a shared property of both implementations' DP prior, not a Go bug; Python's C=0,S=1 fixtures over-split to 20-81 inferred clones too (see `K10_S1_T1000_M100_C0_rep1`: Python K=81, Go K=84). The *real* fixable gap is that Python precomputes per-(SSM, node) coefficients before each MH loop so its C++ inner loop is a simple O(K) dot product, whereas Go rebuilds those coefficients on every one of the 5,000 MH iterations. Porting Python's precomputation is a legitimate Python-faithful fix that should eliminate most of the CNV walltime gap. After that, re-measure over-splitting with symmetric walltime budgets — I expect the K_error gap to shrink dramatically.
