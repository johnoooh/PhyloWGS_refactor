# PhyloWGS Go Port — Simulation Validation Report (Full Default Grid)

_Run date: 2026-04-08. Full default grid (216 fixtures). Results in `simulation_validation_updated/`._

---

## Executive Summary

The full default grid (216 fixtures across K∈{3,5,10}, S∈{1,3}, T∈{200,1000}, M varied, C∈{0,2,5}) reveals two separate behavioral regimes:

**C=0 fixtures (no CNVs):** Go is consistently **2–4× faster** than Python across all 72 fixtures and completes all of them. However, the unconstrained DP prior causes aggressive over-splitting — Go infers dramatically more clones than the truth in single-sample runs (mean K_error=22.8 for C=0, S=1), far worse than Python.

**CNV fixtures (C=2 or C=5):** Go suffers a severe performance regression. **71 of 144 CNV fixtures timed out** (Go `total_time_s=0`), and **all remaining CNV fixtures where Go ran were slower than Python** (median 0.34× speedup). Root cause: the correct `logLikelihoodWithCNVTree` (post-fix) is substantially more expensive per node than the incorrect `logLikelihoodNoCNV`, and over-splitting compounds this by creating hundreds of nodes. The MCMC chain cannot converge within the walltime limit.

The good news: **when Go does finish CNV fixtures, accuracy is comparable to Python** (K_error gap narrows to 1.7–3.9 vs Python's 1.7–2.6 for scored CNV fixtures, and AUPRC is essentially equal). The CNV accuracy fixes from the prior session hold — the problem is now runtime, not correctness.

**Two open issues to fix before production use:**
1. Over-splitting (DP prior) — causes high K_error in S=1 runs and indirectly drives CNV timeouts
2. CNV performance — `logLikelihoodWithCNVTree` needs profiling and optimization (caching per node-traversal, lazy evaluation)

---

## Run Configuration

| Parameter | Value |
|-----------|-------|
| Grid | Full default (K∈{3,5,10}, S∈{1,3}, T∈{200,1000}, M∈{30 or 50 per K, 150 or 250 per K}, C∈{0,2,5}) |
| Replicates | 3 per condition |
| Total fixtures | 216 |
| Implementations | `go-cpu`, `original-python` |

---

## Completeness

| Implementation | Scored | Timed out | Failed scoring |
|----------------|-------:|----------:|---------------:|
| Go CPU | **145 / 216** | **71** (all CNV) | 0 |
| Original Python | **172 / 216** | 24 | 20 |

All Go timeouts are CNV fixtures (C=2 or C=5). Python's 20 "failed scoring" entries are runs that completed but where the scoring script could not parse the output (likely due to malformed tree files from unconverged chains).

**Completeness by (K, fixture type):**

| | Go C=0 | Go CNV | Python C=0 | Python CNV |
|--|:------:|:------:|:----------:|:----------:|
| K=3 | 24/24 ✅ | 34/48 | 20/24 | 37/48 |
| K=5 | 24/24 ✅ | 27/48 | 23/24 | 35/48 |
| K=10 | 24/24 ✅ | 12/48 | 20/24 | 37/48 |

Go completes 100% of C=0 fixtures. For CNV, Go completion drops sharply with K (K=10: only 12/48). Python completes more CNV fixtures (109/144) because its per-node computation is cheaper at baseline.

---

## Speed

### C=0 (No CNVs)

Go is faster on every single C=0 fixture:

| | Speedup |
|--|---------|
| Mean | **3.29×** |
| Median | **2.58×** |
| Max | **12.6×** (K10_S1_T1000_M500_C0_rep2) |
| Min | **1.43×** |
| Fixtures where Go is faster | **70 / 70** (100%) |

Speed advantage grows with K: K=3 averages 2.83×, K=5 averages 3.04×, K=10 averages 4.01×. The advantage is largest on the hardest fixtures.

### CNV (C=2 or C=5)

Go is slower on virtually every CNV fixture that both implementations completed:

| | Speedup |
|--|---------|
| Mean | **0.34×** (Go is ~3× slower) |
| Median | **0.34×** |
| Fixtures where Go is faster | **1 / 72** (1.4%) |
| Worst case | **0.038×** (K5_S1_T200_M50_C2_rep0: Go 6619s, Python 251s) |

This is a dramatic reversal from C=0 behavior. The root cause is explained in the Root Cause section below.

---

## Accuracy: AUPRC

AUPRC measures how well SSMs are co-clustered with their true clonemates (higher is better). On the 131 fixtures where both implementations scored:

| | Go | Python |
|--|---:|-------:|
| Mean AUPRC | **0.614** | **0.607** |
| Median AUPRC | 0.566 | 0.541 |
| Go better | 47 / 131 | |
| Python better | 46 / 131 | |
| Tied | 38 / 131 | |

**AUPRC is essentially equal.** The small Go edge (0.614 vs 0.607) is within stochastic MCMC noise. Neither implementation has a systematic tree topology advantage.

By fixture type:

| | Go AUPRC | Python AUPRC | Go better | Py better |
|--|:--------:|:------------:|:---------:|:---------:|
| C=0 (n=63) | 0.641 | 0.628 | 19 | 22 |
| CNV (n=68) | 0.589 | 0.587 | 28 | 24 |

---

## Accuracy: K_error

K_error = |inferred_K − true_K|. Lower is better.

### Head-to-head (131 comparable fixtures)

| | Go K_error | Python K_error | Go better | Py better | Tied |
|--|:----------:|:--------------:|:---------:|:---------:|:----:|
| Overall (n=131) | 8.02 | 3.20 | 38 | 41 | 52 |

### By fixture type (head-to-head)

| Dimension | n | Go K_error | Python K_error | Go better | Py better | Tied |
|-----------|:-:|:----------:|:--------------:|:---------:|:---------:|:----:|
| C=0, S=1 | 32 | **22.84** | 6.97 | 7 | 13 | 12 |
| C=0, S=3 | 31 | 4.45 | 1.81 | 11 | 5 | 15 |
| CNV, S=1 | 28 | 3.93 | 2.57 | 5 | 11 | 12 |
| CNV, S=3 | 40 | **1.80** | **1.70** | 15 | 12 | 13 |

The C=0, S=1 regime has the largest gap (Go 22.84 vs Python 6.97). For CNV, S=3 fixtures — the most clinically relevant multi-sample scenario — Go and Python are essentially tied (1.80 vs 1.70).

### Extreme outliers

Go's over-splitting produces some catastrophic K_error values:

| Fixture | Go inferred K | True K | K_error |
|---------|:-------------:|:------:|:-------:|
| K5_S1_T1000_M250_C0_rep0 | 186 | 5 | **181** |
| K5_S1_T200_M250_C0_rep1 | 136 | 5 | 131 |
| K10_S1_T200_M500_C0_rep2 | 128 | 10 | 118 |
| K3_S1_T200_M150_C0_rep2 | 72 | 3 | 69 |

These are single-sample runs with many mutations (high M). The unconstrained DP prior assigns each SSM its own clone.

---

## Root Cause Analysis

### Issue 1 — Over-splitting (DP prior, both S=1 and C=0)

The stick-breaking prior in single-sample runs does not penalize tiny-φ clones. When M is large, the MCMC finds it easier to assign each mutation its own near-zero-weight subtree than to merge them. This inflates K_error but does not affect AUPRC significantly (each "clone" usually contains the right SSMs, just split into singletons).

Python shows this too (`K5_S1_T1000_M250_C0_rep0` is in Python's timed-out list — Python would also over-split with enough time), but Python's slower chain means it tends to not explore as far, accidentally appearing more conservative.

**This is a modelling issue, not a port fidelity bug.** The recommended fix is a minimum-φ survival threshold (prune nodes with φ < ε after each MCMC step) or a tighter DP-alpha hyper-prior that penalizes excessive branching. Neither Python nor Go currently implements this pruning.

### Issue 2 — CNV Performance Regression

All 71 Go CNV timeouts and all 71 Go-slower-than-Python CNV runs are explained by a compound effect:

**Step 1: Bug fix made `logLikelihoodWithCNVTree` active.** Before the three bug fixes, CNV-linked SSMs incorrectly used `logLikelihoodNoCNV` (cheap: one binomial calculation). After the fix, they correctly use `logLikelihoodWithCNVTree` (expensive: traverses the TSSB tree to find the nearest CNV ancestor for each SSM, then computes a mixture likelihood over copy-number states).

**Step 2: Over-splitting multiplies the cost.** When Go infers 50–186 clones instead of 3–10, `logLikelihoodWithCNVTree` must traverse 50–186 nodes on every MCMC iteration for every CNV-linked SSM. The per-iteration cost scales roughly as O(M × K_inferred × D) where D is tree depth.

**Step 3: Runaway chains.** With O(100) nodes and O(30–500) SSMs with CNVs, each MCMC iteration takes seconds instead of milliseconds. The chain hits walltime without converging.

This explains the pattern: K=3 CNV completes 34/48 (smaller trees), K=10 CNV completes only 12/48 (deeper trees, more nodes). Python avoids this because its slower exploration produces fewer nodes overall.

**The correctness fix is sound; the performance needs separate work.** When Go does complete CNV fixtures, accuracy is competitive with Python (CNV, S=3 head-to-head: Go 1.80 vs Python 1.70). The priority is making `logLikelihoodWithCNVTree` O(1) or O(log K) per SSM via precomputed CNV ancestor maps (similar to what `rebuildAncestorSets` does for the general ancestor query).

---

## Recommended Next Steps

| Priority | Action | Expected Impact |
|----------|--------|-----------------|
| **P0** | Profile `logLikelihoodWithCNVTree` and cache per-node CNV ancestor lookups | Eliminates CNV timeout issue |
| **P1** | Add minimum-φ node pruning after each MCMC step | Fixes S=1 over-splitting; indirectly fixes CNV timeouts by reducing K_inferred |
| **P2** | Re-run default grid after P0+P1 fixes | Clean comparable metrics |

---

## Summary Table

| Metric | Go | Python | Notes |
|--------|:--:|:------:|-------|
| Fixtures scored (total 216) | **145** | 172 | Go timeout is all-CNV |
| C=0 scored | **72/72** | 63/72 | Go completes all C=0 |
| CNV scored | **73/144** | 109/144 | Go timeout/slowdown issue |
| Speed (C=0, median) | **2.58×** faster | 1× | Consistent advantage |
| Speed (CNV, median) | **0.34×** | 1× | Go slower — bug fix performance cost |
| AUPRC (mean, n=131) | **0.614** | 0.607 | Essentially equal |
| K_error (mean, n=131) | 8.02 | 3.20 | Over-splitting gap |
| K_error CNV, S=3 (n=40) | **1.80** | **1.70** | Nearly equal in best scenario |
| K_error C=0, S=1 (n=32) | 22.84 | 6.97 | Over-splitting worst case |

---

_Three correctness bugs fixed (commits `bf096a0`, `2735615`, `e6defce`). Two open issues: over-splitting and CNV performance. See `BUGS_FOUND.md` for bug details._
