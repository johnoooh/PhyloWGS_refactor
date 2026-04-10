# RESOLVED: K-Inflation Was a Comparison Error, Not a Bug

**Date:** 2026-04-09  
**Author:** Investigation session  
**Status:** RESOLVED — No code fix needed  

## Executive Summary

The reported K-inflation in the Go port was caused by comparing **raw MCMC node
counts** (Go) against **post-processed population counts** (Python). When
comparing like-for-like, the Go port matches the original Python 2 code almost
perfectly.

## The Comparison Error

### What was reported

> "Go produces K=67 while Python produces K=10 for K3_S1_T200_M150_C0_rep2"

### What actually happened

| Metric | Go port | Original Python 2 |
|--------|---------|-------------------|
| **Raw total nodes (K)** | median 288 | median 348 |
| **Non-empty populations (K_pop)** | **median 66** | **median 70** |
| **Post-processed populations** | N/A (not implemented) | **median 10** |

The "K=10" Python number came from `tree_summaries.json.gz`, which applies two
rounds of post-processing:

1. **`remove_empty_nodes()`** — removes all nodes without SSM/CNV data,
   collapsing empty internal nodes (reduces K from ~350 to ~70)
2. **`munge_results.py`** → `remove_superclones()` + additional merging
   (further reduces to ~10-15)

The Go port reported `num_populations` as the count of non-empty nodes from the
live MCMC tree — equivalent to Python's raw K_pop (~70), not the post-processed
count (~10).

## Evidence

### Original Python 2 HPC run (K3_S1_T200_M150_C0_rep2, seed 1173227769)

Raw TSSB state extracted from `trees.zip` pickles:

```
Burnin trajectory (raw K from pickled TSSB):
  iter -500:  K=36,   K_pop=18
  iter -400:  K=172,  K_pop=52
  iter -250:  K=346,  K_pop=70
  iter -150:  K=383,  K_pop=71
  iter -50:   K=353,  K_pop=70

Sampling trajectory:
  iter 0:     K=438,  K_pop=71
  iter 500:   K=537,  K_pop=71
  iter 999:   K=159,  K_pop=61

Final state (state.last.pickle):
  K=159, K_pop=61, mh_std=6400, alpha_decay=0.64
```

**The original Python 2 code produces K=150-537 in the raw MCMC state.** The
K=10 in tree_summaries is purely a post-processing artifact.

### Go port (same fixture, 4 chains)

```
Chain 0: num_populations median=66, mean=64.2
Chain 1: num_populations median=67, mean=65.8
Chain 2: num_populations median=66, mean=62.3
Chain 3: num_populations median=67, mean=65.9
```

### Apples-to-apples comparison

```
Go non-empty populations:  median=66
Python raw K_pop:          median=70
Ratio (Go/Python):         0.94  ← essentially identical
```

## Earlier Analysis Was Incorrect

The earlier "root cause analysis" document (`2026-04-09-dirichlet-perturbation-analysis.md`)
identified the missing `+0.0001` Dirichlet perturbation as the cause of
M-dependent over-splitting. While the perturbation difference between Go and
C++ is real, it was evaluated against the wrong baseline:

- **Claimed:** Go K=67, Python K=10, therefore Go has a bug
- **Actual:** Go K=67, Python K=70, therefore Go matches Python

The +0.0001 perturbation in `util.cpp:dirichlet_sample` is present in the C++
`mh.o` binary used by ALL implementations (Python 2, Python 3, and Go via C++
proxy). Both Go and Python experience the same MH dynamics because they both
call the same `mh.o` binary.

## Hypotheses Tested This Session

| # | Hypothesis | Result |
|---|-----------|--------|
| 26 | numpy PRNG differs between 1.x and 2.x | **FALSE** — identical sequences for all distributions (rand, beta, gamma) |
| 27 | Python 3 port introduced K inflation | **FALSE** — Python 2 has same raw K |
| 28 | K=10 in HPC results is raw MCMC state | **FALSE** — it's post-processed by remove_empty_nodes + munge_results |
| 29 | Go K diverges from Python K | **FALSE** — Go median=66, Python median=70, ratio=0.94 |

## What Needs to Happen

### No code changes needed for K inflation

The Go port correctly reproduces the original Python 2 MCMC behavior. The raw
node counts match. The "K inflation" was never a real discrepancy.

### Go port needs post-processing

To produce final output comparable to Python's `tree_summaries.json.gz`, the
Go port needs to implement:

1. `remove_empty_nodes()` equivalent — prune empty leaves and bypass empty
   internal nodes
2. Optional: `remove_superclones()` equivalent — merge parent-child pairs
   meeting specific criteria

These are **output post-processing** steps, not MCMC algorithm changes.

### The +0.0001 perturbation

The Go `dirichletSample()` is missing the `+0.0001` additive perturbation
that C++ has. However, since the Go port uses C++'s `mh.o` for the actual MH
step (which includes the perturbation internally), this only matters if Go's
Dirichlet sampling is used for something other than the MH proposal — which
it currently is not. This should still be added for correctness if Go ever
uses its own MH implementation.

## File References

- Go tree output: `PhyloWGS_refactor/main.go:3199-3214`
- Python remove_empty_nodes: `phylowgs/util2.py:128-155`
- Python result_generator: `phylowgs/pwgsresults/result_generator.py:37-41`
- Python result_munger: `phylowgs/pwgsresults/result_munger.py:23-78`
- HPC comparison script: `simulation_validation4/hpc_head_to_head.py:37`
- HPC results: `simulation_validation4/results/{go-cpu,original-python}/`
