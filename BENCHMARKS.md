# PhyloWGS Go Rewrite — Benchmark Results

**Machine:** AMD Ryzen 7 9800X3D (8C/16T) | RTX 5080 16GB VRAM | WSL2 (Windows 11)
**Date:** 2026-03-11/12

---

## Implementations

| Label | Description | Branch |
|---|---|---|
| `original-python` | Original PhyloWGS Python 2/3 | — |
| `optimized-python` | Python 3 refactor (this repo, `main`) | main |
| `go-cpu` | Pure Go rewrite, CPU only | go/main |
| `go-gpu` | Go + CUDA kernel (RTX 5080, sm_86 fallback) | go/feature/cuda-likelihood |
| `go-cpu-opt` | Go CPU with MH caching + perf optimizations | go/feature/parallel-traversal |

---

## Small-Scale Benchmark (11 SSMs, 1 CNV)

**Data:** Synthetic test data (`phylowgs-optimized/ssm_data.txt`)
**Config:** B=1000, s=2500, j=4 (production-equivalent)

| Implementation | Wall time | Samp/s/chain | Peak RAM | LLH |
|---|---|---|---|---|
| original-python | 5m 49s | ~0.07 | ~80 MB | -1748041.14 |
| go-cpu | 0m 38s | ~1.1 | ~11 MB | -1748041.14 |
| go-gpu | 0m 40s | ~1.0 | ~45 MB | -1748041.14 |

**Speedup:** Go CPU **9.2×** faster than original Python. GPU slower than CPU at this scale (kernel launch overhead dominates with only 11 SSMs).

---

## Large-Scale Benchmark (13,862 SSMs, 309 CNVs)

**Data:** TCGA SKCM melanoma (real patient data, 1 timepoint)
**Config:** B=5, s=10, j=2 (smoke config — original Python too slow for full run)

### Pre-optimization

| Implementation | Wall time | Samp/s/chain | Peak RAM | Notes |
|---|---|---|---|---|
| original-python | ~40s/iter | ~0.04 | ~138 MB | Extrapolated — too slow to benchmark |
| optimized-python | 96s (B=5,s=10,j=1) | 0.10 | 138 MB | Single chain |
| go-cpu | 35s (B=5,s=10,j=2) | 0.28 | **27 MB** | |
| go-gpu | 39s (B=5,s=10,j=2) | 0.28 | 143 MB | No GPU benefit at this scale |

**Speedup vs original-python:** Go CPU **~7×**. Optimized Python **~2.5×**.

### Post-optimization (`go/feature/parallel-traversal`)

**Config:** B=10, s=20 (comparable)

| | go-cpu (before) | go-cpu-opt (after) | Speedup |
|---|---|---|---|
| j=1 | 0.26 samp/s/chain | **0.56 samp/s/chain** | **2.15×** |
| j=4 wall time | 76.6s | **42.4s** | **1.81×** |

**Speedup vs original-python (post-opt):** Go CPU **~14×**

---

## GPU Scaling Analysis

GPU acceleration via CUDA kernel (`go/feature/cuda-likelihood`):

| Dataset size | Go CPU | Go GPU | GPU benefit |
|---|---|---|---|
| 11 SSMs | 9.2× vs Python | 8.8× vs Python | ❌ Slower than CPU |
| 13,862 SSMs | ~7× vs Python | ~7× vs Python | ❌ No gain |

**Bottleneck at large scale:** Tree traversal (MCMC moves, node creation, stick sampling) — CPU-bound, not parallelizable on GPU. Likelihood kernel is no longer the bottleneck once MH caching is applied.

**When GPU would help:** Dense multi-timepoint data (5+ timepoints) × many SSMs, where likelihood computation outweighs tree ops. Not typical for single-sample PhyloWGS runs.

---

## Profiling Results (pprof, 13,862 SSMs pre-optimization)

```
42%   math.archLog              ← log() calls in binomial likelihood
31%   logBinomialLikelihood
16%   logLikelihoodNoCNV        ← per-SSM likelihood computation
 6%   paramPost / metropolis    ← MH phi sampling
 0.6% resampleAssignments       ← NOT the bottleneck (assumed wrong)
```

**Key finding:** `resampleAssignments` was only 0.6% of CPU. The real issue was `metropolis()` recomputing `oldLLH` on every proposal (~5,000/iteration), despite it only changing on MH accept (~3–8% acceptance rate). ~4,700 redundant full-tree evaluations per iteration.

---

## Optimization Summary

| Optimization | Rationale | Impact |
|---|---|---|
| Cache `oldLLH` in MH loop | Only recompute on accept; frozen between rejected proposals | **Main win** |
| Cache mixture weights | `getMixture()` is O(N_nodes), called in inner SSM loop | Medium |
| Pre-compute ancestry sets | `isAncestorOf()` was O(depth) per call → O(1) map lookup | Medium |
| Parallel `resampleAssignments` | Tried — goroutine overhead (20% `runtime.futex`) outweighed gains | Not shipped |

**Note:** The `oldLLH` caching bug is also present in the original Python PhyloWGS (same structural pattern in `tssb.py`). The fix is equally applicable there.

---

## Memory Scaling

| Implementation | 11 SSMs | 13,862 SSMs | Scaling |
|---|---|---|---|
| optimized-python | ~80 MB | 138 MB | ~1.7× |
| go-cpu | 11 MB | **27 MB** | ~2.5× |
| go-gpu | 45 MB | 143 MB | ~3.2× (CUDA context overhead) |

Go CPU is the most memory-efficient implementation across all scales.

---

## Algorithm Fidelity Score (vs original PhyloWGS)

Evaluated by automated code review against original `tssb.py`, `data.py`, `evolve.py`:

| Version | Score | Notes |
|---|---|---|
| optimized-python | 100/100 | Mathematically equivalent to original |
| go-cpu (initial) | 79.5/100 | CNV likelihood wrong, missing MCMC moves |
| go-cpu (post-fix) | **89.5/100** | CNV tree traversal, `resampleStickOrders()`, dynamic nodes implemented |

Remaining gap (10.5 pts): `munge_results.py`/`posterior_trees.py` require pickle format not supported by Go output; minor MCMC move schedule differences.

---

## HPC Recommendation (SLURM)

For MSK Juno/Iris HPC cluster:

```bash
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=4        # one task per chain (-j 4)
#SBATCH --mem=2G          # go-cpu uses ~27MB at 13k SSMs; headroom for larger datasets

./phylowgs-go-cpu --no-gpu -j 4 -B 1000 -s 2500 ssm_data.txt cnv_data.txt
```

Chain-level parallelism (`-j N`) scales cleanly. Intra-chain goroutine parallelism does not benefit at typical MSK-IMPACT scale due to goroutine overhead.

For very large datasets (50k+ SSMs, 5+ timepoints): GPU acceleration expected to pay off. Re-benchmark with `go-gpu` at that scale.
