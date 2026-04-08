# PhyloWGS Go Port — Gap Analysis Report

_Prepared: 2026-04-06_
_Reviewed files: `phylowgs/` (original Python/C++) and `PhyloWGS_refactor/` (Go port)_
_Post-review note: All findings independently verified against source code before final report._

---

## Executive Summary

The Go port (`PhyloWGS_refactor/main.go`) is a high-quality, substantially faithful reimplementation of the original PhyloWGS MCMC sampler, achieving ~10x wall-clock speedup at production scale. Core algorithms — TSSB stick-breaking, slice sampling, Metropolis-Hastings, hyperparameter resampling, and likelihood computation — all match the original.

The main gaps are in **output completeness** (only the single best tree is saved, not the full posterior) and **operational features** (chain merging, run resumption). One genuine low-severity code quality fix is warranted (`logFactorial`). One dead-code cleanup is needed.

| # | Issue | Severity | Category | Action |
|---|-------|----------|----------|--------|
| 1 | No full posterior sample output | **HIGH** | Missing feature | Implement |
| 2 | `logFactorial` uses O(n) loop vs `math.Lgamma` | **LOW** | Correctness/perf | Fix |
| 3 | Dead `comparePaths` function | **LOW** | Code hygiene | Remove |
| 4 | No chain inclusion-factor filtering | **MEDIUM** | Missing feature | Defer |
| 5 | No run resumption / checkpointing | **MEDIUM** | Missing feature | Defer |
| 6 | Slice sampler iteration hard-cap at 100 | **LOW** | Minor diff | Note only |
| 7 | Likelihood caching absent from slice sampler | **LOW** | Performance only | Note only |
| 8 | Post-processing pipeline not ported | **INFO** | Known scope gap | Document |

---

## Retracted Finding: MH Proposal Alpha

An earlier draft of this report flagged a mismatch between the Dirichlet alpha used for *sampling* (`std·π + 1`) and the alpha used in *correction terms* (`std·π`). After reading `mh.cpp` directly, this is confirmed to be **present in the original C++ code** as well:

```cpp
// mh.cpp: sample_cons_params() — uses +1
alpha[i] = conf.MH_STD*pi[i]+1;
dirichlet_sample(NNODES, alpha, pi_new, rand);

// mh.cpp: mh_loop() correction terms — no +1
theta[i] = conf.MH_STD * pi_new[i];   // forward: q(old|new)
a += gsl_ran_dirichlet_lnpdf(conf.NNODES, theta, pi);
theta[i] = conf.MH_STD * pi[i];        // backward: q(new|old)
a -= gsl_ran_dirichlet_lnpdf(conf.NNODES, theta, pi_new);
```

The Go port reproduces this exactly. This is a pre-existing characteristic of the original algorithm, not a porting bug. It is **not** changed here.

---

## Issue 1 · Missing Feature: Full Posterior Sample Output (HIGH)

### Description

The original PhyloWGS stores **every post-burnin sampled tree** in `trees.zip` as pickled Python objects. This archive is the primary deliverable of the MCMC run — it feeds all downstream tools:

| Downstream tool | What it needs |
|-----------------|---------------|
| `write_results.py` | Full posterior → `*.summ.json.gz`, `*.muts.json.gz`, `*.mutass.zip` |
| `posterior_trees.py` | Full posterior → tree summaries |
| `witness/` | Formatted JSON from `write_results.py` |

The Go port currently saves only:

- `summary.json` — aggregate run statistics (NumChains, BestLLH, time)
- `chain_N_samples.txt` — per-chain LLH trace (Iteration, LLH, NumNodes)
- `best_tree.json` — snapshot of the **single best tree** from the single best chain

The `TreeSample` struct is:
```go
type TreeSample struct {
    Iteration int
    LLH       float64
    NumNodes  int
}
```
It does not persist tree topology or mutation assignments. The `FinalTree *TSSB` in `ChainResult` holds the final tree state per chain, but only the highest-LLH chain's final state is written to `best_tree.json`.

### Impact

Users cannot perform posterior analysis — the essential step between the MCMC run and biological interpretation — without the full sample set.

### Recommended Fix

For each post-burnin sample, serialise the full `summarizePops()` output as a JSON entry, appended line-by-line to a newline-delimited file (`all_trees.ndjson`). This is:
- Go-native and human-readable
- Sufficient for a Go-native or Python downstream pipeline
- Avoids the Python pickle dependency

`TreeSample` should be extended to carry (or reference) the tree snapshot, and `writeResults` should write one `all_trees.ndjson` per chain, or one merged file with a `chain_id` field.

---

## Issue 2 · Bug: `logFactorial` Uses O(n) Loop Instead of `math.Lgamma` (LOW)

### Description

```go
// main.go, line 158
func logFactorial(n int) float64 {
    if n <= 1 { return 0 }
    result := 0.0
    for i := 2; i <= n; i++ {
        result += math.Log(float64(i))
    }
    return result
}
```

This sums `log(2) + log(3) + ... + log(n)` in a loop. For read depths in the tens-of-thousands (realistic for whole-genome sequencing), this is O(d) per data point and accumulates floating-point rounding error.

Go's standard library provides `math.Lgamma(n+1)` which uses the Lanczos approximation — the same numerical method as `scipy.special.gammaln` used by the Python original.

### Impact

`LogBinNormConst` is precomputed once at load time (not in the MCMC hot path), so performance impact is negligible. The accuracy divergence from the Python result is small but non-zero and grows with read depth.

### Fix

```go
func logFactorial(n int) float64 {
    if n <= 1 { return 0 }
    lgamma, _ := math.Lgamma(float64(n + 1))
    return lgamma
}
```

This is a two-line change with no behavioural side effects.

---

## Issue 3 · Code Hygiene: `comparePaths` is Dead Code (LOW)

### Description

`main.go` contains two path comparison functions:

- `pathLT` (line 1363) — used in `resampleAssignments`, correct, matches Python
- `comparePaths` (line 1391) — defined but **never called anywhere**

```go
// comparePaths compares two paths lexicographically
// Returns -1 if path1 < path2, 0 if equal, 1 if path1 > path2
func comparePaths(path1, path2 []int) int { ... }
```

Having two similar-looking but behaviourally distinct comparison functions creates a confusion risk for future contributors.

### Fix

Remove `comparePaths`. The active function is `pathLT`.

---

## Issue 4 · Missing Feature: Chain Inclusion-Factor Filtering (MEDIUM) — Deferred

### Description

Python's `multievolve.py` implements principled chain merging: each chain's log-sum-likelihood is computed, and only chains within a configurable factor of the best are included in the merged posterior.

```python
# Include chain if logSumLH >= chain_inclusion_factor * bestLogSumLH
# Default factor: 1.1
```

The Go port runs all chains and aggregates results without filtering. Chains that get stuck in low-probability regions contribute equally to the output.

### Impact

In practice, when chains diverge significantly (e.g., a local optimum in tree-space), including inferior chains dilutes the posterior estimate. Default `chain_inclusion_factor=1.1` from the original is conservative and would typically include most well-mixing chains.

### Recommendation

Implement after Issue 1 (full posterior output), since the filtering logic acts on the per-chain log-sum-LLH which requires the full LLH trace (already stored in `SampleLLH`). Add a `--chain-inclusion-factor` flag defaulting to `1.1`.

---

## Issue 5 · Missing Feature: Run Resumption / Checkpointing (MEDIUM) — Deferred

### Description

The original `evolve.py` saves `state.last.pickle` periodically. On restart it detects this file and resumes from the exact saved state (RNG state, iteration count, full tree).

The Go port has no equivalent. A killed job must start from scratch.

### Impact

Critical in Slurm/HPC environments with wall-time limits. Production runs (3500+ iterations × 6 chains) may exceed typical queue limits.

### Recommendation

After each sample (or every N=100 iterations), serialise each chain's state to `chain_N_state.json`. Fields: RNG seed + iteration state, iteration counter, tree topology, all Pi/Params vectors, hyperparameters (dp_alpha, dp_gamma, alpha_decay), mhStd.

---

## Issue 6 · Minor: Slice Sampler Hard-Capped at 100 Iterations (LOW)

### Description

```go
for iter := 0; iter < 100; iter++ { ... }   // Go
```
```python
while True: ...  # Python, terminates on epsilon-shrink
```

In practice slice sampling converges in <20 steps for typical PhyloWGS inputs. The cap is safe. Python terminates on `abs(max_u - min_u) < machine_epsilon` which also rarely requires more than ~20 steps.

This is acceptable as-is. Consider logging a warning if the cap is hit, to detect degenerate inputs:

```go
if iter == 99 {
    log.Printf("Warning: slice sampler hit 100-iteration cap for datum %d", n)
}
```

---

## Issue 7 · Minor: No Likelihood Cache in Slice Sampler (LOW, Performance Only)

### Description

Python caches node log-likelihoods within each datum's slice loop (`llhmap[node] = llh`). Go recomputes every time.

For non-CNV SSMs (O(1) likelihood evaluation), this is negligible. For CNV SSMs, `logLikelihoodWithCNVTree` traverses the entire tree on every proposal — potentially O(nodes × iterations) per datum.

Adding a `map[*Node]float64` cache per datum (cleared before each datum) would reduce redundant CNV likelihood evaluations. This is a pure performance optimisation with no correctness implications.

---

## Issue 8 · Known Scope Gap: Post-Processing Pipeline Not Ported (INFO)

### Description

The following Python components have no Go equivalent and are not expected to be ported:

| Component | Purpose |
|-----------|---------|
| `write_results.py` | Loads posterior samples, filters clones, writes output JSON |
| `pwgsresults/result_generator.py` | Population summary extraction |
| `pwgsresults/result_munger.py` | Subclone filtering, multiprimary removal |
| `pwgsresults/json_writer.py` | Output formatting for web viewer |
| `posterior_trees.py` | Posterior tree distribution summary |
| `cc.py` | Correlation clustering for mutation grouping |
| `witness/` | Web-based tree visualisation |

These tools remain usable as Python post-processing steps once the Go port produces compatible output (see Issue 1).

### Recommendation

Document the intended two-stage workflow in `README.md`:
1. Run the Go binary → produces `all_trees.ndjson` (pending Issue 1)
2. Run `write_results.py` → produces final JSON files for downstream / visualisation

---

## Confirmed Correct Components

The following were reviewed and found to match the original faithfully:

| Component | Status |
|-----------|--------|
| MCMC iteration order | ✓ Match |
| Node spawn / kill (Pi conservation) | ✓ Match |
| `findOrCreateNode` path construction (root index included at depth=0) | ✓ Match |
| `pathLT` comparison (string %03d encoding, reversed s2>s1) | ✓ Match |
| MH proposal: sampling uses `std·π + 1`, correction uses `std·π` | ✓ Faithful reproduction of C++ original |
| `resampleSticks` (reverse-order Beta posteriors, root forced 0.999999) | ✓ Match |
| `resampleStickOrders` (weighted ordering, Pi return on kill) | ✓ Match |
| `cullTree` (trailing-empty-child removal, Pi return) | ✓ Match |
| `resampleHypers` (dp_alpha, alpha_decay, dp_gamma with correct bounds) | ✓ Match |
| Non-CNV likelihood: `mu = (1-φ)·mu_r + φ·mu_v` | ✓ Match |
| CNV likelihood (four-case maternal/paternal timing model) | ✓ Match |
| Data loading (SSM + CNV TSV, LogBinNormConst precomputation) | ✓ Match |
| Hyperparameter initialisation (dp_alpha=25, dp_gamma=1, alpha_decay=0.25) | ✓ Match |

---

## Summary of Recommended Actions

| Priority | Issue | Action |
|----------|-------|--------|
| **P1** | No full posterior sample output (#1) | Implement `all_trees.ndjson` per chain |
| **P2** | `logFactorial` loop vs `math.Lgamma` (#2) | Two-line fix; improves accuracy |
| **P2** | Dead `comparePaths` function (#3) | Remove |
| **P3** | Chain inclusion-factor filtering (#4) | Design + implement after #1 |
| **P3** | Run resumption / checkpointing (#5) | Design + implement; HPC requirement |
| **P4** | Slice sampler cap warning (#6) | Optional log warning |
| **P4** | Slice sampler LLH cache (#7) | Optional performance improvement |
| **P4** | Document two-stage workflow (#8) | README update |
