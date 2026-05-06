# CNV Implementation Gap Analysis: Go Port vs Python Original

_Date: 2026-04-06_
_Reviewed files:_
- _Python: `phylowgs/data.py`, `phylowgs/util2.py`, `phylowgs/params.py`, `phylowgs/evolve.py`, `phylowgs/mh.hpp`, `phylowgs/mh.cpp`_
- _Go: `PhyloWGS_refactor/main.go`_
_Every finding below was verified by reading both source files directly._

---

## Executive Summary

The Go port faithfully implements the CNV likelihood model. The four-scenario `computeNGenomes` tree traversal (the most mathematically complex piece) is correct and produces identical results to the Python original. The MH and slice sampling paths for CNV-bearing SSMs are wired up correctly. Three issues were found: one minor functional discrepancy (edge-case scenario weighting), one missing output annotation (no impact on MCMC), and one block of dead code.

---

## Issue 1 — Genotype Scenario Weighting: Go/Python vs C++ MH `[P2 — Minor]`

### Background

When an SSM has an associated CNV, the likelihood is computed by marginalising over up to four genotype scenarios representing different maternal/paternal phasing assumptions and SSM-before-CNV vs CNV-before-SSM timing. Scenarios where the variant genome count would be zero (`nv = 0`) are excluded as they are unobservable.

### The Discrepancy

The original Python codebase uses **two separate implementations** for this likelihood:

**Python-side (slice sampler / `data.py:55`):**
```python
poss_n_genomes = [x for x in poss_n_genomes if x[1] > 0]   # filter nv=0
for (nr, nv) in poss_n_genomes:
    ll.append(log_binomial_likelihood(...) + log(1.0/len(poss_n_genomes)) + ...)
llh = logsumexp(ll)
```
Uniform prior over the **remaining valid scenarios** (dynamic denominator).

**C++ MH (`mh.hpp:106–160`):**
```cpp
ll[0] = log_binomial_likelihood(...) + log(0.25) + ...;  // fixed 1/4
// ... repeat for ll[1], ll[2], ll[3] ...
llh = logsumexp(ll, 4);  // all 4, nv=0 ones contribute log(1e-99)
```
Fixed `log(0.25)` for all four scenarios. Zero-nv scenarios contribute `log(1e-99)`.

These two formulations give the same result when all four scenarios are valid (both normalise to 1/4 each). They differ only when one or more scenarios yield `nv = 0`, which happens when total copy number is 0 or 1.

**Go port (`main.go:905, 1583`):**
```go
prior := math.Log(1.0 / float64(len(validPairs)))  // dynamic, after filtering nv<=0
lls[i] = logBinomialLikelihood(...) + prior + ...
llh += logsumexp(lls)
```

The Go port consistently uses the Python-side dynamic normalization for **both** slice sampling and MH. In the original, the MH step uses C++ which uses a fixed `1/4`.

### Impact

Only affects CNV SSMs where the total copy number collapses to 0 or 1 (e.g., homozygous deletion CNVs where `major_cn + minor_cn ≤ 1`). In that case:
- Go/Python: effective prior is `1/k` where `k < 4` valid scenarios remain
- C++ original MH: effective prior is `1/4` but `log(1e-99)` terms drag the sum

The difference is numerically small and bounded. The Go behaviour is arguably more statistically principled (the excluded scenarios are genuinely unobservable).

### Recommendation

**No action required** for correctness. The discrepancy only surfaces with very low copy number CNVs and is a pre-existing inconsistency between the Python and C++ sides of the original. Document in code comments.

---

## Issue 2 — `physical_cnvs` Column Not Parsed `[P3 — Output only]`

### Python behaviour (`util2.py:40–51, 75–94`)

The Python parser reads a fifth tab-separated column `physical_cnvs` from `cnv_data.txt`:
```
chrom=1,start=141500000,end=148899999,major_cn=2,minor_cn=1,cell_prev=0.0|0.718
```

This genomic annotation data (chromosome, coordinates, copy numbers per physical segment, cell prevalence) is stored as `cnv_logical_physical_mapping` and embedded in the `trees.zip` output file:

```python
# evolve.py:104
tree_writer.add_extra_file('cnv_logical_physical_mapping.json', json.dumps(cnv_logical_physical_mapping))
```

It is then used only in `pwgsresults/result_generator.py` to annotate the `mutlist` output JSON with physical genomic locations.

### Go behaviour

The Go parser (`main.go:383–457`) reads columns 0–3 (CNV ID, variant reads, total reads, SSM links) but does not read or store the `physical_cnvs` column. No Go output includes physical CNV genomic annotation.

### Impact on MCMC

**None.** The `physical_cnvs` data plays no role in likelihood computation, MCMC sampling, or parameter estimation. It is purely an output annotation.

### Impact on output

The Go port's `best_tree.json` `mut_assignments` section lists CNV IDs but does not annotate them with genomic coordinates, copy numbers per segment, or cell prevalence ranges. If downstream tools that consume the original Python output expect the `cnv_logical_physical_mapping.json` sidecar file, they will not find it.

### Recommendation

**P3 — Low priority.** If Go output needs to be consumed by the same downstream tools as Python output, add parsing and passthrough of the `physical_cnvs` column into the output JSON. This is a one-field read and JSON write; no algorithmic change required.

---

## Issue 3 — Dead Code: `logLikelihood` / `logLikelihoodWithCNV` Fallback `[P3 — Dead code]`

### What exists

`main.go:697–703` defines `(ssm *SSM) logLikelihood(phi []float64)` which dispatches to `(ssm *SSM) logLikelihoodWithCNV(phi []float64)` for CNV-bearing SSMs. The latter (`main.go:928–993`) implements a **simplified two-node model** that does not perform tree traversal:

```go
// Case 1: Variant on major allele
nr1 := (1-p)*2 + p*float64(minorCN)
nv1 := p * float64(majorCN)
```

This simplified model computes likelihood from a single phi value without accounting for the contributions from all TSSB tree nodes, which is incorrect for multi-node trees.

### Is it called?

**No.** Searching all call sites confirms that the MCMC paths exclusively use:
- `logLikelihoodWithCNVTree(ssm, t)` for slice sampling (`resampleAssignments`, `completeDataLogLikelihood`)
- `logLikelihoodWithCNVTreeMH(ssm, t, newState)` for MH (`paramPost`)
- `logLikelihoodNoCNV(params)` for non-CNV SSMs

The `ssm.logLikelihood()` function is unreachable in normal execution. It would only be reachable if a future caller invokes it directly without a TSSB context.

### Risk

If future code inadvertently calls `ssm.logLikelihood(phi)` for a CNV SSM (e.g., in a new hyperparameter sampler or utility function), it would silently produce incorrect likelihoods without error.

### Recommendation

Either remove `logLikelihood` and `logLikelihoodWithCNV`, or add a `panic("logLikelihoodWithCNV: requires tree context; call logLikelihoodWithCNVTree instead")` guard so the error is surfaced immediately.

---

## Verified Correct Implementations

The following CNV components were verified correct against the Python source:

### `computeNGenomes` — Four-Scenario Tree Traversal

**Python (`data.py:65–126`) vs Go (`main.go:761–875`)**

All four cases match exactly:

| Case | Condition | Python | Go |
|------|-----------|--------|----|
| No SSM, no CNV | `!ssmInAncestors && mrCNV == nil` | `nr += pi*2` all scenarios | identical |
| SSM, no CNV | `ssmInAncestors && mrCNV == nil` | `nr += pi, nv += pi` all scenarios | identical |
| No SSM, CNV | `!ssmInAncestors && mrCNV != nil` | `nr += pi*(cp+cm)` all scenarios | identical |
| SSM + CNV, SSM first | `ssmInAncestors && cnvNode in ancestors(ssmNode)` | `nr1 += pi*mr_cnv[1], nv1 += pi*mr_cnv[2]` (swapped for s2) | `nr1 += pi*cm, nv1 += pi*cp` where `cm=parts[1], cp=parts[2]` → same indices |
| SSM + CNV, CNV first | `ssmInAncestors && otherwise` | `nr += pi*max(0,cp+cm-1), nv += pi*min(1,cp+cm)` | identical |

**Return value:** Python returns 4 pairs when `len(cnv)==1 && ssm.node==cnv.node`; Go uses identical condition. ✓

**Note on field naming:** Python stores CNV copy numbers as `(cnv, cp, cm)` where `cp = tok[1]` and `cm = tok[2]` from the `ssms` column. Go stores them as `MaternalCN = parts[1]` and `PaternalCN = parts[2]`. The names are swapped between the two (Python's `cp` ≡ Go's `MaternalCN`; Python's `cm` ≡ Go's `PaternalCN`), but both index `parts[1]` / `mr_cnv[1]` for `nr1` and `parts[2]` / `mr_cnv[2]` for `nv1`, so the computation is identical. The `// FIXED` comments in `main.go:841–845` have a misleading description ("mr_cnv[1]=minor") but the code is correct.

### CNV Likelihood Formula

**Python (`data.py:54`):**
```python
mu = (nr * mu_r + nv*(1-mu_r)) / (nr + nv)
```

**Go (`main.go:914`):**
```go
mu := (nr*ssm.MuR + nv*(1-ssm.MuR)) / total
```

Algebraically identical. ✓

### CNV Datum Likelihood (simple binomial)

**Python (`data.py:60`):**
```python
mu = (1 - phi) * mu_r + phi*mu_v   # mu_r=0.999, mu_v=0.5
```

**Go (`main.go:1540`):**
```go
mu := (1-p)*0.999 + p*0.5
```

Identical. ✓

### `findMostRecentCNV`

**Python (`data.py:128–134`):** Iterates `nd.get_ancestors()[::-1]` (deepest ancestor first) looking for a node matching any of the SSM's CNVs.

**Go (`main.go:731–747`):** Iterates `nd.getAncestors()` from `len-1` downward (deepest first), same logic. ✓

### MH State Handling for CNV SSMs

**Python (`data.py:69`):**
```python
pi = nd.pi1[tp] if new_state else nd.pi[tp]
```

**Go (`main.go:776–780`, `logLikelihoodWithCNVTreeMH:1568`):**
```go
if newState {
    pi = nd.Pi1[tp]
} else {
    pi = nd.Pi[tp]
}
// ... passed via computeNGenomes(ssm, tssb, tp, newState)
```

Both slice sampling and MH correctly toggle between current and proposed state. ✓

### MCMC Integration

Python's per-iteration CNV flow:
1. `tssb.resample_assignments()` — reassigns SSMs and CNVs to nodes
2. `set_node_height / set_path_from_root_to_node / map_datum_to_node` — update tree metadata
3. `metropolis(...)` — MH updates for node phi values

Go's equivalent (`main.go:2155–2170`):
1. `tssb.resampleAssignments(rng)` — ✓
2. `tssb.setNodeHeights() / tssb.setNodePaths() / tssb.mapDatumToNode()` — ✓
3. `tssb.metropolis(...)` — ✓

All three steps include CNV data handling. ✓

---

## Summary Table

| Finding | Severity | Affects MCMC? | Affects Output? | Action |
|---------|----------|--------------|----------------|--------|
| Scenario weighting (Go = Python dynamic, original C++ = fixed 1/4) | P2 | Marginal edge case | No | Document only |
| `physical_cnvs` column not parsed | P3 | No | Yes (annotations) | Parse if downstream tools need it |
| Dead code `logLikelihoodWithCNV` simplified fallback | P3 | No (unreachable) | No | Add panic guard or remove |
| `computeNGenomes` 4-scenario traversal | ✅ Correct | — | — | No action |
| Likelihood formula | ✅ Correct | — | — | No action |
| MH state (pi1 vs pi) | ✅ Correct | — | — | No action |
| findMostRecentCNV ancestor walk | ✅ Correct | — | — | No action |
| MCMC integration (resample → metadata → MH) | ✅ Correct | — | — | No action |
| CNV datum simple binomial | ✅ Correct | — | — | No action |
