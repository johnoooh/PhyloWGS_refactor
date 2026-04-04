# PhyloWGS Go Port — Code Review vs Original Python 2

_Reviewed 2026-04-04 by Sevry (manual line-by-line comparison)_

## Executive Summary

The Go port is a **substantially faithful reimplementation** of the original Python MCMC. After the overnight fix session (commits f925127, 0c58226), most critical bugs were resolved. One **medium-severity bug** remains in the slice sampler path comparison, plus one **unknown-severity divergence** in the MH proposal (Go reimplements C++ `mh.o` from scratch). Several minor differences exist but are acceptable.

**Bugs found:** 1 medium, 1 unknown (MH), 3 low
**Recommendation:** Fix the path bug, then validate.

---

## Component-by-Component Analysis

### 1. MCMC Iteration Loop Order
**MATCH** ✓

Python (`evolve.py:do_mcmc`):
```
resample_assignments → cull_tree → [set heights/paths/datum map] → metropolis → 
[adapt MH std] → resample_sticks → resample_stick_orders → resample_hypers → LLH
```

Go (`runChain`):
```
resampleAssignments → cullTree → [setNodeHeights/setNodePaths/mapDatumToNode] → 
metropolis → [adapt] → resampleSticks → resampleStickOrders → resampleHypers → LLH
```

Identical ordering. Python also assigns sequential node IDs between cull_tree and metropolis; Go does this inside metropolis(). Equivalent.

### 2. Node Initialization (`alleles.__init__` vs `spawnChild`)
**MATCH** ✓ (after overnight fix)

Python:
```python
self.pi = rand(1) * parent.pi
parent.pi = parent.pi - self.pi
self.params = self.pi
```

Go:
```go
frac := rng.Float64()
childNode.Pi[i] = frac * parent.Node.Pi[i]
parent.Node.Pi[i] -= childNode.Pi[i]
childNode.Params[i] = childNode.Pi[i]
```

Pi conservation is correct in both.

### 3. Node Kill (`alleles.kill` vs `killNode`)
**MATCH** ✓ (after overnight fix)

Python: `self._parent.pi += self.pi`
Go: `parent.Pi[i] += child.Pi[i]`

### 4. Slice Sampler (`resample_assignments` vs `resampleAssignments`)

#### 4a. Path construction — **BUG (MEDIUM)**

**`currentIndices`** is built by walking from Root to the assigned node:
```go
ancestors := oldNode.getAncestors()
current := t.Root
for _, anc := range ancestors[1:] { // skip root
    for i, child := range current.Children {
        if child.Node == anc { currentIndices = append(currentIndices, i) ... }
    }
}
```
This **includes** the root-level index (e.g., `[0, 2, 1]` where 0 = root→child0).

**`newPath`** from `findOrCreateNode` at depth=0:
```go
// At root (depth 0):
node, path := descend(root.Children[index], u, depth+1)
return node, path  // Does NOT prepend index=0
```
This **omits** the root-level index (e.g., `[2, 1]`).

**Python:** Both `indices` and `find_node`'s returned path include the root-level index — `path.insert(0, index)` is outside the if/else and always executes.

**Impact:** The `pathLT` comparison compares paths of different depths. Since root always has one child (index=0), `currentIndices` = `[0, ...]` while `newPath` = `[...]`. String comparison of `"000..."` vs `"..."` causes systematic bias: `pathLT` almost always returns -1 (path2 < path1), so the slice sampler always sets `minU = newU` and never `maxU = newU`. The upper bracket rarely shrinks.

**Fix:** Prepend root index in `findOrCreateNode` at depth=0:
```go
node, path := descend(root.Children[index], u, depth+1)
return node, append([]int{index}, path...)  // Prepend root index
```

#### 4b. Iteration limit — **LOW**
Python: `while True` (no limit). Go: `for iter := 0; iter < 100; iter++`. Slice sampler converges much faster than 100 iterations in practice.

#### 4c. Likelihood caching — **LOW (performance only)**
Python caches `llhmap[node] = llh` and reuses. Go recomputes each time. Correct but slower.

#### 4d. Data movement during evaluation — **MATCH** ✓
Python moves data before evaluation, restores on rejection. Go only moves on acceptance for non-CNV (correct since non-CNV LLH depends only on phi, not data assignment). Go does temporary move for CNV SSMs (correct).

### 5. `find_node` vs `findOrCreateNode`
**MINOR DIFF** — Go has `maxChildCreations = 20` safeguard. Python has no limit. Acceptable safety measure.

### 6. `boundBeta` — **LOW**
Python: `(1-eps)*(beta(a,b)-0.5) + 0.5` where eps ≈ 2.2e-16. Effectively beta(a,b) with infinitesimal bounding.
Go: `gamma(a)/(gamma(a)+gamma(b))` with clamp at 1e-6. The 1e-6 vs ~1e-16 boundary difference is negligible for MCMC.

### 7. `resample_sticks` vs `resampleSticks`
**MATCH** ✓

Both iterate children in reverse, compute post_alpha = 1 + child_data, post_beta = dp_gamma + data_down. Both force root main to 1e-30 and root sticks to 0.999999. Minor: Go avoids a wasted RNG call at depth=0 (goes directly to 1e-30 instead of sampling then overriding).

### 8. `resample_stick_orders` vs `resampleStickOrders`
**MATCH** ✓

Both sample children weighted by sticks, create new children if u falls beyond existing sticks, kill unrepresented children (with pi return in Go), reset sticks to 0, then immediately call resample_sticks. Go has `maxCreations = 10` safeguard.

### 9. `cull_tree` vs `cullTree`
**MATCH** ✓ (after overnight fix)

Python: `trim_zeros(counts, 'b')` to find last non-zero, then `child['node'].kill()` for trailing empty children.
Go: Reverse scan for last non-zero, then `killNode()` for trailing empty children.

### 10. `resample_hypers` vs `resampleHypers`
**MATCH** ✓ (after overnight fix added dp_gamma)

All three hyperparameters (dp_alpha, alpha_decay, dp_gamma) resampled via slice sampling with matching bounds:
- dp_alpha: [1.0, 50.0]
- alpha_decay: [0.05, 0.80]  
- dp_gamma: [1.0, 10.0]

Go uses 100-iteration limit; Python uses while True with exception on zero-shrink. Acceptable.

### 11. Metropolis-Hastings — **UNKNOWN SEVERITY**

**Critical structural difference:** Python calls external C++ binary `mh.o` via subprocess for MH updates. Go reimplements MH in pure Go.

Python's `params.py:metropolis()`:
```python
sp.check_call(['mh.o', MH_ITR, MH_STD, ...])
update_tree_params(tssb, FNAME_C_PARAMS)
```

Without the C++ source (`mh.cpp`), we cannot verify:
1. Whether the Go Dirichlet proposal matches the C++ proposal
2. Whether `alpha[i] = std*pi[i] + 1` vs `alpha[i] = std*pi[i]` is correct
3. The exact params-from-pi computation in C++

Go's implementation is a **reasonable** tree-structured Dirichlet MH sampler with:
- Joint pi proposal via Dirichlet(std*pi + 1)
- Leaf-to-root params accumulation: phi(n) = pi(n) + Σ phi(children)
- MH correction terms for proposal asymmetry
- Acceptance ratio: newLLH - oldLLH + correction

This should produce valid MCMC samples but mixing behavior may differ from C++ implementation.

### 12. Likelihood: Non-CNV SSMs
**MATCH** ✓

Python: `mu = (1-phi)*mu_r + phi*mu_v`
Go: `mu := (1-p)*ssm.MuR + p*ssm.MuV`

### 13. Likelihood: CNV SSMs (`compute_n_genomes`)
**MATCH** ✓ (after overnight fix of maternal/paternal swap)

All four cases match:
1. No variant, no CNV: nr += pi*2
2. Has variant, no CNV: nr += pi, nv += pi
3. No variant, has CNV: nr += pi*(cp+cm)
4. Has variant AND CNV: Complex timing cases with 2 or 4 (nr,nv) pairs

### 14. `complete_data_log_likelihood`
**MATCH** ✓

Python: `Σ [count * log(weight) + Σ datum_llh]` per node
Go: Same structure, plus CNV datum likelihoods (added in overnight fix)

### 15. `remove_empty_nodes` / `removeEmptyNodes`
**MATCH** ✓

Both handle empty leaves (remove) and empty internal nodes (reparent children to grandparent). Root never removed.

### 16. `getMixture` / `get_mixture`
**MATCH** ✓

Both compute weights via mass * main, propagate (1-main) * edge_weights to children.

### 17. Data Loading
**MATCH** ✓

SSM: id, gene, a, d, mu_r, mu_v parsed correctly.
CNV: id, a, d, ssm_refs with maternal_cn/paternal_cn parsed correctly.

---

## Bugs Summary

| # | Component | Severity | Description |
|---|-----------|----------|-------------|
| 1 | findOrCreateNode depth-0 | **MEDIUM** | Root index not prepended to path, causing asymmetric pathLT comparison in slice sampler |
| 2 | MH proposal | **UNKNOWN** | Go reimplements C++ mh.o; proposal distribution may differ |
| 3 | Slice sampler limit | LOW | 100 iteration cap vs Python's unlimited |
| 4 | boundBeta boundary | LOW | 1e-6 vs ~1e-16 clamping |
| 5 | LLH caching | LOW | No caching in slice sampler (performance only) |

---

## Acceptable Differences

1. **Go uses goroutines for parallel chains** — Python uses subprocess/threading
2. **Go uses Marsaglia-Tsang for gamma variates** — Python uses numpy
3. **Go has maxChildCreations safeguards** — prevents infinite loops
4. **Go outputs best_tree.json** — Python outputs trees.zip with pickled states
5. **Random seeds differ** — different RNG implementations, results are stochastic anyway
