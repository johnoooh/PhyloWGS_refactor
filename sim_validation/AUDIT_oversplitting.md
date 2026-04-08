# Over-Splitting Audit — Go vs Python Side-by-Side

**Date:** 2026-04-08
**Branch:** go-port
**Reference Python:** /sessions/dreamy-pensive-knuth/mnt/phylowgs_cowork/phylowgs/

This document records every concrete difference between the Go port and
the upstream Python on the suspect functions identified in
docs/plans/2026-04-08-oversplitting-investigation.md.

## Format

For each function:
- Python file:line range
- Go file:line range
- Differences observed (line-by-line diff in prose)
- Severity: BLOCKER (will measurably affect over-splitting) /
            SUSPICIOUS (semantic divergence, may matter) /
            BENIGN (cosmetic / equivalent rewrite)
- Whether to fix in this plan

## Functions audited

### 1. `resampleHypers` — DP hyperparameter slice samplers

- **Python:** `phylowgs/tssb.py:245-316` (`TSSB.resample_hypers`)
- **Go:** `PhyloWGS_refactor/main.go:2352-2447` (`(*TSSB).resampleHypers`)

#### 1a. Hyperparameter bounds — BENIGN

Python class constants at `tssb.py:10-15`:
```
min_dp_alpha = 1.0, max_dp_alpha = 50.0
min_dp_gamma = 1.0, max_dp_gamma = 10.0
min_alpha_decay = 0.05, max_alpha_decay = 0.80
```
Go hard-codes the same values at `main.go:2354-2355, 2391-2392, 2415-2416`. Bit-identical. BENIGN.

#### 1b. Sampling order — BENIGN

Both order the three slice samplers as: dp_alpha → alpha_decay → dp_gamma. Matches.

#### 1c. 100-iteration cap in Go, `while True` in Python — SUSPICIOUS

- Python uses `while True` (`tssb.py:261, 278, 305`) and explicitly raises `"Slice sampler shrank to zero!"` if `new_param == current_param` (the else branch at `tssb.py:270-271, 287-288, 314-315`).
- Go caps the inner loop at 100 iterations (`main.go:2376, 2396, 2434`). If no proposal exceeds `llhSlice` in 100 tries, Go silently exits **without updating the hyperparameter**, leaving it at the previous iteration's value. No error, no warning.

**Reasoning about magnitude:** A standard shrinking-interval slice sampler should accept within O(log2(range/tolerance)) ≈ 10-20 iterations in healthy cases because each rejection halves the search interval. 100 iters is >5× the typical need. However, this behavior does silently differ from Python in two cases:
  1. LLH surface is pathological and the slice sampler legitimately needs >100 iters (rare).
  2. Shrink-to-zero (`new_param == current_param`): Python raises (hard stop), Go silently keeps shrinking via `upper = newAlpha`, then eventually hits 100 and returns without updating.

**Severity: SUSPICIOUS.** Not the smoking gun I hoped, but worth empirically checking in Phase 3 (is `dp_alpha` pinned on one side of the interval in Go's traces?).

#### 1d. Missing depth-0 term in `dpAlphaLLH` — **BLOCKER**

This is the first real divergence.

Python `dp_alpha_llh` at `tssb.py:247-255`:
```python
def descend(dp_alpha, root, depth=0):
    llh = betapdfln(root['main'], 1.0,
                    (alpha_decay ** depth) * dp_alpha) if self.min_depth <= depth else 0.0
    for child in root['children']:
        llh += descend(dp_alpha, child, depth + 1)
    return llh
```
With `self.min_depth == 0` (the default set in `TSSB.__init__` at `tssb.py:17-23`, and `evolve.py:85` never overrides it), the condition `self.min_depth <= depth` is `True` for all `depth >= 0` including the root at `depth=0`. So **Python evaluates `betapdfln(root.main, 1.0, dp_alpha)` at the root**.

Go `dpAlphaLLH` at `main.go:2357-2370`:
```go
descend = func(root *TSSBNode, depth int) float64 {
    llh := 0.0
    if depth >= 1 {
        llh = betaPDFLn(root.Main, 1.0, math.Pow(t.AlphaDecay, float64(depth))*alpha)
    }
    for _, child := range root.Children {
        llh += descend(child, depth+1)
    }
    return llh
}
```
Go explicitly gates on `depth >= 1`, **excluding** the root's `main` stick. This is not a translation of Python's condition: Python's gate is `min_depth <= depth` (always True with default `min_depth=0`), Go's is `depth >= 1`. These are not equivalent; Go's behavior corresponds to `min_depth == 1`, which is **not** the PhyloWGS default.

**Impact on dp_alpha sampling:** The dp_alpha slice sampler compares `new_llh(alpha') > log(rand()) + llh(alpha_current)`. Go drops one beta term from both sides of this comparison. If that term were a constant in `dp_alpha`, the drop would cancel. **It is not constant.** The depth-0 term is `betapdfln(root.main, 1.0, alpha_decay^0 * dp_alpha) = betapdfln(root.main, 1.0, dp_alpha)`, which depends directly on `dp_alpha`. So Go's sampler is targeting a subtly different posterior than Python's.

**Impact on alpha_decay sampling:** The depth-0 term is `betapdfln(root.main, 1.0, dp_alpha)` — it does not depend on `alpha_decay` (since `alpha_decay^0 = 1`). So the term is constant across alpha_decay proposals and cancels out of the slice comparison. The alpha_decay sampler is equivalent.

**Impact on dp_gamma sampling:** `dp_gamma_llh` doesn't touch `main` at all (only sticks), so this divergence is irrelevant to the dp_gamma sampler. Equivalent.

**Magnitude:** In a tree with N nodes, the dp_alpha llh has N terms from Python (including root) and N-1 from Go. For small trees this is ~10-20% of the llh signal; for large over-split trees (N=30-80) it's ~1-3%. On the root's `main` stick, `betapdfln(x, 1.0, dp_alpha)` is concave in `dp_alpha` with a peak that depends on `x = root.Main`. Dropping a concave term that depends on `dp_alpha` biases the sampler; the direction of bias depends on where root.Main sits. In PhyloWGS the root holds no data, so root.Main is driven by the prior alone and should sit close to the beta mean.

**Severity: BLOCKER.** This is a semantic divergence from Python, by the plan's definition ("will measurably affect over-splitting"). The magnitude may be small but it's present on *every* iteration of *every* chain, and it's exactly the kind of slow bias that could compound into pinned hyperparameters and excess K growth. Also note: per the plan's Phase 1.8 decision rule, **one BLOCKER skips Phases 2-4 and goes straight to Phase 5 TDD fix**. I will re-check whether this is *the* cause after writing the fix, but either way it needs fixing.

**Fix sketch:** Change Go's `dpAlphaLLH` to include the depth-0 term whenever `t.MinDepth <= 0` (Go should track `min_depth` as a TSSB field mirroring Python, with default 0). Minimal fix: unconditionally include `depth == 0` in Go, matching Python's default. Long-term fix: add `MinDepth int` to TSSB and make both places (`dpAlphaLLH` and `TSSB.Root.Main` initialization) consult it.

#### 1e. Closure capture of `alpha_decay` in Go — BENIGN (but fragile)

Python passes `(dp_alpha, alpha_decay)` explicitly into `dp_alpha_llh`, then calls `dp_alpha_llh(self.dp_alpha, new_alpha_decay)` inside the alpha_decay sampler (line 280). Go's `dpAlphaLLH` closure captures `t.AlphaDecay` by reference from the receiver, so the alpha_decay sampler has to mutate-and-restore:
```go
oldDecay := t.AlphaDecay
t.AlphaDecay = newDecay
newLLH := dpAlphaLLH(t.DPAlpha)
t.AlphaDecay = oldDecay
```
Functionally equivalent in single-threaded code. Fragile if the TSSB were ever read concurrently, but `resampleHypers` is called on a single chain's TSSB from a single goroutine, so this is fine. BENIGN.

#### 1f. Shrink-toward-current comparison — BENIGN

Python uses `elif new_* < self.*: lower = new_*; elif new_* > self.*: upper = new_*; else: raise`. Go uses `if new_* < t.*: lower = new_*; else: upper = new_*`. The only observable difference is the missing shrink-to-zero raise, already covered in 1c. BENIGN.

#### 1.R Summary

| Sub-item | Severity | Fix in this plan |
|---|---|---|
| 1a bounds | BENIGN | no |
| 1b order | BENIGN | no |
| 1c 100-iter cap | SUSPICIOUS | investigate empirically (Phase 3) |
| **1d missing depth-0 term in dpAlphaLLH** | **BLOCKER** | **yes** |
| 1e closure capture | BENIGN | no |
| 1f shrink comparison | BENIGN | no |

Phase 1.8 rule: one BLOCKER found → skip Phases 2-4 after audit is complete and proceed to Phase 5. **I will still finish the rest of Phase 1 audits first** because there may be additional BLOCKERs worth batching into the same fix commit, and the audit document is cheap.

---

### 2. `resampleStickOrders` — child reorder + lazy spawn

- **Python:** `phylowgs/tssb.py:193-243` (`TSSB.resample_stick_orders`)
- **Go:** `PhyloWGS_refactor/main.go:2071-2248` (`(*TSSB).resampleStickOrders`)

The algorithm: for each non-leaf node, look at its data-bearing children (the `represented` set); for each, draw `u ~ Uniform(0,1)` and walk the stick-breaking mass to find which child (or uncreated region) u falls in. If u falls in an uncreated region, spawn a new child and retry; if u falls in an existing child, append it to `new_order` and move on. At the end, the parent's `children` list is permuted into `new_order` and the rest are killed.

#### 2a. Leftover mass computation — **BLOCKER**

This is the second concrete bug and is likely the single largest over-splitting driver.

**Python `tssb.py:208-210`:**
```python
sub_indices = filter(lambda i: i not in new_order, range(root['sticks'].shape[0]))
sub_weights = hstack([all_weights[sub_indices], 1.0 - sum(all_weights)])
sub_weights = sub_weights / sum(sub_weights)
```
The "uncreated region" bucket is `1.0 - sum(all_weights)` — the **true** remaining stick mass that has never been allocated to any existing child. After normalization, `sub_weights` sums to 1 over `len(sub_indices)+1` buckets.

**Go `main.go:2133-2152`:**
```go
subWeights := make([]float64, len(subIndices)+1)
totalUsed := 0.0
for i, idx := range subIndices {
    subWeights[i] = allWeights[idx]
    totalUsed += allWeights[idx]
}
subWeights[len(subIndices)] = 1.0 - totalUsed // leftover
```
Go computes `totalUsed` as `sum(allWeights[i] for i in subIndices)` — the mass of the **not-yet-ordered** children. The "leftover" bucket is then `1 - totalUsed`, which equals `sum(already_ordered_children) + true_unallocated`. **Go's leftover bucket is inflated by the mass of children that have already been placed into `new_order`**, then normalized as if it were all genuinely unallocated.

**Quantitative verification:** With 3 children with stick weights [0.4, 0.3, 0.2] (sum 0.9, true unallocated = 0.1) and child 0 already in `new_order`:

| | sub_weights before norm | P(fall into leftover bucket) |
|---|---|---|
| Python | [0.3, 0.2, 0.1] (sum 0.6) | 0.1/0.6 = **16.7%** |
| Go | [0.3, 0.2, 0.5] (sum 1.0) | **50.0%** |

Go is 3× more likely to spawn a new child on this step. The bug only manifests when the outer loop `for len(represented) > 0` has iterated at least once (i.e. a parent with ≥2 data-bearing children), so it compounds as the tree grows. **This explains why Go over-splits more than Python on fixtures with many true subclones.**

**Severity: BLOCKER.** Direct, quantifiable bias toward spawning new children, with trivial fix: change the Go line `subWeights[len(subIndices)] = 1.0 - totalUsed` to `subWeights[len(subIndices)] = 1.0 - sum(allWeights)` (precomputed before the for loop). Then normalization by sumW (which will no longer be trivially 1) produces Python-equivalent probabilities.

**Relationship to the `maxCreations = 10` cap in 2b below:** the cap is almost certainly there because the original port hit runaway spawning in internal testing, which is precisely what this bug causes. Fixing 2a will likely make 2b unnecessary, but we should retain the cap for now (it acts as a safety net) and remove it in a follow-up after Phase 6 validates the primary fix.

#### 2b. `maxCreations = 10` cap, no Python equivalent — SUSPICIOUS

Python's inner spawn loop (`tssb.py:207`) is `while True`. Go caps at 10 creations per `u` draw (`main.go:2109`). When hit, Go force-picks an existing child rather than allowing the spawn loop to continue. This is a one-directional divergence that **reduces** spawn count in Go. In isolation this would cause Go to under-split, but it's masked by 2a (which biases toward spawning) and serves as a safety net against runaway spawning when 2a fires. Severity: SUSPICIOUS, do not change in this plan — removing it without fixing 2a would cause infinite loops.

#### 2c. Kill-by-index loop semantics — BENIGN

Python `tssb.py:233-237` iterates `k not in new_order` and deletes from `root['children'][k]` while iterating — Python list semantics mean this could misbehave, but in practice it uses `filter` to build the list of indices first and then calls `kill()` then `del`. Go `main.go:2229-2238` builds `inNewOrder` set first, then iterates original indices and kills each. Functionally equivalent. BENIGN.

#### 2d. Immediate `resampleSticks` at end — BENIGN

Both Python `tssb.py:243` and Go `main.go:2247` call `resampleSticks` immediately after the reorder pass. Matches.

#### 2e. `resampleSticks` depth-0 main override — BENIGN (cosmetic divergence)

Tangential to stick-orders but worth noting since it was in the Python source I read:

Python `tssb.py:186-187` does `root['main'] = boundbeta(...) if min_depth <= depth else 0.0; if depth == 0: root['main'] = 1e-30`. So at depth 0 Python always calls `boundbeta` (consuming an RNG draw) then overwrites with 1e-30.

Go `main.go:1959-1963` branches on `depth >= 1` and only calls `boundBeta` then; at depth 0 it directly sets `1e-30` without consuming an RNG draw.

Net effect on `Main` value: identical (1e-30 at depth 0, boundBeta at depth > 0). Only difference is a 1-draw offset in the RNG stream at each `resampleSticks` call. Since Go's RNG stream is independent of Python's anyway (`math/rand` vs numpy), this is observationally unobservable. BENIGN.

#### 2.R Summary

| Sub-item | Severity | Fix in this plan |
|---|---|---|
| **2a leftover mass bucket** | **BLOCKER** | **yes** |
| 2b maxCreations=10 cap | SUSPICIOUS | keep as safety net, follow-up |
| 2c kill-by-index | BENIGN | no |
| 2d resampleSticks call | BENIGN | no |
| 2e depth-0 Main override | BENIGN | no |

---

### 3. `resampleAssignments` + `findOrCreateNode`

- **Python:** `phylowgs/tssb.py:82-150` (`TSSB.resample_assignments`) + `tssb.py:340-377` (`TSSB.find_node`)
- **Go:** `PhyloWGS_refactor/main.go:1268-1417` (`(*TSSB).resampleAssignments`) + `main.go:2253-2350` (`(*TSSB).findOrCreateNode`)

#### 3a. Slice sampler 100-iter cap — SUSPICIOUS

Python `tssb.py:115` uses `while True`. Go `main.go:1307` caps at 100 iterations. Same pattern as 1c. Python has an explicit escape at `tssb.py:132-137` when `abs(max_u - min_u) < epsilon` (reverts to old state and prints a warning). Go has the same epsilon escape at `main.go:1346-1348` but silently. Since in Go the move is either never made (non-CNV) or restored before the check (CNV), the "keep current assignment" behavior on collapse/timeout matches Python. Severity: SUSPICIOUS on parity grounds; no direct over-splitting link. No fix this plan.

#### 3b. `maxChildCreations = 20` cap in `findOrCreateNode` — SUSPICIOUS

Python `tssb.py:353-361` has `while not root['children'] or (1.0 - prod(1.0 - root['sticks'])) < u:` — unbounded. Go caps at 20 (`main.go:2272, 2291`). Same band-aid pattern as 2b. When the cap is hit, Go's subsequent index-finding loop picks `index = len(edges) - 1` (the last child), which differs from Python's behavior (which guarantees edge > u on exit). In practice this cap is rarely hit because 20 creations at one node would already be pathological. SUSPICIOUS, keep as safety net, no fix this plan.

#### 3c. No `llhmap` cache in Go — BENIGN (performance)

Python `tssb.py:99, 125-129` caches per-node LLH across slice-sampler iterations to avoid recomputation when the sampler proposes the same node multiple times. Go recomputes every iteration. This is a performance optimization, not a correctness divergence. BENIGN, no fix.

#### 3d. SSM and CNV datum assignments resampled in separate passes — BENIGN

Python treats SSMs and CNVs as rows of a single `self.data` array and resamples them interleaved in a single loop. Go separates SSMs (`t.Data`) and CNVs (`t.CNVData`) into two sequential passes (`main.go:1277-1361` for SSMs, `main.go:1366-1417` for CNVs). Because each slice-sampler update targets its own conditional posterior, the sequential-vs-interleaved ordering doesn't change the stationary distribution. The RNG stream is different, but Go's RNG stream is independent of Python's anyway. BENIGN.

#### 3e. Root-node empty convention — BENIGN

Python `tssb.py:118-120` and Go `main.go:1314-1317` both redirect the root to `root.children[0]` with `new_path = [0]` when `find_node` returns root. Since both set `root.Main = 1e-30` in `resampleSticks`, the root is essentially never returned except at `max_depth`, and the redirect is equivalent. BENIGN.

#### 3f. Data move semantics during slice sampling — BENIGN

Python `tssb.py:121-124` moves the datum from old_node to new_node before computing LLH, then reverts on rejection. Go does the move only for CNV SSMs (`main.go:1323-1331`) where it matters for the CNV-aware tree likelihood, and leaves non-CNV SSMs alone (non-CNV LLH is a pure function of `newNode.Params`, which isn't affected by datum membership). Both approaches produce the same LLH value; Go's is more efficient. BENIGN.

#### 3g. Path-based slice bracket shrinking — BENIGN

Python `tssb.py:142-146` computes `path_lt(indices, new_path)` and shrinks `min_u`/`max_u` based on sign. Go `main.go:1354-1359` does the same via `pathLT`. Would need to confirm `pathLT` matches Python's `path_lt` exactly; it's reportedly verified in earlier reports. BENIGN assuming prior verification stands. (Sanity check: both treat `pathComp >= 0` identically — Go's `else` branch covers both `== 0` and `> 0`, matching Python's `elif path_comp >= 0`. Good.)

#### 3.R Summary

| Sub-item | Severity | Fix in this plan |
|---|---|---|
| 3a 100-iter cap | SUSPICIOUS | no |
| 3b maxChildCreations=20 in findOrCreateNode | SUSPICIOUS | no (safety net) |
| 3c no llhmap cache | BENIGN (perf) | no |
| 3d SSM/CNV separate passes | BENIGN | no |
| 3e root-empty redirect | BENIGN | no |
| 3f data-move semantics | BENIGN | no |
| 3g path-based shrinking | BENIGN | no |

**No BLOCKER found in resampleAssignments or findOrCreateNode.** All divergences are either performance optimizations or cap/safety-net patterns whose removal would require fixing BLOCKER 2a first.

---

### 4. `cullTree` and `killNode`

- **Python:** `phylowgs/tssb.py:152-166` (`TSSB.cull_tree`) + `phylowgs/alleles.py:45-50` (`alleles.kill`)
- **Go:** `PhyloWGS_refactor/main.go:1993-2029` (`(*TSSB).cullTree`) + `main.go:1977-1991` (`killNode`)

#### 4a. `cullTree` trim-trailing semantics — BENIGN

Python uses `trim_zeros(counts, 'b')` which trims trailing zeros only (NOT arbitrary internal empty children). Go replicates with a reverse loop that finds the highest index with non-zero data and sets `keep = i + 1`. Traced through several cases: all-zero, mixed-trailing, single-survivor. Semantically equivalent. Previously flagged in Report 3 as "already verified" — reconfirmed. BENIGN.

#### 4b. Recursion order — BENIGN

Both Python and Go recurse into children BEFORE computing `keep`, so empty deepest descendants die first and a surviving grandchild's data is still counted in the parent's total (preventing a data-bearing subtree from being killed from above). Matches. BENIGN.

#### 4c. `killNode` — pi return to parent — BENIGN

Python `alleles.kill()` at `alleles.py:45-50`:
```python
def kill(self):
    if self._parent is not None:
        self._parent._children.remove(self)
    self._parent.pi = self._parent.pi + self.pi
    self._parent   = None
    self._children = None
```

Go `killNode` at `main.go:1977-1991`:
```go
func killNode(child *Node, parent *Node) {
    for i := range parent.Pi {
        parent.Pi[i] += child.Pi[i]
    }
    for i, c := range parent.Children {
        if c == child {
            parent.Children = append(parent.Children[:i], parent.Children[i+1:]...)
            break
        }
    }
}
```

Both return the dead child's full `pi` to the parent and remove the child from the parent's children slice. Go does not null out `child._parent` and `child._children` since Go's GC handles the freed node, and no subsequent code paths should see a reference to the dead `child`. BENIGN.

#### 4d. Sticks truncation on cull — BENIGN

Python `tssb.py:161-162` truncates `sticks` and `children` to `keep`. Go `main.go:2016-2018` does the same with an extra length check. Note that after cullTree, both Python and Go rely on `resampleSticks` to repopulate any missing stick entries on the next iteration, so even a tiny off-by-one here would be corrected on the next loop. BENIGN.

#### 4.R Summary

| Sub-item | Severity | Fix in this plan |
|---|---|---|
| 4a trim-trailing | BENIGN | no |
| 4b recursion order | BENIGN | no |
| 4c killNode pi return | BENIGN | no |
| 4d sticks truncation | BENIGN | no |

**No BLOCKER in cullTree or killNode.**

---

### 5. Tree initialization — `newTSSB` / evolve.py startup

- **Python:** `phylowgs/tssb.py:17-46` (`TSSB.__init__`) + `phylowgs/evolve.py:84-98` (hack block) + `phylowgs/alleles.py:17-37` (`alleles.__init__`)
- **Go:** `PhyloWGS_refactor/main.go:589-661` (`newTSSB`) + `main.go:2031-2066` (`spawnChild`)

#### 5a. Initial `root.Main` value — BENIGN (transient)

- Python: root main is initialized to `boundbeta(1.0, dp_alpha)` in `TSSB.__init__` (tssb.py:30) and then overridden to `1e-30` on the first `resample_sticks` call (tssb.py:187). So Python's root main is non-trivial for a brief window before MCMC starts.
- Go: root main is initialized directly to `1e-30` in `newTSSB` (main.go:602) and stays there.

Net effect after the first MCMC iteration: identical (both `1e-30`). The "before-first-iteration" state is only observable if the code computes a likelihood before running any MCMC moves, which neither implementation does. BENIGN.

#### 5b. Initial child's `pi` — **SUSPICIOUS** (broader scope than init alone)

Python `alleles.__init__` at `alleles.py:35-37`:
```python
self.pi = rand(1)*parent.pi
parent.pi = parent.pi - self.pi
self.params = self.pi
```

`rand(1)` returns a length-1 array; broadcasting `(1,) * (ntps,)` produces a length-`ntps` vector where **every element is the same random scalar** times the corresponding parent pi. Verified with numpy REPL: `np.random.rand(1) * np.ones(5)` → `[0.374, 0.374, 0.374, 0.374, 0.374]`. So Python uses **one** random fraction per spawn, applied uniformly across all time-points.

Go `spawnChild` at `main.go:2045-2050` (and the duplicated init block in `newTSSB` at `main.go:632-637`):
```go
for i := range childNode.Pi {
    frac := rng.Float64()
    childNode.Pi[i] = frac * parent.Node.Pi[i]
    parent.Node.Pi[i] -= childNode.Pi[i]
    childNode.Params[i] = childNode.Pi[i]
}
```

Go draws a **fresh** `rng.Float64()` **per time-point**, so the child's pi vector is **not proportional** to the parent's pi — instead, each dimension has an independent random scaling factor.

**Scope of divergence:** this applies every time `spawnChild` is called, which is: (a) the one-time initial child in `newTSSB`, (b) every new child in `resampleStickOrders`, (c) every new child in `findOrCreateNode`. So the bug is not limited to initialization; it fires every time the tree grows.

**Impact on over-splitting:** indirect but real. Node params in PhyloWGS are used by the MH move to score phi-proposal acceptance, and more importantly by the CNV-aware tree-likelihood. If newly spawned children's pi vectors don't inherit the parent's across-sample structure proportionally, they get placed at arbitrary points in the pi-simplex and the MH move has to do extra work to pull them toward plausible values. This is a fidelity issue regardless of over-splitting impact.

**Severity: SUSPICIOUS (BLOCKER-adjacent).** The magnitude for over-splitting is unclear without empirical testing, but per the "match Python's logic" directive this is a divergence and should be fixed. I'm marking it SUSPICIOUS (not BLOCKER) because it doesn't directly inflate node counts — it affects the *content* of newly spawned nodes, not their *creation rate*. But since fixing it is a trivial 2-line change and since we're already writing a Phase 5 fix commit for BLOCKER 2a and BLOCKER 1d, I recommend including this fix in the same commit.

**Fix sketch:** In `spawnChild` and the duplicated `newTSSB` block, draw `frac` **once** outside the for-loop, apply the same value across all `i`. Matches Python's broadcast semantics exactly.

#### 5c. Initial child stick value 0.999999 — BENIGN

- Python evolve.py:89: `boundbeta(1, dp_gamma) if depth!=0 else .999999` — with depth=0, uses 0.999999.
- Go main.go:646: `0.999999` directly.
Match. BENIGN.

#### 5d. Initial child `Main` — BENIGN

- Python evolve.py:91: `boundbeta(1.0, (alpha_decay**(depth+1))*dp_alpha) if min_depth <= depth+1 else 0.0`, with depth=0 and min_depth=0: `boundbeta(1.0, alpha_decay * dp_alpha)`.
- Go main.go:641: `boundBeta(1.0, alphaDecay*dpAlpha, rng)`.
Match. BENIGN.

#### 5e. Initial data assignment — BENIGN

- Python evolve.py:95-98 moves all data from root to the new child.
- Go main.go:649-653 directly assigns all SSMs to the new child (skipping the root-first-then-move indirection).
Net effect: identical. BENIGN.

#### 5f. CNV data assignment at init — BENIGN (Go-specific detail)

Go explicitly assigns all CNVs to the initial child at `main.go:656-658`. Python treats CNVs as rows of the same `data` array and assigns them with the rest. Net effect: identical initial placement. BENIGN.

#### 5.R Summary

| Sub-item | Severity | Fix in this plan |
|---|---|---|
| 5a initial root.Main | BENIGN | no |
| **5b child pi spawn (per-dim vs broadcast scalar)** | **SUSPICIOUS (include in fix)** | **yes** |
| 5c initial stick | BENIGN | no |
| 5d initial child Main | BENIGN | no |
| 5e data assignment | BENIGN | no |
| 5f CNV placement | BENIGN | no |

---

### 6. MCMC loop ordering — `runChain` vs `do_mcmc`

- **Python:** `phylowgs/evolve.py:159-213` (`do_mcmc` iteration body)
- **Go:** `PhyloWGS_refactor/main.go:2519-2587` (`runChain` iteration body)

#### 6a. Call order — matches

| Step | Python (evolve.py) | Go (main.go) | Notes |
|---|---|---|---|
| 1 | `resample_assignments()` | `resampleAssignments(rng)` | |
| 2 | `cull_tree()` | `cullTree()` | |
| 3 | `get_mixture()` + assign `node.id` | — | Go assigns node ids lazily elsewhere |
| 4 | `set_node_height` | `setNodeHeights()` | |
| 5 | `set_path_from_root_to_node` | `setNodePaths()` | |
| 6 | `map_datum_to_node` | `mapDatumToNode()` | |
| 7 | — | `precomputeMHStates()` | Go-only optimization from commit f969d36; no semantic effect |
| 8 | `metropolis(...)` | `metropolis(mhIters, mhStd, rng)` | |
| 9 | MH std adapt (×2 / /2) | MH std adapt (×2 / /2) | same thresholds 0.08/0.5/0.99 |
| 10 | — | `ssm.MHStateValid = false` for all SSMs | Invalidates the Go-only MH cache |
| 11 | `resample_sticks()` | `resampleSticks(rng)` | |
| 12 | `resample_stick_orders()` | `resampleStickOrders(rng)` | (internally re-calls resampleSticks in both) |
| 13 | `resample_hypers()` | `resampleHypers(rng)` | |
| 14 | `complete_data_log_likelihood()` | `completeDataLogLikelihood()` | |

Ordering matches 1:1. The `precomputeMHStates` + invalidation around metropolis is a pure caching optimization and has no effect on the posterior. The missing explicit `get_mixture + node.id` step in Go is only used by Python for logging/output purposes and does not affect sampling. BENIGN.

#### 6b. Burnin semantics — BENIGN

Go uses `for iter := -burnin; iter < samples; iter++` and discards samples with `iter < 0`. Python uses `range(start_iter, num_samples)` with `start_iter = -burnin - 1 + 1 = -burnin`. Equivalent burnin counting. BENIGN.

#### 6.R Summary

No BLOCKERs. MCMC loop ordering matches Python.

---

## Decision

### Tally

| Severity | Count | Items |
|---|---:|---|
| BLOCKER | 2 | 1d (dpAlphaLLH skips root main term), 2a (leftover-mass bucket wrong in resampleStickOrders) |
| SUSPICIOUS | 5 | 1c (100-iter cap in resampleHypers), 2b (maxCreations=10 in stickOrders), 3a (100-iter cap in resampleAssignments), 3b (maxChildCreations=20 in findOrCreateNode), 5b (spawnChild per-dim random vs Python broadcast scalar) |
| BENIGN | 17 | 1a 1b 1e 1f 2c 2d 2e 3c 3d 3e 3f 3g 4a 4b 4c 4d 5a 5c 5d 5e 5f 6a 6b |

### Decision rule from the plan

> If there is even one BLOCKER, the next step is "Phase 5: TDD fix" — skip Phases 2-4.

Two BLOCKERs found. Skipping Phase 2 (diagnostic instrumentation), Phase 3 (comparative run), Phase 4 (targeted experiment). Proceeding directly to **Phase 5: TDD fix**.

### What to fix in the Phase 5 commit

Three changes, all small, all strictly tightening Go's fidelity to Python:

1. **BLOCKER 1d** — `dpAlphaLLH` in `main.go:2357-2370`: replace `if depth >= 1` with a gate based on `t.MinDepth`. Since `t.MinDepth == 0` is always set in `newTSSB`, the condition simplifies to unconditionally including the root's `main` term, matching Python's default. Add a TSSB field `MinDepth int` if not already present (it IS present per `main.go:618`) and gate on `t.MinDepth <= depth`.

2. **BLOCKER 2a** — leftover mass bucket in `resampleStickOrders` at `main.go:2133-2142`: change the computation from `1.0 - totalUsed` (sum of subIndices' weights) to `1.0 - sumAllWeights` (precomputed sum of all stick weights). Then the weights no longer pre-normalize to 1 trivially, and the existing `sumW` normalization step will do the right thing.

3. **SUSPICIOUS 5b** — `spawnChild` at `main.go:2045-2050` and the duplicated init block in `newTSSB` at `main.go:632-637`: hoist `frac := rng.Float64()` out of the for-loop so one scalar is applied uniformly across all `ntps` dimensions. Matches Python's `rand(1) * parent.pi` broadcast semantics.

### What NOT to fix in the Phase 5 commit (deferred)

- **1c, 3a** (100-iter caps on slice samplers): leave in place. Not confirmed as a cause, and removing them requires confidence that the slice samplers will always converge — best done after 1d/2a are fixed and the trajectory traces in a subsequent session can empirically verify.
- **2b, 3b** (maxCreations safety caps on spawn loops): leave in place. These are band-aids for bug 2a; removing them without an independent verification pass could expose infinite loops if 2a isn't fully fixed. Remove in a follow-up PR after Phase 6 validates the main fix.
- **3c** (missing llhmap cache): performance optimization only, no correctness impact.

### Tests to write in Phase 5

1. **`TestDPAlphaLLHIncludesRootMain`** (for 1d): Build a minimal TSSB with a known root.Main, compute `dpAlphaLLH(alpha)` for two different `alpha` values. Independently compute `betaPDFLn(root.Main, 1.0, alpha)` for the depth-0 contribution and assert it's included in the sum. Before the fix, the test should FAIL because Go's function excludes the depth-0 term.

2. **`TestResampleStickOrdersLeftoverMass`** (for 2a): Construct a parent with 3 children with known stick weights (e.g. [0.4, 0.3, 0.2]) and put child 0 into a simulated `new_order`. Compute both Python's and Go's leftover-bucket probability. Assert Go matches Python (0.167, not 0.5). This may be easier as a pure-function unit test of the leftover-weight computation extracted out of `resampleStickOrders`.

3. **`TestSpawnChildBroadcastScalar`** (for 5b): Build a parent node with known pi vector (e.g. [0.3, 0.6, 0.9] for ntps=3). Seed the RNG, call `spawnChild`, assert `child.Pi[i] / parent.Pi[i]` is the same constant for all `i` (within float tolerance). Before the fix, this ratio will differ across dimensions.

All three tests can live in a new file `oversplit_test.go` (rather than `mh_states_test.go`) per the Phase 5 plan instructions.







