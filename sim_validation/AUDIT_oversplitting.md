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

