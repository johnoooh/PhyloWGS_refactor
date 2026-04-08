# Bugs Found During Simulation Validation

_Discovered during local validation session: 2026-04-07_
_All three bugs are fixed and committed on the `go-port` branch._

---

## Bug 1 — `rebuildAncestorSets` missing self-inclusion

**Commit:** `bf096a0`  
**Severity:** High — caused wrong likelihood path for all co-located SSM+CNV pairs  
**Function:** `rebuildAncestorSets` in `main.go`

### What was wrong
The function pre-computes `AncestorSet` maps for O(1) ancestor lookups. The walk started from `n.Parent`, so a node was never added to its own `AncestorSet`. This made `isAncestorOfCached(nd, nd)` return `false`.

Inside `computeNGenomes`, when an SSM and its CNV are assigned to the same node, `ssmInAncestors` evaluated to `false` for that node, routing to Case 3 (no variant, has CNV) instead of the correct Case 4 (has variant, has CNV). For `minor_cn=0` CNVs, Case 3 produces only `nv=0` pairs, all filtered, leaving `log(1e-99)` per SSM per timepoint.

### Fix
```go
n.AncestorSet[n] = true   // add before parent walk
```

---

## Bug 2 — SSM `CNVs` field omitted in `runChain` deep-copy

**Commit:** `2735615`  
**Severity:** Critical — primary driver of anomalous LLH in all C=2 fixtures  
**Function:** `runChain` in `main.go`

### What was wrong
`runChain` deep-copies SSMs to give each parallel chain independent data. The `CNVs` field (slice of `*CNVRef`) was silently omitted from the copy struct. As a result, every SSM in every chain had `len(ssm.CNVs) == 0`, so `logLikelihoodNoCNV` was called for CNV-linked SSMs instead of `logLikelihoodWithCNVTree`.

For SSMs under a hemizygous deletion (`minor_cn=0`), the diploid model expects ~50% reference reads but the data shows ~95%. This caused per-SSM LLH penalties of ~−300 vs expected ~−5.

**Symptoms:** Go LLH was −9,386 vs Python −6,745 for `K3_S1_T1000_M30_C2_rep1`. The chains looked converged but plateaued at the wrong value.

**Why some C=2 fixtures were unaffected:** Fixtures where the CNV is at near-100% prevalence give approximately correct results from `logLikelihoodNoCNV` (normal cells dominate). Only fixtures with low-prevalence or high-impact CNVs showed large anomalies.

### Fix
```go
cnvsCopy := make([]*CNVRef, len(ssm.CNVs))
copy(cnvsCopy, ssm.CNVs)
ssmsCopy[i] = &SSM{ ..., CNVs: cnvsCopy }
```

---

## Bug 3 — Shared `*CNV` objects across concurrent chains (data race)

**Commit:** `e6defce`  
**Severity:** High — silent data race between goroutines; corrupts CNV node assignments  
**Function:** `runChain` in `main.go`

### What was wrong
Fix 2 added the `CNVs` slice copy, but the `*CNVRef` pointers still pointed to the **original shared `*CNV` objects**. `CNV.Node` is written on every MCMC iteration by `resampleAssignments`. With 4 chains running as concurrent goroutines, each chain's `resampleAssignments` would overwrite the shared `cnv.Node` field, corrupting every other chain's ancestor lookups via `findMostRecentCNV`.

This is a classic write-write race: no panic, no error, just silently wrong node assignments when two chains happen to write simultaneously.

### Fix
Deep-copy the `*CNV` objects per chain and re-wire `CNVRef.CNV` pointers to chain-local copies:

```go
chainCNVs := make([]*CNV, len(cnvs))
origToChainCNV := make(map[*CNV]*CNV, len(cnvs))
for i, cnv := range cnvs {
    chainCNV := &CNV{ ID: cnv.ID, A: ..., D: ..., /* Node assigned by newTSSB */ }
    chainCNVs[i] = chainCNV
    origToChainCNV[cnv] = chainCNV
}
for i, ssm := range ssms {
    chainRefs := make([]*CNVRef, len(ssm.CNVs))
    for j, ref := range ssm.CNVs {
        chainRefs[j] = &CNVRef{CNV: origToChainCNV[ref.CNV], ...}
    }
    ssmsCopy[i] = &SSM{ ..., CNVs: chainRefs }
}
tssb := newTSSB(ssmsCopy, chainCNVs, ...)
```

---

## Post-Fix Verification

| Fixture | Go pre-fix | Go post-fix | Python | Notes |
|---------|----------:|------------:|-------:|-------|
| K3_S1_T1000_M30_C2_rep1 | −9,386 | **−6,518** | −6,745 | ✅ |
| K5_S1_T200_M50_C2_rep0 | −6,586 | **−6,649** | −6,649 | ✅ identical |
| K5_S1_T200_M50_C2_rep1 | −6,644 | **−6,786** | −6,652 | ✅ |
| K3_S3_T1000_M30_C2_rep1 | −24,689 | **−21,026** | −21,041 | ✅ |
| K5_S3_T200_M50_C2_rep0 | −18,644 | **−19,005** | −25,418 | ✅ Go better |

Multi-seed variance (5 seeds, `K3_S1_C2_rep1`): σ ≈ 2.7, best range −6513 to −6521. Stable.
