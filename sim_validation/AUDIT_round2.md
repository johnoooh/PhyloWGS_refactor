# Round 2 AUDIT — Localizing the Dominant K-Inflation Source

**Date:** 2026-04-08
**Branch:** `go-port`
**Instrumentation commit:** `c3fb58b`
**Trace archive:** `sim_validation/round2_trace_K3C2.ndjson`

## Method

Built the Go binary at `c3fb58b` with the new `--trace` flag (Tasks 2.2 +
2.3) and ran a short MCMC on the K3_S1_T200_M30_C2_rep0 fixture:

```
phylowgs-go --no-gpu -B 50 -s 100 -j 1 -r 42 -i 500 \
    --trace trace.ndjson -O out \
    ssm_data.txt cnv_data.txt
```

150 iters total (50 burnin + 100 samples). The NDJSON log records, per
iter: K before/after `cullTree`, spawns inside `findOrCreateNode`, spawns
inside `resampleStickOrders`, DP hyperparams, and best llh.

## Headline numbers

```
N iters:                  150
fixture true K:             3
SSMs:                      30

K_before_cull   min=22   max=243   mean=95.6
K_after_cull    min= 7   max=103   mean=44.9
kills_in_cull   min= 1   max=179   mean=50.7

spawns_in_find_node       min=6 max=165 mean=52.6  TOTAL=7896  (79.0%)
spawns_in_stick_orders    min=0 max= 53 mean=14.0  TOTAL=2094  (21.0%)
TOTAL spawns                                        9990

dp_alpha     min=3.35   max=49.99  mean=22.33
alpha_decay  min=0.256  max=0.800  mean=0.523
best_llh     min=-841.3 max=-335.9 mean=-408.0
```

On a 30-SSM fixture with true K=3, the sampler averages **~95 nodes in
the tree before cullTree each iter** and still carries **~45 nodes after
culling**. Spawns average 66 per iter (53 from `findOrCreateNode`, 14 from
`resampleStickOrders`). The sampler is treading water: cullTree kills ~51
per iter, findOrCreateNode re-creates ~53 per iter, and stick-order
reshuffling adds another ~14.

## Attribution

**79% of all child spawns come from `findOrCreateNode`** (7896 out of
9990 over 150 iters). Only 21% come from `resampleStickOrders`. The
audit-and-fix work in Rounds 0+1 corrected three Python-fidelity bugs in
`resampleStickOrders` and `dpAlphaLLH`, which together address *at most*
21% of the spawn events. That explains the Report 5 outcome (2/6 fixtures
with slightly lower K, no meaningful symptom reduction): the fixes landed
on the minority contributor.

## Why `findOrCreateNode` is dominant

`findOrCreateNode` is called once per datum per slice-sampler iteration
inside `resampleAssignments`. With 30 SSMs and a loose slice distribution,
each iter issues 30+ descents. Each descent's inner loop (main.go
lines 2343-2375, depth>0 branch) lazily spawns children until either (a)
`u` falls into the existing stick space or (b) a local safeguard cap
`maxChildCreations = 20` is hit.

Observed: mean spawns_in_find_node = 53/iter ≈ 1.77 spawns per datum per
iter. Peak was 165/iter, which is within reach of 30 data × 5 spawns/datum.
K_before_cull peaked at 243, roughly 30 × 8 — consistent with the cap
being hit at multiple depths simultaneously.

## Candidates for the bug inside `findOrCreateNode`

Two suspicious code paths relative to Python's `tssb.py:find_node` (340-377):

1. **The `maxChildCreations = 20` safeguard (main.go:2346)** is NOT in
   Python. Python's `find_node` loops without a cap until `u < edges[-1]`
   (tssb.py:365-368). When the Go cap triggers (line 2365:
   `if u < edge || creations >= maxChildCreations`), the descent proceeds
   with a truncated stick set and the new child is appended without `u`
   actually falling into its range. This can leave orphan nodes: the
   datum may be assigned into a child whose stick does not actually cover
   the slice sampler's `u`, which then gets re-walked next iter and more
   children get spawned on top of the existing ones.

2. **The depth==0 branch (main.go:2409-2420)** unconditionally hard-codes
   `index := 0` and does NOT do lazy child creation at root. Python's
   `find_node` uses the same code for both depth==0 and depth>0 (single
   while loop at tssb.py:358). This Go divergence means Go never creates
   new root-level children via find_node, only via the initial `newTSSB`
   spawn and via `resampleStickOrders`. But the data then all flows into
   the single depth-1 sub-tree, which then gets flooded with deeper
   children in the depth>0 inner loop — exactly matching the observed
   flat structure `0 -> [1..30]` post-`removeEmptyNodes` (from Task 1.3).

   The flat post-cull structure is strong corroborating evidence: deep
   intermediate nodes get reparented to root by `removeEmptyNodes`, and
   the reported fan-out of 30 siblings is exactly the set of depth-N
   leaves holding one SSM each.

## Null hypotheses tested and rejected

- **"DP hyper slice sampler is bad"**: alpha_decay ranges 0.256-0.800
  (mean 0.52), dp_alpha ranges 3.35-49.99 (mean 22.3). These are within
  reasonable ranges for a well-behaved sampler. Not a stuck slice.
- **"cullTree is broken"**: it kills 51 nodes/iter on average, including
  up to 179 at once. It is actively working. The problem is spawn rate,
  not cull rate.
- **"resampleStickOrders is the driver"**: accounts for only 21% of
  spawns. Even if we drove this to zero, expected K would still inflate
  from the 79% that find_node contributes.

## Decision

**Proceed to Phase 3:** line-by-line audit of `findOrCreateNode` against
Python's `tssb.py:find_node`, with the two suspicious paths above as
the primary targets. Items to confirm or reject against Python source:

1. Python has no child-creation cap → Go's `maxChildCreations = 20` is
   an unconditional divergence. Quantify what happens when it triggers
   (silent orphan append? silent fallback? wrong datum assignment?).

2. Python treats depth==0 and depth>0 in the same while loop → Go's
   special-case root branch (main.go:2409-2420) is a divergence.
   Confirm whether Python's `find_node` at depth 0 also goes through
   the lazy-create loop, and whether that loop can create NEW root
   children.

3. Python passes `min_depth` and stops early if `depth < min_depth`.
   Confirm the Go path honors the same min_depth behavior (default 0,
   so currently a no-op, but important if min_depth is ever raised).

Phase 4 (TDD fix + slice re-run + Report 6) follows Phase 3.

---

## Phase 3 Results — Both hypotheses falsified, paradigm shift

**Date:** 2026-04-08 (same day, later)

### Falsification 1: the `maxChildCreations = 20` cap

Added a `CapHitsFindNode` diagnostic counter that increments every time
the safeguard fires *before* `u` falls into the stick space. Re-ran the
K3_C2 fixture short chain (50 burnin + 100 samples) and a long chain
(500 burnin + 500 samples). Result:

```
cap_hits_find_node:   total = 0   across all 150 + 1000 iters
```

The Go-only cap **never fires in practice**. It cannot be the source of
K-inflation — it is dead code on this workload. Hypothesis rejected.

### Falsification 2: the depth==0 special case

Re-read `tssb.py:find_node` (340-377). At depth 0, Python also
hard-codes `index = 0` and descends into `root['children'][0]`. Python's
`find_node` *does not* create new root children either. The Go
special-case branch (main.go:2409-2420) is structurally equivalent to
Python's depth-0 behavior. Hypothesis rejected.

### Full structural diff of `findOrCreateNode` vs `find_node`

Line-by-line against `phylowgs/tssb.py:340-377` confirmed equivalence in:

- depth-0 index hard-coding
- `u` rescaling on the remaining stick mass at each depth
- lazy `spawnChild` append when `u > edges[-1]`
- the `boundbeta(1, alpha_decay^depth * dp_alpha)` stick draw
- no min_depth handling in either (both default to 0)
- no gamma-prior or DP-prior sampling inside this function

The **only** divergence is the `maxChildCreations=20` cap, which the
counter just proved is unreachable on this fixture. No structural bug
remains to fix in `findOrCreateNode`.

### Head-to-head on K3_C2 (4000 Go trees vs 1000 Python trees)

```
                        Python (orig)        Go (current head)
N trees:                1000                 4000 (4 chains × 1000)
range (K):              4 – 16               3 – 15
median K:               13                   12
mean K:                 11.95                11.42
mode K:                 13 (279/1000)        12 (1228/4000)
best-llh tree K:        5                    4
best-llh value:         −315.8               −307.85   ← Go better
```

**Go and Python produce essentially the same K distribution on K3_C2.**
Go is marginally *tighter* (median 12 vs 13) and finds a slightly better
likelihood. Both samplers over-split the posterior relative to the true
K=3 — because PhyloWGS's DP prior with dp_alpha ∈ [1,50] routinely
admits 10+ populations on 30-SSM data. This is a property of the model,
not a bug in either port.

### The Round 2 premise was wrong

The round 1 / round 2 framing of "Go over-splits relative to Python" on
the K3_C2 fixture is **not supported by the data**. Both implementations
land in the same K=10–14 band. The high per-iter `K_before_cull ≈ 80`
observed by the tracer is a normal burnin/slice-sampler dynamic shared
by both ports; it is not evidence of a fidelity bug because `cullTree`
immediately reclaims ~half the spawned nodes and the retained `K_after`
aligns with Python.

### What the Round 1 / Report 5 "K=30" number actually measured

The original "Go produces K≈30" observation came from the
post-`removeEmptyNodes` flat structure reported in Task 1.3. That
function reparents deep intermediate-empty nodes to the root, producing
an artificial fan-out of leaf populations. Report 5's K count was
measuring *reporting-layer* flatness, not sampler-layer fanout. The
sampler itself carries a hierarchical K ≈ 12, matching Python.

### Decision

**Do not proceed to Phase 4 without user direction.** The premise for
the Phase 4 fix (K-inflation is a Go bug in `findOrCreateNode`) is
falsified. Options:

1. **Abandon the over-splitting investigation.** Go matches Python on
   K3_C2 within sampling noise. Close out round 2 with a negative
   result and revise the sim-validation report.
2. **Broaden the head-to-head.** Run Go + Python on K5/K8 fixtures to
   confirm parity holds beyond K3_C2. If a real divergence appears on
   other fixtures, restart the audit on that evidence.
3. **Investigate the `removeEmptyNodes` flat-structure reporting.** If
   the user cares about the post-processing flatness (K=30 in Report
   5), the fix belongs in the reporting layer, not the sampler.

Phase 3 closes here pending user direction.

