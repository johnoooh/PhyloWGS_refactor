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
