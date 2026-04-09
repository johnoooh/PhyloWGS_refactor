# Local Trace Baseline Report — Round 2 Phase 3

**Date:** 2026-04-09 (overnight, before HPC validation results)
**Branch:** `go-port` @ `1128460`
**Purpose:** Build a calibration dataset of Go-port spawn dynamics on a
4-fixture local slice so the HPC validation results can be interpreted
in the morning.

## TL;DR

Three independent experiments (4 fixtures × HPC-default chain length;
4 fixtures × 5× chain length; 1 fixture × 4 parallel chains) all
**confirm the Phase 3 finding**: Go and Python produce effective K
distributions within ±1 of each other across the entire local slice.
The Round 2 over-splitting premise is falsified more broadly than the
single K3_C2 fixture used in Phase 3.

A new finding from this baseline: **the Go sampler internally carries a
much larger node tree than the effective population count it reports.**
On K3_C2 the sampler median is 49 nodes per iteration but only 12 of
those are non-empty populations. The other 37 are empty intermediate
spine nodes from the stick-breaking process. This is the source of the
"K_before_cull ≈ 100" trace numbers observed in Phase 3 — they were
counting spine nodes, not effective populations.

## Experiment matrix and runtimes

```
expA_baseline    HPC-parity (B=500/S=1000, 1 chain)       4 fixtures   121s
expB_long        Long convergence (B=2000/S=5000, 1 ch.)  4 fixtures   572s
expC_multichain  4-chain variance K3_C2 (B=500/S=1000)    1 fixture     90s
                                                            TOTAL      783s
```

All runs used `phylowgs-go` built locally from `go-port @ 1128460`
(includes Round 1 fidelity fixes + the `--trace` flag + the
`CapHitsFindNode` diagnostic counter).

## Headline result: effective K matches Python within ±1 across all 4 fixtures

Effective populations = `num_populations` field on each tree record =
count of nodes that contain at least one SSM or CNV. This is the
metric that Python's `tree_summaries.json.gz` also reports.

```
                       Python ref      Go expA      Go expB      Go expC
Fixture                median K        median K     median K     median K (avg of 4 chains)
─────────────────────────────────────────────────────────────────────────
K3_S1_T200_M30_C2_rep0    13              12           12          11.6
K3_S1_T200_M30_C0_rep0     3               2            2          —
K5_S1_T200_M50_C2_rep0     4               3            3          —
K5_S1_T200_M50_C0_rep0     6               7            6          —
```

Three of four fixtures show Go median *one less* than Python (slight
under-splitting). One fixture (K5_C0) shows Go median one more (slight
over-splitting). The arithmetic distance is symmetric across the slice
and within sampling noise. **No evidence of systematic Go-vs-Python
divergence.**

Best-LLH on K3_C2 with the longer chain (expB) reached **−325.0**,
much closer to Python's −315.8 than the short-chain expA value of
−347.4. Burn-in length matters; HPC defaults at B=500 are below the
sweet spot for K3_C2.

## Sampler vs effective K ratio — the source of the "K=30" mystery

For each fixture, the per-iter trace counts ALL nodes after `cullTree`
including empty intermediates. The reported `num_populations` only
counts non-empty nodes. The ratio quantifies how bushy the
stick-breaking spine is:

```
                       Sampler nodes   Effective K   Empty intermediate
Fixture                (median tree)   (median tree) ratio
─────────────────────────────────────────────────────────────────
K3_S1_T200_M30_C2_rep0      49             12         4.1×
K3_S1_T200_M30_C0_rep0       3              2         1.5×
K5_S1_T200_M50_C2_rep0       4              3         1.3×
K5_S1_T200_M50_C0_rep0      15              7         2.1×
```

**The "high churn" pattern is fixture-specific, not universal.** It is
very pronounced on K3_C2 (4.1×) and moderate on K5_C0 (2.1×). The other
two fixtures show near-1× ratios.

Plausible mechanism: K3_C2 has 30 SSMs distributed across deep
hierarchical structure with CNV constraints. The stick-breaking sampler
reaches deep nodes by spawning empty intermediate "rungs" along the way,
then leaves them in place because `cullTree` only removes *trailing*
empty children, not interior ones. `removeEmptyNodes` (the post-
processing pass that runs before writing tree NDJSON) removes the
intermediates from the reported tree. Whether Python's sampler also
carries these intermediates is unknown without instrumenting Python.

This explains the original "Go produces K=30" claim from
`SIM_VALIDATION_REPORT_5`: that number was almost certainly counting
sampler-side spine nodes rather than effective populations.

## Spawn dynamics — what's actually happening per iter

```
Fixture          spawns_find/iter    spawns_stick/iter   cap_hits/iter
                 (sample-phase mean) (sample-phase mean) (sample-phase total)
─────────────────────────────────────────────────────────────────────
K3_C2 expA           67.0                 18.9                 0
K3_C0 expA            2.3                  0.2                 0
K5_C2 expA            4.4                  0.5                 0
K5_C0 expA           21.4                  5.3                 0
```

Two confirmations:

1. **`cap_hits_find_node` = 0 across every iter of every experiment.**
   The non-Python `maxChildCreations = 20` safeguard never fires on any
   of these 4 fixtures, even on the fixtures with high spawn rates.
   Phase 3's central finding (the cap is dead code on this workload) is
   reconfirmed at scale.

2. **`spawns_find` >> `spawns_stick` everywhere.** Most child creation
   happens inside `findOrCreateNode` during the slice-sampler descent,
   not in `resampleStickOrders`. The 79% ratio from the original
   audit holds across all four fixtures.

## Convergence — does longer burnin help?

Comparing expA (B=500, S=1000) to expB (B=2000, S=5000) on K3_C2:

```
                  expA        expB        Δ
K_after_cull       47.7 ±20    44.2 ±14   −3.5  (similar mean, lower variance)
spawns_find        56.8 ±41    43.1 ±23   −13.7 (much less spawn churn)
best_llh         −347.4      −325.0      +22.4  (substantially better)
effective K (med)    12          12        0
```

**Effective K is stable across chain length** — the sampler converges
to its preferred K quickly. What longer burnin buys you is (a) better
likelihoods inside that K and (b) lower variance in the sampler-side
node count. K3_C0 / K5_C2 / K5_C0 all show essentially identical
results between expA and expB — they converge in well under 500 burnin.

## Multi-chain variance (Exp C)

4 independent chains on K3_C2 with HPC parameters (B=500/S=1000):

```
chain 0: median total=50  effective=12  best K=5 llh=−334.8
chain 1: median total=40  effective=11  best K=5 llh=−330.6
chain 2: median total=49  effective=12  best K=4 llh=−330.4
chain 3: median total=50  effective=12  best K=4 llh=−323.2
```

**Effective K is consistent across chains** (median 11-12). Sampler-side
total node count is also stable (40-50). Best-LLH varies modestly
(−323 to −335) — chain 3 explored a slightly better region. This is
exactly what good MCMC mixing looks like.

## What this tells us about tomorrow's HPC results

When HPC results land, expect:

1. **Effective K (`num_populations`) within ±1 of Python** for each
   fixture. If a fixture diverges by more than 1-2, that's a real
   signal worth investigating.

2. **`cap_hits_find_node = 0`** unless the HPC fixture set includes
   something dramatically larger than K3_C2. A non-zero cap_hits would
   contradict every local observation and warrant immediate
   investigation.

3. **Best-LLH should be in the same ballpark as Python**, possibly
   slightly worse on K3_C2-style fixtures if the HPC defaults
   (B=500/S=1000) under-burn for high-K hierarchies. K3_C2 with B=500
   reached −347 locally vs Python's −315; with B=2000 it reached −325.

4. **Per-fixture variability** is normal: the 4 fixtures here show
   wildly different spawn rates (K3_C2: 67/iter; K3_C0: 2.3/iter).
   Don't be alarmed by high spawn counts on dense fixtures unless
   they correlate with reported K divergence.

## Decision matrix for tomorrow

Bring the HPC `chain_*_trees.ndjson` files back and compute effective K
medians per fixture. Then:

| HPC outcome                                          | Action             |
|-------------------------------------------------------|--------------------|
| All fixtures within ±1 of Python                     | **Close Round 2.** Write Report 6 as a negative result. |
| 1 fixture diverges by 2-3                            | Re-run *that fixture* with `--trace` to localize. |
| Multiple fixtures diverge systematically             | Restart audit with HPC trace data as evidence. |
| Any `cap_hits_find_node > 0` on any chain            | Investigate immediately — contradicts local baseline. |

## Files

```
sim_validation/local_trace_runs/
  ├── run_local_baseline.sh        # the runner script
  ├── analyze_traces.py            # trace + posterior analyzer
  ├── explore_churn.py             # window/spearman explorer
  ├── REPORT.md                    # full analyzer output
  ├── CHURN.txt                    # full explorer output
  ├── POPSTATS.txt                 # sampler vs effective K table
  ├── BASELINE_REPORT.md           # this document
  ├── expA_baseline/<fixture>/{trace.ndjson, chain_0_trees.ndjson, ...}
  ├── expB_long/<fixture>/{trace.ndjson, chain_0_trees.ndjson, ...}
  └── expC_multichain/K3_S1_T200_M30_C2_rep0/{trace.ndjson.chain[0-3], chain_[0-3]_trees.ndjson, ...}
```

Total wall time: 783s (13.0 minutes). Total disk: ~6 MB of NDJSON.
