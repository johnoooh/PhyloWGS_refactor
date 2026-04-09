# Simulation Validation Report 6 — HPC Head-to-Head

**Date:** 2026-04-09
**Branch:** `go-port` @ `110eea5` (Round 1 fixes + trace instrumentation)
**HPC data:** `simulation_validation4/` (216 fixtures × {go-cpu, original-python})

## TL;DR

The Phase 3 conclusion ("Go matches Python, close Round 2 as negative")
was **premature**. It was based on M=30/M=50 fixtures where parity
holds. The full 216-fixture HPC slice reveals severe, **M-dependent
over-splitting** in the Go port:

```
SSM count (M)    Δmedian Go−Python    Worst case
─────────────    ─────────────────    ──────────
     30              −1.2             Δ = −3  (Go tighter)
     50              −1.0             Δ = −2  (parity)
    100              −0.4             Δ = −4 / +3  (parity)
    150              +6.8             Δ = +69  (over-split)
    250             +12.9             Δ = +132 (catastrophic)
    500             +20.6             Δ = +126 (catastrophic)
```

34 of 171 paired fixtures have |Δmedian| > 2. All 27 fixtures with
|Δ| > 5 are M ≥ 150. The over-splitting is **not** correlated with
CNV count, true K, or S (subclone structure). **SSM count is the sole
driver.**

The Round 2 investigation must restart with an M≥150 fixture as the
canonical test case.

## Methodology

- **171 paired fixtures** where both Go and Python completed
  (Go: 207/216, Python: 171/216). Missing runs are mostly large
  K10/M500 fixtures that hit the 2hr slurm wall time.
- Effective K = `num_populations` from Go's `chain_*_trees.ndjson`
  (count of non-empty populations) and Python's `tree_summaries.json.gz`.
- Δmedian = Go median K − Python median K. Positive = Go over-splits.
- Analysis script: `simulation_validation4/hpc_head_to_head.py`

## Results

### Go completion: 207/216 (96%)

9 Go failures, all M500 or high-M/high-C combinations.

### Python completion: 171/216 (79%)

45 Python failures, scattered across K/M/C variants.

### Overall statistics

```
N paired fixtures:          171
Δmedian range:              −4 to +132
|Δmedian| mean:             7.02
|Δmedian| median:           1.0
Fixtures with |Δmedian| > 2:  34 (20%)
Fixtures with |Δmedian| > 5:  27 (16%)
```

The median |Δ| = 1.0 means the *typical* fixture matches. But the
*tail* is enormous — 20% of fixtures diverge by more than 2, and the
worst cases diverge by 130+.

### The M-dependence is the defining feature

```
SSMs (M)   N     Δmed mean   Δmed median   Δmed range    |Δ|>2
───────────────────────────────────────────────────────────────
  30       27      −1.22       −1.0        [−3, +0]       2
  50       32      −0.95       −1.0        [−2, +1]       0
 100       33      −0.36       −1.0        [−4, +3]       4
 150       30      +6.83       −1.0        [−2, +69]      8
 250       25     +12.88        0.0        [−1, +132]     8
 500       24     +20.58       +5.0        [−2, +126]    12
```

At M ≤ 100, Go is at parity or slightly tighter than Python (mean
Δ < 0). At M ≥ 150, Go systematically over-splits. The magnitude
grows roughly linearly with M.

### Not CNV-correlated

```
CNVs (C)   N     Δmed mean   |Δ|>2
──────────────────────────────────
  0        61      +7.38       9
  2        52      +5.05      14
  5        58      +4.02      11
```

The overnight local baseline (M=30 only) incorrectly attributed high
churn to CNV presence. At scale, CNV count has no predictive value —
the effect is purely M-driven.

### Not K-correlated (after controlling for M)

```
True K    N     Δmed mean   |Δ|>2
────────────────────────────────
  3       57      +3.02      10
  5       57      +5.11       8
 10       57      +8.46      16
```

K=10 shows the highest mean Δ, but K=10 also has the highest M values
(M=100 and M=500 variants). After controlling for M, the effect of
true K is negligible.

### LLH: Go finds better likelihoods, even when over-splitting

```
Best-LLH comparison: Go wins 155, Python wins 15, ties 1
Mean ΔLLH: +413.58 (Go higher)
```

Go achieves better data likelihood on 91% of fixtures. This includes
the over-split cases: the highest-Δ fixtures show Go LLH hundreds of
nats above Python. This is expected — more populations = more free
parameters = better data fit — and is precisely why the DP prior
exists (to penalize unnecessary populations). The fact that Go isn't
being regularized by the prior at high M is the bug.

### Top 10 outliers

```
Fixture                           Go_med  Py_med    Δ     Go_bestLLH  Py_bestLLH   ΔLLH
K5_S1_T200_M250_C0_rep1            137       5    +132    −1886.85    −1986.15    +99.29
K10_S1_T200_M500_C0_rep2           133       7    +126    −2360.03    −2516.25   +156.23
K5_S1_T200_M250_C2_rep0            102      15     +87    −1242.14    −1462.89   +220.75
K10_S1_T1000_M500_C0_rep2           77       4     +73    −3114.09    −3191.96    +77.87
K10_S1_T1000_M500_C0_rep1           74       4     +70    −2783.03    −2915.25   +132.21
K3_S1_T200_M150_C2_rep1             72       3     +69    −2079.80    −2167.45    +87.65
K10_S1_T200_M500_C5_rep1            74      13     +61    −2357.95    −2604.05   +246.10
K3_S1_T200_M150_C0_rep2             67      10     +57     −487.01     −778.61   +291.61
K5_S1_T1000_M250_C5_rep0            61       5     +56    −1513.85    −1600.93    +87.08
K10_S1_T1000_M500_C5_rep2           47      10     +37    −3505.32    −3527.31    +21.99
```

The worst fixture (K5_M250_C0_rep1) has Go median 137 populations vs
Python median 5. That is a 27× over-split. Go finds 99 nats better
likelihood — the extra 132 populations are fitting data, not air.

## Root cause hypothesis

**Mechanism:** With M SSMs, the slice sampler in `resampleAssignments`
calls `findOrCreateNode` approximately M times per MCMC iteration.
Each call walks the stick-breaking spine and may create new children.
If there is a per-call bias in Go vs Python — even a tiny one — it
compounds: M calls × per-call bias = O(M) excess nodes per iter.

This explains:
- Why M=30 and M=50 show parity (30 × ε ≈ 0)
- Why M=250 shows catastrophic divergence (250 × ε >> 0)
- Why the effect is purely M-correlated, not K/C/T-correlated

**Suspects (now re-prioritized for high-M):**

1. The `maxChildCreations=20` cap in `findOrCreateNode` was proven
   dead code on M=30. It may fire on M=250 where the per-iter total
   spawn count is ~250×higher. **MUST re-test cap_hits on M≥150.**

2. The `for iter < 100` slice-walk cap in `resampleAssignments` (Go)
   vs `while True` in Python. If the cap causes Go to accept suboptimal
   nodes more often, this becomes a per-datum bias that compounds.

3. Subtle floating-point differences in stick-product accumulation
   across many children. At M=250 with K≈100+ populations, the product
   `(1-ν₁)(1-ν₂)...(1-νₖ)` involves many more terms and numerical
   edge cases.

## Correction to Phase 3

Phase 3 concluded: "Go and Python produce the same K distribution on
K3_C2. The Round 2 premise is falsified." This was based on:

- K3_S1_T200_M30_C2_rep0 (M=30): Go median=12, Python median=13 ✓
- Overnight local baseline: 4 fixtures, all M=30 or M=50 ✓

The conclusion was true *for those fixtures* but not generalizable.
The Phase 3 test set was under-powered — it did not include any M≥150
fixture, which is where the divergence lives. The HPC slice at 216
fixtures has now revealed the full picture.

## Next steps

1. **Pick a canonical high-M fixture** for the next audit round. Good
   candidates: `K3_S1_T200_M150_C0_rep2` (Δ=+57, small K, no CNV,
   isolates the M effect) or `K5_S1_T200_M250_C0_rep1` (the Δ=+132
   worst case).

2. **Re-run the canonical fixture locally with `--trace`** to get
   per-iter spawn dynamics at M≥150. Key question: does
   `cap_hits_find_node` go non-zero?

3. **Restart the Phase 3 audit from `findOrCreateNode`** with the
   high-M fixture. The structural diff was correct at M=30, but
   behavioral divergence may only manifest when the function is called
   150+ times per iter.

4. **Test whether removing the `maxChildCreations=20` cap fixes the
   divergence.** This is the simplest possible intervention and was
   already identified as a non-Python divergence. Even though it didn't
   fire at M=30, the dynamics are different at M=250.

5. **Test whether removing the `for iter < 100` slice-walk cap fixes
   the divergence.** This is the other non-Python divergence in the
   slice sampler.

## Files

- `simulation_validation4/hpc_head_to_head.py` — analysis script
- `simulation_validation4/results/go-cpu/` — 207 Go fixtures
- `simulation_validation4/results/original-python/` — 171 Python fixtures
- This report: `PhyloWGS_refactor/sim_validation/SIM_VALIDATION_REPORT_6.md`
