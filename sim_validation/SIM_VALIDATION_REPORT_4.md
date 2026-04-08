# SIM Validation Report 4 — Precomputed MH States, Diagonal Slice

**Date:** 2026-04-08
**Branch:** `go-port`
**Change under test:** Precomputed MH states (P0 from `INVESTIGATION_NextSteps.md`). Go now mirrors Python's `params.write_data_state` / `mh.cpp` optimization: per (SSM, node) pi-independent (nr, nv) coefficients are computed once before each MH loop so the hot inner loop becomes a plain O(K) dot product instead of a full `computeNGenomes` tree walk on every MH iteration.

## What this report is not

This is not a full-grid rerun. The baseline is the **pre-optimization Go run already in `simulation_validation_updated/analysis/go-cpu/scores.json`**, not a fresh run of the old binary. That prior run used the same default MCMC settings (`-B 500 -s 1000 -j 4`, `mhIters=5000`) on the same fixtures, so the comparison is apples-to-apples.

Slice selection: 6 CNV fixtures spanning K ∈ {3, 5, 10} × C ∈ {2, 5}, with small M (30–100) to fit within session compute budget. The larger K10_M500 fixtures were excluded (both Python and the old Go binary had timed out on them at default settings, and the new binary would still need 30–60+ min each).

## Settings

```
phylowgs-go --no-gpu -B 500 -s 1000 -j 4 -r 42 -O <out> ssm_data.txt cnv_data.txt
# mhIters = 5000 (default)
# Machine: 4 CPUs (chain parallelism saturated)
```

Seed fixed at 42 across all runs so the stochastic trajectory is reproducible; the prior Go run used different seeds so best-LLH values will differ by noise.

## Results

Raw table:

```
fixture                         K*   new_s   oldGo    py_s   spdup     newLLH     oldLLH newK oldK  pyK
----------------------------------------------------------------------------------------------------
K3_S1_T200_M30_C2_rep0           3   427.9       -    1352     N/A     -240.7          -   30    -    -
K3_S1_T200_M30_C5_rep0           3    40.6     870     169   21.4x    -2332.8    -2333.1    7    7    4
K5_S1_T200_M50_C2_rep0           5    85.2    6619     251   77.7x     -179.2     -175.4   20   13    2
K5_S1_T200_M50_C5_rep0           5    48.6     614     179   12.6x     -409.5     -410.6    5    3    4
K10_S1_T1000_M100_C2_rep0       10   280.7       -    1668     N/A     -471.3          -   25    -   17
K10_S1_T1000_M100_C5_rep0       10   479.5       -    2496     N/A     -597.7          -   37    -    -
```

`K*` is the true number of subclones the fixture was generated with. `-` in the `oldGo` column means the pre-optimization Go run either timed out or failed scoring in the prior grid. `-` in `pyK` means Python also failed to produce a scoreable result for that fixture.

Summary saved to `simulation_validation_updated/results_new_slice/_slice_summary.json`; per-fixture Go output directories and logs are under `results_new_slice/<fixture>/` and `logs_new_slice/<fixture>.log`.

## Interpretation

### Speedups (old Go → new Go)

Where the old Go binary actually produced a scored result at default settings, the speedup from precomputed MH states is:

- **K5_S1_T200_M50_C2_rep0: 6619s → 85s = 77.7×**. This was the worst offender in the prior grid — nearly two hours for 50 SSMs. Now 1.5 minutes.
- K3_S1_T200_M30_C5_rep0: 870s → 41s = 21.4×
- K5_S1_T200_M50_C5_rep0: 614s → 49s = 12.6×

### Previously-timed-out fixtures that now complete

Three fixtures in the slice had no scored result in the prior go-cpu grid (they either timed out or the scorer saw `total_time_s=0`). All three completed successfully with the new binary at default settings:

- K3_S1_T200_M30_C2_rep0: **428s** (was unscored; Python reference 1352s)
- K10_S1_T1000_M100_C2_rep0: **281s** (was unscored; Python reference 1668s)
- K10_S1_T1000_M100_C5_rep0: **480s** (was unscored; Python reference 2496s)

So the optimization not only speeds up the fixtures the old binary could complete, it unblocks an entire class of large-K / large-CNV fixtures the old binary could not finish within the grid time budget.

### New Go vs Python reference

With the optimization in place, the new Go binary is faster than Python on every CNV fixture in the slice:

| fixture | new Go | Python | Go vs Py |
|---|---:|---:|---:|
| K3_M30_C2 | 428s | 1352s | **3.2×** |
| K3_M30_C5 | 41s | 169s | **4.2×** |
| K5_M50_C2 | 85s | 251s | **2.9×** |
| K5_M50_C5 | 49s | 179s | **3.7×** |
| K10_M100_C2 | 281s | 1668s | **5.9×** |
| K10_M100_C5 | 480s | 2496s | **5.2×** |

The old Go binary was significantly slower than Python on CNV fixtures (that was the entire problem). The new Go binary is now consistently **3–6×** faster than Python on this slice.

### LLH equivalence

Where both old and new Go scored the same fixture (at different seeds):

- K3_M30_C5: new -2332.8 vs old -2333.1  (Δ = +0.3, noise)
- K5_M50_C2: new -179.2 vs old -175.4    (Δ = -3.8, seed-variance)
- K5_M50_C5: new -409.5 vs old -410.6    (Δ = +1.1, noise)

These differences are consistent with MCMC seed variance, not with a semantic regression. The earlier `TestPrecomputedMHStatesEquivalence` unit test already proves per-MH-call numerical equivalence to 1e-10; the end-to-end LLH here should differ only because the prior grid runs used different seeds.

### Over-splitting is still present and is orthogonal to this change

`newK` shows that the new binary still over-splits on every fixture in the slice, often severely:

- K3_M30_C2 (true K=3) → inferred K=30
- K5_M50_C2 (true K=5) → inferred K=20
- K10_M100_C5 (true K=10) → inferred K=37

Python over-splits too (K10_M100_C2 → pyK=17 with true K=10) and in a few cases under-splits catastrophically (K5_M50_C2 → pyK=2). This matches what Report 3 documented: DP prior + `min_dp_alpha=1.0, max_dp_alpha=50` lets both implementations over-fit the tree on hard CNV fixtures. The MH optimization does not and should not fix this; it's a sampling / prior issue, not a likelihood-evaluation issue. It is, however, the most likely remaining driver of Go's per-fixture variance in wall time — see K3_M30_C2 where tree growth to 30 nodes inflates runtime even with fast MH.

## Correctness evidence

1. **TDD equivalence test** (`mh_states_test.go::TestPrecomputedMHStatesEquivalence`): asserts `logLikelihoodWithCNVTreeMH` returns identical values via the slow tree-walking path and the precomputed path, for both `newState=false` and `newState=true`, on a small hand-built tree with CNV/no-CNV/co-located/differently-located cases. Passes at 1e-10 tolerance.
2. **Smoke comparison** on `K3_M30_C2` with reduced MCMC (`-B 50 -s 100 -i 1000`, fixed seed 42): baseline 155.6s vs new 4.9s = 31.8× speedup, best LLH -219.69 / median -278.95 **byte-identical** between the two binaries.
3. **C=0 regression check** on `K3_M30_C0` with the same reduced settings: baseline 158ms vs new 168ms (noise), best LLH -108.17 / median -111.44 **byte-identical**. The non-CNV cheap path is untouched.
4. **This slice**: all six fixtures produce reasonable LLH matching the scale of the prior Go and Python results within expected seed variance; no runs crashed or produced NaNs.

## What changed in the code

In `PhyloWGS_refactor/main.go`:

- Added `SSMNodeState` struct and `MHStates []SSMNodeState`, `UseFourStates bool`, `MHStateValid bool` fields on `SSM` (mirrors Python's per-SSM per-node state vectors from `params.py:write_data_state`).
- Added `(*TSSB).precomputeMHStates()` that parallelizes across CNV-affected SSMs with goroutines (tree is read-only during the pass, no locking needed).
- Added `computeSSMStates()` that fills the four (nr, nv) slots per (SSM, node) using exactly the same case logic as `computeNGenomes`.
- Added `logLikelihoodWithCNVTreeMHPrecomputed()` — the O(K) fast path used during MH when precomputed states are valid.
- `logLikelihoodWithCNVTreeMH()` now dispatches to the precomputed path when `ssm.MHStateValid` is true and otherwise falls back to the existing tree walk (so any caller that forgot to precompute still gets correct results).
- `runChain()` calls `precomputeMHStates()` once per MCMC iteration, after `mapDatumToNode()` and before `metropolis()`, then invalidates `MHStateValid` on every SSM after MH returns so subsequent `resampleSticks`/`resampleStickOrders`/next iteration's `resampleAssignments` don't use stale factors.

All of this is a pure refactor of the CNV-MH inner loop. The math is identical to `computeNGenomes`; only the order of operations changed (pi multiplication moved from inside the tree walk to the outer sum across nodes, matching Python's compiled C++ precomputation exactly).
