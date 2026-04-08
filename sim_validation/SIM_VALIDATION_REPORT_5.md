# SIM Validation Report 5 — Over-Splitting Python-Fidelity Fixes, Diagonal Slice

**Date:** 2026-04-08
**Branch:** `go-port`
**Change under test:** Three Python-fidelity corrections to TSSB MCMC sampling identified in the Phase 1 side-by-side audit (`sim_validation/AUDIT_oversplitting.md`) and implemented as a TDD cycle in commit `1ed4e3c`:

1. **BLOCKER 1d** — `dpAlphaLLH` root main-stick term: Python `tssb.py:247-255` evaluates `betapdfln(root.main, 1.0, dp_alpha)` at depth 0 whenever `min_depth <= depth` (default `min_depth=0`). Go's closure gated on `if depth >= 1`, silently dropping the term that the slice sampler depends on. Fix: extract as `TSSB.dpAlphaLLH` method with the correct `t.MinDepth <= depth` gate.

2. **BLOCKER 2a** — `resample_stick_orders` leftover bucket: Python `tssb.py:208-210` computes `hstack([all_weights[sub_indices], 1.0 - sum(all_weights)])`. The unallocated bucket is `1 - sum(ALL stick weights)`. Go computed `1 - sum(NOT-yet-ordered weights)`, inflating P(spawn) by a factor of ~3× per inner iteration on representative inputs. Fix: extract `stickOrderSubWeights` as a pure helper implementing Python's exact formula.

3. **SUSPICIOUS 5b** — `spawnChild` / `newTSSB` broadcast proportionality: Python `alleles.py:35-37` writes `self.pi = rand(1)*parent.pi`; `rand(1)` is a length-1 numpy array, and `(1,) * (ntps,)` broadcasting applies ONE scalar uniformly across all time-point dims. Go drew `rng.Float64()` inside the per-dim loop, producing independent fractions that broke proportionality. Fix: hoist the draw above the loop in both `spawnChild` and the duplicated `newTSSB` init block.

All three fixes were covered by `oversplit_test.go` (5 test cases) before being applied. Full `go test ./...` passes with no regressions.

## Settings

```
phylowgs-go-fixed --no-gpu -B 500 -s 1000 -j 4 -r 42 -O <out> ssm_data.txt cnv_data.txt
# mhIters = 5000 (default)
# Machine: 4 CPUs (chain parallelism saturated)
```

Same fixtures, same seed (42), and same MCMC settings as `SIM_VALIDATION_REPORT_4.md` so comparison is apples-to-apples.

## Results

```
fixture                         K*  prev_K  new_K  prev_s  new_s  prev_llh   new_llh    Δllh
-------------------------------------------------------------------------------------------
K3_S1_T200_M30_C2_rep0           3      30     30   427.9  425.0   -240.69   -241.66   -0.97
K3_S1_T200_M30_C5_rep0           3       7      6    40.6   42.5  -2332.80  -2333.47   -0.67
K5_S1_T200_M50_C2_rep0           5      20     19    85.2   84.5   -179.17   -176.97   +2.20
K5_S1_T200_M50_C5_rep0           5       5      5    48.6   48.9   -409.50   -413.22   -3.72
K10_S1_T1000_M100_C2_rep0       10      25     25   280.7  267.9   -471.31   -442.25  +29.06
K10_S1_T1000_M100_C5_rep0       10      37     40   479.5  429.0   -597.70   -609.04  -11.34
```

`K*` is the true number of subclones the fixture was generated with. `prev_*` columns are from `results_new_slice/_slice_summary.json` (the Report 4 baseline, which is the post-MH-precompute Go binary). `new_*` columns are from `results_phase6_slice/`.

Raw per-fixture JSON summaries and chain logs are under `simulation_validation_updated/results_phase6_slice/<fixture>/` and `logs_phase6_slice/<fixture>.log`. Per-run `best_tree.json` files contain the `num_populations` counts used for `new_K`.

## Interpretation

**Success criterion from the plan:** "at least 3 of 6 fixtures show lower K_inferred without a major LLH regression."

**Observed:**
- **Lower K:** 2 fixtures (K3_C5 7→6, K5_C2 20→19).
- **Unchanged K:** 3 fixtures (K3_C2 30→30, K5_C5 5→5, K10_C2 25→25).
- **Higher K:** 1 fixture (K10_C5 37→40).

**This does not meet the success criterion.** The fixes are real and Python-faithful, but they did not materially reduce the over-splitting symptom.

LLH deltas are all within ~30 nats, consistent with seed noise; there is no regression and no improvement beyond what a different random trajectory would produce. Runtime is essentially identical, as expected — none of the fixes touches the hot inner loop.

## What this tells us

The three audited bugs are **real Python divergences** and the fixes are **correct**, but the audit was **incomplete** as a theory of the K-inflation symptom. Evidence:

1. The `dpAlphaLLH` bug biases only the hyperparameter slice sampler's acceptance region; its effect on the tree structure is indirect (via `dp_alpha`), and the chain can self-correct at slice-sampling steps.
2. The `stickOrderSubWeights` bug was quantitatively large in isolation (50% vs 16.7% P(spawn) on the canonical input), but its effect only manifests when the outer `resample_stick_orders` loop actually reaches the "spawn new child" branch. Instrumenting the rate of spawn events in that branch would tell us whether it's even a hot path at the settings we're running.
3. The broadcast-scalar bug (5b) changes per-node variance of child `pi`, not the number of children the sampler retains. Its effect on K is at most second-order.

Put together, the audit found bugs that matter for Python-fidelity correctness but are not the dominant driver of over-splitting in practice. The real driver is still at large.

## What's still left

The next investigation should:

1. **Instrument the Go binary** to count (per MCMC iteration) (a) how many times `resample_stick_orders` reaches the spawn-new-child branch, (b) how many new children are created, and (c) how many of those survive `cull_tree`. This produces the quantitative trace needed to localize the symptom.
2. **Audit functions not covered in Phase 1**: `resample_sticks`, `resample_sigma`, and `find_node` in particular. The audit explicitly scoped out `find_node`'s lazy child creation loop; that is now the most likely candidate.
3. **Run a Python-vs-Go trace diff** on the smallest fixture (K3_C2, 30 SSMs) with identical seeds, comparing per-iteration `(K, best_llh, dp_alpha, alpha_decay, dp_gamma)`. Python's trajectory on this fixture terminates with a much smaller K than Go's; if we can identify the iteration where the trajectories first diverge, we have a pinpoint debug target.
4. **Consider the possibility that the bug is not in Go at all** but in the fixture scoring path (how `best_tree.json` is written or how `num_populations` is counted in the downstream analysis). This should be a quick sanity check before the next long investigation. Specifically: compare `num_populations` in `best_tree.json` to the Python reference's equivalent output on a matched run.

## Correctness evidence

- `go test ./...` clean on commit `1ed4e3c` (new `oversplit_test.go` added, no pre-existing tests regressed).
- `TestPrecomputedMHStatesEquivalence` (the MH-precompute pin test from `f969d36`) still passes — the refactor did not touch the hot inner loop.
- Each fix has an exact Python-source citation in the commit message and in `AUDIT_oversplitting.md`.
- RED-green-GREEN cycle was followed per TDD discipline: tests failed with compile errors on the extracted symbols before the fixes were applied, then failed semantically on the broadcast assertion, then passed after each fix was made in turn.

## What changed in the code

```
 main.go          | 156 ++++++++++++++++++++++++++++++++++++++++++++----
 oversplit_test.go | 200 ++++++++++++++++++++++++++++++++++++++++++++++++++++++
 2 files changed, 356 insertions(+), 50 deletions(-)
```

Commit: `1ed4e3c` on branch `go-port`. Files touched:
- `PhyloWGS_refactor/main.go`: `dpAlphaLLH` extracted, `resampleStickOrders` rerouted through `stickOrderSubWeights`, broadcast scalars hoisted in `spawnChild` and `newTSSB`.
- `PhyloWGS_refactor/oversplit_test.go`: new TDD tests.

## Decision

Per the plan's decision table: "partial success → next plan". The branch is NOT ready to merge as a fix for over-splitting. The fixes should stay (they are correct Python-fidelity corrections), but the investigation must continue in a new plan focused on the untouched `find_node` / `resample_sticks` code paths and the Python-vs-Go trace diff described above.
