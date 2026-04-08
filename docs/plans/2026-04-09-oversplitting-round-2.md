# Over-Splitting Investigation Round 2 — Trace-Driven Localization

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Localize the dominant driver of K-inflation in the Go PhyloWGS port by (a) ruling out the scoring path as a false source, (b) instrumenting per-iteration spawn/retain counts, (c) producing a Python-vs-Go per-iteration trace diff on the smallest fixture, and (d) auditing the two hot paths scoped out of Round 1 (`findOrCreateNode`, `resampleSticks`).

**Architecture:** This is an investigation plan, not a feature plan. Most "tasks" produce evidence (logs, diffs, audit docs), not shipped code. Any Go code changes are **strictly Python-fidelity fixes**, introduced via TDD. No custom regularization, no φ-pruning, no `min_ssms` thresholds, no DP-prior bound tweaks.

**Tech Stack:**
- Go 1.x (PhyloWGS_refactor/main.go, monolithic)
- Python 2.7 (phylowgs/, reference implementation)
- fixture: `simulation_validation_updated/fixtures/K3_S1_T200_M30_C2_rep0` (smallest; baseline K=3, Go infers K=30)

---

## Background — what Round 1 showed

`sim_validation/SIM_VALIDATION_REPORT_5.md` (commit `bf924da`) ran the fixed binary from commit `1ed4e3c` on the 6-fixture diagonal slice. Results:

```
fixture                     K*  prev_K  new_K    Δllh
K3_S1_T200_M30_C2_rep0       3     30     30   -0.97
K3_S1_T200_M30_C5_rep0       3      7      6   -0.67
K5_S1_T200_M50_C2_rep0       5     20     19   +2.20
K5_S1_T200_M50_C5_rep0       5      5      5   -3.72
K10_S1_T1000_M100_C2_rep0   10     25     25  +29.06
K10_S1_T1000_M100_C5_rep0   10     37     40  -11.34
```

Two fixtures moved by 1 (noise); four didn't meaningfully change. The three bugs from Round 1 (`dpAlphaLLH` root term, `stickOrderSubWeights` leftover bucket, `spawnChild` broadcast scalar) are **real Python divergences and the fixes stay**, but they are not the dominant driver.

## What Round 1 did NOT investigate (the scope of this round)

1. **`findOrCreateNode` (main.go:2300-~2380)** — lazy child creation during slice-sampled assignment walks. Spot-checked during audit but flagged as the most likely next candidate. Contains a non-Python safeguard `maxChildCreations := 20` that may silently emit orphan children when hit.
2. **`resampleSticks` (main.go:1937-1968)** — inspected briefly in Round 1 under assumption of equivalence; re-audit with fresh eyes.
3. **The scoring path** — `best_tree.json` writer at `main.go:2970-3015` and `summarizePops`. `num_populations` is computed as "non-empty pops" (nodes with ≥1 SSM or CNV). Python (`phylowgs/pwgsresults/result_generator.py:39`) uses `load_trees_and_metadata(remove_empty_vertices=True)` then `len(mut_assignments)`. These may not be identical; must be verified before trusting the reported `K_inferred`.
4. **No trace-level empirical evidence** of which MCMC operation actually accumulates the extra nodes on a real fixture.

## Guiding principle

**Evidence before code.** Phase 1-3 gather evidence; Phase 4 only writes code if Phase 1-3 localize a Python divergence. Any fix must follow the Round 1 TDD cadence (red-green, pinned test, commit before moving on).

---

## Phase 1 — Rule out the scoring path (fast sanity check)

### Task 1.1: Inspect Go's `num_populations` computation

**Files:**
- Read: `PhyloWGS_refactor/main.go:2970-3015` (summarizePops)
- Read: `PhyloWGS_refactor/main.go:2810-2830` (best_tree.json writer)

**Step 1:** Read the exact Go code that computes `num_populations`, `populations`, and `mut_assignments` for `best_tree.json`.

**Step 2:** Record in a scratch file `/tmp/round2_scoring_notes.md`:
- How Go counts non-empty populations
- Whether Go's tree-to-pops walk includes empty intermediate nodes in the tree structure (even if not counted)
- Whether `num_populations` and `len(mut_assignments)` would produce the same number

### Task 1.2: Inspect Python's scoring path

**Files:**
- Read: `phylowgs/pwgsresults/result_generator.py:37-90` (_summarize_all_pops, _summarize_pops)
- Read: `phylowgs/util2.py` — the `TreeReader.load_trees_and_metadata` method and `remove_empty_vertices` parameter

**Step 1:** Follow the `remove_empty_vertices=True` code path. Determine exactly what gets removed:
- Intermediate empty nodes collapsed?
- Leaf empty nodes deleted?
- How are children reparented?

**Step 2:** Append findings to `/tmp/round2_scoring_notes.md`.

### Task 1.3: Construct a small synthetic tree and compare counts

**Step 1:** Write a 30-line Python script `/tmp/compare_counts.py` that:
- Loads the `best_tree.json` from `results_phase6_slice/K3_S1_T200_M30_C2_rep0/`
- Counts populations three ways: `d['num_populations']`, `len(d['mut_assignments'])`, `sum(1 for p in d['populations'].values() if p['num_ssms']>0 or p['num_cnvs']>0)`
- Prints all three and the tree structure

**Step 2:** Run it: `python3 /tmp/compare_counts.py`

**Expected outputs to interpret:** If all three numbers agree at 30, the scoring path is consistent and the symptom is real. If they disagree, investigate which one Python's scoring would produce and fix the Go reporter.

**Step 3:** If the scoring path IS the bug, stop this plan and go straight to Phase 4 with a narrow fix. Otherwise continue to Phase 2.

### Task 1.4: Decision point — scoring or sampling?

**Step 1:** Write a 2-3 sentence conclusion in `/tmp/round2_scoring_notes.md`.

**Step 2:** If scoring: skip to Phase 4 with a dedicated "fix scoring" task. If sampling: continue to Phase 2.

---

## Phase 2 — Instrument the Go binary (spawn/retain trace)

### Task 2.1: Design the trace

**Files:**
- Scratch: `/tmp/round2_trace_design.md`

**Step 1:** Write down the trace schema. Per MCMC iteration, record:
- `iter`
- `K_before_cull` = number of nodes before `cullTree`
- `K_after_cull` = number of nodes after `cullTree`
- `spawns_in_find_node` = count of new children created inside `findOrCreateNode` during `resampleAssignments`
- `spawns_in_stick_orders` = count of new children created inside `resampleStickOrders`
- `kills_in_cull` = count of children removed by `cullTree`
- `dp_alpha`, `dp_gamma`, `alpha_decay`, `best_llh_this_iter`

**Step 2:** Decide the output format: one-line-per-iter JSON to a trace file `<out>/trace.ndjson` when a new CLI flag `--trace` is set.

### Task 2.2: Add a `--trace` flag

**Files:**
- Modify: `PhyloWGS_refactor/main.go` — CLI parsing block (search for `flag.String` or `flag.Bool`)

**Step 1:** Write a failing test in `PhyloWGS_refactor/trace_test.go` that runs `main.go` with `--trace` on a 2-iter mini run and asserts the trace file exists and has ≥2 lines of valid JSON.

**Step 2:** Run: `go test -run TestTraceFlag ./...`
Expected: FAIL (flag not wired).

**Step 3:** Add a `traceFile string` CLI flag. In `runChain`, if non-empty, open/append to the file at start and close at end. Write nothing yet — just wire the plumbing.

**Step 4:** Re-run test. Expected: PASS (empty file is acceptable for this step).

**Step 5:** Commit: `feat: add --trace flag for per-iteration diagnostics`

### Task 2.3: Add counters inside the hot paths

**Files:**
- Modify: `PhyloWGS_refactor/main.go` — `findOrCreateNode`, `resampleStickOrders`, `cullTree`, `runChain`

**Step 1:** Write a failing integration test `TestTraceFieldsPopulated` in `trace_test.go` that runs a 3-iter mini run and asserts trace lines contain the 8 schema fields from Task 2.1.

**Step 2:** Run test. Expected: FAIL (fields missing).

**Step 3:** Add a `TraceRecord` struct on the TSSB state or threaded through the iter loop. Increment counters at:
- Every `root.Children = append(..., newChild)` inside `findOrCreateNode` (`main.go:2327, 2346`) → `spawns_in_find_node++`
- Every new-child branch inside `resampleStickOrders` → `spawns_in_stick_orders++`
- `cullTree`: compute `kills_in_cull` as `before - after` node count.
- `runChain` per-iter tail: snapshot counts, write the JSON line, reset counters.

**Step 4:** Run test. Expected: PASS.

**Step 5:** Commit: `feat: per-iter spawn/cull counters in --trace output`

### Task 2.4: Run traced mini-chain on K3_C2

**Step 1:**
```bash
/tmp/phylowgs-go-fixed-traced --no-gpu -B 50 -s 100 -j 1 -r 42 --trace \
    -O /tmp/round2_trace_k3c2 \
    fixtures/K3_S1_T200_M30_C2_rep0/ssm_data.txt \
    fixtures/K3_S1_T200_M30_C2_rep0/cnv_data.txt
```

(Single chain, 50+100=150 iterations, seed 42.)

**Step 2:** Summarize the trace with a small Python one-liner:
```bash
python3 -c "
import json
lines = [json.loads(l) for l in open('/tmp/round2_trace_k3c2/trace.ndjson')]
for l in lines[::10]:
    print(f\"iter={l['iter']:3d} K_pre={l['K_before_cull']:3d} K_post={l['K_after_cull']:3d} \
find_spawns={l['spawns_in_find_node']:3d} so_spawns={l['spawns_in_stick_orders']:3d} \
kills={l['kills_in_cull']:3d} dpa={l['dp_alpha']:.2f} llh={l['best_llh_this_iter']:.1f}\")
"
```

**Step 3:** Save the output to `sim_validation/TRACE_round2_k3c2.txt`.

### Task 2.5: Localize the dominant spawn source

**Step 1:** Look at `TRACE_round2_k3c2.txt`. Compute:
- mean spawns per iter from `findOrCreateNode`
- mean spawns per iter from `resampleStickOrders`
- mean kills per iter from `cullTree`
- net K drift = mean(K_after_cull[t+1] - K_after_cull[t])

**Step 2:** Write a 1-paragraph finding in `sim_validation/AUDIT_round2.md` (create it). Identify which path dominates K inflation. Mark it "PATH X" for Phase 3.

---

## Phase 3 — Audit the dominant path + Python trace diff

### Task 3.1: Line-by-line audit of PATH X

**Files:**
- Create: `sim_validation/AUDIT_round2.md` (already started in 2.5)
- Read: relevant Go function + relevant Python function side-by-side

**Step 1:** For PATH X, produce the same scaffold used in `AUDIT_oversplitting.md`:
- Python section with exact line citations
- Go section with exact line citations
- Divergence list with severity (BLOCKER / SUSPICIOUS / BENIGN) and quantitative reasoning

**Known candidate divergence to verify:** `findOrCreateNode` has `maxChildCreations := 20` (main.go:2319) — Python `find_node` has no such cap. When the cap is hit, Go falls through to the cumulative-edges picker even though `u > edge`. Quantify: how often is the cap hit on K3_C2 during a traced run?

**Step 2:** Instrument and count cap-hits. Either extend trace or run a quick ad hoc print. Record in audit doc.

### Task 3.2: Python trace instrumentation

**Files:**
- Modify: `phylowgs/tssb.py` — add a temporary `print >> sys.stderr` in `find_node`, `resample_stick_orders`, and `resample_hypers` at the same counter locations
- Modify: `phylowgs/evolve.py` — emit one JSON line per MCMC iter with the same 8 fields as Go's trace

**Step 1:** Make the Python changes in a dirty local copy only (do NOT commit to the Python side). These are throwaway diagnostics.

**Step 2:** Run Python PhyloWGS on the same fixture with the same seed (or closest achievable):
```bash
cd /sessions/dreamy-pensive-knuth/mnt/phylowgs_cowork/phylowgs
python2 evolve.py -B 50 -s 100 -r 42 \
    ../simulation_validation_updated/fixtures/K3_S1_T200_M30_C2_rep0/ssm_data.txt \
    ../simulation_validation_updated/fixtures/K3_S1_T200_M30_C2_rep0/cnv_data.txt \
    2> /tmp/py_trace_k3c2.ndjson
```

**Step 3:** Sanity-check the Python trace exists and is parseable. If Python 2 is unavailable in the shell, note this in the audit and skip to Task 3.4 (compare summary statistics instead of per-iter).

### Task 3.3: Per-iter diff

**Files:**
- Create: `/tmp/diff_traces.py`

**Step 1:** Write a Python 3 script that loads both traces and prints:
- First iter where `K_after_cull` diverges by more than 2
- Cumulative spawn counts (find_node vs stick_orders) at that iter
- Cumulative kill counts at that iter
- `dp_alpha`, `dp_gamma`, `alpha_decay` at that iter

**Step 2:** Run it. Save output to `sim_validation/TRACE_DIFF_round2.txt`.

**Step 3:** Interpret: which counter first diverges materially? That is the proximate cause.

### Task 3.4: Decision point

**Step 1:** Write a 2-paragraph conclusion in `sim_validation/AUDIT_round2.md`:
- Which function is the dominant K-inflation source
- Whether it's a Python divergence (fixable) or a legitimate difference due to seed / implementation-specific RNG ordering
- Concrete hypothesis for Phase 4

**Step 2:** If no Python divergence is found, STOP and escalate to the user. Do not invent a fix.

---

## Phase 4 — Fix the localized divergence (TDD)

**Precondition:** Phase 3 produced a concrete Python divergence with quantitative evidence.

### Task 4.1: Write the failing test

**Files:**
- Modify: `PhyloWGS_refactor/oversplit_test.go` (add new test group) OR create `PhyloWGS_refactor/round2_test.go`

**Step 1:** Write a test that pins the exact Python-equivalent behavior. Pattern: construct a small TSSB state deterministically, drive the affected function with a fixed seed, assert the spawn/retain outcome matches Python.

**Step 2:** Run: `go test -run <TestName> ./...`
Expected: FAIL with the exact divergence.

### Task 4.2: Apply the Python-fidelity fix

**Step 1:** Minimal code change in main.go that restores Python behavior. No regularization.

**Step 2:** Re-run the test. Expected: PASS.

**Step 3:** Run full suite: `go test ./...`. Expected: no regressions (especially `TestPrecomputedMHStatesEquivalence` and all 5 Round 1 oversplit tests).

**Step 4:** Commit: `fix: <specific divergence>`

### Task 4.3: Re-run the slice

**Step 1:** Rebuild: `cd PhyloWGS_refactor && go build -o /tmp/phylowgs-go-round2 .`

**Step 2:** Re-run `/tmp/run_slice.sh` after editing the `BIN=` line and `OUT_BASE=/.../results_round2_slice` / `LOG_BASE=/.../logs_round2_slice`.

**Step 3:** Parse results and build a comparison table against Report 5's numbers.

### Task 4.4: Write SIM_VALIDATION_REPORT_6

**Files:**
- Create: `PhyloWGS_refactor/sim_validation/SIM_VALIDATION_REPORT_6.md`

Use `SIM_VALIDATION_REPORT_5.md` as a template. Report honestly: if K finally drops toward true K, celebrate and invoke `superpowers:finishing-a-development-branch`. If not, document what remains and escalate.

### Task 4.5: Commit + finish branch

```bash
git add sim_validation/SIM_VALIDATION_REPORT_6.md
git commit -m "docs: SIM_VALIDATION_REPORT_6 - round 2 fix verified on slice"
```

Then, only if K dropped meaningfully on ≥3/6 fixtures, invoke `superpowers:finishing-a-development-branch`.

---

## Things to absolutely not do

- Do not introduce φ-pruning, post-MCMC tree simplification, `min_ssms` thresholds, subclone-size filters, or any other custom regularization that Python doesn't have. The user has been explicit about this three times now across investigations.
- Do not change the DP prior bounds (`min_dp_alpha=1.0`, `max_dp_alpha=50.0`, `min_dp_gamma=1.0`, `max_dp_gamma=10.0`, `min_alpha_decay=0.05`, `max_alpha_decay=0.80`).
- Do not touch `logLikelihoodWithCNVTreeMHPrecomputed`, `precomputeMHStates`, or `computeNGenomes` (the MH inner loop from `f969d36`).
- Do not leave temporary instrumentation (Task 2.3 counters, Task 3.2 Python prints) in committed code. Either gate everything behind `--trace` or revert before commit.
- Do not skip TDD in Phase 4. Iron Law from `superpowers:test-driven-development`.
- Do not commit without running `go test .` per `superpowers:verification-before-completion`.
- Do not modify the Python reference (`phylowgs/`) in a way that gets committed. Any diagnostic prints are throwaway.

## Decision points summary

| After phase | Decide |
|---|---|
| 1 | Scoring bug found → narrow fix task; otherwise continue |
| 2 | PATH X identified (find_node or stick_orders or cull) |
| 3 | Divergence found → Phase 4; no divergence → escalate |
| 4 | Success → finish branch; partial → Round 3 plan; regression → revert |
