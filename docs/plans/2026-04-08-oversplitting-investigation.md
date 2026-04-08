# Over-Splitting / K Inflation Investigation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use `superpowers:executing-plans` to implement this plan task-by-task. The first half of this plan is investigation/diagnostics — execute the audit and instrumentation tasks even when no code change is the obvious next step. Each phase ends with a decision checkpoint that determines the next phase's tasks.

**Goal:** Identify and fix the cause of over-splitting (K_inferred ≫ K_true) in the Go PhyloWGS port, while keeping Go's logic semantically equivalent to the upstream Python implementation. The bar is **match Python's tree growth dynamics** on the same fixture, not introduce custom regularization.

**Architecture:** Investigation-first. We don't yet know whether the divergence is in (a) the slice samplers for `dp_alpha`/`alpha_decay`/`dp_gamma`, (b) `resampleAssignments`, (c) `resampleStickOrders` (which handles splits), (d) `cullTree`, (e) tree initialization, or (f) the MCMC loop ordering. Each is plausible. We instrument both implementations, run a controlled comparison on a small fixture, find the divergence, then fix it via TDD.

**Tech Stack:** Go 1.x, Python 2.7 (PhyloWGS upstream, run via singularity image), Python 3 for analysis/plotting (matplotlib, json), bash for orchestration.

---

## Background — Read this first

You are taking over from a previous session. Here is the state of the world:

### What's been done

The Go port at `PhyloWGS_refactor/main.go` is a single-file Go translation of the upstream Python PhyloWGS at `phylowgs/`. Several rounds of bug fixes have already landed (see `sim_validation/SIM_VALIDATION_REPORT.md` through `SIM_VALIDATION_REPORT_4.md`). The most recent commit, `f969d36 perf(mh): precompute per-(SSM, node) (nr, nv) factors before MH loop`, made the CNV-MH inner loop into an O(K) dot product mirroring Python's `params.write_data_state` + `mh.cpp`. That landed a 12-78× speedup on CNV fixtures and Go is now 3-6× faster than Python end-to-end on CNV cases. **That optimization is complete and is not what this plan is about.**

### What's left — over-splitting

Both the Go port and upstream Python over-split clones on hard fixtures (K_inferred can be 5-10× the true K). Sometimes Go over-splits *more* than Python on the same fixture; sometimes they're equal. Examples from prior runs:

| fixture | true K | Go K_inferred | Python K_inferred |
|---|---:|---:|---:|
| K10_S1_T1000_M100_C0_rep1 | 10 | 84 | 81 |
| K3_S1_T200_M30_C5_rep0 | 3 | 7 | 4 |
| K5_S1_T200_M50_C2_rep0 | 5 | 20 | 2 |

The K5_S1_T200_M50_C2 case is interesting: Go over-splits to K=20 while Python under-splits to K=2 on the same fixture. That kind of asymmetry strongly suggests at least one of the resamplers diverges.

### What is in scope vs out of scope

**In scope:**
- Reading both Python and Go side-by-side and finding any divergence.
- Adding diagnostic instrumentation (iteration-level traces of `dp_alpha`, `dp_gamma`, `alpha_decay`, `num_nodes`, `llh`) to both implementations on a feature branch.
- Running a controlled comparison on a small, fast fixture.
- Fixing whatever divergence the comparison reveals, via TDD.
- Re-running the diagonal slice from `SIM_VALIDATION_REPORT_4.md` and confirming K_inferred drops on at least one fixture.

**Out of scope:**
- Inventing custom regularization (φ-pruning, minimum subclone size thresholds, etc.). These were already explicitly retracted in `sim_validation/INVESTIGATION_NextSteps.md` because Python doesn't do them and the user has been explicit: **the Go port must match Python's logic**.
- Tuning the DP prior bounds away from Python's defaults (`min_dp_alpha=1.0, max_dp_alpha=50.0, min_dp_gamma=1.0, max_dp_gamma=10.0, min_alpha_decay=0.05, max_alpha_decay=0.80`). If the bounds turn out to be wrong, the right fix is upstream-Python first.
- Touching the MH inner loop or `computeNGenomes` semantics. Those landed in `f969d36` and are now correct + fast.

### Known suspects

These came out of the previous session's investigation and code-reading. Phase 1 of this plan audits each one against Python:

1. **`resampleHypers` slice sampler (`main.go:2352-2447`).** Go's slice samplers cap the inner loop at 100 iterations (`for iter := 0; iter < 100; iter++`). Python uses `while True` (`tssb.py:261, 278, 305`) and explicitly raises if the sampler shrinks to zero. The cap means Go silently exits without updating the hyperparameter when the sampler doesn't accept in 100 tries — leaving the hyperparameter stuck at its previous value. This is likely benign in steady-state but could matter on early iterations or pathological trees.
2. **`resampleStickOrders` (`main.go:1814` onwards).** This is where new clones get spawned (the "stick-breaking" half of the resampler). Any divergence here directly inflates K. Especially the depth bookkeeping and the `findOrCreateNode` call pattern.
3. **`resampleAssignments` (`main.go:1248`).** Where each datum is reassigned to a node sampled from the conditional posterior. Bugs here can push data into spurious nodes.
4. **`cullTree` (`main.go:1736`).** Currently only trims trailing empty children — matches Python's `cull_tree` (`tssb.py:152-166`). But the bookkeeping around `kill_node` could differ.
5. **`newTSSB` initialization (`main.go:569`).** Both implementations start with all data assigned to a single child of the root, but the random `pi` allocation in `newTSSB` may not match `evolve.py:71-95` exactly.
6. **`runChain` MCMC loop ordering (`main.go:2453-2583`).** Verified in the previous session to match `evolve.py:171-211` (assignments → cull → MH → sticks → stickOrders → hypers → llh). Probably fine but worth re-checking after Phase 1 lands traces.

### File map

You will be reading from and (later) writing to these files:

| File | Purpose |
|---|---|
| `PhyloWGS_refactor/main.go` | The entire Go port. |
| `PhyloWGS_refactor/mh_states_test.go` | Existing TDD test file (the only one). New tests for this plan go here or in a new `oversplit_test.go`. |
| `PhyloWGS_refactor/sim_validation/` | Scripts and reports from the validation sweep. |
| `PhyloWGS_refactor/sim_validation/SIM_VALIDATION_REPORT_4.md` | Most recent slice report — defines the diagonal slice you'll re-validate against. |
| `PhyloWGS_refactor/sim_validation/INVESTIGATION_NextSteps.md` | Previous investigation document. Section "Why φ-pruning is wrong" is critical context. |
| `phylowgs/tssb.py` | Upstream Python TSSB implementation. The reference. |
| `phylowgs/evolve.py` | Upstream Python MCMC driver loop. |
| `phylowgs/mh.cpp`, `phylowgs/mh.hpp`, `phylowgs/params.py` | Compiled MH inner loop and per-SSM precomputation. **You should not need to modify these.** |
| `simulation_validation_updated/fixtures/` | Simulated fixtures used for validation. |
| `simulation_validation_updated/results_new_slice/` | Output of the diagonal slice from `SIM_VALIDATION_REPORT_4.md`. |
| `simulation_validation_updated/analysis/go-cpu/scores.json` | Pre-optimization Go grid scores. |
| `simulation_validation_updated/analysis/original-python/scores.json` | Python grid scores from the same default-grid run. |

### Environment setup commands

```bash
# Go binary location used in slice runs
ls -la /tmp/phylowgs-go-new

# Go compiler
export PATH=/tmp/go/bin:$PATH
which go && go version

# Python 2 reference is run via singularity. Check whether the SIF is available:
ls $NXF_SINGULARITY_CACHEDIR/*phylowgs* 2>&1 || echo "no SIF cached"
# If no SIF, you'll fall back to comparing against the recorded
# original-python scores in simulation_validation_updated/analysis/original-python/

# Python 3 for analysis
python3 -c 'import json, matplotlib; print("ok")'

# Working dir (use absolute paths in commands)
cd /sessions/dreamy-pensive-knuth/mnt/phylowgs_cowork/PhyloWGS_refactor
```

### Skills to use

- `superpowers:systematic-debugging` — for Phase 1 audit and Phase 3 hypothesis verification
- `superpowers:test-driven-development` — for any code change in Phase 5
- `superpowers:verification-before-completion` — before claiming any phase complete
- `superpowers:finishing-a-development-branch` — at the very end

---

## Phase 1 — Side-by-side audit of suspects (no code changes)

**Goal of phase:** For each suspect listed in the background, read Go and Python in parallel and write down every concrete difference. No edits. Output is a text artifact: `sim_validation/AUDIT_oversplitting.md`.

### Task 1.1: Create the audit document

**Files:**
- Create: `PhyloWGS_refactor/sim_validation/AUDIT_oversplitting.md`

**Step 1: Create with template**

Write the file with this skeleton:

```markdown
# Over-Splitting Audit — Go vs Python Side-by-Side

**Date:** <today>
**Branch:** <git rev-parse --abbrev-ref HEAD>
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

(filled in as Phase 1 progresses)
```

**Step 2: Commit**

```bash
git add PhyloWGS_refactor/sim_validation/AUDIT_oversplitting.md
git commit -m "docs: scaffold over-splitting audit document"
```

### Task 1.2: Audit `resampleHypers` slice samplers

**Files:**
- Read: `PhyloWGS_refactor/main.go:2352-2447`
- Read: `phylowgs/tssb.py:245-316`
- Append to: `PhyloWGS_refactor/sim_validation/AUDIT_oversplitting.md`

**Step 1: Read both functions in full**

Open both ranges. Compare:
- Bounds (`min_dp_alpha`, `max_dp_alpha`, `min_dp_gamma`, `max_dp_gamma`, `min_alpha_decay`, `max_alpha_decay`).
- Loop structure (Python `while True` vs Go `for iter := 0; iter < 100; iter++`).
- Slice level computation (`llh_s = log(rand()) + dp_alpha_llh(...)`).
- Order of resampling (alpha, then decay, then gamma).
- How `alpha_decay` is passed into the `dp_alpha_llh` closure (Python passes both args explicitly, Go closure captures `t.AlphaDecay` from the receiver and the `alpha_decay` slice sampler has to mutate-test-restore).
- The boundary condition `elif new_dp_alpha == self.dp_alpha: raise Exception("Slice sampler shrank to zero!")` in Python — does Go handle this case?

**Step 2: Document each difference in `AUDIT_oversplitting.md`**

For each difference, write:
- One paragraph describing it precisely (with line numbers in both files).
- A severity rating.
- Whether to fix it.

**Known seed observation to capture:**

> Go's slice samplers cap the inner loop at 100 iterations (`main.go:2376, 2396, 2434`). Python uses `while True` (`tssb.py:261, 278, 305`). When the slice sampler doesn't accept within 100 tries, Go exits silently with the hyperparameter unchanged, while Python would either eventually accept (if the LLH is well-behaved) or raise an exception. **Severity: SUSPICIOUS — could leave dp_alpha pinned high if the early-iteration LLH is rough.**

**Step 3: Commit**

```bash
git add PhyloWGS_refactor/sim_validation/AUDIT_oversplitting.md
git commit -m "audit: resampleHypers slice samplers Go vs Python"
```

### Task 1.3: Audit `resampleStickOrders`

**Files:**
- Read: `PhyloWGS_refactor/main.go:1814` to the end of the function
- Read: `phylowgs/tssb.py` — find `def resample_stick_orders` (use grep)
- Append to: `PhyloWGS_refactor/sim_validation/AUDIT_oversplitting.md`

**Step 1: Locate Python function**

```bash
grep -n "def resample_stick_orders" /sessions/dreamy-pensive-knuth/mnt/phylowgs_cowork/phylowgs/tssb.py
```

Read the entire function in both languages.

**Step 2: Compare**

Especially focus on:
- The branch where a child is *spawned* (look for code that creates a new node — this is where K grows).
- Reordering logic (which children move where in the stick order).
- Whether empty children are pruned during reordering or only by `cull_tree`.
- The recursion termination (max depth, min depth).
- Whether `findOrCreateNode` is called and with what `u` parameter — `u` controls how deeply the new node sits.

**Step 3: Document differences with severity ratings, then commit**

```bash
git add PhyloWGS_refactor/sim_validation/AUDIT_oversplitting.md
git commit -m "audit: resampleStickOrders Go vs Python"
```

### Task 1.4: Audit `resampleAssignments`

**Files:**
- Read: `PhyloWGS_refactor/main.go:1248-1404`
- Read: `phylowgs/tssb.py` — find `def resample_assignments`
- Append to: `PhyloWGS_refactor/sim_validation/AUDIT_oversplitting.md`

**Step 1: Read both**

Compare:
- The conditional posterior over node assignments (the `path` u-sampling).
- How CNV-aware likelihoods are wired in (Go's `logLikelihoodWithCNVTree` vs Python's path).
- Slice sampler shrinkage on `u` (look for similar `while True` vs capped loops).
- The "spawn new node if u falls in unallocated mass" branch — this is the second place where K can grow.
- Whether the data point is removed from its old node before the slice sample, and whether the new node is committed atomically.

**Step 2: Document, commit**

```bash
git add PhyloWGS_refactor/sim_validation/AUDIT_oversplitting.md
git commit -m "audit: resampleAssignments Go vs Python"
```

### Task 1.5: Audit `cullTree` and `killNode`

**Files:**
- Read: `PhyloWGS_refactor/main.go:1722-1775`
- Read: `phylowgs/tssb.py:152-200` (grep `def cull_tree`, `def kill`)
- Append to: `AUDIT_oversplitting.md`

**Step 1: Read both**

This was previously believed to match Python (Report 3 confirmed that both only trim *trailing* empty children, not arbitrary internal empty nodes). Confirm one more time. Also check:
- What `killNode` removes from (`Data` slice, `Children` slice, `Sticks` slice).
- Whether stick weights are renormalized after a cull.
- Whether the parent's `Pi` is updated to reclaim mass.

**Step 2: Document, commit**

### Task 1.6: Audit tree initialization (`newTSSB` vs `evolve.py` startup)

**Files:**
- Read: `PhyloWGS_refactor/main.go:569-641`
- Read: `phylowgs/evolve.py:60-110`
- Append to: `AUDIT_oversplitting.md`

**Step 1: Read both**

Verify:
- Initial child count (Python: 1; Go: 1, set in `newTSSB`).
- Initial `pi` allocation (Python uses `boundbeta(1, dp_gamma)` for sticks and `boundbeta(1.0, alpha_decay*dp_alpha)` for `main`; Go uses `rng.Float64()` for the initial child's pi fraction — **this is a possible divergence**).
- The Dirichlet vs uniform choice for initial node parameters.

**Step 2: Document — flag the rng.Float64() vs boundbeta divergence as SUSPICIOUS**

**Step 3: Commit**

### Task 1.7: Re-confirm MCMC loop ordering

**Files:**
- Read: `PhyloWGS_refactor/main.go:2519-2573` (the iter loop in `runChain`)
- Read: `phylowgs/evolve.py:163-220`

**Step 1: Diff the call sequences**

Both should be: `resample_assignments → cull_tree → set_node_height → set_path_from_root_to_node → map_datum_to_node → metropolis → resample_sticks → resample_stick_orders → resample_hypers → complete_data_log_likelihood`. Confirm. Note any extra calls in either implementation (e.g., the new precompute call we just landed in Go).

**Step 2: Document, commit**

### Task 1.8: Phase 1 decision checkpoint

**Files:**
- Edit: `PhyloWGS_refactor/sim_validation/AUDIT_oversplitting.md` — add a "Decision" section at the bottom.

**Step 1: Re-read every audit entry**

Tally up: how many BLOCKERs, how many SUSPICIOUS, how many BENIGN.

**Step 2: Write the Decision section**

If there is even one BLOCKER, the next step is "Phase 5: TDD fix" — skip Phases 2-4.

If there are only SUSPICIOUS items, the next step is "Phase 2: Diagnostic instrumentation" — we need empirical evidence to rank them.

If everything is BENIGN, this is more serious than expected and the fix is probably in the prior bounds themselves. In that case, escalate by re-reading the Python paper's appendix on the DP prior and asking the user.

**Step 3: Commit**

```bash
git add PhyloWGS_refactor/sim_validation/AUDIT_oversplitting.md
git commit -m "audit: phase 1 complete — decision checkpoint"
```

---

## Phase 2 — Diagnostic instrumentation (only if Phase 1 found no BLOCKERs)

**Goal of phase:** Make both Go and Python emit per-iteration traces of `(dp_alpha, dp_gamma, alpha_decay, num_nodes, num_data_in_root, llh)` to a JSONL file. Don't change any sampling logic. The traces let us *see* the divergence empirically.

### Task 2.1: Add Go trace output behind a flag

**Files:**
- Modify: `PhyloWGS_refactor/main.go` — add `--trace <path>` flag in `main()`, plumb a trace writer through `runChain`.
- Test: `PhyloWGS_refactor/oversplit_test.go` — new test verifying trace JSONL is well-formed and one record is emitted per iteration.

**Step 1: Write the failing test**

```go
package main

import (
    "encoding/json"
    "math/rand"
    "os"
    "strings"
    "testing"
)

func TestTraceWriterEmitsOnePerIteration(t *testing.T) {
    var buf strings.Builder
    tw := newTraceWriter(&buf)

    rng := rand.New(rand.NewSource(7))
    ssms, cnvs := makeTinyDataset(t, rng)
    tssb := newTSSB(ssms, cnvs, 25.0, 1.0, 0.25, rng)

    // Drive 5 iters of fake "after-iter" emit
    for i := 0; i < 5; i++ {
        tw.emit(i, tssb, -1234.5)
    }
    tw.close()

    lines := strings.Split(strings.TrimRight(buf.String(), "\n"), "\n")
    if len(lines) != 5 {
        t.Fatalf("want 5 trace lines, got %d", len(lines))
    }
    for i, ln := range lines {
        var rec map[string]interface{}
        if err := json.Unmarshal([]byte(ln), &rec); err != nil {
            t.Fatalf("line %d not valid json: %v", i, err)
        }
        for _, k := range []string{"iter", "dp_alpha", "dp_gamma", "alpha_decay",
            "num_nodes", "llh"} {
            if _, ok := rec[k]; !ok {
                t.Errorf("line %d missing key %q", i, k)
            }
        }
    }
}

// makeTinyDataset is defined in mh_states_test.go style — reuse buildTestTSSB
// helper if available, or write a 3-SSM no-CNV version here.
```

**Step 2: Run test, verify it fails**

```bash
cd PhyloWGS_refactor && go test -run TestTraceWriterEmitsOnePerIteration .
```
Expected: FAIL with `newTraceWriter undefined` (or similar).

**Step 3: Implement minimal `traceWriter`**

In `main.go` near other utility types, add:

```go
type traceWriter struct {
    w   io.Writer
    enc *json.Encoder
}

func newTraceWriter(w io.Writer) *traceWriter {
    return &traceWriter{w: w, enc: json.NewEncoder(w)}
}

func (tw *traceWriter) emit(iter int, t *TSSB, llh float64) {
    if tw == nil {
        return
    }
    nodes := t.getNodes()
    rec := map[string]interface{}{
        "iter":        iter,
        "dp_alpha":    t.DPAlpha,
        "dp_gamma":    t.DPGamma,
        "alpha_decay": t.AlphaDecay,
        "num_nodes":   len(nodes),
        "llh":         llh,
    }
    _ = tw.enc.Encode(rec)
}

func (tw *traceWriter) close() { /* no-op for now */ }
```

(Add `"io"` to imports if needed.)

**Step 4: Run test, verify it passes**

```bash
go test -run TestTraceWriterEmitsOnePerIteration .
```

**Step 5: Wire into runChain**

Add a `traceWriter *traceWriter` parameter to `runChain`. Inside the iter loop, after `llh := tssb.completeDataLogLikelihood()`, call `traceWriter.emit(iter, tssb, llh)`.

In `main()`, add a `--trace` flag, open the file, wrap in a buffered writer, pass into `runChain`, close in a `defer`.

**Step 6: Smoke test the wiring**

```bash
go build -o /tmp/phylowgs-go-trace .
mkdir -p /tmp/trace-out
/tmp/phylowgs-go-trace --no-gpu -B 20 -s 50 -j 1 -i 500 -r 42 \
    --trace /tmp/trace-out/trace.jsonl \
    -O /tmp/trace-out \
    /sessions/dreamy-pensive-knuth/mnt/phylowgs_cowork/simulation_validation_updated/fixtures/K3_S1_T200_M30_C0_rep0/ssm_data.txt \
    /sessions/dreamy-pensive-knuth/mnt/phylowgs_cowork/simulation_validation_updated/fixtures/K3_S1_T200_M30_C0_rep0/cnv_data.txt
wc -l /tmp/trace-out/trace.jsonl
head -3 /tmp/trace-out/trace.jsonl
```
Expected: 70 lines (20 burnin + 50 sample), each a JSON object with the keys above.

**Step 7: Commit**

```bash
git add PhyloWGS_refactor/main.go PhyloWGS_refactor/oversplit_test.go
git commit -m "feat: --trace flag emits per-iteration hyperparameter trace JSONL"
```

### Task 2.2: Add Python trace output

**Files:**
- Modify: `phylowgs/evolve.py` — add a similar `--trace <path>` argument.

**Step 1: Find argparse / sys.argv handling in evolve.py**

```bash
grep -n "argparse\|sys.argv" /sessions/dreamy-pensive-knuth/mnt/phylowgs_cowork/phylowgs/evolve.py | head
```

**Step 2: Add a `--trace` flag**

After the existing argparse setup, add:
```python
parser.add_argument('--trace', dest='trace_path', default=None,
                    help='Write per-iteration JSONL trace of hyperparameters and num_nodes')
```

In the iteration loop (around `evolve.py:213`), after `last_llh = tssb.complete_data_log_likelihood()`, add:
```python
if state.get('trace_path'):
    import json
    with open(state['trace_path'], 'a') as f:
        weights, nodes = tssb.get_mixture()
        rec = {
            'iter': iteration,
            'dp_alpha': float(tssb.dp_alpha),
            'dp_gamma': float(tssb.dp_gamma),
            'alpha_decay': float(tssb.alpha_decay),
            'num_nodes': int(len(nodes)),
            'llh': float(last_llh),
        }
        f.write(json.dumps(rec) + '\n')
```

Plumb `state['trace_path'] = args.trace_path` near where other args are stored.

**Step 3: Smoke test through singularity (if SIF available) or skip if not**

If the SIF is available, run a quick test the same way as Task 2.1's smoke test. If not, defer Python trace generation to when the user has access to the HPC environment and proceed to Task 2.3 with Go-only data.

**Step 4: Commit (in the phylowgs/ submodule or as a patch in PhyloWGS_refactor/)**

If `phylowgs/` is a git submodule, you'll need to commit there separately. If it isn't, save the patch to `PhyloWGS_refactor/sim_validation/python_trace.patch` and commit that:

```bash
cd /sessions/dreamy-pensive-knuth/mnt/phylowgs_cowork
git diff phylowgs/evolve.py > PhyloWGS_refactor/sim_validation/python_trace.patch
cd PhyloWGS_refactor
git add sim_validation/python_trace.patch
git commit -m "patch: add --trace flag to Python evolve.py for diagnostic comparison"
```

---

## Phase 3 — Comparative diagnostic run

**Goal of phase:** Run both implementations on the same small fixture with traces enabled, then plot/inspect the trajectories.

### Task 3.1: Pick the diagnostic fixture

Use `K3_S1_T200_M30_C0_rep0` (no CNVs, 30 SSMs, true K=3). It's the smallest fixture in the slice, runs in seconds, and the previous slice run already showed Go infers K=30 here — extreme over-splitting.

### Task 3.2: Run Go with trace enabled

**Step 1: Run**

```bash
mkdir -p /tmp/diag/go
/tmp/phylowgs-go-trace --no-gpu -B 200 -s 500 -j 1 -i 2000 -r 42 \
    --trace /tmp/diag/go/trace.jsonl \
    -O /tmp/diag/go \
    /sessions/dreamy-pensive-knuth/mnt/phylowgs_cowork/simulation_validation_updated/fixtures/K3_S1_T200_M30_C0_rep0/ssm_data.txt \
    /sessions/dreamy-pensive-knuth/mnt/phylowgs_cowork/simulation_validation_updated/fixtures/K3_S1_T200_M30_C0_rep0/cnv_data.txt
wc -l /tmp/diag/go/trace.jsonl
```
Expected: 700 lines (200 burnin + 500 sample). One chain only (`-j 1`) so the trace is not interleaved.

**Step 2: Inspect head and tail**

```bash
head -5 /tmp/diag/go/trace.jsonl
tail -5 /tmp/diag/go/trace.jsonl
```

Note the `num_nodes` trajectory. Does it grow monotonically? Does `dp_alpha` get pinned at 50 (the max)?

### Task 3.3: Run Python with trace enabled (if SIF available)

**Step 1: Run inside singularity**

```bash
SIF=$(ls $NXF_SINGULARITY_CACHEDIR/*phylowgs* 2>/dev/null | head -1)
if [[ -n "$SIF" ]]; then
    mkdir -p /tmp/diag/py
    singularity exec --bind /sessions/dreamy-pensive-knuth/mnt/phylowgs_cowork "$SIF" \
        python2 /usr/bin/phylowgs/evolve.py \
        -B 200 -s 500 \
        --trace /tmp/diag/py/trace.jsonl \
        /sessions/dreamy-pensive-knuth/mnt/phylowgs_cowork/simulation_validation_updated/fixtures/K3_S1_T200_M30_C0_rep0/ssm_data.txt \
        /sessions/dreamy-pensive-knuth/mnt/phylowgs_cowork/simulation_validation_updated/fixtures/K3_S1_T200_M30_C0_rep0/cnv_data.txt
    wc -l /tmp/diag/py/trace.jsonl
else
    echo "No SIF — falling back to recorded analysis/original-python data"
fi
```

If no SIF, you can still proceed to Task 3.4 with Go-only data plus the recorded inferred_K from `simulation_validation_updated/analysis/original-python/scores.json`. Note the limitation in your write-up.

### Task 3.4: Plot trajectories side by side

**Files:**
- Create: `PhyloWGS_refactor/sim_validation/plot_traces.py`

**Step 1: Write the plot script**

```python
#!/usr/bin/env python3
"""Plot Go vs Python iteration-level traces from --trace JSONL output."""
import json, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def load(path):
    with open(path) as f:
        return [json.loads(line) for line in f if line.strip()]

go = load('/tmp/diag/go/trace.jsonl')
try:
    py = load('/tmp/diag/py/trace.jsonl')
except FileNotFoundError:
    py = None

fig, axes = plt.subplots(5, 1, figsize=(10, 12), sharex=True)
for ax, key, label in zip(axes,
        ['num_nodes', 'dp_alpha', 'dp_gamma', 'alpha_decay', 'llh'],
        ['# nodes', 'dp_alpha', 'dp_gamma', 'alpha_decay', 'log-likelihood']):
    ax.plot([r['iter'] for r in go], [r[key] for r in go], label='Go', color='C0')
    if py:
        ax.plot([r['iter'] for r in py], [r[key] for r in py], label='Python', color='C1')
    ax.set_ylabel(label)
    ax.legend(loc='best')
axes[-1].set_xlabel('iteration')
plt.tight_layout()
plt.savefig('/tmp/diag/traces.png', dpi=120)
print('wrote /tmp/diag/traces.png')
```

**Step 2: Run**

```bash
python3 PhyloWGS_refactor/sim_validation/plot_traces.py
```

**Step 3: Look at the plot**

The user will need to view this manually (it's in `/tmp/diag/traces.png`). Save a copy to the workspace folder so the user can see it:

```bash
cp /tmp/diag/traces.png /sessions/dreamy-pensive-knuth/mnt/phylowgs_cowork/simulation_validation_updated/diag_traces_K3_M30_C0.png
```

**Step 4: Write findings**

Append to `sim_validation/AUDIT_oversplitting.md` a "Phase 3 findings" section. Specifically answer:
- Does Go's `num_nodes` grow without bound while Python's stabilizes? Or do both grow?
- Where does Go's `dp_alpha` settle? At the upper bound (50)? Mid-range? Compared to Python?
- Is the LLH still improving when over-splitting starts, or does the model start fitting noise?

**Step 5: Commit**

```bash
git add PhyloWGS_refactor/sim_validation/plot_traces.py PhyloWGS_refactor/sim_validation/AUDIT_oversplitting.md
git commit -m "diag: comparative trace plots for over-splitting investigation"
```

### Task 3.5: Phase 3 decision checkpoint

Re-read the Phase 1 audit and the Phase 3 findings together. The combination should now point to a single most-likely cause:

| Pattern | Likely cause |
|---|---|
| Go's `dp_alpha` pinned at 50, Python's free | Slice sampler bug (likely the 100-iter cap) |
| Go's `num_nodes` grows but Python's doesn't, hypers similar | Bug in `resampleStickOrders` or `resampleAssignments` |
| Both `num_nodes` grow similarly | The bug is in the prior/likelihood balance — likely in `paramPost` or `completeDataLogLikelihood` |
| Hypers diverge but `num_nodes` doesn't | Cosmetic — defer |

Write the chosen hypothesis as a one-paragraph "Hypothesis" section in `AUDIT_oversplitting.md`. Commit.

---

## Phase 4 — Confirm hypothesis with a targeted experiment

**Goal of phase:** Before changing production code, run a single targeted experiment that would falsify the hypothesis. If the hypothesis survives, proceed to Phase 5. If not, return to Phase 1 with new information.

**This is intentionally vague** because the experiment depends on the hypothesis. Examples:

- **Hypothesis: 100-iter cap is the bug.** Experiment: temporarily change Go's slice samplers to `for iter := 0; iter < 10000; iter++`, rerun the diagnostic fixture, see if `num_nodes` trajectory matches Python better.
- **Hypothesis: `resampleStickOrders` spawns extra children.** Experiment: add a counter that increments every time a node is spawned. Compare counts between Go and (instrumented) Python.
- **Hypothesis: initial pi allocation is wrong.** Experiment: replace `rng.Float64()` in `newTSSB` with `boundBeta(1, dpGamma)`, rerun, see if early-iteration `num_nodes` differs.

Do the experiment on a throwaway branch (`git checkout -b experiment/<hypothesis>`). Don't merge. Just document the result in `AUDIT_oversplitting.md`.

If the experiment **confirms** the hypothesis: discard the experimental branch, proceed to Phase 5 with TDD.
If the experiment **refutes** it: discard the experimental branch, return to Phase 1 with new clues.

---

## Phase 5 — TDD fix (only after Phase 4 confirms)

**Goal of phase:** Land the fix as a reviewed, tested commit on `go-port` (not on the experimental branch).

### Task 5.1: Write the failing test

The test must fail without the fix and pass with it. Use the diagnostic fixture data or a synthetic minimal case. Examples:

- For the 100-iter slice cap fix: a unit test that constructs a TSSB with hand-set hyperparameters in a region where the slice sampler needs >100 iters to accept, calls `resampleHypers`, and asserts the hyperparameter changed.
- For a `resampleStickOrders` spawn-bug fix: a property-based test that runs N iterations on a fixed-seed tree and asserts `num_nodes` stays below some bound that the buggy version exceeds.

Add the test to `oversplit_test.go` (created in Phase 2) or `mh_states_test.go`.

### Task 5.2: Run the test, verify it fails

```bash
cd PhyloWGS_refactor && go test -run TestOverSplitting <details> .
```
Expected: FAIL.

### Task 5.3: Implement the minimal fix

Edit `main.go`. Keep the change as small as possible — this is a fidelity fix, not a refactor.

### Task 5.4: Run the test, verify it passes

```bash
go test -run TestOverSplitting <details> .
```
Expected: PASS.

### Task 5.5: Run the full test suite

```bash
go test .
```
Expected: all tests pass, including the existing `TestPrecomputedMHStatesEquivalence`.

### Task 5.6: Commit

```bash
git add PhyloWGS_refactor/main.go PhyloWGS_refactor/oversplit_test.go
git commit -m "fix(<area>): <one-line summary of what diverged from Python>

<2-3 sentence explanation of root cause>

<2-3 sentence explanation of the fix>

Test: TestOverSplitting<...> covers the divergence by <how>."
```

---

## Phase 6 — Validate against the slice

**Goal of phase:** Re-run the 6-fixture diagonal slice from `SIM_VALIDATION_REPORT_4.md` with the fixed binary and compare K_inferred and runtime to the recorded baselines.

### Task 6.1: Rebuild the binary

```bash
cd PhyloWGS_refactor && go build -o /tmp/phylowgs-go-fixed .
```

### Task 6.2: Re-run the slice

Reuse the runner script from Report 4's session — it lives at `/tmp/run_slice.sh` in the previous session and may not exist in this one. If not, recreate it:

```bash
cat > /tmp/run_slice.sh <<'EOF'
#!/usr/bin/env bash
set -u
BIN=/tmp/phylowgs-go-fixed
FX_DIR=/sessions/dreamy-pensive-knuth/mnt/phylowgs_cowork/simulation_validation_updated/fixtures
OUT_BASE=/sessions/dreamy-pensive-knuth/mnt/phylowgs_cowork/simulation_validation_updated/results_phase6_slice
LOG_BASE=/sessions/dreamy-pensive-knuth/mnt/phylowgs_cowork/simulation_validation_updated/logs_phase6_slice
mkdir -p "$OUT_BASE" "$LOG_BASE"
FIXTURES=(
  K3_S1_T200_M30_C2_rep0
  K3_S1_T200_M30_C5_rep0
  K5_S1_T200_M50_C2_rep0
  K5_S1_T200_M50_C5_rep0
  K10_S1_T1000_M100_C2_rep0
  K10_S1_T1000_M100_C5_rep0
)
for fx in "${FIXTURES[@]}"; do
  out="$OUT_BASE/$fx"
  log="$LOG_BASE/$fx.log"
  rm -rf "$out"; mkdir -p "$out"
  echo "=== $fx $(date +%T) ===" | tee -a "$LOG_BASE/_overall.log"
  start=$(date +%s)
  "$BIN" --no-gpu -B 500 -s 1000 -j 4 -r 42 -O "$out" \
      "$FX_DIR/$fx/ssm_data.txt" "$FX_DIR/$fx/cnv_data.txt" > "$log" 2>&1
  rc=$?
  end=$(date +%s)
  echo "  exit=$rc elapsed=$((end-start))s" | tee -a "$LOG_BASE/_overall.log"
  grep -E "(Total time|Best LLH|Median LLH)" "$log" | tee -a "$LOG_BASE/_overall.log"
done
EOF
chmod +x /tmp/run_slice.sh
nohup /tmp/run_slice.sh > /tmp/run_slice_stdout.log 2>&1 &
echo "pid=$!"
```

This will take ~25-90 minutes depending on whether the fix changes per-iteration cost.

### Task 6.3: Parse and compare

Reuse the parsing script style from Report 4 (`SIM_VALIDATION_REPORT_4.md`). Build a table:

```
fixture                      true_K  prev_K  new_K  prev_s  new_s  best_llh
```

Where `prev_*` is from `simulation_validation_updated/results_new_slice/_slice_summary.json` (the post-MH-precompute slice), and `new_*` is from this Phase 6 run.

The success criterion is **at least 3 of 6 fixtures show lower K_inferred without a major LLH regression** (LLH may move ±10 nats due to seed; that's noise).

### Task 6.4: Write SIM_VALIDATION_REPORT_5

**Files:**
- Create: `PhyloWGS_refactor/sim_validation/SIM_VALIDATION_REPORT_5.md`

Use the structure of `SIM_VALIDATION_REPORT_4.md` as a template. Sections: change-under-test, settings, results table, interpretation, correctness evidence, what changed in the code, what's still left.

### Task 6.5: Commit + finish branch

```bash
git add PhyloWGS_refactor/sim_validation/SIM_VALIDATION_REPORT_5.md
git commit -m "docs: SIM_VALIDATION_REPORT_5 — over-splitting fix verified on slice"
```

Then invoke the `superpowers:finishing-a-development-branch` skill to choose merge / PR / keep-as-is.

---

## Decision points summary

| After phase | Decide |
|---|---|
| 1 | Skip to Phase 5 (BLOCKER found) or proceed to Phase 2 |
| 3 | Pick a single hypothesis, optionally collapse Phase 4 if confidence is high |
| 4 | Confirmed → Phase 5; Refuted → Phase 1 with new clues |
| 6 | Success → finish branch; partial success → next plan; regression → revert |

## Things to absolutely not do

- Do not introduce φ-pruning, post-MCMC tree simplification, minimum subclone size thresholds, or any other custom regularization that Python doesn't have. The user has been explicit about this twice.
- Do not change the DP prior bounds (`min_dp_alpha=1.0`, `max_dp_alpha=50.0`, `min_dp_gamma=1.0`, `max_dp_gamma=10.0`, `min_alpha_decay=0.05`, `max_alpha_decay=0.80`). These match Python and should remain matched. If they appear to be the actual root cause, escalate to the user before changing them.
- Do not touch the MH inner loop (`logLikelihoodWithCNVTreeMHPrecomputed`, `precomputeMHStates`, `computeNGenomes`). That work landed in `f969d36` and is verified.
- Do not skip TDD on the fix in Phase 5. The Iron Law from `superpowers:test-driven-development` is non-negotiable.
- Do not commit without running `go test .` and seeing it pass first (`superpowers:verification-before-completion`).
