# Round 2 Phase 3 Status — K-inflation premise falsified

**Date:** 2026-04-08
**Branch:** `go-port` @ `7faba3d`
**Plan being executed:** `docs/plans/2026-04-09-oversplitting-round-2.md`

## TL;DR

The Round 2 plan was built on the premise that the Go port over-splits
relative to Python (Round 1 / Report 5 reported "Go produces K≈30 vs
true K=3"). Phase 3 investigation **falsified that premise** on the K3
fixture used as the canonical comparison point. Go and Python produce
essentially identical K distributions on `K3_S1_T200_M30_C2_rep0`. The
remaining K=30 number from Report 5 appears to come from the
*reporting-layer* `removeEmptyNodes` flattening, not from the sampler.

Phase 4 (TDD fix in `findOrCreateNode`) is **paused** pending broader
HPC validation across the slice and a user decision on direction.

## What was done in this session

### Phase 1 — instrumentation

Added a `--trace` flag that writes per-iter NDJSON diagnostics
(`K_before_cull`, `K_after_cull`, spawn counts by site, kill counts,
hyperparameters, llh) for one or all chains. Backed by tests in
`trace_test.go`. See commits `5fb3249`, `c3fb58b`.

### Phase 2 — first audit pass

Ran K3_C2 short chain with the new tracer. Found 79% of all spawns
happen inside `findOrCreateNode` (vs `resampleStickOrders`). Two
suspicious divergences flagged for Phase 3 investigation:

1. The Go-only `maxChildCreations = 20` safeguard
2. The depth==0 special case in `findOrCreateNode`

Documented in `sim_validation/AUDIT_round2.md`.

### Phase 3 — both hypotheses falsified

**Falsification 1 — the cap.** Added a `CapHitsFindNode` diagnostic
counter that increments every time the safeguard fires *before* `u`
falls into the stick space. Re-ran K3_C2 short (150 iters) and long
(1000 iters) chains. Result across 1150 total iters: `cap_hits = 0`.
The cap never fires on this workload. It cannot be the source of
K-inflation. Hypothesis rejected.

**Falsification 2 — the depth==0 branch.** Re-read Python's
`tssb.py:find_node` (lines 340-377). Python also hardcodes `index = 0`
at depth 0 and descends into `root['children'][0]`. Python's `find_node`
*does not* create new root children either. The Go branch at
`main.go:2409-2420` is structurally equivalent. Hypothesis rejected.

**Full structural diff.** Line-by-line comparison of `findOrCreateNode`
against `find_node` confirmed equivalence in:

- depth-0 index hard-coding
- `u` rescaling on remaining stick mass at each depth
- lazy `spawnChild` append when `u > edges[-1]`
- the `boundbeta(1, alpha_decay^depth * dp_alpha)` stick draw
- absence of min_depth handling in either (both default to 0)
- absence of gamma-prior or DP-prior sampling inside this function

The only structural divergence remaining is the now-proven-dead cap.

### Head-to-head on K3_S1_T200_M30_C2_rep0

```
                       Python (orig)        Go (current head)
N trees:               1000                 4000 (4 chains × 1000)
range (K):             4 – 16               3 – 15
median K:              13                   12
mean K:                11.95                11.42
mode K:                13 (279/1000)        12 (1228/4000)
best-llh tree K:       5                    4
best-llh value:        −315.8               −307.85   ← Go better
```

Go and Python produce essentially the same K distribution. Go is
marginally **tighter** (median 12 vs 13) and finds a slightly better
likelihood. Both samplers over-split the posterior relative to the true
K=3 — because PhyloWGS's DP prior with `dp_alpha ∈ [1, 50]` admits
10+ populations on 30-SSM data. **This is a model property, not a
fidelity bug in either implementation.**

### Where the original "K=30" claim came from

The Round 1 / Report 5 K=30 observation used the post-`removeEmptyNodes`
flat structure. That post-processing function reparents deep
intermediate-empty nodes to the root, producing an artificial leaf
fan-out. Report 5 was measuring *reporting-layer* flatness, not sampler
fanout. The sampler itself carries a hierarchical K ≈ 12, matching
Python.

## Round 1 fidelity fixes (still valid)

The Round 1 fixes from `1ed4e3c` are independently correct and remain
in place:

- Two BLOCKER fixes in `findOrCreateNode` / spawn flow
- The pi broadcast scalar bug in `spawnChild`
- Unconditional root term in `dpAlphaLLH`
- Hyperparameter slice-sampler bounds matching Python
  (`dp_alpha [1,50]`, `alpha_decay [0.05, 0.80]`, `dp_gamma [1, 10]`)

These were necessary for the Go port to reach the K=12 ballpark in the
first place. Without them, the gap was much wider.

## Next steps

Two parallel tracks:

1. **HPC validation slice (in progress).** John kicked off
   `submit_sim_validation.sh` on `go-port @ 7faba3d` against the full
   K3 / K5 / K8 fixture set with default parameters (B=500, S=1000,
   4 chains). When results arrive: per-fixture head-to-head Go vs
   Python K distribution comparison.

2. **Local trace baseline (in progress).** Run a small set of fixtures
   locally with `--trace` enabled to build a per-iter spawn-dynamics
   reference *before* the HPC results land. This gives us calibration
   data for distinguishing "normal burnin slice-walk dynamics" from a
   real divergence when interpreting tomorrow's HPC numbers.

After both arrive, decide between:

- **(a)** Close Round 2 as a negative result. Write Report 6, revise
  the over-splitting framing in prior reports, keep diagnostics for
  future use.
- **(b)** Restart audit on a *specific* fixture if the slice surfaces a
  real Go vs Python divergence.
- **(c)** Re-scope to investigate the `removeEmptyNodes` reporting-layer
  flattening separately, which is what actually drove the K=30 claim.

## Artifacts

- `sim_validation/AUDIT_round2.md` — full Phase 3 audit with line refs
- `sim_validation/round2_trace_K3C2.ndjson` — first traced run
- `trace_test.go` — tests for the `--trace` flag and counters
- `main.go` — `TraceCounters` struct, three increment sites,
  `CapHitsFindNode` diagnostic, `runChain` extended with `traceFile`
  parameter, `--trace` CLI flag
- 17 commits on `go-port` ahead of `origin/go-port` (now pushed by John)
