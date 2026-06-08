# RESOLVED / CORRECTION (2026-06-04)

**Conclusion: there is NO sampler bug and NO regression.** The Go port is a
verified-faithful line-by-line port of PhyloWGS. The "K=10 under-recovery vs
Python" described in the original handoff below was a **validation-harness
artifact**, not a defect in the sampler. The investigation was closed without
any change to the statistical/mathematical behavior of the sampler.

The phantom regression was produced by **three compounding harness/methodology
errors**:

1. **Wrong Python baseline (mismatched seeds).** The Go run analyzed
   (`results/go-cpu/go-cpu/`, seed set **B**, intended dir
   `sim_validation_juno_529`) was compared against an original-Python run from a
   **different seed set A** (`junoresults/results/original-python/`). The
   correctly-matched Python baseline is **`junoresultsmulti`** (seed set **B**,
   216/216 identical seeds). Comparing across seed sets on a stochastic,
   multimodal sampler is not a valid apples-to-apples comparison.

2. **Inconsistent population counting (Go-merged vs Python-raw).** Go's reported
   `npop` (from `best_tree.json` `num_populations`) is counted **after**
   `removeSmallNodes` merging — a population is "small" if its mutation fraction
   is `< max(1%, 2/M)`. The Python `npop` was counted as **raw** non-empty
   populations of the best-LLH tree, with **no** merge step. Counted
   consistently, Go (~5.3) ≈ Python (~5.7) at K=10. The headline gap was a
   counting-convention mismatch, not premature tree truncation.

3. **Multimodal mixing / best-of-N-chains seed luck.** High-K fixtures are
   multimodal; roughly 40% of chains get stuck in a worse mode. A best-of-N
   chains result is therefore luck-sensitive on seeds. Python exhibits the same
   property. Both implementations badly under-recover K=10 — this is inherent
   TSSB difficulty at high K, not a Go-specific failure.

**Note on the cited commit:** the original handoff (below) cites the May 29
sweep as running at `da953bc`. **That commit does not exist in the repo.** The
May 29 sweep was at **`d985693`** (the tip of `go-port` referenced as the
worktree base). Treat `da953bc` as a typo/stale reference.

**Disposition of the original "where to look" leads:** `dirichletSample`,
`resampleSticks`, and the MH proposal were all re-verified against the
Python/C++ references this session and match line-by-line (the MH proposal was
checked against C++ `mh.cpp`, including the asymmetric `+1` in
`sample_cons_params`, the GSL `dirichlet_lnpdf` formula, and the dirichlet
correction term). No divergence was found. See
`docs/analysis/CODE_REVIEW.md` (2026-06-04 note) and
`sim_validation/BUGS_FOUND.md` (harness-bug entries).

> **Everything below this line is the ORIGINAL handoff, retained for history.
> It is SUPERSEDED by the correction above. Do NOT re-chase the "regression" —
> it was never real.**

---

# Handoff — K=10 under-recovery investigation (go-port)  *(SUPERSEDED — see correction above)*

**Date:** 2026-06-03
**Branch:** `claude/go-port-k10-investigation` (off `origin/go-port` @ `d985693`)
**Worktree:** `/Users/orgeraj/hla/Phylowgs/PhyloWGS_refactor/.claude/worktrees/go-port-k10`

## The question to answer

On the May 29 simulation sweep (go-port @ `da953bc`, tip of the `go-port` branch),
the Go port systematically **under-recovers populations on K=10 ground-truth fixtures**
compared to the Python reference run on the same inputs.

This is **not** explained by sample budget — Go ran 8× more trees per fixture
(8 chains × 2,500 = 20,000) and still under-recovered. K=3 and K=5 are fine.
So this is a sampler/mixing issue specific to high-K trees.

## Numbers (paired comparison, 165 fixtures)

`results/go-cpu/go-cpu/` vs `junoresults/results/original-python/`.

**Population-recovery accuracy** (|recovered_npop − ground_truth_K|):

| K  | n  | go median err | py median err | go closer | py closer | tie |
|----|----|---------------|---------------|-----------|-----------|-----|
| 3  | 61 | 1.0           | 1.0           | 27        | 15        | 19  |
| 5  | 58 | 1.0           | 1.0           | 18        | 23        | 17  |
| **10** | **46** | **5.0**       | **4.0**       | **15**    | **23**    | 8   |

K=3 actually better in Go. K=5 a wash. K=10 worse — that's the regression.

Across all 165 paired fixtures, Go ends up with fewer populations than Python
in **102/165** cases (vs more in 26, equal in 37). Combined with the K=10
numbers above, the picture is "Go truncates trees prematurely at high K."

Possibly related (but distinct, not in scope of this investigation):
- 7 silent SLURM kills on K=10/M=500 fixtures (heaviest cells).
- 10 catastrophic-LLH outliers ≤ −25,000, heavily skewed to **C=5**
  (high CNV burden). Both Go and Python have such outliers — neither
  port is robust on C=5. Mentioned for context only.

## Where to look (file + line pointers)

All in `main.go` (single-file Go port).

| Function           | main.go line | Notes |
|--------------------|--------------|-------|
| `dirichletSample`  | 301          | Recent churn here. `4e427b4` added pseudocount; `d4c39ea` removed it; `01c0114` restored it. Verify it matches the C++ `util.cpp` reference. |
| `resampleSticks`   | 1982         | `be7a235` "preserve Python tssb.py:187 depth-0 pin"; `99207e8` "honor TSSB.MinDepth instead of hard-coded 1". This is the most likely site for "tree never grows tall enough." |
| `resampleHypers`   | 2441         | `oversplit_test.go` BLOCKER 1d documents the depth-0 pin for `dp_alpha_llh`. Cross-check the alpha-sampling path. |

The corresponding Python references for parity checking:
- `tssb.py:187` (depth-0 pin in `resample_sticks`)
- `tssb.py:247-255` (`dp_alpha_llh` includes betapdfln at depth 0)
- `util.cpp` dirichlet sample with pseudocount

## What already exists (use, don't rewrite)

- **`sim_validation/BUGS_FOUND.md`** — fully written history of three prior CNV bugs (rebuildAncestorSets, CNVs deep-copy omission, *CNV data race). Read this first — same format your finding should land in.
- **`sim_validation/compare_implementations.py`** — existing Go-vs-Python comparison script. Extend rather than rewrite.
- **`sim_validation/score_results.py`**, `plot_results.py`, `generate_fixtures.py` — fixture generation and scoring.
- **`oversplit_test.go`** — already pins the dp_alpha_llh depth-0 behavior. Use as a template for any new sampler-parity tests you write.
- **`docs/analysis/CODE_REVIEW.md`** — current branch status block (last updated 2026-05-20).
- **`CHANGELOG.md`** — branch-level changelog, append findings here.

## Suggested investigation plan

1. **Read the prior context** (15 min):
   - `sim_validation/BUGS_FOUND.md` (whole file)
   - `docs/analysis/CODE_REVIEW.md`
   - `CHANGELOG.md` "branch" section
   - Commits touching the suspect functions:
     - `01c0114 fix: restore dirichletSample pseudocount`
     - `be7a235 fix: preserve Python tssb.py:187 depth-0 pin`
     - `99207e8 fix: resampleSticks honors TSSB.MinDepth`
     - `4e427b4 fix: add pseudocount to dirichletSample` / `d4c39ea` (reverted)
     - `bf2f74f fix: best tree from best-LLH snapshot + VAF-based orphan reassignment`

2. **Reproduce on one K=10 fixture** (1–2 hr):
   - Pick a fixture where Go under-recovered badly. From the paired data, candidates:
     - `K10_S3_T200_M100_C5_rep1`: go npop=7 (truth=10), py npop=9
     - `K10_S3_T1000_M100_C2_rep0`: go npop=5, py npop=7
     - `K10_S3_T200_M500_C2_rep1`: go npop=3, py npop=6
   - Run go-cpu and the Python reference on it locally, single chain, deterministic seed.
   - Eyeball the per-iteration `NumNodes` column in `chain_*_samples.txt`. Does the Go chain ever explore trees with ≥K populations? If it never visits that region of state space, the issue is upstream of likelihood evaluation (stick-breaking proposes too narrow a prior).

3. **Sampler-parity tests** (half day):
   - Compare `resampleSticks` against `tssb.py` line-by-line in the depth-0 path.
   - Compare `dirichletSample` against `util.cpp` — check that the pseudocount matches and is applied consistently.
   - If you find divergence, write a test in the style of `oversplit_test.go` BLOCKER 1d that pins the correct behavior.

4. **Fix + verify** (half day):
   - Re-run on the chosen K=10 fixture. Population should recover toward 10.
   - Re-run a small batch of K=10 fixtures (say, 10 random reps) to confirm the pattern shifts in aggregate.
   - Update `sim_validation/BUGS_FOUND.md` with a new entry in the existing format.
   - Update `CHANGELOG.md`.

5. **Full re-validation** (out of scope, plan only):
   - Don't re-run the full 216-fixture sweep from this session. Once a fix lands, sketch what re-validation sweep is needed (probably just the 46 K=10 fixtures and the 16 missing/empty ones).

## Data locations

| Path | Contents |
|------|----------|
| `/Users/orgeraj/hla/Phylowgs/phylowgs_cowork/results/go-cpu/go-cpu/` | May 29 go-port run output (216 fixtures, 209 non-empty). Per fixture: `best_tree.json`, `chain_0..7_samples.txt`, `summary.json`, `timing.json`. |
| `/Users/orgeraj/hla/Phylowgs/phylowgs_cowork/junoresults/results/original-python/` | Apr 21 Python reference output (1,152 fixtures, 830 non-empty). Per fixture: `mutass.zip`, `tree_summaries.json.gz`, `timing.json`. |
| `/Users/orgeraj/hla/Phylowgs/phylowgs_cowork/junoresults/fixtures/` | Source fixtures (SSM/CNV input files). |
| `/Users/orgeraj/hla/Phylowgs/phylowgs_cowork/sim_validation_juno_529/` | **EMPTY** — output path the May 29 orchestrator was supposed to write to. Real output landed in `results/go-cpu/go-cpu/` instead. Path-config bug, not in scope but worth noting in the eventual writeup. |

Python reference's `tree_summaries.json.gz` structure:
```python
import gzip, json
with gzip.open("tree_summaries.json.gz") as f: d = json.load(f)
# d["trees"]  -> dict keyed by tree_idx str
# d["trees"]["0"]["llh"], ["populations"] (populations is a dict, len = num_pops)
```

## Caveats / things to NOT spend time on

- **Don't chase the C=5 LLH outliers in this task.** Both implementations have them; that's a separate investigation about CNV likelihood stability, not the K=10 issue.
- **Don't worry about the 7 missing K=10/M=500 fixtures.** Almost certainly SLURM walltime/OOM kills (silent — no exit code, no files). Re-submit on juno with bigger budget later, separately.
- **Don't trust the headline "Go LLH wins 112/165"** if it shows up in any docs. That comparison was Go's 20k trees vs Python's 2.5k — 8× sample budget difference. Median Δ(llh) was only +241; the wins are inflated by sample count. Population recovery is the apples-to-apples comparison and that's where Go loses.

## Starting commands

```bash
cd /Users/orgeraj/hla/Phylowgs/PhyloWGS_refactor/.claude/worktrees/go-port-k10

# build & run baseline
go build -o phylowgs-go ./...
# (or `go test ./...` to make sure tests pass before changing anything)

# read prior bug history
$EDITOR sim_validation/BUGS_FOUND.md docs/analysis/CODE_REVIEW.md CHANGELOG.md

# inspect suspect functions
$EDITOR main.go +301   # dirichletSample
$EDITOR main.go +1982  # resampleSticks
$EDITOR main.go +2441  # resampleHypers
```
