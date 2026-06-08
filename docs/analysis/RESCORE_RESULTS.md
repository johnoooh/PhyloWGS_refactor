# Re-scored Go-vs-Python comparison (matched seed-B, fixed harness)

_Date: 2026-06-04. Produced after fixing the validation harness (consistent
population counting + matched-seed guard). Supersedes the HANDOFF.md "K=10
under-recovery" numbers, which were an artifact of (1) a mismatched Python
baseline and (2) inconsistent population counting._

## Method

- **Matched data**: all results below are on the SAME seed set B (216/216
  fixtures share `truth.json` `params.seed`), so every Go-vs-Python pair is the
  same simulated instance.
  - Python: `junoresultsmulti/results/original-python` (seed-B).
  - Go (April build): `junoresultsmulti/results/go-cpu` (seed-B, same sweep).
  - Go (May/current build, d985693): `sim_validation_juno_529/.../results/go-cpu` (seed-B).
- **Consistent counting**: population counts use the *merged* rule for BOTH
  implementations — a population is dropped if its size `< max(0.01*M, 2)`,
  mirroring Go's `removeSmallNodes` threshold `max(1%, 2/M)`. (Go's
  `best_tree.json` is already merged; the harness now applies the same merge to
  the Python tree, instead of counting Python's raw nodes.)
- Scored with the fixed `sim_validation/score_results.py`.

## Table 1 — Population recovery (merged counts)

`MedErr` = median |inferred_npop − true_K|. `[raw mean go/py]` shows the
pre-merge means to expose the counting artifact.

**GO April vs matched Python**

| K  | n  | Go MedErr | Py MedErr | Go mean | Py mean | GO closer | PY closer | tie | raw mean go/py |
|----|----|-----------|-----------|---------|---------|-----------|-----------|-----|----------------|
| 3  | 56 | 1.0       | 1.0       | 3.1     | 3.2     | 12        | 8         | 36  | 3.1 / 3.5      |
| 5  | 59 | 2.0       | 2.0       | 3.9     | 3.9     | 9         | 19        | 31  | 3.9 / 4.2      |
| 10 | 56 | 4.0       | 5.0       | 6.4     | 5.7     | 19        | 13        | 24  | 6.4 / 8.3      |

**GO May/current vs matched Python**

| K  | n  | Go MedErr | Py MedErr | Go mean | Py mean | GO closer | PY closer | tie | raw mean go/py |
|----|----|-----------|-----------|---------|---------|-----------|-----------|-----|----------------|
| 3  | 58 | 1.0       | 1.0       | 3.2     | 3.4     | 13        | 4         | 41  | 3.2 / 3.6      |
| 5  | 61 | 2.0       | 2.0       | 3.9     | 4.0     | 6         | 12        | 43  | 3.9 / 4.3      |
| 10 | 56 | 5.0       | 5.0       | 5.3     | 5.7     | 9         | 15        | 32  | 5.3 / 8.3      |

**Key point:** the headline "Go under-recovers at K=10" disappears under
consistent counting. Python's K=10 mean drops from **8.3 raw → 5.7 merged** — that
~2.6-population collapse *is* the artifact (Python's raw tree is full of tiny
nodes Go merges away). Counted the same way, Go ≈ Python at every K, dominated by
ties. The April↔May swing at K=10 (Go mean 6.4 → 5.3) is multimodal seed luck
(proven separately: ~40% of chains stick in a worse mode; not a code change).

## Table 2 — Co-clustering AUPRC (mutation-pair clustering accuracy; higher = better)

**GO April vs matched Python**

| K  | n  | Go mean | Py mean | Go median | Py median | Go wins | Py wins |
|----|----|---------|---------|-----------|-----------|---------|---------|
| 3  | 56 | 0.688   | 0.756   | 0.705     | 0.811     | 4       | 32      |
| 5  | 59 | 0.531   | 0.608   | 0.512     | 0.586     | 12      | 39      |
| 10 | 56 | 0.420   | 0.456   | 0.325     | 0.395     | 11      | 40      |

**GO May/current vs matched Python**

| K  | n  | Go mean | Py mean | Go median | Py median | Go wins | Py wins |
|----|----|---------|---------|-----------|-----------|---------|---------|
| 3  | 58 | 0.764   | 0.761   | 0.817     | 0.823     | 12      | 16      |
| 5  | 61 | 0.606   | 0.616   | 0.600     | 0.622     | 14      | 22      |
| 10 | 56 | 0.455   | 0.456   | 0.412     | 0.395     | 27      | 20      |

**Key point:** the April Go build was clearly *behind* Python on AUPRC (Py won
32–40 of ~56 at each K). The **May/current build is at PARITY** with Python
across all K (means within 0.001–0.010), and at K=10 it wins more head-to-head
(27 vs 20; higher median). This is an **improvement, not a regression** — driven
by the May 7 post-processing fixes (`bf2f74f`: best-LLH snapshot selection +
VAF-based per-mutation reassignment matching Python's `_move_muts_to_best_node`).

## Conclusions

1. **No regression, no sampler bug.** Population recovery: Go ≈ Python at all K
   once counted consistently. AUPRC: current Go is at parity with Python (a gain
   over the April build).
2. The original "K=10 Go under-recovery" was a harness artifact (mismatched
   baseline + inconsistent counting), now fixed.
3. Residual K=10 npop variance is inherent multimodal mixing, shared with Python.

## Caveats / next

- The "May/current" Go numbers come from the `juno_529` run, which was a single
  8-chain run subject to multimodal seed luck on the hard cells. A clean
  cluster re-run (`submit_matched_validation_lsf.sh`) with **equal, generous
  chain counts for both implementations** will give a fair, low-variance,
  current-code number and is expected to confirm/tighten these results.
- AUPRC here is co-clustering AUPRC over mutation pairs (the metric
  `score_results.py` computes). It does not depend on the merge rule.
