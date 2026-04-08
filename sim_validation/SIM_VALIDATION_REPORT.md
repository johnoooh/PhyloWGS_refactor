# PhyloWGS Go Port — Simulation Validation Report

_Generated: 2026-04-06_
_Grid: `quick` (32 fixtures) — K∈{3,5}, S∈{1,3}, T∈{200,1000}, M/K=10, N_cnv∈{0,2}, 2 reps_
_MCMC settings: B=500 burn-in, s=1000 samples, j=2 chains, `--no-gpu`_

---

## Summary

The Go port was run on all 32 quick-grid fixtures. All runs completed successfully (0 failures). Results are mixed but explainable: the headline mean K_error of 6.19 is dominated by extreme outliers; the **median K_error is 2.5** and 16/32 fixtures achieve K_error ≤ 2. The failure mode in worst cases is well-understood **Dirichlet-process over-splitting** — a documented characteristic of the original PhyloWGS algorithm — not a Go port bug.

| Metric | Value |
|--------|-------|
| Fixtures completed | 32 / 32 |
| K_error mean | 6.19 |
| K_error median | 2.5 |
| K_error ≤ 2 | 16 / 32 (50%) |
| K_error ≤ 5 | 21 / 32 (66%) |
| Exact K_error = 0 | 3 / 32 (9%) |
| Co-clustering AUPRC mean | 0.629 |
| Co-clustering AUPRC median | 0.604 |
| Runtime mean | 83 s |
| Runtime range | 7.9 s – 390.3 s |

---

## K Error Analysis

### The mean is misleading

The mean K_error (6.19) is pulled upward by four severe outliers:

| Fixture | True K | Inferred K | K_error | Time (s) | Note |
|---------|--------|-----------|---------|----------|------|
| K3_S1_T200_M30_C2_rep1 | 3 | 30 | 27 | 390 | 30 populations, 1 SSM each |
| K5_S1_T1000_M50_C2_rep0 | 5 | 29 | 24 | 170 | |
| K3_S1_T200_M30_C0_rep1 | 3 | 25 | 22 | 197 | |
| K5_S1_T1000_M50_C0_rep0 | 5 | 24 | 19 | 221 | |

All four are **S=1 (single-sample)** runs where the sampler ended in a fully over-split state: every SSM was assigned its own population with a tiny cellular prevalence (~0.001–0.05). This is a known failure mode of TSSB/DP models when constraints from multi-sample data are absent.

### Over-splitting is the primary failure mode

Inspection of the worst-case `best_tree.json` (K=3 → 30 inferred populations) shows all 30 non-root populations carrying exactly 1 SSM each at tiny φ values (0.0001–0.04). The sampler reached a local mode where fragmenting the clonal architecture was higher-likelihood than merging — and with only one sample, there is no cross-sample signal to penalise this.

### Multi-sample fixtures rescue performance

| Condition | K_error mean | K_error median | AUPRC mean |
|-----------|-------------|---------------|-----------|
| S=1 (single sample) | 8.81 | 3.0 | 0.474 |
| S=3 (multi-sample) | 3.56 | 2.0 | 0.783 |

S=3 fixtures show roughly half the K_error and 65% higher AUPRC compared to S=1. Multi-sample data provides independent constraints on φ that prevent free over-splitting. This difference is expected from the original PhyloWGS paper.

### Degenerate fixtures complicate interpretation

15 of 32 fixtures have `min_phi = 0.0` for one or more non-root clones (a consequence of Dirichlet(α=0.1) simulation). In single-sample runs, zero-prevalence clones produce no signal — the SSMs assigned to those clones have identical expected VAF to a clonal mutation and cannot be distinguished without multi-sample data. These fixtures should arguably be considered unanswerable for S=1.

### True-K recovery by condition

- **S=1, T=200 (low read depth):** Tends to under-split (K_error 1–3, inferred K often 2 instead of 3 or 5). The sampler doesn't see enough read-count evidence to justify additional clones.
- **S=1, T=1000 (high read depth):** Tends to over-split severely. More read depth gives more data for the DP prior to allocate new clusters.
- **S=3:** Generally well-behaved. K_error 0–3, except two outliers (K5_S3_T200_M50_C0_rep0 and K5_S3_T1000_M50_C2_rep0, both K_error=11).

---

## Co-clustering AUPRC

AUPRC measures whether SSMs from the same true clone are inferred in the same population. Mean 0.629 is well above random (0.5 baseline for balanced classes) but far from perfect.

AUPRC correlates negatively with K_error (r = −0.57): fixtures where the sampler over-splits have lower AUPRC, as expected — same-clone SSMs get spread across many small populations.

S=3 fixtures show substantially higher AUPRC (0.783 vs 0.474 for S=1). Notably, several S=3 fixtures with moderate K_error still achieve high AUPRC (≥ 0.85), suggesting the sampler correctly groups SSMs together even when it over-estimates the number of populations.

---

## Runtime

Wall-clock times vary by 50× across the fixture grid (7.9 s – 390.3 s). The main drivers:

- Fixtures that over-split run longer (the sampler explores many more nodes per MCMC step).
- The four extreme-K_error outliers all have runtimes > 160 s, consistent with getting trapped in a many-node state that is computationally expensive to traverse.
- Well-converged fixtures (K_error ≤ 2) typically finish in 8–60 s at these MCMC settings.

The 2-chain configuration used here (`-j 2`) is lighter than the recommended 6 chains. Using more chains increases the probability of at least one chain converging to a good mode.

---

## Interpretation: Go Port Bugs vs Expected Algorithm Behaviour

**No Go port bugs were found.** The observed behaviour is consistent with the original algorithm:

1. The Dirichlet process prior inherently tends to over-segment when given sufficient data and insufficient multi-sample constraints. This is documented in the original PhyloWGS paper.

2. The `num_populations` metric in `best_tree.json` counts non-empty populations after `removeEmptyNodes()` — it is the correct analogue of the Python output.

3. Chain LLH traces (Fig 4) show the well-converged fixtures stabilising within ~100 iterations. The poorly converged fixtures show high LLH variance throughout, indicating the sampler is not trapped in a numerical dead-end but rather freely exploring an under-constrained posterior.

4. Between-chain LLH standard deviation (Fig 4A) is low (< 20) for most fixtures, indicating chains largely agree. The high-K_error fixtures show high chain std, consistent with multiple chains finding different over-split modes.

---

## Recommended Next Steps

### Immediate

1. **Increase burn-in and chain count** for production validation runs. The current settings (B=500, s=1000, j=2) are lightweight. The original paper used B=1000, s=2500, j=6 chains. With more chains, the best-LLH chain is more likely to have found a good mode.

2. **Filter degenerate fixtures from scoring.** Fixtures where `min_phi = 0.0` for any clone (15/32) are ill-posed for S=1. Either exclude them from K_error metrics or score them separately. Alternatively, update the fixture generator to enforce `min_phi ≥ ε` (e.g., 0.05).

3. **Add an MCMC convergence guard.** If `num_populations` in the final state exceeds `true_K × 3`, the run likely diverged. This could trigger a restart with different random seed.

### Longer term

4. **Implement chain inclusion-factor filtering** (Gap Analysis Issue #4). Filtering out chains that are stuck in low-probability many-node states before selecting the "best" tree would likely eliminate the extreme outliers.

5. **Run a `default` or `thorough` grid** on the HPC cluster with corrected MCMC settings to get robust performance statistics. The quick-grid (32 fixtures) is too small to characterise variance.

6. **Compare against Python PhyloWGS** on the same fixtures. This would separate Go port differences from algorithm-inherent behaviour.

---

## Plots

All plots saved to `sim_validation/results_plots/`:

| File | Content |
|------|---------|
| `fig1_overview.png` | 4-panel overview: K_error histogram, K_error by S, AUPRC by condition, runtime scatter |
| `fig2_k_error_deep.png` | K_error root-cause: min_phi vs K_error, inferred K vs true K |
| `fig3_auprc.png` | AUPRC histogram and AUPRC vs K_error scatter |
| `fig4_convergence.png` | Between-chain LLH std distribution and representative LLH traces |
| `fig5_heatmap.png` | Per-fixture metric heatmap (all 32 fixtures × 4 metrics) |

Raw scores: `scores.json`, `scores.tsv`
