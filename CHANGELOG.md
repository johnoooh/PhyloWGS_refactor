# Changelog

All notable changes to the Go port (`phylowgs-go`) are recorded here.
Format: reverse-chronological, grouped by date. Each entry references the
landing commit short SHA on `go-port`.

## 2026-05-07

### Fixed — best tree selection and orphan reassignment ([`459b2d2`](https://github.com/johnoooh/PhyloWGS_refactor/commit/459b2d2))

Two post-processing bugs that diverged from the Python reference and
contributed to a measurable AUPRC gap.

- **`writeResults` now writes the actual best-LLH snapshot.** Previously,
  it located the chain holding the best sample LLH but then dumped that
  chain's `FinalTree` (the last MCMC state). Now uses
  `Trees[bestSampleIdx].Snapshot` directly. Added `pickBestSample`
  helper with explicit error paths for empty / mismatched chains, plus
  a runtime assertion that the selected sample's LLH matches.
- **`removeSmallNodes` now does VAF-based per-mutation reassignment.**
  Previously dumped every orphan into the highest-CP population. Now
  matches Python's `_move_muts_to_best_node` (`result_munger.py:247-269`):
  per mutation, compute `implied_phi = 2 * (D - A) / D`, find the
  non-root population with closest mean cellular prevalence, send it
  there. SSMs and CNVs go through identical math; both have
  `ref_reads`/`total_reads` in Python's mutlist.
  - Sorted-key iteration → ties resolve to lowest pidx (matches Python's
    insertion-order dict iteration).
  - Mutations with no read data are silently dropped.
- **Adaptive small-node threshold preserved**: `max(1%, 2/M)`.
  Validated improvement over flat 0.01 in earlier benchmarks.

Tests: 9 new in `main_postprocess_test.go` covering low/high VAF
reassignment, CNV path, deterministic tie-breaking (20-trial run),
missing-data drop, zero-total-reads drop, and the extracted
`pickBestSample` helper (mid-chain best, multi-chain, included-only,
ties, empty, length-mismatch error paths).
