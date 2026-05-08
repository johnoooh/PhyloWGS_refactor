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

### Changed — `dirichletSample` no longer adds a pseudocount

- `dirichletSample` no longer adds the C++-style 0.0001 pseudocount.
  The pseudocount masks rather than fixes -Inf likelihoods on near-zero
  pi components and biases acceptance ratios. Validation showed our
  downstream binomial likelihood already clamps mu to [1e-15, 1-1e-15],
  so a zero pi cannot produce -Inf in the LLH path. If a numerical issue
  arises downstream, fix it at the LLH site, not in the proposal.
- Added `TestDirichletSample_NormalizesToOne` in `main_postprocess_test.go`
  to lock the contract: 1000 trials with alpha=[1,1,1,1] must produce
  non-negative components summing to 1 within 1e-12.

### Added

- **`trees.zip` and `mutass.zip` outputs after MCMC**
  ([`c4e0160`](https://github.com/johnoooh/PhyloWGS_refactor/commit/c4e0160),
  [`f1cf897`](https://github.com/johnoooh/PhyloWGS_refactor/commit/f1cf897),
  [`8aeeee7`](https://github.com/johnoooh/PhyloWGS_refactor/commit/8aeeee7)).
  Mirrors Python's `multievolve.py` / `write_results.py` file layout
  (JSON content, not pickle). Aggregated across all chains that pass
  the `-I` inclusion filter and indexed by global tree number.
  - `trees.zip` entries are named `tree_<idx>_<llh>` and contain the
    per-sample tree snapshot.
  - `mutass.zip` entries are named `<idx>.json` and contain the
    mutation-assignment payload for that sample.
- **`-D / --dataset-name` flag**
  ([`be3265c`](https://github.com/johnoooh/PhyloWGS_refactor/commit/be3265c)).
  Embedded in mutass.zip entries so downstream tooling can identify
  which run produced a given archive.
- **`phylowgs-go posterior-trees [-n N] [-no-pdf] <output-dir>`
  subcommand**
  ([`3ce059f`](https://github.com/johnoooh/PhyloWGS_refactor/commit/3ce059f),
  [`28a7a15`](https://github.com/johnoooh/PhyloWGS_refactor/commit/28a7a15),
  [`9962dbd`](https://github.com/johnoooh/PhyloWGS_refactor/commit/9962dbd),
  [`9f025dc`](https://github.com/johnoooh/PhyloWGS_refactor/commit/9f025dc),
  [`f73b03b`](https://github.com/johnoooh/PhyloWGS_refactor/commit/f73b03b),
  [`109c813`](https://github.com/johnoooh/PhyloWGS_refactor/commit/109c813),
  [`94307be`](https://github.com/johnoooh/PhyloWGS_refactor/commit/94307be)).
  Full port of `posterior_trees.py`: reads `trees.zip`, groups samples
  by topology signature, ranks groups by posterior probability, and
  writes per-group standalone-LaTeX summaries (with PDFs, when
  `pdflatex` is on PATH).

### Fixed

- **SSM `mu_r` / `mu_v` default to 0 when columns are absent**
  ([`5ee4450`](https://github.com/johnoooh/PhyloWGS_refactor/commit/5ee4450)).
  Matches Python `util2.py`. Previously the Go parser hard-coded
  defaults of 0.999 / 0.5, which silently changed likelihoods on
  header-less or short-form inputs.
- **`resampleSticks` honors `TSSB.MinDepth`**
  ([`21056b0`](https://github.com/johnoooh/PhyloWGS_refactor/commit/21056b0))
  instead of a hard-coded `depth >= 1`. Matches Python `tssb.py:186`.
  Latent bug; default configuration unaffected.
