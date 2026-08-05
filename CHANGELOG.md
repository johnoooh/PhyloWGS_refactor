# Changelog

All notable changes to the Go port (`phylowgs-go`) are recorded here.
Format: reverse-chronological, grouped by date. Each entry references the
landing commit short SHA on `go-port`.

## 2026-08-05 — branch `go-port`

Closes out the confirmed findings from the 2026-06-24 equivalence audit
(`docs/analysis/2026-06-24-go-python-equivalence-audit.md`). Each item
below states whether it changes inference output or only output
labeling/reporting — only M1 changes inference.

### Fixed — M5, b1, M2 (schema/output only)
[`be78801`](https://github.com/johnoooh/PhyloWGS_refactor/commit/be78801)

- `mutlist.json`'s `physical_cnvs[].cell_prev` is now parsed to
  `[]float64` (was a raw `"0.0|0.718"` string), matching `util2.py:70`.
- `mutlist.json` now carries a top-level `dataset_name`, matching
  `json_writer.py`.
- The GPU branch of `completeDataLogLikelihood` was silently dropping
  all CNV-datum log-likelihood contributions (unreachable on default
  stub-CUDA builds, but real on `-tags cuda`); extracted
  `cnvDatumLogLikelihood` and now call it from both branches.

### Fixed — rounding-mode bug (output only)
[`29f3374`](https://github.com/johnoooh/PhyloWGS_refactor/commit/29f3374)

`removeSmallNodes`'s pruning threshold used `math.Round`
(half-away-from-zero) where Python 3's `int(round(...))` is
round-half-to-even, causing permanent divergence at specific mutation
counts (e.g. M=250: Python gives 2, Go gave 3) independent of the
(intentionally-kept) adaptive-threshold question. Added
`roundHalfEven`. Note: Go's adaptive `minFrac=max(1%, 2/M)` threshold
itself was investigated and left as-is — Python's fixed-1% rule
degenerates to a no-op below M≈100 under banker's rounding, so it isn't
a considered design choice to replicate, and Go's version measurably
beats it on K-error accuracy.

### Fixed — M7: posterior frequency aggregation keyed by mutation content (output only)
[`4c4403b`](https://github.com/johnoooh/PhyloWGS_refactor/commit/4c4403b)

`aggregateFreqsByNode` pooled cellular-prevalence stats by node ID, but
node IDs are per-sample DFS ordinals with no cross-sample identity —
whenever two MCMC samples in the same posterior-topology group had
siblings whose phi rank swapped between samples, the pooled mean/std in
`posterior_trees/*.tex` was silently an order statistic rather than
each mutation set's own posterior. Now keyed by the sorted-joined
SSM+CNV name set (`nodeMutKey`), matching `posterior_trees.py:88-107`.

### Fixed — M4 + M6: population renumbering and topology-signature timing (output only)
[`54e4f3a`](https://github.com/johnoooh/PhyloWGS_refactor/commit/54e4f3a)

- **M4**: `postProcessSummary` renumbered surviving populations by
  ascending old index after removing empty nodes, rather than Python's
  post-removal phi-descending order (`result_generator.py:81-86`).
  Fixed by walking the reparented tree in phi-descending preorder.
  Deleted `removeEmptyNodes` as confirmed-dead code — it would have
  corrupted the live MCMC tree if ever wired into the hot path, since
  it has no deep-copy support.
- **M6**: the topology signature was always computed on an
  already-empty-node-cleaned snapshot, so its (already-faithful)
  genotype-inheritance-break logic could never fire — topologies Python
  keeps separate were silently over-merged. Fixed by computing and
  stamping a `topology_signature` field in `snapshotTree`, before
  `postProcessSummary` runs.

### Fixed — M1: MH kernel now matches the reference's mh.hpp, not data.py, for co-located CNV-SSMs (**inference output changes**)
[`70b7187`](https://github.com/johnoooh/PhyloWGS_refactor/commit/70b7187)

The reference pipeline runs `data.py` for assignment resampling /
complete-data LLH and a separate C++ `mh.hpp`/`mh.o` kernel for MH phi
sampling. Go had collapsed all three call sites onto `data.py`'s
semantics — correct for two of three, wrong for MH sampling. For an SSM
co-located on its own zero-copy-allele CNV (LOH `minor_cn=0`, or
homozygous deletion), this biased phi estimates via a systematically
different prior/channel-mix than the reference actually uses. Ported
the reference's static per-iteration state collapse
(`params.py:213-219`) into `computeSSMStates`, and replaced the MH
kernel's dynamic `nv>0`-filtered reduction with the reference's
unconditional four-term, `nr+nv>0`-guarded, fixed-`log(0.25)`-prior
reduction. Assignment resampling and complete-data LLH are unchanged.
Verified against an independent hand-derived LOH fixture (not derived
from any Go code path), an algebraic invariance check for the common
non-co-located case (1e-12 tolerance), and a full multi-chain
LOH-bearing MCMC run clean under `-race`. **Not yet verified**: a full
`sim_validation/` run against the actual Python/C++ pipeline on HPC —
no local Python/GSL/`mh.o` build environment was available this
session; recommended before treating this as validated for a
publication-facing run.

## 2026-05-08 — branch `claude/go-port-gaps-and-correctness`

Closes the remaining behavioral gaps with the Python `multievolve.py` /
`posterior_trees.py` pipeline and fixes three latent or silent-divergence
bugs in the MCMC layer.

### Added — `trees.zip` and `mutass.zip` outputs after MCMC

Mirrors Python's `multievolve.py` / `write_results.py` archive layout
(JSON payloads instead of pickle). Aggregated across all chains that
pass the `-I` inclusion filter and indexed by global tree number.

- `trees.zip` entries are named `tree_<idx>_<llh>` and contain the
  per-sample tree snapshot
  ([`d7aa848`](https://github.com/johnoooh/PhyloWGS_refactor/commit/d7aa848),
  [`b508915`](https://github.com/johnoooh/PhyloWGS_refactor/commit/b508915)).
- `mutass.zip` entries are named `<idx>.json` and contain the
  mutation-assignment payload
  ([`77f83c7`](https://github.com/johnoooh/PhyloWGS_refactor/commit/77f83c7)).
- New `-D / --dataset-name` flag is embedded in `mutass.zip` entries so
  downstream tooling can identify which run produced a given archive
  ([`d599e49`](https://github.com/johnoooh/PhyloWGS_refactor/commit/d599e49)).

### Added — `posterior-trees` subcommand

```
phylowgs-go posterior-trees [-n N] [-no-pdf] <output-dir>
```

Full port of `posterior_trees.py`: reads `trees.zip`, groups samples by
topology signature, ranks groups by posterior probability, and writes
per-group standalone-LaTeX summaries to
`<output-dir>/posterior_trees/tree_<rank>_<probability>.tex`. PDFs are
produced alongside when `pdflatex` is on `PATH`.

Implementation landed in pieces:
[`67f387d`](https://github.com/johnoooh/PhyloWGS_refactor/commit/67f387d) `TreeArchiveReader`,
[`c5370ca`](https://github.com/johnoooh/PhyloWGS_refactor/commit/c5370ca) topology signatures,
[`5c16082`](https://github.com/johnoooh/PhyloWGS_refactor/commit/5c16082) grouping + ranking,
[`6604333`](https://github.com/johnoooh/PhyloWGS_refactor/commit/6604333) LaTeX writer,
[`549bc52`](https://github.com/johnoooh/PhyloWGS_refactor/commit/549bc52) per-group phi aggregation,
[`7a43968`](https://github.com/johnoooh/PhyloWGS_refactor/commit/7a43968) `runPosteriorTrees` driver,
[`6e061ae`](https://github.com/johnoooh/PhyloWGS_refactor/commit/6e061ae) subcommand wiring.

### Reverted — `dirichletSample` keeps the `0.0001` pseudocount

A prior commit on this branch (`d4c39ea`) removed the pseudocount on the
argument that the binomial-likelihood `mu` clamp made it unnecessary.
Code review caught that the clamp protects `mu` in the SSM binomial
path, but `piNew` from `dirichletSample` is also passed directly to
`dirichletLogPDF` in the MH acceptance step, which does
`(alpha[i]-1)·log(x[i])` and returns `-Inf` whenever `x[i] <= 0`. A
zero-valued draw therefore stalls the sampler at the affected node and
biases toward fewer-node trees — exactly the failure mode the original
fix in `4e427b4` ("add pseudocount to dirichletSample to match C++
util.cpp") was added to prevent. The pseudocount is restored.

Locked by `TestDirichletSample_PseudocountFloor`: 1000 trials with
`alpha=[0.0001, 1000, 1000]` must produce strictly positive components
(min `~5e-8` given the pseudocount). `TestDirichletSample_NormalizesToOne`
still verifies the simplex constraint.

### Fixed — SSM `mu_r` / `mu_v` default to 0 when columns are absent ([`37ccfb6`](https://github.com/johnoooh/PhyloWGS_refactor/commit/37ccfb6))

Matches Python `util2.py:63` (`mu_r=mu_v=0` when columns absent) and
the `Datum` class defaults in `data.py:7`. Previously the Go parser
hard-coded defaults of `0.999` / `0.5`, which silently changed
likelihoods on header-less or short-form inputs without warning.

### Fixed — `resampleSticks` honors `TSSB.MinDepth` *and* preserves the depth-0 pin ([`99207e8`](https://github.com/johnoooh/PhyloWGS_refactor/commit/99207e8))

The MinDepth gate was previously hard-coded to `depth >= 1`. The
initial fix changed it to `t.MinDepth <= depth` to match Python
`tssb.py:186`, but missed `tssb.py:187`, which unconditionally pins
`root['main'] = 1e-30` at depth 0 ("shankar"). With the default
`MinDepth=0` the old `depth >= 1` accidentally produced the same end
state, so the regression was only visible against the Python reference.
Both lines of the Python contract are now ported: the conditional
`boundBeta` resample plus the unconditional depth-0 pin.

Locked by `TestResampleSticks_RootStickIsPinnedAtDepthZero`: after
`resampleSticks` from any seed and starting value, `root.Main` must
equal `1e-30` within `1e-40`.

### Docs

- README updated with the `-I`, `-D` flags, the `mu_r` / `mu_v` default
  change, the new `Outputs` section listing every file produced by a
  run, and a `Posterior-tree summaries` section documenting the new
  subcommand
  ([`2143b0a`](https://github.com/johnoooh/PhyloWGS_refactor/commit/2143b0a)).
- This CHANGELOG entry.

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
