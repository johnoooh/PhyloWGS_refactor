# Go-native PhyloWGS results: drop Python from the PhyloWGS step

Status: for approval (2026-11-02)

## Intent

Make the Go port a drop-in replacement for PhyloWGS inside `mskcc/neoantigen-pipeline`.
Python leaves the PhyloWGS sampling and result-writing steps entirely: no
`multievolve`, no python2 `write_results.py`, no second container. Nothing downstream
of PhyloWGS (`NEOANTIGENUTILS_NEOANTIGENINPUT` and the final neoantigen tree JSON it
builds) may be able to tell which engine ran.

Why: Go is 2–16× faster, it now survives straggler chains (per-chain flush,
`-chain-timeout`), and a single engine/container is simpler to maintain than the
current Go → python2 adapter (`bin/write_results_go.py`) chain.

## What exists today

- Pipeline branches `feat/phylowgs-go-comparison` / `feat/phylowgs-go-on-rc-2.0.0`:
  `params.phylowgs_engine = python | go | both`. The Go path is
  `PHYLOWGS_MULTIEVOLVEGO` (writes JSON `trees.zip` + `mutlist.json`) →
  `PHYLOWGS_WRITERESULTSGO` (python2 container; `write_results_go.py` loads Go's JSON
  trees and calls Python's own `ResultMunger` + `JsonWriter`) → `summ/muts/mutass`.
- Python's result writing (`pwgsresults/` in `mskcc/phylowgs`, ~500 lines):
  1. `ResultGenerator` — per-tree `{llh, structure, populations}` + `mut_assignments`
     and the mutation list. Go's `summarizePops` / `writeMutList` already match this
     (audited; that is what makes the adapter work).
  2. `ResultMunger` — `remove_small_nodes(min_ssms)` (moves SSMs/CNVs of small nodes to
     the node with the closest implied phi, renumbers), `remove_superclones()`,
     `remove_multiprimary_trees(max_multiprimary)` (deletes and renumbers trees, raises
     above the threshold).
  3. `IndexCalculator` — `linearity_index`, `branching_index`, `clustering_index`
     per tree.
  4. `JsonWriter` — `<name>.summ.json.gz` `{dataset_name, params, trees,
     tree_densities}` where `tree_densities` = `scipy.stats.gaussian_kde` over
     (clustering, branching/(branching+linearity+ε)) with 1-D and all-zero fallbacks;
     `<name>.muts.json.gz` = mutlist + `dataset_name`; `<name>.mutass.zip` with one
     `<tidx>.json` = `{mut_assignments, dataset_name}` per tree.
- Downstream (`neoantigen-utils` `generate_input.py`) reads only: `summ.trees[*].llh`,
  `.structure`, `.populations[*].cellular_prevalence[0]`; `mutass/<tidx>.json`
  `.mut_assignments[clone].ssms`; `muts.ssms[id].name`. The indices and densities are
  not read by the pipeline but are part of the published files, so they are in scope.

## Goals

1. `phylowgs-go` writes `<prefix>.summ.json.gz`, `<prefix>.muts.json.gz`,
   `<prefix>.mutass.zip` itself, with the same CLI options as `write_results.py`
   (`--min-ssms`, `--include-multiprimary`, `--keep-superclones`, `--max-multiprimary`,
   `--include-ssm-names`) and the same defaults.
2. **Format equivalence**, defined as: on the same Go trees, the Go-native files and the
   files from `write_results_go.py` (the Python reference) are equal after decoding —
   same keys, same nesting, same JSON types (int vs float: `1` vs `1.0` matters), same
   tree indices and node numbering, same SSM/CNV lists (order-insensitive where Python's
   order comes from dict iteration). Floats equal exactly, except `tree_densities`
   (relative tolerance 1e-9; LAPACK vs Go linear algebra) — see "Densities" below.
3. Same failure behaviour: >= `max_multiprimary` multiprimary trees → non-zero exit with
   the same message text.
4. Pipeline: `--phylowgs_engine go` runs two Go-only processes from one image —
   `PHYLOWGS_MULTIEVOLVEGO` (sampling, writes `trees.zip` + `mutlist.json`) and
   `PHYLOWGS_WRITERESULTSGO` (now `phylowgs-go write-results`, no python2) — emitting
   `summ/muts/mutass` named `<meta.id>.*` (no `-go` suffix, dataset name = `meta.id`).
   It is a drop-in for `PHYLOWGS_MULTIEVOLVE` + `PHYLOWGS_WRITERESULTS`. Keeping writing
   as its own step is deliberate (claude's call): a writer failure (multiprimary raise,
   munger assert) or a filter change then re-runs in seconds with `-resume` instead of
   discarding a multi-day sampling task's outputs.
5. End-to-end check: on a few real samples, `--phylowgs_engine go` produces a final
   neoantigen `*_input.json` with the same structure (keys, types, nesting) as a Python
   run. Content differs (different MCMC draws) and is not compared.

## Non-goals

- Picking the same trees as Python. MCMC differs; posterior similarity is covered by
  the existing `sim_validation/` / NeoQual comparisons.
- Byte-identical files: Python 2 dict key order is hash order, gzip headers carry
  mtimes. Equality is on decoded content.
- Changing `PHYLOWGS_PARSECNVS` / `PHYLOWGS_CREATEINPUT`: they stay Python (decided
  2026-11-02); "no Python" applies to sampling and result writing only.
- Changing anything in `neoantigen-utils`.

## Anti-goals (fails even if it "works")

- Any downstream field changes type or name, or tree/node numbering differs from what
  Python's munger would produce on the same trees.
- Silently writing an empty or partial `summ`/`mutass` (e.g. all trees filtered) where
  Python would raise.
- Removing the `python` engine before the Go path has passed the end-to-end check.

## Decided

- Python reference = `write_results_go.py` (pipeline repo `bin/`) running in the
  existing `phylowgs:v1.5-msk` container on Go's `trees.zip` + `mutlist.json`. It uses
  Python's own munger/writer, so matching it on the same trees means matching Python.
- Filter settings follow the latest release candidate, 2.0.0 (`feature/release_candidate_2.0.0`
  `conf/modules.config`): `--max-multiprimary 1.0 --include-multiprimary`, for both
  engines (decided 2026-11-02). `feat/phylowgs-go-on-rc-2.0.0` lost `--include-multiprimary`
  in the merge; restore it.
- `summ.params` is `{}`: the pipeline never passes a params file to `multievolve`, so
  production files carry `{}` too.
- Pipeline work goes on `feat/phylowgs-go-on-rc-2.0.0` (the RC-based branch that was
  run on Iris), not the older `feat/phylowgs-go-comparison` (claude's call).
- Go branch: start from `go-port` with `fix/per-chain-flush` merged first.
- Keep `python | go | both` during validation; `both` keeps emitting `go_*` channels.
  Removal of the Python engine is a later, separate change.
- Per-sample seed: the Go module passes `-r` derived from `meta.id` (stable hash) so
  samples do not share RNG streams (the shared-seed artefact from the K comparison).
- Image: built from `go-port` after `fix/per-chain-flush` is merged; stays on the
  personal Docker Hub repo `docker.io/orgeraj/phylowgs-go` for now (decided 2026-11-02),
  under a new tag so the existing `go-port` tag is untouched. Moving to `ghcr` is later.

## Design notes for the builder

- **Where:** new file(s) in `PhyloWGS_refactor` (e.g. `results_writer.go`), called at
  the end of `writeResults` after `trees.zip`. Keep `trees.zip`/`mutlist.json`/
  `best_tree.json` as they are (the comparison tooling uses them).
- **Munger port:** port `ResultMunger` line-for-line, including its quirks: sorted
  iteration orders, `int(round(min_ssms * n_ssms))` with Python 2 rounding (half away
  from zero), `np.mean` over `cellular_prevalence`, `implied_phi = min(2*(d-a)/d, 1)`
  with `ref_reads`/`total_reads` averaged across samples, strict `<` tie-breaking in
  best-node search over Python dict iteration order of `populations` (int keys —
  CPython 2 small-int dict order is ascending for these sizes; verify against the
  reference), superclone weighted-mean CP, multiprimary tree renumbering. The asserts in
  the munger become errors.
- **JSON types, everywhere Go writes JSON (including `trees.zip` and `mutlist.json`):**
  the reference reads Go's own files with `json.load`, so Go's types flow through it;
  Go must therefore emit Python types at the source, or the reference would inherit
  Go's `1` for `1.0`. Python writes floats with a decimal point or exponent (`1.0`, `1e-05`)
  and ints bare. Go's `encoding/json` writes `float64(1)` as `1`. Use a float type that
  marshals Python-style (`repr` rules: shortest round-trip, `.0` for integral values,
  exponent outside 1e-4 ≤ |x| < 1e16) and keep `num_ssms`/`num_cnvs`/indices-of-nodes as
  ints. Known current mismatch: root `cellular_prevalence` (Go `1`, Python `[1.0]`).
- **Indices:** `IndexCalculator` port; same asserts (as errors).
- **Densities:** port scipy (container's version — record it) `gaussian_kde` for 2-D
  and 1-D: Scott's factor `n^(-1/(d+4))`, `np.cov` (ddof=1) × factor², inverse and
  `sqrt(det(2π·cov))·n` normaliser, sum of `exp(-0.5·energy)`. Fallback order and
  triggers as Python: singular matrix or invalid float op (`np.seterr(invalid='raise')`
  → e.g. negative det under sqrt) → 1-D → zeros. Near-singular cases may diverge from
  LAPACK; the test set must include all-linear trees (BI=0), all-identical trees, and a
  single tree.
- **`--include-ssm-names`:** Python semantics — names only with the flag; with the
  flag, an empty name is written as `""` (Go's `mutlist.json` currently omits it). The
  no-flag golden case is the reference output with `name` removed (the reference
  ignores the flag).
- **mutass.zip:** one `<tidx>.json` per remaining tree, deflate, `{mut_assignments,
  dataset_name}`; node keys as strings.
- **Pipeline module:** `PHYLOWGS_WRITERESULTSGO` switches to the Go image and runs
  `phylowgs-go write-results <args> <prefix> trees.zip mutlist.json <summ> <muts>
  <mutass>` (same argument order as `write_results.py`). `bin/write_results_go.py`
  leaves the pipeline; it stays as the test reference (copied into the Go repo's test
  harness). Per-sample seed is passed via `PHYLOWGS_MULTIEVOLVEGO`'s `ext.args`.

## Verification

1. Unit tests (Go) per munger operation on small hand-built trees, each with a
   recorded kill (mutate the port, test fails).
2. Golden-file equivalence: synthetic fixtures only (no patient data in the repo):
   run Go → `trees.zip`; run `write_results_go.py` in the `phylowgs:v1.5-msk` container
   with each flag combination; commit decoded reference outputs; Go test compares.
   Cases: default filters, `--include-multiprimary`, `--keep-superclones`, `--min-ssms`
   as count and fraction, multiprimary-threshold raise, all-linear trees, single tree,
   CNV-only nodes.
3. Real-data equivalence on Iris (outputs stay on Iris; access only via `hpcgate`,
   code via git push/clone): ~20 NeoQual samples across size classes, using the Go
   `trees.zip` from the running NeoQual-1500 run (`go_runs/nq1500`, binary cf49940);
   Go-native vs reference on the same trees; report counts of mismatching files/fields
   only.
4. End-to-end: `--phylowgs_engine go` pipeline run on 3–5 samples whose Python
   pipeline outputs already exist (NeoQual `MASTER_RESULTS`); compare the final
   `*_input.json` structure (keys, types, nesting) against those.

## Left to the builder

File layout inside the Go repo, the float-marshalling approach, test-harness scripts,
and whether `write_results_go.py` stays in the pipeline repo after goal 4.
