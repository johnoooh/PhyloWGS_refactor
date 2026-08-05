# Go–Python Equivalence Audit — 2026-06-24

Multi-agent adversarial audit of whether the Go port (`go-port` branch:
`main.go`, `posterior_trees.go`, `archive.go`) faithfully reproduces the
reference Python/C++ implementation (`main` branch: `tssb.py`, `alleles.py`,
`node.py`, `params.py`, `util.py`/`util2.py`, `write_results.py`, `evolve.py`,
`multievolve.py`, `posterior_trees.py`, `pwgsresults/*`; original C++
`mh.hpp`/`mh.cpp`, `util.cpp`).

Method: 12 Go↔reference module pairs reviewed in parallel at high reasoning
effort, every claimed divergence independently re-verified by a skeptic agent
defaulting to *refute*, then synthesized. 34 candidates → **13 confirmed, 21
refuted**. Companion to `CODE_REVIEW.md`.

## 1. Bottom line

The **core MCMC sampler is faithful**. TSSB construction & stick-breaking, tree
mutation/ancestor bookkeeping, datum-assignment resampling, hyperparameter
resampling, math primitives/RNG, and the chain drivers were checked line-by-line
and confirmed equivalent. **No critical divergences survived verification** —
nothing changes the sampled distribution, likelihood, or accept/reject logic on
the default (CPU, stub-CUDA) build for ordinary inputs.

Confirmed divergences cluster in two non-sampler areas:
1. one genuine **edge-case** likelihood difference for an SSM co-located on its
   own zero-copy-allele CNV, where the Python reference (`data.py`) and the
   original C++ kernel (`mh.hpp`) **disagree with each other** — Go matches Python;
2. **post-MCMC summarization / output serialization** differences affecting
   reported subclone counts, population labeling, and JSON schema — not inference.

This resolves the long-standing "MH proposal divergence is the largest unknown"
flag from `CODE_REVIEW.md`: see M1.

## 2. Confirmed divergences

### Critical
None.

### Major

**M1 — Co-located CNV-SSM likelihood: Go/Python drop `nv==0` timing states + use
`log(1/n_valid)` prior; C++ `mh.hpp` keeps all 4 with fixed `log(0.25)`.**
- Go: `main.go` fast MH path (`logLikelihoodWithCNVTreeMHPrecomputed`) and slow
  path (`logLikelihoodWithCNVTree`) filter channels on `nv>0`, prior
  `math.Log(1.0/n)`, logsumexp over survivors.
- Ref: `phylowgs/mh.hpp` `log_complete_ll` keeps 4 terms each guarded
  `if(nr+nv>0) else log(1e-99)`, fixed `log(0.25)` prior; `params.py:159-166`
  applies static state corrections before invoking the compiled `mh.o`.
- Consequence: for an SSM co-located on its own CNV (`len(poss_n_genomes)==4`)
  where an aggregated timing state has `nv_k==0` but `nr_k>0` (LOH `minor_cn=0`
  or homozygous deletion `total_cn=0`), Go follows Python's `data.py`
  drop+renormalize while the C++ kernel keeps the term. ≈0.29 nats per affected
  datum per timepoint; the channel mix changes between old/new proposed states so
  it does **not** cancel in `newLLH − oldLLH` → shifts MH acceptance for that
  SSM's φ. **Go faithfully matches Python; Python and C++ genuinely disagree.**
- Status: **not fixed — reference-authority decision required.** The reference
  *pipeline* compiles and runs `mh.o`, so for true pipeline parity the Go MH
  kernel should mirror `mh.hpp` (4 terms keyed on `nr+nv>0`, fixed `log(0.25)`,
  plus the `params.py:159-166` static pre-collapse) rather than the dynamic
  `data.py` filter. Also: the comment at `main.go:1834` ("Mirrors Python's
  mh.cpp log_complete_ll inner loop") is misleading — it mirrors `data.py`.

**M2 — GPU complete-data LLH path dropped all CNV-datum contributions.** ✅ FIXED
- Was: `completeDataLogLikelihood` GPU branch returned `gpuLLH` before the
  CNV-datum block, omitting each CNV datum's mixture-weight + binomial terms when
  `useGPU==true`. Unreachable on default builds (stub CUDA forces `useGPU=false`).
- Fix: extracted `cnvDatumLogLikelihood(weights, nodes)` and called it in both
  the GPU and CPU branches so the totals agree. CPU output unchanged by
  construction (pure extraction).

**M3 — Adaptive min-SSM pruning threshold replaces Python's fixed 1%.**
- Go: `minFrac = max(0.01, 2/totalSSMs)`, `threshold = round(minFrac*total)` floored ≥1.
- Ref: `result_munger.py:191-205` uses `int(round(0.01*num_ssms))`, no floor.
- Consequence: for M<200 SSMs Go prunes more subclones → fewer reported
  populations; converges at M≥200. Output-summarization only; documented as
  intentional in the Go comment.
- Status: **not fixed — intentional deviation.** Drop the adaptive bump + ≥1
  floor if strict reference parity is wanted; otherwise document in user docs.

**M4 — Empty internal-node removal does not re-sort reparented children by φ.**
- Go: `postProcessSummary` reparents grandchildren and renumbers by ascending old
  index; the faithful `removeEmptyNodes` (`main.go:3064`) is **dead code**.
- Ref: `util2.py:158-186` removes empty nodes on the live tree, then
  `result_generator.py:81-86` numbers populations by φ-sorted preorder.
- Consequence: with an empty internal node having data-bearing descendants, Go
  and Python assign different population indices / `structure` ordering (CPs and
  mutation sets identical). Summary labeling only.
- Status: **not fixed — moderate scope.** Number populations *after* empty-node
  removal using Python's φ-descending DFS (or wire the existing
  `removeEmptyNodes` into the live tree before `summarizePops`).

**M5 — `mutlist.json physical_cnvs[].cell_prev` was a raw string, not floats.** ✅ FIXED
- Was: `PhysicalCNV.CellPrev string` stored `"0.0|0.718"` raw.
- Ref: `util2.py:71` → `[float(C) for C in cell_prev.split('|')]` → JSON array.
- Fix: `CellPrev []float64`; `parsePhysicalCNVs` splits on `|` and parses floats.

**M6 — Posterior aggregation: empty intermediate nodes removed before signature.**
- Go: `posterior_trees.go` `treeSignature` consumes `postProcessSummary`-cleaned
  snapshots; Ref `posterior_trees.py:27-50` descends the raw `tssb.root`.
- Consequence: for `P(a)->E(empty)->C(b)`, Python signature `_a_;_b_;` vs Go
  `_a_;_a_b_;` → changes topology grouping → posterior probabilities and top-tree
  ranking. Posterior-trees report only, not the sampler.
- Status: **not fixed.** Compute `treeSignature` on the raw pre-`postProcessSummary`
  snapshot, replicating Python's genotype-inheritance break at empty nodes.

**M7 — Posterior frequency aggregation keys by node ID, not mutation content.**
- Go: `posterior_trees.go` `aggregateFreqsByNode` keyed by node ID; Ref
  `posterior_trees.py:88-107` keys by `;`-joined `datum.name`.
- Consequence: same-signature trees whose siblings sort differently by φ get
  mismatched node-ID→mutation mappings → wrong mean±std in the posterior table.
  Reporting only.
- Status: **not fixed.** Key by the node's mutation-name set, matching Python.

### Minor

- **m1 — move-to-best initial delta:** Go uses `+Inf` (always assigns nearest);
  Python uses `1.0` with strict `<` (leaves `best_node=None` → crash). Differs
  only at exactly `phi_delta==1.0`. Go is arguably more robust. No action.
- **m2 — fraction-threshold rounding:** Go `math.Round` (half-away) vs Python3
  banker's rounding. One-SSM difference only at exact-half parity (verified at
  M=250: Go 3 vs Py 2). Fix with round-half-to-even only if exact parity needed.
- **m3 — `mutlist.json` always emits SSM `name`:** Python gates behind
  `include_ssm_names` (default off, treated as possibly sensitive); Go always
  emits (with `omitempty`). **Not fixed** — defaulting off risks breaking
  downstream viewers; flagged for an explicit decision. Note the data-sensitivity
  angle if names are ever gene/patient-linked.
- **m4 — posterior probability not rounded to 4 dp before ranking/display:**
  cosmetic (`.tex` filename + displayed string); ranking effect only at
  ns>10000 (never hit at typical 1000–5000 samples). Fix: round before keying.

### Benign

- **b1 — `mutlist.json` missing top-level `dataset_name`.** ✅ FIXED — the mutass
  archive writer already included it (`archive.go:102`); `writeMutList` now
  takes `datasetName` and emits the key.

## 3. Module-by-module faithfulness

| Module | Verdict | Notes |
|---|---|---|
| TSSB construction & stick-breaking | Faithful | dirichlet pseudocount, depth-0 pin verified |
| Tree mutation & ancestor bookkeeping | Faithful (line-by-line) | `cullTree`/`trim_zeros`, `killNode` π-return, ancestor self-inclusion |
| Datum assignment resampling | Faithful | remove-then-add ordering matches `tssb.py:124-128` |
| Hyperparameter resampling | Faithful | |
| Math primitives & RNG | Faithful | |
| MCMC chain drivers & deep-copy | Faithful | SSM.CNVs + per-chain `*CNV` deep-copy preserved |
| SSM likelihood & CNV genome-count | Faithful to Python; diverges from C++ | M1 edge case |
| Metropolis-Hastings params | Mostly faithful; edge-case only | common case identical (`log(2)+log(0.25)=log(0.5)`) |
| Complete-data LLH | CPU faithful; GPU fixed | M2 (now fixed) |
| Post-MCMC summarization & pruning | Divergent (labels/counts) | M3, M4, m1, m2 — not sampling |
| Data loading & output serialization | Schema fixes applied | M5, b1 fixed; m3 open |
| Posterior tree aggregation | Divergent (reporting) | M6, M7, m4 — not sampling |

## 4. Fixes applied this pass

CPU/inference behavior unchanged; all three are output-schema or
unreachable-path corrections, each pure-extraction or additive.

- **M5** — `cell_prev` parsed to `[]float64` (`main.go` `PhysicalCNV`, `parsePhysicalCNVs`).
- **b1** — `dataset_name` added to `mutlist.json` (`writeMutList` signature + caller).
- **M2** — CNV-datum contribution extracted to `cnvDatumLogLikelihood` and added to the GPU branch.
- Tests: `equivalence_fixes_test.go` (`TestParsePhysicalCNVs_CellPrevFloats`,
  `TestWriteMutList_CellPrevArrayAndDatasetName`). Full suite green.

## 5. Residual unknowns — recommend runtime numeric cross-checks

1. **MH proposal density (highest priority).** This audit confirmed the
   *likelihood* terms but did not exhaustively trace the proposal/Hastings
   correction against `mh.cpp`. Feed identical fixed states/seeds/proposals
   through Go `metropolis` and compiled `mh.o`; compare per-step
   `log q(old|new) − log q(new|old)` and the acceptance ratio. A proposal
   asymmetry would be a silent critical bug static reading can miss.
2. **M1 magnitude in practice.** Run a fixture with an LOH (`minor_cn=0`) or
   homozygous-deletion CNV co-located with an SSM through both pipelines; diff
   per-SSM φ posteriors to decide fix-vs-document.
3. **GPU path (M2).** On a CUDA build, validate `completeDataLogLikelihood` on a
   CNV-bearing dataset against the CPU path (should now agree up to FP ordering).
4. **End-to-end determinism on matched seeds.** Per the May-29 seed/baseline
   lesson, compare Go vs Python complete-data LLH trajectories step-by-step on
   the *same seed* for the first N iterations, not just final population counts.
5. **Rounding edge cases (m2, m4)** — only if exact output parity is a release
   requirement.
