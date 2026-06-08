# Tree Selection & Reporting Notes

**Scope:** This is a written recommendation only. It proposes **zero changes to the
sampler.** Nothing here touches likelihoods, MH proposals, acceptance ratios,
stick-breaking, or hyperparameter sampling. Every option below operates strictly
on the *post-sampling reporting layer* (which already-sampled tree we surface to
the user), which is identical math regardless of which tree we print.

**Motivating context:** The earlier "Go under-recovers K at high K vs Python"
finding was traced to a *validation-harness* artifact (wrong baseline seed set,
inconsistent population counting, and multimodal best-of-N luck), not a sampler
bug. The Go sampler is a verified-faithful line-by-line port. See
`MEMORY.md` / `sim_validation/BUGS_FOUND.md`. This doc addresses the third of
those three issues — multimodal mixing making the *reported* tree luck-sensitive
— purely at the reporting layer.

---

## 1. How Go currently picks the reported tree

The reported tree is `best_tree.json`, produced in `writeResults`
(`main.go:2836`):

1. **Collect snapshots.** During sampling, every retained post-burnin sample is
   captured by `snapshotTree` (`main.go:3979`), which serializes the live TSSB
   state (after `summarizePops` + `postProcessSummary` to drop empty pops) into a
   self-contained JSON line. These are written per chain as
   `chain_<id>_trees.ndjson` and into `trees.zip`/`mutass.zip`. The per-sample LLH
   is stored in lockstep in `SampleLLH`.

2. **Filter chains.** `filterChainsByInclusion` keeps chains whose
   `logsumexp(tree LLHs)` is within `chainInclusionFactor` of the best chain —
   a faithful port of Python `multievolve.determine_chains_to_merge`
   (`multievolve.py:215`). Default factor 1.1.

3. **Pick the single global argmax.** `pickBestSample` (`main.go:2809`) scans all
   samples in the included chains and returns the **one sample with the highest
   LLH** (first-seen wins on ties, matching Python's sequential argmax).

4. **Post-process that one tree.** The winning snapshot is run through
   `removeSmallNodes` (`main.go:3644`) with adaptive threshold
   `max(1%, 2/M)` of total SSMs, then `removeSuperclones`, then written as
   `best_tree.json`. **The `num_populations` reported in `best_tree.json` is
   counted AFTER this merge step.**

### Why this is luck-sensitive on multimodal fixtures

- **It is a point estimate from a single MCMC draw.** The reported topology is
  whatever one sample happened to have the maximum instantaneous LLH. On a
  unimodal posterior this is fine — all the high-LLH samples agree. On a
  **multimodal** posterior (typical at high K), the chains split across two or
  more basins, and the global LLH-argmax can land in any of them depending on
  seeds. Two runs with different seeds can report structurally different trees
  with near-identical LLH.

- **Best-of-N-chains amplifies seed dependence.** Roughly ~40% of chains on the
  high-K fixtures get stuck in a worse mode (established this session). Because we
  report the *single* best sample across the included chains, the reported tree is
  governed by "did *any* of my N chains escape into the good mode?" — a function
  of the seed set, not of the sampler's correctness. Add more lucky chains and the
  reported K creeps up; an unlucky seed set reports lower K. This is exactly the
  variance the validation harness mistook for under-recovery.

- **LLH-argmax ≠ posterior mode.** The single most-likely *individual sample* is
  not necessarily representative of where the posterior mass actually concentrates.
  A narrow, slightly-higher-LLH spike can outrank a broad basin that holds far more
  posterior probability. This is a generic property of reporting a max instead of a
  density summary — it is not a Go-specific defect.

This luck-sensitivity is **a property of the reporting choice (single best-LLH),
not of the sampler.** Python has the identical property if you reduce its output
to a single best-LLH tree.

---

## 2. What Python / PhyloWGS does

The upstream PhyloWGS pipeline does **not** present a single best-LLH tree as its
primary output. Two relevant pieces:

- **`posterior_trees.py` (`compute_lineages`, line 17).** Reads the merged
  `trees.zip`, computes a **topology signature** for every MCMC sample
  (`descend` + `_sort` + `sort_and_merge`, lines 31–46), groups samples by that
  signature into `G`, and ranks groups by **posterior probability**
  `len(group)/num_samples` (lines 53–61, via a max-heap on `-prob`). It then emits
  one figure per topology group, **ordered by posterior probability, explicitly
  not by likelihood** (see the comment at lines 59–60). Each group's clonal
  frequencies are **aggregated (mean ± sd) across all samples in the group**
  (lines 82–99, 225–228). So the canonical Python answer to "what tree?" is *a
  ranked list of topologies with their posterior support*, with per-node
  frequencies pooled over the whole group — a density summary, not a single draw.

- **`write_results.py` / `munge_results.py`.** The JSON the witness viewer
  consumes summarizes the **whole sampled population** of trees (with
  `remove_small_nodes` applied per the same `min_ssms` threshold), again a
  distribution rather than a lone argmax.

**The Go port already has this capability.** `posterior_trees.go` is a faithful
port: `treeSignature` (line 42) mirrors `descend`/`_sort`/`sort_and_merge`
exactly (the empty-trailing-token and empty-intermediate-node subtleties are
reproduced and commented), `groupTreesBySignature` + `rankGroups` (lines 193,
203) reproduce the posterior-probability ranking and tie-break, and
`aggregateFreqsByNode` (line 474) pools clonal frequencies across a group. It is
exposed as the `posterior-trees` subcommand (`main.go:4136`,
`posteriorTreesSubcommand` at `posterior_trees.go:453`), operating on the
`trees.zip` that `writeResults` already emits.

So Go *produces* the same posterior-grouped output Python does — it just doesn't
*headline* it. The headline (`best_tree.json`) is the single-best-LLH point
estimate, which is the source of the variance.

---

## 3. Faithful options to reduce reporting variance

None of these change the sampler. They change *which already-sampled tree(s) we
surface*. Listed with pros/cons.

### (a) Headline a posterior / consensus summary instead of single best-LLH

Make the primary reported artifact the top posterior-probability topology group
from `posterior-trees` (the Python-canonical answer), or a consensus over the
included chains, rather than `best_tree.json`.

- **Pros:** Directly matches Python's actual primary output. Far less seed-/luck-
  sensitive: a topology supported by, say, 35% of post-burnin samples is stable
  across seed sets in a way a single argmax draw is not. The machinery already
  exists (`posterior_trees.go`) and is verified faithful, so this is a reporting
  wiring change, not new math. Node frequencies come with mean ± sd, conveying
  uncertainty.
- **Cons:** "Top posterior topology" still depends on whether chains explored the
  good mode at all — if 60% of mass is stuck in a bad basin, the modal topology
  is the bad one. (Mitigated by combining with (b).) The posterior-trees signature
  groups on *topology + mutation membership*, which can fragment when the same
  biological clone tree differs by a few mutation reassignments; many tiny groups
  can dilute the top group's probability. A point estimate is also simpler for
  some downstream consumers that expect one tree; this changes the contract of
  `best_tree.json`.

### (b) Run more chains

Increase chain count (and/or iterations) so the probability that *no* chain
escapes the bad mode shrinks.

- **Pros:** Trivially faithful — identical sampler, just more independent runs.
  Reduces the best-of-N luck directly: if a single chain escapes with prob ~0.6,
  P(all of N stuck) = 0.4^N falls off fast. Helps *every* reporting option
  including the current one. Embarrassingly parallel.
- **Cons:** Linear compute cost. Does not *fix* the variance of the reporting rule
  — it just buys down the failure probability; an unlucky small-N run still
  misreports. Diminishing returns once the good mode is reliably hit. Doesn't make
  Go match a *specific* Python run; it makes both converge to the same answer in
  expectation.

### (c) Report the distribution of npop across chains, not a point estimate

Alongside (or instead of) a single `num_populations`, report the per-chain
distribution — e.g. min/median/max and the histogram of non-empty population
counts across included chains (and optionally across post-burnin samples).

- **Pros:** Honest about multimodality and the cheapest to add — it is pure
  aggregation over data already in `chain_*_trees.ndjson` / `SampleLLH`. Makes the
  Go-vs-Python comparison *valid by construction*: you compare distributions on
  matched seed sets instead of comparing two single luck-draws. Surfaces the
  ~40%-stuck phenomenon explicitly instead of hiding it behind one number.
- **Cons:** Not a "the tree" answer — it's diagnostic, so it complements rather
  than replaces (a). Requires the consumer/validation harness to compare
  distributions, not scalars. Care needed to count populations *consistently*
  between Go and Python (raw non-empty vs post-`removeSmallNodes` — see the known
  counting-convention mismatch); whatever convention is chosen must be applied
  identically on both sides.

### Cross-cutting fix that costs nothing: count populations consistently

Independent of (a)–(c): the headline `num_populations` is currently counted
**after** `removeSmallNodes` merging on the Go side, while the earlier Python
baseline counted **raw non-empty** populations. Any comparison must apply one
convention to both. This is a measurement fix, not a sampler change, and it alone
closed most of the apparent Go-vs-Python gap this session (Go ~5.3 ≈ Python ~5.7
at K=10 when matched).

---

## 4. Recommendation

**Primary recommendation: (a) + (c) together, with (b) as a tunable knob.**

1. **Surface the posterior-grouped summary as the primary answer** (option a),
   using the already-faithful `posterior-trees` subcommand. Keep `best_tree.json`
   for backward compatibility but document it as "single highest-LLH draw
   (luck-sensitive at high K); see posterior_trees/ for the supported-topology
   ranking." This aligns Go's *headline* with what Python actually headlines and
   removes single-draw variance from the primary output. **No sampler change.**

2. **Report the npop distribution** (option c) as a diagnostic in `summary.json`
   (e.g. per-chain non-empty population counts + median), computed from existing
   snapshots. This makes the Go-vs-Python comparison valid on matched seeds and
   exposes multimodality directly. **No sampler change.**

3. **Treat chain count (option b) as the variance knob** for production runs on
   the cluster, not as a code change — more chains buy down best-of-N luck for
   whichever reporting rule is used.

Crucially, fix the **population-counting convention** so both engines are counted
identically before any comparison is drawn.

### What validation evidence would confirm this (to gather on the cluster)

Run on the matched seed set (e.g. Go `sim_validation_juno_529` seed-set B vs
Python `junoresultsmulti` seed-set B, 216/216 same seeds):

- **Consistency check:** For each K fixture, count non-empty populations with the
  **same** convention on both engines (raw non-empty *and*, separately,
  post-`removeSmallNodes`). Expect Go ≈ Python within noise on both
  conventions (confirming no sampler-driven gap).
- **Posterior-group agreement:** Run `phylowgs-go posterior-trees` and Python
  `posterior_trees.py` on the same merged `trees.zip` and compare the ranked
  topology lists and their posterior probabilities. Expect matching top groups
  and probabilities within Monte-Carlo error (this directly validates the
  `posterior_trees.go` port end-to-end and confirms (a) gives a stable answer).
- **Variance reduction:** For a fixed fixture, report the spread of npop /
  top-topology across seed sets under (i) single best-LLH and (ii) top posterior
  group, at N and 2N chains. Expect (ii) to show materially lower seed-to-seed
  spread than (i), and the spread to shrink as N grows — quantifying how much
  reporting variance (a)+(b) remove without touching the sampler.

A passing result is: matched-seed Go and Python agree on npop distributions and
top posterior topologies, and the posterior-summary headline is visibly more
stable across seeds than the single best-LLH headline — confirming the original
"under-recovery" was reporting variance, now reduced, with the sampler untouched.
