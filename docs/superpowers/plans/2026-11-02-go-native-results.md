# Go-native PhyloWGS results Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `phylowgs-go write-results` writes `summ.json.gz` / `muts.json.gz` / `mutass.zip` identical (decoded) to Python's `write_results.py`, and the neoantigen pipeline's Go engine uses it instead of the python2 adapter.

**Architecture:** A faithful, typed Go port of `pwgsresults` (`ResultMunger`, `IndexCalculator`, `JsonWriter` incl. scipy `gaussian_kde`) in new files of `package main`, exposed as a subcommand that reads Go's existing `trees.zip` + `mutlist.json`. Go's JSON outputs switch to Python number typing (`1.0` not `1`) so the Python reference (`write_results_go.py`) inherits the right types. Equivalence is proven with committed golden files produced by the Python reference in its own container, then on real Go trees on Iris, then end-to-end in the pipeline.

**Tech Stack:** Go 1.21 stdlib only (no new modules); python2 reference in `ghcr.io/mskcc-omics-workflows/phylowgs:v1.5-msk`; Nextflow DSL2.

**Spec:** `docs/superpowers/specs/2026-11-02-go-native-results-design.md`

## Global Constraints

- Go repo worktree: branch `feat/go-native-results` (= `go-port` + `fix/per-chain-flush` + spec). Pipeline: new branch `feat/phylowgs-go-native` off `origin/feat/phylowgs-go-on-rc-2.0.0` in `~/Documents/GitHub/neoantigen-pipeline`.
- `go.mod` stays dependency-free. New files ≤ 300 lines, functions ≤ 50 lines (user rule; `main.go` is legacy, only touch what a task names).
- Production filters, both engines: `--max-multiprimary 1.0 --include-multiprimary` (plus `--include-ssm-names`, added by the module script).
- Defaults identical to `write_results.py`: `--min-ssms 0.01`, `--max-multiprimary 0.8`, `--include-multiprimary` off, `--keep-superclones` off, `--include-ssm-names` off.
- `summ.params` is always `{}`.
- Python reference = `write_results_go.py` (copied from the pipeline repo `bin/`) run in `ghcr.io/mskcc-omics-workflows/phylowgs:v1.5-msk` (`docker run --platform linux/amd64`; Iris fallback: `/data1/core005/ccs/orgeraj/go_runs/m1xcheck/phylowgs-1.5.2-msk.sif`).
- Equality = decoded JSON: same keys, same JSON number class (int literal vs float literal), exact values; `tree_densities` values relative tolerance 1e-9; `mut_assignments` SSM/CNV lists compared as sorted lists.
- No patient data in either repo: committed fixtures are synthetic only. Real-data outputs stay on Iris; print counts only. Iris only via `hpcgate` (`run`, `submit`); if a submit is blocked, hand the exact command to the user.
- Image: `docker.io/orgeraj/phylowgs-go:go-native-<shortsha>`; the existing `go-port` tag is not touched.
- Commits end with `Claude-Session: https://claude.ai/code/session_01LRKp1yBzdDy6GKf6o3dCoS`; before every commit grep the staged diff for `P-[0-9]{7}|C-[0-9A-Z]{6}|s_C_`.

## Review Focus

1. Best-node search finds no node (every `|mean(cp) − implied_phi| ≥ 1`) — Python dies with `KeyError: None`; Go must return an error, not drop the mutations. Test in Task 4.
2. `--min-ssms` fraction that rounds to 0 (0.01 × 40 SSMs) keeps every node, and `--min-ssms 1` is a count, not 100 %. Test in Task 4.
3. A tree with ≤ 1 SSM in total — Python's index division raises `ZeroDivisionError`; Go must error with the tree index. Test in Task 5.
4. `trees.zip` with no `tree_*` entries or a missing `mutlist.json` — error, never empty output files. Test in Task 3.
5. CNV-only node next to the clonal node (`num_ssms == 0`) — superclone check must skip it (Python's division guard), and CNVs moved by small-node removal must land in `cnvs`. Golden case `cnvonly` (Task 2) + test in Task 4.

---

### Task 1: Python number typing in Go's JSON

**Files:**
- Create: `pyjson.go`, `pyjson_test.go`
- Modify: `main.go` — `summarizePops` (`popSummary.CellPrev`, `llh`), `snapshotTree` (re-type after `postProcessSummary`), the `best_tree.json` marshal in `writeResults`, `writeMutList` (`ExpectedRefIn*`, `Name`)

**Interfaces:**
- Produces: `type PyFloat float64` with `MarshalJSON() ([]byte, error)`; `func pyFloats(xs []float64) []PyFloat`; `func retypeSummary(generic map[string]interface{})` (wraps `llh` and every `populations[*].cellular_prevalence` element in `PyFloat`, in place).

- [ ] **Step 1: Write the failing tests** in `pyjson_test.go`:
  - `TestPyFloat_Marshal`: table — `1.0→"1.0"`, `0→"0.0"`, `-1300→"-1300.0"`, `0.5→"0.5"`, `1e-05→"1e-05"`, `0.0001→"0.0001"`, `1e16→"1e+16"`, `123456789012345.0→"123456789012345.0"`, `0.1+0.2→"0.30000000000000004"`, `NaN→"NaN"`, `+Inf→"Infinity"`, `-Inf→"-Infinity"` (Python 2.7 `json.dumps(float)` = `repr`: shortest round-trip; exponent form iff exponent < −4 or ≥ 16, written `1e-05`/`1e+16`; otherwise fixed with at least one decimal).
  - `TestSnapshotTree_PythonTypes`: build a 1-node-plus-root TSSB (reuse the helpers in `cnv_treedata_test.go`), `snapshotTree(...)` with `llh = -1300`; assert the raw JSON contains `"cellular_prevalence":[1.0]` for population `"0"` and `"llh":-1300.0`, and `"num_ssms":` values are bare ints.
  - `TestWriteMutList_PythonTypes`: an SSM with `MuR = 1`, `Name = ""` → `mutlist.json` has `"expected_ref_in_ref": 1.0` and `"name": ""` (no `omitempty`).
- [ ] **Step 2:** `go test -run 'PyFloat|PythonTypes' ./...` → FAIL (undefined `PyFloat`).
- [ ] **Step 3: Implement** `PyFloat.MarshalJSON` with `strconv.FormatFloat(f, 'r'-equivalent)`: use `'g', -1` to get shortest digits, then decide form from the decimal exponent per the rule above; append `.0` when fixed form has no `.`. Apply `pyFloats` in `summarizePops`, call `retypeSummary` in `snapshotTree` right before its final `json.Marshal`, and on the `best_tree.json` `generic` map before marshalling. `writeMutList`: `ExpectedRefInRef/InVariant PyFloat`, `Name string \`json:"name"\``.
- [ ] **Step 4:** `go test ./...` → all PASS (existing byte-equality tests such as `chain_flush_test.go` compare Go to Go and must still pass), `go vet ./...` and `gofmt -l .` clean.
- [ ] **Step 5: Commit** `feat: write Python number types (1.0 not 1) in trees.zip, best_tree.json, mutlist.json`.

### Task 2: Golden fixtures, Python references, and the comparator

**Files:**
- Create: `testdata/native_results/make_cases.py` (python3, builds hand-made cases), `testdata/native_results/regen.sh`, `testdata/native_results/write_results_go.py` (verbatim copy of the pipeline repo `bin/write_results_go.py`), `testdata/native_results/<case>/{trees.zip,mutlist.json}`, `testdata/native_results/<case>/<flags>/expected/{summ.json,muts.json,mutass.json}` or `expected/error.txt`, `testdata/native_results/SCIPY_VERSION`
- Create: `results_compare_test.go`

**Interfaces:**
- Produces: `func comparePythonJSON(t *testing.T, label string, got, want []byte)` — decodes both with `json.Decoder.UseNumber()`, walks recursively, reports the first ≤ 20 differences by JSON path: missing/extra keys, type differences, number class (literal contains `.`, `e`, `E`, `NaN`, `Infinity` ⇒ float, else int), exact value inequality; under a path whose first key is `tree_densities` uses relative tolerance 1e-9; any list under `.../mut_assignments/<node>/ssms|cnvs` is sorted before comparing. `mutass.json` is the decoded `mutass.zip` as `{"<tidx>": <entry>}`.
- Produces: fixture layout above; `regen.sh` is the only way expected files change.

- [ ] **Step 1: Write the comparator's own test** `TestComparePythonJSON_Kills` (uses a fake `testing.TB` recorder): `{"a":1}` vs `{"a":1.0}` → reports number-class diff at `a`; `{"tree_densities":{"0":1.0000000001}}` vs `{"tree_densities":{"0":1.0}}` → no diff; same pair under key `x` → diff; `mut_assignments/1/ssms ["s1","s0"]` vs `["s0","s1"]` → no diff; missing key → diff. Run → FAIL (undefined), implement, run → PASS.
- [ ] **Step 2: Build the cases.**
  - `sim`: synthetic sample from `sim_validation/generate_fixtures.py` (K=4, ~200 SSMs, ~30 CNVs, one sample); run the Task 1 binary `go run . -B 50 -s 200 -j 2 -r 7 -O <tmp> -D sim <ssm> <cnv>`; copy `trees.zip`, `mutlist.json`.
  - `make_cases.py` writes Go-format `trees.zip` entries (`tree_<idx>_<llh>`, JSON `{llh, structure, populations, mut_assignments}`, Python-typed floats) + `mutlist.json` for: `linear` (5 linear trees, branching index 0 → KDE 1-D path), `identical` (3 identical trees → zeros path), `single` (1 tree), `superclone` (clonal node 1 with 10 SSMs, cp 0.95; only child 2 with 40 SSMs, cp 0.9; must merge), `cnvonly` (node with 0 SSMs + 2 CNVs as child of clonal, plus a node below `--min-ssms` whose SSMs/CNVs move), `multiprimary` (4 of 5 trees with two root children).
- [ ] **Step 3: Write `regen.sh`**: for each case × flag set, run `write_results_go.py <flags> --include-ssm-names <case> trees.zip mutlist.json …` in the reference container, gunzip/unzip, dump each file with `json.dumps(obj, sort_keys=True, indent=1)` **inside the py2 container** (keeps Python's float literals), write `expected/`; on non-zero exit write stderr's last line to `expected/error.txt`. Flag sets: `default` (none), `prod` (`--max-multiprimary 1.0 --include-multiprimary`), `keepsc` (`--keep-superclones`), `mincount` (`--min-ssms 5`), `nonames` (default flags; expected = `default` with every `ssms/*/name` removed — the reference ignores the flag). Run `multiprimary` with `default` (expects the raise) and `prod`. All other cases run `default` + `prod`; `sim` runs all five. Record `python2 -c 'import scipy,numpy;print scipy.__version__, numpy.__version__'` into `SCIPY_VERSION`.
- [ ] **Step 4: Run `regen.sh`**; check `expected/` exists for every case/flag pair, `multiprimary/default/expected/error.txt` contains `of trees are multiprimary`, and `linear`/`identical` show the expected density shapes (`identical` all `0.0`).
- [ ] **Step 5: Commit** fixtures, script, comparator (`test: Python-reference golden files for write-results`).

### Task 3: Load Go results

**Files:**
- Create: `results_model.go`, `results_model_test.go`

**Interfaces:**
- Produces:
  ```go
  type Population struct { CellPrev []float64; NumSSMs, NumCNVs int }
  type MutAss struct { SSMs, CNVs []string }
  type TreeSummary struct { LLH float64; Structure map[int][]int; Populations map[int]*Population }
  type MutListSSM struct { RefReads, TotalReads []int; ExpRefInRef, ExpRefInVariant float64; Name string }
  type MutList struct { SSMs map[string]*MutListSSM; CNVs map[string]json.RawMessage } // CNVs passed through untouched
  type ResultSet struct { Trees map[int]*TreeSummary; MutAss map[int]map[int]*MutAss; Muts *MutList }
  func loadGoResults(treesZip, mutlistPath string) (*ResultSet, error)
  ```
  CNV entries are pass-through `RawMessage` except `ref_reads`/`total_reads` needed for best-node search: keep a parallel `CNVReads map[string][2][]int` on `MutList`.

- [ ] **Step 1: Failing tests:** `TestLoadGoResults_Sim` (loads `testdata/native_results/sim`, tree count = number of `tree_*` entries, every tree has population 0, `MutAss` keys = trees keys); `TestLoadGoResults_NoTrees` (zip with only a `burnin_0_…` entry → error containing `no trees`); `TestLoadGoResults_MissingMutlist` → error.
- [ ] **Step 2:** run → FAIL. **Step 3:** implement with `newTreeArchiveReader` (`archive.go`); tree idx from the entry name; `dataset_name` key in `mutlist.json` dropped. **Step 4:** run → PASS. **Step 5: Commit** `feat: load Go trees.zip + mutlist.json into a typed result set`.

### Task 4: Port ResultMunger

**Files:**
- Create: `results_munger.go`, `results_munger_test.go`

**Interfaces:**
- Consumes: Task 3 types.
- Produces: `func (rs *ResultSet) RemoveSmallNodes(minSSMs float64) error`, `func (rs *ResultSet) RemoveSuperclones() error`, `func (rs *ResultSet) RemoveMultiprimaryTrees(maxProp float64) error`.

Port `pwgsresults/result_munger.py` (mskcc/phylowgs) operation-for-operation. Decisions the source does not make obvious:
- Iterate tree indices and node indices ascending everywhere Python iterates a dict (CPython 2 int-key dicts of this size iterate ascending); best-node search: ascending node order, strict `<`, start `lowest = 1`, skip node 0; no node found → error `tree %d: no node within phi 1 of mutation %s`.
- `min_ssms`: `>= 1` → `int(min_ssms)`; else `int(math.Round(min_ssms * float64(len(Muts.SSMs))))` (Python 2 `round` = half away from zero = `math.Round`).
- implied phi: `min(2*(mean(total) − mean(ref))/mean(total), 1.0)` from `Muts` SSM or CNV reads; `mean` = arithmetic mean as float64 summed in index order.
- `_find_small_nodes`' phi-ordering assert, the superclone `clonal_idx == 1` / `child_idx == 2` asserts, and the renumbering asserts become errors naming the tree.
- Superclone: skip if root has ≠ 1 child, clonal has ≠ 1 child, child `num_ssms == 0`, `float(clone)/child > 0.33`, or `mean(|cp_clone − cp_child|) > 0.1`; new clone cp elementwise `(cp1*n1 + cp2*n2)/(n1+n2)`; child's muts appended to node 1 (SSMs then CNVs, in the child's order).
- Multiprimary raise text, exactly: `fmt.Sprintf("%d%% of trees are multiprimary (%d of %d), so not enough to report good posterior.", int(100*frac), n, total)` when `frac >= maxProp`.
- After every removal, `num_ssms`/`num_cnvs` of non-root nodes are recounted from `MutAss` (Python `_correct_mut_counts`).

- [ ] **Step 1: Failing tests** (hand-built `ResultSet`s, assert full resulting `Structure`, `Populations` keys/counts/cp, `MutAss`):
  - `TestRemoveSmallNodes_MovesToClosestPhiAndRenumbers` — 0→1(cp .9, 20 SSMs)→2(cp .5, 1 SSM)→3(cp .3, 10 SSMs), min 5: node 2 removed, 3 renumbered to 2 and reparented to 1, the moved SSM (ref 60 / total 100 → implied phi 0.8) lands on node 1.
  - `TestRemoveSmallNodes_TieGoesToLowerIndex` — same shape with cp .75 / .5 / .25 and the moved SSM at ref 75 / total 100 (implied phi 0.5): deltas are exactly 0.25 and 0.25 → node 1.
  - `TestRemoveSmallNodes_FractionRoundsToZero` (40 SSMs, 0.01 → threshold 0, nothing removed) and `TestRemoveSmallNodes_CountOfOne` (`1.0` = count 1).
  - `TestRemoveSmallNodes_NoBestNode` → error (Review Focus 1): only node 1 survives with cp 0.0 and the moved SSM has implied phi 1.0.
  - `TestRemoveSmallNodes_CNVsMoveToCNVs` (Review Focus 5).
  - `TestRemoveSuperclones_Merges` (counts 10/40, cp .95/.9 → one node, cp 0.91, 50 SSMs) and `_SkipsCNVOnlyChild`, `_SkipsRatioAbove033`, `_SkipsMultiprimary`.
  - `TestRemoveMultiprimary_RenumbersTrees` (trees 0..4, 1 and 3 multiprimary, max 0.8 → keys 0,1,2 = old 0,2,4) and `_RaisesAtThreshold` (exact message for 4 of 5, max 0.8 → `80% … (4 of 5)`).
- [ ] **Step 2:** run → FAIL. **Step 3:** implement. **Step 4:** run → PASS; then kill check: flip the strict `<` to `<=` → `…TieGoesToLowerIndex` fails; revert. **Step 5: Commit** `feat: Go port of pwgsresults ResultMunger`.

### Task 5: Tree indices and densities

**Files:**
- Create: `results_indices.go`, `results_indices_test.go`

**Interfaces:**
- Produces: `func treeIndices(t *TreeSummary) (lin, br, cl float64, err error)`; `func treeDensities(trees map[int]*TreeSummary, idx map[int][3]float64) (map[int]float64, string)` — second value `"2d" | "1d" | "zeros"` (logged by the writer).

Indices: port `index_calculator.py` (relations from root by preorder; linearity = 2·Σ anc→desc nᵢnⱼ / (N(N−1)); branching = Σ over ordered cousin pairs; clustering = Σ nᵢ(nᵢ−1) / (N(N−1))); N ≤ 1 → error `tree %d: fewer than 2 SSMs` (Review Focus 3); the `0 ≤ x ≤ 1` and sum-to-1 (`np.isclose`: rtol 1e-5, atol 1e-8) asserts → errors.

Densities: port scipy `gaussian_kde` at the version in `SCIPY_VERSION` (read that version's `scipy/stats/kde.py` before coding). Points sorted by tree index; X = clustering, Y = branching / (branching + linearity + 1e-4). Bandwidth factor `n^(-1/(d+4))`; covariance = `np.cov` (ddof 1) × factor²; for d = 2, invert and take the determinant through a 2×2 LU with partial pivoting exactly as LAPACK `getrf` (swap rows iff |a21| > |a11|; a zero pivot ⇒ singular); normaliser `sqrt(det(2π·cov)) · n` (or the weighted form of that version — use whichever its source has); evaluation sums `exp(−energy)` over data points in index order. Fall back 2-D → 1-D on singular matrix, NaN/negative under the square root, or any NaN/Inf intermediate (`np.seterr(invalid='raise')`); 1-D → all zeros on the same conditions.

- [ ] **Step 1: Failing tests:** `TestTreeIndices_HandComputed` — 0→1(4 SSMs), 1→2(2), 1→3(2): `lin == 2*(16.0/56)`, `br == 8.0/56`, `cl == 16.0/56` exactly (Python's operation order); `TestTreeIndices_OneSSM` → error; `TestTreeDensities_Paths` — load `linear`, `identical`, `sim` fixtures, compute indices, assert mode `1d`, `zeros`, `2d`, and the densities match `expected/summ.json` `tree_densities` within 1e-9 relative.
- [ ] **Step 2:** run → FAIL. **Step 3:** implement. **Step 4:** run → PASS. **Step 5: Commit** `feat: Go port of tree indices and scipy gaussian_kde tree densities`.

### Task 6: Writer and `write-results` subcommand

**Files:**
- Create: `results_writer.go`, `results_writer_test.go`
- Modify: `main.go` `main()` subcommand dispatch (next to `posterior-trees`)

**Interfaces:**
- Consumes: Tasks 1, 3–5.
- Produces: `func writePythonResults(rs *ResultSet, dataset string, includeSSMNames bool, summPath, mutsPath, mutassPath string) error`; `func writeResultsSubcommand(args []string) error`. CLI, same order and flags as `write_results.py`:
  `phylowgs-go write-results [--include-ssm-names] [--min-ssms F] [--include-multiprimary] [--keep-superclones] [--max-multiprimary F] DATASET TREES_ZIP MUTLIST_JSON SUMM_OUT MUTS_OUT MUTASS_OUT`
  (Go `flag` needs flags before positionals; that matches the module's argument order.) `--max-multiprimary` outside [0, 1] → error.

Order of operations = `write_results.py` `main()`. Output: `summ` = gzip of `{dataset_name, params: {}, trees: {"<tidx>": {llh, structure: {"<node>": [ints]}, populations: {"<node>": {cellular_prevalence: [PyFloat], num_ssms, num_cnvs}}, linearity_index, branching_index, clustering_index}}, tree_densities: {"<tidx>": PyFloat}}`; `muts` = gzip of `{ssms, cnvs, dataset_name}` with `name` only when `--include-ssm-names`; `mutass` = deflate zip, one `<tidx>.json` = `{mut_assignments: {"<node>": {ssms, cnvs}}, dataset_name}` per remaining tree. Every float through `PyFloat`.

- [ ] **Step 1: Failing test** `TestWriteResults_MatchesPythonReference`: for every `testdata/native_results/<case>/<flags>/expected`, run `writeResultsSubcommand` with that flag set into `t.TempDir()`, decode outputs (gunzip; unzip mutass into `{"<tidx>": …}`), `comparePythonJSON` against each expected file; for `error.txt` cases assert an error whose message ends with the file's text. Plus `TestWriteResultsSubcommand_BadArgs` (5 positionals → usage error; `--max-multiprimary 1.5` → error).
- [ ] **Step 2:** run → FAIL. **Step 3:** implement + dispatch. **Step 4:** `go test ./...` → PASS; kill check: change `0.33` to `0.34` in the superclone rule → the `superclone`/`sim` golden fails; revert. `go vet`, `gofmt -l .` clean.
- [ ] **Step 5: Commit** `feat: phylowgs-go write-results, Python-identical summ/muts/mutass`.

### Task 7: Real-data equivalence on Iris

**Files:**
- Create: `sim_validation/native_results/compare_results.py` (python3; same rules as `comparePythonJSON`, plus `--loose-int-float`; prints per sample only `ok` or counts of differing paths by top-level key — never values or keys under `ssms`/`cnvs`), `sim_validation/native_results/iris_check.sh` (sbatch script)

**Interfaces:**
- Consumes: Task 6 binary built on Iris from the pushed branch; `/data1/core005/ccs/orgeraj/go_runs/posterior_cmp/seed_rerun/tasks.tsv` + `runs/k*/chains/{trees.zip,mutlist.json}` (20 real samples, all size classes).

- [ ] **Step 1:** Push the branch; on Iris `hpcgate run "git -C <clone> fetch && git -C <clone> checkout feat/go-native-results && cd <clone> && go build -o /data1/core005/ccs/orgeraj/go_runs/bin/phylowgs-go-native-<sha> ."` (clone: `/data1/core005/ccs/orgeraj/go_runs/gosrc_flush`).
- [ ] **Step 2:** `iris_check.sh`: for each of the 20 runs, `phylowgs-go-native write-results --max-multiprimary 1.0 --include-multiprimary --include-ssm-names k<k> …` and the reference `write_results_go.py` with the same flags in the SIF, into `/data1/core005/ccs/orgeraj/go_runs/native_check/k<k>/{go,py}`, then `compare_results.py --loose-int-float` (these `trees.zip` files predate Task 1, so the reference inherits Go's old int typing; types are covered by Task 6 and Task 9). `#SBATCH --partition=cmobic_cpu`, 16G, 12 h, output to `native_check/check.log`.
- [ ] **Step 3:** Submit with `hpcgate submit --local sim_validation/native_results/iris_check.sh -D /data1/core005/ccs/orgeraj/go_runs/native_check -p cmobic_cpu -t 12:00:00 --mem 16G -c 2 -J native_check` (hand to the user if blocked). Pass = 20/20 `ok`; any difference → fix in Tasks 4–6 with a new golden case reproducing it synthetically.
- [ ] **Step 4: Commit** scripts (`test: Iris real-data check for write-results`).

### Task 8: Image

**Files:**
- Modify: `Dockerfile` only if the binary is not already built from the repo root (check first).

- [ ] **Step 1:** `docker buildx build --platform linux/amd64,linux/arm64 -t docker.io/orgeraj/phylowgs-go:go-native-<shortsha> --push .` (ask the user before pushing — it publishes; if the Docker Hub token is expired, hand them `docker login`).
- [ ] **Step 2:** Verify: `docker run --rm --platform linux/amd64 docker.io/orgeraj/phylowgs-go:go-native-<shortsha> phylowgs-go write-results 2>&1 | head -3` prints the usage line; record the image digest.

### Task 9: Pipeline wiring

**Files (pipeline repo, branch `feat/phylowgs-go-native`):**
- Modify: `modules/local/phylowgs/writeresultsgo/main.nf` (container → `docker.io/orgeraj/phylowgs-go:go-native-<shortsha>`; script `phylowgs-go write-results ${args} --include-ssm-names ${prefix} ${trees} ${mutlist} ${prefix}.summ.json.gz ${prefix}.muts.json.gz ${prefix}.mutass.zip`; versions line `phylowgs_go: \$PHYLOWGS_GO_TAG`)
- Modify: `modules/local/phylowgs/multievolvego/main.nf` (same new image tag)
- Modify: `conf/modules.config`: `PHYLOWGS_WRITERESULTS` and `PHYLOWGS_WRITERESULTSGO` `ext.args = '--max-multiprimary 1.0 --include-multiprimary'`; `PHYLOWGS_WRITERESULTSGO` `ext.prefix = { params.phylowgs_engine == 'both' ? "${meta.id}-go" : "${meta.id}" }`; `PHYLOWGS_MULTIEVOLVEGO` `ext.args = { "-B ${params.phylo_burnin_samples} -s ${params.phylo_mcmc_samples} -j ${params.phylo_num_chains} -r ${Math.floorMod(meta.id.hashCode(), 1000000000) + 1}" }`
- Delete: `bin/write_results_go.py`
- Modify: `conf/test_go.config` header comment (image name)

- [ ] **Step 1:** `nextflow run . -profile test_go,docker -stub-run --outdir <scratch>/stub_go` → completes; `<scratch>/stub_go` has `<id>.summ.json.gz` (no `-go`). Same with `--phylowgs_engine both` → both `<id>.*` and `<id>-go.*`.
- [ ] **Step 2:** `nextflow config -profile test_go` shows the `-r <n>` arg resolved per sample (or print it from a stub task's `.command.sh`): two different sample ids give different seeds, neither 0.
- [ ] **Step 3: Commit** (`feat(phylowgs-go): native Go write-results, per-sample seed, RC filters`) and push the branch; no PR yet.

### Task 10: End-to-end check on Iris

**Files:**
- Create: `sim_validation/native_results/compare_input_json.py` (Go repo; python3; reports key paths with list indices collapsed to `[]` and JSON types that differ between two `*_input.json` files — never values; counts per path)

- [ ] **Step 1:** Reuse the launch script and samplesheet pattern of `/data1/core005/ccs/orgeraj/go_runs/test1` (read it via `hpcgate run "cat …"`, script only). Pick 3–5 of its samples spanning small → large mutation counts whose Python results exist under NeoQual `MASTER_RESULTS`. Launch with `--phylowgs_engine go`, `-r feat/phylowgs-go-native`, `--short_queue cmobic_short`, fresh launch dir, no storeDir, via `hpcgate submit`.
- [ ] **Step 2:** For each sample, `compare_input_json.py <go>/…_input.json <MASTER_RESULTS>/…_input.json` → pass = no type/key-path differences (missing optional paths that differ only because a tree had fewer clones are reported separately as counts, not failures).
- [ ] **Step 3:** Commit the script; write the result (counts only) to the vault project log.
