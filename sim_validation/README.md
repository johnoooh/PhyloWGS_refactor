# PhyloWGS Simulation Validation Pipeline

## Purpose

Validates the Go port of PhyloWGS against known ground truth by:
1. Generating simulated tumor phylogenies with known tree structure, mutations, and cellular prevalences
2. Running the Go port MCMC sampler on the simulated data
3. Scoring the inferred tree against the ground truth
4. Optionally comparing against the original Python PhyloWGS

## Quick Start (No SLURM)

### Prerequisites

- **Go 1.21+** (to build the Go binary)
- **Python 3.10+** with `numpy` and `matplotlib`
- *Optional*: Singularity + PhyloWGS container image (for Python comparison)

### Step 1: Build the Go binary

```bash
cd /path/to/PhyloWGS_refactor   # root of the go-port branch
go build -o phylowgs-go .
```

### Step 2: Generate fixtures

```bash
# Single fixture
python sim_validation/generate_fixtures.py \
    --outdir /tmp/fixtures/K3_test \
    -K 3 -S 1 -T 1000 -M 30 --seed 42

# Grid of fixtures (quick = 32 fixtures)
python sim_validation/generate_fixtures.py \
    --outdir /tmp/fixtures --grid quick
```

Each fixture directory contains:
- `ssm_data.txt` — PhyloWGS SSM input (tab-separated: id, gene, a, d, mu_r, mu_v)
- `cnv_data.txt` — PhyloWGS CNV input (empty header for SSM-only, populated for CNV fixtures)
- `truth.json` — Ground truth: parents vector, phi, eta, SSM assignments, clusters, CNV info

### Step 3: Run the Go port on each fixture

```bash
# Single fixture
./phylowgs-go --no-gpu -B 500 -s 1000 -j 4 \
    -O /tmp/results/K3_test \
    /tmp/fixtures/K3_test/ssm_data.txt \
    /tmp/fixtures/K3_test/cnv_data.txt

# Loop over all fixtures
for dir in /tmp/fixtures/K*; do
    name=$(basename "$dir")
    mkdir -p /tmp/results/"$name"
    ./phylowgs-go --no-gpu -B 500 -s 1000 -j 4 \
        -O /tmp/results/"$name" \
        "$dir/ssm_data.txt" \
        "$dir/cnv_data.txt"
done
```

**Go port CLI flags:**
| Flag | Default | Description |
|------|---------|-------------|
| `-B` | 1000 | Burn-in iterations |
| `-s` | 2500 | MCMC sample iterations |
| `-j` | 6 | Parallel chains |
| `-O` | output | Output directory |
| `--no-gpu` | false | Disable CUDA (use for CPU-only) |
| `-i` | 5000 | Metropolis-Hastings iterations per MCMC step |
| `-r` | 0 | Random seed (0 = use time) |

**Go port outputs per run:**
- `summary.json` — timing, best LLH, chain info
- `best_tree.json` — **final tree**: populations, cellular prevalences, SSM assignments, tree structure
- `chain_N_samples.txt` — per-chain MCMC trace (Iteration, LLH, NumNodes)

### Step 4: Score results

```bash
# Single fixture
python sim_validation/score_results.py \
    --fixture-dir /tmp/fixtures/K3_test \
    --result-dir /tmp/results/K3_test

# All fixtures
python sim_validation/score_results.py --all \
    --fixture-base /tmp/fixtures \
    --result-base /tmp/results \
    --outdir /tmp/analysis
```

### Step 5: Generate plots

```bash
python sim_validation/plot_results.py \
    --analysis-dir /tmp/analysis \
    --result-base /tmp/results \
    --fixture-base /tmp/fixtures
```

Plots saved to `/tmp/analysis/plots/`.

## Simulation Model

### Generative Process (SSM-only)

1. **Random tree** — K clones + root (node 0). Parents assigned stochastically (mu=0.75 linear extension probability)
2. **Subclone proportions** — eta ~ Dirichlet(alpha=0.1) per sample, shape (K+1, S)
3. **Ancestral matrix** — Z[i,j] = 1 if node i is ancestor of node j
4. **Cellular prevalences** — phi = Z @ eta, shape (K+1, S)
5. **SSM assignments** — M mutations randomly assigned to clones 1..K (≥1 per clone)
6. **Read counts** — ref_reads ~ Binomial(T, 1 - 0.5 * phi[clone, sample])

### CNV Extension (when `--N-cnv > 0`)

Adds copy number aberrations to the tree:

1. **CNV events** assigned to clones with sampled copy number states: (1,2), (2,1), (0,1), (1,0), (0,2), (2,0), (1,3)
2. **SSM-CNV overlap** — ~30% of SSMs affected per CNV (max 1 CNV per SSM)
3. **Timing classification** per SSM-CNV pair:
   - `ssm_before_cnv` — variant gets amplified/deleted by the CNV
   - `cnv_before_ssm` — variant on post-CNV background
   - `same_node` — ambiguous timing (4 genotype scenarios)
   - `independent` — neither is ancestor of the other
4. **CNV-aware read counts** — mirrors Go port's `computeNGenomes` (main.go:756-863) exactly, computing (nr, nv) pairs per tree node and averaging over valid genotype scenarios

### Parameter Grids

| Grid | K | S | T | M/K | N_cnv | Reps | Total |
|------|---|---|---|-----|-------|------|-------|
| `quick` | 3,5 | 1,3 | 200,1000 | 10 | 0,2 | 2 | 32 |
| `default` | 3,5,10 | 1,3 | 200,1000 | 10,50 | 0,2,5 | 3 | 216 |
| `thorough` | 3,5,10,20 | 1,3,10 | 50,200,1000 | 10,50 | 0,2,5,10 | 4 | 1536 |

### Fixture File Formats

**ssm_data.txt** (tab-separated):
```
id    gene    a              d              mu_r   mu_v
s0    SIM_0   900,850,800    1000,1000,1000 0.999  0.5
```
- `a` = reference reads per sample (comma-separated)
- `d` = total reads per sample
- `mu_r` = P(ref read | ref allele) = 0.999
- `mu_v` = P(ref read | var allele) = 0.5

**cnv_data.txt** (tab-separated):
```
cnv   a            d              ssms                physical_cnvs
c0    8500,9000    10000,10000    s2,1,2;s4,0,1      chrom=1,start=...
```
- `ssms` = semicolon-separated `ssm_id,maternal_cn,paternal_cn` triplets

**truth.json**:
```json
{
  "parents": [0, 1, 0],
  "phi": [[1.0, 1.0], [0.5, 0.4], [0.3, 0.2], [0.2, 0.1]],
  "eta": [[...], ...],
  "assignments": [1, 2, 1, 3, ...],
  "clusters": [["s0", "s2"], ["s1"], ["s3"]],
  "cnvs": [...],
  "params": {"K": 3, "S": 2, "T": 1000, "M": 30, ...}
}
```

## Scoring Metrics

| Metric | Source | Description |
|--------|--------|-------------|
| `K_error` | `best_tree.json` | \|inferred populations - true K\| |
| `cocluster_auprc` | `best_tree.json` | Co-clustering AUPRC: do SSMs from the same clone cluster together? |
| `best_llh` | `summary.json` | Best log-likelihood across all chains |
| `median_llh` | chain files | Median LLH of post-burnin samples |
| `llh_chain_std` | chain files | Std dev of final LLH across chains (convergence) |
| `total_time_s` | `timing.json` or `summary.json` | Wall-clock runtime |

## Go Port Output: `best_tree.json`

The Go port extracts the final tree from the best chain after MCMC, runs `removeEmptyNodes()` (matching Python's `util2.py:128-155`), and writes:

```json
{
  "chain_id": 0,
  "llh": -1234.5,
  "num_populations": 3,
  "structure": {"0": [1, 2], "1": [3]},
  "populations": {
    "0": {"cellular_prevalence": [1.0], "num_ssms": 0, "num_cnvs": 0},
    "1": {"cellular_prevalence": [0.5], "num_ssms": 15, "num_cnvs": 0},
    "2": {"cellular_prevalence": [0.3], "num_ssms": 10, "num_cnvs": 0},
    "3": {"cellular_prevalence": [0.1], "num_ssms": 5, "num_cnvs": 0}
  },
  "mut_assignments": {
    "1": {"ssms": ["s0", "s1", ...], "cnvs": []},
    "2": {"ssms": ["s5", "s6", ...], "cnvs": []}
  }
}
```

Population 0 is always root (normal cells). Children sorted by decreasing mean phi.

## Comparing Against Original Python

If you have a singularity container with the original PhyloWGS:

```bash
# Run original Python on a fixture
singularity exec phylowgs.sif python2 /usr/bin/phylowgs/evolve.py \
    -B 500 -s 1000 \
    ssm_data.txt cnv_data.txt

# Extract results
singularity exec phylowgs.sif python2 /usr/bin/phylowgs/write_results.py \
    fixture_name trees.zip \
    tree_summaries.json.gz mutlist.json.gz mutass.zip \
    --include-ssm-names
```

Without singularity, you need the original morrislab/phylowgs repo with Python 2 + numpy + scipy + GSL (for mh.o compilation).

The scoring script handles both output formats automatically.

## SLURM Usage

For HPC with SLURM, use the wrapper scripts:

```bash
# One-time setup
bash sim_validation/setup_validation.sh --workdir $SCRATCH/sim_val --grid quick

# Submit all jobs (Go + Python)
bash sim_validation/submit_sim_validation.sh \
    --workdir $SCRATCH/sim_val \
    --sif /path/to/phylowgs.sif \
    --mem 8000
```
