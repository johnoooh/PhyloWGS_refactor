# PhyloWGS-Go

Pure Go implementation of PhyloWGS - a Bayesian MCMC sampler for cancer phylogenetics.

## Overview

PhyloWGS infers the evolutionary history of tumors from bulk sequencing data. This Go implementation provides **~10x speedup** over the original Python/C++ version through:

- True parallel MCMC chains via goroutines (no GIL)
- Pure Go implementation (no CGo, no external dependencies)
- In-memory tree operations (no disk I/O between iterations)
- Single static binary deployment

## Installation

### CPU-only build (recommended for small datasets)

```bash
# Requires Go 1.18+
go build -o phylowgs-go .
```

### CUDA-accelerated build (for large datasets)

For datasets with hundreds or thousands of SSMs, GPU acceleration can provide significant speedups.

**Requirements:**
- CUDA toolkit 11.7+ (12.x recommended for newer GPUs)
- NVIDIA GPU with compute capability 8.6+ (Ampere/Ada/Blackwell)
- GCC/G++ for linking

```bash
# 1. Build the CUDA library
cd cuda && ./build.sh && cd ..

# 2. Build Go binary with CUDA support
go build -tags cuda -o phylowgs-go .

# 3. Set library path at runtime
export LD_LIBRARY_PATH=$PWD/cuda:/usr/local/cuda/lib64:$LD_LIBRARY_PATH
./phylowgs-go [options] ssm_data.txt cnv_data.txt
```

**Note:** For small datasets (<100 SSMs), the CPU-only build is typically faster due to GPU launch overhead.

## Usage

```bash
./phylowgs-go [options] ssm_data.txt cnv_data.txt
```

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `-s` | 2500 | Number of MCMC samples |
| `-B` | 1000 | Burn-in samples |
| `-j` | 6 | Number of parallel chains |
| `-O` | output | Output directory |
| `-i` | 5000 | MH iterations per MCMC iteration |
| `-r` | 0 | Random seed (0 = use time) |
| `-I` | 1.0 | Chain inclusion factor (fraction of best-LLH chains kept in `trees.zip`/`mutass.zip`) |
| `-D` | "" | Dataset name embedded in `mutass.zip` entries |
| `--no-gpu` | false | Disable GPU acceleration (CUDA build only) |

### Example

```bash
# Quick test run
./phylowgs-go -B 20 -s 50 -j 4 -O quick_test ssm_data.txt cnv_data.txt

# Production run (match original PhyloWGS defaults)
./phylowgs-go -B 1000 -s 2500 -j 6 -O results ssm_data.txt cnv_data.txt
```

## Input Format

### SSM Data (ssm_data.txt)

Tab-separated file with columns:
- `id`: SSM identifier (e.g., s0, s1, ...)
- `gene`: Gene name
- `a`: Variant read counts per sample (comma-separated)
- `d`: Total read depth per sample (comma-separated)
- `mu_r`: Reference allele probability (default: 0 if column absent — matches Python `util2.py`)
- `mu_v`: Variant allele probability (default: 0 if column absent — matches Python `util2.py`)

```tsv
id	gene	a	d	mu_r	mu_v
s0	BRCA1	66023,50883,62757	126755,100469,121941	0.999	0.5
s1	TP53	71532,50933,64719	135031,97485,120826	0.999	0.5
```

### CNV Data (cnv_data.txt)

Tab-separated file with columns:
- `cnv`: CNV identifier
- `a`: Variant reads (comma-separated per sample)
- `d`: Total depth (comma-separated per sample)
- `ssms`: Affected SSMs (`<ssm_id>,<maternal_cn>,<paternal_cn>;...`)
- `physical_cnvs`: Physical CNV details

```tsv
cnv	a	d	ssms	physical_cnvs
c0	66023,50883,62757	126755,100469,121941	s2,1,2;s4,0,1	chrom=1,start=1234,end=5678,...
```

## Outputs

After a run, `<output-dir>` contains:

- `best_tree.json` — single-best tree (highest-LLH sample across all
  included chains).
- `mutlist.json` — mutation reference table.
- `summary.json` — per-chain trace and global summary (see below).
- `chain_<N>_trees.ndjson` — per-chain newline-delimited JSON snapshots.
- `chain_<N>_samples.txt` — per-chain LLH trace.
- `trees.zip` — every included sample's tree snapshot, named
  `tree_<idx>_<llh>`. Matches Python `multievolve.py` `trees.zip`
  layout, with JSON payloads instead of pickle.
- `mutass.zip` — every included sample's mutation assignments, named
  `<idx>.json`. Includes the dataset name supplied via `-D`.

`summary.json` example:

```json
{
  "num_chains": 6,
  "total_time": "38.694s",
  "chain_times": ["38.6s", "38.5s", ...],
  "best_llh": -1748041.14,
  "trees_per_chain": [2500, 2500, ...]
}
```

`chain_<N>_samples.txt` example:

```tsv
Iteration	LLH	NumNodes
0	-1748042.31	3
1	-1748041.89	3
...
```

## Posterior-tree summaries

After a run, group MCMC trees by topology and emit per-group LaTeX
summaries:

```bash
phylowgs-go posterior-trees [-n 5] [-no-pdf] <output-dir>
```

Reads `<output-dir>/trees.zip`, groups samples by topology signature,
ranks groups by posterior probability, and writes
`<output-dir>/posterior_trees/tree_<rank>_<probability>.tex` for the
top groups (ranked by descending posterior probability). If `pdflatex`
is on PATH, the corresponding PDFs are produced alongside.

Flags:

- `-n N` — number of top groups to emit (default 5).
- `-no-pdf` — skip `pdflatex` invocation; emit `.tex` only.

## Performance

On AMD Ryzen 9800X3D (8 cores / 16 threads):

| Config | Samples | Chains | Time (Go) | Time (Python) | Speedup |
|--------|---------|--------|-----------|---------------|---------|
| Production | 2500 | 6 | 39s | 381s | 9.8x |
| Medium | 1000 | 6 | 17s | 170s | 10.2x |
| Fast | 100 | 4 | 1.7s | 17s | 10.0x |

Throughput: **60-66 samples/sec/chain** (vs 8-12 for Python)

## Algorithm

PhyloWGS uses a Tree-Structured Stick-Breaking Process (TSSB) as the prior over phylogenetic trees:

1. **Resample assignments**: Slice-sample each SSM's node assignment
2. **Metropolis-Hastings**: Sample cellular prevalence parameters
3. **Resample sticks**: Update tree branch weights
4. **Resample hyperparameters**: Update dp_alpha, dp_gamma, alpha_decay
5. **Cull tree**: Remove empty nodes

The tree structure encodes clonal evolution, with nodes representing clones and branches representing descent relationships.

## Limitations

- CNV handling is simplified (full genome counting not yet implemented)
- No pickle output (for Python interop, use JSON export)
- Single-machine only (no distributed MCMC)

## License

MIT

## References

- Deshwar et al. (2015) "PhyloWGS: Reconstructing subclonal composition and evolution from whole-genome sequencing of tumors"
- Jiao et al. (2014) "Inferring clonal evolution of tumors from single nucleotide somatic mutations"
