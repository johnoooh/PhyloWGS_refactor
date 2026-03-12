# PhyloWGS SLURM Benchmark Pipeline

End-to-end pipeline: MAF + FACETS → PhyloWGS runs across all implementations → comparison report.

## Quick Start

```bash
# 1. One-time setup (on login node)
bash setup.sh --workdir /path/to/workdir

# 2. Prepare your CSV
cat samples.csv
# sample_id,maf,facets
# SAMPLE001,/path/to/SAMPLE001.maf,/path/to/SAMPLE001_hisens.cncf.txt
# SAMPLE002,/path/to/SAMPLE002.maf,/path/to/SAMPLE002_hisens.cncf.txt

# 3. Submit all jobs
bash submit_benchmark.sh samples.csv --workdir /path/to/workdir

# 4. Monitor
squeue -u $USER

# 5. Results land in:
# /path/to/workdir/analysis/report.html
# /path/to/workdir/analysis/summary.tsv
# /path/to/workdir/analysis/concordance.tsv
```

## What gets submitted

For each sample:
1. **Convert job** (`cmobic_cpu`) — MAF + FACETS → `ssm_data.txt` + `cnv_data.txt`
2. **Run jobs** (one per implementation, after convert):
   - `optimized-python` → `cmobic_cpu`
   - `go-cpu` → `cmobic_cpu`
   - `go-cpu-opt` → `cmobic_cpu`
   - `go-gpu` → `gpu` partition, `--gres=gpu:1`
3. **Analysis job** (`cmobic_cpu`) — depends on ALL run jobs completing

## Implementations

| Name | Branch | Partition | Notes |
|---|---|---|---|
| `optimized-python` | `main` | cmobic_cpu | Python 3 refactor |
| `go-cpu` | `go/main` | cmobic_cpu | Pure Go, no opts |
| `go-cpu-opt` | `go/feature/parallel-traversal` | cmobic_cpu | MH caching, 2.15× |
| `go-gpu` | `go/feature/cuda-likelihood` | gpu | CUDA kernel |

## Options

```
--burnin N      MCMC burn-in (default: 1000)
--samples N     MCMC samples (default: 2500)
--chains N      Parallel chains (default: 4)
--time HH:MM:SS Wall time limit (default: 4:00:00)
--mem MB        Memory per job (default: 8000)
--dry-run       Print jobs without submitting
```

## Input CSV Format

```
sample_id,maf,facets
SAMPLE001,/abs/path/to.maf,/abs/path/to_hisens.cncf.txt
```

- `maf` — standard MAF file (with t_ref_count/t_alt_count columns)
- `facets` — FACETS `*_hisens.cncf.txt` output

## Analysis Outputs

- **`report.html`** — speedup table + wall times + LLH concordance per sample
- **`summary.tsv`** — `sample_id, impl, status, wall_seconds, best_llh, median_llh`
- **`concordance.tsv`** — best LLH comparison across implementations (threshold: <5% rel diff)

## File Structure

```
workdir/
├── env.sh                    # PATH setup (sourced by all SLURM jobs)
├── implementations.tsv       # Implementation manifest
├── convert_inputs.py         # MAF + FACETS converter
├── impls/
│   ├── optimized-python/     # Cloned + uv venv
│   ├── go-cpu/               # Cloned + compiled binary
│   ├── go-cpu-opt/           # Cloned + compiled binary
│   └── go-gpu/               # Cloned + CUDA binary
├── inputs/
│   └── SAMPLE001/
│       ├── ssm_data.txt
│       └── cnv_data.txt
├── results/
│   └── SAMPLE001/
│       ├── optimized-python/
│       ├── go-cpu/
│       ├── go-cpu-opt/
│       └── go-gpu/
├── logs/                     # SLURM stdout/stderr
└── analysis/
    ├── report.html
    ├── summary.tsv
    └── concordance.tsv
```
