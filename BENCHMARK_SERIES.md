# PhyloWGS Optimized — Benchmark Series

**Date:** 2026-03-10  
**Host:** DESKTOP-BLN7CUN (WSL2, Linux 5.10)  
**Dataset:** 11 SSMs, 5 samples (`ssm_data.txt` / `cnv_data.txt`)  
**Constraint:** max 6 chains simultaneously  

## Results

| Run | Label | Burnin | MCMC samples | Chains | Total iters | Wall time (s) | Samp/s/chain |
|-----|-------|--------|-------------|--------|-------------|---------------|--------------|
| 1 ★ | production-6c | 1000 | 2500 | 6 | 21,000 | 381 | 9.16 |
| 2 | medium | 500 | 1000 | 6 | 9,000 | 170 | 8.82 |
| 3 | light | 200 | 500 | 4 | 2,800 | 59 | 11.86 |
| 4 | lighter | 100 | 250 | 4 | 1,400 | 46 | 7.61 |
| 5 | fast | 50 | 100 | 4 | 600 | 17 | 8.82 |
| 6 | minimal | 20 | 50 | 4 | 280 | 9 | 7.78 |

★ = production configuration

## Key Observations

- **Roughly linear** scaling: production (21k iters) takes ~42× longer than minimal (280 iters), ratio = 381/9 ≈ 42×  
- **Throughput** is consistent at ~8–12 samples/sec/chain across all configs  
- **6-chain overhead** for medium vs light (both 4 chains but different chain counts for medium): chain management adds modest overhead  
- Startup cost is ~6–8s regardless of iteration count (visible from minimal run = 9s for only 280 iters)

## Scaling Estimate for Larger Datasets

On 11 SSMs, production config = 6.35 min. For N SSMs, expect roughly O(N²) scaling due to tree topology space — production on 100 SSMs may take hours.

## Files

- `benchmark_runs/results.json` — machine-readable timing data
- `runtime_plot.png` — runtime vs iteration count plot (log-log)
- `benchmark_runs/run_*/` — per-run output directories

## Notes

- Original PhyloWGS comparison pending (requires Docker/Singularity — scheduled for SLURM tomorrow)
- All runs used `ssm_data.txt` and `cnv_data.txt` from the optimized repo root
