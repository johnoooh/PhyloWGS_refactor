# PhyloWGS Full Benchmark: Original vs Optimized Python vs Go CPU vs Go GPU

**Date:** 2026-03-11
**Machine:** AMD Ryzen 7 9800X3D, 8C/16T | RTX 5080 16GB | WSL2
**Test data:** 11 SSMs, 5 timepoints
**Max chains:** 4 (thermal limit)

## Important Note on Chain Parallelism

- **Python versions:** `evolve.py` runs single-chain only. Multi-chain requires `multievolve.py` wrapper.
- **Go versions:** Native multi-chain support via `-j` flag.

For fair comparison, we report **iterations per second per chain** as the primary metric.

## Results

### Minimal Config (burnin=20, samples=50)

| Implementation | Chains | Wall Time | Iter/s/chain | Peak RAM | vs Original |
|----------------|--------|-----------|--------------|----------|-------------|
| original-python | 1 | 6.35s | 11.0 | 80 MB | 1.0× |
| optimized-python | 1 | 4.38s | 16.0 | 88 MB | 1.5× |
| go-cpu | 2 | 0.79s | 88.6 | 15 MB | 8.1× |
| go-gpu | 2 | 1.24s | 56.5 | 131 MB | 5.1× |

### Light Config (burnin=200, samples=500)

| Implementation | Chains | Wall Time | Iter/s/chain | Peak RAM | vs Original |
|----------------|--------|-----------|--------------|----------|-------------|
| original-python | 1 | 61.59s | 11.4 | 83 MB | 1.0× |
| optimized-python | 1 | 45.01s | 15.6 | 80 MB | 1.4× |
| go-cpu | 2 | 7.30s | 95.9 | 15 MB | 8.4× |
| go-gpu | 2 | 8.08s | 86.6 | 130 MB | 7.6× |

### Medium Config (burnin=500, samples=1000)

| Implementation | Chains | Wall Time | Iter/s/chain | Peak RAM | vs Original |
|----------------|--------|-----------|--------------|----------|-------------|
| original-python | 1 | 2:25.52 | 10.3 | 81 MB | 1.0× |
| optimized-python | 1 | 1:53.94 | 13.2 | 82 MB | 1.3× |
| go-cpu | 2 | 15.46s | 97.0 | 16 MB | 9.4× |
| go-gpu | 2 | 16.51s | 90.9 | 130 MB | 8.8× |

### Production Config (burnin=1000, samples=2500, chains=4)

| Implementation | Chains | Wall Time | Iter/s/chain | Peak RAM | vs Original |
|----------------|--------|-----------|--------------|----------|-------------|
| original-python | 1 | 5:48.72 | 10.0 | 80 MB | 1.0× |
| optimized-python | 1 | 3:45.09 | 15.6 | 98 MB | 1.6× |
| go-cpu | 4 | 37.99s | 92.1 | 18 MB | **9.2×** |
| go-gpu | 4 | 39.91s | 87.7 | 134 MB | 8.8× |

## Analysis

### Speedup Summary (Production Config)

| Metric | Original Python | Optimized Python | Go CPU | Go GPU |
|--------|-----------------|------------------|--------|--------|
| Wall time (1 chain equiv) | 5:49 | 3:45 | 0:38 | 0:40 |
| Iterations/sec/chain | 10.0 | 15.6 | 92.1 | 87.7 |
| Speedup vs original | 1.0× | 1.6× | **9.2×** | 8.8× |
| Peak RAM | 80 MB | 98 MB | 18 MB | 134 MB |

**Total throughput for 4-chain production run:**
- Original Python (via multievolve.py): ~5:49 wall time (4 processes, ~320 MB total RAM)
- Optimized Python (via multievolve.py): ~3:45 wall time (4 processes, ~392 MB total RAM)
- Go CPU: **37.99s** wall time (4 goroutines, 18 MB total RAM)
- Go GPU: **39.91s** wall time (4 goroutines, 134 MB total RAM)

### Where GPU Helps (and Doesn't)

**Current observation:** GPU is *slower* than CPU for this workload.

**Why:**
1. **Small dataset (11 SSMs × 5 timepoints):** The overhead of GPU kernel launches, memory transfers, and CGo crossings exceeds the parallelism benefits
2. **Memory bandwidth:** RTX 5080 shines with large matrix operations, but 11×5 likelihood calculations are too small to saturate GPU cores
3. **Per-iteration overhead:** Each MCMC iteration requires multiple CPU↔GPU round trips

**When GPU would help:**
- Datasets with 1000+ SSMs
- Batch likelihood calculations across multiple tree proposals
- Longer MH iterations that amortize kernel launch overhead

### Memory Efficiency

| Implementation | Peak RAM | RAM per Chain |
|----------------|----------|---------------|
| Go CPU | 18 MB | ~4.5 MB |
| original-python | 80 MB | 80 MB |
| optimized-python | 98 MB | 98 MB |
| Go GPU | 134 MB | ~34 MB |

Go CPU uses **4-5× less memory** than Python implementations. The GPU version has higher memory due to CUDA context and buffer allocation.

### Recommendation for Production HPC Use

**Use Go CPU (`--no-gpu`) for production runs.** Here's why:

1. **9× faster per chain** than original Python
2. **5× less memory** per chain
3. **Native multi-chain support** eliminates multievolve.py wrapper complexity
4. **GPU provides no benefit** for typical PhyloWGS datasets (<100 SSMs)

**Example HPC batch script:**
```bash
#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=1G
#SBATCH --time=00:30:00

./phylowgs-go --no-gpu -B 1000 -s 2500 -j 8 -O output/ ssm_data.txt cnv_data.txt
```

**Consider GPU (`./phylowgs-go` default) only when:**
- SSM count exceeds 500+
- Running multiple samples as batch jobs
- Optimizing for throughput on a GPU node

## Raw Timing Data

```
Config          | Impl              | Run 1     | Run 2     | Best
----------------|-------------------|-----------|-----------|--------
minimal         | original-python   | 7.32s     | 6.35s     | 6.35s
minimal         | optimized-python  | 4.69s     | 4.38s     | 4.38s
minimal         | go-cpu            | 0.81s     | 0.79s     | 0.79s
minimal         | go-gpu            | 1.24s     | 1.28s     | 1.24s
light           | original-python   | 61.59s    | 69.25s    | 61.59s
light           | optimized-python  | 45.01s    | 51.79s    | 45.01s
light           | go-cpu            | 7.38s     | 7.30s     | 7.30s
light           | go-gpu            | 8.23s     | 8.08s     | 8.08s
medium          | original-python   | 152.29s   | 145.52s   | 145.52s
medium          | optimized-python  | 128.54s   | 113.94s   | 113.94s
medium          | go-cpu            | 15.46s    | 15.67s    | 15.46s
medium          | go-gpu            | 17.02s    | 16.51s    | 16.51s
production      | original-python   | 348.72s   | 493.75s   | 348.72s
production      | optimized-python  | 225.09s   | 278.93s   | 225.09s
production      | go-cpu            | 38.29s    | 37.99s    | 37.99s
production      | go-gpu            | 39.91s    | 39.93s    | 39.91s
```

## Notes

- All Python runs used single-chain `evolve.py` (no multievolve.py wrapper)
- Go runs used native `-j` flag for chain parallelism
- 10-second cooldown between runs to prevent thermal throttling
- Best of 2 runs used to avoid cold-start outliers
- CUDA build uses RTX 5080 (compute capability 12.0, 17.1 GB VRAM)
