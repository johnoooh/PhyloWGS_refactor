# PhyloWGS: Original (Python 2) vs Optimized (Python 3) Benchmark

Generated: 2026-03-10

## Setup

- **Original:** Python 2.7.18, `/home/john/.openclaw/workspace/phylowgs/`  
  Deps: numpy 1.16.6, scipy 1.2.3, ete2 2.3.10, pyvcf 0.6.8  
  Compiled: `g++ -o mh.o mh.cpp util.cpp $(gsl-config --cflags --libs)`
- **Optimized:** Python 3.x, `/home/john/.openclaw/workspace/phylowgs-optimized/`
- **Test data:** `ssm_data.txt` + `cnv_data.txt` (same dataset for both)
- **Host:** DESKTOP-BLN7CUN (WSL2, Linux 5.10)

## Status

✅ **Original ran successfully** — all 3 benchmark configs completed with exit code 0.

## Results

| Config | B (burnin) | s (samples) | Chains | Orig Py2 (s) | Opt Py3 (s) | Speedup |
|--------|-----------|-------------|--------|--------------|-------------|---------|
| minimal | 20 | 50 | 4 | 12.3 | 9 | **1.37×** |
| fast | 50 | 100 | 4 | 19.1 | 17 | **1.12×** |
| production-6c | 1000 | 2500 | 6 | 519.6 | 381 | **1.36×** |

## Notes

- Speedup = Original time / Optimized time (>1.0 = optimized is faster)
- All runs used the same flags and dataset
- Chain 2 in production-6c ran notably slower on the original (~8m40s total vs 6m21s for optimized), likely due to Python 2 GIL + subprocess overhead and slower numpy 1.16 ops
- The optimized version (Py3) is consistently **~1.1–1.4× faster** across all configs
- At production scale (B=1000, s=2500), the optimized version saves **~139 seconds** (~2m19s) per run

## Conclusion

The Python 3 port delivers a real and consistent speedup over the original Python 2 codebase — approximately **1.3× faster** at production scale. The gains are likely from:
1. Python 3 interpreter improvements
2. Newer numpy/scipy (Py3 can use modern optimized builds)
3. Any algorithmic changes in the port
