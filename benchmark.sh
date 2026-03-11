#!/bin/bash
# Benchmark script for PhyloWGS Go vs Python comparison

set -e

# Paths
WORKSPACE="/home/john/.openclaw/workspace"
GO_BIN="$WORKSPACE/phylowgs-go/phylowgs-go"
PY_DIR="$WORKSPACE/phylowgs-optimized"
SSM_FILE="$PY_DIR/ssm_data.txt"
CNV_FILE="$PY_DIR/cnv_data.txt"
RESULTS_DIR="$WORKSPACE/phylowgs-go/benchmark_results"

# Add conda to path
export PATH=~/miniconda3/bin:$PATH

mkdir -p "$RESULTS_DIR"

# Benchmark configurations
# Format: name burnin samples chains
configs=(
    "production-6c 1000 2500 6"
    "medium 500 1000 6"
    "light 200 500 4"
    "lighter 100 250 4"
    "fast 50 100 4"
    "minimal 20 50 4"
)

echo "==============================================="
echo "PhyloWGS Go vs Python Benchmark"
echo "==============================================="
echo "Machine: $(hostname)"
echo "CPUs: $(nproc)"
echo "Go version: $(go version)"
echo "Python: $(python3 --version)"
echo "Dataset: 11 SSMs, 5 samples"
echo "==============================================="
echo ""

# Build Go binary with optimizations
echo "Building Go binary with optimizations..."
cd "$WORKSPACE/phylowgs-go"
go build -ldflags="-s -w" -o phylowgs-go .
echo ""

# Results arrays
declare -a go_times
declare -a py_times

for config in "${configs[@]}"; do
    read -r name burnin samples chains <<< "$config"
    
    echo "=== Configuration: $name ==="
    echo "  Burnin: $burnin, Samples: $samples, Chains: $chains"
    
    # Run Go implementation
    echo "  Running Go implementation..."
    GO_OUT="$RESULTS_DIR/go_${name}"
    rm -rf "$GO_OUT"
    GO_START=$(date +%s.%N)
    "$GO_BIN" -B "$burnin" -s "$samples" -j "$chains" -O "$GO_OUT" "$SSM_FILE" "$CNV_FILE" 2>&1 | tail -10
    GO_END=$(date +%s.%N)
    GO_TIME=$(echo "$GO_END - $GO_START" | bc)
    echo "  Go time: ${GO_TIME}s"
    go_times+=("$GO_TIME")
    
    # Run Python implementation
    echo "  Running Python implementation..."
    PY_OUT="$RESULTS_DIR/py_${name}"
    rm -rf "$PY_OUT"
    cd "$PY_DIR"
    PY_START=$(date +%s.%N)
    python3 multievolve.py --ssms "$SSM_FILE" --cnvs "$CNV_FILE" \
        -B "$burnin" -s "$samples" -n "$chains" -O "$PY_OUT" 2>&1 | tail -20
    PY_END=$(date +%s.%N)
    PY_TIME=$(echo "$PY_END - $PY_START" | bc)
    echo "  Python time: ${PY_TIME}s"
    py_times+=("$PY_TIME")
    
    # Calculate speedup
    SPEEDUP=$(echo "scale=2; $PY_TIME / $GO_TIME" | bc)
    echo "  Speedup: ${SPEEDUP}x"
    echo ""
done

# Print summary table
echo "==============================================="
echo "                    SUMMARY"
echo "==============================================="
echo ""
printf "%-15s %12s %12s %12s\n" "Config" "Go (s)" "Python (s)" "Speedup"
printf "%-15s %12s %12s %12s\n" "------" "-------" "----------" "-------"

i=0
for config in "${configs[@]}"; do
    read -r name burnin samples chains <<< "$config"
    go_t="${go_times[$i]}"
    py_t="${py_times[$i]}"
    speedup=$(echo "scale=2; $py_t / $go_t" | bc)
    printf "%-15s %12.2f %12.2f %12.2fx\n" "$name" "$go_t" "$py_t" "$speedup"
    i=$((i+1))
done

echo ""
echo "Benchmark complete. Results in $RESULTS_DIR"
