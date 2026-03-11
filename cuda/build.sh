#!/bin/bash
# build.sh - Build CUDA likelihood library for PhyloWGS-Go
# 
# Usage: ./build.sh
# 
# Outputs: libcuda_bridge.so

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Find CUDA
if [ -z "$CUDA_HOME" ]; then
    if [ -d "/usr/local/cuda" ]; then
        export CUDA_HOME="/usr/local/cuda"
    elif [ -d "/usr/local/cuda-11.7" ]; then
        export CUDA_HOME="/usr/local/cuda-11.7"
    else
        echo "Error: CUDA not found. Set CUDA_HOME environment variable."
        exit 1
    fi
fi

NVCC="$CUDA_HOME/bin/nvcc"
if [ ! -x "$NVCC" ]; then
    echo "Error: nvcc not found at $NVCC"
    exit 1
fi

echo "Using CUDA from: $CUDA_HOME"
$NVCC --version | head -4

# Detect GPU compute capability
GPU_CC=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | head -1 | tr -d '.')
if [ -z "$GPU_CC" ]; then
    GPU_CC="89"  # Default to Ada Lovelace
    echo "Warning: Could not detect GPU compute capability, using $GPU_CC"
else
    echo "Detected GPU compute capability: ${GPU_CC:0:1}.${GPU_CC:1:1}"
fi

# Determine best architecture to compile for
# CUDA 11.7 supports up to sm_89 (Ada Lovelace)
# For newer GPUs (Blackwell sm_120), use PTX forward compatibility
NVCC_VERSION=$($NVCC --version | grep -oP 'release \K[0-9]+\.[0-9]+')
NVCC_MAJOR=$(echo $NVCC_VERSION | cut -d. -f1)
NVCC_MINOR=$(echo $NVCC_VERSION | cut -d. -f2)

echo "NVCC version: $NVCC_VERSION"

# Build architecture flags
# We compile to the highest arch we support + PTX for forward compatibility
if [ "$NVCC_MAJOR" -ge 12 ]; then
    # CUDA 12+ supports newer architectures
    if [ "$GPU_CC" -ge 120 ]; then
        # Blackwell
        ARCH_FLAGS="-arch=compute_120 -code=sm_120,compute_120"
    elif [ "$GPU_CC" -ge 100 ]; then
        # Hopper
        ARCH_FLAGS="-arch=compute_100 -code=sm_100,compute_100"
    elif [ "$GPU_CC" -ge 90 ]; then
        # Hopper/Grace
        ARCH_FLAGS="-arch=compute_90 -code=sm_90,compute_90"
    else
        ARCH_FLAGS="-arch=compute_86 -code=sm_86,compute_86"
    fi
else
    # CUDA 11.x - max is compute_87 (Jetson Orin) but use 86 for broader compat
    # Generate PTX for compute_86 - driver will JIT compile for actual GPU
    ARCH_FLAGS="-arch=compute_86 -code=compute_86"
    echo "Note: Using PTX forward compatibility (compute_86) for GPU cc=$GPU_CC with CUDA $NVCC_VERSION"
    echo "      The driver will JIT compile PTX to native code at runtime."
fi

echo "Using architecture flags: $ARCH_FLAGS"

# Compile CUDA kernel
echo "Compiling CUDA kernel..."
$NVCC -O3 \
    $ARCH_FLAGS \
    -Xcompiler -fPIC \
    -Xptxas -v \
    -c likelihood_kernel.cu -o likelihood_kernel.o \
    2>&1 | grep -E "(registers|spill|ptxas)" || true

# Link shared library
echo "Linking shared library..."
$NVCC -shared \
    -o libcuda_bridge.so \
    likelihood_kernel.o \
    -lcudart

# Check result
if [ -f libcuda_bridge.so ]; then
    SIZE=$(ls -lh libcuda_bridge.so | awk '{print $5}')
    echo "Success: libcuda_bridge.so ($SIZE)"
    
    # Verify symbols
    echo "Exported symbols:"
    nm -D libcuda_bridge.so | grep -E "cuda_(init|compute|cleanup|available|upload)" | head -10
else
    echo "Error: Failed to create libcuda_bridge.so"
    exit 1
fi

# Cleanup intermediate files
rm -f likelihood_kernel.o

echo ""
echo "Build complete. To use:"
echo "  export LD_LIBRARY_PATH=$SCRIPT_DIR:\$LD_LIBRARY_PATH"
echo "  cd $SCRIPT_DIR/.. && go build"
