#!/usr/bin/env bash
# setup.sh — One-time setup on HPC login node
# Clones all PhyloWGS branches, builds uv envs and Go binaries
#
# Usage: bash setup.sh [--workdir /path/to/workdir]
#
# Requires: git, uv, Go (or will download Go), CUDA toolkit (for GPU build)

set -euo pipefail

REPO="https://github.com/johnoooh/PhyloWGS_refactor.git"
WORKDIR="${PHYLOWGS_WORKDIR:-$(pwd)/phylowgs_benchmark}"
GO_VERSION="1.22.0"
GO_INSTALL_DIR="$WORKDIR/go_install"

while [[ $# -gt 0 ]]; do
    case $1 in
        --workdir) WORKDIR="$2"; shift 2 ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

echo "=== PhyloWGS Benchmark Setup ==="
echo "Workdir: $WORKDIR"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

# ── Go installation ──────────────────────────────────────────────────────────
if ! command -v go &>/dev/null; then
    echo "[go] Not found — installing Go $GO_VERSION to $GO_INSTALL_DIR"
    mkdir -p "$GO_INSTALL_DIR"
    curl -sSL "https://go.dev/dl/go${GO_VERSION}.linux-amd64.tar.gz" \
        | tar -C "$GO_INSTALL_DIR" -xz
    export PATH="$GO_INSTALL_DIR/go/bin:$PATH"
    echo "export PATH=\"$GO_INSTALL_DIR/go/bin:\$PATH\"" >> "$WORKDIR/env.sh"
else
    GO_BIN=$(command -v go)
    echo "[go] Found at $GO_BIN ($(go version))"
    echo "export PATH=\"$(dirname $GO_BIN):\$PATH\"" >> "$WORKDIR/env.sh"
fi

# ── uv installation ──────────────────────────────────────────────────────────
if ! command -v uv &>/dev/null; then
    echo "[uv] Not found — installing"
    curl -LsSf https://astral.sh/uv/install.sh | sh
    export PATH="$HOME/.cargo/bin:$PATH"
    echo "export PATH=\"$HOME/.cargo/bin:\$PATH\"" >> "$WORKDIR/env.sh"
fi

# ── Define implementations ───────────────────────────────────────────────────
declare -A BRANCHES=(
    ["optimized-python"]="main"
    ["go-cpu"]="go/main"
    ["go-cpu-opt"]="go/feature/parallel-traversal"
    ["go-gpu"]="go/feature/cuda-likelihood"
)

# ── Clone each branch ────────────────────────────────────────────────────────
for impl in "${!BRANCHES[@]}"; do
    branch="${BRANCHES[$impl]}"
    dest="$WORKDIR/impls/$impl"

    if [[ -d "$dest" ]]; then
        echo "[$impl] Already cloned — pulling"
        git -C "$dest" pull --rebase --autostash 2>&1 | tail -2
    else
        echo "[$impl] Cloning branch $branch"
        git clone --depth=1 --branch "$branch" "$REPO" "$dest" 2>&1 | tail -3
    fi
done

# ── Python env: optimized-python ────────────────────────────────────────────
echo "[optimized-python] Building uv environment"
cd "$WORKDIR/impls/optimized-python"
uv venv --python 3.10 .venv
uv pip install --quiet numpy scipy 2>&1 | tail -2

# ── Go builds ───────────────────────────────────────────────────────────────
source "$WORKDIR/env.sh" 2>/dev/null || true

for impl in go-cpu go-cpu-opt; do
    echo "[$impl] Building CPU binary"
    cd "$WORKDIR/impls/$impl"
    go build -o "phylowgs-cpu" . 2>&1
    echo "[$impl] Built: $(ls -lh phylowgs-cpu | awk '{print $5, $9}')"
done

# ── CUDA build: go-gpu ───────────────────────────────────────────────────────
echo "[go-gpu] Building CUDA binary"
cd "$WORKDIR/impls/go-gpu"
if command -v nvcc &>/dev/null; then
    cd cuda && bash build.sh 2>&1 | tail -5 && cd ..
    CGO_ENABLED=1 LD_LIBRARY_PATH="$WORKDIR/impls/go-gpu/cuda" \
        go build -tags cuda -o phylowgs-gpu . 2>&1
    echo "[go-gpu] CUDA build OK"
else
    echo "[go-gpu] WARNING: nvcc not found — falling back to CPU build"
    go build -o phylowgs-gpu . 2>&1
fi

# ── Input conversion scripts ─────────────────────────────────────────────────
# Copy conversion helpers into workdir
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cp "$SCRIPT_DIR/convert_inputs.py" "$WORKDIR/convert_inputs.py" 2>/dev/null || \
    echo "[warn] convert_inputs.py not found — copy manually"

# ── Write implementation manifest ────────────────────────────────────────────
cat > "$WORKDIR/implementations.tsv" << EOF
name	type	binary	branch	partition	extra_args
optimized-python	python	python	main	cmobic_cpu	
go-cpu	go	phylowgs-cpu	go/main	cmobic_cpu	--no-gpu
go-cpu-opt	go	phylowgs-cpu	go/feature/parallel-traversal	cmobic_cpu	--no-gpu
go-gpu	go	phylowgs-gpu	go/feature/cuda-likelihood	gpu	
EOF

echo ""
echo "=== Setup complete ==="
echo "Workdir: $WORKDIR"
echo "Next: bash submit_benchmark.sh <samples.csv> --workdir $WORKDIR"
