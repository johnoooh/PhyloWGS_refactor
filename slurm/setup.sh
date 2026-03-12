#!/usr/bin/env bash
# setup.sh — One-time setup on HPC login node
# Clones all PhyloWGS implementations, builds uv envs and Go binaries
#
# Usage: bash setup.sh [--workdir /path/to/workdir]

set -euo pipefail

REPO="https://github.com/johnoooh/PhyloWGS_refactor.git"
ORIGINAL_REPO="https://github.com/morrislab/phylowgs.git"
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
> env.sh  # reset env file

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
else
    echo "[uv] Found at $(command -v uv)"
    echo "export PATH=\"$(dirname $(command -v uv)):\$PATH\"" >> "$WORKDIR/env.sh"
fi

mkdir -p "$WORKDIR/impls"

# ── Clone: original PhyloWGS (morrislab) ────────────────────────────────────
echo "[original-python] Cloning from morrislab/phylowgs"
if [[ -d "$WORKDIR/impls/original-python" ]]; then
    git -C "$WORKDIR/impls/original-python" pull --rebase --autostash 2>&1 | tail -2
else
    git clone --depth=1 "$ORIGINAL_REPO" "$WORKDIR/impls/original-python" 2>&1 | tail -3
fi
echo "[original-python] Building uv environment (Python 2 compat via Python 3)"
cd "$WORKDIR/impls/original-python"
uv venv --python 3.10 .venv
# Original PhyloWGS deps
uv pip install --quiet numpy scipy ete3 pyvcf3 2>&1 | tail -2

# ── Clone: optimized Python (refactor main branch) ───────────────────────────
echo "[optimized-python] Cloning from PhyloWGS_refactor:main"
if [[ -d "$WORKDIR/impls/optimized-python" ]]; then
    git -C "$WORKDIR/impls/optimized-python" pull --rebase --autostash 2>&1 | tail -2
else
    git clone --depth=1 --branch main "$REPO" "$WORKDIR/impls/optimized-python" 2>&1 | tail -3
fi
echo "[optimized-python] Building uv environment"
cd "$WORKDIR/impls/optimized-python"
uv venv --python 3.10 .venv
uv pip install --quiet numpy scipy 2>&1 | tail -2

# ── Clone: Go implementation (go/main — single source for cpu + gpu) ─────────
echo "[go] Cloning from PhyloWGS_refactor:go/main"
if [[ -d "$WORKDIR/impls/go-src" ]]; then
    git -C "$WORKDIR/impls/go-src" pull --rebase --autostash 2>&1 | tail -2
else
    git clone --depth=1 --branch go/main "$REPO" "$WORKDIR/impls/go-src" 2>&1 | tail -3
fi

source "$WORKDIR/env.sh" 2>/dev/null || true

# ── Build: Go CPU binary ──────────────────────────────────────────────────────
echo "[go-cpu] Building CPU binary"
cd "$WORKDIR/impls/go-src"
go build -o "$WORKDIR/impls/go-cpu/phylowgs-cpu" . 2>&1
mkdir -p "$WORKDIR/impls/go-cpu"
go build -o "$WORKDIR/impls/go-cpu/phylowgs-cpu" . 2>&1
echo "[go-cpu] Built: $(ls -lh $WORKDIR/impls/go-cpu/phylowgs-cpu | awk '{print $5, $9}')"

# ── Build: Go GPU (CUDA) binary ───────────────────────────────────────────────
mkdir -p "$WORKDIR/impls/go-gpu"
echo "[go-gpu] Building CUDA binary"
if command -v nvcc &>/dev/null; then
    cd "$WORKDIR/impls/go-src/cuda"
    bash build.sh 2>&1 | tail -5
    cd "$WORKDIR/impls/go-src"
    CGO_ENABLED=1 LD_LIBRARY_PATH="$WORKDIR/impls/go-src/cuda" \
        go build -tags cuda -o "$WORKDIR/impls/go-gpu/phylowgs-gpu" . 2>&1
    # Copy shared lib next to binary
    cp cuda/liblikelihood.so "$WORKDIR/impls/go-gpu/" 2>/dev/null || true
    echo "[go-gpu] CUDA build OK: $(ls -lh $WORKDIR/impls/go-gpu/phylowgs-gpu | awk '{print $5, $9}')"
else
    echo "[go-gpu] WARNING: nvcc not found — GPU implementation unavailable"
    echo "CUDA_UNAVAILABLE" > "$WORKDIR/impls/go-gpu/.skip"
fi

# ── Copy conversion helper ────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cp "$SCRIPT_DIR/convert_inputs.py" "$WORKDIR/convert_inputs.py"

# ── Write implementation manifest ────────────────────────────────────────────
cat > "$WORKDIR/implementations.tsv" << EOF
name	type	binary_or_script	partition	extra_args	notes
original-python	python	evolve.py	cmobic_cpu		Original morrislab/phylowgs
optimized-python	python	evolve.py	cmobic_cpu		Python 3 refactor (this repo main)
go-cpu	go	phylowgs-cpu	cmobic_cpu	--no-gpu	Go rewrite, CPU only
go-gpu	go	phylowgs-gpu	gpu		Go rewrite, CUDA kernel
EOF

echo ""
echo "=== Setup complete ==="
echo "Workdir: $WORKDIR"
echo ""
echo "Implementations:"
column -t "$WORKDIR/implementations.tsv"
echo ""
echo "Next: bash submit_benchmark.sh samples.csv --workdir $WORKDIR"
