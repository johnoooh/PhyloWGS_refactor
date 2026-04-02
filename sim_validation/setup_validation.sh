#!/usr/bin/env bash
# setup_validation.sh — One-time setup for simulation validation on HPC
#
# Generates fixtures locally, builds the Go binary, and prepares the work tree.
#
# Usage:
#   bash setup_validation.sh [--workdir /path/to/workdir] [--grid quick|default|thorough]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

WORKDIR="${SIM_VALIDATION_WORKDIR:-$(pwd)/sim_validation_workdir}"
GRID="default"
BASE_SEED=42
GO_VERSION="1.22.0"

while [[ $# -gt 0 ]]; do
    case $1 in
        --workdir) WORKDIR="$2"; shift 2 ;;
        --grid)    GRID="$2";    shift 2 ;;
        --seed)    BASE_SEED="$2"; shift 2 ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

echo "=== Simulation Validation Setup ==="
echo "Workdir:  $WORKDIR"
echo "Grid:     $GRID"
echo "Seed:     $BASE_SEED"
echo ""

mkdir -p "$WORKDIR"
cd "$WORKDIR"

# ── env.sh for SLURM jobs to source ─────────────────────────────────────────
> env.sh

# ── Go installation ──────────────────────────────────────────────────────────
GO_INSTALL_DIR="$WORKDIR/go_install"
need_go_install=false

if ! command -v go &>/dev/null; then
    need_go_install=true
else
    CURRENT_GO=$(go version | grep -oP 'go\K[0-9]+\.[0-9]+' | head -1)
    CURRENT_MAJOR=$(echo "$CURRENT_GO" | cut -d. -f1)
    CURRENT_MINOR=$(echo "$CURRENT_GO" | cut -d. -f2)
    if [[ "$CURRENT_MAJOR" -lt 1 ]] || \
       [[ "$CURRENT_MAJOR" -eq 1 && "$CURRENT_MINOR" -lt 21 ]]; then
        need_go_install=true
    else
        echo "[go] Found $(go version)"
        echo "export PATH=\"$(dirname "$(command -v go)"):\$PATH\"" >> env.sh
    fi
fi

if [[ "$need_go_install" == true ]]; then
    echo "[go] Installing Go $GO_VERSION"
    mkdir -p "$GO_INSTALL_DIR"
    curl -sSL "https://go.dev/dl/go${GO_VERSION}.linux-amd64.tar.gz" \
        | tar -C "$GO_INSTALL_DIR" -xz
    export PATH="$GO_INSTALL_DIR/go/bin:$PATH"
    echo "export PATH=\"$GO_INSTALL_DIR/go/bin:\$PATH\"" >> env.sh
fi

# ── Python environment (numpy only) ─────────────────────────────────────────
if ! command -v uv &>/dev/null; then
    echo "[uv] Installing uv"
    curl -LsSf https://astral.sh/uv/install.sh | sh
    export PATH="$HOME/.cargo/bin:$PATH"
    echo "export PATH=\"\$HOME/.cargo/bin:\$PATH\"" >> env.sh
else
    echo "[uv] Found at $(command -v uv)"
    echo "export PATH=\"$(dirname "$(command -v uv)"):\$PATH\"" >> env.sh
fi

echo "[python] Creating venv"
uv venv "$WORKDIR/.venv"
source "$WORKDIR/.venv/bin/activate"
uv pip install --quiet numpy 2>&1 | tail -2
echo "source \"$WORKDIR/.venv/bin/activate\"" >> env.sh

# ── Clone / update Go source ────────────────────────────────────────────────
REPO_URL="https://github.com/johnoooh/PhyloWGS_refactor.git"
GO_SRC="$WORKDIR/go-src"

if [[ -d "$GO_SRC" ]]; then
    echo "[go-src] Updating"
    git -C "$GO_SRC" pull --rebase --autostash 2>&1 | tail -2
else
    echo "[go-src] Cloning go-port branch"
    git clone --depth=1 --branch go-port "$REPO_URL" "$GO_SRC" 2>&1 | tail -3
fi

# ── Build Go binary ─────────────────────────────────────────────────────────
echo "[go] Building phylowgs-go (CPU)"
mkdir -p "$WORKDIR/bin"
cd "$GO_SRC"
source "$WORKDIR/env.sh" 2>/dev/null || true
go build -o "$WORKDIR/bin/phylowgs-go" . 2>&1
echo "[go] Built: $(ls -lh "$WORKDIR/bin/phylowgs-go" | awk '{print $5}')"

# ── Build GPU binary if CUDA available ───────────────────────────────────────
if command -v nvcc &>/dev/null; then
    echo "[go-gpu] Building CUDA binary"
    cd "$GO_SRC/cuda"
    bash build.sh 2>&1 | tail -3
    cd "$GO_SRC"
    CGO_ENABLED=1 LD_LIBRARY_PATH="$GO_SRC/cuda" \
        go build -tags cuda -o "$WORKDIR/bin/phylowgs-go-gpu" . 2>&1
    cp cuda/liblikelihood.so "$WORKDIR/bin/" 2>/dev/null || true
    echo "[go-gpu] Built OK"
else
    echo "[go-gpu] nvcc not found — GPU binary skipped"
fi

# ── Copy scripts ─────────────────────────────────────────────────────────────
cp "$SCRIPT_DIR/generate_fixtures.py" "$WORKDIR/"
cp "$SCRIPT_DIR/score_results.py"     "$WORKDIR/"

# ── Generate fixtures ────────────────────────────────────────────────────────
echo ""
echo "[fixtures] Generating simulation fixtures (grid=$GRID)"
python3 "$WORKDIR/generate_fixtures.py" \
    --outdir "$WORKDIR/fixtures" \
    --grid "$GRID" \
    --base-seed "$BASE_SEED"

# ── Create directories ───────────────────────────────────────────────────────
mkdir -p "$WORKDIR/results" "$WORKDIR/logs" "$WORKDIR/analysis"

echo ""
echo "=== Setup complete ==="
echo "Workdir:  $WORKDIR"
echo "Binary:   $WORKDIR/bin/phylowgs-go"
echo "Fixtures: $WORKDIR/fixtures/"
echo ""
echo "Next: bash submit_sim_validation.sh --workdir $WORKDIR"
