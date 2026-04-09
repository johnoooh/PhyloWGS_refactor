#!/usr/bin/env bash
# Local trace baseline for Round 2 Phase 3 validation.
# Runs the Go port with --trace on a small slice while HPC validation runs.
#
# Three experiments:
#   A: HPC-parity baseline (4 fixtures × 1 chain × B=500/S=1000)
#   B: Convergence/burnin (4 fixtures × 1 chain × B=2000/S=5000)
#   C: Multi-chain variance (K3_C2 × 4 chains × B=500/S=1000)
#
# Output structure:
#   sim_validation/local_trace_runs/expA_baseline/<fixture>/{trace.ndjson, chain_0_trees.ndjson, ...}
#   sim_validation/local_trace_runs/expB_long/<fixture>/{trace.ndjson, ...}
#   sim_validation/local_trace_runs/expC_multichain/{trace.ndjson.0, trace.ndjson.1, ...}

set -euo pipefail

REPO_ROOT="/sessions/dreamy-pensive-knuth/mnt/phylowgs_cowork"
BIN="$REPO_ROOT/PhyloWGS_refactor/phylowgs-go"
FIXTURES_ROOT="$REPO_ROOT/simulation_validation/fixtures"
OUT_ROOT="$REPO_ROOT/PhyloWGS_refactor/sim_validation/local_trace_runs"

FIXTURES=(
  K3_S1_T200_M30_C2_rep0
  K3_S1_T200_M30_C0_rep0
  K5_S1_T200_M50_C2_rep0
  K5_S1_T200_M50_C0_rep0
)

START_TOTAL=$(date +%s)
echo "==> START $(date)"
echo "==> Binary: $BIN"
echo "==> Binary mtime: $(stat -c '%y' "$BIN")"

run_one() {
  local label="$1" fixture="$2" chains="$3" burnin="$4" samples="$5"
  local out_dir="$OUT_ROOT/$label/$fixture"
  mkdir -p "$out_dir"
  local ssm="$FIXTURES_ROOT/$fixture/ssm_data.txt"
  local cnv="$FIXTURES_ROOT/$fixture/cnv_data.txt"
  echo
  echo "  --- $label / $fixture (chains=$chains, B=$burnin, S=$samples)"
  local t0=$(date +%s)
  "$BIN" --no-gpu \
        -B "$burnin" -s "$samples" -j "$chains" \
        --trace "$out_dir/trace.ndjson" \
        -O "$out_dir" \
        "$ssm" "$cnv" \
    > "$out_dir/stdout.log" 2>&1
  local t1=$(date +%s)
  echo "      done in $((t1 - t0))s"
}

# ─────────── Experiment A: HPC parity ───────────
echo
echo "=== Experiment A: HPC-parity baseline (B=500, S=1000, 1 chain) ==="
for f in "${FIXTURES[@]}"; do
  run_one "expA_baseline" "$f" 1 500 1000
done

# ─────────── Experiment B: long convergence ───────────
echo
echo "=== Experiment B: long convergence (B=2000, S=5000, 1 chain) ==="
for f in "${FIXTURES[@]}"; do
  run_one "expB_long" "$f" 1 2000 5000
done

# ─────────── Experiment C: multi-chain variance ───────────
echo
echo "=== Experiment C: multi-chain variance K3_C2 (B=500, S=1000, 4 chains) ==="
run_one "expC_multichain" "K3_S1_T200_M30_C2_rep0" 4 500 1000

END_TOTAL=$(date +%s)
echo
echo "==> END $(date)  total wall=$((END_TOTAL - START_TOTAL))s"
