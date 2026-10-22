#!/bin/bash
# One sample of the per-sample-seed Go rerun; array index = row i of tasks.tsv.
# Same binary and flags as test5 (-B 1000 -s 2500 -j 15); only -r differs.
#SBATCH --job-name=go_seedrerun
#SBATCH --partition=cmobic_cpu
#SBATCH --cpus-per-task=12
set -euo pipefail
OUT=${OUT:?set OUT to the rerun dir}
BIN=${BIN:-/data1/core005/ccs/orgeraj/go_runs/bin/phylowgs-go-e2b2df4}
IFS=$'\t' read -r I K CLS NSSM SEED SSM CNV < <(awk -F'\t' -v i=$SLURM_ARRAY_TASK_ID '$1==i' $OUT/tasks.tsv)
[ -n "${K:-}" ] || { echo "no row $SLURM_ARRAY_TASK_ID in tasks.tsv" >&2; exit 1; }
export GOMEMLIMIT=${GOMEMLIMIT:-8GiB}
D=$OUT/runs/k$K
mkdir -p $D && cd $D
echo "k=$K cls=$CLS n_ssm=$NSSM seed=$SEED"
$BIN -B 1000 -s 2500 -j 15 -r $SEED -O chains -D seedrerun_k$K $SSM $CNV
