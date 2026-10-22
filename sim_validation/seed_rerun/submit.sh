#!/bin/bash
# Submit the per-sample-seed Go rerun, one array per size class (test5 limits, +margin).
# Usage: submit.sh <rerun_out_dir>   (tasks.tsv from select_samples.py must exist)
set -euo pipefail
OUT=$(realpath "$1")
HERE=$(dirname "$(realpath "$0")")
mkdir -p $OUT/logs
for spec in "A 10G 3:00:00" "B 18G 6:00:00" "C 34G 15:00:00" "D 34G 24:00:00"; do
    read -r CLS MEM TIME <<< "$spec"
    IDX=$(awk -F'\t' -v c=$CLS '$3==c{print $1}' $OUT/tasks.tsv | paste -sd,)
    [ -n "$IDX" ] || continue
    sbatch --array=$IDX --mem=$MEM --time=$TIME --export=ALL,OUT=$OUT \
        --output=$OUT/logs/%x_%a.log $HERE/go_task.sh
done
