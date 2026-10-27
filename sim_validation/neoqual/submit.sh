#!/bin/bash
# Submit Go, Python-conversion and comparison arrays per size class.
# Go arrays are throttled so at most ~50 Go jobs run at once. cmp waits for the
# whole Go and conv arrays (afterany): a Go task that timed out on one straggler
# still leaves 14 finished chains, which cmp accepts.
# Usage: submit.sh <work_dir>   (manifest.json + tasks_<cls>.txt from mk_manifest.py)
set -euo pipefail
W=$(realpath "$1")
T=$(dirname "$(realpath "$0")")/task.sh
mkdir -p $W/logs
#        cls  go-throttle go-mem go-time
for spec in "S 10 10G 6:00:00" "M 15 16G 12:00:00" "L 15 32G 24:00:00" "X 10 48G 72:00:00"; do
    read -r CLS THR MEM TIME <<< "$spec"
    n=$(grep -c . $W/tasks_$CLS.txt || true)
    [ "$n" -gt 0 ] || continue
    A="--array=0-$((n-1))"
    go=$(sbatch --parsable $A%$THR -J nq_go_$CLS -c 12 --mem=$MEM -t $TIME \
        --export=ALL,W=$W,CLS=$CLS,STEP=go -o $W/logs/go_${CLS}_%a.log $T)
    cv=$(sbatch --parsable $A%20 -J nq_conv_$CLS -c 1 --mem=8G -t 4:00:00 \
        --export=ALL,W=$W,CLS=$CLS,STEP=conv -o $W/logs/conv_${CLS}_%a.log $T)
    cm=$(sbatch --parsable $A%20 -J nq_cmp_$CLS -c 1 --mem=16G -t 8:00:00 \
        --dependency=afterany:$go:$cv \
        --export=ALL,W=$W,CLS=$CLS,STEP=cmp -o $W/logs/cmp_${CLS}_%a.log $T)
    echo "$CLS n=$n go=$go conv=$cv cmp=$cm"
done
