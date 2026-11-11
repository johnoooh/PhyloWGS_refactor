#!/bin/bash
# Submit Go, Python-conversion and comparison arrays for chunk tags, from the Mac
# through hpcgate (task.sh reads STEP/CLS from the job name nq_<step>_<tag>).
# Go arrays are throttled per class (S 10, M 15, L 15, X 10: ~50 Go jobs at once),
# and a class's chunks are chained (go of chunk i waits on go of chunk i-1), so
# the throttle holds across chunks. cmp waits for its chunk's whole Go and conv
# arrays (afterany): a Go task that timed out on one straggler still leaves 14
# finished chains, which cmp accepts.
# Usage: submit.sh TAG...   e.g. submit.sh S1 S2 M1 L1
# Prints: TAG n=<tasks> go=<id> conv=<id> cmp=<id>
set -euo pipefail
W=/data1/core005/ccs/orgeraj/go_runs/nq1500
T=/data1/core005/ccs/orgeraj/go_runs/posterior_cmp/src/sim_validation/neoqual/task.sh
jid() { grep -oE 'Submitted batch job [0-9]+' | grep -oE '[0-9]+$'; }
PREV_S= PREV_M= PREV_L= PREV_X=
for TAG in "$@"; do
    CLS=${TAG:0:1}
    case $CLS in
        S) THR=10 MEM=10G TIME=6:00:00 PREV=$PREV_S ;;
        M) THR=15 MEM=16G TIME=12:00:00 PREV=$PREV_M ;;
        L) THR=15 MEM=32G TIME=24:00:00 PREV=$PREV_L ;;
        X) THR=10 MEM=48G TIME=72:00:00 PREV=$PREV_X ;;
        *) echo "bad tag $TAG" >&2; exit 1 ;;
    esac
    n=$(hpcgate run "grep -c . $W/tasks_$TAG.txt" | tail -1)
    [ "$n" -gt 0 ] || { echo "$TAG: no tasks" >&2; exit 1; }
    A="0-$((n-1))"
    go=$(hpcgate submit $T -D $W -p cmobic_cpu -a "$A%$THR" -J nq_go_$TAG -c 12 --mem $MEM -t $TIME \
        ${PREV:+-d afterany:$PREV} | jid)
    cv=$(hpcgate submit $T -D $W -p cmobic_cpu -a "$A%20" -J nq_conv_$TAG -c 1 --mem 8G -t 4:00:00 | jid)
    cm=$(hpcgate submit $T -D $W -p cmobic_cpu -a "$A%20" -J nq_cmp_$TAG -c 1 --mem 16G -t 8:00:00 \
        -d afterany:$go:$cv | jid)
    [ -n "$go" ] && [ -n "$cv" ] && [ -n "$cm" ] || { echo "$TAG: a submission returned no job id" >&2; exit 1; }
    case $CLS in S) PREV_S=$go ;; M) PREV_M=$go ;; L) PREV_L=$go ;; X) PREV_X=$go ;; esac
    echo "$TAG n=$n go=$go conv=$cv cmp=$cm"
done
