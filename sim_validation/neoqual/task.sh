#!/bin/bash
# One step for one NeoQual sample. Array index = line of tasks_<CLS>.txt, where
# CLS is a chunk tag: size class + chunk number (S0 = pilot chunk of class S).
#   STEP=go    run Go (per-sample seed, same flags as test5); chains are written
#              as each finishes, so a straggler hitting the limit loses only itself.
#              Reused Go runs (go_src != new) are only checked for 15 chains.
#              Records the Slurm task id in go_time/<k>.job for runtime.py.
#   STEP=conv  write_results.py on each of the 15 Python chains (no tree filters)
#   STEP=cmp   runtime.py -> out/runtime_<k>.json, then pool_compare.py -> out/<k>.json
#              (skipped if present; needs >= 14 Go chains and all 15 Python chains)
# STEP and CLS come from the env or, failing that, the job name nq_<STEP>_<CLS>
# (for submitters without --export, e.g. hpcgate submit -J nq_go_S0).
# Env: W (work dir, default below); optional BIN, SIF, PY, SRC.
#SBATCH --partition=cmobic_cpu
#SBATCH --output=/data1/core005/ccs/orgeraj/go_runs/nq1500/logs/%x_%a.log
set -euo pipefail
W=${W:-/data1/core005/ccs/orgeraj/go_runs/nq1500}
if [ -z "${STEP:-}" ] && [[ "${SLURM_JOB_NAME:-}" =~ ^nq_(go|conv|cmp)_([SMLX][0-9]+)$ ]]; then
    STEP=${BASH_REMATCH[1]}; CLS=${BASH_REMATCH[2]}
fi
: "${CLS:?set CLS or name the job nq_<step>_<cls>}" "${STEP:?set STEP or name the job nq_<step>_<cls>}"
BIN=${BIN:-/data1/core005/ccs/orgeraj/go_runs/bin/phylowgs-go-native-3b381c0}
SIF=${SIF:-/data1/core005/ccs/orgeraj/go_runs/m1xcheck/phylowgs-1.5.2-msk.sif}
PY=${PY:-/home/orgeraj/miniforge3/bin/python3}
# Not derived from $0: sbatch runs a spooled copy of this script.
SRC=${SRC:-/data1/core005/ccs/orgeraj/go_runs/posterior_cmp/src/sim_validation}
k=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" $W/tasks_$CLS.txt)
[ -n "$k" ] || { echo "no task line $SLURM_ARRAY_TASK_ID in tasks_$CLS.txt" >&2; exit 1; }
get() { $PY -c "import json,sys;e=json.load(open('$W/manifest.json'))[$k];v=e['$1'];print('\n'.join(v) if isinstance(v,list) else v)"; }
echo "k=$k step=$STEP host=$(hostname)"

case $STEP in
go)
    D=$(get go_dir)
    if [ "$(get go_src)" != new ]; then
        n=$(ls $D/chains/chain_*_trees.ndjson 2>/dev/null | wc -l)
        [ "$n" -eq 15 ] || { echo "reused Go run $D has $n/15 chains" >&2; exit 1; }
        echo "reused Go chains ($(get go_src))"; exit 0
    fi
    mkdir -p $D $W/go_time && cd $D
    [ -f DONE ] && { echo "already done"; exit 0; }
    echo "${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}" > $W/go_time/$k.job
    # GOMEMLIMIT a little under the job's memory so GC works before the OOM killer.
    export GOMEMLIMIT=$(( SLURM_MEM_PER_NODE * 85 / 100 ))MiB
    $BIN -B 1000 -s 2500 -j 15 -r $(get seed) -O chains -D nq_$k "$(get ssm)" "$(get cnv)"
    test -s chains/summary.json && touch DONE
    ;;
conv)
    D=$(get py_dir); mapfile -t SRCS < <(get py_src)
    for c in $(seq 0 14); do
        O=$D/chain_$c; mkdir -p $O; [ -f $O/DONE ] && continue
        (cd $O && singularity exec --bind $W --bind /data1/core006/ccs_pipelines/NeoQual_test:/data1/core006/ccs_pipelines/NeoQual_test:ro \
            $SIF python2 /usr/bin/phylowgs/write_results.py --include-multiprimary --keep-superclones --min-ssms 0 \
            p$k "${SRCS[$c]}" tree_summaries.json.gz mutlist.json.gz mutass.zip \
         && test -s tree_summaries.json.gz && test -s mutass.zip && touch DONE) \
         || echo "chain $c conversion failed" >&2
    done
    n=$(ls $D/chain_*/DONE | wc -l); [ "$n" -eq 15 ] || { echo "only $n/15 chains converted" >&2; exit 1; }
    ;;
cmp)
    G=$(get go_dir)/chains; D=$(get py_dir); mkdir -p $W/out
    ng=$(ls $G/chain_*_trees.ndjson 2>/dev/null | wc -l)
    np=$(ls $D/chain_*/DONE 2>/dev/null | wc -l)
    [ "$ng" -ge 14 ] && [ "$np" -eq 15 ] || { echo "not ready: go $ng chains, py $np/15" >&2; exit 1; }
    $PY $SRC/neoqual/runtime.py $W $k
    [ -s $W/out/$k.json ] && { echo "pool_compare already done"; exit 0; }
    $PY $SRC/pool_compare.py --go $G/chain_*_trees.ndjson --py $D/chain_* --pool 7 --draws 20 \
        --trees-per-chain 250 --max-ssms 1500 --out $W/out/$k.json
    ;;
*) echo "unknown STEP=$STEP" >&2; exit 1 ;;
esac
