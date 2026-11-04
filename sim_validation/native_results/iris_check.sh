#!/bin/bash
#SBATCH --partition=cmobic_cpu
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --job-name=native_check
#SBATCH --output=/data1/core005/ccs/orgeraj/go_runs/native_check/check.log
#
# Real-data equivalence: `phylowgs-go write-results` vs write_results_go.py
# (python2 reference in the SIF) on the 20 seed_rerun samples.
# Patient data: outputs stay on Iris; this script prints only k, counts and
# ok/FAIL/densities_undefined -- never sample names, mutation ids or values.
set -uo pipefail

GO_RUNS=/data1/core005/ccs/orgeraj/go_runs
TASKS=$GO_RUNS/posterior_cmp/seed_rerun/tasks.tsv
RUNS=$GO_RUNS/posterior_cmp/seed_rerun/runs
SIF=$GO_RUNS/m1xcheck/phylowgs-1.5.2-msk.sif
SRC=${SRC:-$GO_RUNS/gosrc_flush}              # checkout of feat/go-native-results
BIN=${BIN:-$(ls -t $GO_RUNS/bin/phylowgs-go-native-* | head -1)}
PY3=/home/orgeraj/miniforge3/bin/python3
OUT=$GO_RUNS/native_check
FLAGS="--max-multiprimary 1.0 --include-multiprimary --include-ssm-names"
REF=$SRC/testdata/native_results
CMP=$SRC/sim_validation/native_results/compare_results.py

[ -x "$BIN" ] || { echo "no binary: $BIN"; exit 2; }
echo "binary=$(basename "$BIN") src=$(git -C "$SRC" rev-parse --short HEAD)"
mkdir -p "$OUT"
n_ok=0; n_dens=0; n_fail=0

# tasks.tsv: header row, then i k cls n_ssm seed ssm cnv (tab-separated)
for k in $(awk -F'\t' 'NR>1 {print $2}' "$TASKS"); do
  chains=$RUNS/k$k/chains
  if [ ! -s "$chains/trees.zip" ] || [ ! -s "$chains/mutlist.json" ]; then
    echo "k$k FAIL missing_inputs"; n_fail=$((n_fail+1)); continue
  fi
  w=$OUT/k$k; rm -rf "$w"; mkdir -p "$w/go" "$w/py"

  "$BIN" write-results $FLAGS "k$k" "$chains/trees.zip" "$chains/mutlist.json" \
      "$w/go/summ.json.gz" "$w/go/muts.json.gz" "$w/go/mutass.zip" >"$w/go/log" 2>&1
  go_rc=$?
  singularity exec --bind "$GO_RUNS" "$SIF" python2 "$REF/write_results_go.py" $FLAGS \
      "k$k" "$chains/trees.zip" "$chains/mutlist.json" \
      "$w/py/summ.json.gz" "$w/py/muts.json.gz" "$w/py/mutass.zip" >"$w/py/log" 2>&1
  py_rc=$?
  if [ $go_rc -ne 0 ] || [ $py_rc -ne 0 ]; then
    echo "k$k FAIL engine_rc go=$go_rc py=$py_rc"; n_fail=$((n_fail+1)); continue
  fi

  # decode both to pretty JSON with the reference's own dumper (same as goldens)
  dec_ok=1
  for e in go py; do
    singularity exec --bind "$GO_RUNS" "$SIF" python2 "$REF/dump_expected.py" "$w/$e" "$w/$e" \
        >"$w/$e/dump.log" 2>&1 || dec_ok=0
  done
  if [ $dec_ok -eq 0 ]; then
    echo "k$k FAIL decode"; n_fail=$((n_fail+1)); continue
  fi

  res=$("$PY3" -I "$CMP" --loose-int-float "k$k" "$w/go" "$w/py"); rc=$?
  echo "$res"
  case "$res" in
    *densities_undefined) n_dens=$((n_dens+1)) ;;
    *" ok")               [ $rc -eq 0 ] && n_ok=$((n_ok+1)) || n_fail=$((n_fail+1)) ;;
    *)                    n_fail=$((n_fail+1)) ;;
  esac
done

echo "SUMMARY ok=$n_ok densities_undefined=$n_dens fail=$n_fail"
[ $n_fail -eq 0 ]
