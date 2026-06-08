#!/usr/bin/env bash
# benchmark_matched.sh — MATCHED-WORK wall-time benchmark: Go vs Python
# ============================================================================
# WHY THIS SCRIPT EXISTS (read before running)
# ----------------------------------------------------------------------------
# A prior speed comparison was unfair: it timed the Go sampler running 8 chains
# x 2500 samples against a Python run with FEWER chains/samples. That is not a
# like-for-like measurement — it inflates or deflates the apparent speedup
# depending on which side did more work.
#
# This script guarantees MATCHED WORK:
#   - SAME burn-in (-B)        for both impls
#   - SAME samples (-s)        for both impls
#   - SAME number of chains    (Go -j  ==  Python -n)
#   - SAME MH iterations       (Go -i 5000; Python evolve.py is fixed at 5000)
#   - SAME input fixtures      (one shared fixtures dir, one base-seed)
#
# It then records wall_seconds for each impl per fixture into a tidy CSV and a
# JSON sidecar so the speed comparison is honest and reproducible.
#
# This script does NOT change any Go or Python source and does NOT alter the
# sampler's statistical behavior — it only measures wall time.
#
# ----------------------------------------------------------------------------
# USAGE
# ----------------------------------------------------------------------------
#   bash benchmark_matched.sh \
#       --workdir /path/to/matched_sweep \
#       --sif     /path/to/phylowgs.sif \
#       [--fixtures NAME1,NAME2] [--max-fixtures N] \
#       [--burnin 1000] [--samples 2500] [--chains 8] [--mh-iters 5000] \
#       [--repeats 1] [--out bench_matched.csv] [--dry-run]
#
# Options:
#   --workdir DIR      Workdir from setup_validation.sh (needs bin/phylowgs-go
#                      and fixtures/). REQUIRED.
#   --sif PATH         PhyloWGS singularity image for Python. REQUIRED unless
#                      --go-only.
#   --fixtures LIST    Comma-separated fixture names to time. Default: pick the
#                      smallest few from the manifest (fast, representative).
#   --max-fixtures N   If --fixtures not given, cap auto-selected count (def 3).
#   --burnin N         Burn-in for BOTH impls (default: 1000).
#   --samples N        Samples for BOTH impls (default: 2500).
#   --chains N         Chains for BOTH impls: Go -j and Python -n (default: 8).
#   --mh-iters N       MH iters per MCMC iter for Go -i (default: 5000, the
#                      Python evolve.py fixed value — keep at 5000 for a match).
#   --repeats N        Repeat each timing N times for variance (default: 1).
#   --bind-paths P     Singularity --bind paths (default: auto from workdir).
#   --phylowgs-dir D   Path to multievolve.py inside the container
#                      (default: /usr/bin/phylowgs).
#   --out FILE         CSV output path (default: $WORKDIR/bench_matched.csv).
#                      A .json sidecar is written alongside.
#   --go-only          Time only Go (skips the Python half; not a comparison).
#   --dry-run          Print the commands; run nothing.
#
# OUTPUT (CSV columns):
#   fixture,impl,burnin,samples,chains,mh_iters,repeat,wall_seconds,exit_code
# Plus a JSON array sidecar with the same records for easy programmatic use.
# ============================================================================

set -euo pipefail

# ── Defaults ─────────────────────────────────────────────────────────────────
WORKDIR=""
PHYLOWGS_SIF=""
FIXTURES_CSV=""
MAX_FIXTURES=3
BURNIN=1000
SAMPLES=2500
CHAINS=8
MH_ITERS=5000
REPEATS=1
BIND_PATHS=""
PHYLOWGS_DIR="/usr/bin/phylowgs"
OUT_CSV=""
GO_ONLY=false
DRY_RUN=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --workdir)      WORKDIR="$2";      shift 2 ;;
        --sif)          PHYLOWGS_SIF="$2"; shift 2 ;;
        --fixtures)     FIXTURES_CSV="$2"; shift 2 ;;
        --max-fixtures) MAX_FIXTURES="$2"; shift 2 ;;
        --burnin)       BURNIN="$2";       shift 2 ;;
        --samples)      SAMPLES="$2";      shift 2 ;;
        --chains)       CHAINS="$2";       shift 2 ;;
        --mh-iters)     MH_ITERS="$2";     shift 2 ;;
        --repeats)      REPEATS="$2";      shift 2 ;;
        --bind-paths)   BIND_PATHS="$2";   shift 2 ;;
        --phylowgs-dir) PHYLOWGS_DIR="$2"; shift 2 ;;
        --out)          OUT_CSV="$2";      shift 2 ;;
        --go-only)      GO_ONLY=true;      shift ;;
        --dry-run)      DRY_RUN=true;      shift ;;
        -h|--help)      sed -n '2,75p' "${BASH_SOURCE[0]}"; exit 0 ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

# ── Validate ─────────────────────────────────────────────────────────────────
[[ -z "$WORKDIR" ]] && { echo "ERROR: --workdir is required."; exit 1; }
[[ ! -d "$WORKDIR" ]] && { echo "ERROR: workdir not found: $WORKDIR — run setup_validation.sh first"; exit 1; }

BINARY="$WORKDIR/bin/phylowgs-go"
[[ ! -x "$BINARY" ]] && { echo "ERROR: Go binary not found/executable: $BINARY"; exit 1; }

FIXTURES_DIR="$WORKDIR/fixtures"
MANIFEST="$FIXTURES_DIR/manifest.json"
[[ ! -f "$MANIFEST" ]] && { echo "ERROR: manifest not found: $MANIFEST"; exit 1; }

if [[ "$GO_ONLY" == false ]]; then
    [[ -z "$PHYLOWGS_SIF" ]]     && { echo "ERROR: --sif required (or pass --go-only)."; exit 1; }
    [[ ! -f "$PHYLOWGS_SIF" ]]   && { echo "ERROR: SIF not found: $PHYLOWGS_SIF"; exit 1; }
fi

[[ -z "$OUT_CSV" ]] && OUT_CSV="$WORKDIR/bench_matched.csv"
OUT_JSON="${OUT_CSV%.csv}.json"

# Auto-detect bind paths (mirror submit_sim_validation_lsf.sh logic).
if [[ -z "$BIND_PATHS" && "$GO_ONLY" == false ]]; then
    _top_mount=$(echo "$WORKDIR" | cut -d/ -f1-2)
    BIND_PATHS="$_top_mount"
    _sif_mount=$(echo "$(readlink -f "$PHYLOWGS_SIF")" | cut -d/ -f1-2)
    if [[ "$_sif_mount" != "$_top_mount" ]]; then
        BIND_PATHS="${BIND_PATHS},${_sif_mount}"
    fi
fi

# ── Resolve fixture list ─────────────────────────────────────────────────────
declare -a FIXTURES
if [[ -n "$FIXTURES_CSV" ]]; then
    IFS=',' read -r -a FIXTURES <<< "$FIXTURES_CSV"
else
    # Auto-pick the smallest fixtures (by M in the name) for a fast, fair run.
    while IFS= read -r name; do
        FIXTURES+=("$name")
    done < <(python3 -c "
import json, re
with open('$MANIFEST') as f:
    m = json.load(f)
def mcount(e):
    g = re.search(r'_M(\d+)_', e['name'])
    return int(g.group(1)) if g else 10**9
for e in sorted(m, key=mcount)[:$MAX_FIXTURES]:
    print(e['name'])
")
fi
[[ ${#FIXTURES[@]} -eq 0 ]] && { echo "ERROR: no fixtures selected."; exit 1; }

echo "=== MATCHED-WORK Timing Benchmark: Go vs Python ==="
echo "Workdir:    $WORKDIR"
echo "Binary:     $BINARY"
echo "Fixtures:   ${FIXTURES[*]}"
echo "MCMC:       B=$BURNIN s=$SAMPLES chains=$CHAINS mh_iters=$MH_ITERS  (IDENTICAL both impls)"
echo "Repeats:    $REPEATS"
echo "Impls:      $([[ "$GO_ONLY" == true ]] && echo "go-cpu" || echo "go-cpu original-python")"
echo "CSV out:    $OUT_CSV"
echo "JSON out:   $OUT_JSON"
echo ""

if [[ "$DRY_RUN" != true ]]; then
    echo "fixture,impl,burnin,samples,chains,mh_iters,repeat,wall_seconds,exit_code" > "$OUT_CSV"
    : > "$OUT_JSON.tmp"
fi

# ── Timing helper: runs a command, prints elapsed integer seconds ────────────
time_run() {
    # echoes "ELAPSED EXITCODE"
    local start end ec
    start=$(date +%s)
    if "$@" > /dev/null 2>&1; then ec=0; else ec=$?; fi
    end=$(date +%s)
    echo "$((end - start)) $ec"
}

record() {
    # fixture impl repeat wall ec
    local fx="$1" impl="$2" rep="$3" wall="$4" ec="$5"
    if [[ "$DRY_RUN" == true ]]; then return; fi
    echo "${fx},${impl},${BURNIN},${SAMPLES},${CHAINS},${MH_ITERS},${rep},${wall},${ec}" >> "$OUT_CSV"
    cat >> "$OUT_JSON.tmp" << JSON
{"fixture":"${fx}","impl":"${impl}","burnin":${BURNIN},"samples":${SAMPLES},"chains":${CHAINS},"mh_iters":${MH_ITERS},"repeat":${rep},"wall_seconds":${wall},"exit_code":${ec}}
JSON
}

# Make the installed Go/uv env available if present.
source "$WORKDIR/env.sh" 2>/dev/null || true

for fx in "${FIXTURES[@]}"; do
    FX_DIR="$FIXTURES_DIR/$fx"
    SSM="$FX_DIR/ssm_data.txt"
    CNV="$FX_DIR/cnv_data.txt"
    [[ ! -f "$SSM" ]] && { echo "  SKIP $fx: missing $SSM"; continue; }

    for ((rep=1; rep<=REPEATS; rep++)); do
        echo "--- fixture=$fx repeat=$rep ---"

        # ── Go (matched work) ────────────────────────────────────────────────
        GO_OUT="$WORKDIR/bench_results/go-cpu/${fx}_rep${rep}"
        echo "  [go-cpu]  -B $BURNIN -s $SAMPLES -j $CHAINS -i $MH_ITERS"
        if [[ "$DRY_RUN" == true ]]; then
            echo "    [dry-run] $BINARY --no-gpu -B $BURNIN -s $SAMPLES -j $CHAINS -i $MH_ITERS -O $GO_OUT $SSM $CNV"
        else
            rm -rf "$GO_OUT"; mkdir -p "$GO_OUT"
            read -r gw ge <<< "$(time_run "$BINARY" --no-gpu \
                -B "$BURNIN" -s "$SAMPLES" -j "$CHAINS" -i "$MH_ITERS" \
                -O "$GO_OUT" "$SSM" "$CNV")"
            echo "    go-cpu wall=${gw}s exit=${ge}"
            record "$fx" "go-cpu" "$rep" "$gw" "$ge"
        fi

        # ── Python (IDENTICAL work) ──────────────────────────────────────────
        if [[ "$GO_ONLY" == false ]]; then
            PY_OUT="$WORKDIR/bench_results/original-python/${fx}_rep${rep}"
            echo "  [python]  -B $BURNIN -s $SAMPLES -n $CHAINS (evolve.py MH=5000 fixed)"
            if [[ "$DRY_RUN" == true ]]; then
                echo "    [dry-run] singularity exec --bind $BIND_PATHS $PHYLOWGS_SIF \\"
                echo "                python2 $PHYLOWGS_DIR/multievolve.py -n $CHAINS -B $BURNIN -s $SAMPLES \\"
                echo "                -O $PY_OUT/chains --ssms $SSM --cnvs $CNV"
            else
                rm -rf "$PY_OUT"; mkdir -p "$PY_OUT"
                read -r pw pe <<< "$(time_run singularity exec --bind "$BIND_PATHS" "$PHYLOWGS_SIF" \
                    python2 "$PHYLOWGS_DIR/multievolve.py" \
                    -n "$CHAINS" -B "$BURNIN" -s "$SAMPLES" \
                    -O "$PY_OUT/chains" --ssms "$SSM" --cnvs "$CNV")"
                echo "    python wall=${pw}s exit=${pe}"
                record "$fx" "original-python" "$rep" "$pw" "$pe"
            fi
        fi
    done
done

# ── Finalize JSON (wrap NDJSON tmp into a JSON array) ────────────────────────
if [[ "$DRY_RUN" != true ]]; then
    python3 - "$OUT_JSON.tmp" "$OUT_JSON" << 'PYEOF'
import json, sys
src, dst = sys.argv[1], sys.argv[2]
recs = []
with open(src) as f:
    for line in f:
        line = line.strip()
        if line:
            recs.append(json.loads(line))
with open(dst, "w") as f:
    json.dump(recs, f, indent=2)
PYEOF
    rm -f "$OUT_JSON.tmp"

    echo ""
    echo "=== Matched-work benchmark complete ==="
    echo "CSV:  $OUT_CSV"
    echo "JSON: $OUT_JSON"
    echo ""
    echo "Quick speedup summary (mean wall per impl):"
    python3 - "$OUT_CSV" << 'PYEOF'
import csv, sys, statistics
from collections import defaultdict
rows = list(csv.DictReader(open(sys.argv[1])))
by = defaultdict(list)
for r in rows:
    by[r["impl"]].append(float(r["wall_seconds"]))
means = {k: statistics.mean(v) for k, v in by.items()}
for k in sorted(means):
    print(f"  {k:16s} mean wall = {means[k]:.1f}s  (n={len(by[k])})")
if "go-cpu" in means and "original-python" in means and means["go-cpu"] > 0:
    print(f"  speedup (py/go) = {means['original-python']/means['go-cpu']:.2f}x")
PYEOF
else
    echo ""
    echo "[dry-run] No jobs run; no output written."
fi
