#!/usr/bin/env bash
# submit_sim_validation.sh — Submit SLURM array jobs to run Go port AND original
# Python on simulated fixtures, then score and compare.
#
# Runs two implementations per fixture:
#   - go-cpu: Go port (from setup_validation.sh binary)
#   - original-python: Original PhyloWGS via singularity container
#
# Jobs are submitted as Slurm array jobs grouped by M-tier:
#   - Each tier gets its own wall time (derived from P95 × 2.33 × 1.05 safety)
#   - Tiers under 2h are submitted to cmobic_preempt to reduce fairshare
#   - Longer tiers go to cmobic_cpu (or --partition override)
#
# Usage:
#   bash submit_sim_validation.sh --workdir /path/to/workdir --sif /path/to/phylowgs.sif [options]
#
# Options:
#   --workdir DIR       Workdir from setup_validation.sh
#   --sif PATH          PhyloWGS singularity image (.sif) for original Python
#   --burnin N          MCMC burn-in (default: 1000)
#   --samples N         MCMC samples (default: 2500)
#   --chains N          Parallel chains for Go (default: 8)
#   --mem MB            Memory in MB (default: 8000)
#   --partition NAME    SLURM partition for long jobs (default: cmobic_cpu)
#   --preempt NAME      SLURM partition for short jobs (default: cmobic_preempt)
#   --preempt-cutoff S  Max seconds for preemptable (default: 7200 = 2h)
#   --bind-paths PATHS  Singularity --bind paths (default: auto-detect from workdir)
#   --phylowgs-dir DIR  Path to evolve.py inside container (default: /opt/phylowgs)
#   --go-only           Only run Go port (skip original Python)
#   --dry-run           Print jobs without submitting

set -euo pipefail

WORKDIR="${SIM_VALIDATION_WORKDIR:-$(pwd)/sim_validation_workdir}"
BURNIN=1000
SAMPLES=2500
CHAINS=8
MEM_MB=8000
PARTITION="cmobic_cpu"
PREEMPT_PARTITION="cmobic_preempt"
PREEMPT_CUTOFF=7200   # 2 hours in seconds
PY_CHAINS=8
DRY_RUN=false
GO_ONLY=false
PHYLOWGS_SIF=""
BIND_PATHS=""
PHYLOWGS_DIR="/usr/bin/phylowgs"
CHAIN_INCLUSION_FACTOR="1.1"

# ── Wall-time tiers based on mutation count (M) ─────────────────────────────
# Derived from P95 of observed runtimes at B=500/s=1000, scaled by 2.33×
# (for new B=1000/s=2500) with 1.05× safety margin, rounded up to 15 min.
# Python tiers include +20 min headroom for write_results.py post-processing.
#
# Returns: "HH:MM:SS SECONDS" (two space-separated tokens)
get_go_tier() {
    local m_val
    if [[ "$1" =~ _M([0-9]+)_ ]]; then m_val="${BASH_REMATCH[1]}"; else echo "8:15:00 29700"; return; fi
    if   (( m_val <= 30 ));  then echo "0:30:00 1800"
    elif (( m_val <= 50 ));  then echo "0:45:00 2700"
    elif (( m_val <= 100 )); then echo "1:15:00 4500"
    elif (( m_val <= 150 )); then echo "4:30:00 16200"
    elif (( m_val <= 250 )); then echo "8:15:00 29700"
    else                          echo "7:45:00 27900"
    fi
}

get_py_tier() {
    # Python now runs multievolve.py with 8 chains. Each chain runs evolve.py
    # sequentially within multievolve's subprocess pool. Wall time is roughly
    # 8× single-chain evolve.py + write_results.py overhead.
    local m_val
    if [[ "$1" =~ _M([0-9]+)_ ]]; then m_val="${BASH_REMATCH[1]}"; else echo "72:00:00 259200"; return; fi
    if   (( m_val <= 30 ));  then echo "10:00:00 36000"
    elif (( m_val <= 50 ));  then echo "22:00:00 79200"
    elif (( m_val <= 100 )); then echo "44:00:00 158400"
    elif (( m_val <= 150 )); then echo "50:00:00 180000"
    elif (( m_val <= 250 )); then echo "56:00:00 201600"
    else                          echo "72:00:00 259200"
    fi
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --workdir)        WORKDIR="$2";           shift 2 ;;
        --sif)            PHYLOWGS_SIF="$2";      shift 2 ;;
        --burnin)         BURNIN="$2";            shift 2 ;;
        --samples)        SAMPLES="$2";           shift 2 ;;
        --chains)         CHAINS="$2";            shift 2 ;;
        --mem)            MEM_MB="$2";            shift 2 ;;
        --partition)      PARTITION="$2";         shift 2 ;;
        --preempt)        PREEMPT_PARTITION="$2"; shift 2 ;;
        --preempt-cutoff) PREEMPT_CUTOFF="$2";    shift 2 ;;
        --bind-paths)     BIND_PATHS="$2";        shift 2 ;;
        --phylowgs-dir)   PHYLOWGS_DIR="$2";      shift 2 ;;
        --chain-inclusion-factor) CHAIN_INCLUSION_FACTOR="$2"; shift 2 ;;
        --go-only)        GO_ONLY=true;           shift ;;
        --dry-run)        DRY_RUN=true;           shift ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

# ── Validate workdir ─────────────────────────────────────────────────────────
[[ ! -d "$WORKDIR" ]]           && { echo "Workdir not found: $WORKDIR — run setup_validation.sh first"; exit 1; }
[[ ! -f "$WORKDIR/env.sh" ]]    && { echo "env.sh not found — run setup_validation.sh first"; exit 1; }
[[ ! -d "$WORKDIR/fixtures" ]]  && { echo "No fixtures — run setup_validation.sh first"; exit 1; }

BINARY="$WORKDIR/bin/phylowgs-go"
[[ ! -x "$BINARY" ]] && { echo "Go binary not found: $BINARY"; exit 1; }

# Find SIF if not provided
if [[ "$GO_ONLY" == false && -z "$PHYLOWGS_SIF" ]]; then
    if [[ -n "${NXF_SINGULARITY_CACHEDIR:-}" ]] && ls "${NXF_SINGULARITY_CACHEDIR}"/*phylowgs* &>/dev/null 2>&1; then
        PHYLOWGS_SIF=$(ls "${NXF_SINGULARITY_CACHEDIR}"/*phylowgs* | head -1)
    fi
    if [[ -z "$PHYLOWGS_SIF" ]]; then
        echo "WARNING: No --sif provided and no PhyloWGS SIF found. Use --go-only or provide --sif."
        echo "  Falling back to --go-only mode."
        GO_ONLY=true
    fi
fi

if [[ "$GO_ONLY" == false ]]; then
    [[ ! -f "$PHYLOWGS_SIF" ]] && { echo "SIF not found: $PHYLOWGS_SIF"; exit 1; }
    echo "SIF image: $PHYLOWGS_SIF"
    echo "PhyloWGS dir (in container): $PHYLOWGS_DIR"
fi

# Auto-detect bind paths from workdir if not specified
if [[ -z "$BIND_PATHS" ]]; then
    # Bind the top-level mount of the workdir (e.g. /scratch, /data1, /juno)
    _top_mount=$(echo "$WORKDIR" | cut -d/ -f1-2)
    BIND_PATHS="$_top_mount"
    # Also bind the SIF parent if it's on a different mount
    if [[ "$GO_ONLY" == false ]]; then
        _sif_mount=$(echo "$(readlink -f "$PHYLOWGS_SIF")" | cut -d/ -f1-2)
        if [[ "$_sif_mount" != "$_top_mount" ]]; then
            BIND_PATHS="${BIND_PATHS},${_sif_mount}"
        fi
    fi
    echo "Bind paths (auto): $BIND_PATHS"
fi

FIXTURES_DIR="$WORKDIR/fixtures"
RESULTS_DIR="$WORKDIR/results"
LOGS_DIR="$WORKDIR/logs"
ARRAY_DIR="$WORKDIR/array_lists"
mkdir -p "$RESULTS_DIR" "$LOGS_DIR" "$ARRAY_DIR"

# ── Read manifest ────────────────────────────────────────────────────────────
MANIFEST="$FIXTURES_DIR/manifest.json"
[[ ! -f "$MANIFEST" ]] && { echo "manifest.json not found in fixtures/"; exit 1; }

FIXTURE_NAMES=()
while IFS= read -r name; do
    FIXTURE_NAMES+=("$name")
done < <(python3 -c "
import json, sys
with open('$MANIFEST') as f:
    manifest = json.load(f)
for entry in manifest:
    print(entry['name'])
")

N_FIXTURES="${#FIXTURE_NAMES[@]}"
IMPLS="go-cpu"
[[ "$GO_ONLY" == false ]] && IMPLS="go-cpu original-python"

echo "=== Simulation Validation Job Submission (Array Jobs) ==="
echo "Workdir:         $WORKDIR"
echo "Fixtures:        $N_FIXTURES"
echo "Implementations: $IMPLS"
echo "Config:          B=$BURNIN s=$SAMPLES j=$CHAINS (Go); B=$BURNIN s=$SAMPLES n=$PY_CHAINS (Python multievolve)"
echo "Partitions:      $PARTITION (long), $PREEMPT_PARTITION (< ${PREEMPT_CUTOFF}s)"
echo "Memory:          ${MEM_MB}MB"
echo ""

# ── Helper: group fixtures by time tier for an impl ─────────────────────────
# Populates associative arrays: TIER_COUNTS, TIER_TIMES, TIER_PARTS
# and writes per-tier fixture list files.
declare -A GO_TIER_COUNTS GO_TIER_TIMES GO_TIER_PARTS
declare -A PY_TIER_COUNTS PY_TIER_TIMES PY_TIER_PARTS

group_fixtures() {
    local impl="$1"          # "go-cpu" or "original-python"
    local tier_fn="$2"       # "get_go_tier" or "get_py_tier"
    local skip_file="$3"     # file to check for skip (e.g., "summary.json")
    local skip_file2="${4:-}" # optional second skip file
    local prefix="$5"        # file prefix: "go" or "py"
    # $6 = assoc-array name for counts, $7 = times, $8 = parts
    # Using eval instead of nameref (local -n) for Bash 4.2 compat.
    local _c_name="$6" _t_name="$7" _p_name="$8"

    for name in "${FIXTURE_NAMES[@]}"; do
        result_dir="$RESULTS_DIR/${impl}/$name"
        mkdir -p "$result_dir"

        # Skip if already completed
        if [[ -f "$result_dir/$skip_file" ]]; then
            continue
        fi
        if [[ -n "$skip_file2" && -f "$result_dir/$skip_file2" ]]; then
            continue
        fi

        # Get tier: "HH:MM:SS SECONDS"
        read -r tier_time tier_secs <<< "$($tier_fn "$name")"
        tier_key="${tier_time//:/}"  # e.g., "01500" for 0:15:00

        # Pick partition
        local part="$PARTITION"
        if (( tier_secs < PREEMPT_CUTOFF )); then
            part="$PREEMPT_PARTITION"
        fi

        tier_file="$ARRAY_DIR/${prefix}_fixtures_${tier_key}.txt"
        local _existing
        eval "_existing=\${${_c_name}[\$tier_key]+x}"
        if [[ -z "$_existing" ]]; then
            > "$tier_file"
            eval "${_c_name}[\$tier_key]=0"
            eval "${_t_name}[\$tier_key]=\$tier_time"
            eval "${_p_name}[\$tier_key]=\$part"
        fi
        echo "$name" >> "$tier_file"
        local _cur
        eval "_cur=\${${_c_name}[\$tier_key]}"
        eval "${_c_name}[\$tier_key]=$(( _cur + 1 ))"
    done
}

# ── Group Go fixtures ───────────────────────────────────────────────────────
group_fixtures "go-cpu" "get_go_tier" "summary.json" "" "go" \
    GO_TIER_COUNTS GO_TIER_TIMES GO_TIER_PARTS

# ── Group Python fixtures ───────────────────────────────────────────────────
if [[ "$GO_ONLY" == false ]]; then
    group_fixtures "original-python" "get_py_tier" "tree_summaries.json.gz" "" "py" \
        PY_TIER_COUNTS PY_TIER_TIMES PY_TIER_PARTS
fi

# ── Summary ─────────────────────────────────────────────────────────────────
echo "--- Array Job Summary ---"
GO_TOTAL=0
for tier in $(echo "${!GO_TIER_COUNTS[@]}" | tr ' ' '\n' | sort); do
    c="${GO_TIER_COUNTS[$tier]}"
    echo "  Go array:      $c fixtures @ ${GO_TIER_TIMES[$tier]} → ${GO_TIER_PARTS[$tier]}"
    GO_TOTAL=$((GO_TOTAL + c))
done
[[ $GO_TOTAL -eq 0 ]] && echo "  Go:            all fixtures already completed"

if [[ "$GO_ONLY" == false ]]; then
    PY_TOTAL=0
    for tier in $(echo "${!PY_TIER_COUNTS[@]}" | tr ' ' '\n' | sort); do
        c="${PY_TIER_COUNTS[$tier]}"
        echo "  Python array:  $c fixtures @ ${PY_TIER_TIMES[$tier]} → ${PY_TIER_PARTS[$tier]}"
        PY_TOTAL=$((PY_TOTAL + c))
    done
    [[ $PY_TOTAL -eq 0 ]] && echo "  Python:        all fixtures already completed"
fi
echo ""

# ── Submit Go array jobs (one per time tier) ────────────────────────────────
declare -a ALL_JIDS

for tier in $(echo "${!GO_TIER_COUNTS[@]}" | tr ' ' '\n' | sort); do
    count="${GO_TIER_COUNTS[$tier]}"
    tier_time="${GO_TIER_TIMES[$tier]}"
    tier_part="${GO_TIER_PARTS[$tier]}"
    tier_file="$ARRAY_DIR/go_fixtures_${tier}.txt"
    max_idx=$(( count - 1 ))

    GO_SCRIPT=$(cat << 'HEREDOC_END'
#!/usr/bin/env bash
set -euo pipefail
source "WORKDIR_PLACEHOLDER/env.sh" 2>/dev/null || true

FIXTURE_LIST="TIER_FILE_PLACEHOLDER"
FIXTURE_NAME=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$FIXTURE_LIST")

FIXTURE_DIR="FIXTURES_DIR_PLACEHOLDER/$FIXTURE_NAME"
RESULT_DIR="RESULTS_DIR_PLACEHOLDER/go-cpu/$FIXTURE_NAME"
mkdir -p "$RESULT_DIR"

if [[ -f "$RESULT_DIR/summary.json" ]]; then
    echo "SKIP: $FIXTURE_NAME already completed"
    exit 0
fi

echo "START: $(date) | impl=go-cpu | fixture=$FIXTURE_NAME | host=$(hostname) | task=$SLURM_ARRAY_TASK_ID"
START=$(date +%s)

"BINARY_PLACEHOLDER" \
    --no-gpu \
    -B BURNIN_PLACEHOLDER \
    -s SAMPLES_PLACEHOLDER \
    -j CHAINS_PLACEHOLDER \
    -I CHAIN_INCLUSION_FACTOR_PLACEHOLDER \
    -O "$RESULT_DIR" \
    "$FIXTURE_DIR/ssm_data.txt" \
    "$FIXTURE_DIR/cnv_data.txt"

END=$(date +%s)
ELAPSED=$((END - START))
echo "{\"fixture\":\"$FIXTURE_NAME\",\"impl\":\"go-cpu\",\"wall_seconds\":$ELAPSED,\"exit_code\":0}" \
    > "$RESULT_DIR/timing.json"
echo "END: $(date) | elapsed ${ELAPSED}s"
HEREDOC_END
)

    GO_SCRIPT="${GO_SCRIPT//WORKDIR_PLACEHOLDER/$WORKDIR}"
    GO_SCRIPT="${GO_SCRIPT//TIER_FILE_PLACEHOLDER/$tier_file}"
    GO_SCRIPT="${GO_SCRIPT//FIXTURES_DIR_PLACEHOLDER/$FIXTURES_DIR}"
    GO_SCRIPT="${GO_SCRIPT//RESULTS_DIR_PLACEHOLDER/$RESULTS_DIR}"
    GO_SCRIPT="${GO_SCRIPT//BINARY_PLACEHOLDER/$BINARY}"
    GO_SCRIPT="${GO_SCRIPT//BURNIN_PLACEHOLDER/$BURNIN}"
    GO_SCRIPT="${GO_SCRIPT//SAMPLES_PLACEHOLDER/$SAMPLES}"
    GO_SCRIPT="${GO_SCRIPT//CHAINS_PLACEHOLDER/$CHAINS}"
    GO_SCRIPT="${GO_SCRIPT//CHAIN_INCLUSION_FACTOR_PLACEHOLDER/$CHAIN_INCLUSION_FACTOR}"

    if [[ "$DRY_RUN" == true ]]; then
        echo "  [dry-run] go-cpu (${tier_time}, ${tier_part}): --array=0-${max_idx} (${count} tasks)"
    else
        GO_JID=$(echo "$GO_SCRIPT" | sbatch --parsable \
            --job-name="simval_go_${tier}" \
            --partition="$tier_part" \
            --cpus-per-task="$CHAINS" \
            --mem="${MEM_MB}M" \
            --time="$tier_time" \
            --array="0-${max_idx}" \
            --output="${LOGS_DIR}/go_${tier}_%a.out" \
            --error="${LOGS_DIR}/go_${tier}_%a.err")
        echo "  [submit] go-cpu (${tier_time}, ${tier_part}): job $GO_JID (${count} tasks)"
        ALL_JIDS+=("$GO_JID")
    fi
done

# ── Submit Python array jobs (one per time tier) ────────────────────────────
if [[ "$GO_ONLY" == false ]]; then
    for tier in $(echo "${!PY_TIER_COUNTS[@]}" | tr ' ' '\n' | sort); do
        count="${PY_TIER_COUNTS[$tier]}"
        tier_time="${PY_TIER_TIMES[$tier]}"
        tier_part="${PY_TIER_PARTS[$tier]}"
        tier_file="$ARRAY_DIR/py_fixtures_${tier}.txt"
        max_idx=$(( count - 1 ))

        PY_SCRIPT=$(cat << 'HEREDOC_END'
#!/usr/bin/env bash
set -euo pipefail
source "WORKDIR_PLACEHOLDER/env.sh" 2>/dev/null || true

FIXTURE_LIST="TIER_FILE_PLACEHOLDER"
FIXTURE_NAME=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$FIXTURE_LIST")

FIXTURE_DIR="FIXTURES_DIR_PLACEHOLDER/$FIXTURE_NAME"
RESULT_DIR="RESULTS_DIR_PLACEHOLDER/original-python/$FIXTURE_NAME"
mkdir -p "$RESULT_DIR"

if [[ -f "$RESULT_DIR/tree_summaries.json.gz" ]]; then
    echo "SKIP: $FIXTURE_NAME already completed"
    exit 0
fi

echo "START: $(date) | impl=original-python | fixture=$FIXTURE_NAME | host=$(hostname) | task=$SLURM_ARRAY_TASK_ID"
START=$(date +%s)

cd "$RESULT_DIR"
singularity exec --bind BIND_PLACEHOLDER "SIF_PLACEHOLDER" \
    python2 PHYLOWGSDIR_PLACEHOLDER/multievolve.py \
    -n PY_CHAINS_PLACEHOLDER \
    -B BURNIN_PLACEHOLDER -s SAMPLES_PLACEHOLDER \
    -O "$RESULT_DIR/chains" \
    --ssms "$FIXTURE_DIR/ssm_data.txt" \
    --cnvs "$FIXTURE_DIR/cnv_data.txt"

# write_results.py on the merged output
if [[ -f "$RESULT_DIR/chains/trees.zip" ]]; then
    singularity exec --bind BIND_PLACEHOLDER "SIF_PLACEHOLDER" \
        python2 PHYLOWGSDIR_PLACEHOLDER/write_results.py \
        "$FIXTURE_NAME" \
        "$RESULT_DIR/chains/trees.zip" \
        "$RESULT_DIR/tree_summaries.json.gz" \
        "$RESULT_DIR/mutlist.json.gz" \
        "$RESULT_DIR/mutass.zip" \
        --include-ssm-names 2>/dev/null || true
fi

END=$(date +%s)
ELAPSED=$((END - START))
echo "{\"fixture\":\"$FIXTURE_NAME\",\"impl\":\"original-python\",\"wall_seconds\":$ELAPSED,\"exit_code\":0}" \
    > "$RESULT_DIR/timing.json"
echo "END: $(date) | elapsed ${ELAPSED}s"
HEREDOC_END
)

        PY_SCRIPT="${PY_SCRIPT//WORKDIR_PLACEHOLDER/$WORKDIR}"
        PY_SCRIPT="${PY_SCRIPT//TIER_FILE_PLACEHOLDER/$tier_file}"
        PY_SCRIPT="${PY_SCRIPT//FIXTURES_DIR_PLACEHOLDER/$FIXTURES_DIR}"
        PY_SCRIPT="${PY_SCRIPT//RESULTS_DIR_PLACEHOLDER/$RESULTS_DIR}"
        PY_SCRIPT="${PY_SCRIPT//SIF_PLACEHOLDER/$PHYLOWGS_SIF}"
        PY_SCRIPT="${PY_SCRIPT//BIND_PLACEHOLDER/$BIND_PATHS}"
        PY_SCRIPT="${PY_SCRIPT//PHYLOWGSDIR_PLACEHOLDER/$PHYLOWGS_DIR}"
        PY_SCRIPT="${PY_SCRIPT//BURNIN_PLACEHOLDER/$BURNIN}"
        PY_SCRIPT="${PY_SCRIPT//SAMPLES_PLACEHOLDER/$SAMPLES}"
        PY_SCRIPT="${PY_SCRIPT//PY_CHAINS_PLACEHOLDER/$PY_CHAINS}"

        if [[ "$DRY_RUN" == true ]]; then
            echo "  [dry-run] python (${tier_time}, ${tier_part}): --array=0-${max_idx} (${count} tasks)"
        else
            PY_JID=$(echo "$PY_SCRIPT" | sbatch --parsable \
                --job-name="simval_py_${tier}" \
                --partition="$tier_part" \
                --cpus-per-task="$PY_CHAINS" \
                --mem="${MEM_MB}M" \
                --time="$tier_time" \
                --array="0-${max_idx}" \
                --output="${LOGS_DIR}/py_${tier}_%a.out" \
                --error="${LOGS_DIR}/py_${tier}_%a.err")
            echo "  [submit] python (${tier_time}, ${tier_part}): job $PY_JID (${count} tasks)"
            ALL_JIDS+=("$PY_JID")
        fi
    done
fi

# ── Submit scoring job (depends on all array jobs) ──────────────────────────
if [[ ${#ALL_JIDS[@]} -gt 0 ]]; then
    ALL_DEP=$(IFS=':'; echo "${ALL_JIDS[*]}")

    SCORE_SCRIPT=$(cat << EOF
#!/usr/bin/env bash
#SBATCH --job-name=simval_score
#SBATCH --partition=${PARTITION}
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0:30:00
#SBATCH --output=${LOGS_DIR}/scoring.out
#SBATCH --error=${LOGS_DIR}/scoring.err

set -euo pipefail
source "$WORKDIR/env.sh" 2>/dev/null || true

echo "START scoring: \$(date)"

python3 "$WORKDIR/score_results.py" \
    --all \
    --fixture-base "$FIXTURES_DIR" \
    --result-base "$RESULTS_DIR/go-cpu" \
    --outdir "$WORKDIR/analysis/go-cpu"

if [[ -d "$RESULTS_DIR/original-python" ]]; then
    python3 "$WORKDIR/score_results.py" \
        --all \
        --fixture-base "$FIXTURES_DIR" \
        --result-base "$RESULTS_DIR/original-python" \
        --outdir "$WORKDIR/analysis/original-python"
fi

echo "START plotting: \$(date)"

python3 "$WORKDIR/plot_results.py" \
    --analysis-dir "$WORKDIR/analysis/go-cpu" \
    --result-base "$RESULTS_DIR/go-cpu" \
    --fixture-base "$FIXTURES_DIR" \
    --outdir "$WORKDIR/analysis/go-cpu/plots"

if [[ -d "$WORKDIR/analysis/original-python" ]]; then
    python3 "$WORKDIR/plot_results.py" \
        --analysis-dir "$WORKDIR/analysis/original-python" \
        --result-base "$RESULTS_DIR/original-python" \
        --fixture-base "$FIXTURES_DIR" \
        --outdir "$WORKDIR/analysis/original-python/plots"
fi

if [[ -f "$WORKDIR/analysis/go-cpu/scores.json" && -f "$WORKDIR/analysis/original-python/scores.json" ]]; then
    python3 "$WORKDIR/compare_implementations.py" \
        --go-scores "$WORKDIR/analysis/go-cpu/scores.json" \
        --py-scores "$WORKDIR/analysis/original-python/scores.json" \
        --outdir "$WORKDIR/analysis/comparison"
fi

echo "END scoring + plotting: \$(date)"
EOF
)

    if [[ "$DRY_RUN" == true ]]; then
        echo ""
        echo "  [dry-run] scoring job (depends on ${#ALL_JIDS[@]} array jobs)"
    else
        SCORE_JID=$(echo "$SCORE_SCRIPT" | sbatch --parsable --dependency="afterany:${ALL_DEP}")
        echo ""
        echo "=== All jobs submitted ==="
        echo "  Array jobs:  ${#ALL_JIDS[@]}"
        echo "  Score job:   $SCORE_JID (runs after all arrays complete)"
        echo ""
        echo "Monitor:  squeue -u \$USER | grep simval"
        echo "Results:  $WORKDIR/analysis/"
        echo "Plots:    $WORKDIR/analysis/go-cpu/plots/"
        [[ "$GO_ONLY" == false ]] && echo "Compare:  $WORKDIR/analysis/comparison/"
    fi
elif [[ "$DRY_RUN" == true ]]; then
    echo ""
    echo "  [dry-run] All fixtures already completed — scoring only"
else
    echo ""
    echo "All fixtures already completed. Running scoring + plotting directly:"
    source "$WORKDIR/env.sh" 2>/dev/null || true
    python3 "$WORKDIR/score_results.py" \
        --all \
        --fixture-base "$FIXTURES_DIR" \
        --result-base "$RESULTS_DIR/go-cpu" \
        --outdir "$WORKDIR/analysis/go-cpu"
    python3 "$WORKDIR/plot_results.py" \
        --analysis-dir "$WORKDIR/analysis/go-cpu" \
        --result-base "$RESULTS_DIR/go-cpu" \
        --fixture-base "$FIXTURES_DIR" \
        --outdir "$WORKDIR/analysis/go-cpu/plots"
fi
