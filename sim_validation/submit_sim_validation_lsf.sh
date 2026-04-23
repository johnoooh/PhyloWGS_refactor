#!/usr/bin/env bash
# submit_sim_validation_lsf.sh — Submit LSF array jobs to run Go port AND original
# Python on simulated fixtures, then score and compare.
#
# LSF port of submit_sim_validation.sh (Slurm version).
#
# Runs two implementations per fixture:
#   - go-cpu: Go port (from setup_validation.sh binary)
#   - original-python: Original PhyloWGS via singularity container
#
# Jobs are submitted as LSF job arrays grouped by M-tier:
#   - Each tier gets its own wall time (derived from P95 × 2.33 × 1.05 safety)
#   - All jobs go to the "general" queue by default
#
# Usage:
#   bash submit_sim_validation_lsf.sh --workdir /path/to/workdir --sif /path/to/phylowgs.sif [options]
#
# Options:
#   --workdir DIR       Workdir from setup_validation.sh
#   --sif PATH          PhyloWGS singularity image (.sif) for original Python
#   --burnin N          MCMC burn-in (default: 1000)
#   --samples N         MCMC samples (default: 2500)
#   --chains N          Parallel chains for Go (default: 8)
#   --mem GB            Memory in GB (default: 8)
#   --queue NAME        LSF queue (default: general)
#   --go-only           Only run Go port (skip original Python)
#   --dry-run           Print jobs without submitting

set -euo pipefail

WORKDIR="${SIM_VALIDATION_WORKDIR:-$(pwd)/sim_validation_workdir}"
BURNIN=1000
SAMPLES=2500
CHAINS=8
MEM_GB=8
QUEUE="general"
PY_CHAINS=8
DRY_RUN=false
GO_ONLY=false
PHYLOWGS_SIF=""
BIND_PATHS=""
PHYLOWGS_DIR="/usr/bin/phylowgs"

# ── Wall-time tiers based on mutation count (M) ─────────────────────────────
# Derived from P95 of observed runtimes at B=500/s=1000, scaled by 2.33×
# (for new B=1000/s=2500) with 1.05× safety margin, rounded up to 15 min.
# Then doubled for LSF server (slower I/O than Slurm cluster).
# Python tiers include +20 min headroom for write_results.py post-processing.
#
# Returns: "HH:MM SECONDS" (two space-separated tokens)
# Note: LSF -W uses HH:MM format (no seconds).
get_go_tier() {
    local m_val
    if [[ "$1" =~ _M([0-9]+)_ ]]; then m_val="${BASH_REMATCH[1]}"; else echo "16:30 59400"; return; fi
    if   (( m_val <= 30 ));  then echo "1:00 3600"
    elif (( m_val <= 50 ));  then echo "1:30 5400"
    elif (( m_val <= 100 )); then echo "2:30 9000"
    elif (( m_val <= 150 )); then echo "9:00 32400"
    elif (( m_val <= 250 )); then echo "16:30 59400"
    else                          echo "15:30 55800"
    fi
}

get_py_tier() {
    # Python now runs multievolve.py with 8 chains. LSF tiers are 2× Slurm
    # (slower I/O), and 8× single-chain for multievolve overhead.
    local m_val
    if [[ "$1" =~ _M([0-9]+)_ ]]; then m_val="${BASH_REMATCH[1]}"; else echo "144:00 518400"; return; fi
    if   (( m_val <= 30 ));  then echo "20:00 72000"
    elif (( m_val <= 50 ));  then echo "44:00 158400"
    elif (( m_val <= 100 )); then echo "88:00 316800"
    elif (( m_val <= 150 )); then echo "100:00 360000"
    elif (( m_val <= 250 )); then echo "112:00 403200"
    else                          echo "144:00 518400"
    fi
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --workdir)    WORKDIR="$2";       shift 2 ;;
        --sif)        PHYLOWGS_SIF="$2";  shift 2 ;;
        --burnin)     BURNIN="$2";        shift 2 ;;
        --samples)    SAMPLES="$2";       shift 2 ;;
        --chains)     CHAINS="$2";        shift 2 ;;
        --mem)        MEM_GB="$2";        shift 2 ;;
        --queue)          QUEUE="$2";         shift 2 ;;
        --bind-paths)     BIND_PATHS="$2";    shift 2 ;;
        --phylowgs-dir)   PHYLOWGS_DIR="$2";  shift 2 ;;
        --go-only)        GO_ONLY=true;       shift ;;
        --dry-run)        DRY_RUN=true;       shift ;;
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
    _top_mount=$(echo "$WORKDIR" | cut -d/ -f1-2)
    BIND_PATHS="$_top_mount"
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

echo "=== Simulation Validation Job Submission — LSF (Array Jobs) ==="
echo "Workdir:         $WORKDIR"
echo "Fixtures:        $N_FIXTURES"
echo "Implementations: $IMPLS"
echo "Config:          B=$BURNIN s=$SAMPLES j=$CHAINS (Go); B=$BURNIN s=$SAMPLES (Python)"
echo "Queue:           $QUEUE"
echo "Memory:          ${MEM_GB}GB"
echo ""

# ── Helper: group fixtures by time tier for an impl ─────────────────────────
declare -A GO_TIER_COUNTS GO_TIER_TIMES
declare -A PY_TIER_COUNTS PY_TIER_TIMES

group_fixtures() {
    local impl="$1"          # "go-cpu" or "original-python"
    local tier_fn="$2"       # "get_go_tier" or "get_py_tier"
    local skip_file="$3"     # file to check for skip (e.g., "summary.json")
    local skip_file2="${4:-}" # optional second skip file
    local prefix="$5"        # file prefix: "go" or "py"
    # $6 = assoc-array name for counts, $7 = times
    # Using eval instead of nameref (local -n) for Bash 4.2 compat.
    local _c_name="$6" _t_name="$7"

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

        # Get tier: "HH:MM SECONDS"
        read -r tier_time tier_secs <<< "$($tier_fn "$name")"
        tier_key="${tier_time//:/}"  # e.g., "0130" for 1:30

        tier_file="$ARRAY_DIR/${prefix}_fixtures_${tier_key}.txt"
        local _existing
        eval "_existing=\${${_c_name}[\$tier_key]+x}"
        if [[ -z "$_existing" ]]; then
            > "$tier_file"
            eval "${_c_name}[\$tier_key]=0"
            eval "${_t_name}[\$tier_key]=\$tier_time"
        fi
        echo "$name" >> "$tier_file"
        local _cur
        eval "_cur=\${${_c_name}[\$tier_key]}"
        eval "${_c_name}[\$tier_key]=$(( _cur + 1 ))"
    done
}

# ── Group Go fixtures ───────────────────────────────────────────────────────
group_fixtures "go-cpu" "get_go_tier" "summary.json" "" "go" \
    GO_TIER_COUNTS GO_TIER_TIMES

# ── Group Python fixtures ───────────────────────────────────────────────────
if [[ "$GO_ONLY" == false ]]; then
    group_fixtures "original-python" "get_py_tier" "tree_summaries.json.gz" "" "py" \
        PY_TIER_COUNTS PY_TIER_TIMES
fi

# ── Summary ─────────────────────────────────────────────────────────────────
echo "--- Array Job Summary ---"
GO_TOTAL=0
for tier in $(echo "${!GO_TIER_COUNTS[@]}" | tr ' ' '\n' | sort); do
    c="${GO_TIER_COUNTS[$tier]}"
    echo "  Go array:      $c fixtures @ ${GO_TIER_TIMES[$tier]}"
    GO_TOTAL=$((GO_TOTAL + c))
done
[[ $GO_TOTAL -eq 0 ]] && echo "  Go:            all fixtures already completed"

if [[ "$GO_ONLY" == false ]]; then
    PY_TOTAL=0
    for tier in $(echo "${!PY_TIER_COUNTS[@]}" | tr ' ' '\n' | sort); do
        c="${PY_TIER_COUNTS[$tier]}"
        echo "  Python array:  $c fixtures @ ${PY_TIER_TIMES[$tier]}"
        PY_TOTAL=$((PY_TOTAL + c))
    done
    [[ $PY_TOTAL -eq 0 ]] && echo "  Python:        all fixtures already completed"
fi
echo ""

# ── Submit Go array jobs (one per time tier) ────────────────────────────────
declare -a ALL_JOB_NAMES

for tier in $(echo "${!GO_TIER_COUNTS[@]}" | tr ' ' '\n' | sort); do
    count="${GO_TIER_COUNTS[$tier]}"
    tier_time="${GO_TIER_TIMES[$tier]}"
    tier_file="$ARRAY_DIR/go_fixtures_${tier}.txt"
    job_name="simval_go_${tier}"

    # Write the job script to a file (LSF bsub reads from file or stdin)
    JOB_SCRIPT="$ARRAY_DIR/go_job_${tier}.sh"
    cat > "$JOB_SCRIPT" << HEREDOC_END
#!/usr/bin/env bash
set -euo pipefail
source "$WORKDIR/env.sh" 2>/dev/null || true

# LSF array index is 1-based (\$LSB_JOBINDEX)
FIXTURE_LIST="$tier_file"
FIXTURE_NAME=\$(sed -n "\${LSB_JOBINDEX}p" "\$FIXTURE_LIST")

FIXTURE_DIR="$FIXTURES_DIR/\$FIXTURE_NAME"
RESULT_DIR="$RESULTS_DIR/go-cpu/\$FIXTURE_NAME"
mkdir -p "\$RESULT_DIR"

if [[ -f "\$RESULT_DIR/summary.json" ]]; then
    echo "SKIP: \$FIXTURE_NAME already completed"
    exit 0
fi

echo "START: \$(date) | impl=go-cpu | fixture=\$FIXTURE_NAME | host=\$(hostname) | index=\$LSB_JOBINDEX"
START=\$(date +%s)

"$BINARY" \\
    --no-gpu \\
    -B $BURNIN \\
    -s $SAMPLES \\
    -j $CHAINS \\
    -O "\$RESULT_DIR" \\
    "\$FIXTURE_DIR/ssm_data.txt" \\
    "\$FIXTURE_DIR/cnv_data.txt"

END=\$(date +%s)
ELAPSED=\$((END - START))
echo "{\"fixture\":\"\$FIXTURE_NAME\",\"impl\":\"go-cpu\",\"wall_seconds\":\$ELAPSED,\"exit_code\":0}" \\
    > "\$RESULT_DIR/timing.json"
echo "END: \$(date) | elapsed \${ELAPSED}s"
HEREDOC_END
    chmod +x "$JOB_SCRIPT"

    if [[ "$DRY_RUN" == true ]]; then
        echo "  [dry-run] go-cpu (${tier_time}): -J '${job_name}[1-${count}]' (${count} tasks)"
    else
        # LSF arrays are 1-based: [1-N]
        BSUB_OUT=$(bsub \
            -J "${job_name}[1-${count}]" \
            -q "$QUEUE" \
            -n "$CHAINS" \
            -R "rusage[mem=${MEM_GB}]" \
            -W "$tier_time" \
            -o "${LOGS_DIR}/go_${tier}_%I.out" \
            -e "${LOGS_DIR}/go_${tier}_%I.err" \
            < "$JOB_SCRIPT" 2>&1)
        echo "  [submit] go-cpu (${tier_time}): ${BSUB_OUT} (${count} tasks)"
        ALL_JOB_NAMES+=("$job_name")
    fi
done

# ── Submit Python array jobs (one per time tier) ────────────────────────────
if [[ "$GO_ONLY" == false ]]; then
    for tier in $(echo "${!PY_TIER_COUNTS[@]}" | tr ' ' '\n' | sort); do
        count="${PY_TIER_COUNTS[$tier]}"
        tier_time="${PY_TIER_TIMES[$tier]}"
        tier_file="$ARRAY_DIR/py_fixtures_${tier}.txt"
        job_name="simval_py_${tier}"

        JOB_SCRIPT="$ARRAY_DIR/py_job_${tier}.sh"
        cat > "$JOB_SCRIPT" << HEREDOC_END
#!/usr/bin/env bash
set -euo pipefail
source "$WORKDIR/env.sh" 2>/dev/null || true

FIXTURE_LIST="$tier_file"
FIXTURE_NAME=\$(sed -n "\${LSB_JOBINDEX}p" "\$FIXTURE_LIST")

FIXTURE_DIR="$FIXTURES_DIR/\$FIXTURE_NAME"
RESULT_DIR="$RESULTS_DIR/original-python/\$FIXTURE_NAME"
mkdir -p "\$RESULT_DIR"

if [[ -f "\$RESULT_DIR/tree_summaries.json.gz" ]]; then
    echo "SKIP: \$FIXTURE_NAME already completed"
    exit 0
fi

echo "START: \$(date) | impl=original-python | fixture=\$FIXTURE_NAME | host=\$(hostname) | index=\$LSB_JOBINDEX"
START=\$(date +%s)

cd "\$RESULT_DIR"
singularity exec --bind $BIND_PATHS "$PHYLOWGS_SIF" \\
    python2 $PHYLOWGS_DIR/multievolve.py \\
    -n $PY_CHAINS \\
    -B $BURNIN -s $SAMPLES \\
    -O "\$RESULT_DIR/chains" \\
    --ssms "\$FIXTURE_DIR/ssm_data.txt" \\
    --cnvs "\$FIXTURE_DIR/cnv_data.txt"

if [[ -f "\$RESULT_DIR/chains/trees.zip" ]]; then
    singularity exec --bind $BIND_PATHS "$PHYLOWGS_SIF" \\
        python2 $PHYLOWGS_DIR/write_results.py \\
        "\$FIXTURE_NAME" \\
        "\$RESULT_DIR/chains/trees.zip" \\
        "\$RESULT_DIR/tree_summaries.json.gz" \\
        "\$RESULT_DIR/mutlist.json.gz" \\
        "\$RESULT_DIR/mutass.zip" \\
        --include-ssm-names 2>/dev/null || true
fi

END=\$(date +%s)
ELAPSED=\$((END - START))
echo "{\"fixture\":\"\$FIXTURE_NAME\",\"impl\":\"original-python\",\"wall_seconds\":\$ELAPSED,\"exit_code\":0}" \\
    > "\$RESULT_DIR/timing.json"
echo "END: \$(date) | elapsed \${ELAPSED}s"
HEREDOC_END
        chmod +x "$JOB_SCRIPT"

        if [[ "$DRY_RUN" == true ]]; then
            echo "  [dry-run] python (${tier_time}): -J '${job_name}[1-${count}]' (${count} tasks)"
        else
            BSUB_OUT=$(bsub \
                -J "${job_name}[1-${count}]" \
                -q "$QUEUE" \
                -n "$PY_CHAINS" \
                -R "rusage[mem=${MEM_GB}]" \
                -W "$tier_time" \
                -o "${LOGS_DIR}/py_${tier}_%I.out" \
                -e "${LOGS_DIR}/py_${tier}_%I.err" \
                < "$JOB_SCRIPT" 2>&1)
            echo "  [submit] python (${tier_time}): ${BSUB_OUT} (${count} tasks)"
            ALL_JOB_NAMES+=("$job_name")
        fi
    done
fi

# ── Submit scoring job (depends on all array jobs) ──────────────────────────
if [[ ${#ALL_JOB_NAMES[@]} -gt 0 ]]; then
    # Build LSF dependency expression: done(name1) && done(name2) && ...
    DEP_EXPR=""
    for jn in "${ALL_JOB_NAMES[@]}"; do
        if [[ -z "$DEP_EXPR" ]]; then
            DEP_EXPR="done(${jn})"
        else
            DEP_EXPR="${DEP_EXPR} && done(${jn})"
        fi
    done

    SCORE_SCRIPT="$ARRAY_DIR/score_job.sh"
    cat > "$SCORE_SCRIPT" << EOF
#!/usr/bin/env bash
set -euo pipefail
source "$WORKDIR/env.sh" 2>/dev/null || true

echo "START scoring: \$(date)"

python3 "$WORKDIR/score_results.py" \\
    --all \\
    --fixture-base "$FIXTURES_DIR" \\
    --result-base "$RESULTS_DIR/go-cpu" \\
    --outdir "$WORKDIR/analysis/go-cpu"

if [[ -d "$RESULTS_DIR/original-python" ]]; then
    python3 "$WORKDIR/score_results.py" \\
        --all \\
        --fixture-base "$FIXTURES_DIR" \\
        --result-base "$RESULTS_DIR/original-python" \\
        --outdir "$WORKDIR/analysis/original-python"
fi

echo "START plotting: \$(date)"

python3 "$WORKDIR/plot_results.py" \\
    --analysis-dir "$WORKDIR/analysis/go-cpu" \\
    --result-base "$RESULTS_DIR/go-cpu" \\
    --fixture-base "$FIXTURES_DIR" \\
    --outdir "$WORKDIR/analysis/go-cpu/plots"

if [[ -d "$WORKDIR/analysis/original-python" ]]; then
    python3 "$WORKDIR/plot_results.py" \\
        --analysis-dir "$WORKDIR/analysis/original-python" \\
        --result-base "$RESULTS_DIR/original-python" \\
        --fixture-base "$FIXTURES_DIR" \\
        --outdir "$WORKDIR/analysis/original-python/plots"
fi

if [[ -f "$WORKDIR/analysis/go-cpu/scores.json" && -f "$WORKDIR/analysis/original-python/scores.json" ]]; then
    python3 "$WORKDIR/compare_implementations.py" \\
        --go-scores "$WORKDIR/analysis/go-cpu/scores.json" \\
        --py-scores "$WORKDIR/analysis/original-python/scores.json" \\
        --outdir "$WORKDIR/analysis/comparison"
fi

echo "END scoring + plotting: \$(date)"
EOF
    chmod +x "$SCORE_SCRIPT"

    if [[ "$DRY_RUN" == true ]]; then
        echo ""
        echo "  [dry-run] scoring job (depends on: $DEP_EXPR)"
    else
        BSUB_OUT=$(bsub \
            -J "simval_score" \
            -q "$QUEUE" \
            -n 1 \
            -R "rusage[mem=4]" \
            -W "0:30" \
            -w "$DEP_EXPR" \
            -o "${LOGS_DIR}/scoring.out" \
            -e "${LOGS_DIR}/scoring.err" \
            < "$SCORE_SCRIPT" 2>&1)
        echo ""
        echo "=== All jobs submitted ==="
        echo "  Array jobs:  ${#ALL_JOB_NAMES[@]}"
        echo "  Score job:   ${BSUB_OUT} (runs after all arrays complete)"
        echo ""
        echo "Monitor:  bjobs -w | grep simval"
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
