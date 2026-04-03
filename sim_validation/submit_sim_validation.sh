#!/usr/bin/env bash
# submit_sim_validation.sh — Submit SLURM jobs to run Go port AND original Python
# on simulated fixtures, then score and compare.
#
# Runs two implementations per fixture:
#   - go-cpu: Go port (from setup_validation.sh binary)
#   - original-python: Original PhyloWGS via singularity container
#
# Usage:
#   bash submit_sim_validation.sh --workdir /path/to/workdir --sif /path/to/phylowgs.sif [options]
#
# Options:
#   --workdir DIR       Workdir from setup_validation.sh
#   --sif PATH          PhyloWGS singularity image (.sif) for original Python
#   --burnin N          MCMC burn-in (default: 500)
#   --samples N         MCMC samples (default: 1000)
#   --chains N          Parallel chains for Go (default: 4)
#   --time HH:MM:SS     Wall time per job (default: 2:00:00)
#   --mem MB            Memory in MB (default: 8000)
#   --partition NAME    SLURM partition (default: cmobic_cpu)
#   --go-only           Only run Go port (skip original Python)
#   --dry-run           Print jobs without submitting

set -euo pipefail

WORKDIR="${SIM_VALIDATION_WORKDIR:-$(pwd)/sim_validation_workdir}"
BURNIN=500
SAMPLES=1000
CHAINS=4
TIME_LIMIT="2:00:00"
MEM_MB=8000
PARTITION="cmobic_cpu"
DRY_RUN=false
GO_ONLY=false
PHYLOWGS_SIF=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --workdir)    WORKDIR="$2";       shift 2 ;;
        --sif)        PHYLOWGS_SIF="$2";  shift 2 ;;
        --burnin)     BURNIN="$2";        shift 2 ;;
        --samples)    SAMPLES="$2";       shift 2 ;;
        --chains)     CHAINS="$2";        shift 2 ;;
        --time)       TIME_LIMIT="$2";    shift 2 ;;
        --mem)        MEM_MB="$2";        shift 2 ;;
        --partition)  PARTITION="$2";     shift 2 ;;
        --go-only)    GO_ONLY=true;       shift ;;
        --dry-run)    DRY_RUN=true;       shift ;;
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
fi

FIXTURES_DIR="$WORKDIR/fixtures"
RESULTS_DIR="$WORKDIR/results"
LOGS_DIR="$WORKDIR/logs"
mkdir -p "$RESULTS_DIR" "$LOGS_DIR"

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

echo "=== Simulation Validation Job Submission ==="
echo "Workdir:         $WORKDIR"
echo "Fixtures:        $N_FIXTURES"
echo "Implementations: $IMPLS"
echo "Config:          B=$BURNIN s=$SAMPLES j=$CHAINS (Go); B=$BURNIN s=$SAMPLES (Python)"
echo "Partition:       $PARTITION"
echo "Time:            $TIME_LIMIT"
echo "Memory:          ${MEM_MB}MB"
echo ""

# ── Submit jobs ──────────────────────────────────────────────────────────────
declare -a RUN_JIDS

for name in "${FIXTURE_NAMES[@]}"; do
    fixture_dir="$FIXTURES_DIR/$name"

    for impl in $IMPLS; do
        result_dir="$RESULTS_DIR/${impl}/$name"
        mkdir -p "$result_dir"

        # Skip if already completed
        if [[ -f "$result_dir/summary.json" ]] || [[ -f "$result_dir/mcmc_samples.txt" ]]; then
            echo "  [skip] $impl × $name (already completed)"
            continue
        fi

        SCRIPT=""
        case "$impl" in

          go-cpu)
            SCRIPT=$(cat << EOF
#!/usr/bin/env bash
#SBATCH --job-name=simval_go_${name}
#SBATCH --partition=${PARTITION}
#SBATCH --cpus-per-task=${CHAINS}
#SBATCH --mem=${MEM_MB}M
#SBATCH --time=${TIME_LIMIT}
#SBATCH --output=${LOGS_DIR}/go-cpu_${name}.out
#SBATCH --error=${LOGS_DIR}/go-cpu_${name}.err

set -euo pipefail
source "$WORKDIR/env.sh" 2>/dev/null || true

echo "START: \$(date) | impl=go-cpu | fixture=$name | host=\$(hostname)"
START=\$(date +%s)

"$BINARY" \
    --no-gpu \
    -B $BURNIN \
    -s $SAMPLES \
    -j $CHAINS \
    -O "$result_dir" \
    "$fixture_dir/ssm_data.txt" \
    "$fixture_dir/cnv_data.txt"

END=\$(date +%s)
ELAPSED=\$((END - START))
echo "{\"fixture\":\"$name\",\"impl\":\"go-cpu\",\"wall_seconds\":\$ELAPSED,\"exit_code\":0}" \
    > "$result_dir/timing.json"
echo "END: \$(date) | elapsed \${ELAPSED}s"
EOF
)
            ;;

          original-python)
            SCRIPT=$(cat << EOF
#!/usr/bin/env bash
#SBATCH --job-name=simval_py_${name}
#SBATCH --partition=${PARTITION}
#SBATCH --cpus-per-task=1
#SBATCH --mem=${MEM_MB}M
#SBATCH --time=${TIME_LIMIT}
#SBATCH --output=${LOGS_DIR}/original-python_${name}.out
#SBATCH --error=${LOGS_DIR}/original-python_${name}.err

set -euo pipefail
source "$WORKDIR/env.sh" 2>/dev/null || true

echo "START: \$(date) | impl=original-python | fixture=$name | host=\$(hostname)"
START=\$(date +%s)

# Run evolve.py inside the container (Python 2, original morrislab code)
cd "$result_dir"
singularity exec --bind /data1,$WORKDIR "$PHYLOWGS_SIF" \
    python2 /usr/bin/phylowgs/evolve.py \
    -B $BURNIN -s $SAMPLES \
    "$fixture_dir/ssm_data.txt" \
    "$fixture_dir/cnv_data.txt"

# Run write_results.py to generate JSON output (tree summaries + mutation assignments)
if [[ -f "$result_dir/trees.zip" ]]; then
    singularity exec --bind /data1,$WORKDIR "$PHYLOWGS_SIF" \
        python2 /usr/bin/phylowgs/write_results.py \
        "$name" \
        "$result_dir/trees.zip" \
        "$result_dir/tree_summaries.json.gz" \
        "$result_dir/mutlist.json.gz" \
        "$result_dir/mutass.zip" \
        --include-ssm-names 2>/dev/null || true
fi

END=\$(date +%s)
ELAPSED=\$((END - START))
echo "{\"fixture\":\"$name\",\"impl\":\"original-python\",\"wall_seconds\":\$ELAPSED,\"exit_code\":0}" \
    > "$result_dir/timing.json"
echo "END: \$(date) | elapsed \${ELAPSED}s"
EOF
)
            ;;
        esac

        if [[ "$DRY_RUN" == true ]]; then
            echo "  [dry-run] $impl × $name"
        else
            JID=$(echo "$SCRIPT" | sbatch --parsable)
            echo "  [submit] $impl × $name → job $JID"
            RUN_JIDS+=("$JID")
        fi
    done
done

# ── Submit scoring job (depends on all run jobs) ─────────────────────────────
if [[ ${#RUN_JIDS[@]} -gt 0 ]]; then
    ALL_DEP=$(IFS=':'; echo "${RUN_JIDS[*]}")

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

# Score Go port results
python3 "$WORKDIR/score_results.py" \
    --all \
    --fixture-base "$FIXTURES_DIR" \
    --result-base "$RESULTS_DIR/go-cpu" \
    --outdir "$WORKDIR/analysis/go-cpu"

# Score original Python results (if present)
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

# Compare implementations
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
        echo "  [dry-run] scoring job (depends on ${#RUN_JIDS[@]} run jobs)"
    else
        SCORE_JID=$(echo "$SCORE_SCRIPT" | sbatch --parsable --dependency="afterany:${ALL_DEP}")
        echo ""
        echo "=== All jobs submitted ==="
        echo "  Run jobs:    ${#RUN_JIDS[@]}"
        echo "  Score job:   $SCORE_JID (runs after all complete, including failures)"
        echo ""
        echo "Monitor:  squeue -u \$USER -n simval_%"
        echo "Results:  $WORKDIR/analysis/"
        echo "Plots:    $WORKDIR/analysis/go-cpu/plots/"
        [[ "$GO_ONLY" == false ]] && echo "Compare:  $WORKDIR/analysis/comparison/"
    fi
elif [[ "$DRY_RUN" == true ]]; then
    echo ""
    echo "  [dry-run] $N_FIXTURES fixtures × $(echo $IMPLS | wc -w) impls would be submitted"
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
