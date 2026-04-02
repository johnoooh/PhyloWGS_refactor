#!/usr/bin/env bash
# submit_sim_validation.sh — Submit SLURM jobs to run Go port on simulated fixtures
#
# Reads fixtures from setup_validation.sh output, submits one job per fixture,
# then a final scoring/analysis job.
#
# Usage:
#   bash submit_sim_validation.sh --workdir /path/to/sim_validation_workdir [options]
#
# Options:
#   --workdir DIR       Workdir from setup_validation.sh
#   --burnin N          MCMC burn-in (default: 500)
#   --samples N         MCMC samples (default: 1000)
#   --chains N          Parallel chains (default: 4)
#   --time HH:MM:SS     Wall time per job (default: 2:00:00)
#   --mem MB            Memory in MB (default: 4000)
#   --partition NAME    SLURM partition (default: cmobic_cpu)
#   --dry-run           Print jobs without submitting

set -euo pipefail

WORKDIR="${SIM_VALIDATION_WORKDIR:-$(pwd)/sim_validation_workdir}"
BURNIN=500
SAMPLES=1000
CHAINS=4
TIME_LIMIT="2:00:00"
MEM_MB=4000
PARTITION="cmobic_cpu"
DRY_RUN=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --workdir)    WORKDIR="$2";    shift 2 ;;
        --burnin)     BURNIN="$2";     shift 2 ;;
        --samples)    SAMPLES="$2";    shift 2 ;;
        --chains)     CHAINS="$2";     shift 2 ;;
        --time)       TIME_LIMIT="$2"; shift 2 ;;
        --mem)        MEM_MB="$2";     shift 2 ;;
        --partition)  PARTITION="$2";  shift 2 ;;
        --dry-run)    DRY_RUN=true;    shift ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

# ── Validate workdir ─────────────────────────────────────────────────────────
[[ ! -d "$WORKDIR" ]]           && { echo "Workdir not found: $WORKDIR — run setup_validation.sh first"; exit 1; }
[[ ! -f "$WORKDIR/env.sh" ]]    && { echo "env.sh not found — run setup_validation.sh first"; exit 1; }
[[ ! -d "$WORKDIR/fixtures" ]]  && { echo "No fixtures — run setup_validation.sh first"; exit 1; }

BINARY="$WORKDIR/bin/phylowgs-go"
[[ ! -x "$BINARY" ]] && { echo "Go binary not found: $BINARY"; exit 1; }

FIXTURES_DIR="$WORKDIR/fixtures"
RESULTS_DIR="$WORKDIR/results"
LOGS_DIR="$WORKDIR/logs"
mkdir -p "$RESULTS_DIR" "$LOGS_DIR"

# ── Read manifest ────────────────────────────────────────────────────────────
MANIFEST="$FIXTURES_DIR/manifest.json"
[[ ! -f "$MANIFEST" ]] && { echo "manifest.json not found in fixtures/"; exit 1; }

# Extract fixture names from manifest
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

echo "=== Simulation Validation Job Submission ==="
echo "Workdir:   $WORKDIR"
echo "Fixtures:  $N_FIXTURES"
echo "Config:    B=$BURNIN s=$SAMPLES j=$CHAINS"
echo "Partition: $PARTITION"
echo "Time:      $TIME_LIMIT"
echo ""

# ── Submit one job per fixture ───────────────────────────────────────────────
declare -a RUN_JIDS

for name in "${FIXTURE_NAMES[@]}"; do
    fixture_dir="$FIXTURES_DIR/$name"
    result_dir="$RESULTS_DIR/$name"
    mkdir -p "$result_dir"

    # Skip if already completed
    if [[ -f "$result_dir/summary.json" ]]; then
        echo "  [skip] $name (already completed)"
        continue
    fi

    SCRIPT=$(cat << EOF
#!/usr/bin/env bash
#SBATCH --job-name=simval_${name}
#SBATCH --partition=${PARTITION}
#SBATCH --cpus-per-task=${CHAINS}
#SBATCH --mem=${MEM_MB}M
#SBATCH --time=${TIME_LIMIT}
#SBATCH --output=${LOGS_DIR}/${name}.out
#SBATCH --error=${LOGS_DIR}/${name}.err

set -euo pipefail
source "$WORKDIR/env.sh" 2>/dev/null || true

echo "START: \$(date) | fixture=$name | host=\$(hostname)"
echo "Binary: $BINARY"
echo "Config: B=$BURNIN s=$SAMPLES j=$CHAINS"
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

echo "{\"fixture\":\"$name\",\"wall_seconds\":\$ELAPSED,\"exit_code\":0}" \
    > "$result_dir/timing.json"

echo "END: \$(date) | elapsed \${ELAPSED}s"
EOF
)

    if [[ "$DRY_RUN" == true ]]; then
        echo "  [dry-run] $name"
    else
        JID=$(echo "$SCRIPT" | sbatch --parsable)
        echo "  [submit] $name → job $JID"
        RUN_JIDS+=("$JID")
    fi
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

python3 "$WORKDIR/score_results.py" \
    --all \
    --fixture-base "$FIXTURES_DIR" \
    --result-base "$RESULTS_DIR" \
    --outdir "$WORKDIR/analysis"

echo "START plotting: \$(date)"

python3 "$WORKDIR/plot_results.py" \
    --analysis-dir "$WORKDIR/analysis" \
    --result-base "$RESULTS_DIR" \
    --fixture-base "$FIXTURES_DIR"

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
        echo "Results:  $WORKDIR/analysis/scores.tsv"
        echo "Plots:    $WORKDIR/analysis/plots/"
    fi
elif [[ "$DRY_RUN" == true ]]; then
    echo ""
    echo "  [dry-run] $N_FIXTURES fixtures would be submitted"
else
    echo ""
    echo "All fixtures already completed. Running scoring + plotting directly:"
    source "$WORKDIR/env.sh" 2>/dev/null || true
    python3 "$WORKDIR/score_results.py" \
        --all \
        --fixture-base "$FIXTURES_DIR" \
        --result-base "$RESULTS_DIR" \
        --outdir "$WORKDIR/analysis"
    python3 "$WORKDIR/plot_results.py" \
        --analysis-dir "$WORKDIR/analysis" \
        --result-base "$RESULTS_DIR" \
        --fixture-base "$FIXTURES_DIR"
fi
