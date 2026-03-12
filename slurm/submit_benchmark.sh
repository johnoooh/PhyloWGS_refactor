#!/usr/bin/env bash
# submit_benchmark.sh — Submit PhyloWGS benchmark SLURM jobs for all implementations
#
# Usage:
#   bash submit_benchmark.sh samples.csv --workdir /path/to/workdir [options]
#
# CSV format (header required):
#   sample_id,maf,facets
#
# Options:
#   --workdir DIR       Workdir created by setup.sh (default: ./phylowgs_benchmark)
#   --burnin N          MCMC burn-in iterations (default: 1000)
#   --samples N         MCMC samples (default: 2500)
#   --chains N          Parallel chains (default: 4)
#   --time HH:MM:SS     Wall time limit (default: 4:00:00)
#   --mem MB            Memory per job in MB (default: 8000)
#   --dry-run           Print jobs without submitting

set -euo pipefail

# ── Defaults ─────────────────────────────────────────────────────────────────
WORKDIR="${PHYLOWGS_WORKDIR:-$(pwd)/phylowgs_benchmark}"
BURNIN=1000
SAMPLES=2500
CHAINS=4
TIME_LIMIT="4:00:00"
MEM_MB=8000
DRY_RUN=false
CSV_FILE=""

# ── Parse args ────────────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case $1 in
        --workdir)  WORKDIR="$2";    shift 2 ;;
        --burnin)   BURNIN="$2";     shift 2 ;;
        --samples)  SAMPLES="$2";    shift 2 ;;
        --chains)   CHAINS="$2";     shift 2 ;;
        --time)     TIME_LIMIT="$2"; shift 2 ;;
        --mem)      MEM_MB="$2";     shift 2 ;;
        --dry-run)  DRY_RUN=true;    shift ;;
        *.csv|*.tsv) CSV_FILE="$1";  shift ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

[[ -z "$CSV_FILE" ]] && { echo "Usage: $0 samples.csv [options]"; exit 1; }
[[ ! -f "$CSV_FILE" ]] && { echo "CSV not found: $CSV_FILE"; exit 1; }
[[ ! -d "$WORKDIR" ]] && { echo "Workdir not found: $WORKDIR — run setup.sh first"; exit 1; }

RESULTS_DIR="$WORKDIR/results"
LOGS_DIR="$WORKDIR/logs"
INPUTS_DIR="$WORKDIR/inputs"
mkdir -p "$RESULTS_DIR" "$LOGS_DIR" "$INPUTS_DIR"

# ── Read samples CSV ──────────────────────────────────────────────────────────
declare -a SAMPLE_IDS MAFS FACETS_FILES
while IFS=',' read -r sample_id maf facets; do
    [[ "$sample_id" == "sample_id" ]] && continue  # skip header
    [[ -z "$sample_id" ]] && continue
    SAMPLE_IDS+=("$sample_id")
    MAFS+=("$maf")
    FACETS_FILES+=("$facets")
done < "$CSV_FILE"

N_SAMPLES="${#SAMPLE_IDS[@]}"
echo "=== PhyloWGS Benchmark Submission ==="
echo "Samples: $N_SAMPLES | B=$BURNIN s=$SAMPLES j=$CHAINS"
echo "Workdir: $WORKDIR"
echo ""

# ── Read implementations manifest ────────────────────────────────────────────
declare -a IMPL_NAMES IMPL_TYPES IMPL_BINARIES IMPL_PARTITIONS IMPL_EXTRA_ARGS
while IFS=$'\t' read -r name type binary branch partition extra_args; do
    [[ "$name" == "name" ]] && continue
    IMPL_NAMES+=("$name")
    IMPL_TYPES+=("$type")
    IMPL_BINARIES+=("$binary")
    IMPL_PARTITIONS+=("$partition")
    IMPL_EXTRA_ARGS+=("$extra_args")
done < "$WORKDIR/implementations.tsv"

# ── Step 1: Input conversion jobs (one per sample) ───────────────────────────
declare -a CONVERT_JIDS
for i in "${!SAMPLE_IDS[@]}"; do
    sid="${SAMPLE_IDS[$i]}"
    maf="${MAFS[$i]}"
    facets="${FACETS_FILES[$i]}"
    outdir="$INPUTS_DIR/$sid"
    mkdir -p "$outdir"

    CONVERT_SCRIPT=$(cat << EOF
#!/usr/bin/env bash
#SBATCH --job-name=phylo_convert_${sid}
#SBATCH --partition=cmobic_cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0:30:00
#SBATCH --output=${LOGS_DIR}/convert_${sid}.out
#SBATCH --error=${LOGS_DIR}/convert_${sid}.err

set -euo pipefail
source "$WORKDIR/env.sh" 2>/dev/null || true
cd "$WORKDIR"

python3 convert_inputs.py \
    --sample-id "$sid" \
    --maf "$maf" \
    --facets "$facets" \
    --outdir "$outdir"

echo "DONE: $sid"
EOF
)
    if [[ "$DRY_RUN" == true ]]; then
        echo "[dry-run] convert: $sid"
        CONVERT_JIDS+=("DRY_$i")
    else
        JID=$(echo "$CONVERT_SCRIPT" | sbatch --parsable)
        echo "  [convert] $sid → job $JID"
        CONVERT_JIDS+=("$JID")
    fi
done

# Dependency string for all conversion jobs
CONVERT_DEP=$(IFS=':'; echo "${CONVERT_JIDS[*]}")

# ── Step 2: PhyloWGS run jobs (per sample × per implementation) ───────────────
declare -a ALL_RUN_JIDS
for j in "${!IMPL_NAMES[@]}"; do
    impl="${IMPL_NAMES[$j]}"
    itype="${IMPL_TYPES[$j]}"
    binary="${IMPL_BINARIES[$j]}"
    partition="${IMPL_PARTITIONS[$j]}"
    extra="${IMPL_EXTRA_ARGS[$j]}"

    impl_dir="$WORKDIR/impls/$impl"

    declare -a IMPL_JIDS
    for i in "${!SAMPLE_IDS[@]}"; do
        sid="${SAMPLE_IDS[$i]}"
        input_dir="$INPUTS_DIR/$sid"
        out_dir="$RESULTS_DIR/$sid/$impl"
        mkdir -p "$out_dir"

        # GPU-specific SLURM flags
        GPU_FLAGS=""
        if [[ "$partition" == "gpu" ]]; then
            GPU_FLAGS="#SBATCH --gres=gpu:1"
        fi

        # Build the run command
        if [[ "$itype" == "python" ]]; then
            RUN_CMD="source \"$impl_dir/.venv/bin/activate\"
python3 \"$impl_dir/evolve.py\" \
    -B $BURNIN -s $SAMPLES \
    -O \"$out_dir\" \
    \"$input_dir/ssm_data.txt\" \
    \"$input_dir/cnv_data.txt\""
        else
            # Go binary
            LD_PREFIX=""
            if [[ "$partition" == "gpu" ]]; then
                LD_PREFIX="LD_LIBRARY_PATH=\"$impl_dir/cuda\""
            fi
            RUN_CMD="${LD_PREFIX} \"$impl_dir/$binary\" \
    -B $BURNIN -s $SAMPLES -j $CHAINS \
    $extra \
    -O \"$out_dir\" \
    \"$input_dir/ssm_data.txt\" \
    \"$input_dir/cnv_data.txt\""
        fi

        RUN_SCRIPT=$(cat << EOF
#!/usr/bin/env bash
#SBATCH --job-name=phylo_${impl}_${sid}
#SBATCH --partition=${partition}
#SBATCH --cpus-per-task=${CHAINS}
#SBATCH --mem=${MEM_MB}M
#SBATCH --time=${TIME_LIMIT}
#SBATCH --output=${LOGS_DIR}/${impl}_${sid}.out
#SBATCH --error=${LOGS_DIR}/${impl}_${sid}.err
${GPU_FLAGS}

set -euo pipefail
source "$WORKDIR/env.sh" 2>/dev/null || true

echo "START: \$(date)"
echo "Impl: $impl | Sample: $sid"

START_TIME=\$(date +%s)

$RUN_CMD

END_TIME=\$(date +%s)
ELAPSED=\$((END_TIME - START_TIME))

# Write timing file for analysis
echo "{\"sample_id\": \"$sid\", \"impl\": \"$impl\", \"wall_seconds\": \$ELAPSED, \"exit_code\": 0}" \
    > "$out_dir/timing.json"

echo "END: \$(date) | Elapsed: \${ELAPSED}s"
EOF
)
        if [[ "$DRY_RUN" == true ]]; then
            echo "[dry-run] run: $impl × $sid"
            IMPL_JIDS+=("DRY_${j}_${i}")
        else
            JID=$(echo "$RUN_SCRIPT" | sbatch --parsable --dependency="afterok:${CONVERT_JIDS[$i]}")
            echo "  [run] $impl × $sid → job $JID (after convert ${CONVERT_JIDS[$i]})"
            IMPL_JIDS+=("$JID")
            ALL_RUN_JIDS+=("$JID")
        fi
    done
done

# ── Step 3: Analysis job (depends on all run jobs) ────────────────────────────
ALL_RUN_DEP=$(IFS=':'; echo "${ALL_RUN_JIDS[*]}")

ANALYSIS_SCRIPT=$(cat << EOF
#!/usr/bin/env bash
#SBATCH --job-name=phylo_analyze
#SBATCH --partition=cmobic_cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --output=${LOGS_DIR}/analysis.out
#SBATCH --error=${LOGS_DIR}/analysis.err

set -euo pipefail
source "$WORKDIR/env.sh" 2>/dev/null || true
source "$WORKDIR/impls/optimized-python/.venv/bin/activate"

echo "START analysis: \$(date)"
python3 "$WORKDIR/impls/go-cpu-opt/slurm/analyze_benchmark.py" \
    --results-dir "$RESULTS_DIR" \
    --outdir "$WORKDIR/analysis" \
    --implementations $(IFS=' '; echo "${IMPL_NAMES[*]}")

echo "END analysis: \$(date)"
echo "Report: $WORKDIR/analysis/report.html"
EOF
)

if [[ "$DRY_RUN" == true ]]; then
    echo ""
    echo "[dry-run] analysis job would depend on: ${ALL_RUN_JIDS[*]:-none}"
else
    ANA_JID=$(echo "$ANALYSIS_SCRIPT" | sbatch --parsable --dependency="afterok:${ALL_RUN_DEP}")
    echo ""
    echo "=== All jobs submitted ==="
    echo "  Input conversion: ${#CONVERT_JIDS[@]} jobs"
    echo "  PhyloWGS runs:    ${#ALL_RUN_JIDS[@]} jobs (${N_SAMPLES} samples × ${#IMPL_NAMES[@]} impls)"
    echo "  Analysis:         job $ANA_JID (depends on all runs)"
    echo ""
    echo "Monitor:"
    echo "  squeue -u \$USER"
    echo "  sacct -j $ANA_JID --format=JobID,State,Elapsed,MaxRSS"
fi
