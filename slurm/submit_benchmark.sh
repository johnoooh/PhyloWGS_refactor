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
#   --chains N          Parallel chains for Go impls (default: 4)
#   --time HH:MM:SS     Wall time limit CPU jobs (default: 6:00:00)
#   --gpu-time HH:MM:SS Wall time limit GPU job (default: 4:00:00)
#   --mem MB            Memory per job in MB (default: 8000)
#   --py2 CMD           Python 2 command for original impl (default: conda run -n phylo_py2 python)
#   --dry-run           Print jobs without submitting

set -euo pipefail

WORKDIR="${PHYLOWGS_WORKDIR:-$(pwd)/phylowgs_benchmark}"
BURNIN=1000
SAMPLES=2500
CHAINS=4
TIME_LIMIT="6:00:00"
GPU_TIME_LIMIT="4:00:00"
MEM_MB=8000
DRY_RUN=false
CSV_FILE=""
PY2_CMD="conda run -n phylo_py2 python"
PHYLOWGS_SIF=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --workdir)   WORKDIR="$2";       shift 2 ;;
        --burnin)    BURNIN="$2";        shift 2 ;;
        --samples)   SAMPLES="$2";       shift 2 ;;
        --chains)    CHAINS="$2";        shift 2 ;;
        --time)      TIME_LIMIT="$2";    shift 2 ;;
        --gpu-time)  GPU_TIME_LIMIT="$2";shift 2 ;;
        --mem)       MEM_MB="$2";        shift 2 ;;
        --py2)       PY2_CMD="$2";       shift 2 ;;
        --sif)       PHYLOWGS_SIF="$2";  shift 2 ;;
        --dry-run)   DRY_RUN=true;       shift ;;
        *.csv|*.tsv) CSV_FILE="$1";      shift ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

[[ -z "$CSV_FILE" ]]   && { echo "Usage: $0 samples.csv [options]"; exit 1; }
[[ ! -f "$CSV_FILE" ]] && { echo "CSV not found: $CSV_FILE"; exit 1; }
[[ ! -d "$WORKDIR" ]]  && { echo "Workdir not found: $WORKDIR — run setup.sh first"; exit 1; }

# Find Singularity image
if [[ -z "$PHYLOWGS_SIF" ]]; then
    # Check workdir first, then NXF cache
    if [[ -f "$WORKDIR/phylowgs_v1.5-msk.sif" ]]; then
        PHYLOWGS_SIF="$WORKDIR/phylowgs_v1.5-msk.sif"
    elif [[ -n "${NXF_SINGULARITY_CACHEDIR:-}" ]] && ls "${NXF_SINGULARITY_CACHEDIR}"/*phylowgs* &>/dev/null; then
        PHYLOWGS_SIF=$(ls "${NXF_SINGULARITY_CACHEDIR}"/*phylowgs* | head -1)
    else
        echo "ERROR: PhyloWGS Singularity image not found. Provide --sif /path/to/phylowgs.sif"; exit 1
    fi
fi
[[ ! -f "$PHYLOWGS_SIF" ]] && { echo "SIF not found: $PHYLOWGS_SIF"; exit 1; }
echo "SIF image: $PHYLOWGS_SIF"

RESULTS_DIR="$WORKDIR/results"
LOGS_DIR="$WORKDIR/logs"
INPUTS_DIR="$WORKDIR/inputs"
mkdir -p "$RESULTS_DIR" "$LOGS_DIR" "$INPUTS_DIR"

# ── Read samples CSV ──────────────────────────────────────────────────────────
declare -a SAMPLE_IDS MAFS FACETS_FILES
while IFS=',' read -r sample_id maf facets; do
    [[ "$sample_id" == "sample_id" ]] && continue
    [[ -z "$sample_id" ]] && continue
    SAMPLE_IDS+=("$sample_id")
    MAFS+=("$maf")
    FACETS_FILES+=("$facets")
done < "$CSV_FILE"

N_SAMPLES="${#SAMPLE_IDS[@]}"
echo "=== PhyloWGS Benchmark Submission ==="
echo "Samples:  $N_SAMPLES"
echo "Config:   B=$BURNIN s=$SAMPLES j=$CHAINS"
echo "Workdir:  $WORKDIR"
echo ""

# ── Step 1: Input conversion (one per sample) ─────────────────────────────────
declare -a CONVERT_JIDS
for i in "${!SAMPLE_IDS[@]}"; do
    sid="${SAMPLE_IDS[$i]}"
    maf="${MAFS[$i]}"
    facets="${FACETS_FILES[$i]}"
    outdir="$INPUTS_DIR/$sid"
    mkdir -p "$outdir"

    SCRIPT=$(cat << EOF
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

# Step 1: Convert FACETS cncf to PhyloWGS CNV format
python3 "$WORKDIR/facets_to_phylowgs_cnv.py" "$facets" -o "$outdir/cnv_for_phylowgs.txt"

# Step 2: Run create_phylowgs_inputs.py via container Python 2 (file deployed by setup.sh)
singularity exec --bind /data1 "$PHYLOWGS_SIF" \
    python2 "$WORKDIR/create_phylowgs_inputs.py" \
    --cnvs S1="$outdir/cnv_for_phylowgs.txt" \
    --output-cnvs "$outdir/cnv_data.txt" \
    --output-variants "$outdir/ssm_data.txt" \
    --output-params "$outdir/params.json" \
    --vcf-type S1=maf \
    S1="$maf"
echo "DONE convert: $sid"
EOF
)
    if [[ "$DRY_RUN" == true ]]; then
        echo "[dry-run] convert: $sid"
        CONVERT_JIDS+=("DRY_$i")
    else
        JID=$(echo "$SCRIPT" | sbatch --parsable)
        echo "  [convert] $sid → job $JID"
        CONVERT_JIDS+=("$JID")
    fi
done

# ── Step 2: Run jobs (4 implementations × N samples) ─────────────────────────
declare -a ALL_RUN_JIDS

run_job() {
    local impl="$1"
    local sid="$2"
    local convert_jid="$3"
    local out_dir="$RESULTS_DIR/$sid/$impl"
    mkdir -p "$out_dir"

    local script_body=""

    case "$impl" in

      original-python)
        script_body=$(cat << EOF
#SBATCH --job-name=phylo_orig_${sid}
#SBATCH --partition=cmobic_cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=${MEM_MB}M
#SBATCH --time=${TIME_LIMIT}
#SBATCH --output=${LOGS_DIR}/original-python_${sid}.out
#SBATCH --error=${LOGS_DIR}/original-python_${sid}.err

set -euo pipefail
source "$WORKDIR/env.sh" 2>/dev/null || true

echo "START: \$(date) | $impl | $sid"
START=\$(date +%s)

singularity exec --bind /data1 "$PHYLOWGS_SIF" \
    python2 /opt/phylowgs/evolve.py \
    -B $BURNIN -s $SAMPLES \
    -O "$out_dir" \
    "$INPUTS_DIR/$sid/ssm_data.txt" \
    "$INPUTS_DIR/$sid/cnv_data.txt"

END=\$(date +%s)
echo "{\"sample_id\":\"$sid\",\"impl\":\"$impl\",\"wall_seconds\":\$((END-START)),\"exit_code\":0}" \
    > "$out_dir/timing.json"
echo "END: \$(date) | elapsed \$((END-START))s"
EOF
)
        ;;

      optimized-python)
        script_body=$(cat << EOF
#SBATCH --job-name=phylo_opt_${sid}
#SBATCH --partition=cmobic_cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=${MEM_MB}M
#SBATCH --time=${TIME_LIMIT}
#SBATCH --output=${LOGS_DIR}/optimized-python_${sid}.out
#SBATCH --error=${LOGS_DIR}/optimized-python_${sid}.err

set -euo pipefail
source "$WORKDIR/env.sh" 2>/dev/null || true
source "$WORKDIR/impls/optimized-python/.venv/bin/activate"

# Extract mh.o from container if not already present
if [[ ! -f "$WORKDIR/impls/optimized-python/mh.o" ]]; then
    MH_SRC=\$(singularity exec "$PHYLOWGS_SIF" find / -name "mh.o" 2>/dev/null | head -1)
    singularity exec "$PHYLOWGS_SIF" cat "\$MH_SRC" > "$WORKDIR/impls/optimized-python/mh.o"
    chmod +x "$WORKDIR/impls/optimized-python/mh.o"
fi

echo "START: \$(date) | $impl | $sid"
START=\$(date +%s)

cd "$WORKDIR/impls/optimized-python"
python3 evolve.py \
    -B $BURNIN -s $SAMPLES \
    -O "$out_dir" \
    "$INPUTS_DIR/$sid/ssm_data.txt" \
    "$INPUTS_DIR/$sid/cnv_data.txt"

END=\$(date +%s)
echo "{\"sample_id\":\"$sid\",\"impl\":\"$impl\",\"wall_seconds\":\$((END-START)),\"exit_code\":0}" \
    > "$out_dir/timing.json"
echo "END: \$(date) | elapsed \$((END-START))s"
EOF
)
        ;;

      go-cpu)
        script_body=$(cat << EOF
#SBATCH --job-name=phylo_gcpu_${sid}
#SBATCH --partition=cmobic_cpu
#SBATCH --cpus-per-task=${CHAINS}
#SBATCH --mem=${MEM_MB}M
#SBATCH --time=${TIME_LIMIT}
#SBATCH --output=${LOGS_DIR}/go-cpu_${sid}.out
#SBATCH --error=${LOGS_DIR}/go-cpu_${sid}.err

set -euo pipefail
source "$WORKDIR/env.sh" 2>/dev/null || true

echo "START: \$(date) | $impl | $sid"
START=\$(date +%s)

"$WORKDIR/impls/go-cpu/phylowgs-cpu" \
    --no-gpu -B $BURNIN -s $SAMPLES -j $CHAINS \
    -O "$out_dir" \
    "$INPUTS_DIR/$sid/ssm_data.txt" \
    "$INPUTS_DIR/$sid/cnv_data.txt"

END=\$(date +%s)
echo "{\"sample_id\":\"$sid\",\"impl\":\"$impl\",\"wall_seconds\":\$((END-START)),\"exit_code\":0}" \
    > "$out_dir/timing.json"
echo "END: \$(date) | elapsed \$((END-START))s"
EOF
)
        ;;

      go-gpu)
        # Skip if CUDA unavailable
        if [[ -f "$WORKDIR/impls/go-gpu/.skip" ]]; then
            echo "  [skip] go-gpu: CUDA not available on this machine"
            return 0
        fi
        script_body=$(cat << EOF
#SBATCH --job-name=phylo_ggpu_${sid}
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=${CHAINS}
#SBATCH --mem=${MEM_MB}M
#SBATCH --time=${GPU_TIME_LIMIT}
#SBATCH --gres=gpu:1
#SBATCH --output=${LOGS_DIR}/go-gpu_${sid}.out
#SBATCH --error=${LOGS_DIR}/go-gpu_${sid}.err

set -euo pipefail
source "$WORKDIR/env.sh" 2>/dev/null || true
command -v module &>/dev/null && module load cuda/12 2>/dev/null || true

echo "START: \$(date) | $impl | $sid"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null || true
START=\$(date +%s)

LD_LIBRARY_PATH="$WORKDIR/impls/go-gpu" \
"$WORKDIR/impls/go-gpu/phylowgs-gpu" \
    -B $BURNIN -s $SAMPLES -j $CHAINS \
    -O "$out_dir" \
    "$INPUTS_DIR/$sid/ssm_data.txt" \
    "$INPUTS_DIR/$sid/cnv_data.txt"

END=\$(date +%s)
echo "{\"sample_id\":\"$sid\",\"impl\":\"$impl\",\"wall_seconds\":\$((END-START)),\"exit_code\":0}" \
    > "$out_dir/timing.json"
echo "END: \$(date) | elapsed \$((END-START))s"
EOF
)
        ;;
    esac

    local full_script="#!/usr/bin/env bash
${script_body}"

    if [[ "$DRY_RUN" == true ]]; then
        echo "[dry-run] run: $impl × $sid"
    else
        local JID
        JID=$(echo "$full_script" | sbatch --parsable --dependency="afterok:${convert_jid}")
        echo "  [run] $impl × $sid → job $JID"
        ALL_RUN_JIDS+=("$JID")
    fi
}

for i in "${!SAMPLE_IDS[@]}"; do
    sid="${SAMPLE_IDS[$i]}"
    conv_jid="${CONVERT_JIDS[$i]}"
    for impl in original-python optimized-python go-cpu go-gpu; do
        run_job "$impl" "$sid" "$conv_jid"
    done
done

# ── Step 3: Analysis job ───────────────────────────────────────────────────────
ALL_RUN_DEP=$(IFS=':'; echo "${ALL_RUN_JIDS[*]}")

ANALYSIS=$(cat << EOF
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
python3 "$WORKDIR/impls/go-cpu/slurm/analyze_benchmark.py" \
    --results-dir "$RESULTS_DIR" \
    --outdir "$WORKDIR/analysis" \
    --implementations original-python optimized-python go-cpu go-gpu

echo "END analysis: \$(date)"
echo "Report: $WORKDIR/analysis/report.html"
EOF
)

if [[ "$DRY_RUN" == true ]]; then
    echo ""
    echo "[dry-run] analysis depends on: ${ALL_RUN_JIDS[*]:-none}"
else
    ANA_JID=$(echo "$ANALYSIS" | sbatch --parsable --dependency="afterok:${ALL_RUN_DEP}")
    echo ""
    echo "=== All jobs submitted ==="
    echo "  Convert jobs: ${#CONVERT_JIDS[@]}"
    echo "  Run jobs:     ${#ALL_RUN_JIDS[@]} (${N_SAMPLES} samples × 4 impls)"
    echo "  Analysis:     job $ANA_JID (depends on all)"
    echo ""
    echo "Monitor:  squeue -u \$USER"
    echo "Results:  $WORKDIR/analysis/report.html"
fi
