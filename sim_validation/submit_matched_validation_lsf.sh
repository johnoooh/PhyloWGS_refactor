#!/usr/bin/env bash
# submit_matched_validation_lsf.sh — MATCHED Go-vs-Python validation sweep (juno/LSF)
# ============================================================================
# WHY THIS SCRIPT EXISTS (read before running)
# ----------------------------------------------------------------------------
# A previous validation effort produced a bogus "Go under-recovers at K=10"
# finding. That finding was a HARNESS artifact, not a sampler bug. The root
# cause: Go and Python results came from TWO SEPARATE SWEEPS that used
# DIFFERENT base-seeds, so the two implementations were never run on the same
# simulated fixtures. Comparing them was apples-to-oranges.
#
# This script makes that mistake STRUCTURALLY IMPOSSIBLE:
#
#   1. Fixtures are generated EXACTLY ONCE, via setup_validation.sh, into a
#      single shared $WORKDIR/fixtures dir, from a single --base-seed.
#   2. BOTH implementations (go-cpu AND original-python) are then submitted
#      against that ONE fixtures dir in ONE sweep. --go-only is intentionally
#      NOT exposed here — you cannot accidentally run a Go-only half-sweep.
#   3. The shared fixtures dir is fingerprinted (manifest sha) and the
#      fingerprint is written into the workdir, so a later re-run that points
#      at a different fixtures set is detected loudly.
#
# Net effect: every Go result has a Python result on the IDENTICAL fixture,
# generated from the IDENTICAL seed. The comparison is matched by construction.
#
# This is a thin, opinionated wrapper around the proven entrypoints:
#   - setup_validation.sh            (one-time fixture generation + Go build)
#   - submit_sim_validation_lsf.sh   (LSF array submission, both impls)
# It does NOT re-implement submission logic or touch any Go/Python source.
#
# ----------------------------------------------------------------------------
# USAGE
# ----------------------------------------------------------------------------
#   bash submit_matched_validation_lsf.sh \
#       --workdir /juno/work/.../matched_sweep \
#       --sif     /path/to/phylowgs.sif \
#       [--grid default] [--base-seed 42] \
#       [--burnin 1000] [--samples 2500] [--chains 8] \
#       [--queue general] [--mem 8] [--dry-run] [--skip-setup]
#
# Options:
#   --workdir DIR     Shared workdir (created if absent). Holds the ONE
#                     fixtures/ dir both impls run against. REQUIRED.
#   --sif PATH        PhyloWGS singularity image for original Python. REQUIRED
#                     (omit only if you truly want Go-only, which defeats the
#                     purpose — the script will refuse unless --allow-go-only).
#   --grid NAME       Fixture grid: quick|default|thorough (default: default).
#   --base-seed N     SINGLE base seed for fixture generation (default: 42).
#                     Both impls share these fixtures; the seed is recorded.
#   --burnin N        MCMC burn-in, applied to BOTH impls (default: 1000).
#   --samples N       MCMC samples, applied to BOTH impls (default: 2500).
#   --chains N        Parallel chains, applied to BOTH impls (default: 8).
#   --queue NAME      LSF queue (default: general).
#   --mem GB          Memory per job in GB (default: 8).
#   --skip-setup      Reuse an existing $WORKDIR/fixtures (still validates the
#                     fingerprint). Use for resubmits after a partial run.
#   --allow-go-only   Escape hatch: permit a Go-only submission (NOT matched;
#                     use only for a quick Go smoke test). Off by default.
#   --dry-run         Print what would happen; submit nothing.
#
# After completion, compare results in:
#   $WORKDIR/analysis/comparison/     (matched Go-vs-Python report)
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Defaults ─────────────────────────────────────────────────────────────────
WORKDIR=""
PHYLOWGS_SIF=""
GRID="default"
BASE_SEED=42
BURNIN=1000
SAMPLES=2500
CHAINS=8
QUEUE="general"
MEM_GB=8
SKIP_SETUP=false
ALLOW_GO_ONLY=false
DRY_RUN=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --workdir)       WORKDIR="$2";       shift 2 ;;
        --sif)           PHYLOWGS_SIF="$2";  shift 2 ;;
        --grid)          GRID="$2";          shift 2 ;;
        --base-seed)     BASE_SEED="$2";     shift 2 ;;
        --seed)          BASE_SEED="$2";     shift 2 ;;  # alias
        --burnin)        BURNIN="$2";        shift 2 ;;
        --samples)       SAMPLES="$2";       shift 2 ;;
        --chains)        CHAINS="$2";        shift 2 ;;
        --queue)         QUEUE="$2";         shift 2 ;;
        --mem)           MEM_GB="$2";        shift 2 ;;
        --skip-setup)    SKIP_SETUP=true;    shift ;;
        --allow-go-only) ALLOW_GO_ONLY=true; shift ;;
        --dry-run)       DRY_RUN=true;       shift ;;
        -h|--help)
            sed -n '2,70p' "${BASH_SOURCE[0]}"; exit 0 ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

# ── Validate required inputs ─────────────────────────────────────────────────
[[ -z "$WORKDIR" ]] && { echo "ERROR: --workdir is required."; exit 1; }

if [[ -z "$PHYLOWGS_SIF" && "$ALLOW_GO_ONLY" == false ]]; then
    echo "ERROR: --sif is required for a MATCHED sweep (Go + Python on the same"
    echo "       fixtures). A matched validation needs the Python baseline."
    echo "       If you really want a Go-only smoke test, pass --allow-go-only."
    exit 1
fi
if [[ -n "$PHYLOWGS_SIF" && ! -f "$PHYLOWGS_SIF" ]]; then
    echo "ERROR: SIF not found: $PHYLOWGS_SIF"; exit 1
fi

SETUP_SH="$SCRIPT_DIR/setup_validation.sh"
SUBMIT_SH="$SCRIPT_DIR/submit_sim_validation_lsf.sh"
[[ ! -f "$SETUP_SH" ]]  && { echo "ERROR: missing $SETUP_SH"; exit 1; }
[[ ! -f "$SUBMIT_SH" ]] && { echo "ERROR: missing $SUBMIT_SH"; exit 1; }

FIXTURES_DIR="$WORKDIR/fixtures"
MANIFEST="$FIXTURES_DIR/manifest.json"
SEED_RECORD="$WORKDIR/.matched_seed"
FINGERPRINT_RECORD="$WORKDIR/.fixtures_fingerprint"

echo "=== MATCHED Validation Sweep (LSF) ==="
echo "Workdir:     $WORKDIR"
echo "Grid:        $GRID"
echo "Base seed:   $BASE_SEED  (single shared seed — both impls)"
echo "MCMC:        B=$BURNIN s=$SAMPLES chains=$CHAINS  (identical for Go + Python)"
echo "SIF:         ${PHYLOWGS_SIF:-<none — go-only escape hatch>}"
echo "Queue/mem:   $QUEUE / ${MEM_GB}GB"
echo ""

# ── Step 1: generate the ONE shared fixtures dir (or reuse + verify) ─────────
sha_of_manifest() {
    # Portable manifest fingerprint (Linux: sha256sum; macOS: shasum -a 256).
    if command -v sha256sum &>/dev/null; then
        sha256sum "$MANIFEST" | awk '{print $1}'
    elif command -v shasum &>/dev/null; then
        shasum -a 256 "$MANIFEST" | awk '{print $1}'
    else
        # Last resort: byte count (still detects a different fixture set).
        wc -c < "$MANIFEST" | tr -d ' '
    fi
}

if [[ "$SKIP_SETUP" == true ]]; then
    echo "[fixtures] --skip-setup: reusing existing $FIXTURES_DIR"
    [[ ! -f "$MANIFEST" ]] && { echo "ERROR: --skip-setup but no manifest at $MANIFEST"; exit 1; }
else
    if [[ -f "$MANIFEST" ]]; then
        echo "[fixtures] WARNING: $FIXTURES_DIR already exists."
        echo "           A matched sweep must use ONE consistent fixtures set."
        echo "           Keeping the existing fixtures (will verify fingerprint)."
        echo "           To regenerate from scratch, delete $WORKDIR and rerun."
    else
        echo "[fixtures] Generating the single shared fixtures dir via setup_validation.sh"
        if [[ "$DRY_RUN" == true ]]; then
            echo "  [dry-run] bash $SETUP_SH --workdir $WORKDIR --grid $GRID --seed $BASE_SEED"
        else
            bash "$SETUP_SH" --workdir "$WORKDIR" --grid "$GRID" --seed "$BASE_SEED"
        fi
    fi
fi

# Record / verify the single base-seed used for these fixtures.
if [[ "$DRY_RUN" != true ]]; then
    if [[ -f "$SEED_RECORD" ]]; then
        RECORDED_SEED="$(cat "$SEED_RECORD")"
        if [[ "$RECORDED_SEED" != "$BASE_SEED" ]]; then
            echo "ERROR: this workdir's fixtures were built with base-seed"
            echo "       $RECORDED_SEED, but you passed --base-seed $BASE_SEED."
            echo "       Mixing seeds is exactly the bug this script prevents."
            echo "       Use the recorded seed, or use a fresh --workdir."
            exit 1
        fi
    else
        echo "$BASE_SEED" > "$SEED_RECORD"
    fi

    # Fingerprint the manifest so an accidental fixtures swap is caught.
    if [[ -f "$MANIFEST" ]]; then
        CUR_FP="$(sha_of_manifest)"
        if [[ -f "$FINGERPRINT_RECORD" ]]; then
            OLD_FP="$(cat "$FINGERPRINT_RECORD")"
            if [[ "$CUR_FP" != "$OLD_FP" ]]; then
                echo "ERROR: fixtures manifest fingerprint changed for this workdir."
                echo "       Recorded: $OLD_FP"
                echo "       Current:  $CUR_FP"
                echo "       The fixtures set was modified — a matched comparison"
                echo "       across runs is no longer guaranteed. Use a fresh workdir."
                exit 1
            fi
        else
            echo "$CUR_FP" > "$FINGERPRINT_RECORD"
        fi
        echo "[fixtures] OK — shared manifest fingerprint: $CUR_FP"
    fi
fi

# ── Step 2: submit BOTH impls against the ONE shared fixtures dir ────────────
echo ""
echo "[submit] Submitting Go AND Python on the SAME fixtures (single sweep)"

SUBMIT_ARGS=(
    --workdir "$WORKDIR"
    --burnin  "$BURNIN"
    --samples "$SAMPLES"
    --chains  "$CHAINS"
    --queue   "$QUEUE"
    --mem     "$MEM_GB"
)

if [[ -n "$PHYLOWGS_SIF" ]]; then
    SUBMIT_ARGS+=(--sif "$PHYLOWGS_SIF")
elif [[ "$ALLOW_GO_ONLY" == true ]]; then
    echo "[submit] WARNING: --allow-go-only set; this is NOT a matched sweep."
    SUBMIT_ARGS+=(--go-only)
fi

[[ "$DRY_RUN" == true ]] && SUBMIT_ARGS+=(--dry-run)

echo "  bash $SUBMIT_SH ${SUBMIT_ARGS[*]}"
echo ""
bash "$SUBMIT_SH" "${SUBMIT_ARGS[@]}"

echo ""
echo "=== Matched sweep dispatched ==="
echo "Both impls ran on:  $FIXTURES_DIR  (seed $BASE_SEED)"
echo "Monitor:            bjobs -w | grep simval"
echo "Matched comparison: $WORKDIR/analysis/comparison/"
