#!/usr/bin/env bash
# =============================================================================
# submit_full_suite.sh — HPC analog of scripts/Studio/run_batch_tests.sh
#
# Submits the same default curated test suite that the local Studio batch
# runner uses, but on the cluster. Each test is fanned out as its own sbatch
# job and they all land in a single shared parent folder under $SCRATCH for
# easy bulk download. This is a thin wrapper around submit_batch.sh: all the
# resource-sizing, build, and BATCH_OUT_DIR plumbing lives there.
#
# Usage (run from project root):
#
#   ./scripts/HPC/submit_full_suite.sh                        # default suite
#   ./scripts/HPC/submit_full_suite.sh --tag mylabel
#   ./scripts/HPC/submit_full_suite.sh --skip-1d              # 2D only
#   ./scripts/HPC/submit_full_suite.sh --skip-hires           # exclude *_hires geometries
#   ./scripts/HPC/submit_full_suite.sh --tests "g1:e1,g2:e2"  # override list
#   ./scripts/HPC/submit_full_suite.sh --dry-run              # print resolved list, do not submit
#
# Append extra sbatch flags after --:
#   ./scripts/HPC/submit_full_suite.sh --tag overnight -- --time=0-08:00:00
#
# Design note: the test list mirrors scripts/Studio/run_batch_tests.sh::DEFAULT_TESTS
# byte-for-byte. Keep both in sync when adding new canonical tests so the
# local-vs-HPC workflows produce comparable results.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
SUBMIT_BATCH="$SCRIPT_DIR/submit_batch.sh"

if [[ ! -x "$SUBMIT_BATCH" ]]; then
    echo "ERROR: submit_batch.sh not executable at $SUBMIT_BATCH"
    exit 1
fi

# -----------------------------------------------------------------------------
# Default test list — mirrors scripts/Studio/run_batch_tests.sh::DEFAULT_TESTS.
# Format: "<geometry>:<experiment>" (no .opts extension)
# -----------------------------------------------------------------------------
DEFAULT_TESTS=(
    # ---- 1D tests at h=0.95, T=-20 ----
    "1D_ice_slab:1day_T-20_h0.95"
    "1D_ice_slab_hires:1day_T-20_h0.95"
    "1D_single_ice:1day_T-20_h0.95"
    "1D_single_ice_hires:1day_T-20_h0.95"
    "1D_ice_sed_pair:1day_T-20_h0.95"
    "1D_ice_sed_pair_hires:1day_T-20_h0.95"
    "1D_separated_grains:1day_T-20_h0.95"
    "1D_separated_grains_hires:1day_T-20_h0.95"
    "1D_touching_grains:1day_T-20_h0.95"
    "1D_touching_grains_hires:1day_T-20_h0.95"

    # ---- 2D tests at saturation (multi-grain physics — Ostwald ripening, sintering) ----
    "2D_ice_sed_pair:1day_T-20_h1.00"
    "2D_ice_sed_pair_hires:1day_T-20_h1.00"
    "2D_touching_grains:1day_T-20_h1.00"
    "2D_touching_grains_hires:1day_T-20_h1.00"
    "2D_separated_grains:1day_T-20_h1.00"
    "2D_separated_grains_hires:1day_T-20_h1.00"
    "2D_single_sed:1day_T-20_h1.00"

    # ---- 2D tests at undersaturation (bulk sublimation driving force) ----
    "2D_single_ice:1day_T-20_h0.95"
    "2D_single_ice_hires:1day_T-20_h0.95"
    "2D_ice_slab:1day_T-20_h0.95"
    "2D_ice_slab_hires:1day_T-20_h0.95"
)

# -----------------------------------------------------------------------------
# Arg parsing
# -----------------------------------------------------------------------------
tag=""
skip_1d=false
skip_hires=false
custom_tests=""
dry_run=false
sbatch_extra=()

usage() { sed -n '2,30p' "$0"; exit 1; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        --tag)         tag="$2"; shift 2 ;;
        --tests)       custom_tests="$2"; shift 2 ;;
        --skip-1d)     skip_1d=true; shift ;;
        --skip-hires)  skip_hires=true; shift ;;
        --dry-run)     dry_run=true; shift ;;
        --)            shift; sbatch_extra=("$@"); break ;;
        -h|--help)     usage ;;
        *)             echo "Unknown argument: $1"; usage ;;
    esac
done

# Build the final list
TESTS=()
if [[ -n "$custom_tests" ]]; then
    IFS=',' read -ra TESTS <<< "$custom_tests"
    for i in "${!TESTS[@]}"; do
        s="${TESTS[$i]}"
        s="${s#"${s%%[![:space:]]*}"}"
        s="${s%"${s##*[![:space:]]}"}"
        TESTS[$i]="$s"
    done
else
    for t in "${DEFAULT_TESTS[@]}"; do
        $skip_1d    && [[ "$t" == 1D_*       ]] && continue
        $skip_hires && [[ "$t" == *_hires:*  ]] && continue
        TESTS+=("$t")
    done
fi

if [[ ${#TESTS[@]} -eq 0 ]]; then
    echo "❌ No tests selected after filtering."
    exit 1
fi

# Render the comma-separated list submit_batch.sh expects
joined=""
for t in "${TESTS[@]}"; do
    joined+="${joined:+,}$t"
done

echo "============================================================"
echo "  Permafrost full-suite HPC submission"
echo "  Tag        : ${tag:-<none>}"
echo "  Tests      : ${#TESTS[@]}"
echo "  Filters    : skip_1d=$skip_1d  skip_hires=$skip_hires"
echo "  Extra sbatch flags: ${sbatch_extra[*]:-<none>}"
echo "============================================================"
for t in "${TESTS[@]}"; do echo "   - $t"; done
echo "------------------------------------------------------------"

if $dry_run; then
    echo "DRY RUN — not invoking submit_batch.sh."
    exit 0
fi

# Hand off to submit_batch.sh — it owns the build, BATCH_OUT_DIR placement,
# per-geometry resource sizing, and SLURM submission.
cd "$PROJECT_ROOT"
if [[ ${#sbatch_extra[@]} -gt 0 ]]; then
    exec "$SUBMIT_BATCH" \
        ${tag:+--tag "$tag"} \
        --tests "$joined" \
        -- "${sbatch_extra[@]}"
else
    exec "$SUBMIT_BATCH" \
        ${tag:+--tag "$tag"} \
        --tests "$joined"
fi
