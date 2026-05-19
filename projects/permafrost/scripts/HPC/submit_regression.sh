#!/usr/bin/env bash
# =============================================================================
# submit_regression.sh — submit the PenaltyWeight / k_pen=0 regression sweep.
#
# Submits every job in parallel via sbatch. Each job uses the standard
# `run_permafrost.sh <geom> <exp> <tag>` interface and inherits the #SBATCH
# directives at the top of that script.
#
# Run from project root:
#   ./scripts/HPC/submit_regression.sh
#
# Override the tag (appended to the run folder name) with the first arg:
#   ./scripts/HPC/submit_regression.sh kpen0_check
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
RUN_SCRIPT="$SCRIPT_DIR/run_permafrost.sh"

TAG="${1:-kpen0_regression}"

# Quick regression: 1 day, -20°C, saturated. All geometries below should run
# cleanly with k_pen=0 and the new PenaltyWeight() (active only deep in solid).
REGRESSION_EXP="1day_T-20_h1.00"
REGRESSION_GEOMS=(
    "1D_ice_sed_pair"
    "1D_ice_slab"
    "1D_separated_grains"
    "1D_single_ice"
    "1D_touching_grains"
    "2D_single_ice"
    "2D_single_sed"
    "2D_ice_sed_pair"
    "2D_touching_grains"
)

# Diagnostic: the wide-separation Ostwald ripening run (30 day, -5°C).
DIAGNOSTIC_EXP="30day_T-5_h1.00"
DIAGNOSTIC_GEOM="2D_separated_grains"

cd "$PROJECT_ROOT"

if [[ ! -f permafrost ]]; then
    echo "❌ ./permafrost not found. Build first (make all)."
    exit 1
fi

echo "============================================================"
echo "  Permafrost regression sweep"
echo "  Tag        : $TAG"
echo "  Regression : $REGRESSION_EXP  (${#REGRESSION_GEOMS[@]} jobs)"
echo "  Diagnostic : $DIAGNOSTIC_EXP on $DIAGNOSTIC_GEOM"
echo "============================================================"

submit_one() {
    local geom="$1" exp="$2" tag="$3"
    local job_name="${geom}__${exp}"
    echo "→ sbatch -J $job_name : geom=$geom exp=$exp tag=$tag"
    sbatch --job-name="$job_name" "$RUN_SCRIPT" "$geom" "$exp" "$tag"
}

for geom in "${REGRESSION_GEOMS[@]}"; do
    submit_one "$geom" "$REGRESSION_EXP" "$TAG"
done

submit_one "$DIAGNOSTIC_GEOM" "$DIAGNOSTIC_EXP" "wideSep_${TAG}"

echo "============================================================"
echo "  Submitted ${#REGRESSION_GEOMS[@]} regression jobs + 1 diagnostic."
echo "  Check with: squeue -u \$USER"
echo "============================================================"
