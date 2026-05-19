#!/usr/bin/env bash
# =============================================================================
# submit_regression.sh — submit the PenaltyWeight / k_pen=0 regression sweep.
#
# Submits every job in parallel via sbatch. Each job sizes its allocation to
# the geometry's actual problem size (target: TARGET_DOFS_PER_CORE DoFs/core),
# overriding the default #SBATCH directives in run_permafrost.sh so 1D jobs
# don't grab a full 64-core allocation.
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

# Resource-sizing parameters — keep in sync with submit_permafrost.sh.
TARGET_DOFS_PER_CORE=10000          # rule-of-thumb DoFs/core target
MAX_TASKS_PER_NODE=32               # safe minimum across icelake|skylake|cascadelake

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

# Build once on the submission host so the fanned-out jobs don't all race to
# `make clean && make all` in the shared obj/ directory (which produces
# "Stale file handle" / "No space left on device" errors when one job's
# `make clean` deletes obj/*.o while another job is mid-write).
echo ""
echo "--- Building permafrost on submission host ---"
if ! make all; then
    echo "❌ Build failed. Fix the build before submitting jobs."
    exit 1
fi
if [[ ! -x ./permafrost ]]; then
    echo "❌ ./permafrost still missing after make all."
    exit 1
fi
echo "✅ Build complete: $(ls -la ./permafrost | awk '{print $5, $6, $7, $8}')"

# Compute optimal (nprocs, nnodes, ntasks_per_node) for a geometry's grid size.
# Echoes the three values space-separated; caller reads with `read`.
compute_alloc() {
    local geom_file="$1"
    local nx ny nz
    nx=$(awk '$1=="-Nx"{print $2}' "$geom_file" | head -n1); nx=${nx:-1}
    ny=$(awk '$1=="-Ny"{print $2}' "$geom_file" | head -n1); ny=${ny:-1}
    nz=$(awk '$1=="-Nz"{print $2}' "$geom_file" | head -n1); nz=${nz:-1}

    local total_dofs=$((4 * nx * ny * nz))
    local nprocs=$(( (total_dofs + TARGET_DOFS_PER_CORE - 1) / TARGET_DOFS_PER_CORE ))
    (( nprocs < 1 )) && nprocs=1

    # Cap tasks-per-node to nprocs so a 1-core 1D job doesn't reserve 32 slots.
    local tasks_per_node=$nprocs
    (( tasks_per_node > MAX_TASKS_PER_NODE )) && tasks_per_node=$MAX_TASKS_PER_NODE

    local nnodes=$(( (nprocs + tasks_per_node - 1) / tasks_per_node ))
    (( nnodes < 1 )) && nnodes=1

    echo "$nprocs $nnodes $tasks_per_node $total_dofs"
}

submit_one() {
    local geom="$1" exp="$2" tag="$3"
    local geom_file="$PROJECT_ROOT/inputs/geometry/${geom}.opts"
    local job_name="${geom}__${exp}"

    if [[ ! -f "$geom_file" ]]; then
        echo "⚠ Skipping $geom — geometry file not found: $geom_file"
        return
    fi

    local nprocs nnodes tasks_per_node total_dofs
    read -r nprocs nnodes tasks_per_node total_dofs < <(compute_alloc "$geom_file")

    printf "→ %-25s DoFs=%-8d nprocs=%-3d nodes=%-2d tasks/node=%d\n" \
        "$job_name" "$total_dofs" "$nprocs" "$nnodes" "$tasks_per_node"

    # --export=ALL,SKIP_COMPILE=1 forwards the submitter's environment and adds
    # the skip-compile flag — see compile_code() in run_permafrost.sh.
    sbatch --job-name="$job_name" \
           --nodes="$nnodes" \
           --ntasks="$nprocs" \
           --ntasks-per-node="$tasks_per_node" \
           --export=ALL,SKIP_COMPILE=1 \
           "$RUN_SCRIPT" "$geom" "$exp" "$tag"
}

echo "============================================================"
echo "  Permafrost regression sweep"
echo "  Tag        : $TAG"
echo "  Target     : ${TARGET_DOFS_PER_CORE} DoFs/core (max ${MAX_TASKS_PER_NODE}/node)"
echo "  Regression : $REGRESSION_EXP  (${#REGRESSION_GEOMS[@]} jobs)"
echo "  Diagnostic : $DIAGNOSTIC_EXP on $DIAGNOSTIC_GEOM"
echo "============================================================"

for geom in "${REGRESSION_GEOMS[@]}"; do
    submit_one "$geom" "$REGRESSION_EXP" "$TAG"
done

submit_one "$DIAGNOSTIC_GEOM" "$DIAGNOSTIC_EXP" "wideSep_${TAG}"

echo "============================================================"
echo "  Submitted ${#REGRESSION_GEOMS[@]} regression jobs + 1 diagnostic."
echo "  Check with: squeue -u \$USER"
echo "============================================================"
