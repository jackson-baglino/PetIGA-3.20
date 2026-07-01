#!/usr/bin/env bash
# =============================================================================
# submit_mob_sweep.sh — Fan out one SLURM job per mob_sub value, all writing
# into a single shared parent folder under $SCRATCH.
#
# Usage (from project root):
#   ./scripts/HPC/submit_mob_sweep.sh [experiment] [tag] [-- sbatch_overrides...]
#
#   experiment   Experiment opts name (default: 30day_T-5_h1.00_GTphys)
#   tag          Optional label appended to the batch folder name
#
# Examples:
#   ./scripts/HPC/submit_mob_sweep.sh
#   ./scripts/HPC/submit_mob_sweep.sh 5yr_T-5_h1.00_GTphys longrun
#   ./scripts/HPC/submit_mob_sweep.sh 30day_T-5_h1.00_GTphys -- --time=0-06:00:00
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
RUN_SCRIPT="$SCRIPT_DIR/run_permafrost.sh"
INPUTS_DIR="$PROJECT_ROOT/inputs"
GEOMETRY_DIR="$INPUTS_DIR/geometry"
EXPERIMENT_DIR="$INPUTS_DIR/experiment"

GEOM="2D_two_ice_grains_boundary"
EXP="${1:-30day_T-5_h1.00_GTphys}"; shift || true
TAG="${1:-mob_sweep}"; [[ "${1:-}" == "--" ]] && TAG="mob_sweep" || shift || true

# Remaining args after -- go to sbatch
sbatch_extra=()
if [[ "${1:-}" == "--" ]]; then
    shift
    sbatch_extra=("$@")
fi

TARGET_DOFS_PER_CORE=10000
MAX_TASKS_PER_NODE=32

# mob_sub values to sweep (baseline → 100× reduction)
MOB_VALUES=(
    6.5404e-10   # K&P baseline (safety=0.5); reference
    1.3e-10      #  5× reduction
    6.5e-11      # 10× reduction
    1.3e-11      # 50× reduction
    6.5e-12      # 100× reduction
)

# ---------------------------------------------------------------------------
# Compute allocation from geometry opts
# ---------------------------------------------------------------------------
GEOM_OPTS="$GEOMETRY_DIR/${GEOM}.opts"
Nx=$(awk '$1=="-Nx"{print $2}' "$GEOM_OPTS" | head -n1); Nx=${Nx:-1}
Ny=$(awk '$1=="-Ny"{print $2}' "$GEOM_OPTS" | head -n1); Ny=${Ny:-1}
Nz=$(awk '$1=="-Nz"{print $2}' "$GEOM_OPTS" | head -n1); Nz=${Nz:-1}
total_dofs=$(( 4 * Nx * Ny * Nz ))
nprocs=$(( (total_dofs + TARGET_DOFS_PER_CORE - 1) / TARGET_DOFS_PER_CORE ))
(( nprocs < 1 )) && nprocs=1
tasks_per_node=$nprocs
(( tasks_per_node > MAX_TASKS_PER_NODE )) && tasks_per_node=$MAX_TASKS_PER_NODE
nnodes=$(( (nprocs + tasks_per_node - 1) / tasks_per_node ))
(( nnodes < 1 )) && nnodes=1

# ---------------------------------------------------------------------------
# Build once on the submission host
# ---------------------------------------------------------------------------
echo "--- Compiling ---"
cd "$PROJECT_ROOT"
make clean && make all
echo "✅ Build complete."

# ---------------------------------------------------------------------------
# Create shared parent folder
# ---------------------------------------------------------------------------
TS=$(date +%Y-%m-%d__%H.%M.%S)
if [[ -d "${SCRATCH:-}" ]]; then
    BATCH_PARENT="$SCRATCH/permafrost/${GEOM}/${TS}_${TAG}"
else
    BATCH_PARENT="$PROJECT_ROOT/scratch/${GEOM}/${TS}_${TAG}"
fi
mkdir -p "$BATCH_PARENT"

echo ""
echo "============================================================"
echo "  mob_sub sweep — $GEOM × $EXP"
echo "  Batch dir   : $BATCH_PARENT"
echo "  DOFs        : $total_dofs  →  nprocs=$nprocs  nodes=$nnodes"
echo "  Values      : ${MOB_VALUES[*]}"
echo "============================================================"

# ---------------------------------------------------------------------------
# Fan out: one sbatch job per mob_sub value
# ---------------------------------------------------------------------------
for mob in "${MOB_VALUES[@]}"; do
    job_tag="mob_${mob}"
    sbatch \
        --job-name="pf_mob_${mob}" \
        --nodes="$nnodes" \
        --ntasks="$nprocs" \
        --ntasks-per-node="$tasks_per_node" \
        --export=ALL,SKIP_COMPILE=1,BATCH_OUT_DIR="$BATCH_PARENT" \
        "${sbatch_extra[@]}" \
        "$RUN_SCRIPT" "$GEOM" "$EXP" "$job_tag" -mob_sub "$mob"

    echo "  → submitted mob_sub=$mob"
done

echo ""
echo "  All jobs submitted. Check with: squeue -u \$USER"
echo "  Results will appear in: $BATCH_PARENT"
