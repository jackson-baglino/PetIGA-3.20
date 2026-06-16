#!/usr/bin/env bash
# =============================================================================
# submit_permafrost.sh — compute optimal MPI ranks and submit via sbatch
#
# Usage (run from project root):
#   ./scripts/HPC/submit_permafrost.sh <geometry> <experiment> [tag] [sbatch_overrides...]
#
#   geometry    Name (without .opts) of a file in inputs/geometry/
#   experiment  Name (without .opts) of a file in inputs/experiment/
#   tag         Optional label appended to the run folder name
#
# Any extra arguments after the tag (starting with --) are forwarded verbatim
# to sbatch and can override any of the computed or default resource flags.
#
# Examples:
#   ./scripts/HPC/submit_permafrost.sh 2D_multi_grain_test 2day_T-20_h0.95
#
#   ./scripts/HPC/submit_permafrost.sh 2D_multi_grain_test 2day_T-20_h0.95 p2_run \
#       --time=0-12:00:00 --partition=expansion
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

SOLVER_OPTS="$PROJECT_ROOT/inputs/solver.opts"

# Load cost utilities (graceful no-op if missing)
if [[ -f "$SCRIPT_DIR/hpc_cost.sh" ]]; then
    source "$SCRIPT_DIR/hpc_cost.sh"
else
    hpc_cost_pre_submit() { :; }
fi

if [[ "$#" -lt 2 ]]; then
    echo "Usage: $0 <geometry> <experiment> [tag] [extra_sbatch_flags...]"
    echo "  geometry   : name (without .opts) from inputs/geometry/"
    echo "  experiment : name (without .opts) from inputs/experiment/"
    echo "  tag        : optional run-folder label"
    exit 1
fi

geom_name="$1"
exp_name="$2"
shift 2

# Third arg is the optional tag (doesn't start with --); rest are sbatch flags
title=""
if [[ "${1:-}" != "" && "${1:-}" != --* ]]; then
    title="$1"
    shift 1
fi

geom_file="$PROJECT_ROOT/inputs/geometry/${geom_name}.opts"
exp_file="$PROJECT_ROOT/inputs/experiment/${exp_name}.opts"

if [ ! -f "$geom_file" ]; then
    echo "❌ Geometry opts not found: $geom_file"
    exit 1
fi
if [ ! -f "$exp_file" ]; then
    echo "❌ Experiment opts not found: $exp_file"
    exit 1
fi

# ---------------------------------------------------------------------------
# Compute optimal NPROCS — kept in sync with compute_optimal_nprocs() in
# run_permafrost.sh. Formula: ceil(dof * Nx * Ny * Nz / TARGET_DOFS_PER_CORE)
# ---------------------------------------------------------------------------
TARGET_DOFS_PER_CORE=10000
NTASKS_PER_NODE=32   # safe minimum across icelake|skylake|cascadelake

# Read dof from solver.opts (default 4 if absent)
dof=$(awk '$1=="-dof"{print $2}' "$SOLVER_OPTS" 2>/dev/null | head -n1)
[[ -z "${dof:-}" ]] && dof=4

# For -geom_file meshes the DOF grid is in the "# DOF_GRID: nx ny [nz]" comment
if grep -q "^-geom_file" "$geom_file"; then
    read -r Nx Ny Nz <<< "$(awk '$1=="#" && $2=="DOF_GRID:"{print $3, $4, $5}' "$geom_file" | head -n1)"
    grid_src="DOF_GRID comment"
else
    Nx=$(awk '$1=="-Nx"{print $2}' "$geom_file" | head -n1)
    Ny=$(awk '$1=="-Ny"{print $2}' "$geom_file" | head -n1)
    Nz=$(awk '$1=="-Nz"{print $2}' "$geom_file" | head -n1)
    grid_src="-Nx/-Ny/-Nz flags"
fi
[[ -z "${Nx:-}" ]] && Nx=1
[[ -z "${Ny:-}" ]] && Ny=1
[[ -z "${Nz:-}" ]] && Nz=1

total_dofs=$((dof * Nx * Ny * Nz))
NPROCS=$(((total_dofs + TARGET_DOFS_PER_CORE - 1) / TARGET_DOFS_PER_CORE))
(( NPROCS < 1 )) && NPROCS=1

NNODES=$(( (NPROCS + NTASKS_PER_NODE - 1) / NTASKS_PER_NODE ))
(( NNODES < 1 )) && NNODES=1

echo "============================================================"
echo "  Permafrost submission"
echo "  Geometry   : ${geom_name}"
echo "  Experiment : ${exp_name}"
echo "  Title      : ${title}"
echo "------------------------------------------------------------"
echo "  Grid source: ${grid_src}"
echo "  Grid       : Nx=${Nx}, Ny=${Ny}, Nz=${Nz}  (dof=${dof})"
echo "  Total DoFs : ${dof} × ${Nx} × ${Ny} × ${Nz} = ${total_dofs}"
echo "  Target     : ${TARGET_DOFS_PER_CORE} DoFs/core"
echo "  NPROCS     : ${NPROCS}"
echo "  Nodes      : ${NNODES}  (${NTASKS_PER_NODE} tasks/node)"
echo "============================================================"
hpc_cost_pre_submit "${NPROCS}"

# ---------------------------------------------------------------------------
# Submit — --ntasks/--nodes override the #SBATCH defaults in run_permafrost.sh
# ---------------------------------------------------------------------------
run_args=("$geom_name" "$exp_name")
[[ -n "$title" ]] && run_args+=("$title")

sbatch \
    --job-name="${geom_name}__${exp_name}" \
    --nodes="${NNODES}" \
    --ntasks="${NPROCS}" \
    --ntasks-per-node="${NTASKS_PER_NODE}" \
    "$@" \
    "$SCRIPT_DIR/run_permafrost.sh" \
    "${run_args[@]}"
