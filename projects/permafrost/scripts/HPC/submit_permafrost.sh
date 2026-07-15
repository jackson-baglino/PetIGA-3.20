#!/usr/bin/env bash
# =============================================================================
# submit_permafrost.sh — compute optimal MPI ranks and submit via sbatch
#
# Usage (run from project root):
#   ./scripts/HPC/submit_permafrost.sh <geometry> <experiment> [tag] \
#       [sbatch_overrides...] [-- extra_permafrost_opts...]
#
#   geometry    Name (without .opts) of a file in inputs/geometry/
#   experiment  Name (without .opts) of a file in inputs/experiment/
#   tag         Optional label appended to the run folder name
#
# Extra arguments after the tag are split on a literal `--`:
#   - before `--` (or if no `--` is given): forwarded verbatim to sbatch,
#     can override any of the computed or default resource flags.
#     Special flag (consumed here, not passed to sbatch):
#       --half-cores   request half the computed MPI ranks — queues faster
#                      on a busy cluster at ~2x wall time.
#   - after `--`: forwarded verbatim to the permafrost executable itself
#     (appended after the three -options_file flags, so they override
#     anything set in solver.opts/geometry/experiment opts files).
#
# Examples:
#   ./scripts/HPC/submit_permafrost.sh 2D_multi_grain_test 2day_T-20_h0.95
#
#   ./scripts/HPC/submit_permafrost.sh 2D_multi_grain_test 2day_T-20_h0.95 p2_run \
#       --time=0-12:00:00 --partition=expansion
#
#   ./scripts/HPC/submit_permafrost.sh 2D_single_bump_two_grains 21day_T-20_h0.95 \
#       d0GT_1e-8 -- -d0_GT 1.0e-8
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

# Split remaining args on a literal `--`: before it are sbatch flags, after
# it are extra options forwarded to the permafrost executable.
# `--half-cores` is intercepted here (not a real sbatch flag): request half
# the computed MPI ranks so the job queues faster on a busy cluster. The
# runner (run_permafrost.sh) clamps its own rank count to the SLURM
# allocation, so the halved request propagates consistently; wall time
# roughly doubles (weak-scaling regime at the 40k DoFs/core target).
sbatch_flags=()
extra_opts=()
sep_seen=0
half_cores=0
for a in "$@"; do
    if [[ "$sep_seen" -eq 0 && "$a" == "--" ]]; then
        sep_seen=1
        continue
    fi
    if [[ "$sep_seen" -eq 0 ]]; then
        if [[ "$a" == "--half-cores" ]]; then
            half_cores=1
            continue
        fi
        sbatch_flags+=("$a")
    else
        extra_opts+=("$a")
    fi
done

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
# 40000 DoFs/rank (raised from 10000, 2026-07-12): PETSc guidance for
# implicit solves is 20k-100k unknowns/rank — below ~20k, reductions and halo
# exchange dominate; and ASM+ILU weakens as subdomain count grows (more ranks
# -> more BiCGStab iterations AND more comms). Empirically the local axisym
# Molaro runs at 108k DoFs/rank did ~7 s/step at 1.3M DoFs, while the old
# target allocated 260 ranks / 9 HPC nodes to a 62-step job whose cost was
# all queue wait. These runs are step-limited, so wall time is nearly flat
# in rank count; the allocation size is what costs. Keep the three copies of
# this constant in sync (Studio/run, HPC/run, HPC/submit, HPC/submit_batch
# -- submit_batch was omitted from this list and silently kept the old
# 10000 until 2026-07-15).
TARGET_DOFS_PER_CORE=40000
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

if [[ "$half_cores" -eq 1 ]]; then
    NPROCS=$(( (NPROCS + 1) / 2 ))
    (( NPROCS < 1 )) && NPROCS=1
    echo "  --half-cores: requesting ${NPROCS} ranks (half the computed optimum)"
fi

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
if [[ "${#extra_opts[@]}" -gt 0 ]]; then
    # title must be positional $3 even if empty, so extra_opts land at $4+
    run_args+=("$title" "${extra_opts[@]}")
    echo "  Extra opts : ${extra_opts[*]}"
elif [[ -n "$title" ]]; then
    run_args+=("$title")
fi

sbatch \
    --job-name="${geom_name}__${exp_name}" \
    --nodes="${NNODES}" \
    --ntasks="${NPROCS}" \
    --ntasks-per-node="${NTASKS_PER_NODE}" \
    "${sbatch_flags[@]}" \
    "$SCRIPT_DIR/run_permafrost.sh" \
    "${run_args[@]}"
