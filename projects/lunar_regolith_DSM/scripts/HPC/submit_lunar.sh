#!/usr/bin/env bash
# =============================================================================
# submit_lunar.sh — compute optimal MPI ranks and submit via sbatch
#
# Usage (run from project root):
#   ./scripts/HPC/submit_lunar.sh <geometry> <experiment> [tag] \
#       [sbatch_overrides...] [-- extra_solver_opts...]
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
#   - after `--`: forwarded verbatim to the lunar_regolith_dsm executable itself
#     (appended after the three -options_file flags, so they override
#     anything set in solver.opts/geometry/experiment opts files).
#
# Examples:
#   ./scripts/HPC/submit_lunar.sh multigrain_2D_L100um_eps0.37um_test base_T-20_h0.95_1d
#
#   ./scripts/HPC/submit_lunar.sh multigrain_2D_L100um_eps0.37um_test base_T-20_h0.95_1d p2_run \
#       --time=0-12:00:00 --partition=expansion
#
#   ./scripts/HPC/submit_lunar.sh twograins_2D_L100um_eps0.37um_2grain base_T-20_h0.95_1d \
#       d0GT_1e-8 -- -d0_GT 1.0e-8
# =============================================================================


# Geometry and experiment .opts live in per-family / per-campaign
# subdirectories. Accept either a bare name or an explicit
# "subdir/name", so the CLI did not change when the files were regrouped.
resolve_opts() {
    local dir="$1" name="$2" hit
    if [ -f "$dir/${name}.opts" ]; then printf '%s' "$dir/${name}.opts"; return 0; fi
    hit=$(find "$dir" -type f -name "${name}.opts" -print 2>/dev/null)
    if [ "$(printf '%s' "$hit" | grep -c .)" -gt 1 ]; then
        echo "❌ Ambiguous name '${name}':" >&2
        printf '%s\n' "$hit" | sed 's|^|     |' >&2
        return 1
    fi
    [ -n "$hit" ] && printf '%s' "$hit" && return 0
    return 1
}

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
# it are extra options forwarded to the lunar_regolith_dsm executable.
# `--half-cores` is intercepted here (not a real sbatch flag): request half
# the computed MPI ranks so the job queues faster on a busy cluster. The
# runner (run_lunar.sh) clamps its own rank count to the SLURM
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

geom_file="$(resolve_opts "$PROJECT_ROOT/inputs/geometry" "$geom_name")" || geom_file="$PROJECT_ROOT/inputs/geometry/${geom_name}.opts"

exp_file="$(resolve_opts "$PROJECT_ROOT/inputs/experiment" "$exp_name")" || exp_file="$PROJECT_ROOT/inputs/experiment/${exp_name}.opts"

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
# run_lunar.sh. Formula: ceil(dof * Nx * Ny * Nz / TARGET_DOFS_PER_CORE)
# ---------------------------------------------------------------------------
# 40000 DoFs/rank (raised from 10000, 2026-07-12): PETSc guidance for
# implicit solves is 20k-100k unknowns/rank — below ~20k, reductions and halo
# exchange dominate; and ASM+ILU weakens as subdomain count grows (more ranks
# -> more BiCGStab iterations AND more comms). Empirically the local axisym
# Molaro runs at 108k DoFs/rank did ~7 s/step at 1.3M DoFs, while the old
# target allocated 260 ranks / 9 HPC nodes to a 62-step job whose cost was
# all queue wait. These runs are step-limited, so wall time is nearly flat
# TARGET_DOFS_PER_CORE and NTASKS_PER_NODE are sourced from
# scripts/lib/alloc.sh (single source of truth; see rationale there).
source "$PROJECT_ROOT/scripts/lib/alloc.sh"

# Read dof from solver.opts (default 3 if absent: ice / temperature / vapor)
dof=$(awk '$1=="-dof"{print $2}' "$SOLVER_OPTS" 2>/dev/null | head -n1)
[[ -z "${dof:-}" ]] && dof=3

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

read -r NNODES TASKS_PER_NODE <<< "$(plan_alloc "$NPROCS")"

echo "============================================================"
echo "  Lunar regolith DSM submission"
echo "  Geometry   : ${geom_name}"
echo "  Experiment : ${exp_name}"
echo "  Title      : ${title}"
echo "------------------------------------------------------------"
echo "  Grid source: ${grid_src}"
echo "  Grid       : Nx=${Nx}, Ny=${Ny}, Nz=${Nz}  (dof=${dof})"
echo "  Total DoFs : ${dof} × ${Nx} × ${Ny} × ${Nz} = ${total_dofs}"
echo "  Target     : ${TARGET_DOFS_PER_CORE} DoFs/core"
echo "  NPROCS     : ${NPROCS}"
echo "  Nodes      : ${NNODES}  (${TASKS_PER_NODE} tasks/node)"
echo "============================================================"
hpc_cost_pre_submit "${NPROCS}"

# ---------------------------------------------------------------------------
# Submit — --ntasks/--nodes override the #SBATCH defaults in run_lunar.sh
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
    --ntasks-per-node="${TASKS_PER_NODE}" \
    "${sbatch_flags[@]}" \
    "$SCRIPT_DIR/run_lunar.sh" \
    "${run_args[@]}"
