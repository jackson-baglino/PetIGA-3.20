#!/usr/bin/env bash
# =============================================================================
# submit_permafrost.sh — compute optimal MPI ranks and submit via sbatch
#
# Usage (run from project root):
#   ./scripts/HPC/submit_permafrost.sh <opts_file> [tag] [sbatch_overrides...]
#
# The output folder name is auto-derived from the opts file (ic_type + basename).
# An optional [tag] is appended to the run folder name for disambiguation.
#
# Examples:
#   ./scripts/HPC/submit_permafrost.sh \
#       inputs/tests/test_2D_TouchingGrainPair.opts
#
#   ./scripts/HPC/submit_permafrost.sh \
#       inputs/tests/test_2D_TouchingGrainPair.opts sweep_a \
#       --time=0-12:00:00 --partition=expansion
#
# Any extra arguments after the tag are forwarded verbatim to sbatch and
# can override any of the computed or default resource flags.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Load cost utilities (graceful no-op if missing)
if [[ -f "$SCRIPT_DIR/hpc_cost.sh" ]]; then
    # shellcheck source=scripts/HPC/hpc_cost.sh
    source "$SCRIPT_DIR/hpc_cost.sh"
else
    hpc_cost_pre_submit() { :; }
fi

if [[ "$#" -lt 1 ]]; then
    echo "Usage: $0 <opts_file> [tag] [extra_sbatch_flags...]"
    exit 1
fi

opts_file="$1"
shift 1   # remaining args: optional tag then sbatch overrides

# Second arg is the optional tag (doesn't start with --); rest are sbatch flags
title=""
if [[ "${1:-}" != "" && "${1:-}" != --* ]]; then
    title="$1"
    shift 1
fi

# Resolve opts_file relative to project root if not absolute
[[ "$opts_file" != /* ]] && opts_file="$PROJECT_ROOT/$opts_file"

if [ ! -f "$opts_file" ]; then
    echo "❌ Options file not found: $opts_file"
    exit 1
fi

# ---------------------------------------------------------------------------
# Compute optimal NPROCS from the grid in the opts file.
# Keep in sync with compute_optimal_nprocs() in run_permafrost.sh.
# ---------------------------------------------------------------------------
TARGET_DOFS_PER_CORE=10000
# Tasks per node for the target partition (icelake|skylake|cascadelake).
# Skylake nodes have 32 cores; this is the safe minimum across all three.
NTASKS_PER_NODE=32

Nx=$(awk '$1=="-Nx"{print $2}' "$opts_file" | head -n1)
Ny=$(awk '$1=="-Ny"{print $2}' "$opts_file" | head -n1)
Nz=$(awk '$1=="-Nz"{print $2}' "$opts_file" | head -n1)
[[ -z "${Nx:-}" ]] && Nx=1
[[ -z "${Ny:-}" ]] && Ny=1
[[ -z "${Nz:-}" ]] && Nz=1

total_dofs=$((4 * Nx * Ny * Nz))
NPROCS=$(((total_dofs + TARGET_DOFS_PER_CORE - 1) / TARGET_DOFS_PER_CORE))
(( NPROCS < 1 )) && NPROCS=1

NNODES=$(( (NPROCS + NTASKS_PER_NODE - 1) / NTASKS_PER_NODE ))
(( NNODES < 1 )) && NNODES=1

echo "============================================================"
echo "  Permafrost submission"
echo "  Options  : $opts_file"
echo "  Title    : $title"
echo "------------------------------------------------------------"
echo "  Grid     : Nx=${Nx}, Ny=${Ny}, Nz=${Nz}"
echo "  Total DoFs: 4 × Nx × Ny × Nz = ${total_dofs}"
echo "  Target   : ${TARGET_DOFS_PER_CORE} DoFs/core"
echo "  NPROCS   : ${NPROCS}"
echo "  Nodes    : ${NNODES}  (${NTASKS_PER_NODE} tasks/node)"
echo "============================================================"
hpc_cost_pre_submit "${NPROCS}"

# ---------------------------------------------------------------------------
# Submit — command-line flags override #SBATCH directives in run_permafrost.sh
# ---------------------------------------------------------------------------
job_name="$(basename "$opts_file" .opts)"

# Build run_permafrost.sh args: always pass opts_file; pass title only if set
run_args=("$opts_file")
[[ -n "$title" ]] && run_args+=("$title")

sbatch \
    --job-name="$job_name" \
    --nodes="${NNODES}" \
    --ntasks="${NPROCS}" \
    --ntasks-per-node="${NTASKS_PER_NODE}" \
    "$@" \
    "$SCRIPT_DIR/run_permafrost.sh" \
    "${run_args[@]}"
