#!/usr/bin/env bash
# =============================================================================
# submit_permafrost.sh — compute optimal MPI ranks and submit via sbatch
#
# Usage (run from project root):
#   ./scripts/HPC/submit_permafrost.sh <opts_file> <title_prefix> [sbatch_overrides...]
#
# Examples:
#   ./scripts/HPC/submit_permafrost.sh \
#       inputs/tests/test3_EnclosedGrainPair.opts EnclosedGrainPair_
#
#   ./scripts/HPC/submit_permafrost.sh \
#       inputs/tests/test3_EnclosedGrainPair.opts EnclosedGrainPair_ \
#       --time=0-12:00:00 --partition=expansion
#
# Any extra arguments after the title are forwarded verbatim to sbatch and
# can override any of the computed or default resource flags.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

if [[ "$#" -lt 2 ]]; then
    echo "Usage: $0 <opts_file> <title_prefix> [extra_sbatch_flags...]"
    exit 1
fi

opts_file="$1"
title="$2"
shift 2   # remaining args passed verbatim to sbatch

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

# ---------------------------------------------------------------------------
# Submit — command-line flags override #SBATCH directives in run_permafrost.sh
# ---------------------------------------------------------------------------
sbatch \
    --nodes="${NNODES}" \
    --ntasks="${NPROCS}" \
    --ntasks-per-node="${NTASKS_PER_NODE}" \
    "$@" \
    "$SCRIPT_DIR/run_permafrost.sh" \
    "$opts_file" \
    "$title"
