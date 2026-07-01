#!/usr/bin/env bash
# =============================================================================
# submit_mob_sweep.sh — Submit the mob_sub sweep as a single SLURM job
#
# Usage (from project root):
#   ./scripts/HPC/submit_mob_sweep.sh [experiment] [sbatch_overrides...]
#
#   experiment   Experiment opts name (default: 30day_T-5_h1.00_GTphys)
#
# Examples:
#   ./scripts/HPC/submit_mob_sweep.sh
#   ./scripts/HPC/submit_mob_sweep.sh 5yr_T-5_h1.00_GTphys --time=0-72:00:00
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXP="${1:-30day_T-5_h1.00_GTphys}"
shift || true   # remaining args forwarded to sbatch

sbatch "$@" "$SCRIPT_DIR/run_mob_sweep.sh" "$EXP"
