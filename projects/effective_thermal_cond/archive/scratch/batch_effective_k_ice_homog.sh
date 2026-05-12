#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'
shopt -s nullglob   # empty globs expand to nothing

# Helper: source an env file (if present) and export expected variables
load_grains_env() {
  local envfile="$1"
  if [[ -f "$envfile" ]]; then
    echo "Loading parameters from: $envfile"
    # Save current values to restore if file omits some
    local _Nx=${Nx-} _Ny=${Ny-} _Nz=${Nz-} _Lx=${Lx-} _Ly=${Ly-} _Lz=${Lz-} _eps=${eps-}
    # shellcheck disable=SC1090
    set +u; set -a
    source "$envfile"
    set +a; set -u
    # Restore any previously set values that were not defined by the file
    [[ -n ${_Nx-}  ]] && export Nx="${Nx-$_Nx}"
    [[ -n ${_Ny-}  ]] && export Ny="${Ny-$_Ny}"
    [[ -n ${_Nz-}  ]] && export Nz="${Nz-$_Nz}"
    [[ -n ${_Lx-}  ]] && export Lx="${Lx-$_Lx}"
    [[ -n ${_Ly-}  ]] && export Ly="${Ly-$_Ly}"
    [[ -n ${_Lz-}  ]] && export Lz="${Lz-$_Lz}"
    [[ -n ${_eps-} ]] && export eps="${eps-$_eps}"
    echo "  Lx=${Lx:-?} Ly=${Ly:-?} Lz=${Lz:-?} | Nx=${Nx:-?} Ny=${Ny:-?} Nz=${Nz:-?} | eps=${eps:-?}"
    return 0
  fi
  return 1
}

###############################################################################
# Script: batch_effective_k_ice_homog.sh
# Purpose:
#   Run run_effective_k_ice_homog.sh over all input subdirectories of PARENT_DIR.
#   If PARENT_DIR has no subdirectories, run it directly as a single case.
#
# Usage:
#   ./scripts/batch_effective_k_ice_homog.sh [PARENT_DIR]
#
# Notes:
#   - Do NOT pass ENV_FILE here; the runner will auto-discover a grains.env/.env
#     inside INIT_DIR. We also explicitly unset ENV_FILE for safety.
#   - If an argument is given, it overrides the default PARENT_DIR below.
###############################################################################

# Default parent directory (override via CLI arg)
DEFAULT_PARENT_DIR="/Users/jacksonbaglino/SimulationResults/HPC_results/dry_snow_metamorphism/DSMgrains__phi=0.24__Lxmm=1__Lymm=1__seed=22"
PARENT_DIR="${1:-$DEFAULT_PARENT_DIR}"

# Path to runner relative to this script
RUNNER="$(dirname "$0")/run_effective_k_ice_homog.sh"

if [[ ! -x "$RUNNER" ]]; then
  echo "ERROR: Runner not found or not executable: $RUNNER" >&2
  exit 1
fi

if [[ ! -d "$PARENT_DIR" ]]; then
  echo "ERROR: PARENT_DIR is not a directory: $PARENT_DIR" >&2
  exit 1
fi

echo "PARENT_DIR: $PARENT_DIR"
processed=0       

#
# Case 1: Iterate all immediate subdirectories of PARENT_DIR (if any)
for dir in "$PARENT_DIR"/*/; do
  [[ -d "$dir" ]] || continue
  ((processed++))
  echo "────────────────────────────────────────────────────────"
  echo "Processing directory: $dir"

  # Load per-case parameters from grains.env if present
  load_grains_env "$dir/grains.env" || echo "WARNING: grains.env not found in $dir — proceeding without these variables"

  # Detect per-case input subdirectories and run per input if present
  input_dirs=()
  if [[ -d "$dir/inputs" ]]; then
    for id in "$dir/inputs"/*/; do [[ -d "$id" ]] && input_dirs+=("$id"); done
  fi
  # Fallback: any subdirs matching input*/
  if (( ${#input_dirs[@]} == 0 )); then
    for id in "$dir"/input*/; do [[ -d "$id" ]] && input_dirs+=("$id"); done
  fi

  if (( ${#input_dirs[@]} > 0 )); then
    echo "Found ${#input_dirs[@]} input director(ies) under $dir"
    for in_dir in "${input_dirs[@]}"; do
      echo "  → Input: $in_dir"
      # If the input dir has its own grains.env, load it (overrides case-level)
      load_grains_env "$in_dir/grains.env" || true
      env -u ENV_FILE INIT_DIR="$dir" INPUT_DIR="$in_dir" \
          Lx="${Lx:-}" Ly="${Ly:-}" Lz="${Lz:-}" \
          Nx="${Nx:-}" Ny="${Ny:-}" Nz="${Nz:-}" \
          eps="${eps:-}" \
          "$RUNNER"
    done
  else
    # No input subdirectories; run once for this case directory
    env -u ENV_FILE INIT_DIR="$dir" \
        Lx="${Lx:-}" Ly="${Ly:-}" Lz="${Lz:-}" \
        Nx="${Nx:-}" Ny="${Ny:-}" Nz="${Nz:-}" \
        eps="${eps:-}" \
        "$RUNNER"
  fi
done

# Case 2: If no subdirectories were processed, treat PARENT_DIR as a single case
if (( processed == 0 )); then
  echo "────────────────────────────────────────────────────────"
  echo "No subdirectories found; running single case: $PARENT_DIR"
  env -u ENV_FILE INIT_DIR="$PARENT_DIR" "$RUNNER"
fi