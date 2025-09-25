#!/usr/bin/env bash
set -euo pipefail             # Exit on error, undefined variable, or error in a pipeline
IFS=$'\n\t'                   # Set internal field separator to newline and tab
shopt -s nullglob             # Allow for loops over empty globs

# Helper function: source a .env file if it exists
load_grains_env() {
  echo "Loading environment variables from grains.env"

  # Check if grains.env exists in the current directory
  if [[ -f "grains.env" ]]; then
    # Source the grains.env file to load environment variables
    set +u; source "grains.env"; set -u
    echo "Environment variables loaded from grains.env"

    # Check that all required variables are set
    if [[ -z "${Nx:-}" ]]; then
      echo "Error: Nx is not set in grains.env"
      exit 1
    elif [[ -z "${Ny:-}" ]]; then
      echo "Error: Ny is not set in grains.env"
      exit 1
    elif [[ -z "${Nz:-}" ]]; then
      echo "Error: Nz is not set in grains.env"
      exit 1
    elif [[ -z "${Lx:-}" ]]; then
      echo "Error: Lx is not set in grains.env"
      exit 1
    elif [[ -z "${Ly:-}" ]]; then
      echo "Error: Ly is not set in grains.env"
      exit 1
    elif [[ -z "${Lz:-}" ]]; then
      echo "Error: Lz is not set in grains.env"
      exit 1
    elif [[ -z "${eps:-}" ]]; then
      echo "Error: eps is not set in grains.env"
      exit 1
    fi

    # Export required variables
    export Nx Ny Nz Lx Ly Lz eps
    return 0
  fi
  return 1  # grains.env not found
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

# Default parent directory containing subdirectories to process
PARENT_DIR="/Users/jacksonbaglino/SimulationResults/HPC_results/dry_snow_metamorphism/all_results"
OUT_ROOT="${OUT_ROOT:-/Users/jacksonbaglino/SimulationResults/effective_thermal_cond/scratch}"


# Check that the parent directory exists
if [[ ! -d "$PARENT_DIR" ]]; then
  echo "Error: Parent directory '$PARENT_DIR' does not exist."
  exit 1
fi 

echo "Processing parent directory: $PARENT_DIR"
processed=0           # Count of processed directories


# Iterate all immediate subdirectories of PARENT_DIR (if any)
for dir in "$PARENT_DIR"/*; do
  if [[ -d "$dir" ]]; then
    echo "Processing directory: $dir"

    # Extract the base name of the input directory (strip trailing slash)
    base="$(basename "${dir%/}")"
    # Skip if already processed (k_eff.csv present in corresponding output subdirectory)
    if [[ -f "$OUT_ROOT/$base/k_eff.csv" ]]; then
      echo "Skipping $dir (already has $OUT_ROOT/$base/k_eff.csv)"
      continue
    fi

    # Unset ENV_FILE to ensure the runner auto-discovers the environment file
    unset ENV_FILE

    # Load environment variables from grains.env; skip if missing
    pushd "$dir" >/dev/null
    if ! load_grains_env; then
      echo "Skipping $dir (no grains.env)"
      popd >/dev/null
      continue   # skip this directory and move on to the next one
    fi

    # Return to original directory
    popd >/dev/null

    # Echo loaded variables for verification
    echo "Loaded environment variables:"
    echo "  Nx: $Nx"
    echo "  Ny: $Ny"
    echo "  Nz: $Nz"
    echo "  Lx: $Lx"
    echo "  Ly: $Ly"
    echo "  Lz: $Lz"
    echo "  eps: $eps"

    # Skip if the input folder has no sol_*.dat files
    if ! compgen -G "$dir/sol_*.dat" > /dev/null; then
      echo "Skipping $dir (no sol_*.dat files found)"
      continue
    fi

    # Run the effective k ice homogenization script
    ./scripts/run_effective_k_ice_homog.sh "$dir"

    # Increment processed count
    ((processed++))
  fi
done

# If no subdirectories were processed, run the parent directory itself
if (( processed == 0 )); then
  echo "No subdirectories found. Processing parent directory directly: $PARENT_DIR"
  unset ENV_FILE
  pushd "$PARENT_DIR" >/dev/null
  if ! load_grains_env; then
    echo "Error: grains.env not found in $PARENT_DIR"
    popd >/dev/null
    exit 1
  fi
  
  # Return to original directory
  popd >/dev/null

  # Echo loaded variables for verification
  echo "Loaded environment variables:"
  echo "  Nx: $Nx"
  echo "  Ny: $Ny"
  echo "  Nz: $Nz"
  echo "  Lx: $Lx"
  echo "  Ly: $Ly"
  echo "  Lz: $Lz"
  echo "  eps: $eps"
  ./scripts/run_effective_k_ice_homog.sh "$PARENT_DIR"
fi