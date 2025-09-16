#!/bin/bash

###############################################################################
# Script: batch_effective_k_ice_homog.sh
# Purpose:
#   Sequentially run the master run_effective_k_ice_homog.sh script over all
#   input directories found under PARENT_DIR. Each subdirectory should contain
#   initial condition files and (optionally) a .env file.
#
# Usage:
#   ./scripts/batch_effective_k_ice_homog.sh
#
# Notes:
#   - This script loops over all subdirectories in PARENT_DIR.
#   - It unsets ENV_FILE for each run so the master script will auto-discover
#     the correct .env in the INIT_DIR.
#   - Change PARENT_DIR to point to your parent input directory.
###############################################################################

# PARENT_DIR="/Users/jacksonbaglino/SimulationResults/HPC_results/dry_snow_metamorphism/DSM-grainReadFile-70/"
PARENT_DIR="/Users/jacksonbaglino/SimulationResults/HPC_results/dry_snow_metamorphism/DSMgrains__phi=0.24__Lxmm=1__Lymm=1__seed=22"

for dir in "$PARENT_DIR"/*; do
    if [ -d "$dir" ]; then
        echo "Processing directory: $dir"
        unset ENV_FILE
        ENV_FILE="" \
        INIT_DIR="$dir" \
          ./scripts/run_effective_k_ice_homog.sh
    fi
done

#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'
shopt -s nullglob   # so empty globs don't expand to literal *

###############################################################################
# Script: batch_effective_k_ice_homog.sh
# Purpose:
#   Run run_effective_k_ice_homog.sh over all input subdirectories of PARENT_DIR.
#   If PARENT_DIR has no subdirectories, run it directly as a single case.
#
# Usage:
#   ./scripts/batch_effective_k_ice_homog.sh
#
# Notes:
#   - Do NOT pass ENV_FILE here; the runner will auto-discover a grains.env/.env
#     inside INIT_DIR.
#   - Update PARENT_DIR below to either a parent folder containing many cases,
#     or a single case folder.
###############################################################################

# Example parent with many cases:
# PARENT_DIR="/Users/jacksonbaglino/SimulationResults/HPC_results/dry_snow_metamorphism/DSM-grainReadFile-70"
# Current setting (single case folder OK too):
PARENT_DIR="/Users/jacksonbaglino/SimulationResults/HPC_results/dry_snow_metamorphism/DSMgrains__phi=0.24__Lxmm=1__Lymm=1__seed=22"

# Path to runner relative to this script
RUNNER="$(dirname "$0")/run_effective_k_ice_homog.sh"

if [[ ! -x "$RUNNER" ]]; then
  echo "ERROR: Runner not found or not executable: $RUNNER" >&2
  exit 1
fi

processed=0

# Case 1: Iterate all immediate subdirectories of PARENT_DIR (if any)
for dir in "$PARENT_DIR"/*/; do
  [[ -d "$dir" ]] || continue
  ((processed++))
  echo "────────────────────────────────────────────────────────"
  echo "Processing directory: $dir"
  INIT_DIR="$dir" "$RUNNER"
done

# Case 2: If no subdirectories were processed, treat PARENT_DIR as a single case
if (( processed == 0 )); then
  if [[ -d "$PARENT_DIR" ]]; then
    echo "────────────────────────────────────────────────────────"
    echo "No subdirectories found; running single case: $PARENT_DIR"
    INIT_DIR="$PARENT_DIR" "$RUNNER"
  else
    echo "ERROR: PARENT_DIR is not a directory: $PARENT_DIR" >&2
    exit 1
  fi
fi