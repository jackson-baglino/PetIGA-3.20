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
PARENT_DIR="/Users/jacksonbaglino/SimulationResults/dry_snow_metamorphism/scratch"

for dir in "$PARENT_DIR"/*; do
    if [ -d "$dir" ]; then
        echo "Processing directory: $dir"
        unset ENV_FILE
        ENV_FILE="" \
        INIT_DIR="$dir" \
          ./scripts/run_effective_k_ice_homog.sh
    fi
done
