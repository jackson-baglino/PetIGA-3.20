#!/bin/bash

PARENT_DIR="/Users/jacksonbaglino/SimulationResults/HPC_results/dry_snow_metamorphism/30G_microstructures/jbaglino"

for dir in "$PARENT_DIR"/*; do
    if [ -d "$dir" ]; then
        echo "Processing directory: $dir"
        unset ENV_FILE
        INIT_DIR="$dir" ./scripts/run_effective_k_ice_homog.sh
    fi
done
