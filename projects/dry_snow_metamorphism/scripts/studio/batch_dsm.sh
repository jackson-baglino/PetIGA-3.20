#!/bin/zsh
###############################################################################
# Script: batch_dsm.sh
# Model: Dry Snow Metamorphism (DSM)
# Purpose:
#   Sweep over temperature(s), humidity(ies), and input file(s), calling the
#   single-run script (run_dsm.sh) sequentially — one job at a time.
#
# Usage:
#   ./scripts/studio/batch_dsm.sh
#
# Precedence of parameters:
#   - This script exports env vars (temp, humidity, filename) before each call.
#   - The run script should use these exported values if set; otherwise it falls
#     back to its own internal defaults. (We will update run_dsm.sh accordingly.)
#
# Notes:
#   - Adjust the arrays below to choose the sweep set.
#   - Optionally export NUM_PROCS before running to override MPI ranks.
#   - All existing comments remain unchanged by design.
###############################################################################

# ============================
# Batch run configuration
# ============================

# Define arrays of values
temps=(-12 -14 -16 -18 -20 -22 -24 -26 -28 -30)
humidities=(0.98 0.90 0.50 0.05)
input_files=(
  "grainReadFile-2G_Molaro_0p25R1_HIGHRES.dat"
  "grainReadFile-2G_Molaro_0p25R1.dat"
)

# Absolute path to the single-run driver script
# Path to run_dsm.sh
RUN_SCRIPT="/Users/jacksonbaglino/PetIGA-3.20/projects/dry_snow_metamorphism/scripts/studio/run_dsm.sh"

if [[ ! -f "$RUN_SCRIPT" ]]; then
  echo "[ERROR] run script not found: $RUN_SCRIPT"
  exit 1
fi

start_time=$(date +%s)
echo "Starting DSM batch: ${#temps[@]} temps × ${#humidities[@]} humidities × ${#input_files[@]} files"

# Loop over combinations
for temp in "${temps[@]}"; do
  for hum in "${humidities[@]}"; do
    for input_file in "${input_files[@]}"; do
      echo "=================================================="
      echo "Starting simulation with T=${temp}, RH=${hum}, file=${input_file}"
      
      # Export variables to be picked up by run_dsm.sh
      export temp=$temp
      export humidity=$hum
      export filename=$input_file

      # Call the main run script
      per_run_start=$(date +%s)
      zsh "$RUN_SCRIPT"
      per_run_end=$(date +%s)
      echo "Run elapsed: $((per_run_end - per_run_start)) s"
      
      echo "Finished simulation with T=${temp}, RH=${hum}, file=${input_file}"
      echo "=================================================="
      echo ""
    done
  done
done

end_time=$(date +%s)
total_time=$((end_time - start_time))
echo "Total batch runtime: ${total_time} seconds"
