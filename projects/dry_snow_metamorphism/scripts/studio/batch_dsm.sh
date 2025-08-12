#!/bin/zsh

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

# Path to run_dsm.sh
RUN_SCRIPT="/Users/jacksonbaglino/PetIGA-3.20/projects/dry_snow_metamorphism/scripts/studio/run_dsm.sh"

start_time=$(date +%s)

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
      zsh "$RUN_SCRIPT"
      
      echo "Finished simulation with T=${temp}, RH=${hum}, file=${input_file}"
      echo "=================================================="
      echo ""
    done
  done
done

end_time=$(date +%s)
total_time=$((end_time - start_time))
echo "Total batch runtime: ${total_time} seconds"
