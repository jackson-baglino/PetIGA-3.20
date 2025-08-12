#!/bin/bash

# List of temperatures to simulate (you can add more values here)
temperatures=(-80 -75 -70 -65 -60 -55 -50 -45 -40 -35 -30 -25 -20 -15 -10)

# Set the path to your run script
RUN_SCRIPT="/resnick/groups/rubyfu/jbaglino/PetIGA-3.20/projects/dry_snow_metamorphism/scripts/HPC/run_dsm.sh"

# Set a base config name
CONFIG_BASENAME="grainReadFile-35_s1-10"

# Base env file path
ENV_FILE="/resnick/groups/rubyfu/jbaglino/PetIGA-3.20/projects/dry_snow_metamorphism/configs/${CONFIG_BASENAME}.env"

# Loop through each temperature and submit the job
for temp in "${temperatures[@]}"; do
  echo "Submitting job for temperature: $temp°C"

  # Create a temporary .env file specific to this temperature
  # ENV_FILE="/resnick/groups/rubyfu/jbaglino/PetIGA-3.20/projects/dry_snow_metamorphism/configs/${CONFIG_BASENAME}_T${temp}.env"
  # cp "/resnick/groups/rubyfu/jbaglino/PetIGA-3.20/projects/dry_snow_metamorphism/configs/${CONFIG_BASENAME}.env" "$ENV_FILE"

  # Replace temperature value in the env file
  # sed -i "s/^temp=.*/temp=${temp}/" "$ENV_FILE"

  # Submit the job with SLURM, overriding SLURM_JOB_NAME to embed the temperature
  sbatch --job-name="DSM-T=${temp}_hum=0.50" \
         --export=ALL,ENV_FILE_OVERRIDE="$ENV_FILE",temp="$temp" \
         "$RUN_SCRIPT"
done