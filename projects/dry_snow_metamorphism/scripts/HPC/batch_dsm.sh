#!/bin/bash

# List of temperatures to simulate (you can add more values here)
temperatures=(-80)
humidities=(0.98)

# Set the path to your run script
RUN_SCRIPT="/resnick/groups/rubyfu/jbaglino/PetIGA-3.20/projects/dry_snow_metamorphism/scripts/HPC/run_dsm.sh"

# Set base config names
FILE_BASENAMES=("grainReadFile-70_s1-10")
        #         "grainReadFile-70_s1-11"
        #         "grainReadFile-70_s1-12"
        #         "grainReadFile-70_s1-13"
        #         "grainReadFile-70_s1-14"
        #         "grainReadFile-70_s1-15"
        #         "grainReadFile-70_s1-16"
        #         "grainReadFile-70_s1-17"
        #         "grainReadFile-70_s1-18"
        #         "grainReadFile-70_s1-19"
        #         "grainReadFile-70_s1-20"
        #         "grainReadFile-70_s1-21"
        #  ) 

# Loop through each filename, temperature, and humidity and submit the job
for basename in "${FILE_BASENAMES[@]}"; do
  ENV_FILE="/resnick/groups/rubyfu/jbaglino/PetIGA-3.20/projects/dry_snow_metamorphism/configs/${basename}.env"
  INPUT_FILE="${basename}.dat"

  for temp in "${temperatures[@]}"; do
    for hum in "${humidities[@]}"; do
      hum_tag="${hum/./p}"
      echo "Submitting job for file: $basename, T=${temp}C, rel_hum=${hum}"

      sbatch --job-name="DSM-${basename}-T=${temp}_hum=${hum_tag}" \
             --export=ALL,ENV_FILE_OVERRIDE="$ENV_FILE",temp="$temp",humidity="$hum",inputFile="$INPUT_FILE" \
             "$RUN_SCRIPT"
    done
  done
done