###############################################################################
# Script: batch_dsm.sh (HPC)
# Model: Dry Snow Metamorphism (DSM)
# Purpose:
#   Submit a sweep of DSM simulations to SLURM, varying temperature(s),
#   humidity(ies), and input file(s). Calls run_dsm.sh for each combination.
#
# Usage:
#   sbatch batch_dsm.sh
#
# Notes:
#   - Overrides temp, humidity, and inputFile by exporting them to SLURM.
#   - RUN_LABEL is generated inside run_dsm.sh using compact tags (2-digit temp/RH).
#   - Adjust FILE_BASENAMES, temperatures, humidities arrays below to choose sweep set.
###############################################################################


# List of temperatures to simulate (you can add more values here)
temperatures=(-80)
humidities=(0.98)

# Set the path to your run script
RUN_SCRIPT="/resnick/groups/rubyfu/jbaglino/PetIGA-3.20/projects/dry_snow_metamorphism/scripts/HPC/run_dsm.sh"

INPUT_DIR="/resnick/groups/rubyfu/jbaglino/PetIGA-3.20/projects/dry_snow_metamorphism/inputs/porespy2"
CONFIG_DIR="/resnick/groups/rubyfu/jbaglino/PetIGA-3.20/projects/dry_snow_metamorphism/configs/porespy2"

shopt -s nullglob
dat_files=("$INPUT_DIR"/*.dat)
shopt -u nullglob

if [ ${#dat_files[@]} -eq 0 ]; then
  echo "No .dat files found in $INPUT_DIR"
  exit 1
fi

for dat_file in "${dat_files[@]}"; do
  basename=$(basename "$dat_file" .dat)
  env_file="$CONFIG_DIR/$basename.env"
  if [ ! -f "$env_file" ]; then
    echo "Warning: Env file $env_file not found for $dat_file, skipping."
    continue
  fi
  ENV_FILE="$env_file"
  INPUT_FILE="$dat_file"

  echo "Submitting job for basename: $basename, input file: $INPUT_FILE, env file: $ENV_FILE"

  for temp in "${temperatures[@]}"; do
    for hum in "${humidities[@]}"; do
      temp_int=$(printf "%.0f" "$temp")
      if [[ "$temp_int" == -* ]]; then
        temp_tag=${temp_int:0:3}
      else
        temp_tag=${temp_int:0:2}
      fi
      hum_int=$(awk "BEGIN{printf \"%d\", $hum*100}")
      hum_tag=$(printf "%02d" "$hum_int"); hum_tag=${hum_tag:0:2}

      echo "Submitting job for file: $basename, T=${temp_tag}C, RH=${hum_tag}%"

      sbatch --job-name="DSM-${basename}_Tm${temp_tag}_hum${hum_tag}" \
             --export=ALL,ENV_FILE_OVERRIDE="$ENV_FILE",temp="$temp",humidity="$hum",inputFile="$INPUT_FILE" \
             "$RUN_SCRIPT"
    done
  done
done