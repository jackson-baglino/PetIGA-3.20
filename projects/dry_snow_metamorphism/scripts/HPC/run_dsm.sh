#!/bin/bash
#SBATCH -J DSM-T=-80_hum=0.90
#SBATCH -A rubyfu
#SBATCH -t 5-00:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=50
#SBATCH --cpus-per-task=1
#SBATCH -o "output_files/%x.o%j"
#SBATCH -e "output_files/%x.e%j"
#SBATCH --partition=expansion
#SBATCH --mem-per-cpu=1G
#SBATCH --mail-user=jbaglino@caltech.edu
#SBATCH --mail-type=END,FAIL,TIME_LIMIT

##############################################
# USER-MODIFIABLE SIMULATION SETTINGS
##############################################

inputFile="grainReadFile-2_Molaro_tight.dat"
temp=-80.0
humidity=0.90
dim=2
grad_temp0X=0.0
grad_temp0Y=3.0e-5
grad_temp0Z=0.0
delt_t=1.0e-4
t_final=$(echo "28*24*60*60" | bc -l)
n_out=56

##############################################
# ENVIRONMENT AND FILE PATH SETUP
##############################################

BASE_DIR="${PETIGA_DIR}/projects/dry_snow_metamorphism"
input_dir="$BASE_DIR/inputs"
output_dir="/resnick/scratch/jbaglino"
exec_file="${BASE_DIR}/NASAv2"
SETTINGS_FILE="$BASE_DIR/configs/${inputFile%.dat}.env"
inputFile="$input_dir/$inputFile"

# Validate input file
if [[ ! -f "$inputFile" ]]; then
  echo "[ERROR] Input file '$inputFile' not found. Exiting."
  exit 1
fi

# Generate .env file if missing
if [[ ! -f "$SETTINGS_FILE" ]]; then
  echo "[INFO] .env file not found. Generating..."
  python3 "$BASE_DIR/scripts/generate_env_from_input.py" "$inputFile" "$SETTINGS_FILE" || exit 1
fi

# Load .env
echo "[INFO] Loading settings from $SETTINGS_FILE"
source "$SETTINGS_FILE"

# Output folder
folder="$output_dir/${SLURM_JOB_NAME}_${SLURM_JOB_ID:0:9}"
mkdir -p "$folder"
echo "[INFO] Output directory created: $folder"

# Export environment vars for simulation
export Lx Ly Lz Nx Ny Nz eps delt_t t_final n_out dim \
       grad_temp0X grad_temp0Y grad_temp0Z humidity temp inputFile folder

##############################################
# COMPILE EXECUTABLE IF NEEDED
##############################################

if [[ ! -x "$exec_file" ]]; then
  echo "[INFO] Executable not found. Attempting to compile..."
  cd "$BASE_DIR/src" || { echo "[ERROR] Failed to enter src directory."; exit 1; }
  make NASAv2
  [[ $? -ne 0 ]] && { echo "[ERROR] Compilation failed. Exiting."; exit 1; }
else
  echo "[INFO] Using existing executable: $exec_file"
fi

##############################################
# WRITE METADATA
##############################################

# Attempt to extract number of grains from input file name
if [[ "$inputFile" =~ grainReadFile-([0-9]+) ]]; then
  num_grains="${BASH_REMATCH[1]}"
else
  num_grains="unknown"
fi

json_file="$folder/metadata.json"
echo "[INFO] Writing metadata to: $json_file"

cat << EOF > "$json_file"
{
  "folder_name": "$(basename "$folder")",
  "run_time": "$(date -u +"%Y-%m-%dT%H:%M:%SZ")",
  "executed_on": "$(hostname)",
  "user": "$(whoami)",
  "project": "dry_snow_metamorphism",
  "sim_dimension": $dim,
  "input_file": "$(basename "$inputFile")",
  "env_file_used": "${SETTINGS_FILE:-"none"}",
  "temperature_C": $temp,
  "humidity": $humidity,
  "num_grains": "$num_grains",
  "grad_temp": {
    "x": $grad_temp0X,
    "y": $grad_temp0Y,
    "z": $grad_temp0Z
  },
  "domain_size_m": {
    "Lx": $Lx,
    "Ly": $Ly,
    "Lz": $Lz
  },
  "mesh_resolution": {
    "Nx": $Nx,
    "Ny": $Ny,
    "Nz": $Nz
  },
  "interface_width_eps": $eps,
  "delt_t": $delt_t,
  "t_final": $t_final,
  "n_out": $n_out
}
EOF

##############################################
# BACKUP INPUTS
##############################################

cp "$inputFile" "$folder/"
cp "$BASE_DIR/src/NASAv2.c" "$folder/"
cp "$0" "$folder/run_script_copy.sh"

##############################################
# RUN SIMULATION
##############################################

echo "[INFO] Launching simulation at $(date)"
mpiexec "$exec_file" -initial_cond -initial_PFgeom \
  -snes_rtol 1e-3 -snes_stol 1e-3 -snes_max_it 6 \
  -ksp_gmres_restart 150 -ksp_max_it 500 -ksp_converged_maxits 1 \
  -ksp_converged_reason -snes_converged_reason -snes_linesearch_monitor \
  -snes_linesearch_type basic | tee "$folder/outp.txt"

echo "[INFO] Simulation completed at $(date)"