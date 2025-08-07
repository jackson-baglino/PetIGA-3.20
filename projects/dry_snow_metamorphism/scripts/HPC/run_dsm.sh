#!/bin/bash
#SBATCH -J DSM-T=-40_hum=0.50
#SBATCH -A rubyfu
#SBATCH -t 5-00:00:00
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=40
#SBATCH --cpus-per-task=1
#SBATCH -o "output_files/%x.o%j"
#SBATCH -e "output_files/%x.e%j"
#SBATCH --partition=expansion
#SBATCH --mem-per-cpu=1G
#SBATCH --mail-user=jbaglino@caltech.edu
#SBATCH --mail-type=END,FAIL,TIME_LIMIT

# module load mpi
# module load mpich

##############################################
# USER-MODIFIABLE SIMULATION SETTINGS
##############################################

# Set input filename (only the filename, not full path)
# inputFile="grainReadFile-2G_Molaro_0p25R1.dat"
inputFile="grainReadFile-35_s1-10.dat"

readFlag=1  # Set to 1 to read grain file, 0 to generate grains


# Define physical & environmental parameters
temp=-40.0
humidity=0.50
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
exec_file="${BASE_DIR}/dry_snow_metamorphism"
SETTINGS_FILE="$BASE_DIR/configs/${inputFile%.dat}.env"
inputFile="$input_dir/$inputFile"

# Validate input file
if [[ ! -f "$inputFile" ]]; then
  echo "[ERROR] Input file '$inputFile' not found. Exiting."
  exit 1
fi

# Generate .env file if missing
if [[ ! -f "$SETTINGS_FILE" ]]; then
  echo "[INFO] .env file not found for $inputFile. Generating..."
  python3 "$BASE_DIR/scripts/generate_env_from_input.py" "$inputFile" "$SETTINGS_FILE" || exit 1
fi

# Load simulation parameters from .env
echo "[INFO] Loading settings from $SETTINGS_FILE"
source "$SETTINGS_FILE"

# Output folder
folder="$output_dir/${SLURM_JOB_NAME}_${SLURM_JOB_ID:0:9}"
mkdir -p "$folder"
echo "[INFO] Output directory created: $folder"

# Export for simulation
export Lx Ly Lz Nx Ny Nz eps delt_t t_final n_out dim \
       grad_temp0X grad_temp0Y grad_temp0Z humidity temp inputFile folder readFlag

##############################################
# COMPILE EXECUTABLE IF NEEDED
##############################################

if [[ ! -x "$exec_file" ]]; then
  echo "[INFO] Executable not found. Attempting to compile..."
  make dry_snow_metamorphism
  echo "[INFO] Compilation complete."
  [[ $? -ne 0 ]] && { echo "[ERROR] Compilation failed. Exiting."; exit 1; }
  echo "[INFO] Compilation successful."
else
  echo "[INFO] Using existing executable: $exec_file"
fi

##############################################
# RUN SIMULATION
##############################################

echo "[INFO] Starting simulation on $(date)"
echo "[INFO] Input file: $(basename "$inputFile")"
echo "[INFO] Temperature = $temp°C, Humidity = $humidity"
echo "[INFO] Domain: ($Lx x $Ly x $Lz), Grid: ($Nx x $Ny x $Nz)"

# Backup inputs to output folder
cp "$inputFile" "$folder/"
cp "$BASE_DIR/src/dry_snow_metamorphism.c" "$folder/"
cp "$0" "$folder/run_script_copy.sh"

# Run simulation
mpiexec "$exec_file" -initial_cond -initial_PFgeom \
  -snes_rtol 1e-3 -snes_stol 1e-3 -snes_max_it 6 \
  -ksp_gmres_restart 150 -ksp_max_it 500 -ksp_converged_maxits 1 \
  -ksp_converged_reason -snes_converged_reason -snes_linesearch_monitor \
  -snes_linesearch_type basic | tee "$folder/outp.txt"

echo "[INFO] Simulation completed on $(date)"