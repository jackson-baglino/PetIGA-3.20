#!/bin/zsh
###############################################################################
# Script: run_dsm.sh
# Model: Dry Snow Metamorphism (DSM)
# Purpose:
#   Configure and launch a single DSM simulation, manage run metadata, and
#   archive inputs/plots to a timestamped results folder.
#
# Typical Usage:
#   1) Set `filename` to the desired grain input (in $BASE_DIR/inputs), or set
#      `readFlag=0` to generate grains instead of reading a file.
#   2) Adjust physical & numerical params below (t_final, n_out, gradients, etc.).
#   3) Run: `./scripts/studio/run_dsm.sh`
#
# Key Variables:
#   BASE_DIR      Project root for DSM model
#   input_dir     Directory with input grain files (*.dat)
#   output_dir    Parent directory to store run outputs (timestamped subfolder)
#   exec_file     DSM executable to run
#   filename      (When readFlag=1) The grain input file name in input_dir
#   readFlag      1=read grain file; 0=procedural generation
#   t_final       Total simulated time [s]
#   n_out         Number of output frames (approx.)
#   grad_temp0*   Temperature gradient components [K/m]
#   dim           2D or 3D
#   NUM_PROCS     MPI ranks (default set below)
#
# Outputs:
#   - Creates $output_dir/<title>_<timestamp>/
#   - Copies input & env files used, code snapshot, and plotting scripts
#   - Saves outp.txt, simulation_parameters.csv, sim_params.dat, metadata.json
#
# Notes:
#   - This script avoids deleting any user comments.
#   - If `readFlag=1`, ensure `filename` is set to a valid file in $input_dir.
#   - Settings file is loaded from $BASE_DIR/configs/${filename%.dat}.env.
###############################################################################

# =======================================
# Base directory and file setup
# =======================================
BASE_DIR="${PETIGA_DIR}/projects/dry_snow_metamorphism"
input_dir="$BASE_DIR/inputs"
output_dir="/Users/jacksonbaglino/SimulationResults/dry_snow_metamorphism/scratch"
exec_file="${BASE_DIR}/dry_snow_metamorphism"

# =======================================
# Define simulation parameters
# =======================================
# filename="grainReadFile-2G_Molaro_0p25R1_HIGHRES.dat"
# Batch override precedence: if $filename is exported, use it; else use this default
: ${filename:="grainReadFile-2G_Molaro_0p25R1.dat"}
inputFile="$input_dir/$filename"

# --- Basic validation (only when reading a grain file) ---
if [[ "${readFlag:-1}" -eq 1 ]]; then
  if [[ -z "$filename" ]]; then
    echo "[ERROR] readFlag=1 but 'filename' is not set. Set it above (e.g., filename=\"grainReadFile-...dat\")."
    exit 1
  fi
  if [[ ! -f "$inputFile" ]]; then
    echo "[ERROR] Input file not found: $inputFile"
    exit 1
  fi
fi

readFlag=1  # Set to 1 to read grain file, 0 to generate grains

delt_t=1.0e-4
# t_final=$((28 * 24 * 60 * 60))  # 14 days in seconds
t_final=$((12 * 60 * 60))  # 2 hours in seconds
# t_final=10
n_out=100
# Batch override precedence: use exported values if provided; otherwise defaults
: ${humidity:=0.95}
: ${temp:=-2.0}
grad_temp0X=0.0
grad_temp0Y=3.0e-6
grad_temp0Z=0.0
dim=2

if [[ readFlag -eq 1 ]]; then
    echo "[INFO] Reading grain file: $inputFile"
else
    echo "[INFO] Generating grains instead of reading from file."
    Lx=0.5e-3
    Ly=0.5e-3
    Lz=0.5e-3

    Nx=275
    Ny=275
    Nz=275

    eps=9.00e-07
fi

# Build a descriptive run title encoding key parameters for easier indexing
clean_name="${filename#grainReadFile-}"
clean_name="${clean_name%.dat}"

# Create compact two-digit tags for temperature (°C) and humidity (%)
# temp_tag: round to nearest integer, then take sign + first two digits
temp_int=$(printf "%.0f" "$temp")
if [[ "$temp_int" == -* ]]; then
  temp_tag=${temp_int:0:3}   # e.g., -12, -9
else
  temp_tag=${temp_int:0:2}   # e.g., 15 -> 15, 7 -> 7
fi
# hum_tag: integer percent, then first two digits (e.g., 98, 95, 100 -> 10)
hum_int=$(awk "BEGIN{printf \"%d\", $humidity*100}")
hum_tag=$(printf "%02d" "$hum_int")
hum_tag=${hum_tag:0:2}

# Build title using compact tags
# Note: keep existing day count in suffix
ndays=$(awk "BEGIN{printf \"%d\", $t_final/86400}")

title="SCRATCH_DSM${clean_name}_${dim}D_Tm${temp_tag}_hum${hum_tag}_tf${ndays}d_"
SETTINGS_FILE="$BASE_DIR/configs/${filename%.dat}.env"

# MPI ranks used for this run (override here or export before calling)
NUM_PROCS=10  # Number of MPI processes

# =======================================
# Timestamped result folder
# =======================================
timestamp=$(date +%Y-%m-%d__%H.%M.%S)
folder="$output_dir/${title}_$timestamp"
mkdir -p "$folder"

# Copy input file to results folder
cp "$inputFile" "$folder"
cp "$SETTINGS_FILE" "$folder"

# =======================================
# Build and run setup
# =======================================
cd "$BASE_DIR" || exit 1

# Compile the DSM executable (requires PETSc/PetIGA set up)
make dry_snow_metamorphism || {
    echo "[ERROR] Build failed. Please check the Makefile and dependencies."
    exit 1
}

# Generate env file if missing
if [ ! -f "$SETTINGS_FILE" ]; then
    echo "[INFO] .env file not found. Generating from input..."
    python3 scripts/generate_env_from_input.py "$inputFile" "$SETTINGS_FILE"
fi

# Load run-time physical/mesh parameters from the settings .env
set -a
source "$SETTINGS_FILE"
set +a

echo "Settings file loaded: $SETTINGS_FILE"
echo "Lx = $Lx  Ly = $Ly Lz = $Lz"

# Export simulation parameters
export folder input_dir inputFile filename title
export delt_t t_final n_out humidity temp grad_temp0X grad_temp0Y grad_temp0Z dim
export readFlag Lx Ly Lz Nx Ny Nz eps

# =======================================
# Save simulation metadata
# =======================================

# Write simulation parameters to CSV file for easy spreadsheet viewing
write_parameters_to_csv() {
    csv="$folder/simulation_parameters.csv"
    echo "Variable,Value" > "$csv"
    for var in folder inputFile title Lx Ly Lz Nx Ny Nz delt_t t_final n_out humidity temp grad_temp0X grad_temp0Y grad_temp0Z dim eps
    do
        echo "$var,${(P)var}" >> "$csv"
    done
}

# Write simulation parameters to a human-readable .dat file
write_parameters_to_dat() {
    cat << EOF > "$folder/sim_params.dat"
----- SIMULATION PARAMETERS -----
Input file: $inputFile
Dimensions: dim = $dim
eps = $eps
Lx = $Lx, Ly = $Ly, Lz = $Lz
Nx = $Nx, Ny = $Ny, Nz = $Nz
delt_t = $delt_t, t_final = $t_final, n_out = $n_out
humidity = $humidity, temp = $temp
grad_temp0X = $grad_temp0X, grad_temp0Y = $grad_temp0Y, grad_temp0Z = $grad_temp0Z
EOF
}

# Write a fully-resolved .env-style snapshot for reproducibility
write_env_snapshot() {
    snapshot="$folder/resolved_params.env"
    echo "# Auto-generated resolved parameters for this run" > "$snapshot"
    echo "# Generated: $(date -u +"%Y-%m-%dT%H:%M:%SZ")" >> "$snapshot"
    for var in folder inputFile title Lx Ly Lz Nx Ny Nz delt_t t_final n_out humidity temp grad_temp0X grad_temp0Y grad_temp0Z dim eps readFlag; do
        echo "$var=${(P)var}" >> "$snapshot"
    done
}

# Function to write metadata.json
write_metadata_json() {
    json_file="$folder/metadata.json"

    cat << EOF > "$json_file"
{
  "folder_name": "$(basename "$folder")",
  "run_time": "$(date -u +"%Y-%m-%dT%H:%M:%SZ")",
  "executed_on": "$(hostname)",
  "user": "$(whoami)",
  "project": "dry_snow_metamorphism",
  "sim_dimension": $dim,
  "input_file": "$inputFile",
  "env_file_used": "${SETTINGS_FILE:-"none"}",
  "temperature_C": $temp,
  "humidity": $humidity,
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
}

write_parameters_to_csv
write_parameters_to_dat

# Also write machine-readable metadata and a resolved .env snapshot
write_metadata_json
write_env_snapshot

# Append resolved parameters to the copied settings file for a single self-contained record
COPIED_ENV_FILE="$folder/$(basename "$SETTINGS_FILE")"
{
  echo ""
  echo "# ---- Resolved run-time parameters (auto-generated) ----"
  echo "# Note: Values below reflect the actual run configuration (after overrides)."
  for var in folder inputFile title Lx Ly Lz Nx Ny Nz delt_t t_final n_out humidity temp grad_temp0X grad_temp0Y grad_temp0Z dim eps readFlag; do
      echo "$var=${(P)var}"
  done
} >> "$COPIED_ENV_FILE"

# =======================================
# Run the simulation
# =======================================

# Solver tolerances and options are set here
echo "[INFO] Launching DRY SNOW METAMORPHISM simulation..."
mpiexec -np "$NUM_PROCS" "$exec_file" -initial_PFgeom -temp_initial \
  -snes_rtol 1e-3 -snes_stol 1e-6 -snes_max_it 7 \
  -ksp_gmres_restart 150 -ksp_max_it 1000 \
  -ksp_converged_reason -snes_converged_reason -snes_linesearch_monitor \
  -snes_linesearch_type basic | tee "$folder/outp.txt"

# =======================================
# Finalize
# =======================================
echo " "
echo "[INFO] Simulation completed."
cp -r src scripts/studio/run_dsm.sh postprocess/plotDSM.py postprocess/plotSSA.py postprocess/plotPorosity.py "$folder"
echo "Copied files to $folder"

echo " "
echo "Running plotting scripts..."
./scripts/run_plotDSM.sh

echo "✅ Simulation complete. Results stored in:"
echo "$folder"