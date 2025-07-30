#!/bin/zsh

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
<<<<<<< HEAD
filename="grainReadFile-2G_Molaro.dat"
# filename="grainReadFile-18FCC.dat"
# filename="grainReadFile-30G_s1-10.dat"
# filename="grainReadFile-2G_Molaro_tight.dat"
=======
# filename="grainReadFile-2G_Molaro_tight.dat"
filename="grainReadFile-18FCC.dat"
>>>>>>> parent of b85b434 (Refactor simulation scripts and update metadata handling for effective thermal conductivity calculations)
inputFile="$input_dir/$filename"
readFlag=1  # Set to 1 to read grain file, 0 to generate grains

delt_t=1.0e-4
<<<<<<< HEAD
# Define simulation end time in weeks, days, hours, and seconds
weeks=0
days=0
hours=0
seconds=100
t_final=$((weeks * 7 * 24 * 60 * 60 + days * 24 * 60 * 60 + hours * 60 * 60 + seconds))
n_out=10
humidity=0.5  # Relative humidity in fraction (0.0 to 1.0)
temp=-10.0

=======
t_final=1e-4
t_final=$((28 * 24 * 60 * 60))  # 14 days in seconds
n_out=40
humidity=1.00
temp=-20.0
>>>>>>> parent of b85b434 (Refactor simulation scripts and update metadata handling for effective thermal conductivity calculations)
grad_temp0X=0.0
grad_temp0Y=3.0e0
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

clean_name="${filename#grainReadFile-}"
clean_name="${clean_name%.dat}"

title="drysnow_${clean_name}_${dim}D_Tm${temp/-}_hum$(printf "%.0f" "$(echo "$humidity * 100" | bc -l)")_tf$(echo "$t_final / 86400" | bc)d_"
SETTINGS_FILE="$BASE_DIR/configs/${filename%.dat}.env"

NUM_PROCS=12  # Number of MPI processes

# =======================================
# Timestamped result folder
# =======================================
timestamp=$(date +%Y-%m-%d__%H.%M.%S)
folder="$output_dir/${title}_$timestamp"
mkdir -p "$folder"

# Copy input file to results folder
cp "$inputFile" "$folder"

# =======================================
# Build and run setup
# =======================================
cd "$BASE_DIR" || exit 1
make dry_snow_metamorphism || {
    echo "[ERROR] Build failed. Please check the Makefile and dependencies."
    exit 1
}

# Generate env file if missing
if [ ! -f "$SETTINGS_FILE" ]; then
    echo "[INFO] .env file not found. Generating from input..."
    python3 scripts/generate_env_from_input.py "$inputFile" "$SETTINGS_FILE"
fi

# Source the env file
set -a
source "$SETTINGS_FILE"
set +a

# Export simulation parameters
export folder input_dir inputFile filename title
export delt_t t_final n_out humidity temp grad_temp0X grad_temp0Y grad_temp0Z dim
export readFlag Lx Ly Lz Nx Ny Nz eps

# =======================================
# Save simulation metadata
# =======================================
write_parameters_to_csv() {
    csv="$folder/simulation_parameters.csv"
    echo "Variable,Value" > "$csv"
    for var in folder inputFile title Lx Ly Lz Nx Ny Nz delt_t t_final n_out humidity temp grad_temp0X grad_temp0Y grad_temp0Z dim eps
    do
        echo "$var,${(P)var}" >> "$csv"
    done
}

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

# Function to write metadata.json
write_metadata_json() {
    json_file="$folder/metadata.json"

    cat << EOF > "$json_file"
{
  "folder_name": "$(basename "$folder")",
  "run_time": "$(date --iso-8601=seconds)",
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

# =======================================
# Run the simulation
# =======================================
echo "[INFO] Launching DRY SNOW METAMORPHISM simulation..."
mpiexec -np "$NUM_PROCS" "$exec_file" \
  -snes_rtol 1e-3 -snes_stol 1e-6 -snes_max_it 7 \
  -ksp_gmres_restart 150 -ksp_max_it 1000 \
  -ksp_converged_reason -snes_converged_reason -snes_linesearch_monitor \
  -snes_linesearch_type basic | tee "$folder/outp.txt"


# mpiexec -np "$NUM_PROCS" "$exec_file" -initial_PFgeom -temp_initial \
#   -snes_rtol 1e-3 -snes_stol 1e-6 -snes_max_it 7 \
#   -ksp_gmres_restart 150 -ksp_max_it 1000 \
#   -ksp_converged_reason -snes_converged_reason -snes_linesearch_monitor \
#   -snes_linesearch_type basic | tee "$folder/outp.txt"

# For debugging...
# mpiexec -np "$NUM_PROCS" "$exec_file" \
#   -snes_rtol 1e-5 \
#   -snes_stol 1e-10 \
#   -snes_max_it 30 \
#   -snes_monitor_short \
#   -snes_monitor_solution draw \
#   -snes_converged_reason \
#   -ksp_monitor_true_residual \
#   -ksp_monitor_short \
#   -ksp_converged_reason \
#   -ksp_view \
#   -snes_linesearch_monitor \
#   -info \
#   -log_view | tee "$folder/outp.txt"


# =======================================
# Finalize
# =======================================
echo " "
echo "[INFO] Simulation completed."
cp -r src scripts/studio/run_dsm.sh postprocess/plotDSM.py postprocess/plotSSA.py postprocess/plotPorosity.py "$folder"
./scripts/run_plotDSM.sh

echo "✅ Simulation complete. Results stored in:"
echo "$folder"