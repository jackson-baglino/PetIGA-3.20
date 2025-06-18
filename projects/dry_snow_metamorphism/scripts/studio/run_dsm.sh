#!/bin/zsh

# =======================================
# Base directory and file setup
# =======================================
BASE_DIR="${PETIGA_DIR}/projects/dry_snow_metamorphism"
input_dir="$BASE_DIR/inputs"
output_dir="/Users/jacksonbaglino/SimulationResults/DrySed_Metamorphism/NASAv2"
exec_file="${BASE_DIR}/NASAv2"

# =======================================
# Define simulation parameters
# =======================================
filename="grainReadFile-2_Molaro_tight.dat"
inputFile="$input_dir/$filename"

delt_t=1.0e-4
t_final=$((14 * 24 * 60 * 60))  # 14 days in seconds
n_out=40
humidity=1.00
temp=-80.0
grad_temp0X=0.0
grad_temp0Y=3.0e-4
grad_temp0Z=0.0
dim=2

# title="NASAv2_30G-${dim}D_T${temp}_hum${humidity}"
title="drysnow_30g_dim${dim}_Tm${temp/-/m}_hum$(printf "%.0f" "$(echo "$humidity * 100" | bc -l)")_dt${delt_t}_tf$(echo "$t_final / 86400" | bc)d_"
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
make NASAv2

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
export folder input_dir inputFile filename title \
       delt_t t_final n_out humidity temp grad_temp0X grad_temp0Y grad_temp0Z dim \
       Lx Ly Lz Nx Ny Nz eps

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
echo "[INFO] Launching NASAv2 simulation..."
mpiexec -np "$NUM_PROCS" "$exec_file" -initial_PFgeom -temp_initial \
  -snes_rtol 1e-3 -snes_stol 1e-6 -snes_max_it 7 \
  -ksp_gmres_restart 150 -ksp_max_it 1000 \
  -ksp_converged_reason -snes_converged_reason -snes_linesearch_monitor \
  -snes_linesearch_type basic | tee "$folder/outp.txt"

# =======================================
# Finalize
# =======================================
cp "$exec_file.c" run_NASAv2.sh plotNASA.py plotSSA.py plotPorosity.py "$folder"
./run_plotNASAv2.sh "$title"  # Optional if plotting is scripted

echo "✅ Simulation complete. Results stored in:"
echo "$folder"