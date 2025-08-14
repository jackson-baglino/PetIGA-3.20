#!/bin/zsh
###############################################################################
# run_effective_k_ice_full.sh
#
# • No verbose tracing, no pipefail, no trap on every error
# • Compiles with *release* flags (BUILD=release)
# • Runs without PETSc debug monitors
# • Still honours the same environment variables / paths
###############################################################################

start_time=$(date +%s)  # seconds since epoch

# ----------------------------
# 🔹  Simulation parameters
# ----------------------------
# export Nx=1100
# export Ny=1100
# export Nz=1          # 1 for 2-D

# export Nx=134
# export Ny=214

# export Lx=2.03e-3
# export Ly=2.03e-3
# export Lz=2.02e-4    # ignored when dim=2

FLUX_BOTTOM=-0.1
TEMP_TOP=$((273.15-30))

# export eps=$((9.09629658751972e-07))
# export dim=2         # 2 = 2-D, 3 = 3-D

# Initial ice-field mode: "circle" | "layered" | /path/to/file.dat
# INIT_MODE="circle"
# INIT_MODE="layered"
INIT_MODE="FILE"
# INIT_DIR="/Users/jacksonbaglino/PetIGA-3.20/projects/effective_thermal_cond/inputs/"\
# "NASAv2_96G-2D_T-20.0_hum0.70_2025-05-31__18.55.56"
INIT_DIR="/Users/jacksonbaglino/SimulationResults/dry_snow_metamorphism/scratch/DSM2G_Molaro_0p25R1_2D_Tm12_hum5_tf0d__2025-08-07__18.59.47"

# Locate and source .env file (portable across bash/zsh)
if [[ -z "$ENV_FILE" ]]; then
    echo "🔎 Searching for .env in: $INIT_DIR"
    # Use find so we don't rely on shell globs/nullglob behavior
    env_files=( $(find "$INIT_DIR" -maxdepth 1 -type f -name '*.env' 2>/dev/null) )
    env_count=${#env_files[@]}
    if [[ -z "$env_count" || "$env_count" -eq 0 ]]; then
        echo "❌ No .env file found in $INIT_DIR"
        exit 1
    elif [[ "$env_count" -gt 1 ]]; then
        echo "❌ Multiple .env files found in $INIT_DIR:";
        for f in "${env_files[@]}"; do echo "  $f"; done
        echo "   Set ENV_FILE=/path/to.env to choose one."
        exit 1
    else
        ENV_FILE="${env_files[1]}"
    fi
fi

if [[ ! -f "$ENV_FILE" ]]; then
    echo "❌ .env file not found: $ENV_FILE"
    exit 1
fi

# Diagnostic: show which env file is used
echo "✅ Using env file: $ENV_FILE"

# Source the .env file to get variables
source "$ENV_FILE"

# Set dim if not already set in script
if [[ -z "$dim" ]]; then
    if [[ -z "$dim" ]]; then
        dim=2
    fi
fi

# Export variables for simulation
export Nx
export Ny
export Nz
export Lx
export Ly
export Lz
export eps
export TEMP_TOP
export FLUX_BOTTOM
export dim

# Extract base folder name
base_folder=$(basename "$INIT_DIR")

grains=2
dims=2
temp=-20.0
hum=1.00

# Timestamped output folder name
timestamp=$(date +%Y-%m-%d__%H.%M.%S)
OUT_FOLDER="ThermalSim_homog_${grains}G_${dims}D_T${temp}_hum${hum}_$timestamp"

echo "✅ OUT_FOLDER: $OUT_FOLDER"

echo "Loaded parameters from .env:"
echo "  dim=$dim, Nx=$Nx, Ny=$Ny, Nz=$Nz, Lx=$Lx, Ly=$Ly, Lz=$Lz, eps=$eps, TEMP_TOP=$TEMP_TOP"

# Output flags
export OUTPUT_VTK=1
export OUTPUT_BINARY=1
export SOL_INDEX=-1

export OUTPUT_DIR="/Users/jacksonbaglino/PetIGA-3.20/projects/effective_thermal_cond/outputs/homog/$OUT_FOLDER"
mkdir -p "$OUTPUT_DIR"

# MPI ranks (override by exporting NUM_PROCS beforehand)
NUM_PROCS=${NUM_PROCS:-1}
# NUM_PROCS=4

# ----------------------------
# 🔹  Helpers
# ----------------------------
compile_code() {
    echo "Compiling (release)…"
    make BUILD=release effective_k_ice_homog
}

run_simulation() {
    echo " "
    echo "Running with $NUM_PROCS MPI proc(s)…"
    mpiexec -np $NUM_PROCS ./effective_k_ice_homog \
        -init_mode "$INIT_MODE" \
        -init_dir "$INIT_DIR" \
}

collect_outputs() {
    echo "📂 Moving output files to $OUTPUT_DIR"
    mv *.dat  "$OUTPUT_DIR" 2>/dev/null
    mv *.bin  "$OUTPUT_DIR" 2>/dev/null
    mv *.info "$OUTPUT_DIR" 2>/dev/null
    mv *.csv  "$OUTPUT_DIR" 2>/dev/null
}

# ----------------------------
# 🔹  Main
# ----------------------------
compile_code   || { echo "❌ compile failed"; exit 1; }
run_simulation || { echo "❌ simulation failed"; exit 1; }
collect_outputs

echo "📈 Post-processing…"
python3 postprocess/plot_vector_field.py "$OUT_FOLDER" "$Nx" "$Ny" "$Lx" "$Ly"

echo "✅ Finished. Outputs in: $OUTPUT_DIR"

end_time=$(date +%s)

elapsed=$(( end_time - start_time ))

echo "Simulation completed in $elapsed seconds."