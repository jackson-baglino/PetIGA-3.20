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
export Nz=1          # 1 for 2-D

export Nx=400
export Ny=400

export Lx=0.000600
export Ly=0.000600
export Lz=2.02e-4    # ignored when dim=2

export FLUX_BOTTOM=-0.1
export TEMP_TOP=$((273.15-30))

export eps=$((9.09629658751972e-07))
export dim=2         # 2 = 2-D, 3 = 3-D

# Initial ice-field mode: "circle" | "layered" | /path/to/file.dat
# INIT_MODE="circle"
# INIT_MODE="layered"
INIT_MODE="FILE"
# INIT_DIR="/Users/jacksonbaglino/PetIGA-3.20/projects/effective_thermal_cond/inputs/"\
# "NASAv2_96G-2D_T-20.0_hum0.70_2025-05-31__18.55.56"
INIT_DIR="/Users/jacksonbaglino/SimulationResults/dry_snow_metamorphism/scratch/drysnow_18FCC_2D_Tm20.0_hum100_tf28d__2025-06-20__14.40.41"
base_folder=$(basename "$INIT_DIR")

grains=2
dims=2
temp=-80.0
hum=1.00


echo "grains: $grains"
echo "dims: $dims"
echo "temp: $temp"
echo "humidity: $hum"

timestamp=$(date +%Y-%m-%d__%H.%M.%S)
OUT_FOLDER="ThermalSim_homog_${grains}G_${dims}D_T${temp}_hum${hum}_$timestamp"
echo "OUT_FOLDER: $OUT_FOLDER"

echo "✅ OUT_FOLDER: $OUT_FOLDER"

echo "✅ OUT_FOLDER: $OUT_FOLDER"

# Output flags
export OUTPUT_VTK=1
export OUTPUT_BINARY=1
export SOL_INDEX=-1

export OUTPUT_DIR="/Users/jacksonbaglino/SimulationResults/effective_thermal_cond/scratch/$OUT_FOLDER"
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

# Move SSA file from input to output directory
SSA_FILE="$INIT_DIR/SSA_evo.dat"
cp "$SSA_FILE" "$OUTPUT_DIR/SSA_evo.dat" 2>/dev/null || echo "SSA file not found, skipping copy."

end_time=$(date +%s)

elapsed=$(( end_time - start_time ))

echo "Simulation completed in $elapsed seconds."