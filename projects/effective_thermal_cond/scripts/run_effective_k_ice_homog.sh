#!/bin/zsh
###############################################################################
# run_effective_k_ice_full.sh
#
# • No verbose tracing, no pipefail, no trap on every error
# • Compiles with *release* flags (BUILD=release)
# • Runs without PETSc debug monitors
# • Still honours the same environment variables / paths
###############################################################################

# ----------------------------
# 🔹  Simulation parameters
# ----------------------------
export Nx=275
export Ny=275
export Nz=1          # 1 for 2-D

export Lx=0.5e-3
export Ly=0.5e-3
export Lz=2.02e-4    # ignored when dim=2

export FLUX_BOTTOM=-0.1
export TEMP_TOP=$((273.15-30))

export eps=$((9.09629658751972e-07))
export dim=2         # 2 = 2-D, 3 = 3-D

# Initial ice-field mode: "circle" | "layered" | /path/to/file.dat
# INIT_MODE="circle"
# INIT_MODE="layered"
INIT_MODE="/Users/jacksonbaglino/PetIGA-3.20/projects/effective_thermal_cond/inputs/"\
"NASAv2_10G_2D_T-20.0_hum0.70_2025-03-13__14.20.59/sol_00444.dat"
INIT_DIR="/Users/jacksonbaglino/PetIGA-3.20/projects/effective_thermal_cond/inputs/"\
"NASAv2_10G_2D_T-20.0_hum0.70_2025-03-13__14.20.59/"

# Output flags
export OUTPUT_VTK=1
export OUTPUT_BINARY=1

# Time-stamped output directory
timestamp=$(date +%Y-%m-%d__%H.%M.%S)
OUT_FOLDER="ThermalSim_homog_$timestamp"
export OUTPUT_DIR="/Users/jacksonbaglino/PetIGA-3.20/projects/effective_thermal_cond/outputs/homog/$OUT_FOLDER"
mkdir -p "$OUTPUT_DIR"

# MPI ranks (override by exporting NUM_PROCS beforehand)
NUM_PROCS=${NUM_PROCS:-1}

# ----------------------------
# 🔹  Helpers
# ----------------------------
compile_code() {
    echo "🔨 Compiling (release)…"
    make BUILD=release effective_k_ice_homog
}

run_simulation() {
    echo "🏃 Running with $NUM_PROCS MPI proc(s)…"
    mpiexec -np $NUM_PROCS ./effective_k_ice_homog \
        -init_mode "$INIT_MODE" -init_dir "$INIT_DIR"
}

collect_outputs() {
    echo "📂 Moving output files to $OUTPUT_DIR"
    mv *.dat  "$OUTPUT_DIR" 2>/dev/null
    mv *.bin  "$OUTPUT_DIR" 2>/dev/null
    mv *.info "$OUTPUT_DIR" 2>/dev/null
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