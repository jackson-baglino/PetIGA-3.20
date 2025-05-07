#!/bin/zsh
###############################################################################
# run_effective_k_ice_full_debug.sh
#
# Debug version of run_effective_k_ice_full.sh:
# - verbose tracing of every command
# - pipefail and immediate exit on error
# - compiles with debug flags (BUILD=debug)
# - runs with extra PETSc debug monitors & log viewers
#
# Usage:
#   ./run_effective_k_ice_full_debug.sh
###############################################################################

set -euxo pipefail
trap 'echo "❌ Error on line $LINENO"; exit 1' ERR

# =============================
# 🔹 Environment Variables
# =============================
# export Nx=$((2**8))
# export Ny=$((2**8))
export Nx=550
export Ny=550
export Nz=1                    # 1 for 2D

export Lx=1.0
export Ly=1.0
# export Lx=0.5e-3
# export Ly=0.5e-3
export Lz=2.02e-4              # only used if dim=3

# =============================
# 🔹 Boundary Conditions
# =============================
export FLUX_BOTTOM=-0.1
export TEMP_TOP=$((273.15-30))

# export FLUX_BOTTOM=1.0
# export TEMP_TOP=$((273.15-4))

# =============================
# 🔹 Interface Width Calculation
# =============================
# export eps=$(awk "BEGIN {print (($Lx/$Nx < $Ly/$Ny) ? $Lx/$Nx : $Ly/$Ny)}")
# export eps=9.09629658751972e-04
export eps=1e-4
export dim=2                   # 2 for 2D, 3 for 3D

# =============================
# 🔹 Initial Conditions
# =============================
# INIT_MODE can be "circle", "layered", or a path to .dat
# INIT_MODE="layered"
# INIT_MODE="circle"
INIT_MODE="/Users/jacksonbaglino/PetIGA-3.20/projects/effective_thermal_cond/inputs/sol_00000.dat"

# =============================
# 🔹 Output Settings
# =============================
export OUTPUT_VTK=1
export OUTPUT_BINARY=1

# timestamped output directory
timestamp=$(date +%Y-%m-%d__%H.%M.%S)
OUT_FOLDER="ThermalSim_full_debug_$timestamp"
export OUTPUT_DIR="/Users/jacksonbaglino/PetIGA-3.20/projects/effective_thermal_cond/outputs/$OUT_FOLDER"
mkdir -p "$OUTPUT_DIR"

# =============================
# 🔹 MPI Settings
# =============================
NUM_PROCS=${NUM_PROCS:-1}      # override by exporting before calling

# =============================
# 🔹 Functions
# =============================
compile_code() {
    echo "🔨 [DEBUG] Compiling effective_k_ice_full with debug flags..."
    make BUILD=debug effective_k_ice_full
    echo "✅ Compilation (debug) successful."
}

run_simulation() {
    echo "🏃 [DEBUG] Running effective_k_ice_full with $NUM_PROCS MPI proc(s)..."
    echo
    mpiexec -np $NUM_PROCS ./effective_k_ice_full \
        -init_mode "$INIT_MODE" \
        -ksp_monitor \
        -ksp_converged_reason \
        -log_summary
        # add any other -log_view or PETSc debug flags here
}

move_output_files() {
    echo "📂 [DEBUG] Moving output files to $OUTPUT_DIR..."
    mv *.dat  "$OUTPUT_DIR" 2>/dev/null || echo "  ⚠️ No .dat files."
    mv *.bin  "$OUTPUT_DIR" 2>/dev/null || echo "  ⚠️ No .bin files."
    mv *.info "$OUTPUT_DIR" 2>/dev/null || echo "  ⚠️ No .info files."
    echo "✅ Files moved."
}

# =============================
# 🔹 Main Workflow
# =============================
echo -e "\n🛠 [DEBUG] Starting debug workflow for effective_k_ice_full"
echo " interface width (eps) = $eps"
echo -e "\n----------------------------------------\n"

compile_code
run_simulation
move_output_files

echo "Plotting output files..."
python3 postprocess/plot_layered_sys.py $OUT_FOLDER

echo -e "\n🎉 [DEBUG] Simulation (debug) complete. Results in:$OUTPUT_DIR\n"