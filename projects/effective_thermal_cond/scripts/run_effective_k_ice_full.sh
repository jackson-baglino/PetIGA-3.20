#!/bin/zsh
###############################################################################
# run_effective_k_ice_full.sh
#
# This script sets up the environment, compiles the effective_k_ice_full simulation,
# runs the single-file debugging build, and then moves the generated output files.
#
# Optional Post-Processing:
#   At the end you can add commands to process simulation output
#   (e.g., generate plots, convert data).
#
# Usage:
#   ./run_effective_k_ice_full.sh
###############################################################################

# =============================
# 🔹 Environment Variables
# =============================
# export Nx=$((2**8))
# export Ny=$((2**8))
export Nx=550
export Ny=550
export Nz=1                    # Set to 1 for 2D simulations

export Lx=1.0
export Ly=1.0
# export Lx=0.5e-3
# export Ly=0.5e-3
export Lz=2.02e-4              # Only used in 3D mode

# =============================
# 🔹 Boundary Conditions
# =============================
export FLUX_BOTTOM=$((1.0))
export TEMP_TOP=$((273.15-30))

# =============================
# 🔹 Interface Width Calculation
# =============================
export eps=$(awk "BEGIN {print (($Lx/$Nx < $Ly/$Ny) ? $Lx/$Nx : $Ly/$Ny)}")
export dim=2                  # Set 2 for 2D, 3 for 3D

# =============================
# 🔹 Initial Conditions
# =============================
# INIT_MODE can be a filename or one of "circle", "layered", etc.
# INIT_MODE="layered"
# e.g.
# INIT_MODE="circle"
INIT_MODE="/Users/jacksonbaglino/PetIGA-3.20/projects/effective_thermal_cond/inputs/sol_00000.dat"

# =============================
# 🔹 Output Settings
# =============================
export OUTPUT_VTK=1
export OUTPUT_BINARY=1

# timestamped output directory
timestamp=$(date +%Y-%m-%d__%H.%M.%S)
OUT_FOLDER="ThermalSim_full_$timestamp"
export OUTPUT_DIR="/Users/jacksonbaglino/PetIGA-3.20/projects/effective_thermal_cond/outputs/$OUT_FOLDER"
mkdir -p "$OUTPUT_DIR"

# =============================
# 🔹 Simulation Settings
# =============================
NUM_PROCS=1   # Adjust MPI process count as needed

# =============================
# 🔹 Functions
# =============================

compile_code() {
    echo "🔨 Compiling effective_k_ice_full..."
    # only build the full single-file executable
    if make effective_k_ice_full; then
        echo "✅ Compilation successful."
    else
        echo "❌ Compilation failed. Exiting..."
        exit 1
    fi
}

run_simulation() {
    echo "🏃 Running effective_k_ice_full with $NUM_PROCS process(es)..."
    echo
    mpiexec -np $NUM_PROCS ./effective_k_ice_full \
        -init_mode "$INIT_MODE" \
        # add any other PETSc flags here
}

move_output_files() {
    echo "📂 Moving output files to $OUTPUT_DIR..."
    echo
    if [ -d "$OUTPUT_DIR" ]; then
        mv *.dat  "$OUTPUT_DIR" 2>/dev/null || echo "  ⚠️ No .dat files."
        mv *.bin  "$OUTPUT_DIR" 2>/dev/null || echo "  ⚠️ No .bin files."
        mv *.info "$OUTPUT_DIR" 2>/dev/null || echo "  ⚠️ No .info files."
        echo "  ✅ Done."
        echo
    else
        echo "⚠️ Output directory missing!"
    fi
}

###############################################################################
# 🔹 Execution Workflow
###############################################################################
echo -e "\n🌡️ Starting effective_k_ice_full Workflow"
echo "Calculated interface width (eps): $eps"
echo -e "\n-----------------------------------------\n"

compile_code
run_simulation
move_output_files

###############################################################################
# 🔹 Optional Post-Processing Section
###############################################################################
echo "🔍 Running post-processing..."
# e.g.: python postprocess/plot_layered_sys.py "$OUT_FOLDER"
echo "✅ Post-processing complete."

echo -e "\n🎉 Simulation complete: results in $OUTPUT_DIR\n"