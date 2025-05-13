#!/bin/zsh
###############################################################################
# run_effective_k_ice.sh
#
# This script sets up the environment, compiles the effective_k_ice simulation,
# runs the simulation, and then moves the generated output files.
#
# Optional Post-Processing:
#   A section has been outlined at the end of this script where you can add
#   commands to process simulation output (e.g., generate plots, convert data).
#
# Usage:
#   ./run_effective_k_ice.sh
###############################################################################

# =============================
# 🔹 Environment Variables
# =============================
export Nx=$((2**8))
export Ny=$((2**8))
export Nz=1                    # Set to 1 for 2D simulations

export Lx=1.0
export Ly=1.0
# export Lx=0.5e-3
# export Ly=0.5e-3
export Lz=2.02e-4              # Only used in 3D mode

# =============================
# 🔹 Boundary Conditions
# =============================
export FLUX_BOTTOM=$((-1.0))
export TEMP_TOP=$((273.15-30))

# =============================
# 🔹 Interface Width Calculation
# =============================
export eps=$(awk "BEGIN {print (($Lx/$Nx < $Ly/$Ny) ? $Lx/$Nx : $Ly/$Ny)}")
export dim=2                  # Set 2 for 2D, 3 for 3D

# =============================
# 🔹 Initial Conditions
# =============================
# Set the initial mode. Uncomment and modify the path if using a file.
# INIT_MODE="/Users/jacksonbaglino/PetIGA-3.20/projects/effective_thermal_cond/inputs/sol_00444.dat"
# INIT_MODE="/Users/jacksonbaglino/PetIGA-3.20/projects/effective_thermal_cond/inputs/LayeredSystem/NoInclusions1/sol_00778.dat"
# INIT_MODE="circle"
INIT_MODE="layered"

# =============================
# 🔹 Output Settings
# =============================
export OUTPUT_VTK=1
export OUTPUT_BINARY=1

# Create a timestamp and define the output directory (change the path as needed)
timestamp=$(date +%Y-%m-%d__%H.%M.%S)
# export OUTPUT_DIR="/Users/jacksonbaglino/SimulationResults/ThermalConductivity/ThermalSim_$timestamp"
OUT_FOLDER="ThermalSim_$timestamp"
export OUTPUT_DIR="/Users/jacksonbaglino/PetIGA-3.20/projects/effective_thermal_cond/outputs/$OUT_FOLDER"

# Ensure the output directory exists
mkdir -p "$OUTPUT_DIR"

# =============================
# 🔹 Simulation Settings
# =============================
NUM_PROCS=1   # Adjust the number of MPI processes as needed

# =============================
# 🔹 Functions
# =============================

compile_code() {
    echo "🔨 Compiling the code..."
    if make all; then
        echo "✅ Compilation successful."
    else
        echo "❌ Compilation failed. Exiting..."
        exit 1
    fi
}

run_simulation() {
    echo "Running effective_k_ice simulation with $NUM_PROCS process(es)..."
    echo " "
    mpiexec -np $NUM_PROCS ./effective_k_ice -init_mode "$INIT_MODE" # -ksp_monitor # Additional flags can be added here
}

move_output_files() {
    echo "📂 Moving output files..."
    echo " "
    if [ -d "$OUTPUT_DIR" ]; then
        mv *.dat "$OUTPUT_DIR" 2>/dev/null || echo "  ⚠️ No .dat files to move."
        mv *.bin "$OUTPUT_DIR" 2>/dev/null || echo "  ⚠️ No .bin files to move."
        mv *.info "$OUTPUT_DIR" 2>/dev/null || echo " ⚠️ No .info files to move."
        echo "  ✅ Output files moved to $OUTPUT_DIR"
        echo " "
    else
        echo "⚠️ Output directory not found."
    fi
}

###############################################################################
# 🔹 Execution Workflow
###############################################################################
echo -e "\n🌡️ Starting Thermal Conductivity Simulation Workflow\n"
echo "Calculated interface width (eps): $eps"
echo -e "\n-----------------------------------------\n"

compile_code
run_simulation
move_output_files

###############################################################################
# 🔹 Optional Post-Processing Section
#
# Uncomment and modify the commands below to perform any post-processing on
# your simulation output (e.g., plotting, data conversion, analysis).
###############################################################################
echo "🔍 Running post-processing..."
python postprocess/plot_layered_sys.py "$OUT_FOLDER"
echo "✅ Post-processing complete."

echo -e "Simulation complete."