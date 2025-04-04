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
export Nx=128
export Ny=128
export Nz=1                    # Set to 1 for 2D simulations

export Lx=0.5e-3
export Ly=0.5e-3
export Lz=2.02e-4              # Only used in 3D mode

export temp=268.15

# Temperature Gradients (currently unused)
export grad_temp0X=0.01
export grad_temp0Y=0
export grad_temp0Z=0

# =============================
# 🔹 Boundary Conditions
# =============================
export TEMP_BOTTOM=265.15
export FLUX_BOTTOM=1.0

export TEMP_TOP=240.15
export FLUX_TOP=10.0

# =============================
# 🔹 Interface Width Calculation
# =============================
export eps=$(awk "BEGIN {print ($Lx/$Nx < $Ly/$Ny) ? $Lx/$Nx : $Ly/$Ny}")
export dim=2                  # Set 2 for 2D, 3 for 3D

# =============================
# 🔹 Initial Conditions
# =============================
# Set the initial mode. Uncomment and modify the path if using a file.
# INIT_MODE="/path/to/your/inputs/circle_phase_field.dat"
INIT_MODE="circle"

# =============================
# 🔹 Output Settings
# =============================
export OUTPUT_VTK=1
export OUTPUT_BINARY=1

# Create a timestamp and define the output directory (change the path as needed)
timestamp=$(date +%Y-%m-%d__%H.%M.%S)
export OUTPUT_DIR="/Users/jacksonbaglino/SimulationResults/ThermalConductivity/ThermalSim_$timestamp"

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
    make clean
    if make all; then
        echo "✅ Compilation successful."
    else
        echo "❌ Compilation failed. Exiting..."
        exit 1
    fi
}

run_simulation() {
    echo "🚀 Running effective_k_ice simulation with $NUM_PROCS process(es)..."
    mpiexec -np $NUM_PROCS ./effective_k_ice -init_mode "$INIT_MODE" # Additional flags can be added here
}

move_output_files() {
    echo "📂 Moving output files..."
    if [ -d "$OUTPUT_DIR" ]; then
        mv *.bin "$OUTPUT_DIR" 2>/dev/null || echo "⚠️ No .bin files to move."
        mv *.dat "$OUTPUT_DIR" 2>/dev/null || echo "⚠️ No .dat files to move."
        mv *.info "$OUTPUT_DIR" 2>/dev/null || echo "⚠️ No .info files to move."
        echo "✅ Output files moved to $OUTPUT_DIR"
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
# echo "🔍 Running post-processing..."
# python postprocess_script.py "$OUTPUT_DIR"
# echo "✅ Post-processing complete."

echo -e "\n✅ Simulation complete."