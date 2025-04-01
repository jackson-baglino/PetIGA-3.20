#!/bin/bash

# =============================
# 🔹 Environment Variables
# =============================

export Nx=128
export Ny=128
export Nz=1  # Set to 1 for 2D simulations

export Lx=0.5e-3
export Ly=0.5e-3
export Lz=2.02e-4  # Only used in 3D mode

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
export dim=2  # Set 2 for 2D, 3 for 3D

# =============================
# 🔹 Initial Conditions
# =============================

INIT_MODE="/Users/jacksonbaglino/PetIGA-3.20/demo/input/Thermal_IO/circle_phase_field.dat"
# INIT_MODE="circle"

# =============================
# 🔹 Output Settings
# =============================

export OUTPUT_VTK=1
export OUTPUT_BINARY=1

timestamp=$(date +%Y-%m-%d__%H.%M.%S)
export OUTPUT_DIR="/Users/jacksonbaglino/SimulationResults/ThermalConductivity/ThermalSim_$timestamp"

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# =============================
# 🔹 Simulation Settings
# =============================

NUM_PROCS=1  # Adjust as needed

# =============================
# 🔹 Functions
# =============================

compile_code() {
    echo "🔨 Compiling the code..."
    make clean
    if make effective_k_ice; then
        echo "✅ Compilation successful."
    else
        echo "❌ Compilation failed. Exiting..."
        exit 1
    fi
}

run_simulation() {
    echo "🚀 Running effective_k_ice simulation with $NUM_PROCS processes..."
    mpiexec -np $NUM_PROCS ./effective_k_ice -init_mode "$INIT_MODE" # -ksp_view -ksp_monitor -log_view
}

move_output_files() {
    echo "📂 Moving output files..."
    if [ -d "$OUTPUT_DIR" ]; then
        mv *.bin "$OUTPUT_DIR" 2>/dev/null || echo "⚠️ No binary files to move."
        mv *.dat "$OUTPUT_DIR" 2>/dev/null || echo "⚠️ No binary files to move."
        echo "✅ Output files moved to $OUTPUT_DIR"
    fi
}

# =============================
# 🔹 Execution Workflow
# =============================

echo -e "\n🌡️ Starting Thermal Conductivity Simulation Workflow\n"
echo "Calculated interface width (eps): $eps"
echo -e "\n-----------------------------------------\n"

compile_code
run_simulation
move_output_files

echo -e "\n✅ Simulation complete."
