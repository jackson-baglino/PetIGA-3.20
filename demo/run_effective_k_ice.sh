#!/bin/bash

# ========== SET ENVIRONMENT VARIABLES ==========
export Nx=128        # Number of elements in x-direction
export Ny=128        # Number of elements in y-direction
export Nz=1          # Number of elements in z-direction (only used if dim=3)

export Lx=1.0e-3        # Length in x-direction (meters)
export Ly=1.0e-3        # Length in y-direction (meters)
export Lz=1.0e-3        # Length in z-direction (meters) (only used if dim=3)

export temp=268.15   # Initial temperature (Kelvin)

# Temperature gradients (currently unused in the code)
export grad_temp0X=0.01  # Temperature gradient in x-direction (K/m)
export grad_temp0Y=0     # Temperature gradient in y-direction (K/m)
export grad_temp0Z=0     # Temperature gradient in z-direction (K/m)

# ========== BOUNDARY CONDITIONS ==========
# Define whether to use fixed temperature or prescribed flux at boundaries

# Bottom boundary
export TEMP_BOTTOM=265.15  # Fixed temperature (Kelvin) at bottom boundary
export FLUX_BOTTOM=1.0    # Set this if prescribing heat flux instead of temperature

# Top boundary
export TEMP_TOP=240.15     # Fixed temperature (Kelvin) at top boundary
export FLUX_TOP=10.0       # Set this if prescribing heat flux instead of temperature

# Interface width (meters)
export eps=$(awk "BEGIN {print ($Lx/$Nx < $Ly/$Ny) ? $Lx/$Nx : $Ly/$Ny}")

export dim=2   # Set dimension (2 for 2D, 3 for 3D)

# ========== INITIAL CONDTIONS ==========
INIT_MODE="circle"  # Initial condition mode (e.g., "circle", "square", etc.)

# ========== OPTIONAL SETTINGS ==========
export OUTPUT_VTK=1    # Set to 1 to enable VTK output
export OUTPUT_BINARY=1 # Set to 1 to enable binary output

# Set timestamp for output directory
timestamp=$(date +%Y-%m-%d__%H.%M.%S)
export OUTPUT_DIR="/Users/jacksonbaglino/SimulationResults/ThermalConductivity/ThermalSim_$timestamp"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Number of processors for MPI execution
NUM_PROCS=1  # Adjust as needed

# ========== FUNCTION: COMPILE THE CODE ==========
compile_code() {
    echo "Compiling the code..."
    make clean
    if make effective_k_ice; then
        echo "Compilation successful."
    else
        echo "Compilation failed. Exiting..."
        exit 1
    fi
}

# ========== FUNCTION: RUN THE SIMULATION ==========
run_simulation() {
    echo "Running effective_k_ice simulation with $NUM_PROCS processes..."
    mpiexec -np $NUM_PROCS ./effective_k_ice -init_mode $INIT_MODE #-ksp_view -ksp_monitor -log_view
}

# ========== FUNCTION: MOVE OUTPUT FILES ==========
move_output_files() {
    if [ -d "$OUTPUT_DIR" ]; then
        # Move binary files only if they exist
        mv *.bin "$OUTPUT_DIR" 2>/dev/null || echo "No binary files to move."
        echo "Output files moved to $OUTPUT_DIR"
    fi
}

# ========== EXECUTE WORKFLOW ==========
echo " "
echo "Starting Thermal Conductivity Simulation Workflow"
echo " "
echo "Calculated interface width (eps): $eps"
echo " "
echo " "

compile_code
run_simulation
move_output_files

echo "Simulation complete."