#!/bin/bash

# ========== SET ENVIRONMENT VARIABLES ==========
export Nx=10        # Number of elements in x-direction
export Ny=10        # Number of elements in y-direction
export Nz=1        # Number of elements in z-direction (only used if dim=3)

export Lx=1.0      # Length in x-direction (meters)
export Ly=1.0      # Length in y-direction (meters)
export Lz=1.0      # Length in z-direction (meters) (only used if dim=3)

export temp=268.15  # Initial temperature (Kelvin)
export grad_temp0X=0.01  # Temperature gradient in x-direction (K/m)
export grad_temp0Y=0   # Temperature gradient in y-direction (K/m)
export grad_temp0Z=0   # Temperature gradient in z-direction (K/m)

export dim=2   # Set dimension (2 for 2D, 3 for 3D)

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
    mpiexec -np $NUM_PROCS ./effective_k_ice -ksp_monitor -log_view
}


# ========== FUNCTION: MOVE OUTPUT FILES ==========
move_output_files() {
    if [ -d "$OUTPUT_DIR" ]; then
        # Move VTK and binary files only if they exist
        # mv *.vtk "$OUTPUT_DIR" 2>/dev/null || echo "No VTK files to move."
        mv *.bin "$OUTPUT_DIR" 2>/dev/null || echo "No binary files to move."
        echo "Output files moved to $OUTPUT_DIR"
    fi
}

# ========== EXECUTE WORKFLOW ==========
echo " "
echo "Starting Thermal Conductivity Simulation Workflow"
echo " "

compile_code
run_simulation
move_output_files

echo "Simulation complete."
