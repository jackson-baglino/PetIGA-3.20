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

# Number of MPI processes
NUM_PROCS=1

echo "Starting Thermal Conductivity Simulation Workflow"

# ========== COMPILATION ==========
echo "Compiling the code..."
make clean && make effective_k_ice
if [ $? -ne 0 ]; then
    echo "Compilation failed. Exiting..."
    exit 1
fi

echo "Compilation successful."

# ========== RUN SIMULATION IN LLDB ==========
echo "Running effective_k_ice simulation with $NUM_PROCS processes in LLDB..."

# ./effective_k_ice -snes_monitor -snes_converged_reason -ksp_monitor -log_view

lldb --batch --one-line "target create ./effective_k_ice" \
     --one-line "settings set -- target.run-args -ksp_monitor" \
     --one-line "run"