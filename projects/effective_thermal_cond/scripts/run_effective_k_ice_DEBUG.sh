#!/bin/zsh
###############################################################################
# run_effective_k_ice_DEBUG.sh
#
# This script is used to debug the effective_k_ice simulation.
# It sets a small mesh for quicker runs, compiles the code, and runs the
# simulation inside LLDB. Additional debugging options and post-processing
# steps are provided to help diagnose issues.
#
# Usage:
#   ./run_effective_k_ice_DEBUG.sh
#
# Note: Shell tracing is enabled (set -x) for detailed output. Remove or
# comment out "set -x" if you prefer less verbose output.
###############################################################################

# Enable shell command tracing for debugging purposes
set -x

###############################################################################
# 🔹 Set Environment Variables
###############################################################################
# Mesh and domain dimensions
export Nx=10                 # Number of elements in x-direction
export Ny=10                 # Number of elements in y-direction
export Nz=1                  # Number of elements in z-direction (only for 3D)

# Domain lengths (in meters)
export Lx=1.0                # Length in x-direction
export Ly=1.0                # Length in y-direction
export Lz=1.0                # Length in z-direction (only for 3D)

# Temperature settings
export temp=268.15           # Initial temperature (Kelvin)
export grad_temp0X=0.01      # Temperature gradient in x-direction (K/m)
export grad_temp0Y=0         # Temperature gradient in y-direction (K/m)
export grad_temp0Z=0         # Temperature gradient in z-direction (K/m)

# Simulation dimension: set to 3 for 3D, 2 for 2D
export dim=3

###############################################################################
# 🔹 Optional Output Settings
###############################################################################
export OUTPUT_VTK=1          # Enable VTK output if supported
export OUTPUT_BINARY=1       # Enable binary output

###############################################################################
# 🔹 MPI and Debugging Settings
###############################################################################
NUM_PROCS=1                  # Number of MPI processes to run (adjust as needed)

# Set PETSC_OPTIONS for debugging; these options instruct PETSc to start in
# debugger mode and output additional debug information.
export PETSC_OPTIONS="-start_in_debugger -malloc_debug -malloc_dump -log_view"

###############################################################################
# 🔹 Compilation Section
###############################################################################
echo "Starting Debug Simulation Workflow..."
echo "Compiling the code..."

# Clean and compile the code. Exits if compilation fails.
make clean && make effective_k_ice
if [ $? -ne 0 ]; then
    echo "Compilation failed. Exiting..."
    exit 1
fi
echo "Compilation successful."

###############################################################################
# 🔹 Running the Simulation in LLDB
###############################################################################
echo "Running effective_k_ice simulation with $NUM_PROCS process(es) in LLDB..."

# Run the simulation under LLDB with custom run arguments. You can modify the 
# LLDB commands as needed. Here, we set a breakpoint on the function "KSPSetFromOptions"
# and enable ksp_monitor.
lldb --batch --one-line "target create ./effective_k_ice" \
     --one-line "breakpoint set --name KSPSetFromOptions" \
     --one-line "settings set -- target.run-args -ksp_monitor" \
     --one-line "run"

###############################################################################
# 🔹 Optional Post-Processing for Debugging
###############################################################################
# Uncomment and modify the commands below to perform any additional debugging
# or post-processing on simulation output (e.g., parsing logs, generating plots).
# echo "Running post-processing..."
# python debug_postprocess.py "$OUTPUT_DIR"
# echo "Post-processing complete."

###############################################################################
# 🔹 macOS Crash Log Extraction (Optional)
###############################################################################
if [[ "$OSTYPE" == "darwin"* ]]; then
    echo "Checking for macOS crash logs..."
    log show --predicate 'process == "effective_k_ice"' --info --last 5m | tee crash_log.txt
fi

echo "Simulation completed. Check logs for debugging output."

# Disable shell tracing once done (optional)
set +x