
set -uo pipefail

###############################################################################
# Usage / argument checking
###############################################################################
usage() {
    echo ""
    echo "Usage:"
    echo "  ./run_permafrost.sh <PETSc_options_file> <run_name_prefix>"
    echo ""
    echo "Arguments:"
    echo "  PETSc_options_file   Path to PETSc .opts file used to configure the run"
    echo "  run_name_prefix      Prefix for the output folder name"
    echo ""
    echo "Example:"
    echo "  ./run_permafrost.sh ./inputs/permafrost_test.opts capillary_bridge_"
    echo ""
}

if [ "$#" -lt 2 ]; then
    echo "❌ Error: Missing required command-line arguments."
    usage
    exit 1
fi

# Don't exit immediately on non-zero status; we want to keep going (e.g., run plotting)
# even if the simulation fails.
trap 'echo "❌ Script error on line $LINENO"' ERR

###############################################################################
# Input arguments
#   $1 : PETSc options file (e.g., params.txt)
#   $2 : Title prefix for the run (used in folder name)
###############################################################################
params_file="${1:-params.txt}"
title="${2:-permafrost_TESTING_}"

###############################################################################
# Create output folder based on timestamp and title
################################################################################
create_folder() {
    name="$title$(date +%Y-%m-%d__%H.%M.%S)"
    dir="/Users/jacksonbaglino/SimulationResults/permafrost/scratch"
    folder="$dir/$name"

    if [ ! -d "$dir" ]; then
        mkdir -p "$dir"
    fi

    mkdir -p "$folder"
}

################################################################################
# Copy plotting / helper scripts to the output folder (before the simulation)
################################################################################
stage_output_folder() {
    echo "Staging scripts into output folder..."

    # Copy plotting and runner scripts
    cd ./scripts
    cp Studio/run_permafrost.sh plotpermafrost.py plotSSA.py plotPorosity.py "$folder"
    cd ../

    # Save a copy of the PETSc options file used for this run
    cp "$params_file" "$folder/$(basename "$params_file")"
}

###############################################################################
# Compile simulation code
################################################################################
compile_code() {
    echo "Compiling..."
    make all
}

###############################################################################
# Run the simulation using MPI
###############################################################################
run_simulation() {
    echo "Running simulation..."

    # Export the output folder so the code can see it as an environment variable
    export folder="$folder"

    # Run the simulation. If it fails, capture the exit code but continue.
    sim_exit=0
    set +e
    mpiexec -np 12 ./permafrost \
        -options_file "$params_file" \
        -output_path "$folder" \
        | tee "$folder/outp.txt"

    sim_exit=${PIPESTATUS[0]}
    set -e

    if [ "$sim_exit" -ne 0 ]; then
        echo "⚠️ Simulation exited with code $sim_exit (continuing to post-processing)"
    fi
}

################################################################################
# Copy relevant scripts to folder and save summary parameters to .dat and CSV
################################################################################
finalize_results() {
    echo "Finalizing results..."

    # Save a full copy of the source code for reproducibility
    if [ -d "./src" ]; then
        cp -r ./src "$folder/"
    else
        echo "⚠️ Warning: src directory not found; source code not copied."
    fi
}

################################################################################
# Run post-processing plotting script
################################################################################
run_plotting() {
    echo "Queuing plotpermafrost.py"
    set +e
    ./scripts/run_plotpermafrost.sh "$name"
    plot_exit=$?
    set -e
    if [ "$plot_exit" -ne 0 ]; then
        echo "⚠️ Plotting script exited with code $plot_exit"
    fi
}

################################################################################
# USER-DEFINED SIMULATION SETTINGS
################################################################################
echo " "
echo "Starting permafrost simulation workflow"
echo " "

delt_t=1.0e-4
t_final=60*24*3600  # 60 days in seconds
n_out=100
t_final=$(echo "$t_final" | bc -l)

if (( $(echo "$t_final <= 0" | bc -l) )); then
    echo "❌ Error: t_final must be > 0. Got $t_final"
    exit 1
fi

compile_code

create_folder

stage_output_folder

run_simulation

finalize_results

run_plotting

# Exit with the simulation exit code if it failed, otherwise 0.
exit ${sim_exit:-0}

echo "-------------------------------------------------------------------------"
echo " "
echo "✅ Done with permafrost simulation!"
echo "-------------------------------------------------------------------------"
echo " "
