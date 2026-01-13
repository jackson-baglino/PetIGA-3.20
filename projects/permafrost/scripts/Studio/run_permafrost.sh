set -euo pipefail
trap 'echo "❌ Error on line $LINENO"; exit 1' ERR

###############################################################################
# Input arguments
#   $1 : PETSc options file (e.g., params.txt)
#   $2 : Title prefix for the run (used in folder name)
###############################################################################
params_file="${1:-params.txt}"
title="${2:-permafrost_}"

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

    mpiexec -np 12 ./permafrost \
        -options_file "$params_file" \
        -output_path "$folder" \
        -snes_rtol 1e-3 \
        -snes_stol 1e-2 \
        -snes_atol 1e-6 \
        -snes_max_it 8 \
        -ksp_gmres_restart 150 \
        -ksp_max_it 1000 \
        -ksp_rtol 1e-10 \
        -ksp_converged_reason \
        -snes_converged_reason \
        -snes_linesearch_monitor \
        -snes_linesearch_type basic | tee "$folder/outp.txt"
}

################################################################################
# Copy relevant scripts to folder and save summary parameters to .dat and CSV
################################################################################
finalize_results() {
    echo "Finalizing results..."
    cd ./scripts
    cp Studio/run_permafrost.sh plotpermafrost.py plotSSA.py plotPorosity.py $folder
    cd ../src
    cp permafrost.c $folder
    cd ../

    # Save a copy of the PETSc options file used for this run
    cp "$params_file" "$folder/$(basename "$params_file")"
}

################################################################################
# Run post-processing plotting script
################################################################################
run_plotting() {
    echo "Queuing plotpermafrost.py"
    ./scripts/run_plotpermafrost.sh $name
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

run_simulation

finalize_results

run_plotting

echo "-------------------------------------------------------------------------"
echo " "
echo "✅ Done with permafrost simulation!"
echo "-------------------------------------------------------------------------"
echo " "
