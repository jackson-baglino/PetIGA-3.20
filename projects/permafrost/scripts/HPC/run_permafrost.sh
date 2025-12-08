#!/bin/bash
#SBATCH -J permafrost_Tm-20_hum98
#SBATCH -A rubyfu
#SBATCH -t 0-4:00:00

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=50
#SBATCH --cpus-per-task=1
#SBATCH --partition=expansion
#SBATCH --mem-per-cpu=1G

#SBATCH -o "output_files/%x.o%j"
#SBATCH -e "output_files/%x.e%j"
#SBATCH --mail-user=jbaglino@caltech.edu
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --constraint='icelake|skylake|cascadelake'

set -euo pipefail
trap 'echo "❌ Error on line $LINENO"; exit 1' ERR

###############################################################################
# Input arguments
#   $1 : PETSc options file (e.g., path/to/permafrost_test1.opts)
#   $2 : Title prefix for the run (used in folder name)
###############################################################################
params_file="${1:-params.txt}"
title="${2:-permafrost_}"

mkdir -p output_files

###############################################################################
# Create output folder based on timestamp and title
###############################################################################
create_folder() {
    local ts
    ts="$(date +%Y-%m-%d__%H.%M.%S)"

    if [[ -n "${SLURM_JOB_ID:-}" ]]; then
        # On Resnick / in a SLURM allocation
        out_root="/resnick/scratch/${USER}/permafrost"
        name="${title}${ts}_job${SLURM_JOB_ID}"
    else
        # Local / interactive run (not via sbatch)
        out_root="${PWD}/scratch"
        name="${title}${ts}_joblocal"
    fi

    folder="${out_root}/${name}"

    echo "Output folder: ${folder}"
    mkdir -p "${folder}"
}

###############################################################################
# Compile simulation code
###############################################################################
compile_code() {
    echo "Compiling..."
    make all
}

###############################################################################
# Run the simulation using MPI (or serial if needed)
###############################################################################
run_simulation() {
    echo "Running simulation..."

    if [[ -n "${SLURM_JOB_ID:-}" ]]; then
        NPROCS="${SLURM_NTASKS:-1}"
    else
        NPROCS=12
    fi

    echo "SLURM_JOB_ID        = ${SLURM_JOB_ID:-none}"
    echo "SLURM_NNODES        = ${SLURM_NNODES:-none}"
    echo "SLURM_NTASKS        = ${SLURM_NTASKS:-none}"
    echo "SLURM_CPUS_PER_TASK = ${SLURM_CPUS_PER_TASK:-none}"
    echo "Using NPROCS        = ${NPROCS}"

    export folder

    if [[ -n "${SLURM_JOB_ID:-}" ]]; then
        ###################################################################
        # Inside a SLURM allocation → use srun
        ###################################################################
        srun -n "${NPROCS}" ./permafrost \
            -options_file "${params_file}" \
            -output_path "${folder}" \
            -snes_rtol 1e-3 \
            -snes_stol 1e-6 \
            -snes_max_it 7 \
            -ksp_gmres_restart 150 \
            -ksp_max_it 1000 \
            -ksp_converged_reason \
            -snes_converged_reason \
            -snes_linesearch_monitor \
            -snes_linesearch_type basic | tee "${folder}/outp.txt"
    else
        ###################################################################
        # Not in SLURM → try mpiexec, otherwise run serial
        ###################################################################
        if command -v mpiexec >/dev/null 2>&1; then
            echo "mpiexec found; running parallel locally with ${NPROCS} ranks."
            mpiexec -n "${NPROCS}" ./permafrost \
                -options_file "${params_file}" \
                -output_path "${folder}" \
                -snes_rtol 1e-3 \
                -snes_stol 1e-6 \
                -snes_max_it 7 \
                -ksp_gmres_restart 150 \
                -ksp_max_it 1000 \
                -ksp_converged_reason \
                -snes_converged_reason \
                -snes_linesearch_monitor \
                -snes_linesearch_type basic | tee "${folder}/outp.txt"
        else
            echo "⚠️  mpiexec not found; running ./permafrost in serial."
            NPROCS=1
            ./permafrost \
                -options_file "${params_file}" \
                -output_path "${folder}" \
                -snes_rtol 1e-3 \
                -snes_stol 1e-6 \
                -snes_max_it 7 \
                -ksp_gmres_restart 150 \
                -ksp_max_it 1000 \
                -ksp_converged_reason \
                -snes_converged_reason \
                -snes_linesearch_monitor \
                -snes_linesearch_type basic | tee "${folder}/outp.txt"
        fi
    fi
}

###############################################################################
# Copy relevant scripts and options to folder
###############################################################################
finalize_results() {
    echo "Finalizing results..."

    # Copy this script into the folder
    cp "$0" "${folder}/run_permafrost.sh" 2>/dev/null || true

    # Copy plotting scripts if they exist
    if [[ -d scripts ]]; then
        cp scripts/run_plotpermafrost.sh scripts/plotpermafrost.py \
           scripts/plotSSA.py scripts/plotPorosity.py "${folder}" 2>/dev/null || true
    fi

    # Copy main C file (for provenance)
    if [[ -f src/permafrost.c ]]; then
        cp src/permafrost.c "${folder}"
    fi

    # Save a copy of the PETSc options file used for this run
    cp "${params_file}" "${folder}/$(basename "${params_file}")"
}

###############################################################################
# Optionally queue plotting
###############################################################################
run_plotting() {
    echo "Queuing plotpermafrost.py (if desired)"
    if [[ -x ./scripts/run_plotpermafrost.sh ]]; then
        ./scripts/run_plotpermafrost.sh "${folder}"
    fi
}

###############################################################################
# Main workflow
###############################################################################
echo
echo "Starting permafrost simulation workflow"
echo

compile_code
create_folder
run_simulation
finalize_results
run_plotting

echo "-------------------------------------------------------------------------"
echo
echo "✅ Done with permafrost simulation!"
echo "Output folder: ${folder}"
echo "-------------------------------------------------------------------------"
echo