#!/bin/bash
#SBATCH -J permafrost_60day_enclosed_pair_hum80
#SBATCH -A rubyfu
#SBATCH -t 0-12:00:00

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=30
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
#   $1 : PETSc options file (e.g., params.txt)
#   $2 : Title prefix for the run (used in folder name)
###############################################################################
params_file="${1:-params.txt}"
title="${2:-permafrost_}"

###############################################################################
# Create output folder based on timestamp and title
###############################################################################
create_folder() {
    local ts
    ts="$(date +%Y-%m-%d__%H.%M.%S)"

    # Base scratch dir on Resnick
    local base_dir="$SCRATCH/permafrost"

    # If not on Resnick, fall back to a local scratch under the repo
    if [[ ! -d "$base_dir" ]]; then
        base_dir="$(pwd)/scratch"
    fi

    # Differentiate SLURM vs local runs in the folder name
    if [[ "${SLURM_JOB_ID:-none}" != "none" ]]; then
        name="${title}${ts}_job${SLURM_JOB_ID}"
    else
        name="${title}${ts}_joblocal"
    fi

    folder="${base_dir}/${name}"
    mkdir -p "$folder"

    echo "Output folder: $folder"
}

###############################################################################
# Compile simulation code
###############################################################################
compile_code() {
    echo "Compiling..."
    make all
}

###############################################################################
# Compute optimal number of MPI ranks from Nx, Ny, Nz and 3 DoFs
# Strategy:
#   total_dofs = 3 * Nx * Ny * Nz
#   target_dofs_per_core ~ 10,000  (tweakable)
#   NPROCS = ceil(total_dofs / target_dofs_per_core)
#   Then cap by SLURM allocation (if present) or hardware (local).
###############################################################################
compute_optimal_nprocs() {
    local Nx Ny Nz
    local TARGET_DOFS_PER_CORE total_dofs

    # Default target ~1e4 DoFs per core
    TARGET_DOFS_PER_CORE=10000

    # Extract Nx, Ny, Nz from params_file; fall back to 1 if missing
    Nx=$(awk '$1=="-Nx"{print $2}' "$params_file" | head -n1)
    Ny=$(awk '$1=="-Ny"{print $2}' "$params_file" | head -n1)
    Nz=$(awk '$1=="-Nz"{print $2}' "$params_file" | head -n1)

    [[ -z "${Nx:-}" ]] && Nx=1
    [[ -z "${Ny:-}" ]] && Ny=1
    [[ -z "${Nz:-}" ]] && Nz=1

    # Integer arithmetic for total DoFs
    total_dofs=$((3 * Nx * Ny * Nz))

    # Ceil division: (total + TARGET - 1) / TARGET
    NPROCS=$(((total_dofs + TARGET_DOFS_PER_CORE - 1) / TARGET_DOFS_PER_CORE))
    (( NPROCS < 1 )) && NPROCS=1

    # Cap NPROCS by SLURM allocation (if in batch job)
    if [[ "${SLURM_JOB_ID:-none}" != "none" ]]; then
        local nnodes ntpn max_procs
        nnodes=${SLURM_NNODES:-1}
        ntpn=${SLURM_NTASKS_PER_NODE:-1}
        max_procs=$((nnodes * ntpn))
        if (( max_procs > 0 && NPROCS > max_procs )); then
            NPROCS=$max_procs
        fi
    else
        # Local run: cap by hardware cores (or 12 if nproc isn't available)
        local hw
        if command -v nproc >/dev/null 2>&1; then
            hw=$(nproc)
        else
            hw=12
        fi
        (( NPROCS > hw )) && NPROCS=$hw
    fi

    echo "------------------------------------------------------------"
    echo "Grid from opts: Nx=${Nx}, Ny=${Ny}, Nz=${Nz}"
    echo "Total DoFs:      3 * Nx * Ny * Nz = ${total_dofs}"
    echo "Target/core:     ${TARGET_DOFS_PER_CORE} DoFs"
    echo "Chosen NPROCS:   ${NPROCS}"
    echo "------------------------------------------------------------"
}

###############################################################################
# Run the simulation using MPI (srun on SLURM, mpiexec locally)
###############################################################################
run_simulation() {
    echo "Running simulation..."

    # Export the output folder so the code can see it as an environment variable
    export folder="$folder"

    if [[ "${SLURM_JOB_ID:-none}" != "none" ]]; then
        echo "SLURM_JOB_ID        = ${SLURM_JOB_ID}"
        echo "SLURM_NNODES        = ${SLURM_NNODES:-unknown}"
        echo "SLURM_NTASKS        = ${SLURM_NTASKS:-unknown}"
        echo "SLURM_CPUS_PER_TASK = ${SLURM_CPUS_PER_TASK:-unknown}"
        echo "Using NPROCS        = ${NPROCS}"

        srun -n "${NPROCS}" ./permafrost \
            -options_file "$params_file" \
            -output_path "$folder" \
            -snes_rtol 1e-3 \
            -snes_stol 1e-6 \
            -snes_max_it 7 \
            -ksp_gmres_restart 150 \
            -ksp_max_it 1000 \
            -ksp_converged_reason \
            -snes_converged_reason \
            -snes_linesearch_monitor \
            -snes_linesearch_type basic | tee "$folder/outp.txt"
    else
        echo "No SLURM environment detected; running locally with mpiexec."
        echo "Using NPROCS = ${NPROCS}"

        mpiexec -n "${NPROCS}" ./permafrost \
            -options_file "$params_file" \
            -output_path "$folder" \
            -snes_rtol 1e-3 \
            -snes_stol 1e-6 \
            -snes_max_it 7 \
            -ksp_gmres_restart 150 \
            -ksp_max_it 1000 \
            -ksp_converged_reason \
            -snes_converged_reason \
            -snes_linesearch_monitor \
            -snes_linesearch_type basic | tee "$folder/outp.txt"
    fi
}

################################################################################
# Copy relevant scripts to folder and save summary parameters
################################################################################
finalize_results() {
    echo "Finalizing results..."
    # Copy some scripts and source for provenance
    cp ./scripts/run_permafrost.sh "$folder" 2>/dev/null || true
    cp ./scripts/plotpermafrost.py ./scripts/plotSSA.py ./scripts/plotPorosity.py "$folder" 2>/dev/null || true
    cp ./src/permafrost.c "$folder" 2>/dev/null || true

    # Save a copy of the PETSc options file used for this run
    cp "$params_file" "$folder/$(basename "$params_file")"
}

################################################################################
# Run post-processing plotting script
################################################################################
run_plotting() {
    echo "Queuing plotpermafrost.py"
    # If you have a separate HPC plotting wrapper, call it here.
    # For now, just call the existing script if it exists.
    if [[ -x ./scripts/run_plotpermafrost.sh ]]; then
        ./scripts/run_plotpermafrost.sh "$name"
    fi
}

################################################################################
# USER-DEFINED SIMULATION SETTINGS (currently unused but kept for reference)
################################################################################
echo " "
echo "Starting permafrost simulation workflow"
echo " "

compile_code
create_folder
compute_optimal_nprocs
run_simulation
finalize_results
run_plotting

echo "-------------------------------------------------------------------------"
echo " "
echo "✅ Done with permafrost simulation!"
echo "-------------------------------------------------------------------------"
echo " "