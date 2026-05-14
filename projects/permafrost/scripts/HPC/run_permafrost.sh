#!/bin/bash
#SBATCH -J permafrost_SlabAndGrains
#SBATCH -A rubyfu
#SBATCH -t 0-24:00:00

#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
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
#   $1 : PETSc options file (e.g., inputs/tests/test3_EnclosedGrainPair.opts)
#   $2 : Title prefix for the run (used in folder name)
###############################################################################
params_file="${1:-inputs/tests/test_2D_TouchingGrainPair.opts}"
title="${2:-}"

###############################################################################
# Resolve project root.
#
# When submitted via `sbatch`, SLURM copies the script to a temp path under
# /var/spool/..., so BASH_SOURCE[0] no longer points to the source tree.
# SLURM exports SLURM_SUBMIT_DIR = the directory where `sbatch` was called,
# which must be the project root.  Fall back to BASH_SOURCE resolution when
# running the script directly (local or interactive HPC session).
###############################################################################
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    PROJECT_ROOT="$SLURM_SUBMIT_DIR"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
fi

if [ ! -f "$PROJECT_ROOT/makefile" ] && [ ! -f "$PROJECT_ROOT/Makefile" ]; then
    echo "❌ Error: Could not find makefile in resolved project root: $PROJECT_ROOT"
    echo "   Run sbatch from the project root directory."
    exit 1
fi

# Load cost utilities (graceful no-op if missing)
_HPC_COST_SH="$PROJECT_ROOT/scripts/HPC/hpc_cost.sh"
if [[ -f "$_HPC_COST_SH" ]]; then
    # shellcheck source=scripts/HPC/hpc_cost.sh
    source "$_HPC_COST_SH"
else
    hpc_cost_post_job() { :; }
fi

EXEC="$PROJECT_ROOT/permafrost"
SRC_DIR="$PROJECT_ROOT/src"
SCRIPTS_DIR="$PROJECT_ROOT/scripts"
INPUTS_DIR="$PROJECT_ROOT/inputs"
UNIVERSAL_OPTS="$INPUTS_DIR/universal.opts"

# Discover mpiexec: prefer the PETSC-bundled binary (always present when
# PETSc downloaded MPICH), fall back to whatever is on PATH.
MPIEXEC=""
if [[ -n "${PETSC_DIR:-}" && -n "${PETSC_ARCH:-}" ]]; then
    _petsc_mpi="$PETSC_DIR/$PETSC_ARCH/bin/mpiexec"
    [[ -x "$_petsc_mpi" ]] && MPIEXEC="$_petsc_mpi"
fi
if [[ -z "$MPIEXEC" ]]; then
    MPIEXEC=$(command -v mpiexec 2>/dev/null || command -v mpirun 2>/dev/null || echo "")
fi
if [[ -z "$MPIEXEC" ]]; then
    echo "❌ Error: mpiexec/mpirun not found."
    echo "   Set PETSC_DIR and PETSC_ARCH, or add mpiexec to PATH."
    exit 1
fi

# Resolve params_file relative to project root if not absolute
if [[ "$params_file" != /* ]]; then
    params_file="$PROJECT_ROOT/$params_file"
fi

if [ ! -f "$params_file" ]; then
    echo "❌ Error: Options file not found: $params_file"
    exit 1
fi

# Global state
folder=""
name=""
sim_exit=0

###############################################################################
# Compile simulation code
###############################################################################
compile_code() {
    echo ""
    echo "--- Compiling ---"
    echo "Project root: $PROJECT_ROOT"

    cd "$PROJECT_ROOT"
    make clean && make all

    if [ ! -f "$EXEC" ]; then
        echo "❌ Executable not found after compilation: $EXEC"
        exit 1
    fi
    echo "✅ Compilation successful."
}

###############################################################################
# derive_ic_subfolder — maps -ic_type to a clean category folder name
###############################################################################
derive_ic_subfolder() {
    local ic
    ic=$(awk '$1=="-ic_type"{print $2}' "$params_file" | head -n1)
    case "${ic:-}" in
        ice_slab)        echo "IceSlab" ;;
        enclosed)        echo "EnclosedGrainPair" ;;
        contact_sed)     echo "ContactSed" ;;
        capillary)       echo "CapillaryBridge" ;;
        slab_and_grains) echo "SlabAndGrains" ;;
        ice_cap)         echo "IceCap" ;;
        single_ice)      echo "SingleIceGrain" ;;
        ice_sed_pair)    echo "IceSedPair" ;;
        *)               echo "Other" ;;
    esac
}

###############################################################################
# Create output folder: $SCRATCH/permafrost/<ic_category>/<opts_name>[_tag]_<ts>
###############################################################################
create_folder() {
    echo ""
    echo "--- Creating output folder ---"

    local subfolder ts opts_name tag job_suffix
    subfolder=$(derive_ic_subfolder)
    opts_name=$(basename "$params_file" .opts)
    ts="$(date +%Y-%m-%d__%H.%M.%S)"
    tag="${title:+_${title}}"

    if [[ "${SLURM_JOB_ID:-none}" != "none" ]]; then
        job_suffix="_job${SLURM_JOB_ID}"
    else
        job_suffix=""
    fi

    name="${opts_name}${tag}_${ts}${job_suffix}"

    # Base scratch dir on Resnick; fall back to local scratch
    local base_dir
    if [[ -d "${SCRATCH:-}" ]]; then
        base_dir="$SCRATCH/permafrost"
    else
        base_dir="$PROJECT_ROOT/scratch"
    fi

    folder="${base_dir}/${subfolder}/${name}"
    mkdir -p "$folder"
    echo "Output folder: $folder"
}

###############################################################################
# Stage output folder — copy opts files, run script, and post-processing
###############################################################################
stage_output_folder() {
    echo ""
    echo "--- Staging output folder ---"

    # Opts files and this run script
    [ -f "$UNIVERSAL_OPTS" ] && cp "$UNIVERSAL_OPTS" "$folder/"
    cp "$params_file" "$folder/$(basename "$params_file")"
    cp "${BASH_SOURCE[0]}" "$folder/run_permafrost.sh"

    # Post-processing scripts — copy the full postprocess/ directory so that
    # results can be reproduced without access to the source tree
    local POSTPROCESS="$PROJECT_ROOT/postprocess"
    if [ -d "$POSTPROCESS" ]; then
        cp -r "$POSTPROCESS" "$folder/postprocess"
        echo "  Copied postprocess/ → $folder/postprocess/"
    else
        echo "⚠️  postprocess/ directory not found at $POSTPROCESS — skipping."
    fi

    echo "✅ Staging complete."
}

###############################################################################
# Compute optimal number of MPI ranks from Nx, Ny, Nz and 4 DoFs
###############################################################################
compute_optimal_nprocs() {
    local Nx Ny Nz
    local TARGET_DOFS_PER_CORE total_dofs

    TARGET_DOFS_PER_CORE=10000

    Nx=$(awk '$1=="-Nx"{print $2}' "$params_file" | head -n1)
    Ny=$(awk '$1=="-Ny"{print $2}' "$params_file" | head -n1)
    Nz=$(awk '$1=="-Nz"{print $2}' "$params_file" | head -n1)

    [[ -z "${Nx:-}" ]] && Nx=1
    [[ -z "${Ny:-}" ]] && Ny=1
    [[ -z "${Nz:-}" ]] && Nz=1

    total_dofs=$((4 * Nx * Ny * Nz))

    NPROCS=$(((total_dofs + TARGET_DOFS_PER_CORE - 1) / TARGET_DOFS_PER_CORE))
    (( NPROCS < 1 )) && NPROCS=1

    if [[ "${SLURM_JOB_ID:-none}" != "none" ]]; then
        local nnodes ntpn max_procs
        nnodes=${SLURM_NNODES:-1}
        ntpn=${SLURM_NTASKS_PER_NODE:-1}
        max_procs=$((nnodes * ntpn))
        (( max_procs > 0 && NPROCS > max_procs )) && NPROCS=$max_procs
    else
        local hw
        hw=$(nproc 2>/dev/null || echo 12)
        (( NPROCS > hw )) && NPROCS=$hw
    fi

    echo "------------------------------------------------------------"
    echo "Grid from opts: Nx=${Nx}, Ny=${Ny}, Nz=${Nz}"
    echo "Total DoFs:     4 * Nx * Ny * Nz = ${total_dofs}"
    echo "Target/core:    ${TARGET_DOFS_PER_CORE} DoFs"
    echo "Chosen NPROCS:  ${NPROCS}"
    echo "------------------------------------------------------------"
}

###############################################################################
# Run the simulation using MPI (srun on SLURM, mpiexec locally)
###############################################################################
run_simulation() {
    echo ""
    echo "--- Running simulation ---"
    echo "Executable   : $EXEC"
    echo "Universal    : $UNIVERSAL_OPTS"
    echo "Options file : $params_file"
    echo "Output path  : $folder"

    export folder

    local universal_arg=""
    [ -f "$UNIVERSAL_OPTS" ] && universal_arg="-options_file $UNIVERSAL_OPTS"

    set +e

    if [[ "${SLURM_JOB_ID:-none}" != "none" ]]; then
        echo "SLURM_JOB_ID        = ${SLURM_JOB_ID}"
        echo "SLURM_NNODES        = ${SLURM_NNODES:-unknown}"
        echo "SLURM_NTASKS        = ${SLURM_NTASKS:-unknown}"
        echo "SLURM_CPUS_PER_TASK = ${SLURM_CPUS_PER_TASK:-unknown}"
        echo "Using NPROCS        = ${NPROCS}"

        srun -n "${NPROCS}" "$EXEC" \
            $universal_arg \
            -options_file "$params_file" \
            -output_path  "$folder" \
            | tee "$folder/outp.txt"
    else
        echo "No SLURM environment detected; running locally with mpiexec."
        echo "Using NPROCS = ${NPROCS}"
        echo "MPIEXEC      = ${MPIEXEC}"

        "$MPIEXEC" -n "${NPROCS}" "$EXEC" \
            $universal_arg \
            -options_file "$params_file" \
            -output_path  "$folder" \
            | tee "$folder/outp.txt"
    fi

    sim_exit=${PIPESTATUS[0]}
    set -e

    if [ "$sim_exit" -ne 0 ]; then
        echo "⚠️  Simulation exited with code $sim_exit (continuing to post-processing)"
    else
        echo "✅ Simulation completed successfully."
    fi
}

###############################################################################
# Copy source files for reproducibility
###############################################################################
copy_source_code() {
    echo ""
    echo "--- Copying source code ---"

    local src_dest="$folder/src"
    mkdir -p "$src_dest"
    local found=0

    if [ -d "$SRC_DIR" ]; then
        for ext in "*.c" "*.h" "*.cpp" "*.hpp"; do
            for f in "$SRC_DIR"/$ext; do
                [ -f "$f" ] && cp "$f" "$src_dest/" && found=1
            done
        done
    fi

    [ -d "$PROJECT_ROOT/include" ] && cp -r "$PROJECT_ROOT/include" "$src_dest/" && found=1

    for mf in makefile Makefile; do
        if [ -f "$PROJECT_ROOT/$mf" ]; then
            cp "$PROJECT_ROOT/$mf" "$src_dest/"
            found=1
            break
        fi
    done

    [ "$found" -eq 0 ] && echo "⚠️  Warning: No source files found." \
                       || echo "✅ Source code copied to $src_dest"
}

# Post-processing is intentionally omitted on HPC.
# All output files are staged in $folder so the user can run
#   bash <run_folder>/postprocess/run_postprocess.sh
# locally after rsyncing the results.

###############################################################################
# Main workflow
###############################################################################
JOB_START_TIME=$(date +%s)

echo ""
echo "========================================================================="
echo "  Permafrost simulation workflow (HPC)"
echo "  Project root : $PROJECT_ROOT"
echo "  Options      : $params_file"
echo "  Title        : $title"
echo "========================================================================="

compile_code
create_folder
stage_output_folder
compute_optimal_nprocs
run_simulation
copy_source_code

# Post-processing is intentionally skipped on HPC.
# All required files (igasol.dat, sol_*.dat, SSA_evo.dat, outp.txt, postprocess/)
# are staged in $folder. After rsyncing to your local machine, run:
#   bash <run_folder>/postprocess/run_postprocess.sh

echo ""
echo "========================================================================="
if [ "${sim_exit:-0}" -ne 0 ]; then
    echo "  ⚠️  Workflow completed with simulation errors (exit code $sim_exit)"
else
    echo "  ✅ Workflow completed successfully"
fi
echo "  Results: $folder"
echo "========================================================================="
echo ""

hpc_cost_post_job \
    "${NPROCS:-1}" \
    "$(( $(date +%s) - JOB_START_TIME ))" \
    "${SLURM_JOB_ID:-}"

exit "${sim_exit:-0}"
