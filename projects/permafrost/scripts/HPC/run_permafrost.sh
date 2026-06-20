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
#   $1   : geometry name (e.g. 2D_touching_grains)
#   $2   : experiment name (e.g. 1day_T-20_h1.00)
#   $3   : Title prefix for the run (used in folder name)
#   $4.. : extra options forwarded verbatim to the permafrost executable,
#          appended after the three -options_file flags so they override
#          anything set in solver.opts/geometry/experiment opts
#          (e.g. -d0_GT 1.0e-8). Populated by submit_permafrost.sh's `--`
#          separator.
###############################################################################
geom_name="${1:-2D_touching_grains}"
exp_name="${2:-1day_T-20_h1.00}"
title="${3:-}"
EXTRA_OPTS=("${@:4}")

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
SOLVER_OPTS="$INPUTS_DIR/solver.opts"
GEOM_OPTS="$INPUTS_DIR/geometry/${geom_name}.opts"
EXP_OPTS="$INPUTS_DIR/experiment/${exp_name}.opts"

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

if [ ! -f "$SOLVER_OPTS" ]; then
    echo "❌ Error: Solver opts not found: $SOLVER_OPTS"
    exit 1
fi
if [ ! -f "$GEOM_OPTS" ]; then
    echo "❌ Error: Geometry opts not found: $GEOM_OPTS"
    echo "   Pass a geometry from inputs/geometry/ as \$1 (without .opts)."
    exit 1
fi
if [ ! -f "$EXP_OPTS" ]; then
    echo "❌ Error: Experiment opts not found: $EXP_OPTS"
    echo "   Pass an experiment from inputs/experiment/ as \$2 (without .opts)."
    exit 1
fi

# For backward-compat with the rest of the script (Nx/Ny/Nz lookups, etc.)
params_file="$GEOM_OPTS"

# Global state
folder=""
name=""
sim_exit=0

###############################################################################
# Compile simulation code
###############################################################################
compile_code() {
    echo ""
    # Skip-compile path: required when multiple parallel jobs share the same
    # project root, since `make clean && make all` in shared obj/ races between
    # jobs and produces "Stale file handle" / "No space left on device" errors
    # when one job's `make clean` deletes obj/*.o while another job is mid-write.
    # Set SKIP_COMPILE=1 (e.g. via `sbatch --export=ALL,SKIP_COMPILE=1`) when the
    # binary has already been built by the submitter before fanning out jobs.
    if [[ "${SKIP_COMPILE:-0}" == "1" ]]; then
        echo "--- Skipping compile (SKIP_COMPILE=1) ---"
        if [ ! -f "$EXEC" ]; then
            echo "❌ SKIP_COMPILE=1 but executable not found: $EXEC"
            echo "   Build first on the submission host, then resubmit."
            exit 1
        fi
        echo "Using existing executable: $EXEC"
        return
    fi

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
# Create output folder:
#   $SCRATCH/permafrost/<geometry>/<timestamp>_<experiment>[_<tag>][_job<id>]
#
# One subfolder per distinct geometry (geom_name itself, not the old
# -ic_type category bucket, which lumped unrelated geometries that happen to
# share an IC type together -- e.g. every multi_grains-based geometry landed
# in the same "Other" folder). The timestamp leads the run-folder name so a
# plain `ls` (alphabetical) sort is also chronological -- no need to parse a
# timestamp out of the middle of a long name to find the most recent run.
###############################################################################
create_folder() {
    echo ""
    echo "--- Creating output folder ---"

    # Batch-mode path: when submit_batch.sh fans out a sweep, every job sets
    # BATCH_OUT_DIR=$SCRATCH/permafrost/batch_<ts>[_<tag>]/, and we just write
    # into a clean `<geom>__<exp>/` subfolder of that shared parent so the
    # whole batch can be downloaded with a single rsync.
    if [[ -n "${BATCH_OUT_DIR:-}" ]]; then
        name="${geom_name}__${exp_name}"
        folder="${BATCH_OUT_DIR}/${name}"
        mkdir -p "$folder"
        echo "Output folder (batch mode): $folder"
        return
    fi

    # Single-run path: one subfolder per geometry, timestamp-led run names.
    local ts tag job_suffix
    ts="$(date +%Y-%m-%d__%H.%M.%S)"
    tag="${title:+_${title}}"

    if [[ "${SLURM_JOB_ID:-none}" != "none" ]]; then
        job_suffix="_job${SLURM_JOB_ID}"
    else
        job_suffix=""
    fi

    name="${ts}_${exp_name}${tag}${job_suffix}"

    # Base scratch dir on Resnick; fall back to local scratch
    local base_dir
    if [[ -d "${SCRATCH:-}" ]]; then
        base_dir="$SCRATCH/permafrost"
    else
        base_dir="$PROJECT_ROOT/scratch"
    fi

    folder="${base_dir}/${geom_name}/${name}"
    mkdir -p "$folder"
    echo "Output folder: $folder"
}

###############################################################################
# Stage output folder — copy opts files, run script, and post-processing
###############################################################################
stage_output_folder() {
    echo ""
    echo "--- Staging output folder ---"

    # Opts files (all three) and this run script
    cp "$SOLVER_OPTS" "$folder/"
    cp "$GEOM_OPTS"   "$folder/"
    cp "$EXP_OPTS"    "$folder/"
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
# Compute optimal number of MPI ranks from Nx, Ny, Nz and dof
###############################################################################
compute_optimal_nprocs() {
    local Nx Ny Nz dof grid_src
    local TARGET_DOFS_PER_CORE total_dofs

    TARGET_DOFS_PER_CORE=10000

    dof=$(awk '$1=="-dof"{print $2}' "$SOLVER_OPTS" | head -n1)
    [[ -z "${dof:-}" ]] && dof=4

    # -geom_file meshes override -Nx/-Ny/-Nz; read the DOF count from the
    # "# DOF_GRID: nx ny [nz]" comment written alongside -geom_file in the
    # geometry opts -- otherwise Nx/Ny default to 1 and NPROCS collapses to 1.
    if grep -q "^-geom_file" "$params_file"; then
        read -r Nx Ny Nz <<< "$(awk '$1=="#" && $2=="DOF_GRID:"{print $3, $4, $5}' "$params_file" | head -n1)"
        grid_src="DOF_GRID comment"
    else
        Nx=$(awk '$1=="-Nx"{print $2}' "$params_file" | head -n1)
        Ny=$(awk '$1=="-Ny"{print $2}' "$params_file" | head -n1)
        Nz=$(awk '$1=="-Nz"{print $2}' "$params_file" | head -n1)
        grid_src="-Nx/-Ny/-Nz flags"
    fi

    [[ -z "${Nx:-}" ]] && Nx=1
    [[ -z "${Ny:-}" ]] && Ny=1
    [[ -z "${Nz:-}" ]] && Nz=1

    total_dofs=$((dof * Nx * Ny * Nz))

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
    echo "Grid source:    ${grid_src}"
    echo "Total DoFs:     ${dof} × Nx × Ny × Nz = ${total_dofs}"
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
    echo "PROJECT_ROOT : $PROJECT_ROOT"
    echo "Executable   : $EXEC"
    echo "Solver       : $SOLVER_OPTS"
    echo "Geometry     : $GEOM_OPTS"
    echo "Experiment   : $EXP_OPTS"
    echo "Output path  : $folder"
    echo "Extra opts   : ${EXTRA_OPTS[*]:-(none)}"

    # Validate the -geom_file .dat so stale/missing mesh files are caught
    # before launching hundreds of MPI ranks.
    if grep -q "^-geom_file" "$GEOM_OPTS"; then
        local rel_dat abs_dat dat_bytes
        rel_dat=$(awk '$1=="-geom_file"{print $2}' "$GEOM_OPTS" | head -n1)
        abs_dat="$PROJECT_ROOT/$rel_dat"
        echo "Geom .dat    : $abs_dat"
        if [ ! -f "$abs_dat" ]; then
            echo "❌ ERROR: geom_file not found: $abs_dat"
            echo "   Run: python3 preprocess/build_geometry_multi_grain.py \\"
            echo "          --out $rel_dat"
            exit 1
        fi
        dat_bytes=$(wc -c < "$abs_dat")
        echo "Geom .dat size: ${dat_bytes} bytes"
        # A 122x122/P=2 mesh is ~371 KB; 64x64/P=2 is ~100 KB.
        # Warn if the file looks too small for the DOF_GRID in this opts.
        local dof_nx dof_ny
        read -r dof_nx dof_ny <<< "$(awk '$1=="#" && $2=="DOF_GRID:"{print $3,$4}' "$GEOM_OPTS" | head -n1)"
        local min_bytes=$(( ${dof_nx:-0} * ${dof_ny:-0} * 8 ))
        if (( dat_bytes < min_bytes )); then
            echo "⚠️  WARNING: geom .dat is ${dat_bytes} bytes but DOF_GRID ${dof_nx}x${dof_ny}"
            echo "   suggests at least ${min_bytes} bytes. The mesh file may be stale."
            echo "   Regenerate with: python3 preprocess/build_geometry_multi_grain.py \\"
            echo "     --out $rel_dat"
        fi
    fi

    export folder

    set +e

    if [[ "${SLURM_JOB_ID:-none}" != "none" ]]; then
        echo "SLURM_JOB_ID        = ${SLURM_JOB_ID}"
        echo "SLURM_NNODES        = ${SLURM_NNODES:-unknown}"
        echo "SLURM_NTASKS        = ${SLURM_NTASKS:-unknown}"
        echo "SLURM_CPUS_PER_TASK = ${SLURM_CPUS_PER_TASK:-unknown}"
        echo "Using NPROCS        = ${NPROCS}"

        srun -n "${NPROCS}" "$EXEC" \
            -options_file "$SOLVER_OPTS" \
            -options_file "$GEOM_OPTS"   \
            -options_file "$EXP_OPTS"    \
            -output_path  "$folder" \
            "${EXTRA_OPTS[@]}" \
            | tee "$folder/outp.txt"
    else
        echo "No SLURM environment detected; running locally with mpiexec."
        echo "Using NPROCS = ${NPROCS}"
        echo "MPIEXEC      = ${MPIEXEC}"

        "$MPIEXEC" -n "${NPROCS}" "$EXEC" \
            -options_file "$SOLVER_OPTS" \
            -options_file "$GEOM_OPTS"   \
            -options_file "$EXP_OPTS"    \
            -output_path  "$folder" \
            "${EXTRA_OPTS[@]}" \
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
