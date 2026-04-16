#!/bin/bash
#SBATCH --job-name=keff_homog
#SBATCH --partition=expansion
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=12:00:00
#SBATCH --output=logs/keff_%j.out
#SBATCH --error=logs/keff_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jbaglino@caltech.edu

# =============================================================================
# run_hpc.slurm — Caltech Resnick HPC SLURM script
# Computes effective thermal conductivity via periodic homogenization.
#
# Usage:
#   sbatch scripts/run_hpc.slurm
#
# Required environment (set in your ~/.bashrc or via sbatch --export):
#   PETSC_DIR   — path to PETSc installation
#   PETSC_ARCH  — PETSc arch string
#   PETIGA_DIR  — path to PetIGA installation
#   INIT_DIR    — directory containing sol_*.dat and igasol.dat (or grains.opts)
#   OUTPUT_DIR  — where to write k_eff.csv and t_vec.dat
#
# Example:
#   INIT_DIR=/path/to/sim OUTPUT_DIR=/path/to/out sbatch scripts/run_hpc.slurm
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# 1. Load modules (adjust module names to match Resnick's available modules)
# ---------------------------------------------------------------------------
module purge
module load gcc/12
module load openmpi/4
# PETSc and PetIGA are typically installed manually on HPC clusters.
# Set PETSC_DIR/PETSC_ARCH/PETIGA_DIR in your environment or uncomment below:
# export PETSC_DIR=/central/scratch/$USER/petsc
# export PETSC_ARCH=arch-linux-c-opt
# export PETIGA_DIR=/central/scratch/$USER/petiga

# Validate required environment variables
for v in PETSC_DIR PETSC_ARCH PETIGA_DIR; do
  if [[ -z "${!v:-}" ]]; then
    echo "ERROR: $v is not set. Export it before submitting with sbatch."
    exit 1
  fi
done

# ---------------------------------------------------------------------------
# 2. I/O directories
# ---------------------------------------------------------------------------
INIT_DIR="${INIT_DIR:-}"
OUTPUT_DIR="${OUTPUT_DIR:-}"

if [[ -z "$INIT_DIR" || ! -d "$INIT_DIR" ]]; then
  echo "ERROR: INIT_DIR must be set to a valid directory."
  echo "Usage: INIT_DIR=/path/to/sim OUTPUT_DIR=/path/to/out sbatch scripts/run_hpc.slurm"
  exit 1
fi

if [[ -z "$OUTPUT_DIR" ]]; then
  OUTPUT_DIR="$SLURM_SUBMIT_DIR/outputs/$(basename "$INIT_DIR")"
fi
mkdir -p "$OUTPUT_DIR"
mkdir -p logs

# ---------------------------------------------------------------------------
# 3. Resolve options file
# ---------------------------------------------------------------------------
if [[ -f "$INIT_DIR/grains.opts" ]]; then
  OPTS_FILE="$INIT_DIR/grains.opts"
elif [[ -f "$INIT_DIR/grains.env" ]]; then
  echo "Generating grains.opts from grains.env ..."
  python3 "$SLURM_SUBMIT_DIR/scripts/gen_opts.py" "$INIT_DIR/grains.env" \
    > "$INIT_DIR/grains.opts"
  OPTS_FILE="$INIT_DIR/grains.opts"
else
  OPTS_FILE="$SLURM_SUBMIT_DIR/inputs/default.opts"
  echo "WARNING: No grains.opts or grains.env found; using $OPTS_FILE"
fi

echo "============================"
echo " Effective k HPC run"
echo " INIT_DIR   : $INIT_DIR"
echo " OUTPUT_DIR : $OUTPUT_DIR"
echo " OPTS_FILE  : $OPTS_FILE"
echo " SLURM_JOB_ID : ${SLURM_JOB_ID:-N/A}"
echo " Nodes/Ranks  : $SLURM_NNODES / $SLURM_NTASKS"
echo "============================"

# ---------------------------------------------------------------------------
# 4. Build
# ---------------------------------------------------------------------------
cd "$SLURM_SUBMIT_DIR"
make effective_k_ice_homog

# ---------------------------------------------------------------------------
# 5. Run
# ---------------------------------------------------------------------------
srun --mpi=pmix ./effective_k_ice_homog \
  -options_file "$OPTS_FILE" \
  -init_dir     "$INIT_DIR" \
  -init_mode    "file" \
  -output_dir   "$OUTPUT_DIR"

echo "Simulation completed. Results in $OUTPUT_DIR"
