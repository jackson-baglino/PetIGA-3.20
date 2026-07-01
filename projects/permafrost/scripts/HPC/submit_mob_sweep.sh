#!/usr/bin/env bash
#SBATCH -J pf_mob_sweep
#SBATCH -A rubyfu
#SBATCH -t 0-12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=11
#SBATCH --cpus-per-task=1
#SBATCH --partition=expansion
#SBATCH --mem-per-cpu=2G
#SBATCH -o "output_files/%x.o%j"
#SBATCH -e "output_files/%x.e%j"
#SBATCH --mail-user=jbaglino@caltech.edu
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --constraint='icelake|skylake|cascadelake'
# =============================================================================
# submit_mob_sweep.sh — Sequential mob_sub sweep, single SLURM job
#
# Runs 5 mob_sub values for 2D_two_ice_grains_boundary, one after another,
# all writing into a single shared parent folder:
#   $SCRATCH/permafrost/2D_two_ice_grains_boundary/mob_sweep_<timestamp>/
#
# Usage (from project root):
#   sbatch scripts/HPC/submit_mob_sweep.sh [experiment]
#
#   experiment   Experiment opts name (default: 30day_T-5_h1.00_GTphys)
#
# Override wall time or nodes via sbatch flags before the script name:
#   sbatch --time=0-24:00:00 scripts/HPC/submit_mob_sweep.sh 5yr_T-5_h1.00_GTphys
# =============================================================================
set -euo pipefail

GEOM="2D_two_ice_grains_boundary"
EXP="${1:-30day_T-5_h1.00_GTphys}"

# mob_sub values to sweep (baseline → 100× reduction)
MOB_VALUES=(
    6.5404e-10   # K&P baseline (safety=0.5); reference
    1.3e-10      #  5× reduction
    6.5e-11      # 10× reduction
    1.3e-11      # 50× reduction
    6.5e-12      # 100× reduction
)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT="${SLURM_SUBMIT_DIR}"
EXEC="$PROJECT_ROOT/permafrost"
INPUTS_DIR="$PROJECT_ROOT/inputs"
SOLVER_OPTS="$INPUTS_DIR/solver.opts"
GEOM_OPTS="$INPUTS_DIR/geometry/${GEOM}.opts"
EXP_OPTS="$INPUTS_DIR/experiment/${EXP}.opts"

BATCH_DIR="$SCRATCH/permafrost/${GEOM}/mob_sweep_$(date +%Y-%m-%d__%H.%M.%S)"
mkdir -p "$BATCH_DIR"

SUMMARY="$BATCH_DIR/summary.txt"
printf "%-14s  %-45s  %-8s  %s\n" "mob_sub" "folder" "status" "wall_s" > "$SUMMARY"
printf "%-14s  %-45s  %-8s  %s\n" \
    "--------------" "---------------------------------------------" "--------" "------" >> "$SUMMARY"

# ---------------------------------------------------------------------------
# Compile once on the compute node
# ---------------------------------------------------------------------------
echo "--- Compiling ---"
cd "$PROJECT_ROOT"
make clean && make all
echo "--- Compilation done ---"
echo ""

NPROCS="${SLURM_NTASKS:-11}"

echo "========================================================================"
echo "  mob_sub sweep — $GEOM × $EXP"
echo "  NPROCS     = $NPROCS"
echo "  Output dir = $BATCH_DIR"
echo "  Values     = ${MOB_VALUES[*]}"
echo "========================================================================"

# ---------------------------------------------------------------------------
# Run each mob_sub value sequentially
# ---------------------------------------------------------------------------
for mob in "${MOB_VALUES[@]}"; do
    run_dir="$BATCH_DIR/mob_${mob}"
    mkdir -p "$run_dir"
    cp "$SOLVER_OPTS" "$run_dir/"
    cp "$GEOM_OPTS"   "$run_dir/"
    cp "$EXP_OPTS"    "$run_dir/"

    echo ""
    echo "=========================================="
    echo "  mob_sub = $mob"
    echo "  → $run_dir"
    echo "=========================================="

    export folder="$run_dir"
    t0=$(date +%s)

    set +e
    srun -n "$NPROCS" "$EXEC"           \
        -options_file "$SOLVER_OPTS"    \
        -options_file "$GEOM_OPTS"      \
        -options_file "$EXP_OPTS"       \
        -output_path  "$run_dir"        \
        -mob_sub      "$mob"            \
        2>&1 | tee "$run_dir/outp.txt"
    exit_code=${PIPESTATUS[0]}
    set -e

    t1=$(date +%s)
    elapsed=$(( t1 - t0 ))
    status="OK"
    (( exit_code != 0 )) && status="FAIL(${exit_code})"

    printf "%-14s  %-45s  %-8s  %s\n" \
        "$mob" "mob_${mob}" "$status" "${elapsed}s" >> "$SUMMARY"

    echo "  → $status  (${elapsed}s)"
done

echo ""
echo "========================================================================"
echo "  Sweep complete. Summary:"
cat "$SUMMARY"
echo "  All results in: $BATCH_DIR"
echo "========================================================================"
