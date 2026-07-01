#!/usr/bin/env bash
# =============================================================================
# run_mob_sweep.sh — Sequential mob_sub sweep for 2D_two_ice_grains_boundary
#
# Runs one simulation per mob_sub value, each into its own timestamped folder
# under $RESULTS_BASE/2D_two_ice_grains_boundary/.
#
# Usage:
#   ./scripts/Studio/run_mob_sweep.sh [experiment] [extra petsc opts...]
#
#   experiment   Experiment opts name (default: 30day_T-5_h1.00_GTphys)
#
# Example:
#   ./scripts/Studio/run_mob_sweep.sh 30day_T-5_h1.00_GTphys
#   ./scripts/Studio/run_mob_sweep.sh 5yr_T-5_h1.00_GTphys
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PETIGA_DIR}/projects/permafrost"
EXEC="$PROJECT_ROOT/permafrost"
INPUTS_DIR="$PROJECT_ROOT/inputs"
RESULTS_BASE="/Users/jacksonbaglino/SimulationResults/permafrost/scratch"
GEOM="2D_two_ice_grains_boundary"
EXP="${1:-30day_T-5_h1.00_GTphys}"

SOLVER_OPTS="$INPUTS_DIR/solver.opts"
GEOM_OPTS="$INPUTS_DIR/geometry/${GEOM}.opts"
EXP_OPTS="$INPUTS_DIR/experiment/${EXP}.opts"

# mob_sub values to sweep (log-spaced: baseline → 100× reduction)
MOB_VALUES=(
    6.5404e-10   # K&P baseline (safety=0.5); reference
    1.3e-10      #  5× reduction
    6.5e-11      # 10× reduction
    1.3e-11      # 50× reduction
    6.5e-12      # 100× reduction
)

# Compute NPROCS from geometry opts
Nx=$(awk '$1=="-Nx"{print $2}' "$GEOM_OPTS" | tail -n1)
Ny=$(awk '$1=="-Ny"{print $2}' "$GEOM_OPTS" | tail -n1)
Nz=$(awk -v default=1 '$1=="-Nz"{v=$2} END{print (v=="")?default:v}' "$GEOM_OPTS")
dof=3
total=$(( Nx * Ny * Nz * dof ))
TARGET=10000
hw=$(sysctl -n hw.logicalcpu 2>/dev/null || nproc 2>/dev/null || echo 12)
cap=$(( hw < 12 ? hw : 12 ))
NPROCS=$(( (total + TARGET - 1) / TARGET ))
(( NPROCS < 1 )) && NPROCS=1
(( NPROCS > cap )) && NPROCS=$cap

echo "========================================================================"
echo "  mob_sub sweep — $GEOM × $EXP"
echo "  NPROCS = $NPROCS  (${total} total DOFs)"
echo "  Values: ${MOB_VALUES[*]}"
echo "========================================================================"

# Compile once before the sweep
echo ""
echo "--- Compiling ---"
make -C "$PROJECT_ROOT" all

SUMMARY="$RESULTS_BASE/$GEOM/mob_sweep_summary.txt"
mkdir -p "$RESULTS_BASE/$GEOM"
printf "%-14s  %-50s  %s\n" "mob_sub" "output folder" "status" > "$SUMMARY"
printf "%-14s  %-50s  %s\n" "--------------" "--------------------------------------------------" "------" >> "$SUMMARY"

for mob in "${MOB_VALUES[@]}"; do
    ts=$(date +%Y-%m-%d__%H.%M.%S)
    tag="mob_${mob}"
    folder="$RESULTS_BASE/$GEOM/${ts}_${EXP}_${tag}"
    mkdir -p "$folder"
    cp "$SOLVER_OPTS" "$folder/"
    cp "$GEOM_OPTS"   "$folder/"
    cp "$EXP_OPTS"    "$folder/"

    echo ""
    echo "=========================================================================="
    echo "  mob_sub = $mob"
    echo "  → $folder"
    echo "=========================================================================="

    export folder
    set +e
    mpiexec -np "$NPROCS" "$EXEC"       \
        -options_file "$SOLVER_OPTS"    \
        -options_file "$GEOM_OPTS"      \
        -options_file "$EXP_OPTS"       \
        -output_path  "$folder"         \
        -mob_sub      "$mob"            \
        2>&1 | tee "$folder/outp.txt"
    exit_code=${PIPESTATUS[0]}
    set -e

    status="OK"
    (( exit_code != 0 )) && status="FAIL(${exit_code})"
    printf "%-14s  %-50s  %s\n" "$mob" "$(basename "$folder")" "$status" >> "$SUMMARY"
    echo "  → $status"
done

echo ""
echo "========================================================================"
echo "  Sweep complete. Summary:"
cat "$SUMMARY"
echo "========================================================================"
