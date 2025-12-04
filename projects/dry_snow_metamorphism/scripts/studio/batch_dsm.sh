#!/bin/zsh
###############################################################################
# Script: batch_dsm.sh
# Model: Dry Snow Metamorphism (DSM)
# Purpose:
#   Iterate over *grain directories* under ./inputs and call run_dsm.sh for
#   each directory (and for each temp/RH/gradient combo defined below).
#   A grain directory must contain:
#     - grains.dat
#     - .env            (Lx, Ly, Lz; Nx/Ny/Nz/eps may be added automatically)
#     - metadata.json
###############################################################################

set -o errexit
set -o nounset
set -o pipefail

# Enable zsh glob qualifiers like (N) and directory slash matches
setopt extended_glob null_glob

# Project roots
: ${PETIGA_DIR:="/Users/jacksonbaglino/PetIGA-3.20"}
BASE_DIR="${PETIGA_DIR}/projects/dry_snow_metamorphism"
INPUT_ROOT="$BASE_DIR/inputs"
RUN_SCRIPT="$BASE_DIR/scripts/studio/run_dsm.sh"

# ==========================
# Sweeps (edit as needed)
# ==========================
# Temperatures (°C)
temps=( -20 -5 )
humidities=( 0.98 )
grad_xs=( 0.0 )
grad_ys=( 3.0e-6 )
grad_zs=( 0.0 )

# MPI/output overrides (optional) — match run_dsm.sh
: ${NUM_PROCS:=12}
: ${output_dir:="/Users/jacksonbaglino/SimulationResults/dry_snow_metamorphism/scratch"}

# ==========================
# Discover grain directories
# ==========================
# Default pattern picks seed-stable names like grains__phi=...__Lxmm=...__seed=...
pattern="$INPUT_ROOT/grains__phi=0.24__Lxmm=2__Lymm=2__seed=21"
# (N) => null glob (no error if empty), (/) => only directories
grain_dirs=( $~pattern(N/) )

if (( ${#grain_dirs[@]} == 0 )); then
  echo "[ERROR] No grain directories found under: $pattern"
  echo "        Generate inputs first (e.g., preprocess/batch_generate_grains.sh)."
  exit 1
fi

echo "Starting DSM batch: ${#temps[@]} temps × ${#humidities[@]} humidities × ${#grad_ys[@]} gradY × ${#grain_dirs[@]} grain-dirs"
echo "================================================================================"

for GRAIN_DIR in ${grain_dirs[@]}; do
  # Basic sanity of required files
  if [[ ! -f "$GRAIN_DIR/grains.dat" || ! -f "$GRAIN_DIR/grains.env" || ! -f "$GRAIN_DIR/metadata.json" ]]; then
    echo "[WARN] Skipping (missing required files): $GRAIN_DIR"
    continue
  fi

  for t in ${temps[@]}; do
    for rh in ${humidities[@]}; do
      for gx in ${grad_xs[@]}; do
        for gy in ${grad_ys[@]}; do
          for gz in ${grad_zs[@]}; do
            # Note: in run_dsm.sh, `-h` flag is for humidity; help is `-?` or `--help`.
            echo "[run] T=$t°C RH=$rh G=($gx,$gy,$gz) dir=$(basename "$GRAIN_DIR")"
            # Build filename relative to input_dir for run_dsm.sh
            rel_path="${GRAIN_DIR#${INPUT_ROOT}/}"
            filename_rel="$rel_path/grains.dat"

            # Export variables expected by run_dsm.sh and pass compatible flags
            NUM_PROCS="$NUM_PROCS" \
            output_dir="$output_dir" \
            input_dir="$INPUT_ROOT" \
            grad_temp0X="$gx" grad_temp0Y="$gy" grad_temp0Z="$gz" \
            "$RUN_SCRIPT" \
              -f "$filename_rel" \
              -t "$t" \
              -h "$rh" \
              -p "$NUM_PROCS"
            echo "[done] T=$t°C RH=$rh G=($gx,$gy,$gz) dir=$(basename "$GRAIN_DIR")"
            echo "--------------------------------------------------------------------------------"
          done
        done
      done
    done
  done

done

echo "[all done] DSM batch complete."
