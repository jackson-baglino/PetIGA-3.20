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
temps=( -25 )
# Relative humidities (0..1)
humidities=( 0.98 )
# Temperature gradient components (K/m)
grad_xs=( 0.0 )
grad_ys=( 3.0e-6 )
grad_zs=( 0.0 )

# MPI/launch overrides (optional)
: ${NTASKS:=8}
: ${LAUNCHER:=mpirun}
: ${OUT_ROOT:="$BASE_DIR/outputs"}

# ==========================
# Discover grain directories
# ==========================
# Default pattern picks seed-stable names like grains__phi=...__Lxmm=...__seed=...
pattern="$INPUT_ROOT/grains__*"
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
            echo "[run] T=$t°C RH=$rh G=($gx,$gy,$gz) dir=$(basename "$GRAIN_DIR")"
            # Pass parameters as flags to the new grain-dir runner
            NTASKS="$NTASKS" LAUNCHER="$LAUNCHER" OUT_ROOT="$OUT_ROOT" \
            "$RUN_SCRIPT" \
              --grain-dir "$GRAIN_DIR" \
              --temp-c "$t" \
              --humidity "$rh" \
              --grad-x "$gx" \
              --grad-y "$gy" \
              --grad-z "$gz"
            echo "[done] T=$t°C RH=$rh G=($gx,$gy,$gz) dir=$(basename "$GRAIN_DIR")"
            echo "--------------------------------------------------------------------------------"
          done
        done
      done
    done
  done

done

echo "[all done] DSM batch complete."
