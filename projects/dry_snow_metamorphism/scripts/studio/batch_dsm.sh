#!/bin/zsh
###############################################################################
# Script: batch_dsm.sh
# Model: Dry Snow Metamorphism (DSM)
# Purpose:
#   Iterate over all input .dat files in ./inputs/porespy/dat and call
#   run_dsm.sh for each file (and for each temp/RH combo defined below).
###############################################################################

set -o errexit
set -o nounset
set -o pipefail

# Enable zsh glob qualifiers like (N) and avoid errors on empty globs
setopt extended_glob null_glob

# Project roots
: ${PETIGA_DIR:="/Users/jacksonbaglino/PetIGA-3.20"}
BASE_DIR="${PETIGA_DIR}/projects/dry_snow_metamorphism"
INPUT_ROOT="$BASE_DIR/inputs"
RUN_SCRIPT="$BASE_DIR/scripts/studio/run_dsm.sh"

# Sweeps (edit as needed)
temps=( -30 )
humidities=( 0.98 )

# Collect all files in ./inputs/porespy/dat (use glob first; fallback to find)
input_glob="$INPUT_ROOT/porespy/dat/*.dat"
input_files=( $input_glob(N) )

if (( ${#input_files[@]} == 0 )); then
  echo "[WARN] Glob found 0 files in: $input_glob"
  if [[ -d "$INPUT_ROOT/porespy/dat" ]]; then
    echo "[INFO] Directory exists. Listing to help debug:"; ls -l "$INPUT_ROOT/porespy/dat" || true
  else
    echo "[ERROR] Directory does not exist: $INPUT_ROOT/porespy/dat"
  fi
  # Fallback: find
  input_files=( ${(f)"$(command -v find >/dev/null 2>&1 && find "$INPUT_ROOT/porespy/dat" -type f -name '*.dat' -print 2>/dev/null)"} )
fi

if (( ${#input_files[@]} == 0 )); then
  echo "[ERROR] No .dat files found under: $INPUT_ROOT/porespy/dat/"
  exit 1
fi

echo "Starting DSM batch: ${#temps[@]} temps × ${#humidities[@]} humidities × ${#input_files[@]} files"
echo "=================================================="

for fpath in ${input_files[@]}; do
  fname=$(basename -- "$fpath")
  for t in ${temps[@]}; do
    for rh in ${humidities[@]}; do
      echo "Starting simulation with T=$t, RH=$rh, file=$fname"
      export temp=$t
      export humidity=$rh
      export filename="$fname"                # basename only
      export FILE_SUBDIR="porespy/dat"        # tell run_dsm.sh which subdir to use
      zsh "$RUN_SCRIPT"
      echo "Finished simulation with T=$t, RH=$rh, file=$fname"
      echo "=================================================="
    done
  done

done

exit 0
