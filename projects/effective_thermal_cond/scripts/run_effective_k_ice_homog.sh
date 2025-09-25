#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# run_effective_k_ice_homog.sh
#
# Goal: Keep this script as simple as possible. The batch script discovers
#       directories and exports env vars (Nx, Ny, Nz, Lx, Ly, Lz, eps, INIT_DIR).
#       This script:
#         0) validates inputs,
#         1) compiles (minimal),
#         2) runs the simulation,
#         3) collects outputs,
#         4) runs post-processing (optional).
#
# Assumptions:
#   - Called as: ./scripts/run_effective_k_ice_homog.sh <init_dir>
#   - The caller (batch script) has already exported: Nx Ny Nz Lx Ly Lz eps
#   - Optional env with defaults here: dim (2), NUM_PROCS (1),
#     TEMP_TOP (273.15-30), FLUX_BOTTOM (-0.1)
#   - Global dry run: set DRY_RUN=1 to print Steps 2–5 and exit without changes.
# -----------------------------------------------------------------------------

# Record start time for total runtime reporting
__t0__=$(date +%s)

# Bash safety + predictable word splitting/globbing
set -euo pipefail
IFS=$'\n\t'
shopt -s nullglob

# --- Step 0: Resolve arguments and minimal defaults --------------------------
if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <INIT_DIR>"
  exit 1
fi

INIT_DIR="$1"   # directory with initial condition files / metadata
if [[ ! -d "$INIT_DIR" ]]; then
  echo "ERROR: INIT_DIR is not a directory: $INIT_DIR"
  exit 1
fi

# Optional knobs (kept here, not in the batch script)
: "${dim:=2}"                  # 2D by default
: "${NUM_PROCS:=1}"            # single-rank default
: "${TEMP_TOP:=$(awk 'BEGIN{print 273.15-30}')}"; # top BC (K)
: "${FLUX_BOTTOM:=-0.1}"  # bottom BC (W/m^2)
: "${OUTPUT_BINARY:=1}"      # whether to output binary files (0/1)
: "${SOL_INDEX:=0}"        # which solid phase to use (0-based index, default 0)

# Resolve the output directory early so the solver can read it via getenv()
OUT_ROOT="${OUT_ROOT:-/Users/jacksonbaglino/SimulationResults/effective_thermal_cond/scratch}"
base="$(basename "$INIT_DIR")"
OUTPUT_DIR="$OUT_ROOT/${base}"
export OUT_ROOT OUTPUT_DIR

# --- Step 1: Validate required env set by the batch script -------------------
missing=()
for v in Nx Ny Nz Lx Ly Lz eps; do
  if [[ -z "${!v:-}" ]]; then missing+=("$v"); fi
done
if (( ${#missing[@]} > 0 )); then
  echo "ERROR: missing required environment variable(s): ${missing[*]}"
  echo "Hint: ensure the batch script sourced grains.env and exported these."
  exit 1
fi

# Compact summary (safe under -u because we just validated)
echo "=== Effective k (homog) – run configuration ==="
echo "INIT_DIR     : $INIT_DIR"
if (( dim == 2 )); then
  echo "Grid         : Nx=$Nx, Ny=$Ny (dim=$dim)"
  echo "Domain (m)   : Lx=$Lx, Ly=$Ly, eps=$eps"
else
  echo "Grid         : Nx=$Nx, Ny=$Ny, Nz=$Nz (dim=$dim)"
  echo "Domain (m)   : Lx=$Lx, Ly=$Ly, Lz=$Lz, eps=$eps"
fi
echo "BCs / Params : TEMP_TOP=$TEMP_TOP, FLUX_BOTTOM=$FLUX_BOTTOM"
echo "MPI          : NUM_PROCS=$NUM_PROCS"

# -----------------------------------------------------------------------------
# Global dry-run: print intended actions for Steps 2–5 and exit.
#   - Set DRY_RUN=1 in environment to enable.
# -----------------------------------------------------------------------------
if [[ "${DRY_RUN:-0}" == "1" ]]; then
  echo "[DRY-RUN] Build : make BUILD=release effective_k_ice_homog"

  if (( NUM_PROCS > 1 )); then
    echo "[DRY-RUN] Run   : mpiexec -np \"$NUM_PROCS\" ./effective_k_ice_homog -init_dir \"$INIT_DIR\""
  else
    echo "[DRY-RUN] Run   : ./effective_k_ice_homog -init_dir \"$INIT_DIR\""
  fi

  echo "[DRY-RUN] Collect -> $OUTPUT_DIR (move: *.dat *.bin *.info *.csv)"

  POST_SCRIPT="${POST_SCRIPT:-postprocess/plot_vector_field.py}"
  echo "[DRY-RUN] Post  : python3 \"$POST_SCRIPT\" \"$OUTPUT_DIR\" \"$Nx\" \"$Ny\" \"$Lx\" \"$Ly\""

  exit 0
fi

# -----------------------------------------------------------------------------
# Step 2: Compile (minimal)
#   - Uses Makefile target `effective_k_ice_homog` with BUILD=release.
#   - `make` is incremental; it will skip if everything is up-to-date.
#   - Set NO_BUILD=1 to skip compilation (useful for quick re-runs).
# -----------------------------------------------------------------------------
if [[ "${NO_BUILD:-0}" != "1" ]]; then
  echo "[Build] make BUILD=release effective_k_ice_homog"
  make BUILD=release effective_k_ice_homog
else
  echo "[Build] Skipping compilation because NO_BUILD=1"
fi

# -----------------------------------------------------------------------------
# Step 3: Run the simulation (minimal)
#   - Use mpiexec only if NUM_PROCS > 1.
#   - Pass only the minimal flag: the INIT_DIR we were given.
#   - Set NO_RUN=1 to skip execution (useful for testing builds).
# -----------------------------------------------------------------------------
# Finish exporting any remaining env vars needed.
export INIT_DIR dim TEMP_TOP FLUX_BOTTOM NUM_PROCS OUTPUT_BINARY SOL_INDEX OUTPUT_DIR

# Ensure output directory exists (solver won't create it)
# Ensure the output root exists (OUTPUT_DIR was computed earlier)
mkdir -p "${OUT_ROOT}"

mkdir -p "$OUTPUT_DIR"

if [[ "${NO_RUN:-0}" == "1" ]]; then
  echo "[Run] Skipping execution because NO_RUN=1"
else
  # Ensure the binary exists before trying to run it.
  if [[ ! -x ./effective_k_ice_homog ]]; then
    echo "ERROR: ./effective_k_ice_homog not found (or not executable). Did the build succeed?"
    exit 1
  fi

  echo "[Run] Launching simulation (NUM_PROCS=$NUM_PROCS)"
  if (( NUM_PROCS > 1 )); then
    mpiexec -np "$NUM_PROCS" ./effective_k_ice_homog -init_dir "$INIT_DIR"
  else
    ./effective_k_ice_homog -init_dir "$INIT_DIR"
  fi
fi

# -----------------------------------------------------------------------------
# Step 4: Collect outputs to a simple OUTPUT_DIR
#   - Predictable root (override with OUT_ROOT env if you like)
#   - Folder name uses only the init-dir basename
#   - Moves: *.dat *.bin *.info *.csv (if present)
# -----------------------------------------------------------------------------
echo "[Collect] Moving outputs to: $OUTPUT_DIR"

# Collect matches and move only if non-empty.
# - Arrays avoid word-splitting problems; `--` protects against leading '-' in filenames.
dat=( *.dat );  (( ${#dat[@]}  ))  && mv -- "${dat[@]}"  "$OUTPUT_DIR"
bin=( *.bin );  (( ${#bin[@]}  ))  && mv -- "${bin[@]}"  "$OUTPUT_DIR"
inf=( *.info ); (( ${#inf[@]}  ))  && mv -- "${inf[@]}"  "$OUTPUT_DIR"
csv=( *.csv );  (( ${#csv[@]}  ))  && mv -- "${csv[@]}"  "$OUTPUT_DIR"

echo "[Collect] Done. Output dir: $OUTPUT_DIR"
echo "[Collect] Contents:"
ls -lh "$OUTPUT_DIR" || true

# -----------------------------------------------------------------------------
# Step 5: Post-process (optional)
#   - Default script: postprocess/plot_vector_field.py
#   - Args: OUTPUT_DIR, Nx, Ny, Lx, Ly (adjust to your script’s signature)
#   - Set NO_POST=1 to skip. Set POST_SCRIPT to override the script path.
# -----------------------------------------------------------------------------
if [[ "${NO_POST:-0}" == "1" ]]; then
  echo "[Post] Skipping post-processing because NO_POST=1"
else
  POST_SCRIPT="${POST_SCRIPT:-postprocess/plot_vector_field.py}"
  if [[ -f "$POST_SCRIPT" ]]; then
    echo "[Post] python3 \"$POST_SCRIPT\" \"$OUTPUT_DIR\" \"$Nx\" \"$Ny\" \"$Lx\" \"$Ly\""
    python3 "$POST_SCRIPT" "$OUTPUT_DIR" "$Nx" "$Ny" "$Lx" "$Ly"
  else
    echo "[Post] WARNING: Post-process script not found: $POST_SCRIPT (skipping)"
  fi
fi

# Copy SSA_evo.dat from input directory to output directory if it exists
if [[ -f "$INIT_DIR/SSA_evo.dat" ]]; then
  cp "$INIT_DIR/SSA_evo.dat" "$OUTPUT_DIR/"
  echo "[Collect] Copied SSA_evo.dat to output directory."
else
  echo "[Collect] SSA_evo.dat not found in INIT_DIR; skipping copy."
fi

# Copy metadata.json from input directory to output directory if it exists
if [[ -f "$INIT_DIR/metadata.json" ]]; then
  cp "$INIT_DIR/metadata.json" "$OUTPUT_DIR/"
  echo "[Collect] Copied metadata.json to output directory."
else
  echo "[Collect] metadata.json not found in INIT_DIR; skipping copy."
fi

# Total runtime reporting
__t1__=$(date +%s)
echo "[Done] Completed in $(( __t1__ - __t0__ )) s"