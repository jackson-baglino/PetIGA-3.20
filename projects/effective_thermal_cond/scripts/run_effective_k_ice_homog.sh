#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# run_effective_k_ice_homog.sh
#
# Goal: Keep this script as simple as possible. This script:
#         0) validates inputs,
#         1) resolves a PETSc options file (grains.opts preferred),
#         2) compiles (minimal),
#         3) runs the simulation,
#         4) collects outputs,
#         5) runs post-processing (optional).
#
# Assumptions:
#   - Called as: ./scripts/run_effective_k_ice_homog.sh <init_dir>
#   - Preferred: INIT_DIR contains grains.opts (see inputs/grains_template.opts
#     and scripts/gen_opts.py for how to create one from grains.env)
#   - Fallback: fall back to inputs/default.opts and use env vars on CLI
#   - Optional env with defaults here: NUM_PROCS (1), INIT_MODE (file)
#   - Global dry run: set DRY_RUN=1 to print Steps 2-5 and exit without changes.
# -----------------------------------------------------------------------------

# Record start time for total runtime reporting
__t0__=$(date +%s)

# Bash safety + predictable word splitting/globbing
set -euo pipefail
IFS=$'\n\t'
shopt -s nullglob

# --- Step 0: Resolve arguments -----------------------------------------------
if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <INIT_DIR>"
  exit 1
fi

INIT_DIR="$1"
if [[ ! -d "$INIT_DIR" ]]; then
  echo "ERROR: INIT_DIR is not a directory: $INIT_DIR"
  exit 1
fi

# Optional knobs
: "${NUM_PROCS:=1}"      # MPI ranks (default: serial)
: "${INIT_MODE:=file}"   # ice field mode: file | circle | layered

# Resolve output directory
OUT_ROOT="${OUT_ROOT:-/Users/jacksonbaglino/SimulationResults/effective_thermal_cond/scratch}"
base="$(basename "$INIT_DIR")"
OUTPUT_DIR="$OUT_ROOT/${base}"
export OUT_ROOT OUTPUT_DIR

# --- Step 1: Resolve options source (grains.opts preferred over env vars) ----

if [[ -f "$INIT_DIR/grains.opts" ]]; then
  OPTS_FILE="$INIT_DIR/grains.opts"
  echo "[Options] Using: $OPTS_FILE"
elif [[ -f "$INIT_DIR/grains.env" ]]; then
  echo "[Options] No grains.opts found; generating from grains.env ..."
  python3 scripts/gen_opts.py "$INIT_DIR/grains.env" > "$INIT_DIR/grains.opts"
  OPTS_FILE="$INIT_DIR/grains.opts"
  echo "[Options] Generated: $OPTS_FILE"
else
  OPTS_FILE="inputs/default.opts"
  echo "[Options] No grains.opts or grains.env in INIT_DIR; using $OPTS_FILE"
fi

echo "=== Effective k (homog) – run configuration ==="
echo "INIT_DIR     : $INIT_DIR"
echo "OUTPUT_DIR   : $OUTPUT_DIR"
echo "OPTS_FILE    : $OPTS_FILE"
echo "INIT_MODE    : $INIT_MODE"
echo "MPI          : NUM_PROCS=$NUM_PROCS"

# -----------------------------------------------------------------------------
# Global dry-run: print intended actions for Steps 2-5 and exit.
#   - Set DRY_RUN=1 in environment to enable.
# -----------------------------------------------------------------------------
if [[ "${DRY_RUN:-0}" == "1" ]]; then
  echo "[DRY-RUN] Build : make effective_k_ice_homog"

  _run_args="-options_file \"$OPTS_FILE\" -init_dir \"$INIT_DIR\" -init_mode \"$INIT_MODE\" -output_dir \"$OUTPUT_DIR\""
  if (( NUM_PROCS > 1 )); then
    echo "[DRY-RUN] Run   : mpiexec -np $NUM_PROCS ./effective_k_ice_homog $_run_args"
  else
    echo "[DRY-RUN] Run   : ./effective_k_ice_homog $_run_args"
  fi

  echo "[DRY-RUN] Collect -> $OUTPUT_DIR (move: *.dat *.bin *.info *.csv)"

  POST_SCRIPT="${POST_SCRIPT:-postprocess/plot_k_eff.py}"
  echo "[DRY-RUN] Post  : python3 \"$POST_SCRIPT\" \"$OUTPUT_DIR\""

  exit 0
fi

# -----------------------------------------------------------------------------
# Step 2: Compile (minimal)
#   - `make` is incremental; it will skip if everything is up-to-date.
#   - Set NO_BUILD=1 to skip compilation (useful for quick re-runs).
# -----------------------------------------------------------------------------
if [[ "${NO_BUILD:-0}" != "1" ]]; then
  echo "[Build] make effective_k_ice_homog"
  make effective_k_ice_homog
else
  echo "[Build] Skipping compilation because NO_BUILD=1"
fi

# -----------------------------------------------------------------------------
# Step 3: Run the simulation
#   - Pass -options_file so all mesh/physics params come from the .opts file.
#   - Also pass -init_dir, -init_mode, and -output_dir on the command line
#     so they can be overridden without editing the .opts file.
#   - Set NO_RUN=1 to skip execution.
# -----------------------------------------------------------------------------
mkdir -p "${OUT_ROOT}"
mkdir -p "$OUTPUT_DIR"

if [[ "${NO_RUN:-0}" == "1" ]]; then
  echo "[Run] Skipping execution because NO_RUN=1"
else
  if [[ ! -x ./effective_k_ice_homog ]]; then
    echo "ERROR: ./effective_k_ice_homog not found or not executable. Did the build succeed?"
    exit 1
  fi

  echo "[Run] Launching simulation (NUM_PROCS=$NUM_PROCS)"
  if (( NUM_PROCS > 1 )); then
    mpiexec -np "$NUM_PROCS" ./effective_k_ice_homog \
      -options_file "$OPTS_FILE" \
      -init_dir     "$INIT_DIR" \
      -init_mode    "$INIT_MODE" \
      -output_dir   "$OUTPUT_DIR"
  else
    ./effective_k_ice_homog \
      -options_file "$OPTS_FILE" \
      -init_dir     "$INIT_DIR" \
      -init_mode    "$INIT_MODE" \
      -output_dir   "$OUTPUT_DIR"
  fi
fi

# -----------------------------------------------------------------------------
# Step 4: Collect outputs
#   - Moves: *.dat *.bin *.info *.csv (if present in CWD)
# -----------------------------------------------------------------------------
echo "[Collect] Moving outputs to: $OUTPUT_DIR"

dat=( *.dat );  (( ${#dat[@]}  )) && mv -- "${dat[@]}"  "$OUTPUT_DIR"
bin=( *.bin );  (( ${#bin[@]}  )) && mv -- "${bin[@]}"  "$OUTPUT_DIR"
inf=( *.info ); (( ${#inf[@]}  )) && mv -- "${inf[@]}"  "$OUTPUT_DIR"
csv=( *.csv );  (( ${#csv[@]}  )) && mv -- "${csv[@]}"  "$OUTPUT_DIR"

echo "[Collect] Done. Contents of $OUTPUT_DIR:"
ls -lh "$OUTPUT_DIR" || true

# Copy ancillary files from INIT_DIR if present
for f in SSA_evo.dat metadata.json; do
  if [[ -f "$INIT_DIR/$f" ]]; then
    cp "$INIT_DIR/$f" "$OUTPUT_DIR/"
    echo "[Collect] Copied $f"
  fi
done

# -----------------------------------------------------------------------------
# Step 5: Post-process (optional)
#   - Default script: postprocess/plot_k_eff.py
#   - Set NO_POST=1 to skip. Set POST_SCRIPT to override the script path.
# -----------------------------------------------------------------------------
if [[ "${NO_POST:-0}" == "1" ]]; then
  echo "[Post] Skipping post-processing because NO_POST=1"
else
  POST_SCRIPT="${POST_SCRIPT:-postprocess/plot_k_eff.py}"
  if [[ -f "$POST_SCRIPT" ]]; then
    echo "[Post] python3 \"$POST_SCRIPT\" \"$OUTPUT_DIR\""
    python3 "$POST_SCRIPT" "$OUTPUT_DIR"
  else
    echo "[Post] WARNING: Post-process script not found: $POST_SCRIPT (skipping)"
  fi
fi

# Total runtime reporting
__t1__=$(date +%s)
echo "[Done] Completed in $(( __t1__ - __t0__ )) s"
