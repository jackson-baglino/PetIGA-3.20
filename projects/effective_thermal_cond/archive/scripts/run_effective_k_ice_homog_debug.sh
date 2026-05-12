#!/bin/zsh
###############################################################################
# Script: run_effective_k_ice_homog_debug.sh
# Purpose:
#   DEBUG-ORIENTED runner for the homogeneous effective thermal conductivity
#   simulation (PetIGA/PETSc). Compiles with debug flags, adds verbose tracing,
#   enables PETSc debug options, captures detailed logs, and runs a quick plot.
#
# Usage:
#   ./scripts/run_effective_k_ice_homog_debug.sh
#   (override any exports below as needed, or export in your shell)
###############################################################################

# --- Shell debug hygiene ------------------------------------------------------
set -Eeuo pipefail
setopt PROMPT_SUBST
# Use prompt escapes for timestamp to avoid relying on $EPOCHREALTIME under set -u
export PS4='[%D{%H:%M:%S}] %N:%i > '
set -x

start_time=$(date +%s)

# ----------------------------
# 🔹  Simulation parameters (EDIT ME)
# ----------------------------
export Nx=${Nx:-150}
export Ny=${Ny:-150}
export Nz=${Nz:-150}          # 1 for 2-D

export Lx=${Lx:-0.20e-03}
export Ly=${Ly:-0.20e-03}
export Lz=${Lz:-0.20e-03}     # ignored when dim=2

export eps=${eps:-9.09629658751972e-07}
export dim=${dim:-2}          # 2 = 2-D, 3 = 3-D
# Required by the executable at startup
export OUTPUT_VTK=${OUTPUT_VTK:-1}
export OUTPUT_BINARY=${OUTPUT_BINARY:-1}
export SOL_INDEX=${SOL_INDEX:-1}

# Boundary / physics knobs (recorded into provenance; solver reads its own opts)
FLUX_BOTTOM=${FLUX_BOTTOM:--0.1}
TEMP_TOP=${TEMP_TOP:-$((273.15-30))}

# Initial ice-field mode: "circle" | "layered" | "FILE"
INIT_MODE=${INIT_MODE:-layered}
INIT_DIR=${INIT_DIR:-}         # required only when INIT_MODE=FILE
ENV_FILE=${ENV_FILE:-}         # optional explicit .env path when INIT_MODE=FILE

# Output roots
OUT_ROOT=${OUT_ROOT:-/Users/jacksonbaglino/SimulationResults/effective_thermal_cond/scratch}

# MPI ranks
NUM_PROCS=${NUM_PROCS:-1}

# PETSc debug options (tweak freely). Applied to ALL ranks.
# Note: keep concise to avoid massive logs; add -info if needed.
PETSC_DEBUG_OPTS=${PETSC_DEBUG_OPTS:--malloc_debug -malloc_test -on_error_abort -snes_converged_reason -ksp_converged_reason}

# --- Project root & build target ---------------------------------------------
BASE_DIR=${BASE_DIR:-${PETIGA_DIR:-/Users/jacksonbaglino/PetIGA-3.20}/projects/effective_thermal_cond}
cd "$BASE_DIR"

# ----------------------------
# 🔹  Optional parameter import from metadata.json / .env
# ----------------------------
METADATA_JSON="${INIT_DIR:+$INIT_DIR/metadata.json}"

if [[ -n "$METADATA_JSON" && -f "$METADATA_JSON" ]]; then
  echo "[info] Using parameters from: $METADATA_JSON"
  json_get() {
    # Usage: json_get dotted.key.path
    python3 - "$METADATA_JSON" "$1" <<'PY'
import json,sys
path = sys.argv[2].split('.')
with open(sys.argv[1]) as f:
    d = json.load(f)
for p in path:
    d = d[p]
print(d)
PY
  }
  : ${dim:=$(json_get sim_dimension)}
  : ${Nx:=$(json_get mesh_resolution.Nx)}
  : ${Ny:=$(json_get mesh_resolution.Ny)}
  : ${Nz:=$(json_get mesh_resolution.Nz)}
  : ${Lx:=$(json_get domain_size_m.Lx)}
  : ${Ly:=$(json_get domain_size_m.Ly)}
  : ${Lz:=$(json_get domain_size_m.Lz)}
  : ${eps:=$(json_get interface_width_eps)}
elif [[ "$INIT_MODE" == "FILE" ]]; then
  if [[ -n "$ENV_FILE" ]]; then
    echo "[info] Using env file (explicit): $ENV_FILE"
    source "$ENV_FILE"
  elif [[ -n "$INIT_DIR" ]]; then
    echo "[info] Searching for single .env in: $INIT_DIR"
    typeset -a env_files
    env_files=( $(find "$INIT_DIR" -maxdepth 1 -type f -name '*.env' 2>/dev/null) )
    if (( ${#env_files[@]} == 1 )); then
      ENV_FILE="${env_files[1]}"; echo "[info] Using env file (found): $ENV_FILE"; source "$ENV_FILE"
    else
      echo "[error] Need exactly one .env in $INIT_DIR (found ${#env_files[@]})."; exit 1
    fi
  else
    echo "[error] INIT_MODE=FILE but neither ENV_FILE nor INIT_DIR provided."; exit 1
  fi
else
  echo "[info] No metadata/.env import needed for procedural mode ($INIT_MODE)."
fi

# Defaults if still empty
: ${dim:=2}
: ${Nx:=150} ; : ${Ny:=150} ; : ${Nz:=150}
: ${Lx:=0.20e-03} ; : ${Ly:=0.20e-03} ; : ${Lz:=0.20e-03}
: ${eps:=9.09629658751972e-07}

# Export variables for simulation
export Nx Ny Nz Lx Ly Lz eps dim TEMP_TOP FLUX_BOTTOM

# ----------------------------
# 🔹  Output folder synthesis (handles layered/circle too)
# ----------------------------
if [[ -n "$INIT_DIR" ]]; then
  base_folder=$(basename "$INIT_DIR")
else
  ts=$(date +%Y%m%dT%H%M%S)
  Lx_mm=$(python3 - "$Lx" <<'PY'
import sys
try:
    v=float(sys.argv[1]); print(f"{v*1e3:g}")
except Exception:
    print("NA")
PY
  )
  
  Ly_mm=$(python3 - "$Ly" <<'PY'
import sys
try:
    v=float(sys.argv[1]); print(f"{v*1e3:g}")
except Exception:
    print("NA")
PY
  )
  mode_tag=${INIT_MODE:-unknown}
  base_folder="homog_${mode_tag}__Lxmm=${Lx_mm}__Lymm=${Ly_mm}__dim=${dim}__${ts}"
fi

OUT_ROOT=${OUT_ROOT%/}
export OUTPUT_DIR="$OUT_ROOT/$base_folder"
mkdir -p "$OUTPUT_DIR"

# record resolved parameters (helps when debugging)
{
  echo "# resolved debug params"
  echo "when=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "pwd=$(pwd)"
  echo "INIT_MODE=$INIT_MODE"
  echo "INIT_DIR=${INIT_DIR:-}"
  echo "ENV_FILE=${ENV_FILE:-}"
  echo "dim=$dim Nx=$Nx Ny=$Ny Nz=$Nz"
  echo "Lx=$Lx Ly=$Ly Lz=$Lz eps=$eps"
  echo "TEMP_TOP=$TEMP_TOP FLUX_BOTTOM=$FLUX_BOTTOM"
  echo "NUM_PROCS=$NUM_PROCS"
  echo "PETSC_DEBUG_OPTS=$PETSC_DEBUG_OPTS"
} > "$OUTPUT_DIR/run.debug.env"

# ----------------------------
# 🔹  Build (DEBUG)
# ----------------------------
compile_code() {
  echo "Compiling (debug)…"
  make BUILD=debug effective_k_ice_homog
}

# ----------------------------
# 🔹  Run with PETSc debug options
# ----------------------------
run_simulation() {
  printf '\nRunning with %s MPI proc(s)…\n' "$NUM_PROCS"
  # Rank-labeled stdout; keep complete logs
  local log_out="$OUTPUT_DIR/run.stdout"
  local log_err="$OUTPUT_DIR/run.stderr"

  # PETSc options can be passed via environment (affects all ranks)
  export PETSC_OPTIONS="$PETSC_DEBUG_OPTS"

  if [[ -n "$INIT_DIR" ]]; then
    mpiexec -l -np $NUM_PROCS ./effective_k_ice_homog \
      -init_mode "$INIT_MODE" \
      -init_dir "$INIT_DIR" \
      1> >(tee "$log_out") 2> >(tee "$log_err" >&2)
  else
    mpiexec -l -np $NUM_PROCS ./effective_k_ice_homog \
      -init_mode "$INIT_MODE" \
      1> >(tee "$log_out") 2> >(tee "$log_err" >&2)
  fi
}

# ----------------------------
# 🔹  Collect outputs
# ----------------------------
collect_outputs() {
  echo "Moving output files to $OUTPUT_DIR"
  mv *.dat  "$OUTPUT_DIR" 2>/dev/null || true
  mv *.bin  "$OUTPUT_DIR" 2>/dev/null || true
  mv *.info "$OUTPUT_DIR" 2>/dev/null || true
  mv *.csv  "$OUTPUT_DIR" 2>/dev/null || true
}

# ----------------------------
# 🔹  Post-process (vector field heatmaps + k_eff)
# ----------------------------
postprocess() {
  echo "Post-processing…"
  python3 postprocess/plot_vector_field.py "$OUTPUT_DIR" "$Nx" "$Ny" "$Lx" "$Ly" || {
    echo "[warn] postprocess/plot_vector_field.py failed" >&2
  }
}

# ----------------------------
# 🔹  Main
# ----------------------------
compile_code
run_simulation
collect_outputs
postprocess

end_time=$(date +%s)
elapsed=$(( end_time - start_time ))

echo "Finished. Outputs in: $OUTPUT_DIR (elapsed ${elapsed}s)"
