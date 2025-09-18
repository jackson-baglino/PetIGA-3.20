#!/bin/zsh
###############################################################################
# Script: run_effective_k_ice_homog.sh
# Purpose:
#   Compile and run the homogeneous effective thermal conductivity simulation
#   (PetIGA/PETSc), collect outputs, and run a quick post-processing plot.
#
# Usage:
#   ./run_effective_k_ice_homog.sh \
#       [with environment variables exported beforehand as needed]
#
# Key env vars (export before running or load via .env):
#   Nx, Ny, Nz           Grid resolution
#   Lx, Ly, Lz           Domain sizes (m)
#   eps                  Interface width
#   dim                  2 for 2D, 3 for 3D (defaults to 2 if unset)
#   TEMP_TOP             Prescribed top temperature (K)
#   FLUX_BOTTOM          Prescribed bottom heat flux (W/m^2)
#   INIT_MODE            "circle" | "layered" | "FILE"
#   INIT_DIR             Directory containing initial condition files (.dat) and a .env
#   ENV_FILE             Explicit path to a .env to source (overrides INIT_DIR search)
#   NUM_PROCS            MPI ranks (default: 1)
#
# Outputs:
#   Creates a timestamped folder under outputs/homog/ and moves *.dat, *.bin,
#   *.info, *.csv there. Then runs postprocess/plot_vector_field.py.
#
# Notes:
#   - This script intentionally avoids set -euo pipefail and heavy tracing.
#   - For INIT_MODE=FILE, either INIT_DIR must point to a directory containing
#     a single .env OR ENV_FILE must be provided explicitly.
###############################################################################
###############################################################################
# run_effective_k_ice_full.sh
#
# • No verbose tracing, no pipefail, no trap on every error
# • Compiles with *release* flags (BUILD=release)
# • Runs without PETSc debug monitors
# • Still honours the same environment variables / paths
###############################################################################

start_time=$(date +%s)  # seconds since epoch

# ----------------------------
# 🔹  Simulation parameters
# ----------------------------
export Nx=571
export Ny=571
export Nz=1          # 1 for 2-D

# export Nx=134
# export Ny=214

export Lx=1.00e-03
export Ly=1.00e-03
export Lz=1.00e-03    # ignored when dim=2

FLUX_BOTTOM=-0.1
TEMP_TOP=$((273.15-30))

export eps=$((9.09629658751972e-07))
export dim=2         # 2 = 2-D, 3 = 3-D

# Initial ice-field mode: "circle" | "layered" | /path/to/file.dat
# INIT_MODE="circle"
# INIT_MODE="layered"
# INIT_MODE="FILE"
# INIT_DIR="/Users/jacksonbaglino/PetIGA-3.20/projects/effective_thermal_cond/inputs/"\
# "NASAv2_96G-2D_T-20.0_hum0.70_2025-05-31__18.55.56"
# INIT_DIR="/Users/jacksonbaglino/SimulationResults/dry_snow_metamorphism/scratch/DSM2G_Molaro_0p25R1_2D_Tm12_hum5_tf0d__2025-08-07__18.59.47"

# Prefer JSON parameters from a DSM run directory if available
METADATA_JSON="$INIT_DIR/metadata.json"

if [[ ! -f "$METADATA_JSON" ]]; then
# Locate and source .env file (portable across bash/zsh)
# Priority:
#  1) If ENV_FILE is set, use it directly
#  2) Else if INIT_MODE=="FILE" and INIT_DIR is set, search for a single .env in INIT_DIR
#  3) Else, if INIT_MODE=="FILE" but neither ENV_FILE nor INIT_DIR is provided, fail with guidance

if [[ -n "$ENV_FILE" ]]; then
    echo "✅ Using env file (explicit): $ENV_FILE"
elif [[ "$INIT_MODE" == "FILE" ]]; then
    if [[ -z "$INIT_DIR" ]]; then
        echo "❌ INIT_MODE is 'FILE' but neither ENV_FILE nor INIT_DIR is set."
        echo "   Provide ENV_FILE=/path/to/.env or INIT_DIR=/path/to/dir containing a single .env."
        exit 1
    fi
    if [[ ! -d "$INIT_DIR" ]]; then
        echo "❌ INIT_DIR is not a directory: $INIT_DIR"
        exit 1
    fi
    echo "🔎 Searching for .env in: $INIT_DIR"
    env_files=( $(find "$INIT_DIR" -maxdepth 1 -type f -name '*.env' 2>/dev/null) )
    env_count=${#env_files[@]}
    if [[ -z "$env_count" || "$env_count" -eq 0 ]]; then
        echo "❌ No .env file found in $INIT_DIR"
        exit 1
    elif [[ "$env_count" -gt 1 ]]; then
        echo "❌ Multiple .env files found in $INIT_DIR:";
        for f in "${env_files[@]}"; do echo "  $f"; done
        echo "   Set ENV_FILE=/path/to.env to choose one."
        exit 1
    else
        ENV_FILE="${env_files[1]}"
        echo "✅ Using env file (discovered): $ENV_FILE"
    fi
else
    echo "ℹ️  No ENV_FILE provided and INIT_MODE != 'FILE'. Proceeding without sourcing .env."
fi

# If an env file was resolved, source it
if [[ -n "$ENV_FILE" ]]; then
    if [[ ! -f "$ENV_FILE" ]]; then
        echo "❌ .env file not found: $ENV_FILE"
        exit 1
    fi
    source "$ENV_FILE"
fi
fi

# If metadata.json is present, read parameters from it (preferred)
if [[ -f "$METADATA_JSON" ]]; then
    echo "🧾 Using parameters from: $METADATA_JSON"

    json_get() {
        # Usage: json_get dotted.key.path
        python3 - "$METADATA_JSON" "$1" <<'PY'
import json,sys
path = sys.argv[2].split('.')
with open(sys.argv[1]) as f:
    d = json.load(f)
for p in path:
    d = d[p]
# Print scalar as-is; for safety cast to str
print(d)
PY
    }

    # Populate variables only if they are not already set in the environment
    : ${dim:=$(json_get sim_dimension)}
    : ${Nx:=$(json_get mesh_resolution.Nx)}
    : ${Ny:=$(json_get mesh_resolution.Ny)}
    : ${Nz:=$(json_get mesh_resolution.Nz)}
    : ${Lx:=$(json_get domain_size_m.Lx)}
    : ${Ly:=$(json_get domain_size_m.Ly)}
    : ${Lz:=$(json_get domain_size_m.Lz)}
    : ${eps:=$(json_get interface_width_eps)}
fi

# Set dim default if not provided anywhere
if [[ -z "$dim" ]]; then
    dim=2
fi

# Export variables for simulation
export Nx
export Ny
export Nz
export Lx
export Ly
export Lz
export eps
export TEMP_TOP
export FLUX_BOTTOM
export dim

# Extract or synthesize base folder name for outputs
if [[ -n "$INIT_DIR" ]]; then
    base_folder=$(basename "$INIT_DIR")
else
    # Build a seed-free, reproducible-ish name for procedural modes
    ts=$(date +%Y%m%dT%H%M%S)
    Lx_mm=$(python3 - <<'PY'
import sys
try:
    v=float(sys.argv[1]); print(f"{v*1e3:g}")
except Exception:
    print("NA")
PY
"$Lx")
    Ly_mm=$(python3 - <<'PY'
import sys
try:
    v=float(sys.argv[1]); print(f"{v*1e3:g}")
except Exception:
    print("NA")
PY
"$Ly")
    mode_tag=${INIT_MODE:-unknown}
    base_folder="homog_${mode_tag}__Lxmm=${Lx_mm}__Lymm=${Ly_mm}__dim=${dim}__${ts}"
fi

# Resolve which sol_XXXXX.dat to load; prefer requested SOL_INDEX, else auto-detect
resolve_sol_index() {
    # If INIT_DIR is set, try to locate sol_*.dat files
    if [[ -n "${INIT_DIR:-}" && -d "$INIT_DIR" ]]; then
        # If SOL_INDEX is set, verify the file exists; otherwise try common fallbacks
        if [[ -n "${SOL_INDEX:-}" ]]; then
            local try="$INIT_DIR/sol_$(printf '%05d' "$SOL_INDEX").dat"
            if [[ -f "$try" ]]; then
                echo "✅ Using SOL_INDEX=$SOL_INDEX ($try)"
                return 0
            fi
            echo "⚠️  Requested SOL_INDEX=$SOL_INDEX not found at $try; attempting auto-detect…"
        fi
        # Try 00000 then 00001, then pick the smallest available index
        if [[ -f "$INIT_DIR/sol_00000.dat" ]]; then
            export SOL_INDEX=0; echo "✅ Auto-detected SOL_INDEX=0"; return 0
        fi
        if [[ -f "$INIT_DIR/sol_00001.dat" ]]; then
            export SOL_INDEX=1; echo "✅ Auto-detected SOL_INDEX=1"; return 0
        fi
        # General scan for sol_*.dat; choose the lowest index present
        local first
        first=$(ls "$INIT_DIR"/sol_*.dat 2>/dev/null | sed -E 's#.*/sol_([0-9]{5})\.dat#\1#' | sort | head -n1 || true)
        if [[ -n "$first" ]]; then
            export SOL_INDEX=$((10#$first))
            echo "✅ Auto-detected SOL_INDEX=$SOL_INDEX"
            return 0
        fi
        echo "❌ No sol_XXXXX.dat files found in $INIT_DIR" >&2
        return 1
    fi
    return 0
}

# If INIT_DIR is provided, default INIT_MODE to FILE unless already set
if [[ -n "${INIT_DIR:-}" && -z "${INIT_MODE:-}" ]]; then
    INIT_MODE="FILE"
fi

# If we are in FILE mode, resolve the solution index to avoid missing-file errors
if [[ "${INIT_MODE:-}" == "FILE" ]]; then
    resolve_sol_index || true
fi

grains=2
dims=2
temp=-20.0
hum=1.00

OUT_ROOT="/Users/jacksonbaglino/SimulationResults/effective_thermal_cond/scratch"

echo "✅ Output root: $OUT_ROOT"

echo "Loaded parameters from .env:"
echo "  dim=$dim, Nx=$Nx, Ny=$Ny, Nz=$Nz, Lx=$Lx, Ly=$Ly, Lz=$Lz, eps=$eps, TEMP_TOP=$TEMP_TOP"

# Output flags
export OUTPUT_VTK=1
export OUTPUT_BINARY=1

export OUTPUT_DIR="$OUT_ROOT/${base_folder}"
mkdir -p "$OUTPUT_DIR"

# MPI ranks (override by exporting NUM_PROCS beforehand)
NUM_PROCS=${NUM_PROCS:-1}
# NUM_PROCS=4

# ----------------------------
# 🔹  Helpers
# ----------------------------
compile_code() {
    echo "Compiling (release)…"
    make BUILD=release effective_k_ice_homog
}

run_simulation() {
    echo " "
    echo "Running with $NUM_PROCS MPI proc(s)…"
    if [[ -n "${INIT_DIR:-}" ]]; then
        if [[ -n "${INIT_MODE:-}" ]]; then
            mpiexec -np $NUM_PROCS ./effective_k_ice_homog \
                -init_mode "$INIT_MODE" \
                -init_dir "$INIT_DIR"
        else
            mpiexec -np $NUM_PROCS ./effective_k_ice_homog \
                -init_dir "$INIT_DIR"
        fi
    else
        if [[ -n "${INIT_MODE:-}" ]]; then
            mpiexec -np $NUM_PROCS ./effective_k_ice_homog \
                -init_mode "$INIT_MODE"
        else
            mpiexec -np $NUM_PROCS ./effective_k_ice_homog
        fi
    fi
}

collect_outputs() {
    echo "📂 Moving output files to $OUTPUT_DIR"
    mv *.dat  "$OUTPUT_DIR" 2>/dev/null
    mv *.bin  "$OUTPUT_DIR" 2>/dev/null
    mv *.info "$OUTPUT_DIR" 2>/dev/null
    mv *.csv  "$OUTPUT_DIR" 2>/dev/null
}

# ----------------------------
# 🔹  Main
# ----------------------------
compile_code   || { echo "❌ compile failed"; exit 1; }
run_simulation || { echo "❌ simulation failed"; exit 1; }
collect_outputs

echo "📈 Post-processing…"
python3 postprocess/plot_vector_field.py "$OUTPUT_DIR" "$Nx" "$Ny" "$Lx" "$Ly"

echo "✅ Finished. Outputs in: $OUTPUT_DIR"

end_time=$(date +%s)

elapsed=$(( end_time - start_time ))

echo "Simulation completed in $elapsed seconds."