#!/bin/zsh
###############################################################################
# Script: run_dsm.sh
# Model: Dry Snow Metamorphism (DSM)
# Purpose:
#   Configure and launch a single DSM simulation, manage run metadata, and
#   archive inputs/plots to a timestamped results folder.
#
# Typical Usage:
#   1) Set `filename` to the desired grain input (in $BASE_DIR/inputs), or set
#      `readFlag=0` to generate grains instead of reading a file.
#   2) Adjust physical & numerical params below (t_final, n_out, gradients, etc.).
#   3) Run: `./scripts/studio/run_dsm.sh`
#
# Key Variables:
#   BASE_DIR      Project root for DSM model
#   input_dir     Directory with input grain files (*.dat)
#   output_dir    Parent directory to store run outputs (timestamped subfolder)
#   exec_file     DSM executable to run
#   filename      (When readFlag=1) The grain input file name in input_dir
#   readFlag      1=read grain file; 0=procedural generation
#   t_final       Total simulated time [s]
#   n_out         Number of output frames (approx.)
#   grad_temp0*   Temperature gradient components [K/m]
#   dim           2D or 3D
#   NUM_PROCS     MPI ranks (default set below)
#
# Outputs:
#   - Creates $output_dir/<title>_<timestamp>/
#   - Copies input & env files used, code snapshot, and plotting scripts
#   - Saves outp.txt, simulation_parameters.csv, sim_params.dat, metadata.json
#
# Notes:
#   - This script avoids deleting any user comments.
#   - If `readFlag=1`, ensure `filename` is set to a valid file in $input_dir.
#   - Settings file is loaded from $BASE_DIR/configs/${filename%.dat}.env.
###############################################################################

# =======================================
# User configuration (override via env or edit here)
# =======================================
# Paths
BASE_DIR="${BASE_DIR:-${PETIGA_DIR}/projects/dry_snow_metamorphism}"
input_dir="${input_dir:-$BASE_DIR/inputs}"
output_dir="${output_dir:-/Users/jacksonbaglino/SimulationResults/dry_snow_metamorphism/scratch}"
exec_file="${exec_file:-$BASE_DIR/dry_snow_metamorphism}"

# Grain selection (relative path under $input_dir)
filename="${filename:-grains__phi=0.24__Lxmm=1__Lymm=1__seed=7/grains.dat}"
readFlag=${readFlag:-1}   # 1=read grains from file; 0=procedural generation (not used here)

# Physics & numerics
delt_t=${delt_t:-1.0e-4}
# t_final=${t_final:-$((28 * 24 * 60 * 60))}  # 28 days in seconds
# n_out=${n_out:-100}
t_final=${t_final:-$((1))}  # 1 second for quick test
n_out=${n_out:-10}
humidity=${humidity:-0.95}
temp=${temp:--2.0}
grad_temp0X=${grad_temp0X:-0.0}
grad_temp0Y=${grad_temp0Y:-3.0e-6}
grad_temp0Z=${grad_temp0Z:-0.0}
dim=${dim:-2}

# Parallel / MPI
NUM_PROCS=${NUM_PROCS:-12}

# Derived
inputFile="$input_dir/$filename"

# =======================================
# Base directory and file setup
# =======================================

# =======================================
# Define simulation parameters
# =======================================
# filename="grainReadFile-2G_Molaro_0p25R1_HIGHRES.dat"
# Batch override precedence: if $filename is exported, use it; else use this default
# : ${filename:="grains__phi=0.24__Lxmm=1__Lymm=1__seed=7/grains.dat"}
# inputFile="$input_dir/$filename"

# --- Basic validation (only when reading a grain file) ---
if [[ "${readFlag:-1}" -eq 1 ]]; then
  if [[ -z "$filename" ]]; then
    echo "[ERROR] readFlag=1 but 'filename' is not set. Set it above (e.g., filename=\"grainReadFile-...dat\")."
    exit 1
  fi
  if [[ ! -f "$inputFile" ]]; then
    echo "[ERROR] Input file not found: $inputFile"
    exit 1
  fi
fi

# readFlag=1  # Set to 1 to read grain file, 0 to generate grains

# delt_t=1.0e-4
# t_final=$((28 * 24 * 60 * 60))  # 28 days in seconds
# # t_final=$((12 * 60 * 60))  # 2 hours in seconds
# # t_final=10
# n_out=100
# # Batch override precedence: use exported values if provided; otherwise defaults
# : ${humidity:=0.95}
# : ${temp:=-2.0}
# grad_temp0X=0.0
# grad_temp0Y=3.0e-6
# grad_temp0Z=0.0
# dim=2

if [[ readFlag -eq 1 ]]; then
    echo "[INFO] Reading grain file: $inputFile"
else
    echo "[INFO] Generating grains instead of reading from file."
fi

# Build a descriptive run title encoding key parameters for easier indexing
clean_name="${filename#grainReadFile-}"
clean_name="${clean_name%.dat}"

# Create compact two-digit tags for temperature (°C) and humidity (%)
# temp_tag: round to nearest integer, then take sign + first two digits
temp_int=$(printf "%.0f" "$temp")
if [[ "$temp_int" == -* ]]; then
  temp_tag=${temp_int:0:3}   # e.g., -12, -9
else
  temp_tag=${temp_int:0:2}   # e.g., 15 -> 15, 7 -> 7
fi
# hum_tag: integer percent, then first two digits (e.g., 98, 95, 100 -> 10)
hum_int=$(awk "BEGIN{printf \"%d\", $humidity*100}")
hum_tag=$(printf "%02d" "$hum_int")
hum_tag=${hum_tag:0:2}

# Build title using compact tags
# Note: keep existing day count in suffix
ndays=$(awk "BEGIN{printf \"%d\", $t_final/86400}")

title="DSM${clean_name}_${dim}D_Tm${temp_tag}_hum${hum_tag}_tf${ndays}d_"
# SETTINGS_FILE="$BASE_DIR/configs/${filename%.dat}.env"
SETTINGS_FILE="$BASE_DIR/inputs/${filename%.dat}.env"

# MPI ranks used for this run (override here or export before calling)
# NUM_PROCS=12  # Number of MPI processes

# =======================================
# Timestamped result folder
# =======================================
timestamp=$(date +%Y-%m-%d__%H.%M.%S)
folder="$output_dir/${title}_$timestamp"
mkdir -p "$folder"

# Copy input file to results folder
cp "$inputFile" "$folder"

# =======================================
# Build and run setup
# =======================================
cd "$BASE_DIR" || exit 1

# Compile the DSM executable (requires PETSc/PetIGA set up)
make dry_snow_metamorphism || {
    echo "[ERROR] Build failed. Please check the Makefile and dependencies."
    exit 1
}

# Generate env file if missing (infer Lx/Ly from metadata under sibling 'meta' dir)
if [ ! -f "$SETTINGS_FILE" ]; then
    echo "[INFO] .env not found at $SETTINGS_FILE; generating it with scripts/generate_env_from_input.py ..."
    python3 scripts/generate_env_from_input.py "$inputFile" "$SETTINGS_FILE" || {
        echo "[ERROR] Failed to generate $SETTINGS_FILE"; exit 1; }
fi

# Copy settings file into the results folder now that it exists
if [ -f "$SETTINGS_FILE" ]; then
    cp "$SETTINGS_FILE" "$folder"
fi

# Load run-time physical/mesh parameters from the settings .env

set -a
source "$SETTINGS_FILE"
set +a

echo "Settings file loaded: $SETTINGS_FILE"
echo "Lx = $Lx  Ly = $Ly Lz = $Lz"
echo "----------------------------------------------"

# Ensure grid variables are defined; if missing, generate and reload
typeset -a MISSING_VARS
MISSING_VARS=()
for v in Nx Ny Nz eps; do
  if [[ -z "${(P)v:-}" ]]; then
    MISSING_VARS+="$v"
  fi
done

if (( ${#MISSING_VARS} > 0 )); then
    echo "[INFO] Missing variables in .env: ${MISSING_VARS[*]}"
    echo "[INFO] Running scripts/generate_env_from_input.py to populate them..."
    python3 scripts/generate_env_from_input.py "$inputFile" "$SETTINGS_FILE" "$Lx" "$Ly" || {
        echo "[ERROR] Failed to generate missing variables in $SETTINGS_FILE"; exit 1; }

    # Re-load and re-validate
    set -a
    source "$SETTINGS_FILE"
    set +a

    for v in Nx Ny Nz eps; do
      if [[ -z "${(P)v:-}" ]]; then
        echo "[ERROR] Variable $v is still undefined in $SETTINGS_FILE after generation.";
        exit 1
      fi
    done
    echo "[OK] .env updated: Nx=$Nx Ny=$Ny Nz=$Nz eps=$eps"
else
    echo "[OK] .env already defines grid variables: Nx=$Nx Ny=$Ny Nz=$Nz eps=$eps"
fi

# Export simulation parameters
export folder input_dir inputFile filename title
export delt_t t_final n_out humidity temp grad_temp0X grad_temp0Y grad_temp0Z dim
export readFlag Lx Ly Lz Nx Ny Nz eps

# =======================================
# Save simulation metadata
# =======================================

# Function to write metadata.json: copy from input folder, then append DSM run info
write_metadata_json() {
    # Require jq for safe JSON merge
    if ! command -v jq >/dev/null 2>&1; then
        echo "[ERROR] jq is required to write/augment metadata.json" >&2
        return 1
    fi

    # Source (from grain input folder) and destination (output folder)
    local grain_dir
    grain_dir="$(dirname "$inputFile")"
    local src_meta="$grain_dir/metadata.json"
    local dst_meta="$folder/metadata.json"

    if [[ ! -f "$src_meta" ]]; then
        echo "[ERROR] Source metadata not found: $src_meta" >&2
        return 1
    fi

    # Copy the original packing metadata
    cp "$src_meta" "$dst_meta"

    # Build a DSM run augmentation object
    local augment
    augment=$(jq -n \
        --arg run_time "$(date -u +"%Y-%m-%dT%H:%M:%SZ")" \
        --arg host     "$(hostname)" \
        --arg user     "$(whoami)" \
        --arg input    "$inputFile" \
        --arg env_path "$SETTINGS_FILE" \
        --argjson dim  "${dim:-null}" \
        --argjson dt   "${delt_t:-null}" \
        --argjson tf   "${t_final:-null}" \
        --argjson nout "${n_out:-null}" \
        --argjson temp "${temp:-null}" \
        --argjson hum  "${humidity:-null}" \
        --argjson gx   "${grad_temp0X:-null}" \
        --argjson gy   "${grad_temp0Y:-null}" \
        --argjson gz   "${grad_temp0Z:-null}" \
        --argjson Lx_v "${Lx:-null}" \
        --argjson Ly_v "${Ly:-null}" \
        --argjson Lz_v "${Lz:-null}" \
        --argjson Nx_v "${Nx:-null}" \
        --argjson Ny_v "${Ny:-null}" \
        --argjson Nz_v "${Nz:-null}" \
        --argjson eps_v "${eps:-null}" \
        '{
           dsm_run: {
             run_time_utc: $run_time,
             executed_on: $host,
             user: $user,
             inputs: { grains_dat: $input, env_file: $env_path },
             parameters: {
               sim_dimension: $dim,
               delt_t: $dt, t_final: $tf, n_out: $nout,
               temperature_C: $temp, humidity: $hum,
               grad_temp: { x: $gx, y: $gy, z: $gz }
             },
             domain_size_m: { Lx: $Lx_v, Ly: $Ly_v, Lz: $Lz_v },
             mesh_resolution: { Nx: $Nx_v, Ny: $Ny_v, Nz: $Nz_v },
             interface_width_eps: $eps_v
           }
         }')

    # Merge augmentation into the copied metadata
    local tmp
    tmp="$(mktemp)"
    jq --argjson add "$augment" '. * $add' "$dst_meta" > "$tmp" && mv "$tmp" "$dst_meta"
    echo "[OK] metadata augmented: $dst_meta"
}

# Write a fully-resolved .env-style snapshot for reproducibility
write_env_snapshot() {
    snapshot="$folder/resolved_params.env"
    echo "# Auto-generated resolved parameters for this run" > "$snapshot"
    echo "# Generated: $(date -u +"%Y-%m-%dT%H:%M:%SZ")" >> "$snapshot"
    for var in folder inputFile title Lx Ly Lz Nx Ny Nz delt_t t_final n_out humidity temp grad_temp0X grad_temp0Y grad_temp0Z dim eps readFlag; do
        echo "$var=${(P)var}" >> "$snapshot"
    done
}

# Also write machine-readable metadata and a resolved .env snapshot
write_metadata_json
write_env_snapshot

# =======================================
# Run the simulation
# =======================================

# Solver tolerances and options are set here
echo "[INFO] Launching DRY SNOW METAMORPHISM simulation..."
mpiexec -np "$NUM_PROCS" "$exec_file" -initial_PFgeom \
  -snes_rtol 1e-3 -snes_stol 1e-6 -snes_max_it 7 \
  -ksp_gmres_restart 150 -ksp_max_it 1000 \
  -ksp_converged_reason -snes_converged_reason -snes_linesearch_monitor \
  -snes_linesearch_type basic | tee "$folder/outp.txt"

# =======================================
# Finalize
# =======================================
echo " "
echo "[INFO] Simulation completed."
cp -r src scripts/studio/run_dsm.sh postprocess/plotDSM.py postprocess/plotSSA.py postprocess/plotPorosity.py "$folder"
echo "Copied files to $folder"

echo " "
echo "Running plotting scripts (if available)..."
if [[ -x ./scripts/run_plotDSM.sh ]]; then
  ./scripts/run_plotDSM.sh
else
  echo "[info] ./scripts/run_plotDSM.sh not found or not executable; skipping."
fi

echo "Simulation complete. Results stored in:"
echo "$folder"