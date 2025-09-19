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

set -euo pipefail

# ------------ helpers ------------
log()  { printf "[INFO] %s\n" "$*"; }
warn() { printf "[WARN] %s\n" "$*"; }
err()  { printf "[ERROR] %s\n" "$*" 1>&2; }
die()  { err "$*"; exit 1; }
require_cmd() { command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"; }

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Options (env vars also supported):
  -f, --filename PATH   Grain input path relative to \$input_dir (default: env/var)
  -r, --read {0|1}      readFlag; 1=read grains, 0=generate (default: 1)
  -t, --temp C          Temperature in °C (default: -2.0)
  -h, --humidity H      Relative humidity [0-1] (default: 0.95)
  -p, --procs N         MPI ranks (default: 12)
  --no-build            Skip rebuilding the executable
  --dry-run             Print actions without running the simulation
  -?, --help            Show this help and exit
EOF
}

# light arg parse (zsh/bash compatible)
NO_BUILD=0
DRY_RUN=0
if [[ ${#@} -gt 0 ]]; then
  while [[ ${#} -gt 0 ]]; do
    case "${1}" in
      -f|--filename) filename="${2:-}"; shift 2;;
      -r|--read)     readFlag="${2:-}"; shift 2;;
      -t|--temp)     temp="${2:-}"; shift 2;;
      -h|--humidity) humidity="${2:-}"; shift 2;;
      -p|--procs)    NUM_PROCS="${2:-}"; shift 2;;
      --no-build)    NO_BUILD=1; shift;;
      --dry-run)     DRY_RUN=1; shift;;
      -\?|--help|-help|-h) usage; exit 0;;
      --) shift; break;;
      *) warn "Ignoring unknown arg: ${1}"; shift;;
    esac
  done
fi

# trap helpful errors
TRAPERR() { err "An error occurred on line $1"; }
set -E
trap 'TRAPERR $LINENO' ERR

# =======================================
# User configuration (override via env or edit here)
# =======================================
# Paths
# Defaults (already in your script)
BASE_DIR="${BASE_DIR:-${PETIGA_DIR}/projects/dry_snow_metamorphism}"
input_dir="${input_dir:-$BASE_DIR/inputs}"
filename="${filename:-grains__phi=0.24__Lxmm=2__Lymm=2__seed=22/grains.dat}"

# --- OVERRIDES FROM BATCH (do not change anything else) ---
# If batch passed a full path inputFile, use it and keep it.
if [[ -n "${inputFile:-}" ]]; then
  input_dir="$(cd "$(dirname "$inputFile")" && pwd)"
  filename="$(basename "$inputFile")"
fi

# Derived inputFile: only set if not already provided by batch
: "${inputFile:="$input_dir/$filename"}"


output_dir="${output_dir:-/Users/jacksonbaglino/SimulationResults/dry_snow_metamorphism/scratch}"
exec_file="${exec_file:-$BASE_DIR/dry_snow_metamorphism}"

# Grain selection (relative path under $input_dir)
filename="${filename:-grains__phi=0.24__Lxmm=2__Lymm=2__seed=22/grains.dat}"
readFlag=${readFlag:-1}   # 1=read grains from file; 0=procedural generation (not used here)

# Physics & numerics
delt_t=${delt_t:-1.0e-4}
t_final=${t_final:-$((28 * 24 * 60 * 60))}  # 28 days in seconds
n_out=${n_out:-100}
# t_final=${t_final:-$((2 * 60 * 60))}  # 2 hours in seconds
# n_out=${n_out:-10}
humidity=${humidity:-0.98}
temp=${temp:--25.0}
grad_temp0X=${grad_temp0X:-0.0}
grad_temp0Y=${grad_temp0Y:-3.0e-6}
grad_temp0Z=${grad_temp0Z:-0.0}
dim=${dim:-2}

# Parallel / MPI
NUM_PROCS=${NUM_PROCS:-12}

#
# Derived (preserve batch-passed absolute inputFile if provided)
if [[ -z "${inputFile:-}" ]]; then
  inputFile="$input_dir/$filename"
fi

[[ -d "$input_dir" ]] || die "input_dir does not exist: $input_dir"
mkdir -p "$output_dir" || die "Could not create output_dir: $output_dir"

# =======================================
# Base directory and file setup
# =======================================
# --- Basic validation (only when reading a grain file) ---
if [[ "${readFlag:-1}" -eq 1 ]]; then
  [[ -n "${inputFile:-}" ]] || die "readFlag=1 but 'inputFile' is unset."
  [[ -f "$inputFile" ]]     || die "Input file not found: $inputFile"
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

if [[ ${readFlag:-1} -eq 1 ]]; then
    log "Reading grain file: $inputFile"
else
    log "Generating grains instead of reading from file."
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
if [[ "$NO_BUILD" -eq 0 ]]; then
  require_cmd make
  log "Building executable..."
  if [[ "$DRY_RUN" -eq 1 ]]; then
    log "(dry-run) make dry_snow_metamorphism"
  else
    make dry_snow_metamorphism || die "Build failed. Please check the Makefile and dependencies."
  fi
else
  warn "Skipping build as requested (--no-build)."
fi

# Generate env file if missing (infer Lx/Ly from metadata under sibling 'meta' dir)
if [ ! -f "$SETTINGS_FILE" ]; then
    log ".env not found at $SETTINGS_FILE; generating via scripts/generate_env_from_input.py ..."
    if [[ "$DRY_RUN" -eq 1 ]]; then
      log "(dry-run) python3 scripts/generate_env_from_input.py \"$inputFile\" \"$SETTINGS_FILE\""
    else
      python3 scripts/generate_env_from_input.py "$inputFile" "$SETTINGS_FILE" || die "Failed to generate $SETTINGS_FILE"
    fi
fi

# Copy settings file into the results folder now that it exists
if [ -f "$SETTINGS_FILE" ]; then
    if [[ "$DRY_RUN" -eq 1 ]]; then
      log "(dry-run) cp \"$SETTINGS_FILE\" \"$folder\""
    else
      cp "$SETTINGS_FILE" "$folder"
    fi
fi

# Load run-time physical/mesh parameters from the settings .env

set -a
source "$SETTINGS_FILE"
set +a

log "Settings file loaded: $SETTINGS_FILE"
log "Lx = $Lx  Ly = $Ly Lz = $Lz"
log "----------------------------------------------"

# Ensure grid variables are defined; if missing, generate and reload
typeset -a MISSING_VARS
MISSING_VARS=()
for v in Nx Ny Nz eps; do
  if [[ -z "${(P)v:-}" ]]; then
    MISSING_VARS+="$v"
  fi
done

if (( ${#MISSING_VARS} > 0 )); then
    warn "Missing variables in .env: ${MISSING_VARS[*]}"
    log  "Populating via scripts/generate_env_from_input.py ..."
    if [[ "$DRY_RUN" -eq 1 ]]; then
      log "(dry-run) python3 scripts/generate_env_from_input.py \"$inputFile\" \"$SETTINGS_FILE\" \"$Lx\" \"$Ly\""
    else
      python3 scripts/generate_env_from_input.py "$inputFile" "$SETTINGS_FILE" "$Lx" "$Ly" || die "Failed to generate missing variables in $SETTINGS_FILE"
    fi

    set -a; source "$SETTINGS_FILE"; set +a

    for v in Nx Ny Nz eps; do
      if [[ -z "${(P)v:-}" ]]; then
        die "Variable $v is still undefined in $SETTINGS_FILE after generation."
      fi
    done
    log "[OK] .env updated: Nx=$Nx Ny=$Ny Nz=$Nz eps=$eps"
else
    log "[OK] .env already defines grid variables: Nx=$Nx Ny=$Ny Nz=$Nz eps=$eps"
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
    require_cmd jq

    local grain_dir src_meta dst_meta augment tmp
    grain_dir="$(dirname "$inputFile")"
    src_meta="$grain_dir/metadata.json"
    dst_meta="$folder/metadata.json"

    [[ -f "$src_meta" ]] || die "Source metadata not found: $src_meta"

    if [[ "$DRY_RUN" -eq 1 ]]; then
      log "(dry-run) cp \"$src_meta\" \"$dst_meta\""
    else
      cp "$src_meta" "$dst_meta"
    fi

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
        '{ dsm_run: { run_time_utc: $run_time, executed_on: $host, user: $user,
           inputs: { grains_dat: $input, env_file: $env_path },
           parameters: { sim_dimension: $dim, delt_t: $dt, t_final: $tf, n_out: $nout,
                         temperature_C: $temp, humidity: $hum, grad_temp: { x: $gx, y: $gy, z: $gz } },
           domain_size_m: { Lx: $Lx_v, Ly: $Ly_v, Lz: $Lz_v },
           mesh_resolution: { Nx: $Nx_v, Ny: $Ny_v, Nz: $Nz_v },
           interface_width_eps: $eps_v } }')

    tmp="$(mktemp)"
    if [[ "$DRY_RUN" -eq 1 ]]; then
      log "(dry-run) jq --slurpfile add <augment> '. * $add' \"$dst_meta\" > \"$tmp\" && mv \"$tmp\" \"$dst_meta\""
      rm -f "$tmp"
    else
      jq --argjson add "$augment" '. * $add' "$dst_meta" > "$tmp" && mv "$tmp" "$dst_meta"
      log "[OK] metadata augmented: $dst_meta"
    fi
}

# Write a fully-resolved .env-style snapshot for reproducibility
write_env_snapshot() {
    snapshot="$folder/resolved_params.env"
    if [[ "$DRY_RUN" -eq 1 ]]; then
      log "(dry-run) write snapshot to $snapshot"
      return 0
    fi
    {
      echo "# Auto-generated resolved parameters for this run"
      echo "# Generated: $(date -u +"%Y-%m-%dT%H:%M:%SZ")"
      for var in folder inputFile title Lx Ly Lz Nx Ny Nz delt_t t_final n_out humidity temp grad_temp0X grad_temp0Y grad_temp0Z dim eps readFlag; do
          echo "$var=${(P)var}"
      done
    } > "$snapshot"
}

# Also write machine-readable metadata and a resolved .env snapshot
write_metadata_json
write_env_snapshot

# =======================================
# Run the simulation
# =======================================

# Solver tolerances and options are set here
log "Launching DRY SNOW METAMORPHISM simulation..."
require_cmd mpiexec
if [[ "$DRY_RUN" -eq 1 ]]; then
  log "(dry-run) mpiexec -np \"$NUM_PROCS\" \"$exec_file\" -initial_PFgeom [..opts..] | tee \"$folder/outp.txt\""
else
  mpiexec -np "$NUM_PROCS" "$exec_file" -initial_PFgeom \
    -snes_rtol 1e-3 -snes_stol 1e-6 -snes_max_it 7 \
    -ksp_gmres_restart 150 -ksp_max_it 1000 \
    -ksp_converged_reason -snes_converged_reason -snes_linesearch_monitor \
    -snes_linesearch_type basic | tee "$folder/outp.txt"
fi

# =======================================
# Finalize
# =======================================
log "Simulation completed."
if [[ "$DRY_RUN" -eq 1 ]]; then
  log "(dry-run) cp -r src scripts/studio/run_dsm.sh postprocess/plotDSM.py postprocess/plotSSA.py postprocess/plotPorosity.py \"$folder\""
else
  cp -r src scripts/studio/run_dsm.sh postprocess/plotDSM.py postprocess/plotSSA.py postprocess/plotPorosity.py "$folder"
fi
log "Copied files to $folder"

log "Running plotting scripts (if available)..."
if [[ -x ./scripts/run_plotDSM.sh ]]; then
  if [[ "$DRY_RUN" -eq 1 ]]; then
    log "(dry-run) ./scripts/run_plotDSM.sh"
  else
    ./scripts/run_plotDSM.sh
  fi
else
  warn "./scripts/run_plotDSM.sh not found or not executable; skipping."
fi

log "Simulation complete. Results stored in:"
printf "%s\n" "$folder"