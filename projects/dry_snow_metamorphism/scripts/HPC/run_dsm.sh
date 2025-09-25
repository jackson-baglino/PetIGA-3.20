#!/bin/bash
#SBATCH -J DSM_Tm-30_hum98
#SBATCH -A rubyfu
#SBATCH -t 0-4:00:00

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --cpus-per-task=1
#SBATCH -o "output_files/%x.o%j"
#SBATCH -e "output_files/%x.e%j"
#SBATCH --partition=expansion
#SBATCH --mem-per-cpu=1G
#SBATCH --mail-user=jbaglino@caltech.edu
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --constraint='icelake|skylake|cascadelake'


###############################################################################
# Script: run_dsm.sh (HPC)
# Model: Dry Snow Metamorphism (DSM)
# Purpose:
#   Configure and launch a single DSM simulation on the cluster, snapshot the
#   fully-resolved parameters used, and write machine-readable metadata.
#
# This mirrors the desktop runner, but uses SLURM-safe launch and paths.
###############################################################################

set -euo pipefail

cleanup_on_err() {
  echo "[ERROR] Aborting due to an error (line $1)." >&2
  if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    scancel "$SLURM_JOB_ID" >/dev/null 2>&1 || true
  fi
}
trap 'cleanup_on_err $LINENO' ERR

# =======================================
# User configuration (edit here or override via environment)
# =======================================
# Paths
: "${PETIGA_DIR:=/Users/jacksonbaglino/PetIGA-3.20}"
BASE_DIR="${BASE_DIR:-${PETIGA_DIR}/projects/dry_snow_metamorphism}"
input_dir="${input_dir:-$BASE_DIR/inputs}"
output_dir="${output_dir:-/resnick/scratch/jbaglino}"
exec_file="${exec_file:-$BASE_DIR/dry_snow_metamorphism}"

# Grain selection (relative path under $input_dir)
filename="${filename:-grains__phi=0.24__Lxmm=3__Lymm=3__seed=21/grains.dat}"
readFlag=${readFlag:-1}   # 1=read grains from file; 0=procedural generation (not used here)

# If batch passed a full path inputFile, prefer it and derive input_dir/filename from it
if [[ -n "${inputFile:-}" ]]; then
  input_dir="$(cd "$(dirname "$inputFile")" && pwd)"
  filename="$(basename "$inputFile")"
fi

# Physics & numerics
delt_t=${delt_t:-1.0e-4}
t_final=${t_final:-$((28 * 24 * 60 * 60))}  # 28 days in seconds
# t_final=${t_final:-$((2 * 60 * 60))}  # 2 hours for quick test
n_out=${n_out:-100}
# n_out=${n_out:-10} # Reduced for quick test
humidity=${humidity:-0.95}
temp=${temp:--30.0}
grad_temp0X=${grad_temp0X:-0.0}
grad_temp0Y=${grad_temp0Y:-3.0e-6}
grad_temp0Z=${grad_temp0Z:-0.0}
dim=${dim:-2}

# Parallel / MPI
# Prefer SLURM-provided task count; fallback to explicit override; final fallback to 40
if [[ -n "${SLURM_NTASKS:-}" ]]; then
  NUM_PROCS="$SLURM_NTASKS"
else
  NUM_PROCS="${NUM_PROCS:-40}"
fi

# Derived (preserve batch-passed absolute inputFile if provided)
: "${inputFile:="$input_dir/$filename"}"

# =======================================
# Define simulation parameters & basic validation
# =======================================
if [[ "${readFlag}" -eq 1 ]]; then
  if [[ -z "${inputFile:-}" ]]; then
    echo "[ERROR] readFlag=1 but 'inputFile' is unset." >&2; exit 1
  fi
  echo "[INFO] Reading grain file: $inputFile"
  if [[ ! -f "$inputFile" ]]; then
    echo "[ERROR] Input file not found: $inputFile" >&2; exit 1
  fi
  echo "[INFO] Using inputFile: $inputFile"
else
  echo "[INFO] Generating grains instead of reading from file."
fi

# Build a descriptive run title encoding key parameters for easier indexing
clean_name="${filename#grainReadFile-}"
clean_name="${clean_name%.dat}"

# Compact two-digit tags for temperature (°C) and humidity (%)
temp_int=$(printf "%.0f" "$temp")
if [[ "$temp_int" == -* ]]; then
  temp_tag=${temp_int:0:3}
else
  temp_tag=${temp_int:0:2}
fi
hum_int=$(awk "BEGIN{printf \"%d\", $humidity*100}")
hum_tag=$(printf "%02d" "$hum_int"); hum_tag=${hum_tag:0:2}

# Inserted code to read phi from metadata.json if it exists
env_dir="$(dirname "$inputFile")"
if [[ -f "$env_dir/metadata.json" ]]; then
  phi=$(jq -r 'porosity_target' "$env_dir/metadata.json")
fi

# Title and settings path
ndays=$(awk "BEGIN{printf \"%d\", $t_final/86400}")
phi_tag=$(awk -v phi="$phi" 'BEGIN{printf "phi%.2f", phi}')
Lx_tag=$(awk -v v="$Lx" 'BEGIN{printf "Lx%.0fmm", v*1000}')
Ly_tag=$(awk -v v="$Ly" 'BEGIN{printf "Ly%.0fmm", v*1000}')
temp_tag=$(printf "Tm%d" "$(printf "%.0f" "$temp")")

title="DSM_${phi_tag}_${Lx_tag}_${Ly_tag}_${temp_tag}_hum${hum_tag}_tf${ndays}d"

# Resolve settings file adjacent to inputFile (prefer grains.env, then basename.env)
if [[ -f "$env_dir/grains.env" ]]; then
  SETTINGS_FILE="$env_dir/grains.env"
elif [[ -f "$env_dir/${filename%.dat}.env" ]]; then
  SETTINGS_FILE="$env_dir/${filename%.dat}.env"
else
  # Target path for generation if none exist yet
  SETTINGS_FILE="$env_dir/grains.env"
fi

# =======================================
# Timestamped result folder
# =======================================
timestamp=$(date +%Y-%m-%d__%H.%M.%S)
folder="$output_dir/${title}_$timestamp"
mkdir -p "$folder"

# Copy input file to results folder (provenance)
cp "$inputFile" "$folder"

# =======================================
# Build and run setup
# =======================================
cd "$BASE_DIR"

# Compile the DSM executable (requires PETSc/PetIGA set up)
make dry_snow_metamorphism || { echo "[ERROR] Build failed." >&2; exit 1; }

# Ensure .env exists (generate if missing)
if [[ ! -f "$SETTINGS_FILE" ]]; then
  echo "[INFO] .env not found at $SETTINGS_FILE; generating via scripts/generate_env_from_input.py ..."
  if command -v srun >/dev/null 2>&1; then
    srun -N1 -n1 --cpus-per-task=1 python3 scripts/generate_env_from_input.py "$inputFile" "$SETTINGS_FILE" || { echo "[ERROR] Failed to generate $SETTINGS_FILE" >&2; exit 1; }
  else
    python3 scripts/generate_env_from_input.py "$inputFile" "$SETTINGS_FILE" || { echo "[ERROR] Failed to generate $SETTINGS_FILE" >&2; exit 1; }
  fi
fi

# Copy settings file into the results folder
cp "$SETTINGS_FILE" "$folder"

# Load run-time physical/mesh parameters from the settings .env
set -a
source "$SETTINGS_FILE"
set +a

echo "Settings file loaded: $SETTINGS_FILE"
echo "Lx = $Lx  Ly = $Ly  Lz = $Lz"
echo "----------------------------------------------"

# Ensure grid variables are defined; if missing, generate and reload
missing=()
for v in Nx Ny Nz eps; do
  eval "val=\${$v:-}"
  [[ -n "$val" ]] || missing+=("$v")
done

if (( ${#missing[@]} > 0 )); then
  echo "[INFO] Missing variables in .env: ${missing[*]}"
  echo "[INFO] Running scripts/generate_env_from_input.py to populate them..."
  if command -v srun >/dev/null 2>&1; then
    srun -N1 -n1 --cpus-per-task=1 python3 scripts/generate_env_from_input.py "$inputFile" "$SETTINGS_FILE" "$Lx" "$Ly" || { echo "[ERROR] Failed to generate missing variables in $SETTINGS_FILE" >&2; exit 1; }
  else
    python3 scripts/generate_env_from_input.py "$inputFile" "$SETTINGS_FILE" "$Lx" "$Ly" || { echo "[ERROR] Failed to generate missing variables in $SETTINGS_FILE" >&2; exit 1; }
  fi
  set -a; source "$SETTINGS_FILE"; set +a
  for v in Nx Ny Nz eps; do
    eval "val=\${$v:-}"
    if [[ -z "$val" ]]; then
      echo "[ERROR] Variable $v is still undefined in $SETTINGS_FILE after generation." >&2; exit 1
    fi
  done
  echo "[OK] .env updated: Nx=$Nx Ny=$Ny Nz=$Nz eps=$eps"
else
  echo "[OK] .env already defines grid variables: Nx=$Nx Ny=$Ny Nz=$Nz eps=$eps"
fi

# Export simulation parameters (for downstream tools / postprocess)
export folder input_dir inputFile filename title
export delt_t t_final n_out humidity temp grad_temp0X grad_temp0Y grad_temp0Z dim
export readFlag Lx Ly Lz Nx Ny Nz eps

# =======================================
# Save simulation metadata (copy + augment)
# =======================================
write_metadata_json() {
  # Require jq for safe JSON merge
  if ! command -v jq >/dev/null 2>&1; then
    echo "[ERROR] jq is required to write/augment metadata.json" >&2
    return 1
  fi
  local grain_dir src_meta dst_meta augment tmp
  grain_dir="$(dirname "$inputFile")"
  src_meta="$grain_dir/metadata.json"
  dst_meta="$folder/metadata.json"
  if [[ ! -f "$src_meta" ]]; then
    echo "[ERROR] Source metadata not found: $src_meta" >&2
    return 1
  fi
  cp "$src_meta" "$dst_meta"
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
    --arg slurm_job   "${SLURM_JOB_ID:-}" \
    --arg slurm_name  "${SLURM_JOB_NAME:-}" \
    --arg slurm_nodes "${SLURM_JOB_NUM_NODES:-}" \
    --arg slurm_ntasks "${SLURM_NTASKS:-}" \
    --arg slurm_cpt   "${SLURM_CPUS_PER_TASK:-}" \
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
         interface_width_eps: $eps_v,
         slurm: { job_id: $slurm_job, job_name: $slurm_name,
                  nodes: $slurm_nodes, ntasks: $slurm_ntasks, cpus_per_task: $slurm_cpt }
       }
     }')
  tmp="$(mktemp)"; jq --argjson add "$augment" '. * $add' "$dst_meta" > "$tmp" && mv "$tmp" "$dst_meta"
  echo "[OK] metadata augmented: $dst_meta"
}

# Write a fully-resolved .env-style snapshot for reproducibility
write_env_snapshot() {
  local snapshot="$folder/resolved_params.env"
  {
    echo "# Auto-generated resolved parameters for this run"
    echo "# Generated: $(date -u +"%Y-%m-%dT%H:%M:%SZ")"
    for var in folder inputFile title Lx Ly Lz Nx Ny Nz delt_t t_final n_out humidity temp grad_temp0X grad_temp0Y grad_temp0Z dim eps readFlag; do
      eval "val=\${$var}"
      echo "$var=$val"
    done
  } > "$snapshot"
}

# Also write machine-readable metadata and a resolved .env snapshot
write_metadata_json
write_env_snapshot

# =======================================
# Run the simulation (SLURM-friendly)
# =======================================

echo "[INFO] Launching DRY SNOW METAMORPHISM simulation..."
# Use srun if available (preferred on SLURM), otherwise mpiexec with -n $NUM_PROCS
if command -v srun >/dev/null 2>&1; then
  set +e
  srun --kill-on-bad-exit=1 "$exec_file" -initial_PFgeom \
    -snes_rtol 1e-3 -snes_stol 1e-6 -snes_max_it 7 \
    -ksp_gmres_restart 150 -ksp_max_it 1000 \
    -ksp_converged_reason -snes_converged_reason -snes_linesearch_monitor \
    -snes_linesearch_type basic 2>&1 | tee "$folder/outp.txt"
  rc=${PIPESTATUS[0]}
  set -e
else
  set +e
  mpiexec -n "$NUM_PROCS" "$exec_file" -initial_PFgeom \
    -snes_rtol 1e-3 -snes_stol 1e-6 -snes_max_it 7 \
    -ksp_gmres_restart 150 -ksp_max_it 1000 \
    -ksp_converged_reason -snes_converged_reason -snes_linesearch_monitor \
    -snes_linesearch_type basic 2>&1 | tee "$folder/outp.txt"
  rc=${PIPESTATUS[0]}
  set -e
fi

# Abort immediately on simulation failure
if [[ ${rc:-0} -ne 0 ]]; then
  echo "[ERROR] Simulation failed with exit code $rc — terminating job." >&2
  # Best-effort cancel of the remaining allocation (if still active)
  if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    scancel "$SLURM_JOB_ID" >/dev/null 2>&1 || true
  fi
  exit "$rc"
fi

# =======================================
# Finalize (only if simulation succeeded)
# =======================================
echo
echo "[INFO] Simulation completed successfully."
cp -r src scripts/studio/run_dsm.sh postprocess/plotDSM.py postprocess/plotSSA.py postprocess/plotPorosity.py scripts/run_plotDSM.sh "$folder" || true

# Run SSA plotting on a single core
if [[ -f "$BASE_DIR/postprocess/plotSSA.py" ]]; then
  echo "[INFO] Running SSA plotting script..."
  if command -v srun >/dev/null 2>&1; then
    srun -N1 -n1 --cpus-per-task=1 python3 "$BASE_DIR/postprocess/plotSSA.py" "$folder" || true
  else
    python3 "$BASE_DIR/postprocess/plotSSA.py" "$folder" || true
  fi
else
  echo "[WARN] plotSSA.py not found at $BASE_DIR/postprocess/plotSSA.py; skipping."
fi

# echo
# echo "Running plotting scripts (if available)..."
# if [[ -x ./scripts/run_plotDSM.sh ]]; then
#   ./scripts/run_plotDSM.sh || true
# else
#   echo "[info] ./scripts/run_plotDSM.sh not found or not executable; skipping."
# fi

# echo "Simulation complete. Results stored in:"
# echo "$folder"