#!/bin/bash
#SBATCH -J DSM_run
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
#   Launch a DSM simulation on the cluster using PETSc .opts files.
#
# Required env vars (set before sbatch or export from batch script):
#   GEOMETRY_OPTS   — path to geometry .opts file
#   EXPERIMENT_OPTS — path to experiment .opts file
#
# Optional:
#   NUM_PROCS       — override MPI rank count (default: SLURM_NTASKS)
#   BASE_DIR        — project root (default: $PETIGA_DIR/projects/dry_snow_metamorphism)
#   output_dir      — parent directory for timestamped result folders
#
# Example sbatch invocation:
#   GEOMETRY_OPTS=inputs/grains__phi=0.24.../grains.opts \
#   EXPERIMENT_OPTS=inputs/experiment/30day_T-20_h1.00.opts \
#   sbatch scripts/HPC/run_dsm.sh
###############################################################################

set -euo pipefail

cleanup_on_err() {
  echo "[ERROR] Aborting due to an error (line $1)." >&2
  [[ -n "${SLURM_JOB_ID:-}" ]] && scancel "$SLURM_JOB_ID" >/dev/null 2>&1 || true
}
trap 'cleanup_on_err $LINENO' ERR

# ---- Paths ----
: "${PETIGA_DIR:=/Users/jacksonbaglino/PetIGA-3.20}"
BASE_DIR="${BASE_DIR:-${PETIGA_DIR}/projects/dry_snow_metamorphism}"
exec_file="${exec_file:-$BASE_DIR/dry_snow_metamorphism}"
output_dir="${output_dir:-/resnick/scratch/jbaglino}"

[[ -n "${GEOMETRY_OPTS:-}" ]]   || { echo "[ERROR] GEOMETRY_OPTS is not set." >&2; exit 1; }
[[ -n "${EXPERIMENT_OPTS:-}" ]] || { echo "[ERROR] EXPERIMENT_OPTS is not set." >&2; exit 1; }
[[ -f "$GEOMETRY_OPTS" ]]   || { echo "[ERROR] Geometry opts not found: $GEOMETRY_OPTS" >&2; exit 1; }
[[ -f "$EXPERIMENT_OPTS" ]] || { echo "[ERROR] Experiment opts not found: $EXPERIMENT_OPTS" >&2; exit 1; }

# ---- MPI ranks ----
if [[ -n "${SLURM_NTASKS:-}" ]]; then
  NUM_PROCS="$SLURM_NTASKS"
else
  NUM_PROCS="${NUM_PROCS:-40}"
fi

# ---- Extract parameters for title/metadata ----
get_opt() {
  local file="$1" flag="$2" default="${3:-}"
  grep -m1 "^-${flag} " "$file" 2>/dev/null | awk '{print $2}' || echo "$default"
}

Lx=$(get_opt "$GEOMETRY_OPTS" Lx "0")
Ly=$(get_opt "$GEOMETRY_OPTS" Ly "0")
dim=$(get_opt "$GEOMETRY_OPTS" dim "2")
grains_file=$(get_opt "$GEOMETRY_OPTS" grains_file "")
grain_dir="$(dirname "${grains_file:-$GEOMETRY_OPTS}")"

temp=$(get_opt "$EXPERIMENT_OPTS" temp "-10")
humidity=$(get_opt "$EXPERIMENT_OPTS" humidity "1.0")
t_final=$(get_opt "$EXPERIMENT_OPTS" t_final "2592000")

phi=""; seed=""
if [[ -f "$grain_dir/metadata.json" ]] && command -v jq >/dev/null 2>&1; then
  phi=$(jq -r '.structure.porosity_target // empty' "$grain_dir/metadata.json" 2>/dev/null || echo "")
  seed=$(jq -r '.generator.seed // empty' "$grain_dir/metadata.json" 2>/dev/null || echo "")
fi

ndays=$(awk "BEGIN{printf \"%d\", $t_final/86400}")
temp_tag=$(printf "Tm%.0f" "$temp")
hum_tag=$(awk "BEGIN{printf \"%02d\", $humidity*100}")

if [[ -n "$phi" && -n "$seed" ]]; then
  phi_tag=$(awk -v p="$phi" 'BEGIN{printf "phi%.2f", p}')
  Lx_tag=$(awk -v v="$Lx" 'BEGIN{printf "Lx%.0fmm", v*1000}')
  Ly_tag=$(awk -v v="$Ly" 'BEGIN{printf "Ly%.0fmm", v*1000}')
  title="DSM_${phi_tag}_${Lx_tag}_${Ly_tag}_seed${seed}_${temp_tag}_hum${hum_tag}_tf${ndays}d"
else
  geom_stem=$(basename "$(dirname "$GEOMETRY_OPTS")")
  title="DSM_${geom_stem}_${dim}D_${temp_tag}_hum${hum_tag}_tf${ndays}d"
fi

# Update SLURM job name
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  scontrol update JobId="$SLURM_JOB_ID" JobName="$title" 2>/dev/null || true
fi

# ---- Timestamped result folder ----
timestamp=$(date +%Y-%m-%d__%H.%M.%S)
folder="$output_dir/${title}_${timestamp}"
mkdir -p "$folder"
echo "[INFO] Output folder: $folder"

# ---- Build ----
cd "$BASE_DIR"
make dry_snow_metamorphism || { echo "[ERROR] Build failed." >&2; exit 1; }

# ---- Copy inputs (provenance) ----
[[ -n "$grains_file" && -f "$grains_file" ]] && cp "$grains_file" "$folder/"
cp "$GEOMETRY_OPTS" "$folder/"
cp "$EXPERIMENT_OPTS" "$folder/"
cp "$BASE_DIR/inputs/solver.opts" "$folder/"

# ---- Metadata ----
write_metadata_json() {
  if ! command -v jq >/dev/null 2>&1; then return 0; fi
  local grain_meta="$grain_dir/metadata.json"
  local dst_meta="$folder/metadata.json"
  [[ -f "$grain_meta" ]] && cp "$grain_meta" "$dst_meta"
  local augment
  augment=$(jq -n \
    --arg run_time "$(date -u +"%Y-%m-%dT%H:%M:%SZ")" \
    --arg host     "$(hostname)" \
    --arg geom_opts "$GEOMETRY_OPTS" \
    --arg exp_opts  "$EXPERIMENT_OPTS" \
    --arg output_dir "$folder" \
    --arg slurm_job "${SLURM_JOB_ID:-}" \
    --arg slurm_ntasks "${SLURM_NTASKS:-}" \
    '{ dsm_run: { run_time_utc: $run_time, executed_on: $host,
       opts_files: { geometry: $geom_opts, experiment: $exp_opts },
       output_dir: $output_dir,
       slurm: { job_id: $slurm_job, ntasks: $slurm_ntasks } } }')
  local tmp
  tmp="$(mktemp)"
  if [[ -f "$dst_meta" ]]; then
    jq --argjson add "$augment" '. * $add' "$dst_meta" > "$tmp" && mv "$tmp" "$dst_meta"
  else
    echo "$augment" > "$dst_meta"
  fi
}
write_metadata_json

# ---- Run simulation ----
echo "[INFO] Launching DSM simulation on $NUM_PROCS ranks..."
if command -v srun >/dev/null 2>&1; then
  set +e
  srun --kill-on-bad-exit=1 "$exec_file" \
    -options_file "$BASE_DIR/inputs/solver.opts" \
    -options_file "$GEOMETRY_OPTS" \
    -options_file "$EXPERIMENT_OPTS" \
    -output_dir "$folder" \
    2>&1 | tee "$folder/outp.txt"
  rc=${PIPESTATUS[0]}
  set -e
else
  set +e
  mpiexec -n "$NUM_PROCS" "$exec_file" \
    -options_file "$BASE_DIR/inputs/solver.opts" \
    -options_file "$GEOMETRY_OPTS" \
    -options_file "$EXPERIMENT_OPTS" \
    -output_dir "$folder" \
    2>&1 | tee "$folder/outp.txt"
  rc=${PIPESTATUS[0]}
  set -e
fi

if [[ ${rc:-0} -ne 0 ]]; then
  echo "[ERROR] Simulation failed with exit code $rc." >&2
  [[ -n "${SLURM_JOB_ID:-}" ]] && scancel "$SLURM_JOB_ID" >/dev/null 2>&1 || true
  exit "$rc"
fi

# ---- Finalize ----
echo "[INFO] Simulation completed. Archiving scripts..."
cp -r src "$folder/" || true
if [[ -f "$BASE_DIR/postprocess/plotSSA.py" ]]; then
  srun -N1 -n1 --cpus-per-task=1 python3 "$BASE_DIR/postprocess/plotSSA.py" "$folder" || true
fi
echo "[INFO] Results stored in: $folder"
