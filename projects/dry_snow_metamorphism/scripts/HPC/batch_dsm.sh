#!/usr/bin/env bash
###############################################################################
# Script: batch_dsm.sh (HPC)
# Model: Dry Snow Metamorphism (DSM)
# Purpose:
#   Submit a sweep of DSM simulations to SLURM, varying temperature(s),
#   humidity(ies), and input file(s). Calls run_dsm.sh for each combination.
#
# Usage:
#   sbatch batch_dsm.sh
#
# Notes:
#   - Overrides temp, humidity, and inputFile by exporting them to SLURM.
#   - RUN_LABEL is generated inside run_dsm.sh using compact tags (2-digit temp/RH).
#   - Adjust INPUT_DIR / CONFIG_DIR, temperatures, humidities below as needed.
###############################################################################

set -euo pipefail

# --- Sweep settings ---
temperatures=(-25)
humidities=(0.98)

# --- Paths ---
RUN_SCRIPT="/resnick/groups/rubyfu/jbaglino/PetIGA-3.20/projects/dry_snow_metamorphism/scripts/HPC/run_dsm.sh"
INPUT_DIR="/resnick/groups/rubyfu/jbaglino/PetIGA-3.20/projects/dry_snow_metamorphism/inputs"
CONFIG_DIR="/resnick/groups/rubyfu/jbaglino/PetIGA-3.20/projects/dry_snow_metamorphism/inputs"

# --- Resource tuning knobs (override by exporting before running) ---
: "${MIN_TASKS_PER_NODE:=35}"       # lower bound for cores/node
: "${MAX_TASKS_PER_NODE:=50}"       # upper bound for cores/node
: "${TASKS_PER_NODE_DEFAULT:=40}"   # preferred packing per node (will be clamped to [MIN,MAX])
: "${MEM_PER_CPU_DEFAULT:=1G}"      # memory request per CPU

shopt -s nullglob
dat_files=("$INPUT_DIR"/grains__phi=*__Lxmm=2__Lymm=2__seed=*/grains.dat)
shopt -u nullglob

if [ ${#dat_files[@]} -eq 0 ]; then
  echo "No .dat files found in $INPUT_DIR"
  exit 1
fi

make clean && make all
for dat_file in "${dat_files[@]}"; do
  basename=$(basename "$dat_file" .dat)
  env_dir="$(dirname "$dat_file")"

  # Prefer the generator’s default name; fall back to basename.env if you ever use that
  if   [ -f "$env_dir/grains.env" ]; then
    env_file="$env_dir/grains.env"
  elif [ -f "$env_dir/$basename.env" ]; then
    env_file="$env_dir/$basename.env"
  else
    echo "Warning: No .env found next to $dat_file (looked for grains.env or $basename.env). Skipping."
    continue
  fi

  ENV_FILE="$env_file"
  INPUT_FILE="$dat_file"

  # --- Derive per-job resources from the .env ---
  set -a; source "$ENV_FILE"; set +a

  # Clamp preferred tasks-per-node into requested band [MIN, MAX]
  pref_tpn="$TASKS_PER_NODE_DEFAULT"
  if (( pref_tpn < MIN_TASKS_PER_NODE )); then pref_tpn="$MIN_TASKS_PER_NODE"; fi
  if (( pref_tpn > MAX_TASKS_PER_NODE )); then pref_tpn="$MAX_TASKS_PER_NODE"; fi

  if [[ -z "${Nx:-}" || -z "${Ny:-}" ]]; then
    echo "Warning: Nx/Ny not found in $ENV_FILE; using defaults (NODES=1, TASKS_PER_NODE=$pref_tpn)"
    NODES=1
    TASKS_PER_NODE="$pref_tpn"
  else
    # Target total cores so that elements/core ~ 10,000 with 3 fields:
    #   target_cores = ceil( 3*Nx*Ny / 10000 )
    target_cores=$(awk -v Nx="$Nx" -v Ny="$Ny" 'BEGIN{
      v=3*Nx*Ny/10000.0; 
      if (v<1) v=1;
      print (v==int(v)?v:int(v)+1)
    }')

    # First pass: choose nodes based on preferred tasks-per-node (clamped)
    NODES=$(awk -v t="$target_cores" -v p="$pref_tpn" 'BEGIN{
      n=int((t + p - 1)/p); if(n<1)n=1; print n
    }')
    TASKS_PER_NODE=$(awk -v t="$target_cores" -v n="$NODES" 'BEGIN{
      p=int((t + n - 1)/n); if(p<1)p=1; print p
    }')

    # If over max, increase nodes until within MAX_TASKS_PER_NODE
    if (( TASKS_PER_NODE > MAX_TASKS_PER_NODE )); then
      NODES=$(awk -v t="$target_cores" -v p="$MAX_TASKS_PER_NODE" 'BEGIN{
        n=int((t + p - 1)/p); if(n<1)n=1; print n
      }')
      TASKS_PER_NODE=$(awk -v t="$target_cores" -v n="$NODES" 'BEGIN{
        p=int((t + n - 1)/n); if(p<1)p=1; print p
      }')
    fi

    # If under min (and we have room), try reducing nodes to raise TPN into band
    if (( TASKS_PER_NODE < MIN_TASKS_PER_NODE )); then
      # Try to find the largest n that keeps p in [MIN,MAX]
      NODES=$(awk -v t="$target_cores" -v minp="$MIN_TASKS_PER_NODE" -v maxp="$MAX_TASKS_PER_NODE" 'BEGIN{
        # start from ceil(t/maxp) (min nodes to stay <= maxp) and try decreasing
        n=int((t + maxp - 1)/maxp); if(n<1)n=1;
        for (; n>=1; n--) {
          p=int((t + n - 1)/n);
          if (p>=minp && p<=maxp) { print n; exit }
        }
        print 1
      }')
      TASKS_PER_NODE=$(awk -v t="$target_cores" -v n="$NODES" 'BEGIN{
        p=int((t + n - 1)/n); if(p<1)p=1; print p
      }')
      # Final clamp to bounds, if still outside we accept nearest feasible
      if (( TASKS_PER_NODE > MAX_TASKS_PER_NODE )); then TASKS_PER_NODE="$MAX_TASKS_PER_NODE"; fi
      if (( TASKS_PER_NODE < MIN_TASKS_PER_NODE )); then TASKS_PER_NODE="$MIN_TASKS_PER_NODE"; fi
    fi
  fi

  CPUS_PER_TASK=1
  MEM_PER_CPU="$MEM_PER_CPU_DEFAULT"

  # For info only (approximate): GRID if Nz provided
  if [[ -n "${Nz:-}" ]]; then
    GRID=$(awk -v Nx="$Nx" -v Ny="$Ny" -v Nz="$Nz" 'BEGIN{print Nx*Ny*Nz}')
  else
    GRID="unknown"
  fi

  per_core=$(awk -v Nx="${Nx:-0}" -v Ny="${Ny:-0}" -v N="$NODES" -v P="$TASKS_PER_NODE" \
             'BEGIN{den=N*P; if(den<=0){print "inf"} else {printf "%.0f", (3*Nx*Ny)/den}}')

  # ---- Pretty diagnostics block -------------------------------------------
  echo "------------------------------------------------------------"
  echo " Input:    $INPUT_FILE"
  echo " Config:   $ENV_FILE"
  echo " Basename: $basename"
  echo
  echo " Parsed from .env:"
  echo "   Nx = ${Nx:-unset}"
  echo "   Ny = ${Ny:-unset}"
  echo "   Nz = ${Nz:-unset}"
  echo
  echo " Resource plan (targets ~10k pts/core, 35–50 cores/node):"
  echo "   target_cores = ${target_cores:-unknown}"
  echo "   NODES        = ${NODES}"
  echo "   TASKS/Node   = ${TASKS_PER_NODE}  (band: ${MIN_TASKS_PER_NODE}-${MAX_TASKS_PER_NODE})"
  echo "   CPUS/Task    = ${CPUS_PER_TASK}"
  echo "   MEM/CPU      = ${MEM_PER_CPU}"
  echo "   GRID         = ${GRID}"
  echo "   Points/Core  ≈ ${per_core}"
  echo "------------------------------------------------------------"
  # -------------------------------------------------------------------------

  for temp in "${temperatures[@]}"; do
    for hum in "${humidities[@]}"; do
      temp_int=$(printf "%.0f" "$temp")
      if [[ "$temp_int" == -* ]]; then temp_tag=${temp_int:0:3}; else temp_tag=${temp_int:0:2}; fi
      hum_int=$(awk "BEGIN{printf \"%d\", $hum*100}")
      hum_tag=$(printf "%02d" "$hum_int"); hum_tag=${hum_tag:0:2}

      # Build descriptive job name and log destinations
      nx_tag=${Nx:-NA}
      ny_tag=${Ny:-NA}
      nz_tag=${Nz:-}
      if [[ -n "$nz_tag" ]]; then
        grid_tag="${nx_tag}x${ny_tag}x${nz_tag}"
      else
        grid_tag="${nx_tag}x${ny_tag}"
      fi
      job_base="DSM-${basename}_Tm${temp_tag}_hum${hum_tag}_grid${grid_tag}"
      log_dir="$env_dir/slurm"
      mkdir -p "$log_dir"

      echo "Submitting job for file: $basename, T=${temp_tag}C, RH=${hum_tag}%"
      sbatch --job-name="$job_base" \
             --output="$log_dir/%x__%j.out" \
             --error="$log_dir/%x__%j.err" \
             --nodes="$NODES" \
             --ntasks-per-node="$TASKS_PER_NODE" \
             --cpus-per-task="$CPUS_PER_TASK" \
             --mem-per-cpu="$MEM_PER_CPU" \
             --export=ALL,ENV_FILE_OVERRIDE="$ENV_FILE",temp="$temp",humidity="$hum",inputFile="$INPUT_FILE",filename="$basename" \
             "$RUN_SCRIPT" \
             --kill-on-bad-exit=1
    done
  done
done