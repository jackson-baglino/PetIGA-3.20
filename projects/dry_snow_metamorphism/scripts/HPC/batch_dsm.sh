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

# Default INPUT_DIR (can be overridden by env var, INPUT_DIR_OVERRIDE, or CLI args)
DEFAULT_INPUT_DIR="$PETIGA_DIR/projects/dry_snow_metamorphism/inputs"
INPUT_DIR="${INPUT_DIR:-$DEFAULT_INPUT_DIR}"

# --- Lightweight CLI parsing / overrides ---
# Usage examples:
#   sbatch batch_dsm.sh --input-dir=/path/to/inputs
#   sbatch batch_dsm.sh /path/to/inputs
#   DSM_TEST=1 sbatch batch_dsm.sh
FORCE_TEST=0
# Allow INPUT_DIR override via arg or env var INPUT_DIR_OVERRIDE
if [[ -n "${INPUT_DIR_OVERRIDE:-}" ]]; then
  INPUT_DIR="$INPUT_DIR_OVERRIDE"
fi
# Parse args
for arg in "$@"; do
  case "$arg" in
    --test) FORCE_TEST=1 ;;
    --input-dir=*) INPUT_DIR="${arg#*=}" ;;
    -i) shift; INPUT_DIR="${1:-$INPUT_DIR}" ;;
    --) shift; break ;;
    -h|--help)
      echo "Usage: sbatch batch_dsm.sh [--test] [--input-dir=/path] [INPUT_DIR]"
      exit 0
      ;;
    *)
      # If it's a path and exists, treat as INPUT_DIR (positional)
      if [[ -d "$arg" ]]; then
        INPUT_DIR="$arg"
      fi
      ;;
  esac
done
if [[ -z "${INPUT_DIR:-}" || ! -d "$INPUT_DIR" ]]; then
  echo "[ERROR] INPUT_DIR not set or does not exist: '${INPUT_DIR:-unset}'"
  exit 1
fi

# --- Sweep settings ---
temperatures=(-24)
humidities=(0.98)

# --- Paths ---
RUN_SCRIPT="/resnick/groups/rubyfu/jbaglino/PetIGA-3.20/projects/dry_snow_metamorphism/scripts/HPC/run_dsm.sh"
CONFIG_DIR="/resnick/groups/rubyfu/jbaglino/PetIGA-3.20/projects/dry_snow_metamorphism/inputs"

# --- Resource tuning knobs (override by exporting before running) ---
: "${MIN_TASKS_PER_NODE:=28}"       # lower bound for cores/node
: "${MAX_TASKS_PER_NODE:=32}"       # upper bound for cores/node
: "${TASKS_PER_NODE_DEFAULT:=32}"   # preferred packing per node (will be clamped to [MIN,MAX])
: "${MEM_PER_CPU_DEFAULT:=1G}"      # memory request per CPU

shopt -s nullglob globstar
# Recursively find any grains.dat under INPUT_DIR
dat_files=("$INPUT_DIR"/grains__phi=*__Lxmm=2__Lymm=2__seed=*/grains.dat)
shopt -u globstar nullglob

if [ ${#dat_files[@]} -eq 0 ]; then
  echo "[ERROR] No grains.dat files found under INPUT_DIR (searched recursively): $INPUT_DIR"
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
      v=3*Nx*Ny/14000.0; 
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

  # --- Optional test mode: force 1 node and 10 cores total ---
  if [[ "${DSM_TEST:-0}" == "1" || "$FORCE_TEST" -eq 1 ]]; then
    NODES=1
    TASKS_PER_NODE=5
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
  if [[ "${DSM_TEST:-0}" == "1" || "$FORCE_TEST" -eq 1 ]]; then
    echo " Mode:     TEST (forcing 1 node, 10 cores)"
  fi
  echo
  echo " Parsed from .env:"
  echo "   Nx = ${Nx:-unset}"
  echo "   Ny = ${Ny:-unset}"
  echo "   Nz = ${Nz:-unset}"
  echo
  echo " Resource plan (targets ~14k pts/core, 35–50 cores/node):"
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

      # Build descriptive job name using physical parameters from the directory name
      dirbase="$(basename "$env_dir")"
      phi_tag="$(echo "$dirbase" | sed -n 's/.*phi=\([0-9.]*\).*/\1/p')"
      Lx_tag="$(echo "$dirbase" | sed -n 's/.*Lxmm=\([0-9.]*\).*/\1/p')"
      Ly_tag="$(echo "$dirbase" | sed -n 's/.*Lymm=\([0-9.]*\).*/\1/p')"
      seed_tag="$(echo "$dirbase" | sed -n 's/.*seed=\([0-9]*\).*/\1/p')"
      # Fallbacks if parsing fails
      [[ -z "$phi_tag"  ]] && phi_tag="NA"
      if [[ -z "$Lx_tag" && -n "${Lx:-}" ]]; then Lx_tag="$(awk -v v="$Lx" 'BEGIN{printf "%.3g", v*1000}')"; fi
      if [[ -z "$Ly_tag" && -n "${Ly:-}" ]]; then Ly_tag="$(awk -v v="$Ly" 'BEGIN{printf "%.3g", v*1000}')"; fi
      [[ -z "$Lx_tag"  ]] && Lx_tag="NA"
      [[ -z "$Ly_tag"  ]] && Ly_tag="NA"
      [[ -z "$seed_tag" ]] && seed_tag="NA"
      job_base="DSM_phi${phi_tag}_Lx${Lx_tag}_Ly${Ly_tag}_seed${seed_tag}_Tm${temp_tag}_hum${hum_tag}"

      # Define okay architectures if not already set
      : "${ARCH_OK:='icelake|skylake|cascadelake&!cascadelake_lowmem'}"

      echo "Submitting job for file: $basename, T=${temp_tag}C, RH=${hum_tag}%"
      sbatch --job-name="$job_base" \
            --nodes="$NODES" \
            --ntasks-per-node="$TASKS_PER_NODE" \
            --cpus-per-task="$CPUS_PER_TASK" \
            --mem-per-cpu="$MEM_PER_CPU" \
            --export=ALL,inputFile="$INPUT_FILE",temp="$temp",humidity="$hum" \
            "$RUN_SCRIPT"
    done
  done
done