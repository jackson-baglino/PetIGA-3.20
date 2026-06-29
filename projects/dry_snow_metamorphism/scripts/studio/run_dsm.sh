#!/bin/zsh
###############################################################################
# Script: run_dsm.sh (studio)
# Model: Dry Snow Metamorphism (DSM)
# Purpose:
#   Configure and launch a single DSM simulation, manage run metadata, and
#   archive inputs/plots to a timestamped results folder.
#
# Parameters are passed via PETSc .opts files, NOT environment variables.
#
# Typical Usage:
#   ./scripts/studio/run_dsm.sh \
#       -g inputs/grains__phi=0.24__Lxmm=2__Lymm=2__seed=7/grains.opts \
#       -e inputs/experiment/30day_T-10_h1.00.opts
#
# Options:
#   -g, --geometry PATH   Path to geometry .opts file (required)
#   -e, --experiment PATH Path to experiment .opts file (required)
#   -p, --procs N         MPI ranks (default: 12)
#   --no-build            Skip rebuilding the executable
#   --dry-run             Print actions without running
#   -?, --help            Show this help and exit
###############################################################################

set -euo pipefail

log()  { printf "[INFO] %s\n" "$*"; }
warn() { printf "[WARN] %s\n" "$*"; }
err()  { printf "[ERROR] %s\n" "$*" 1>&2; }
die()  { err "$*"; exit 1; }
require_cmd() { command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"; }

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Options:
  -g, --geometry PATH   Path to geometry .opts file (required if not using GEOMETRY_OPTS env var)
  -e, --experiment PATH Path to experiment .opts file (required if not using EXPERIMENT_OPTS env var)
  -p, --procs N         MPI ranks (default: 12)
  --no-build            Skip rebuilding the executable
  --dry-run             Print actions without running the simulation
  -?, --help            Show this help and exit

Examples:
  # Standard run
  ./scripts/studio/run_dsm.sh \\
      -g inputs/grains__phi=0.24__Lxmm=2__Lymm=2__seed=7/grains.opts \\
      -e inputs/experiment/30day_T-10_h1.00.opts

  # Quick override via CLI (PETSc merges all sources, later overrides earlier)
  ./scripts/studio/run_dsm.sh -g ... -e ... -- -temp -5.0 -n_out 50
EOF
}

# ---- Argument parsing ----
NO_BUILD=0
DRY_RUN=0
GEOMETRY_OPTS="${GEOMETRY_OPTS:-}"
EXPERIMENT_OPTS="${EXPERIMENT_OPTS:-}"
NUM_PROCS="${NUM_PROCS:-12}"
EXTRA_OPTS=()

if [[ ${#@} -gt 0 ]]; then
  while [[ ${#} -gt 0 ]]; do
    case "${1}" in
      -g|--geometry)   GEOMETRY_OPTS="${2:-}"; shift 2;;
      -e|--experiment) EXPERIMENT_OPTS="${2:-}"; shift 2;;
      -p|--procs)      NUM_PROCS="${2:-}"; shift 2;;
      --no-build)      NO_BUILD=1; shift;;
      --dry-run)       DRY_RUN=1; shift;;
      -\?|--help|-help) usage; exit 0;;
      --) shift; EXTRA_OPTS=("$@"); break;;
      *) warn "Ignoring unknown arg: ${1}"; shift;;
    esac
  done
fi

[[ -n "$GEOMETRY_OPTS" ]]   || die "Geometry .opts file required (-g / --geometry or GEOMETRY_OPTS env var)."
[[ -n "$EXPERIMENT_OPTS" ]] || die "Experiment .opts file required (-e / --experiment or EXPERIMENT_OPTS env var)."
[[ -f "$GEOMETRY_OPTS" ]]   || die "Geometry opts file not found: $GEOMETRY_OPTS"
[[ -f "$EXPERIMENT_OPTS" ]] || die "Experiment opts file not found: $EXPERIMENT_OPTS"

# ---- Paths ----
BASE_DIR="${BASE_DIR:-${PETIGA_DIR}/projects/dry_snow_metamorphism}"
exec_file="${exec_file:-$BASE_DIR/dry_snow_metamorphism}"
output_dir="${output_dir:-/Users/jacksonbaglino/SimulationResults/dry_snow_metamorphism/scratch}"

# ---- Extract parameters from opts files for title/metadata ----
# Helper: extract the value of a -flag from an opts file (first match wins)
get_opt() {
  local file="$1" flag="$2" default="${3:-}"
  grep -m1 "^-${flag} " "$file" 2>/dev/null | awk '{print $2}' || echo "$default"
}

Lx=$(get_opt "$GEOMETRY_OPTS" Lx "0")
Ly=$(get_opt "$GEOMETRY_OPTS" Ly "0")
dim=$(get_opt "$GEOMETRY_OPTS" dim "2")
grains_file=$(get_opt "$GEOMETRY_OPTS" grains_file "")

temp=$(get_opt "$EXPERIMENT_OPTS" temp "-10")
humidity=$(get_opt "$EXPERIMENT_OPTS" humidity "1.0")
t_final=$(get_opt "$EXPERIMENT_OPTS" t_final "2592000")
n_out=$(get_opt "$EXPERIMENT_OPTS" n_out "500")
delt_t=$(get_opt "$EXPERIMENT_OPTS" delt_t "1e-4")

# Try to read phi and seed from metadata.json adjacent to the grains file
grain_dir="$(dirname "${grains_file:-$GEOMETRY_OPTS}")"
phi=""; seed=""
if [[ -f "$grain_dir/metadata.json" ]] && command -v jq >/dev/null 2>&1; then
  phi=$(jq -r '.structure.porosity_target // empty' "$grain_dir/metadata.json" 2>/dev/null || echo "")
  seed=$(jq -r '.generator.seed // empty' "$grain_dir/metadata.json" 2>/dev/null || echo "")
fi

# ---- Build run title ----
ndays=$(awk "BEGIN{printf \"%d\", $t_final/86400}")
temp_tag=$(printf "Tm%.0f" "$temp")
hum_tag=$(awk "BEGIN{printf \"%02d\", $humidity*100}")
exp_stem=$(basename "${EXPERIMENT_OPTS%.opts}")
geom_stem=$(basename "$(dirname "$GEOMETRY_OPTS")")

if [[ -n "$phi" && -n "$seed" ]]; then
  phi_tag=$(awk -v p="$phi" 'BEGIN{printf "phi%.2f", p}')
  Lx_tag=$(awk -v v="$Lx" 'BEGIN{printf "Lx%.0fmm", v*1000}')
  Ly_tag=$(awk -v v="$Ly" 'BEGIN{printf "Ly%.0fmm", v*1000}')
  title="DSM_${phi_tag}_${Lx_tag}_${Ly_tag}_seed${seed}_${temp_tag}_hum${hum_tag}_tf${ndays}d"
else
  title="DSM_${geom_stem}_${dim}D_${temp_tag}_hum${hum_tag}_tf${ndays}d"
fi

# ---- Timestamped result folder ----
timestamp=$(date +%Y-%m-%d__%H.%M.%S)
folder="$output_dir/${title}_${timestamp}"
[[ "$DRY_RUN" -eq 1 ]] || mkdir -p "$folder"

log "Run title: $title"
log "Output folder: $folder"
log "Geometry opts: $GEOMETRY_OPTS"
log "Experiment opts: $EXPERIMENT_OPTS"

# ---- Build ----
cd "$BASE_DIR" || exit 1
if [[ "$NO_BUILD" -eq 0 ]]; then
  require_cmd make
  log "Building executable..."
  if [[ "$DRY_RUN" -eq 1 ]]; then
    log "(dry-run) make dry_snow_metamorphism"
  else
    make dry_snow_metamorphism || die "Build failed. Please check the Makefile and dependencies."
  fi
else
  warn "Skipping build (--no-build)."
fi

# ---- Ensure grains.opts exists alongside grains.dat ----
if [[ -n "$grains_file" && -f "$grains_file" ]]; then
  grain_opts="$grain_dir/grains.opts"
  if [[ ! -f "$grain_opts" ]]; then
    log "grains.opts not found; generating via preprocess/generate_opts_from_input.py ..."
    if [[ "$DRY_RUN" -eq 1 ]]; then
      log "(dry-run) python3 preprocess/generate_opts_from_input.py \"$grains_file\" \"$grain_opts\""
    else
      python3 preprocess/generate_opts_from_input.py "$grains_file" "$grain_opts" || die "Failed to generate $grain_opts"
    fi
  fi
fi

# ---- Copy input files to results folder (provenance) ----
if [[ "$DRY_RUN" -eq 0 ]]; then
  [[ -n "$grains_file" && -f "$grains_file" ]] && cp "$grains_file" "$folder/"
  cp "$GEOMETRY_OPTS" "$folder/"
  cp "$EXPERIMENT_OPTS" "$folder/"
  cp "$BASE_DIR/inputs/solver.opts" "$folder/"
fi

# ---- Write resolved params snapshot ----
write_opts_snapshot() {
  local snapshot="$folder/resolved_params.opts"
  {
    echo "# Resolved opts snapshot for this run"
    echo "# Generated: $(date -u +"%Y-%m-%dT%H:%M:%SZ")"
    echo "# solver.opts:    $BASE_DIR/inputs/solver.opts"
    echo "# geometry opts:  $GEOMETRY_OPTS"
    echo "# experiment opts: $EXPERIMENT_OPTS"
    echo "# output_dir:     $folder"
    echo ""
    echo "-output_dir $folder"
  } > "$snapshot"
}

# ---- Write metadata ----
write_metadata_json() {
  if ! command -v jq >/dev/null 2>&1; then return 0; fi
  local grain_meta="$grain_dir/metadata.json"
  local dst_meta="$folder/metadata.json"
  if [[ -f "$grain_meta" ]]; then cp "$grain_meta" "$dst_meta"; fi
  local augment
  augment=$(jq -n \
    --arg run_time "$(date -u +"%Y-%m-%dT%H:%M:%SZ")" \
    --arg host     "$(hostname)" \
    --arg user     "$(whoami)" \
    --arg geom_opts "$GEOMETRY_OPTS" \
    --arg exp_opts  "$EXPERIMENT_OPTS" \
    --arg solver_opts "$BASE_DIR/inputs/solver.opts" \
    --arg output_dir "$folder" \
    '{ dsm_run: { run_time_utc: $run_time, executed_on: $host, user: $user,
       opts_files: { solver: $solver_opts, geometry: $geom_opts, experiment: $exp_opts },
       output_dir: $output_dir } }')
  local tmp
  tmp="$(mktemp)"
  if [[ -f "$dst_meta" ]]; then
    jq --argjson add "$augment" '. * $add' "$dst_meta" > "$tmp" && mv "$tmp" "$dst_meta"
  else
    echo "$augment" > "$dst_meta"
  fi
  log "Metadata written: $dst_meta"
}

if [[ "$DRY_RUN" -eq 0 ]]; then
  write_opts_snapshot
  write_metadata_json
fi

# ---- Run the simulation ----
log "Launching DSM simulation..."
require_cmd mpiexec

SIM_CMD=(
  mpiexec -np "$NUM_PROCS" "$exec_file"
  -options_file "$BASE_DIR/inputs/solver.opts"
  -options_file "$GEOMETRY_OPTS"
  -options_file "$EXPERIMENT_OPTS"
  -output_dir "$folder"
  "${EXTRA_OPTS[@]}"
)

if [[ "$DRY_RUN" -eq 1 ]]; then
  log "(dry-run) ${SIM_CMD[*]} | tee $folder/outp.txt"
else
  "${SIM_CMD[@]}" | tee "$folder/outp.txt"
fi

# ---- Finalize ----
log "Simulation completed. Archiving source and scripts..."
if [[ "$DRY_RUN" -eq 0 ]]; then
  cp -r src postprocess/plotDSM.py postprocess/plotSSA.py postprocess/plotPorosity.py "$folder/" || true
fi

log "Results stored in: $folder"
