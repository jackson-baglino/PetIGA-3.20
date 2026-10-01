#!/usr/bin/env bash
# =============================================================================
# submit_keff_production.sh — THE submit command for the k_eff sintering
# manuscript runs. Every production run goes through this script, so every run
# gets exactly the same options.
#
#   ./scripts/HPC/submit_keff_production.sh studies/keff_sintering/batch3a_shakedown.txt
#   ./scripts/HPC/submit_keff_production.sh <stage file> --dry-run
#
# The stage file lists RUNS ONLY (geometry:experiment, one per line). All
# options live here, in PRODUCTION_OPTS below. The script refuses a stage file
# that carries its own options (a third field), a geometry that is not one of
# the campaign packings, or an experiment that is not the campaign's 30-day
# snow_T*_h1.00 family -- so a run cannot quietly differ from the others.
#
# k_eff CADENCE: THE SAME OPTION AT EVERY TEMPERATURE (since 2026-09-29).
#   -keff_dlnssa 0.001: sample every 0.1% drop in SSA, from t = 11 tau_sub on
#   (-keff_dlnssa_t0_tau, the baseline); every 5 steps before that (the IC
#   relaxation); never more than 20 tau_sub between samples.
# A step count spent samples evenly in steps, leaving visible corners where
# k_eff bends fastest and wasting samples on the straight late curve. The SSA
# trigger puts them along the k-SSA curve itself, and because temperature acts
# as a time rescaling it samples every temperature at the same states -- so the
# old "-keff_freq 1 at -40 C" rule is gone. On batch 2's every-step reference,
# 0.1% gives 0.02% max interpolation error, 60x below the smallest real kink
# seen (the SSA~15300 event on 3a seed 301, >= 1.3% of the plotted range).
# Samples per run ~308 / 179 / 46 at -5 / -20 / -40 C.
#
# TIME LIMITS PER RUN (2026-10-01). Each job asks for a limit sized to its
# temperature, not a flat 24 h: the scheduler backfills a job into a gap only
# if its limit fits, and with the account over its fair share (sshare,
# 2026-10-01) backfill is how jobs start. Base limits at 200k DoF/core are
# ~2x the wall time predicted from the scaling test (23.4 s/step, 5.2 s per
# k_eff sample at 121 ranks; steps ~ t_final/dtmax + ~70):
#     -5 C 18 h (~8.5 h)   -10 C 12 h (~5.7 h)   -20 C 6 h (~2.7 h)
#     -30/-40 C 4 h (< 1.5 h)
# A larger mesh is scaled by its DoF per rank relative to the target (the
# L/R 56/80 runs are capped at MAX_NODES_PER_JOB, so they carry more), and
# nothing asks for more than 24 h. See time_limit_for() below.
#
# GUARDS. The repo must be committed and pushed: the batch records the commit
# (print_repo_provenance), and a run from an uncommitted tree cannot be
# reproduced. The allocation target is whatever scripts/lib/alloc.sh says
# (60k DoF/core since 2026-09-26) -- printed, not overridden.
#
# RECORD. A PRODUCTION_MANIFEST.txt is written into the batch folder: commit,
# options, cadence rule, and the exact run list. Each run's SLURM .o also
# carries "Extra opts : ..." (run_enceladus.sh), copied into its run folder.
#
# Output: ONE campaign folder for every manuscript run (2026-10-01):
#   /resnick/groups/rubyfu/jbaglino/simulation_outputs/enceladus_DSM/keff_sintering_campaign/
#       <geom>__<exp>/                    one folder per run, names unique
#       stages/<stage>__<timestamp>/      PRODUCTION_MANIFEST.txt, the stage
#                                         file, job ids, and submit_batch's
#                                         inputs/src snapshot for that stage
# A run folder that already exists is refused (submit_batch --parent-dir), so
# resubmitting a stage cannot write over finished results. Download a stage
# with scripts/HPC/fetch_stage.sh <stage file> (one rsync, one 2FA prompt).
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# ---------------------------------------------------------------------------
# THE PRODUCTION OPTIONS. Change these only deliberately, and only between
# stages -- a change mid-campaign splits the manuscript's runs into two sets.
# ---------------------------------------------------------------------------
PRODUCTION_OPTS=(
    -keff 1                    # in-line k_eff (tensor law, the default since 2026-09-23)
    -keff_step0 1              # sample the initial condition too
    -keff_dlnssa 0.001         # sample every 0.1% drop in SSA ...
    -keff_dlnssa_t0_tau 11.05  #   ... from t = 11 tau_sub (1 d at -20 C) on
    -keff_freq 5               #   every 5 steps before that (IC relaxation)
    -keff_max_gap_tau 20       #   and never more than 20 tau_sub apart
    -keff_ksp_type cg          # corrector solve: CG
    -keff_pc_type gamg         #   + algebraic multigrid
    -t_out_log 50              # 50 log-spaced field snapshots (~9 GB/run)
    -t_out_log_t0 60           #   starting at 60 s
)                              # (the t >= 1 s opening frame is the solver default)
OUT_ROOT="/resnick/groups/rubyfu/jbaglino/simulation_outputs"
CAMPAIGN_DIR="$OUT_ROOT/enceladus_DSM/keff_sintering_campaign"
# Campaign packing families: the production matrix and the domain-size
# convergence study (supplement), both one generator recipe (seam +
# percolation gates only, option B 2026-10-01); and the pre-B GATED build,
# run only in batch_rve.txt to measure what its homogeneity gates did.
PACKING_FAMILIES=("inputs/packings/keff_LR40/" "inputs/packings/rve_phi0.325/"
                  "inputs/packings/keff_LR40_gated/")
EXP_PATTERN='^snow_T-?[0-9]+_h1\.00_30d$'

usage() { sed -n '2,57p' "$0"; exit "${1:-0}"; }

# time_limit_for <T in C> <geometry opts file>  ->  "HH:MM:SS"
time_limit_for() {
    local T="$1" g="$2"
    local base
    base=$(awk -v t="$T" 'BEGIN{ if (t >= -7.5) print 18; else if (t >= -15) print 12;
                                 else if (t >= -25) print 6; else print 4 }')
    local nx ny dof
    nx=$(awk '$1=="-Nx"{print $2; exit}' "$g"); ny=$(awk '$1=="-Ny"{print $2; exit}' "$g")
    dof=$(awk '$1=="-dof"{print $2; exit}' "$PROJECT_ROOT/inputs/solver.opts"); dof=${dof:-3}
    awk -v nx="$nx" -v ny="$ny" -v d="$dof" -v b="$base" -v tgt="$TARGET_DOFS_PER_CORE" \
        -v cap="$((MAX_NODES_PER_JOB * MAX_TASKS_PER_NODE))" 'BEGIN{
        N = d * nx * ny; P = int((N + tgt - 1) / tgt); if (cap > 0 && P > cap) P = cap
        f = (N / P) / tgt; if (f < 1) f = 1
        h = int(b * f); if (h < b * f) h++; if (h > 24) h = 24
        printf "%02d:00:00\n", h }'
}
source "$PROJECT_ROOT/scripts/lib/alloc.sh"

stage="" ; dry=0
for arg in "$@"; do
    case "$arg" in
        --dry-run) dry=1 ;;
        -h|--help) usage 0 ;;
        -*) echo "unknown option: $arg" >&2; usage 1 ;;
        *) stage="$arg" ;;
    esac
done
[[ -n "$stage" && -f "$stage" ]] || { echo "❌ stage file not found: ${stage:-<none>}" >&2; usage 1; }
cd "$PROJECT_ROOT"
stage_name="$(basename "$stage" .txt)"
tag="keff_${stage_name#batch}"

# ---------------------------------------------------------------------------
# Validate the stage file and build the per-run specs
# ---------------------------------------------------------------------------
specs=() ; errors=0 ; n=0
while IFS= read -r line; do
    line="${line%%#*}"
    line="${line#"${line%%[![:space:]]*}"}"; line="${line%"${line##*[![:space:]]}"}"
    [[ -z "$line" ]] && continue
    n=$((n + 1))
    IFS=':' read -r geom exp extra <<< "$line"
    if [[ -n "${extra:-}" ]]; then
        echo "❌ line $n carries its own options ('$extra'): options belong in this script" >&2
        errors=$((errors + 1)); continue
    fi
    gfile=$(find inputs/geometry -name "${geom}.opts" -print -quit)
    efile=$(find inputs/experiment -name "${exp}.opts" -print -quit)
    if [[ -z "$gfile" ]]; then echo "❌ geometry not found: $geom" >&2; errors=$((errors+1)); continue; fi
    if [[ -z "$efile" ]]; then echo "❌ experiment not found: $exp" >&2; errors=$((errors+1)); continue; fi
    fam_ok=0
    for fam in "${PACKING_FAMILIES[@]}"; do
        grep -q "^-grains_file ${fam}" "$gfile" && fam_ok=1
    done
    if (( ! fam_ok )); then
        echo "❌ $geom is not a campaign packing (grains_file not under ${PACKING_FAMILIES[*]})" >&2
        errors=$((errors + 1)); continue
    fi
    if ! [[ "$exp" =~ $EXP_PATTERN ]]; then
        echo "❌ $exp is not a campaign experiment (want snow_T<T>_h1.00_30d)" >&2
        errors=$((errors + 1)); continue
    fi
    T=$(awk '$1=="-temp"{print $2; exit}' "$efile")
    TG=$(awk '$1=="-eps_valid_temp"{print $2; exit}' "$gfile")
    if [[ -z "$T" || -z "$TG" ]] || awk -v a="$T" -v b="$TG" 'BEGIN{exit !((a-b)>1 || (b-a)>1)}'; then
        echo "❌ $geom (eps_valid_temp $TG) does not match $exp (temp $T)" >&2
        errors=$((errors + 1)); continue
    fi
    specs+=("${geom}:${exp}:--time $(time_limit_for "$T" "$gfile")")
done < "$stage"

(( errors == 0 )) || { echo "❌ $errors problem(s) in $stage -- nothing submitted" >&2; exit 1; }
(( ${#specs[@]} > 0 )) || { echo "❌ no runs in $stage" >&2; exit 1; }

# ---------------------------------------------------------------------------
# Repo must be committed and pushed
# ---------------------------------------------------------------------------
dirty=$(git status --porcelain -- src include inputs scripts preprocess postprocess makefile 2>/dev/null || true)
head=$(git rev-parse --short HEAD)
upstream_ok=1
if git rev-parse --abbrev-ref --symbolic-full-name '@{u}' >/dev/null 2>&1; then
    git fetch -q 2>/dev/null || true
    [[ "$(git rev-parse HEAD)" == "$(git rev-parse '@{u}')" ]] || upstream_ok=0
fi

echo "============================================================"
echo "  k_eff PRODUCTION submission — $stage_name ($((${#specs[@]})) runs)"
echo "  commit        : $head"
echo "  options       : ${PRODUCTION_OPTS[*]}"
echo "  allocation    : ${TARGET_DOFS_PER_CORE} DoF/core (scripts/lib/alloc.sh)"
echo "  campaign dir  : $CAMPAIGN_DIR"
echo "============================================================"
for s in "${specs[@]}"; do echo "  ${s%%:--time *}   [${s##*--time }]"; done

if [[ -n "$dirty" ]]; then
    echo "❌ uncommitted changes under src/ inputs/ scripts/ pre/postprocess/:" >&2
    echo "$dirty" >&2
    (( dry )) || exit 1
fi
if (( ! upstream_ok )); then
    echo "❌ HEAD ($head) is not the pushed upstream -- pull/push first" >&2
    (( dry )) || exit 1
fi
if (( dry )); then echo "(dry run: nothing submitted)"; exit 0; fi

# ---------------------------------------------------------------------------
# Submit, then write the manifest into the batch folder
# ---------------------------------------------------------------------------
tmp=$(mktemp)
printf '%s\n' "${specs[@]}" > "$tmp"
log=$(mktemp)
stage_dir_name="${stage_name}__$(date +%Y-%m-%d__%H.%M.%S)"
STAGE_DIR_NAME="$stage_dir_name" "$SCRIPT_DIR/submit_batch.sh" --tag "$tag" --tests-file "$tmp" \
    --parent-dir "$CAMPAIGN_DIR" --extra-opts "${PRODUCTION_OPTS[*]}" 2>&1 | tee "$log"
sdir="$CAMPAIGN_DIR/stages/$stage_dir_name"
mkdir -p "$sdir"
{
    echo "k_eff production stage — $stage_name"
    echo "submitted : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "commit    : $(git rev-parse HEAD)"
    echo "stage file: $stage"
    echo "options   : ${PRODUCTION_OPTS[*]}"
    echo "time      : per run, from T and DoF/rank (time_limit_for); shown below"
    echo "alloc     : ${TARGET_DOFS_PER_CORE} DoF/core (+ mem_per_cpu from scripts/lib/alloc.sh)"
    echo "runs (folder name -> job id):"
    for sp in "${specs[@]}"; do
        g="${sp%%:*}"; e="${sp#*:}"; e="${e%%:*}"
        id=$(grep -A3 -F "${g}__${e}" "$log" | grep -oE "Submitted batch job [0-9]+" | head -1 | awk '{print $4}')
        refused=$(grep -F "${g}__${e} already exists" "$log" >/dev/null && echo " (EXISTS -- not resubmitted)")
        echo "  ${g}__${e}  ${id:-none}  [${sp##*--time }]${refused}"
    done
} > "$sdir/PRODUCTION_MANIFEST.txt"
cp "$stage" "$sdir/"
echo "  manifest  : $sdir/PRODUCTION_MANIFEST.txt"
rm -f "$tmp" "$log"
