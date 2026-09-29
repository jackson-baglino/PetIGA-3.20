#!/usr/bin/env bash
# =============================================================================
# submit_scaling_test.sh — how many DoF per core is cheapest for the k_eff runs?
#
#   ./scripts/HPC/submit_scaling_test.sh --run <finished 3a run dir> [--dry-run]
#       [--t-days 10] [--steps 30] [--keff-every 3] [--repeats 2]
#       [--targets "20000 40000 60000 100000 150000 200000"]
#
# Each job RESTARTS a finished production run from its snapshot nearest
# --t-days (so the steps are representative mid-run steps at the capped dt,
# not the tiny ramp-up steps from the initial condition), takes --steps time
# steps with a k_eff sample every --keff-every, writes no field snapshots, and
# stops. -log_view is on. Nothing here is a manuscript run.
#
# One submit_batch.sh call per DoF/core target, because the allocation is set
# per submission (TARGET_DOFS_PER_CORE, scripts/lib/alloc.sh). Each call builds
# once on the login node and its jobs skip compiling, so the calls do not race
# in obj/. --repeats jobs per target measure the node-to-node noise (batch 2
# saw 10x swings in the k_eff time per iteration within one run).
#
# WHY TWO NUMBERS. The run alternates two different solves on the same ranks:
#   phase-field step : 3 coupled fields, 24M DoF, Newton + BiCGStab + ASM/ILU(3)
#   k_eff sample     : 1 scalar field, 8M DoF, CG + GAMG
# At a given rank count the k_eff solve has a third of the DoF per rank, and
# the two preconditioners scale differently (ASM/ILU weakens with more
# subdomains; AMG is limited by its coarse levels and communication), so their
# best rank counts need not agree. The analysis reports both, and the cost per
# production run they imply at each target.
#
# THE RANGE. PETSc's rule of thumb (>= ~10-20k unknowns per rank for good
# parallel efficiency; ~20-100k the usual well-scaling band) is a statement
# about WALL-TIME efficiency. The bill is ranks x wall time, and efficiency only
# falls as ranks are added, so the cheapest point is usually at or above the
# top of the band -- bounded by memory (12% used at 60k in 3a) and the 24 h
# limit (-5 C took 16.6 h at 401 ranks). The targets therefore span 20k (the
# band's floor, where strong scaling should visibly stall) to 200k. They count
# ALL unknowns (3 fields); the k_eff solve has a third of that per rank, so it
# spans ~7k-67k.
#
# Read with: venv_enceladus/bin/python studies/keff_sintering/scaling/analyze_scaling.py <dirs>
# Cost: ~$30 for 6 targets x 2 repeats (each job ~30-40 min; 20k is 1201 ranks).
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
OUT_ROOT="/resnick/groups/rubyfu/jbaglino/simulation_outputs"

run="" ; t_days=10 ; steps=30 ; every=3 ; repeats=2 ; dry=0
targets="20000 40000 60000 100000 150000 200000"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --run) run="$2"; shift 2 ;;
        --t-days) t_days="$2"; shift 2 ;;
        --steps) steps="$2"; shift 2 ;;
        --keff-every) every="$2"; shift 2 ;;
        --repeats) repeats="$2"; shift 2 ;;
        --targets) targets="$2"; shift 2 ;;
        --dry-run) dry=1; shift ;;
        -h|--help) sed -n '2,32p' "$0"; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 1 ;;
    esac
done
[[ -d "$run" && -f "$run/SSA_evo.dat" ]] || { echo "❌ --run must be a finished run dir with SSA_evo.dat" >&2; exit 1; }
cd "$PROJECT_ROOT"

# geometry and experiment from the run folder name: <geom>__<exp>[__<label>]
name="$(basename "$run")"
geom="${name%%__*}"; rest="${name#*__}"; exp="${rest%%__*}"

# snapshot nearest t_days
t_want=$(awk -v d="$t_days" 'BEGIN{print d*86400}')
best="" ; bestd=""
for f in "$run"/sol_*.dat; do
    s=$(basename "$f" .dat); s=$((10#${s#sol_}))
    t=$(awk -v s="$s" '$4==s{print $3; exit}' "$run/SSA_evo.dat")
    [[ -z "$t" ]] && continue
    d=$(awk -v a="$t" -v b="$t_want" 'BEGIN{x=a-b; print (x<0?-x:x)}')
    if [[ -z "$bestd" ]] || awk -v a="$d" -v b="$bestd" 'BEGIN{exit !(a<b)}'; then best="$f"; bestd="$d"; bstep=$s; bt=$t; fi
done
[[ -n "$best" ]] || { echo "❌ no usable sol_*.dat in $run" >&2; exit 1; }
max_steps=$((bstep + steps))

opts="-initial_cond $best -ts_max_steps $max_steps -pf_output 0 -log_view"
opts+=" -keff 1 -keff_step0 0 -keff_freq $every -keff_ksp_type cg -keff_pc_type gamg"

echo "============================================================"
echo "  DoF/core scaling test"
echo "  restart from : $(basename "$best")  (step $bstep, t = $bt s)"
echo "  run          : $name"
echo "  steps        : $steps, k_eff every $every"
echo "  targets      : $targets  (x $repeats repeats)"
echo "============================================================"

for tgt in $targets; do
    tf=$(mktemp)
    for r in $(seq 1 "$repeats"); do
        echo "${geom}:${exp}:--label r${r} ${opts}" >> "$tf"
    done
    if (( dry )); then
        echo "--- TARGET_DOFS_PER_CORE=$tgt"; cat "$tf"; rm -f "$tf"; continue
    fi
    TARGET_DOFS_PER_CORE="$tgt" "$SCRIPT_DIR/submit_batch.sh" --tag "scaling_${tgt}" \
        --tests-file "$tf" --out-root "$OUT_ROOT" -- --time=0-03:00:00
    rm -f "$tf"
done
(( dry )) && echo "(dry run: nothing submitted)"
