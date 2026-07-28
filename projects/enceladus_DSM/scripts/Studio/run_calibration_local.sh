#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run_calibration_local.sh — the Phase-4 arms that are small enough for the Mac.
#
#   ./scripts/Studio/run_calibration_local.sh --dry-run
#   ./scripts/Studio/run_calibration_local.sh          # run + post-process
#   ./scripts/Studio/run_calibration_local.sh molaro   # just the neck arm
#   ./scripts/Studio/run_calibration_local.sh --no-post
#
# After the runs it measures the neck on each and writes the eps-convergence
# figure + table. plot_neck_convergence.py keeps only the NEWEST run per arm,
# so re-running a batch compares it against itself rather than mixing builds.
#
# Split rationale, MEASURED (not guessed -- an earlier guess put a 24 M DOF
# arm on the workstation and exhausted its 64 GB).
#
# Memory: preprocess/estimate_memory.py. The ILU(3) fill dominates, roughly
#   3 x the Jacobian, itself DOF x (2p+1)^2 x dof x 12 B. All arms on <=0.5 mm
#   domains fit except the two finest (ripen_s025 10 GB, pack_s025 21 GB).
#
# Time: benchmarked at 7.5 s/step for 1.5 M DOF on 6 ranks, scaling ~linearly
#   in DOF. Step COUNT is what separates the two families:
#
#   Molaro neck (test A)   t_final 2 h / dtmax 100 s -> ~125 steps including
#     the ramp.  0.6 h for all three arms. LOCAL.
#
#   7-day calibration      t_final 7 d / dtmax 200 s -> ~3075 steps (the
#     adaptive stepper needs ~50 just to climb from dt = 1e-4 to dtmax).
#     55 h serial for the safety sweep alone, before the xi_v arms. HPC --
#     these fit in memory but not in patience.
#
# Runs go through run_enceladus.sh (not the bare binary) so each one stages its
# own solver.opts + geometry .opts -- which is also what effective_thermal_cond
# reads back to match the mesh.
# ---------------------------------------------------------------------------
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/../.."

WHICH="all"
DRY=0
POST=1
for a in "$@"; do
  case "$a" in
    --dry-run) DRY=1 ;;
    --no-post) POST=0 ;;
    -*)        echo "unknown flag: $a" >&2; exit 1 ;;
    *)         WHICH="$a" ;;
  esac
done

RESULTS="${RESULTS_BASE:-$HOME/SimulationResults/enceladus_DSM/scratch}"
PY_ENV="./venv_enceladus/bin/python3"
[[ -x "$PY_ENV" ]] || PY_ENV="python3"

RUN=./scripts/Studio/run_enceladus.sh
cmds=()

# --- A. Molaro neck width: all three eps arms ------------------------------
if [[ "$WHICH" == "all" || "$WHICH" == "molaro" ]]; then
  for e in epsloose epsmid epsstrict; do
    cmds+=("$RUN 2D_molaro_axisym_T-20pair_union_$e molaro_T-20_h1.00_arrh calib_$e")
  done
fi

# The 7-day arms are deliberately NOT here -- see the timing note above.

echo "=== ${#cmds[@]} local calibration jobs ==="
printf '  %s\n' "${cmds[@]}"
echo
echo "  All 7-day arms (calib_ripen_*, calib_pack_*) go to the cluster:"
echo "    ./scripts/HPC/submit_calibration.sh heavy"
echo "  They fit in memory on <=0.5 mm domains but need ~55 h serial here."

if (( POST )); then
  echo
  echo "  then post-process:"
  echo "    $PY_ENV postprocess/plot_neck_convergence.py --runs $RESULTS/*/*/ \\"
  echo "        --validation inputs/validation/molaro2019_fig11_T-20.csv \\"
  echo "        --out $RESULTS/neck_convergence.png"
fi

if (( DRY )); then
  echo
  echo "(dry run — nothing executed)"
  exit 0
fi

echo
joined=""
for c in "${cmds[@]}"; do
  joined+="${joined:+ && }$c"
done
eval "$joined"

if (( POST )); then
  echo
  echo "============================================================"
  echo "  Neck-width convergence"
  echo "============================================================"
  $PY_ENV postprocess/plot_neck_convergence.py \
      --runs "$RESULTS"/*/*/ \
      --validation inputs/validation/molaro2019_fig11_T-20.csv \
      --out "$RESULTS/neck_convergence.png"
fi
