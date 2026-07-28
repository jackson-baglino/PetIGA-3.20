#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run_calibration_local.sh — the Phase-4 arms that are small enough for the Mac.
#
#   ./scripts/Studio/run_calibration_local.sh --dry-run
#   ./scripts/Studio/run_calibration_local.sh
#   ./scripts/Studio/run_calibration_local.sh molaro   # just the neck arm
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

WHICH="${1:-all}"
DRY=0
[[ "${1:-}" == "--dry-run" ]] && { DRY=1; WHICH="${2:-all}"; }

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
