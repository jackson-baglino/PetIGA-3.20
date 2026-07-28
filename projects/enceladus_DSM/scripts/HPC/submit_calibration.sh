#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# submit_calibration.sh — Phase 4: how loose can the discretisation be?
#
# Finds the cheapest (safety, xi_v) that still reproduces the observables of
# the finest run, before committing to the 200-run production matrix.
#
#   ./scripts/HPC/submit_calibration.sh heavy      # ONLY what needs a cluster
#   ./scripts/HPC/submit_calibration.sh            # submit everything
#   ./scripts/HPC/submit_calibration.sh --dry-run  # print the commands only
#   ./scripts/HPC/submit_calibration.sh safety     # just the safety sweep
#   ./scripts/HPC/submit_calibration.sh xiv        # just the xi_v sweep
#   ./scripts/HPC/submit_calibration.sh molaro     # just the neck-width arm
#
# `heavy` omits only the three Molaro geometries, which need ~125 steps and
# finish in 0.6 h on the Mac -- run those with
# scripts/Studio/run_calibration_local.sh.
#
# Every 7-day arm belongs here. On <=0.5 mm domains they all FIT in the
# workstation's memory except ripen_s025 (10 GB) and pack_s025 (21 GB), but at
# a measured 7.5 s/step per 1.5 M DOF on 6 ranks they are ~55 h serial for the
# safety sweep alone. Memory is not the binding constraint; step count is.
#
# Every job is chained with && so a failure stops the sequence rather than
# burning the whole allocation on a broken build.
#
# TWO KNOBS, OPPOSITE DIRECTIONS
#   safety  eps = 2*safety um (the kinetic bound binds at R_feat = 2 um).
#           Cost ~ 1/eps^2, so s100 is 16x cheaper than s025.
#   xi_v    assembly.c solves (1/xi_v)*d(phi_a*rhov)/dt = ..., so 1/xi_v
#           amplifies STORAGE: SMALLER xi_v slows vapour relaxation and allows
#           a larger dt. Small is the cheap direction. The risk is at the small
#           end, where tau_vap = L^2/(xi_v*D_v) stops being negligible against
#           the ice-evolution time.
# ---------------------------------------------------------------------------
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/../.."

WHICH="${1:-all}"
DRY=0
[[ "${1:-}" == "--dry-run" ]] && { DRY=1; WHICH="${2:-all}"; }

SUB=./scripts/HPC/submit_enceladus.sh
EXP=calib_T-20_7d

cmds=()

# --- A. Two-grain sintering: neck width vs time -----------------------------
# Reuses the existing Molaro axisymmetric union geometries, whose neck width is
# held FIXED across eps by preprocess/calibrate_neck_geometry.py -- so any
# change in the measured neck growth is a discretisation effect, not a
# different starting geometry. eps = 8.584e-7 / 6.030e-7 / 3.479e-7.
if [[ "$WHICH" == "all" || "$WHICH" == "molaro" ]]; then   # not in `heavy`
  for e in epsloose epsmid epsstrict; do
    cmds+=("$SUB 2D_molaro_axisym_T-20pair_union_$e molaro_T-20_h1.00_arrh calib_$e")
  done
fi

# --- B/C. Safety sweep at the default xi_v ---------------------------------
if [[ "$WHICH" == "all" || "$WHICH" == "safety" || "$WHICH" == "heavy" ]]; then
  for s in s025 s050 s075 s100; do
    cmds+=("$SUB calib_ripen_$s $EXP calib_${s}")
  done
  for s in s025 s050 s075 s100; do
    cmds+=("$SUB calib_pack_$s $EXP calib_${s}")
  done
fi

# --- B/C. xi_v sweep at safety = 0.5 ---------------------------------------
# 1e-3 is the solver default and is already covered by the safety sweep's
# s050 arms, so it is not repeated here.
if [[ "$WHICH" == "all" || "$WHICH" == "xiv" || "$WHICH" == "heavy" ]]; then
  for x in 1e-4 1e-2 1e-1; do
    tag="xiv${x}"
    cmds+=("$SUB calib_ripen_s050 $EXP $tag -- -xi_v $x")
    cmds+=("$SUB calib_pack_s050  $EXP $tag -- -xi_v $x")
  done
fi

if (( ${#cmds[@]} == 0 )); then
  echo "Nothing selected. Use: heavy | all | safety | xiv | molaro" >&2
  exit 1
fi

echo "=== ${#cmds[@]} calibration jobs ==="
printf '  %s\n' "${cmds[@]}"

if (( DRY )); then
  echo
  echo "(dry run — nothing submitted)"
  exit 0
fi

echo
joined=""
for c in "${cmds[@]}"; do
  joined+="${joined:+ && }$c"
done
eval "$joined"
