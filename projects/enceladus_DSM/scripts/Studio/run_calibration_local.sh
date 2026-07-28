#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run_calibration_local.sh — the Phase-4 arms that are small enough for the Mac.
#
#   ./scripts/Studio/run_calibration_local.sh --dry-run
#   ./scripts/Studio/run_calibration_local.sh
#   ./scripts/Studio/run_calibration_local.sh molaro   # just the neck arm
#
# Split rationale. Cost is (DOF) x (steps), and the two test families sit in
# very different places:
#
#   Molaro neck (test A)   t_final 2 h / dtmax 100 s  ->  ~72 steps minimum.
#     Only 0.4-2.3 M DOF and few steps, so all three eps arms finish locally.
#
#   7-day calibration      t_final 7 d / dtmax 200 s  ->  >=3024 steps, and the
#     adaptive stepper takes more early on while dt climbs from 1e-4. Only the
#     two coarsest ripening arms are sane locally; everything else goes to
#     scripts/HPC/submit_calibration.sh.
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

# --- B. The two coarsest ripening arms -------------------------------------
if [[ "$WHICH" == "all" || "$WHICH" == "ripen" ]]; then
  for s in s100 s075; do
    cmds+=("$RUN calib_ripen_$s calib_T-20_7d calib_$s")
  done
fi

echo "=== ${#cmds[@]} local calibration jobs ==="
printf '  %s\n' "${cmds[@]}"
echo
echo "  Everything else (calib_ripen_s050/s025, all calib_pack_*) needs the"
echo "  cluster: ./scripts/HPC/submit_calibration.sh"

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
