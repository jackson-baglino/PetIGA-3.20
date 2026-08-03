#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# submit_rev_sweep.sh — the REV domain-size sweep.
#
# Six nested, centred windows cut from ONE master packing (3.5 mm, phi = 0.325,
# 100 um grains, tangent contacts). Each window is a strict subset of the next,
# so the microstructure is identical and the ONLY thing changing is how much of
# it is being averaged. k_eff(L) then measures size, not realisation scatter --
# which matters, because two independent 1 mm packings at the same target
# porosity differed by 48% in k_iso in this project's own tests.
#
#   L [mm]     Nx     DOF        cores    grains   phi (measured on the window)
#   0.5       393     463,347      6        9      0.4144
#   1.0       786   1,853,388     24       26      0.3864
#   1.5      1179   4,170,123     53       60      0.3781
#   2.0      1572   7,413,552     93      104      0.3379
#   2.5      1965  11,583,675    145      158      0.3257
#   3.0      2358  16,680,492    209      219      0.3237
#                                 530 core-slots total, submitted smallest first
#
# READ k_eff(L) TOGETHER WITH THAT POROSITY COLUMN. It drifts from 0.414 to
# 0.324 as the window grows, converging on the master's 0.325 -- real, not an
# error, but k_eff responds far more strongly to porosity than to L, so a
# k_eff(L) trend that merely tracks the porosity drift is NOT a size effect.
# The porosity itself is converged by about 2.5 mm, which is a first hint at
# where the REV sits.
#
# The sweep is submitted smallest first so a failure shows up cheaply.
#
# Usage:  ./submit_rev_sweep.sh [tag] [-- extra solver flags]
# ---------------------------------------------------------------------------
set -uo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/../.."
PROJECT_ROOT="$PWD"
SUBMIT="$PROJECT_ROOT/scripts/HPC/submit_enceladus.sh"

TAG="${1:-rev}"; shift || true
[[ "${1:-}" == "--" ]] && shift || true
EXTRA=("$@")

EXP="snow_T-20_rev1day"
SIZES=(0.5 1 1.5 2 2.5 3)

# k_eff every 20 steps: enough points to see the trajectory without the
# corrector solve dominating. The wall_s column in k_eff.csv reports the real
# cost -- retune from it after the first run rather than guessing again.
KEFF=(-keff 1 -keff_freq 20)

[[ -x "$SUBMIT" ]] || { echo "ERROR: $SUBMIT missing" >&2; exit 1; }
[[ -f "$PROJECT_ROOT/inputs/experiment/$EXP.opts" ]] || {
    echo "ERROR: inputs/experiment/$EXP.opts missing" >&2; exit 1; }

echo "=== REV sweep: ${#SIZES[@]} nested windows, one master packing ==="
for L in "${SIZES[@]}"; do
    geom="rev_L${L}mm_T-20_rev"
    if [[ ! -f "$PROJECT_ROOT/inputs/geometry/${geom}.opts" ]]; then
        echo "  MISSING inputs/geometry/${geom}.opts -- skipping" >&2
        continue
    fi
    echo
    echo "--- L = ${L} mm ---"
    "$SUBMIT" "$geom" "$EXP" "$TAG" "${KEFF[@]}" ${EXTRA[@]+"${EXTRA[@]}"}
done

echo
echo "When they finish, collect with:"
echo "  venv_enceladus/bin/python3 studies/snow_thermal/verification/plot_rev.py \\"
echo "      --results \$RESULTS_BASE --tag $TAG"
