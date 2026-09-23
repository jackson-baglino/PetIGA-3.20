#!/usr/bin/env bash
# =============================================================================
# verify_keff_disk.sh
#
# Curved-interface test of -keff_interp tensor: a square array of ice disks,
# whose sharp k_eff is known in closed form (Rayleigh / Perrins et al. 1979).
# Runs the SAME eps ladder under both laws so they can be compared rung by rung:
#
#   arith   first-order biased high; slope predicted with no fitted constants
#   tensor  first-order bias removed; residual O(eps^2)
#
# The laminate ladder (../verification/) cannot test the tensor law's
# off-diagonal terms, because the laminate's normal is exactly the y axis. Here
# the normal takes every direction.
#
# Each rung is -keff_only: one sample from the IC, dim scalar solves, exit. No
# time integration, so it runs locally.
#
# USAGE
#   ./verify_keff_disk.sh [--rungs "50 100 200 400"] [--dry-run]
#
# Writes keff_disk.csv next to this script, then calls gate_keff_disk.py.
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$HERE/../../.." && pwd)"
RUNNER="$PROJECT_ROOT/scripts/Studio/run_enceladus.sh"
PYTHON="$PROJECT_ROOT/venv_enceladus/bin/python"

GEOM="singleice_2D_L1mm_R250um_keff"
EXP="snow_T-20_h1.00_1d"     # -keff_only exits before integrating; conditions are inert
CSV="$HERE/keff_disk.csv"
LOG="$HERE/keff_disk.log"

L=1.0e-3
R=2.5e-4
# Resolution scales with eps, as in the laminate ladder: eps/dy held at 5.12.
EPS_PER_ELEM=5.12

RUNGS="50 100 200 400"
DRY_RUN=0
while [ $# -gt 0 ]; do
    case "$1" in
        --rungs)   RUNGS="$2"; shift 2 ;;
        --dry-run) DRY_RUN=1; shift ;;
        -h|--help) sed -n '2,24p' "${BASH_SOURCE[0]}"; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

[ -x "$RUNNER" ] || { echo "run script not found or not executable: $RUNNER" >&2; exit 1; }
[ -x "$PYTHON" ] || { echo "venv python not found: $PYTHON" >&2; exit 1; }

echo "=== Gate 0: analytic module self-check ==="
"$PYTHON" "$HERE/keff_disk_analytic.py"
echo ""

schedule=""
for denom in $RUNGS; do
    eps=$("$PYTHON" -c "print(f'{$L/$denom:.8e}')")
    # Even Ny keeps the disk centre on a mesh line.
    Ny=$("$PYTHON" -c "n=round($EPS_PER_ELEM*$L/$eps); print(n+(n%2))")
    schedule="$schedule$denom $eps $Ny"$'\n'
done

if [ "$DRY_RUN" -eq 1 ]; then
    echo "=== ladder schedule (dry run -- nothing executed, no files written) ==="
    printf '%s' "$schedule" | while read -r denom eps Ny; do
        echo "  L/$denom :  eps = $eps   Nx = Ny = $Ny   (arith and tensor)"
    done
    exit 0
fi

echo "step,time,k_00,k_01,k_10,k_11,phi_bar,k_iso,ksp_its,ksp_reason,wall_s,eps,Ny,interp,run_dir" > "$CSV"
: > "$LOG"

printf '%s' "$schedule" | while read -r denom eps Ny; do
    for interp in arith tensor; do
        echo "=== rung L/$denom ($interp) :  eps = $eps   Nx = Ny = $Ny ==="
        out=$("$RUNNER" "$GEOM" "$EXP" "disk_${interp}_L$denom" -- \
                -keff 1 -keff_only 1 -keff_interp "$interp" \
                -eps "$eps" -Nx "$Ny" -Ny "$Ny" \
                -keff_ksp_type cg -keff_pc_type gamg 2>&1 | tee -a "$LOG")

        run_dir=$(echo "$out" | sed -n 's/^Output folder: //p' | tail -1)
        [ -n "$run_dir" ] || { echo "could not parse output folder from run script" >&2; exit 1; }
        # Non-default laws write k_eff_<law>.csv (see KeffCreate).
        if [ "$interp" = arith ]; then kcsv="$run_dir/k_eff.csv"
        else kcsv="$run_dir/k_eff_$interp.csv"; fi
        [ -f "$kcsv" ] || { echo "no $(basename "$kcsv") in $run_dir" >&2; exit 1; }

        row=$(tail -n +2 "$kcsv" | tail -1)
        echo "$row,$eps,$Ny,$interp,$run_dir" >> "$CSV"
        echo "    -> $row"
    done
done

echo ""
echo "=== Gates 1-4 ==="
"$PYTHON" "$HERE/gate_keff_disk.py" --csv "$CSV" --L "$L" --R "$R"
