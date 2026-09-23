#!/usr/bin/env bash
# =============================================================================
# verify_keff_sharp_limit.sh
#
# Measure the finite-eps bias of the phase-field homogenization against the
# closed form derived in
#   projects/effective_thermal_cond/docs/calonne_to_phasefield_equivalence.tex
#
# WHAT THIS TESTS, AND WHAT IT DOES NOT
# -------------------------------------
# The note proves that the sharp (Calonne) and phase-field cell problems are
# the SAME problem -- not in a limit, but identically, because the sharp
# coefficient K_star and the interpolation K(phi) are one function evaluated at
# two arguments (the sharp indicator vs. the smoothed phi^eps). That part needs
# no test; it is algebra.
#
# What DOES need a test is the residual finite-eps bias, which the note gives in
# closed form for a laminate. This driver measures it. Three predictions, each
# failing in a diagnostically different way:
#
#   (i)   k_00 (parallel) FLAT across the ladder, at phi*K_i + (1-phi)*K_a.
#         The tangential surface excess is exactly zero at every eps, so a drift
#         here is a REAL DEFECT in the cell solver or the IC -- not eps bias.
#   (ii)  1/k_11 (perpendicular) LINEAR in eps, with intercept <1/K>_sharp and
#         slope (n_Gamma/L)(1/K_a - 1/K_i)ln(K_i/K_a). Both predicted, neither
#         fitted. This is the test of the theory.
#   (iii) k_01, k_10 ~ 0. Off-diagonal isotropy check on the cell solver.
#
# RESOLUTION SCALES WITH EPS -- this is not optional
# --------------------------------------------------
# Ny is set to hold eps/dy fixed at EPS_PER_ELEM, anchored on the geometry
# file's reference point (eps = 2e-5 at Ny = 256, i.e. eps/dy = 5.12). At FIXED
# resolution an eps ladder confounds the interface bias being measured with
# ordinary mesh convergence, and the fitted slope means nothing. If the gate
# reports a bad slope but a good intercept, suspect this first.
#
# COST: each rung is -keff_only, which samples the IC once, solves dim scalar
# Poisson problems, and exits. There is no time integration. Even the 2560^2
# rung is minutes, so this runs locally; it is not an HPC job.
#
# USAGE
#   ./verify_keff_sharp_limit.sh [--rungs "50 64 128 256 512"]
#                                [--interp arith|tensor] [--dry-run]
#
# --interp selects -keff_interp. Under "tensor" (arithmetic along the interface,
# harmonic across) prediction (ii) becomes a FLAT 1/k_11 at the sharp value:
# the first-order bias is what the tensor law exists to remove.
#
# Writes keff_sharp_limit.csv (keff_sharp_limit_tensor.csv under --interp
# tensor) next to this script, then calls plot_keff_sharp_limit.py.
# Re-runnable: the CSV is rewritten from scratch.
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$HERE/../../.." && pwd)"
RUNNER="$PROJECT_ROOT/scripts/Studio/run_enceladus.sh"
PYTHON="$PROJECT_ROOT/venv_enceladus/bin/python"

GEOM="iceslab_2D_L1mm_eps20um_keff"
EXP="snow_T-20_h1.00_1d"     # -keff_only exits before integrating; conditions are inert

# Domain and interface count must match the geometry file and the analytic
# module's assumptions. n_Gamma = 2 because -periodic 1 makes the y=0 seam a
# real second interface.
L=1.0e-3
N_GAMMA=2

# Reference point from the geometry file: eps = 2.0e-5 at Ny = 256 on L = 1 mm.
EPS_PER_ELEM=5.12

RUNGS="50 64 128 256 512"
INTERP="arith"
DRY_RUN=0
while [ $# -gt 0 ]; do
    case "$1" in
        --rungs)   RUNGS="$2"; shift 2 ;;
        --interp)  INTERP="$2"; shift 2 ;;
        --dry-run) DRY_RUN=1; shift ;;
        -h|--help) sed -n '2,52p' "${BASH_SOURCE[0]}"; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

case "$INTERP" in
    arith)  SUFFIX="" ;;
    tensor) SUFFIX="_tensor" ;;
    *) echo "--interp must be arith or tensor, got: $INTERP" >&2; exit 2 ;;
esac
CSV="$HERE/keff_sharp_limit$SUFFIX.csv"
LOG="$HERE/keff_sharp_limit$SUFFIX.log"

[ -x "$RUNNER" ] || { echo "run script not found or not executable: $RUNNER" >&2; exit 1; }
[ -x "$PYTHON" ] || { echo "venv python not found: $PYTHON" >&2; exit 1; }

# -- gate 0: the analytic module must reproduce the two constants the geometry
# file states independently. If this fails the module is wrong, not the solver,
# and every number downstream is meaningless. Check it before burning any runs.
echo "=== Gate 0: analytic module self-check ==="
"$PYTHON" "$HERE/keff_laminate_analytic.py"
echo ""

# Resolve the schedule first, so --dry-run can print it without touching any
# file. Initialising the CSV before this point would let a dry run clobber a
# real ladder, which is the opposite of what --dry-run is for.
schedule=""
for denom in $RUNGS; do
    eps=$("$PYTHON" -c "print(f'{$L/$denom:.8e}')")
    # Hold eps/dy fixed: Ny = EPS_PER_ELEM * L / eps. Rounded to even so the
    # slab's two interfaces sit symmetrically on the grid.
    Ny=$("$PYTHON" -c "n=round($EPS_PER_ELEM*$L/$eps); print(n+(n%2))")
    schedule="$schedule$denom $eps $Ny"$'\n'
done

if [ "$DRY_RUN" -eq 1 ]; then
    echo "=== ladder schedule (dry run -- nothing executed, no files written) ==="
    printf '%s' "$schedule" | while read -r denom eps Ny; do
        echo "  L/$denom :  eps = $eps   Nx = Ny = $Ny"
    done
    exit 0
fi

echo "step,time,k_00,k_01,k_10,k_11,phi_bar,k_iso,ksp_its,ksp_reason,wall_s,eps,Ny,run_dir" > "$CSV"
: > "$LOG"

printf '%s' "$schedule" | while read -r denom eps Ny; do
    echo "=== rung L/$denom ($INTERP) :  eps = $eps   Nx = Ny = $Ny ==="
    # -keff_ksp/-keff_pc: CG+GAMG rather than the direct LU that is only
    # comfortable up to ~256^2 (see the geometry file header).
    out=$("$RUNNER" "$GEOM" "$EXP" "sharplimit_${INTERP}_L$denom" -- \
            -keff 1 -keff_only 1 -keff_interp "$INTERP" \
            -eps "$eps" -Nx "$Ny" -Ny "$Ny" \
            -keff_ksp_type cg -keff_pc_type gamg 2>&1 | tee -a "$LOG")

    run_dir=$(echo "$out" | sed -n 's/^Output folder: //p' | tail -1)
    [ -n "$run_dir" ] || { echo "could not parse output folder from run script" >&2; exit 1; }

    kcsv="$run_dir/k_eff$SUFFIX.csv"
    [ -f "$kcsv" ] || { echo "no $(basename "$kcsv") in $run_dir" >&2; exit 1; }

    # -keff_only writes exactly one sample row after the header.
    row=$(tail -n +2 "$kcsv" | tail -1)
    echo "$row,$eps,$Ny,$run_dir" >> "$CSV"
    echo "    -> $row"
done

echo ""
echo "=== Gates 1-3: measured vs. predicted ==="
"$PYTHON" "$HERE/plot_keff_sharp_limit.py" --csv "$CSV" --L "$L" --n-gamma "$N_GAMMA" --interp "$INTERP"
