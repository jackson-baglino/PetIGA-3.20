#!/usr/bin/env bash
# =============================================================================
# submit_curvature_calib.sh — measure the 2D -> 3D curvature correction.
#
# THE PROBLEM IT SOLVES. Everything in this model is driven by curvature: the
# Gibbs-Thomson term sets the surface vapour pressure, which sets the whole
# sublimation/deposition cycle. A surface in 3D has TWO principal curvatures
# and the driving force goes as 1/R1 + 1/R2; a curve in a plane has ONE. So a
# planar 2D packing is not a scaled 3D packing, it is missing half the driving
# term -- and not uniformly:
#
#   sphere   3D 2/R          planar disc   1/R        factor 2
#   neck     3D 1/r - 1/rho  planar        -1/rho     the +1/r that partially
#                                                     CANCELS in 3D is gone, so
#                                                     a planar neck sees a
#                                                     stronger sink and grows
#                                                     FASTER
#
# Sintering RATES from a planar packing are therefore fast, by an amount no
# planar run can determine from itself. This pair determines it.
#
# HOW. Two runs, identical in every respect except how the solver interprets
# the out-of-plane direction:
#
#   -axisym 1   revolves the profile -> a pair of SPHERES (true 3D curvature)
#   -axisym 0   extrudes the profile -> a pair of CYLINDERS (the planar case)
#
# Same two radii, exactly tangent, same eps, same mesh spacing, same kinetics.
# The axisym domain is a half plane with the grains on the axis, so its planar
# twin doubles Ly and moves the centres to the mid-plane; Ny doubles to hold
# h = eps/sqrt(2). Nothing else differs.
#
# The ratio of neck growth rates is then the correction, and it is measured
# rather than asserted. Extract it with postprocess/neck_width.py on both runs.
#
# WHAT IT DOES AND DOES NOT LICENCE. It gives a factor for the RATE of neck
# growth in a two-grain geometry. It does not make a planar packing into a 3D
# packing: coordination, pore topology and the percolation behaviour of the
# two phases all differ independently of curvature. Use it to state how far
# the timing is off, not to convert 2D results into 3D ones.
#
# These are small (631x198 and 631x396), so they are cheap next to anything
# else in the campaign.
#
# USAGE
#   ./scripts/HPC/submit_curvature_calib.sh [--dry-run] [-- <sbatch flags>]
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$(cd "$SCRIPT_DIR/../.." && pwd)"

AXI="molaro_2D_L387um_eps0.87um_axisym_T-20pair_tangent"
PLA="molaro_2D_L387x242um_eps0.87um_planar_T-20pair_tangent"
EXP="molaro_T-20_h1.00_15h_a1.34e-2"

dry=0; sbatch_extra=()
while [ $# -gt 0 ]; do
    case "$1" in
        --dry-run) dry=1; shift ;;
        --)        shift; sbatch_extra=("$@"); break ;;
        -h|--help) sed -n '2,42p' "${BASH_SOURCE[0]}"; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done
if [ "$dry" -eq 0 ] && ! command -v sbatch >/dev/null 2>&1; then
    echo "ERROR: sbatch not found -- run this on a cluster login node." >&2; exit 1
fi

echo "=== 2D -> 3D curvature calibration: one geometry, two interpretations ==="
for pair in "axisym:$AXI" "planar:$PLA"; do
    tag="${pair%%:*}"; geom="${pair#*:}"
    echo "  $tag  $geom"
    cmd=(./scripts/HPC/submit_enceladus.sh "$geom" "$EXP" "curvcal_$tag")
    [ "${#sbatch_extra[@]}" -gt 0 ] && cmd+=("${sbatch_extra[@]}")
    if [ "$dry" -eq 1 ]; then printf '      $ '; printf '%q ' "${cmd[@]}"; echo
    else "${cmd[@]}"; fi
done
echo ""
echo "When both finish, compare neck growth:"
echo "  venv_enceladus/bin/python postprocess/neck_width.py <axisym_run>"
echo "  venv_enceladus/bin/python postprocess/neck_width.py <planar_run>"
echo "The ratio of dr/dt is the correction. Expect planar to be FASTER."
