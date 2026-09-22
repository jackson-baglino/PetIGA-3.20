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
# RESULT OF THE FIRST PAIR (2026-09-21), and why it needs repeating at finer
# eps. Both ran to t_final cleanly, 1428 steps each. Measured neck WIDTH:
#
#     t [min]      6     41    134    476    900
#     axisym    20.7   27.3   34.2   45.7   53.1  um
#     planar    20.1   24.6   28.6   40.8   48.4  um
#     ratio     1.03   1.11   1.20   1.12   1.10
#
# Two things to take from it, one of which contradicts the prediction written
# above:
#
#   * AXISYM GROWS FASTER, not slower -- by ~10-15%. The argument for planar
#     being faster considered only the neck SINK (a 3D saddle, 1/r - 1/rho,
#     partially cancels where planar has only -1/rho). It ignored the grain
#     SOURCE, which is 2/R for a sphere against 1/R for a cylinder. The
#     stronger source evidently wins here. The sign of this correction is a
#     measurement, not a derivation, which is the whole reason for the pair.
#
#   * BOTH RUNS STAY BELOW THE NECK RESOLUTION FLOOR the whole time. That floor
#     is w = 2R*sqrt(12 eps/R) = 54.9 um at eps = 0.8675 um, and the runs reach
#     53.1 and 48.4 um. neck_width.py says so itself and skips its power-law
#     fit. So the ~1.10-1.15 ratio is indicative, not quotable.
#
# Rerun at eps/4 (0.217 um, Nx = 2524) before quoting a number: that drops the
# floor to ~27 um, which both trajectories clear before the halfway point.
# Still a small run.
#
# SUBMITS THROUGH submit_batch.sh, NOT A LOOP OVER submit_enceladus.sh.
# Each submit_enceladus.sh call runs `make clean && make all` inside the job,
# so N of them race in the shared obj/. run_enceladus.sh:162-165 documents the
# failure exactly: one job's `make clean` deletes obj/*.o while another is
# mid-write. It is not hypothetical -- the first attempt at this pair produced
#
#   rm: cannot remove 'obj/monitoring.o': Stale file handle    (one job)
#   srun: error: tasks 0-3: Bus error (core dumped)            (the other)
#
# the Bus error being the same cause seen from the other side: the executable
# relinked underneath ranks that had it mapped. submit_batch.sh builds ONCE on
# the submission host and passes SKIP_COMPILE=1 to every fanned-out job.
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
echo "  axisym  $AXI"
echo "  planar  $PLA"
echo ""
cmd=(./scripts/HPC/submit_batch.sh --tag curvcal
     --tests "${AXI}:${EXP},${PLA}:${EXP}")
[ "${#sbatch_extra[@]}" -gt 0 ] && cmd+=(-- "${sbatch_extra[@]}")
if [ "$dry" -eq 1 ]; then printf '$ '; printf '%q ' "${cmd[@]}"; echo
else "${cmd[@]}"; fi
echo ""
echo "When both finish, compare neck growth:"
echo "  venv_enceladus/bin/python postprocess/neck_width.py <axisym_run>"
echo "  venv_enceladus/bin/python postprocess/neck_width.py <planar_run>"
echo "The ratio of dr/dt is the correction. Expect planar to be FASTER."
