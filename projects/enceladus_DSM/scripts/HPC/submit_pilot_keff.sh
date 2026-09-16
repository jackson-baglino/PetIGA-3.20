#!/usr/bin/env bash
# =============================================================================
# submit_pilot_keff.sh — the 4-run pilot for the sintering -> k_eff study.
#
# WHAT THIS PILOT IS FOR. It is calibration, not science. Nothing in the main
# campaign should be submitted until its numbers are read, because it sets all
# four of the quantities the campaign design depends on:
#
#   t_final          from where k_eff(t) and SSA(t) flatten
#   L/R_ave(0)       from the coarsening factor g = R_ave(t_f)/R_ave(0).
#                    Coarsening grows R_ave, so L/R_ave FALLS during the run;
#                    starting at the measured REV requirement of 40 and
#                    coarsening 2x ends at 20, below it. Need L/R_ave(0) >= 40g.
#   seeds needed     from the scatter across these 4
#   cost per run     from the run itself
#
# CONFIGURATION, and why each number is what it is:
#
#   L/R_ave = 40      measured REV requirement. At 40 the seed-to-seed scatter
#                     in k_eff is 4.3%; at 10-20 it is 12-13%, which would
#                     swamp any trend. See studies/packing_design/rev_bias.csv.
#   R_ave = 50 um     realistic snow/firn grain. The MESH does not depend on
#                     this -- it is set by L/R_ave and R_feat/R_ave, both
#                     scale-invariant -- only the physical timescale does.
#   phi = 0.325       the study's standard porosity. Solid percolation fails
#                     between 0.40 and 0.45, so 0.40 is the usable ceiling.
#   R_feat = R_ave/25 THE COARSEST DEFENSIBLE MESH, which is what sets cost.
#                     eps = safety*R_feat = 1 um, Nx = 2829, 24M DOF,
#                     ~301 cores/run. See the caveat below.
#   alpha_c = 1e-3    attachment-limited (L* = 144 um > R_ave), the regime
#                     intended for the production runs. Note it does NOT
#                     affect the mesh -- see below -- only the timescale.
#   T = -20 C         centre of the campaign's temperature axis.
#
# TWO THINGS WORTH KNOWING BEFORE READING THE RESULTS
# ---------------------------------------------------
# 1. alpha_c CANNOT COARSEN THE MESH. With --vn_feature, comp_eps.py derives
#    v_n from R_feat, which makes the K&P Eq.(45) kinetic bound exactly equal
#    R_feat -- it is binding at every alpha_c from 1e-5 to 3e-2, and eps comes
#    out as safety*R_feat regardless. The mesh is set by R_feat alone. What
#    alpha_c does control is RUNTIME: the simulated time to reach a given
#    sintering state goes as 1/alpha_c, so alpha_c = 1e-3 costs ~13x the steps
#    of the 1.34e-2 used previously. That is the real price of the
#    attachment-limited regime, and it is in wall-clock, not in memory.
#
# 2. THE COARSE MESH COSTS NECK RESOLUTION. A neck is only trustworthy above
#    r/R = sqrt(12*eps/R). At R_feat = R_ave/25 that floor is r/R = 0.49 --
#    half the grain radius -- so most of the early sintering trajectory is
#    unresolved. That is acceptable for a pilot measuring timescales and cost,
#    and is NOT acceptable for production runs about neck growth. Expect to
#    move to R_feat = R_ave/50 (floor 0.35, 96M DOF, ~1200 cores) for the
#    campaign; the pilot's job includes confirming that.
#
# k_eff IS ALREADY IN-LINE -- there is no separate code to chain
# ------------------------------------------------------------
# The standalone projects/effective_thermal_cond is superseded. The
# homogenization lives in src/keff*.c inside enceladus_dsm and is switched on
# with -keff 1, which samples k_eff on simulated-time intervals DURING the run
# and appends to k_eff.csv. It builds its OWN periodic corrector mesh, so it
# gets the periodic cell homogenization requires whatever the main run does.
#
# If you would rather compute k_eff as a separate step afterwards -- e.g. to
# re-run it at a different eps without repeating the sintering -- the snapshots
# support that:
#
#   ./scripts/HPC/run_enceladus.sh <geom> <exp> -- -keff 1 -keff_replay <rundir>
#
# Both paths write the same k_eff.csv schema. This script uses the in-line
# form, so one job per seed rather than two.
#
# USAGE
#   ./scripts/HPC/submit_pilot_keff.sh [--tag <tag>] [--dry-run] [-- <sbatch flags>]
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT"

EXP="snow_T-20_h1.00_30d"
GEOM_STEM="packing_2D_pilot_phi0.325_Rave50um_LR40_seed%d_L2mm_eps1000nm_perxy_T-20"

# k_eff sampling. -keff_t_interv is in SIMULATED seconds, not steps, so the
# trajectory is evenly sampled in time whatever the adaptive dt does. 5e4 s
# over t_final = 2.592e6 s gives ~52 samples; each costs dim scalar Poisson
# solves on the corrector mesh, which is small beside the nonlinear transient.
# -keff_step0 records t = 0 -- needed as the reference for the bias, NOT as a
# physical baseline (at t = 0 the band is standing in for necks that have not
# grown, so k_eff(0) is an artifact by construction).
KEFF_OPTS="-keff 1 -keff_step0 1 -keff_t_interv 5.0e4 -keff_ksp_type cg -keff_pc_type gamg"

tag="pilot_keff"
dry=0
sbatch_extra=()
while [ $# -gt 0 ]; do
    case "$1" in
        --tag)     tag="$2"; shift 2 ;;
        --dry-run) dry=1; shift ;;
        --)        shift; sbatch_extra=("$@"); break ;;
        -h|--help) sed -n '2,80p' "${BASH_SOURCE[0]}"; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

tests=""
for s in 1 2 3 4; do
    geom=$(printf "$GEOM_STEM" "$s")
    [ -n "$tests" ] && tests="${tests},"
    tests="${tests}${geom}:${EXP}"
done

echo "=== pilot: 4 seeds, phi = 0.325, L/R_ave = 40, R_ave = 50 um, T = -20 C ==="
echo "  alpha_c  1e-3   (attachment-limited; does NOT change the mesh)"
echo "  eps      1 um   (= safety * R_feat, R_feat = R_ave/25)"
echo "  mesh     2829^2 = 24M DOF, ~301 cores/run"
echo "  k_eff    in-line, ~52 samples over 30 simulated days"
echo "  neck floor r/R = 0.49 -- pilot-grade; see the header"
echo ""
echo "  tests: $tests"
echo ""

cmd=(./scripts/HPC/submit_batch.sh --tag "$tag" --tests "$tests"
     --extra-opts "$KEFF_OPTS")
[ "${#sbatch_extra[@]}" -gt 0 ] && cmd+=(-- "${sbatch_extra[@]}")

if [ "$dry" -eq 1 ]; then
    printf '%q ' "${cmd[@]}"; echo
    echo "(dry run: not submitted)"
    exit 0
fi
exec "${cmd[@]}"
