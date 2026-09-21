#!/usr/bin/env bash
# =============================================================================
# submit_rev_check.sh — one run at L/R_ave = 64, to settle the domain size
# before the main campaign commits to it.
#
# WHY. The pilot measured two things that together say 2 mm is not big enough:
#
#   coarsening   g = R_ave(t_f)/R_ave(0) = 1.593-1.602, remarkably consistent
#                across all four seeds, so L/R_ave falls 40 -> 25 during a run
#   divergence   seed-to-seed CV in k_eff grows 1.4% -> 3.0%, i.e. independent
#                realisations agree LESS as the structure coarsens into the box
#
# The REV requirement measured at t = 0 was L/R_ave = 40. Holding that at
# t_final needs L/R_ave(0) >= 40 * 1.59 = 64, which is this run.
#
# WHAT IT DECIDES. Compare k_eff(t_final) here against the L/R_ave = 40
# ensemble mean (0.7570 +/- 0.0229 sd, 1.5% standard error on 4 seeds):
#
#   inside the SEM   -> 40 was adequate after all, and the campaign saves 2.56x
#   outside it       -> 40 is biased, and that is worth knowing before 80 runs
#                       are spent on it rather than after
#
# COST. Nx = 4526, 61.5M DOF, ~769 cores -- 2.56x the pilot per run, because
# DOF goes as (L/eps)^2. One run, not a campaign.
#
# FOUR SEEDS, NOT ONE. The first attempt ran a single L/R_ave = 64 seed. It
# came out +5.3% above the L/R_ave = 40 mean at t=0, which looked like a
# domain-size bias but could not be claimed as one: a single realisation
# cannot be separated from seed scatter, and that packing also carried a
# coordination at the band scale of 3.566, +2.9 sd above the L/R_ave = 40 mean
# -- more contacts being a direct mechanism for higher k_eff, independent of
# domain size. Comparing ENSEMBLES removes both objections.
#
# All four are now gated on --max-void-per-L 0.0335, the size-INDEPENDENT
# criterion, which is also what the L/R_ave = 40 packings satisfy
# (0.0303-0.0328). The earlier run had to loosen the mean-radii gate, so it was
# accepted under a different rule than the set it was compared against; that
# confound is gone. The four sit at 0.0211-0.0251, i.e. relatively cleaner than
# the smaller domains. See studies/packing_design/README.md section 6.
#
# USAGE
#   ./scripts/HPC/submit_rev_check.sh [--dry-run] [-- <sbatch flags>]
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$(cd "$SCRIPT_DIR/../.." && pwd)"

GEOM_STEM="packing_2D_rev_phi0.325_Rave50um_LR64_seed%d_L3.2mm_eps1000nm_perxy_T-20"
SEEDS="1 2 3 4"
EXP="snow_T-20_h1.00_30d"
# Same k_eff sampling as the pilot, so the two are directly comparable.
KEFF_OPTS="-keff 1 -keff_step0 1 -keff_t_interv 5.0e4 -keff_ksp_type cg -keff_pc_type gamg"

dry=0; sbatch_extra=()
while [ $# -gt 0 ]; do
    case "$1" in
        --dry-run) dry=1; shift ;;
        --)        shift; sbatch_extra=("$@"); break ;;
        -h|--help) sed -n '2,34p' "${BASH_SOURCE[0]}"; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

if [ "$dry" -eq 0 ] && ! command -v sbatch >/dev/null 2>&1; then
    echo "ERROR: sbatch not found -- run this on a cluster login node." >&2
    exit 1
fi

echo "=== REV check: L/R_ave = 64, phi = 0.325, T = -20 C, alpha_c = 1e-3 ==="
echo "  4 seeds, 724-805 grains, L = 3.2 mm, eps = 1 um, Nx = 4526"
echo "  61.5M DOF, ~769 cores/run, 3076 core-allocations total"
echo "  compare the ENSEMBLE mean against L/R_ave=40's 0.7570 (SEM 1.5%)"
echo ""
for s in $SEEDS; do
    geom=$(printf "$GEOM_STEM" "$s")
    cmd=(./scripts/HPC/submit_enceladus.sh "$geom" "$EXP" "rev64_s$s")
    [ "${#sbatch_extra[@]}" -gt 0 ] && cmd+=("${sbatch_extra[@]}")
    cmd+=(-- $KEFF_OPTS)
    if [ "$dry" -eq 1 ]; then printf '  $ '; printf '%q ' "${cmd[@]}"; echo
    else "${cmd[@]}"; fi
done
