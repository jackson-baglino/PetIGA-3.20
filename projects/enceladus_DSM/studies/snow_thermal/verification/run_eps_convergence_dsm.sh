#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run_eps_convergence_dsm.sh — does the eps dependence of k_eff survive DSM?
#
# THE QUESTION. The t=0 probe (run_t0_probe.sh) found k_eff differing by 31-39%
# between safety 0.5 and safety 1.0 on an identical microstructure. But t=0 is
# the worst case by construction: the packing is jammed with a 4 um design
# contact gap and the solid barely percolates, so every conduction path runs
# through gaps narrower than the diffuse band itself. At t=0 the band IS the
# conduction path. That says nothing about k_eff once sintering has built solid
# necks, which is what the study actually reports.
#
# So: run the SAME initial condition at both eps for 6 hours of DSM, recording
# k_eff at EVERY accepted step, and watch the ratio
#
#     R(t) = k_iso(t; safety 1.0) / k_iso(t; safety 0.5)
#
# Two outcomes, and they mean different things:
#
#   R(t) -> 1 as necks form
#       The eps dependence is a START-UP TRANSIENT, not a property of the
#       converged physics. safety 1.0 is then usable, and -- more interesting --
#       it means our INITIALISATION is introducing an interface-width dependence
#       that the model does not otherwise have. The union-SDF tanh profile we lay
#       down is not an equilibrium profile of the double well: it is a geometric
#       construction that the Allen-Cahn dynamics must first relax. The relaxation
#       is O(eps)-dependent, so the first few steps are eps-contaminated by
#       construction. That is a fixable initialisation problem, not a physics
#       limit -- see the literature note in the README.
#
#   R(t) stays near 1.35
#       The dependence is structural. safety 1.0 is unusable, and an eps ladder
#       (0.5 -> 0.25 -> 0.125) is needed before any k_eff is trusted.
#
# WHY THIS DOMAIN. 0.5 mm, 28 grains, phi = 0.327, Z = 3.0 -- the calibration
# packing. Small enough that both arms are cheap (Nx 708 and 354, 19 and 5 cores)
# and that every-step field output is affordable, while carrying the same 4 um
# contact gap and skewed throat distribution as the 2 mm study packings.
#
# Both arms go through scripts/Studio/run_enceladus.sh, per CLAUDE.md: that
# stages the .opts and src/ into the output folder and generates the plots.
# Results land in $RESULTS_BASE, NOT in /tmp.
#
# Usage:  ./run_eps_convergence_dsm.sh [tag]
# ---------------------------------------------------------------------------
set -uo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/../../.."
PROJECT_ROOT="$PWD"
RUNNER="$PROJECT_ROOT/scripts/Studio/run_enceladus.sh"

TAG="${1:-epsconv}"
EXP="snow_T-20_6hr"
GEOMS=(phi0.325_0.5mm_T-20_s05 phi0.325_0.5mm_T-20_s10)

[[ -x "$RUNNER" ]] || { echo "ERROR: $RUNNER missing" >&2; exit 1; }

# -keff_freq 1 : a k_eff sample on every accepted step. This is the whole point
#                of the in-line diagnostic -- the old snapshot-driven path could
#                not resolve a start-up transient at all.
# -outp 1      : a field snapshot every step too, so the neck geometry behind any
#                k_eff change can be inspected rather than inferred. Affordable
#                only because the domain is small (~1.5 MB and ~0.4 MB a step).
KEFF_FLAGS=(-keff 1 -keff_freq 1 -outp 1)

echo "=== eps convergence through DSM: 0.5 mm, 6 h, k_eff every step ==="
for g in "${GEOMS[@]}"; do
    [[ -f "$PROJECT_ROOT/inputs/geometry/${g}.opts" ]] || {
        echo "  MISSING inputs/geometry/${g}.opts" >&2; exit 1; }
done

for g in "${GEOMS[@]}"; do
    echo
    echo "--- $g ---"
    "$RUNNER" "$g" "$EXP" "$TAG" -- "${KEFF_FLAGS[@]}"
    rc=$?
    [[ $rc -eq 0 ]] || echo "  run_enceladus.sh returned $rc for $g" >&2
done

echo
echo "Both arms done. Compare with:"
echo "  venv_enceladus/bin/python3 studies/snow_thermal/verification/plot_eps_convergence_dsm.py \\"
echo "      --s05 <results>/phi0.325_0.5mm_T-20_s05/<stamp>_${EXP}_${TAG} \\"
echo "      --s10 <results>/phi0.325_0.5mm_T-20_s10/<stamp>_${EXP}_${TAG}"
