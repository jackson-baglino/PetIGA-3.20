#!/usr/bin/env bash
# =============================================================================
# run_mesh_pair_fits.sh — exponent fits for the mesh-resolution pair
# (batches/mesh_pair.txt): tangent-contact Molaro grains at alpha_c = 1e-3,
# coarse (eps = 0.87 um) vs fine (eps = 0.24 um), against the -20 C data.
#
#     bash studies/sinter_exponent/analysis/run_mesh_pair_fits.sh
#
# Point RUN_DIRS at the run folders (space-separated) if they are not where
# this script looks by default. Outputs land in results/mesh_pair/.
#
# Three prefixes are produced ON PURPOSE, because the fitted exponent is a
# property of the curve PLUS the window:
#
#   meshpair_coarse  default protocol — the per-series resolution floor
#                    sqrt(12*eps/R) is enforced. On the coarse arm that floor
#                    is r/R0 = 0.346 and the whole run stops at 0.309, so the
#                    model has NO fittable window. That is the finding.
#   meshpair_full    floor overridden (--rmin 0), each series over its own
#                    full range. The headline numbers, caveated.
#   meshpair_shared  both series clipped to the overlap, r/R0 = 0.197-0.309.
#                    Apples-to-apples, but leaves the data only 3 points.
#
# NOT produced: an --anchor-neck variant. Anchoring exists to overlay a run
# that opened with a mature neck onto one that did not; these runs start at
# true tangent contact, so their t = 0 is real and re-zeroing the clock would
# throw away the one thing that makes them fittable with t0 = 0.
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
cd "$ROOT"

PY="${PY:-$ROOT/venv_enceladus/bin/python}"
BATCH="${BATCH:-$HOME/SimulationResults/HPC_results/enceladus_DSM/batch_2026-08-11__17.24.16_mesh_pair}"
RES="$(cd "$HERE/.." && pwd)/results/mesh_pair"
DATA20="inputs/validation/molaro2019_fig11_T-20.csv"

mkdir -p "$RES"

# Every arm of the batch that has a neck curve. The fine arm is included
# automatically once its folder is pulled down from the cluster.
if [ -z "${RUN_DIRS:-}" ]; then
    RUN_DIRS="$(ls -1d "$BATCH"/*_a1e-3_*/ 2>/dev/null || true)"
fi

csvs=() labels=()
for d in $RUN_DIRS; do
    [ -f "$d/neck_width.csv" ] || $PY postprocess/neck_width.py "$d" --axisym || continue
    csvs+=("$d/neck_width.csv")
    case "$d" in
        *coarse*) labels+=("model coarse eps0.87um (tangent)") ;;
        *fine*)   labels+=("model fine eps0.24um (tangent)") ;;
        *)        labels+=("$(basename "${d%/}")") ;;
    esac
done
[ "${#csvs[@]}" -gt 0 ] || { echo "no mesh_pair runs found under $BATCH"; exit 1; }
labels+=("Molaro 2019 (-20C)")
LB="$(IFS=,; echo "${labels[*]}")"

FIT="$PY postprocess/fit_neck_growth.py"

echo "== 1/3  default protocol (resolution floor enforced) =="
$FIT "${csvs[@]}" "$DATA20" --labels "$LB" \
    --out "$RES" --prefix meshpair_coarse --demmenie || true

echo "== 2/3  floor overridden, full range of each series =="
$FIT "${csvs[@]}" "$DATA20" --labels "$LB" --rmin 0 \
    --out "$RES" --prefix meshpair_full --demmenie

echo "== 3/3  shared window r/R0 = 0.197-0.309 =="
$FIT "${csvs[@]}" "$DATA20" --labels "$LB" --rmin 0.197 --rmax 0.309 \
    --out "$RES" --prefix meshpair_shared --demmenie

echo
echo "outputs -> $RES"
