#!/usr/bin/env bash
# =============================================================================
# make_rve_packings.sh — packings for the k_eff domain-size convergence study.
#
#   bash studies/keff_sintering/make_rve_packings.sh
#
# phi 0.325, R_ave = 50 um, L/R_ave {20, 30, 40, 56, 80} with {8, 6, 5, 4, 3}
# seeds (L = 1, 1.5, 2, 2.8, 4 mm).
#
# WHY AN ENSEMBLE PER SIZE. A larger periodic cell is a new realization, not a
# zoom-out of the same packing, and cropping nested windows from one big bed
# gives non-periodic windows (a seam again). So convergence is tested the
# statistical-RVE way (Kanit et al. 2003): the seed-mean k_eff curve must stop
# depending on L, and the seed scatter must shrink (~1/L in 2D). More seeds at
# small L, where the scatter is largest.
#
# NO HOMOGENEITY GATES -- void, density CV, half-domain asymmetry are all OFF
# at every size. The first build (2026-09-30) used the production void gate up
# to L/R 40 and a per-L gate above, and coordination at the band jumped
# 3.32 -> 3.47 exactly at that switch (the per-L gate admits larger voids, and
# at fixed porosity larger voids mean denser, better-connected solid
# elsewhere): the study would have confounded size with gate. Every
# homogeneity gate is a statement about the extreme of N quantities, and N
# grows with L, so no fixed gate is size-fair. Ungated, each size is an
# unfiltered sample of the same deposition process.
#   KEPT: the y-seam gate (it removes a generator artifact, not a
#   realization) and solid percolation (k_eff of a disconnected solid is not
#   the quantity studied; at phi 0.325 it rarely binds).
# The GATED production L/R 40 set (keff_LR40) against this ungated L/R 40 set
# then measures what the production recipe does to k_eff, separately.
# Porosity-salted streams; unique seed numbers (blocks 11-15, unused elsewhere).
#
# Output: inputs/packings/rve_phi0.325/phi0.325_Rave50um_LR<N>_seed<S>/
# =============================================================================
set -uo pipefail

PROJ="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT="$PROJ/inputs/packings/rve_phi0.325"
PY="$PROJ/venv_enceladus/bin/python"
JOBS="${JOBS:-8}"
LRS=(20 30 40 56 80)
NSEEDS=(8 6 5 4 3)
BLOCKS=(11 12 15 13 14)
MAX_SEEDS_PER_SLOT=20
mkdir -p "$OUT"

build_slot() {
    local LR="$1" b="$2" k="$3" n="$4"
    local Lx; Lx=$(awk -v r="$LR" 'BEGIN{printf "%.6e", r*50e-6}')
    local extra="--max-void-ratio 99 --max-density-cv 99 --max-asymmetry 99"
    local tries=0 seed=$((100 * b + k))
    while (( tries < MAX_SEEDS_PER_SLOT )); do
        local name="phi0.325_Rave50um_LR${LR}_seed${seed}"
        local dir="$OUT/$name"
        [[ -f "$dir/metadata.json" ]] && { echo "  exists  $name"; return 0; }
        mkdir -p "$dir"
        if "$PY" "$PROJ/preprocess/generate_packing_gravity.py" \
                --Lx "$Lx" --porosity 0.325 --mean-r 50e-6 --sigma-ln 0.5 \
                --periodic xy --seed "$seed" --band-per-mean-r 0.184 \
                --out "$dir" --no-periodic-subdir $extra > "$dir.build.log" 2>&1 \
           && [[ -f "$dir/metadata.json" ]]; then
            mv "$dir.build.log" "$dir/build.log"; echo "  OK      $name"; return 0
        fi
        echo "  reject  $name"
        rmdir "$dir" 2>/dev/null; mv "$dir.build.log" "$OUT/rejected_$name.build.log"
        tries=$((tries + 1))
        seed=$((100 * b + n + k + n * (tries - 1)))
    done
    echo "  FAILED  L/R $LR slot $k"; return 1
}
export -f build_slot
export OUT PY PROJ MAX_SEEDS_PER_SLOT

for i in "${!LRS[@]}"; do
    for k in $(seq 1 "${NSEEDS[$i]}"); do
        echo "${LRS[$i]} ${BLOCKS[$i]} $k ${NSEEDS[$i]}"
    done
done | xargs -P "$JOBS" -n 4 bash -c 'build_slot "$0" "$1" "$2" "$3"'

echo ""; echo "packings in $OUT: $(ls -d "$OUT"/phi*/ 2>/dev/null | wc -l)"
