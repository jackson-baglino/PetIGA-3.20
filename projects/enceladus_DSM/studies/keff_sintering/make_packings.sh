#!/usr/bin/env bash
# =============================================================================
# make_packings.sh — build the production packings for the porosity campaign.
#
#   bash studies/keff_sintering/make_packings.sh            # all of them
#   JOBS=4 bash studies/keff_sintering/make_packings.sh     # fewer in parallel
#
# phi {0.275, 0.325, 0.375, 0.425, 0.475} x 5 packings -- evenly spaced by
# 0.05 (user, 2026-09-28). R_ave = 50 um, L = 2 mm (L/R_ave = 40), periodic xy.
#
# GATES (option B, user 2026-10-01): only the y-seam gate (contact density
# >= 0.76 of the interior -- it removes a generator artifact, not a
# realization; studies/rve_anisotropy/README.md) and solid percolation. The
# homogeneity gates -- largest void, local density CV, half-domain asymmetry
# -- are OFF, exactly as in the convergence set (make_rve_packings.sh). Each
# of them tests an EXTREME over the domain, and at L/R 40 they no longer
# remove bad draws, they remove ordinary ones with one large void: the first
# (gated) build had z_band 3.32 +- 0.07 against 3.48 +- 0.08 ungated at the
# same size, while ungated z_band is converged from L/R 30 (3.47-3.53).
# Ungated, every production packing is an unfiltered sample of the
# deposition process, and the convergence study validates this recipe
# directly. The gated build (seeds 301-905) is kept in
# inputs/packings/keff_LR40_gated/ for provenance: the 3a shakedown ran on it,
# and four of its 0.325 packings are run in batch_rve.txt to measure what
# the gates did to k_eff.
#
# EXCEPT phi 0.475: the solid percolation gate is OFF there too. In 2D the
# solid stops percolating between 0.40 and 0.45; at 0.475 a probe had the ice
# fail to span x in 12 of 16 attempts. Forcing the gate would keep only the
# rare connected realizations -- a biased sample of a near-threshold medium.
# The 0.475 set is the TYPICAL microstructure instead, with percolation
# recorded in metadata.json, and serves to show where 2D stops being a snow
# analogue.
#
# SEEDS ARE UNIQUE. Each porosity has its own block of base seeds, 100*b + k,
# k = 1..5 (0.275 -> b=16, 0.325 -> 17, 0.375 -> 18, 0.425 -> 19,
# 0.475 -> 20). Blocks 1-9 are the dropped and the gated sets, 11-15 the
# convergence set. If a base seed exhausts the generator's 128 retries, the
# next UNUSED seed in that block is tried (100*b + 6, 7, ...), so no two
# packings share a seed number and none is reused across porosities. (The
# generator also salts its stream with the porosity, so they would be
# independent regardless; unique numbers make that obvious from the name.)
#
# Output: inputs/packings/keff_LR40/phi<X>_Rave50um_LR40_seed<N>/ with
# grains.dat, metadata.json, preview.png; the generator log goes beside it as
# build.log. Local, a few minutes at most per packing.
# =============================================================================
set -uo pipefail

PROJ="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT="$PROJ/inputs/packings/keff_LR40"
PY="$PROJ/venv_enceladus/bin/python"
JOBS="${JOBS:-10}"
PHIS=(0.275 0.325 0.375 0.425 0.475)
BLOCKS=(16 17 18 19 20)
N_PER_PHI=5
MAX_SEEDS_PER_SLOT=20          # base seeds to try before giving up on a slot

mkdir -p "$OUT"

# One slot = one packing. Slot k of block i tries seeds 100*i + k, then the
# next unused numbers above 100*i + N_PER_PHI (disjoint per slot: 100*i + 5 +
# k, +5, ...), so parallel slots never pick the same seed.
build_slot() {
    local phi="$1" i="$2" k="$3"
    local extra="--max-void-ratio 99 --max-density-cv 99 --max-asymmetry 99"
    [[ "$phi" == "0.475" ]] && extra="$extra --no-percolation-gate"
    local tries=0 seed=$((100 * i + k))
    while (( tries < MAX_SEEDS_PER_SLOT )); do
        local name="phi${phi}_Rave50um_LR40_seed${seed}"
        local dir="$OUT/$name"
        if [[ -f "$dir/metadata.json" ]]; then
            echo "  exists  $name"; return 0
        fi
        mkdir -p "$dir"
        if "$PY" "$PROJ/preprocess/generate_packing_gravity.py" \
                --Lx 2e-3 --porosity "$phi" --mean-r 50e-6 --sigma-ln 0.5 \
                --periodic xy --seed "$seed" --band-per-mean-r 0.184 \
                --out "$dir" --no-periodic-subdir $extra > "$dir.build.log" 2>&1 \
           && [[ -f "$dir/metadata.json" ]]; then
            mv "$dir.build.log" "$dir/build.log"
            echo "  OK      $name"; return 0
        fi
        echo "  reject  $name (all retries failed; see $name.build.log)"
        rmdir "$dir" 2>/dev/null
        mv "$dir.build.log" "$OUT/rejected_$name.build.log"
        tries=$((tries + 1))
        seed=$((100 * i + N_PER_PHI + k + N_PER_PHI * (tries - 1)))
    done
    echo "  FAILED  phi $phi slot $k after $MAX_SEEDS_PER_SLOT seeds"; return 1
}
export -f build_slot
export OUT PY PROJ N_PER_PHI MAX_SEEDS_PER_SLOT

for i in "${!PHIS[@]}"; do
    for k in $(seq 1 "$N_PER_PHI"); do
        echo "${PHIS[$i]} ${BLOCKS[$i]} $k"
    done
done | xargs -P "$JOBS" -n 3 bash -c 'build_slot "$0" "$1" "$2"'

echo ""
echo "packings in $OUT:"
ls -d "$OUT"/phi*/ 2>/dev/null | wc -l
