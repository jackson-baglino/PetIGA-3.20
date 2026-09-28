#!/usr/bin/env bash
# =============================================================================
# make_packings.sh — build the production packings for the porosity campaign.
#
#   bash studies/keff_sintering/make_packings.sh            # all of them
#   JOBS=4 bash studies/keff_sintering/make_packings.sh     # fewer in parallel
#
# phi {0.275, 0.325, 0.375, 0.425, 0.475} x 5 packings -- evenly spaced by
# 0.05 (user, 2026-09-28; the first set at 0.250/0.300/0.350/0.400 was
# dropped, 0.325 kept). R_ave = 50 um, L = 2 mm (L/R_ave = 40), periodic xy,
# every gate at its default (void and CV envelopes scaled with porosity, solid
# percolation, half-domain asymmetry, y-seam contact density >= 0.76 -- see
# studies/rve_anisotropy/README.md).
#
# EXCEPT phi 0.475: the solid percolation gate is OFF there. In 2D the solid
# stops percolating between 0.40 and 0.45; at 0.475 a probe had the ice fail
# to span x in 12 of 16 attempts. Forcing the gate would keep only the rare
# connected realizations -- a biased sample of a near-threshold medium. The
# 0.475 set is the TYPICAL microstructure instead, with percolation recorded
# in metadata.json, and serves to show where 2D stops being a snow analogue.
#
# SEEDS ARE UNIQUE. Each porosity has its own block of base seeds, 100*b + k,
# k = 1..5 (0.275 -> b=6, 0.325 -> 3, 0.375 -> 7, 0.425 -> 8, 0.475 -> 9;
# blocks 1, 2, 4, 5 belonged to the dropped set). If a base seed exhausts the generator's 128
# retries, the next UNUSED seed in that block is tried (100*i + 6, 7, ...), so
# no two packings share a seed number and none is reused across porosities.
# (The generator also salts its stream with the porosity, so they would be
# independent regardless; unique numbers make that obvious from the name.)
# Within a run the generator's retries use seed + 1000*attempt, which cannot
# collide with another packing's base seed because the last three digits
# differ.
#
# phi 0.325 was rebuilt here even though pilot_LR40/ exists: those predate the
# seam gate (seeds 1 and 4 have seams at 0.33 and 0.36 of the interior) and the
# porosity-salted stream, so the porosity series is built one way throughout.
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
BLOCKS=(6 3 7 8 9)
N_PER_PHI=5
MAX_SEEDS_PER_SLOT=20          # base seeds to try before giving up on a slot

mkdir -p "$OUT"

# One slot = one packing. Slot k of block i tries seeds 100*i + k, then the
# next unused numbers above 100*i + N_PER_PHI (disjoint per slot: 100*i + 5 +
# k, +5, ...), so parallel slots never pick the same seed.
build_slot() {
    local phi="$1" i="$2" k="$3"
    local extra=""
    [[ "$phi" == "0.475" ]] && extra="--no-percolation-gate"
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
