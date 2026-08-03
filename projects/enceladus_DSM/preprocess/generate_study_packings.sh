#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# generate_study_packings.sh — the packings for the snow-thermal k_eff study.
#
# 5 porosities x 5 kept seeds. Each porosity generates a batch of candidates and
# keeps those that RETAIN THE MOST OPEN PORE, via select_packing_seeds.py.
#
# THE RULES THESE PACKINGS FOLLOW (all adopted 2026-08-01..03; the previous
# version of this script used none of them)
# ---------------------------------------------------------------------------
# 1. GRAINS TOUCH.  --contact-gap 0, exact tangency. Each grain's own field
#    phi_k = 0.5 - 0.5*tanh(0.5*(r_k - R_k)/eps) crosses 0.5 at exactly r_k=R_k,
#    so "in contact" is eps-independent and needs no calibration. The previous
#    +4 um gap held the grains APART to keep the pore percolating. That gap was
#    narrower than the diffuse band, so conduction at t=0 ran through the
#    interface representation rather than the microstructure -- k_eff differed by
#    31-39% between eps = 1 and 2 um on an identical packing. It also corrupted
#    the trend being measured: k_eff rises through a run because connectivity
#    grows, and starting disconnected makes part of that rise a numerical gap
#    closing rather than real neck growth.
#
# 2. GRAINS MAY HANG PAST THE EDGE.  Centres span [-R, L+R]; a grain need only
#    be PARTLY inside. These cells represent a window cut from a much larger
#    pack, and dropping the grains that straddle an edge would sever the
#    conduction paths leaving the domain and bias k_eff low. Implemented by
#    emitting each boundary grain's periodic image explicitly into grains.dat,
#    so the solver needs no wrapping and can run with zero-flux walls.
#
# 3. 100 um MEAN RADIUS, matching the Molaro et al. (2019) grain experiments.
#    At 2 mm that is roughly 65-90 grains per cell.
#
# 4. PORE PERCOLATION IS NOT EXPECTED, and is not a failure. In 2D a spanning
#    contact network necessarily cuts the void into cells; only the choice of
#    configuration can mitigate it, which is what the seed ranking does. For
#    Enceladus the vapour budget is small at tiger-stripe conditions anyway.
#    This is a genuine 2D limitation and belongs in the write-up.
#
# Each packing directory gets grains.dat, metadata.json and preview.png. The
# preview draws the domain box and the overhanging grains, so a configuration
# can be checked by eye without running the solver.
#
# Usage:  ./generate_study_packings.sh [out_dir] [Lx]
# ---------------------------------------------------------------------------
set -uo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.."
PY=./venv_enceladus/bin/python3
OUT="${1:-inputs/packings}"
LX="${2:-2.0e-3}"

MEAN_R=100e-6
SIGMA_LN=0.5
GAP=0.0
CANDIDATES="1-12"
KEEP=5
POROSITIES=(0.250 0.2875 0.325 0.3625 0.400)

mkdir -p "$OUT"
echo "=== study packings: ${#POROSITIES[@]} porosities x $KEEP seeds, Lx=$LX, mean r=$MEAN_R ==="

for phi in "${POROSITIES[@]}"; do
    echo
    echo "--- porosity $phi ---"
    $PY preprocess/select_packing_seeds.py \
        --porosity "$phi" --Lx "$LX" --mean-r "$MEAN_R" --sigma-ln "$SIGMA_LN" \
        --contact-gap="$GAP" --seeds "$CANDIDATES" --keep "$KEEP" \
        --out "$OUT" --raster 1024
done

echo
echo "=== summary ==="
$PY - "$OUT" <<'PYEOF'
import json, sys, glob, os
root = sys.argv[1]
rows = []
for d in sorted(glob.glob(os.path.join(root, "phi*"))):
    f = os.path.join(d, "metadata.json")
    if os.path.isfile(f):
        m = json.load(open(f)); m["_dir"] = os.path.basename(d); rows.append(m)
print(f"{'packing':>22} {'N':>4} {'rows':>5} {'img':>4} {'phi':>7} {'Z':>5} "
      f"{'solid%':>7} {'pore%':>6} {'p50[um]':>8}")
for m in rows:
    print(f"{m['_dir']:>22} {m['n_grains']:>4} {m.get('n_grain_rows_emitted','-'):>5} "
          f"{m.get('n_edge_images','-'):>4} {m['porosity_achieved']:>7.4f} "
          f"{m['coordination_number']:>5.2f} "
          f"{m['solid_largest_cluster_frac']*100:>6.1f}% "
          f"{m['pore_largest_cluster_frac']*100:>5.1f}% "
          f"{m['throat_gap_m']['p50']*1e6:>8.1f}")
print(f"\n{len(rows)} packings in {root}")
PYEOF
