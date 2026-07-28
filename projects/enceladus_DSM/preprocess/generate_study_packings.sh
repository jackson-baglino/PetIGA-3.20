#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# generate_study_packings.sh — the 25 packings for the snow-thermal k_eff study.
#
# 5 porosities x 5 seeds, every cell an INDEPENDENT realization. The previous
# study built higher-porosity configurations by adding grains to a lower-porosity
# base, which left the configurations nearly identical and made the porosity
# trend inseparable from the shared microstructure. Nothing is shared here.
#
# Porosity floor: jamming against r + gap/2 costs solid fraction, so
# (1-phi)*(1+gap/2r)^2 must stay under the 2D RCP limit ~0.84. At r = 45 um and
# gap = 4 um that gives phi >~ 0.23, which is why 0.25 is the bottom of the range
# (0.20 fails to jam).
# ---------------------------------------------------------------------------
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.."
PY=./venv_enceladus/bin/python3
OUT="${1:-inputs/packings}"

LX=2.0e-3          # Calonne 2011 REV floor for k_eff
MEAN_R=45e-6
SIGMA_LN=0.5
GAP=4.0e-6         # see generate_packing.py: 0 disconnects the pore space in 2D
POROSITIES=(0.250 0.2875 0.325 0.3625 0.400)
SEEDS=(1 2 3 4 5)

mkdir -p "$OUT"
fail=0
for phi in "${POROSITIES[@]}"; do
  for seed in "${SEEDS[@]}"; do
    name="phi${phi}_seed${seed}"
    if [[ -f "$OUT/$name/metadata.json" ]]; then
      echo "  skip $name (exists)"
      continue
    fi
    if ! $PY preprocess/generate_packing.py \
        --Lx "$LX" --porosity "$phi" --seed "$seed" \
        --mean_r_m "$MEAN_R" --sigma_ln "$SIGMA_LN" --contact-gap "$GAP" \
        --out "$OUT/$name" --raster 1024; then
      echo "  FAILED: $name"
      fail=$((fail+1))
    fi
  done
done

echo
echo "=== summary ==="
$PY - "$OUT" <<'PYEOF'
import json, sys, glob, os
root = sys.argv[1]
rows = []
for d in sorted(glob.glob(os.path.join(root, "phi*"))):
    f = os.path.join(d, "metadata.json")
    if not os.path.isfile(f):
        continue
    m = json.load(open(f))
    m["_dir"] = os.path.basename(d)
    rows.append(m)
print(f"{'packing':>22} {'N':>5} {'phi':>7} {'Z':>5} {'pore%':>6} {'throat p25[um]':>15}  status")
for m in rows:
    t = m["throat_gap_m"].get("p25", 0.0) * 1e6
    st = "OK" if not m["problems"] else m["problems"][0][:40]
    print(f"{m['_dir']:>22} {m['n_grains']:>5} {m['porosity_achieved']:>7.4f} "
          f"{m['coordination_number']:>5.2f} {m['pore_largest_cluster_frac']*100:>6.1f} "
          f"{t:>15.2f}  {st}")
print(f"\n{len(rows)} packings")
PYEOF

exit $(( fail > 0 ))
