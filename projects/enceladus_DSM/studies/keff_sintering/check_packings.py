#!/usr/bin/env python3
"""Summarize the production packings built by make_packings.sh.

    venv_enceladus/bin/python studies/keff_sintering/check_packings.py [dir]

One row per packing: porosity achieved, accepted attempt, y-seam contact ratio,
largest void, coordination at the band, solid percolation, and the contact
fabric (studies/rve_anisotropy/geometry.py) -- then per-porosity means, and a
check that every seed number is unique. Writes packings_summary.csv into the
packing directory.
"""
from __future__ import annotations

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PROJ = HERE.parents[1]
sys.path.insert(0, str(PROJ / "studies/rve_anisotropy"))
sys.path.insert(0, str(PROJ / "preprocess"))
from geometry import analyse                                 # noqa: E402


def main():
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else PROJ / "inputs/packings/keff_LR40"
    rows = []
    for d in sorted(root.glob("phi*/")):
        if not (d / "metadata.json").is_file():
            continue
        m = json.loads((d / "metadata.json").read_text())
        g, _ = analyse(d, "keff_LR40")
        rows.append(dict(
            name=d.name, phi_target=m["porosity_target"], phi=m["porosity_achieved"],
            base_seed=int(d.name.split("seed")[-1]), seed=m["seed"], attempt=m["attempt"],
            seam=m["seam_contact_ratio"], void_per_R=m["max_void_radius_per_mean_r"],
            density_cv=m["local_density_cv"], z_band=m["coordination_at_band"],
            perc_x=m["solid_percolates_x"], perc_y=m["solid_percolates_y"],
            n_grains=m["n_grains"], F_ratio=g["F_ratio"], F_xy=g["F_xy"]))
    if not rows:
        print(f"no packings under {root}")
        return 1

    with open(root / "packings_summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    print(f"{'packing':38s} {'phi':>6} {'att':>3} {'seam':>5} {'void/R':>6} {'z_band':>6} "
          f"{'perc':>4} {'Fxx/Fyy':>7} {'F_xy':>7}")
    for r in rows:
        print(f"{r['name']:38s} {r['phi']:6.4f} {r['attempt']:3d} {r['seam']:5.2f} "
              f"{r['void_per_R']:6.2f} {r['z_band']:6.2f} "
              f"{'xy' if r['perc_x'] and r['perc_y'] else 'NO':>4} "
              f"{r['F_ratio']:7.3f} {r['F_xy']:+7.4f}")

    by = defaultdict(list)
    for r in rows:
        by[r["phi_target"]].append(r)
    print(f"\n{'phi':>6} {'n':>2} {'seam min':>8} {'Fxx/Fyy':>12} {'F_xy mean':>9} {'z_band':>6}")
    for phi, v in sorted(by.items()):
        F = np.array([r["F_ratio"] for r in v]); X = np.array([r["F_xy"] for r in v])
        print(f"{phi:6.3f} {len(v):2d} {min(r['seam'] for r in v):8.2f} "
              f"{F.mean():6.3f}±{F.std(ddof=1) if len(F) > 1 else 0:.3f} {X.mean():+9.4f} "
              f"{np.mean([r['z_band'] for r in v]):6.2f}")

    bases = [r["base_seed"] for r in rows]
    seeds = [r["seed"] for r in rows]
    print(f"\nbase seeds unique: {len(set(bases)) == len(bases)}   "
          f"accepted seeds unique: {len(set(seeds)) == len(seeds)}   ({len(rows)} packings)")
    print(f"wrote {root / 'packings_summary.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
