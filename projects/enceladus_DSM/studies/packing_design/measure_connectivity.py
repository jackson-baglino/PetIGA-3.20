#!/usr/bin/env python3
"""Measure solid/pore connectivity of 2D packings as a function of porosity.

Answers the design question "how much contact should there be between grains"
with data rather than principle, and settles whether a porosity exists at which
both the solid network and the usable pore network connect.

USAGE
    venv_enceladus/bin/python studies/packing_design/measure_connectivity.py \
        --porosities 0.30 0.35 0.40 0.45 0.50 --seeds 1 2 3 \
        --out studies/packing_design

Generation is ~2 s per packing, so the whole sweep is a couple of minutes on a
laptop. No solver involved.

WHAT IS MEASURED, AND WHY NOT THE OBVIOUS THING
-----------------------------------------------
The pore is measured TWICE:

  pore        the sharp complement of the solid. This is the intuitive
              measure and it is misleading twice over -- it uses
              4-connectivity on both phases (which severs both at a diagonal
              pinch) and it counts channels narrower than the diffuse band as
              open, when phi never reaches 0 inside them so they conduct like
              solid and block vapour.

  open pore   pore further than band/2 from any solid, labelled with
              8-connectivity: what the solver can actually use. This is also
              the measure that converges under raster refinement, where the
              sharp one does not.

Report the second. The first is kept only so the difference stays visible.
"""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "preprocess"))
import packing_lib as pl                                   # noqa: E402


def load(d: Path):
    m = json.loads((d / "metadata.json").read_text())
    rows = [l.split() for l in (d / "grains.dat").read_text().splitlines()
            if l and not l.startswith("#")]
    raw = np.array([[float(v) for v in r] for r in rows[1:]])   # row 0 = Lx Ly
    c, r = raw[:, :2], raw[:, 2]
    # grains.dat carries edge images too; the base tile is what tiles the torus
    k = ((c[:, 0] >= 0) & (c[:, 0] < m["Lx"])
         & (c[:, 1] >= 0) & (c[:, 1] < m["Ly"]))
    return m, c[k].copy(), r[k].copy()


def measure(d: Path, band_per_r: float, raster: int) -> dict:
    m, c, r = load(d)
    Lx, Ly = m["Lx"], m["Ly"]
    band = band_per_r * m["mean_r_m_requested"]
    solid = pl.rasterize(c, r, Lx, Ly, raster, raster, True, True)
    op = pl.open_pore(solid, band, Lx / raster, True, True)

    spx, spy, sfr, sn = pl.percolates(solid)
    ppx, ppy, pfr, pn = pl.percolates(~solid, diagonal=True)
    opx, opy, ofr, on = pl.percolates(op, diagonal=True)
    _, _, st = pl.descriptors(c, r, Lx, Ly, True, True, band=band)
    return {
        "name": d.name,
        "porosity": round(1.0 - float(solid.mean()), 5),
        "n_grains": len(c),
        "coordination_number": round(m["coordination_number"], 3),
        "coordination_at_band": round(st.get("coordination_at_band", float("nan")), 3),
        "solid_perc_x": int(spx), "solid_perc_y": int(spy),
        "pore_perc_x": int(ppx), "pore_perc_y": int(ppy),
        "pore_n_clusters": pn, "pore_largest_frac": round(pfr, 4),
        "open_perc_x": int(opx), "open_perc_y": int(opy),
        "open_n_clusters": on, "open_largest_frac": round(ofr, 4),
        "open_area_frac_of_pore": round(float(op.sum()) / float((~solid).sum()), 4),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--porosities", type=float, nargs="+",
                    default=[0.30, 0.35, 0.40, 0.45, 0.50])
    ap.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3])
    ap.add_argument("--Lx", type=float, default=1.0e-4)
    ap.add_argument("--mean-r", dest="mean_r", type=float, default=2.5e-6)
    ap.add_argument("--sigma-ln", dest="sigma_ln", type=float, default=0.5)
    ap.add_argument("--band-per-mean-r", dest="band_per_r", type=float, default=0.166,
                    help="9.2*eps / mean_r; 0.166 is eps = 45 nm at R_ave = 2.5 um")
    ap.add_argument("--raster", type=int, default=2048)
    ap.add_argument("--out", type=Path, default=Path(__file__).parent)
    a = ap.parse_args()

    packdir = a.out / "packings"
    packdir.mkdir(parents=True, exist_ok=True)
    rows = []
    for phi in a.porosities:
        for s in a.seeds:
            d = packdir / f"phi{phi:.2f}_seed{s}"
            if not (d / "metadata.json").is_file():
                # Gates OFF: this sweep is measuring the very tail the gates
                # exist to reject, so leaving them on would censor the answer.
                subprocess.run(
                    [sys.executable, str(ROOT / "preprocess" / "generate_packing_gravity.py"),
                     "--Lx", str(a.Lx), "--porosity", str(phi),
                     "--mean-r", str(a.mean_r), "--sigma-ln", str(a.sigma_ln),
                     "--periodic", "xy", "--seed", str(s), "--no-preview",
                     "--no-periodic-subdir", "--no-percolation-gate",
                     "--max-void-ratio", "99", "--max-density-cv", "99",
                     "--max-asymmetry", "99", "--porosity-tol", "99",
                     "--band-per-mean-r", str(a.band_per_r),
                     "--out", str(d)],
                    check=True, capture_output=True)
            rows.append(measure(d, a.band_per_r, a.raster))
            print(f"  {rows[-1]['name']:18s} phi={rows[-1]['porosity']:.4f} "
                  f"Z_band={rows[-1]['coordination_at_band']:5.2f} "
                  f"solid=({rows[-1]['solid_perc_x']},{rows[-1]['solid_perc_y']}) "
                  f"open=({rows[-1]['open_perc_x']},{rows[-1]['open_perc_y']}) "
                  f"largest={rows[-1]['open_largest_frac']:.3f}")

    csv_path = a.out / "connectivity.csv"
    with csv_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {csv_path}")

    both = [r for r in rows if r["solid_perc_x"] and r["solid_perc_y"]
            and r["open_perc_x"] and r["open_perc_y"]]
    print(f"packings where BOTH phases percolate in both directions: "
          f"{len(both)} / {len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
