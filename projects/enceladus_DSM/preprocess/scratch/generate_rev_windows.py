#!/usr/bin/env python3
"""Cut nested, centred sub-windows out of one master packing, for the REV study.

THE QUESTION
------------
How large must the cell be before k_eff stops depending on its size? The answer
is the representative elementary volume, and it has to be settled before any
absolute k_eff is quoted -- and before the porosity sweep, since a k_eff measured
below the REV is a property of the window, not of the material.

WHY NESTED WINDOWS OF ONE MASTER, rather than an independent packing per size
----------------------------------------------------------------------------
Independent packings differ in microstructure as well as in size, so k_eff(L)
carries realisation scatter on top of the size trend -- and at the small end that
scatter is large (two 1 mm seeds at the same target porosity differed by 48% in
k_iso in this project's own test runs). Nested windows share one microstructure,
so the only thing varying is how much of it is being averaged. The size trend is
then clean, at the cost of describing one realisation; the seed-to-seed spread is
a separate question, answered by repeating this on a second master.

WHY THIS IS NOW STRAIGHTFORWARD
-------------------------------
Cropping a window used to require re-periodising it -- grains straddling the cut
had to be deleted (biasing porosity and severing conduction paths) or wrapped
(creating overlaps at the seam). Neither is needed now:

  * grains may hang past the domain edge, centres spanning [-R, L+R], so a window
    simply keeps every grain whose disk overlaps it, exactly as if the window had
    been cut out of a larger pack -- which it has;
  * the microstructure does not need to be periodic. Periodic BCs on a
    non-periodic window give an apparent conductivity bracketed by the Dirichlet
    and Neumann results and converging to the effective value faster than either
    (Calonne et al. 2011 apply exactly this to tomographic snow images), and the
    corrector builds its own periodic mesh regardless of the solver's walls.

So a window is just: keep the overlapping grains, shift to window coordinates,
re-measure porosity and connectivity on the window itself.

WHAT TO WATCH IN THE RESULT
---------------------------
Porosity is re-measured per window and will drift by a few percent at the small
end -- that is real, not an error, but k_eff responds far more strongly to
porosity than to L, so read the two columns together. A k_eff(L) trend that
tracks the porosity drift is not a size effect.

Usage:
    python3 generate_rev_windows.py --master inputs/packings_rev/master_... \\
        --sizes 0.5 1.0 1.5 2.0 2.5 3.0 --out inputs/packings_rev
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import generate_packing as gp


def read_grains(path: Path):
    rows = [l.split() for l in path.read_text().splitlines()
            if l.strip() and not l.startswith("#")]
    data = [[float(x) for x in r] for r in rows if len(r) == 3]
    hdr = [[float(x) for x in r] for r in rows if len(r) == 2]
    a = np.array(data)
    Lx, Ly = (hdr[0] if hdr else [None, None])
    return a[:, :2], a[:, 2], Lx, Ly


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--master", type=Path, required=True)
    ap.add_argument("--sizes", type=float, nargs="+", required=True,
                    help="window side lengths in mm")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--raster", type=int, default=1024)
    args = ap.parse_args(argv)

    centres, radii, Lm, _ = read_grains(args.master / "grains.dat")
    mmeta = json.loads((args.master / "metadata.json").read_text())
    print(f"master: {args.master.name}  L = {Lm*1e3:.2f} mm  "
          f"{len(radii)} grain rows  target phi = {mmeta['porosity_target']:.4f}")

    # Windows are cut from a master that was jammed periodically, but each window
    # is itself a plain box: no minimum image when measuring its descriptors.
    gp.PERIODIC[0] = False

    print(f"\n{'window':>18} {'L[mm]':>7} {'grains':>7} {'phi':>7} {'Z':>5} "
          f"{'solid%':>7} {'pore%':>6} {'p50[um]':>8}")
    for Lmm in args.sizes:
        L = Lmm * 1e-3
        if L > Lm + 1e-12:
            print(f"  skip {Lmm} mm: larger than the master ({Lm*1e3:.2f} mm)")
            continue
        c0 = 0.5 * (Lm - L)                       # centred window origin

        # Keep every grain whose DISK overlaps the window, not just those whose
        # centre is inside. Dropping the straddlers would cut the conduction
        # paths that leave the window, which is what k_eff is most sensitive to.
        keep = ((centres[:, 0] > c0 - radii) & (centres[:, 0] < c0 + L + radii) &
                (centres[:, 1] > c0 - radii) & (centres[:, 1] < c0 + L + radii))
        wc = centres[keep] - c0
        wr = radii[keep]

        ny = int(round(args.raster))
        solid = gp.rasterize(wc, wr, L, L, args.raster, ny)
        phi = 1.0 - solid.mean()
        sx, sy, s_frac, s_n = gp.percolates(solid)
        px, py, p_frac, p_n = gp.percolates(~solid)
        coord, n_contact, throats, _ = gp.descriptors(wc, wr, L, L, 0.0)

        name = f"rev_L{Lmm:g}mm"
        d = args.out / name
        meta = {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "generator": "preprocess/generate_rev_windows.py",
            "algorithm": f"centred window cut from {args.master.name}",
            "master": str(args.master),
            "master_Lx": float(Lm),
            "periodic": False,
            "Lx": float(L), "Ly": float(L),
            "seed": mmeta.get("seed"),
            "n_grains": int(keep.sum()),
            # Area fraction measured ON THE WINDOW. sum(pi r^2) would overcount,
            # because grains that only partly overlap contribute only partly.
            "porosity_target": mmeta["porosity_target"],
            "porosity_achieved": float(phi),
            "porosity_raster": float(phi),
            "raster_nx": args.raster, "raster_ny": ny,
            "contact_gap_m": 0.0,
            "mean_r_m_realized": float(wr.mean()),
            "sigma_ln": mmeta.get("sigma_ln"),
            "overlap_free": True,
            "solid_percolates_x": bool(sx), "solid_percolates_y": bool(sy),
            "solid_largest_cluster_frac": float(s_frac), "solid_n_clusters": int(s_n),
            "pore_percolates_x": bool(px), "pore_percolates_y": bool(py),
            "pore_largest_cluster_frac": float(p_frac), "pore_n_clusters": int(p_n),
            "coordination_number": float(coord), "n_contacts": int(n_contact),
            "throat_gap_m": throats,
            "problems": [],
            "warnings": [f"window cut from {args.master.name}; not periodic, "
                         f"and the pore is not expected to percolate"],
        }
        gp.write_outputs(str(d), wc, wr, L, L, meta, preview=True)
        print(f"{name:>18} {Lmm:7.2f} {int(keep.sum()):7d} {phi:7.4f} {coord:5.2f} "
              f"{s_frac*100:6.1f}% {p_frac*100:5.1f}% {throats['p50']*1e6:8.1f}")

    print(f"\nwindows written to {args.out}")
    print("Read k_eff(L) together with the porosity column: k_eff responds much "
          "more strongly to porosity\nthan to L, so a trend that tracks porosity "
          "drift is not a size effect.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
