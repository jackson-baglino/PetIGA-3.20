#!/usr/bin/env python3
"""
track_grain_mass.py — Per-region ice mass over time, for Ostwald-ripening runs.

Why this exists
----------------
outp.txt's TOT_ICE and SSA_evo.dat's I-A INTERF are aggregate, whole-domain
integrals. In a two-grain Ostwald-ripening setup, mass moving from the small
grain to the large one leaves TOT_ICE almost exactly flat (conservation) and
I-A INTERF only weakly informative (it mixes both grains' perimeters, which
move in opposite directions) -- so watching either one can make substantial,
real ripening look like "nothing is happening." This script integrates
phi_ice separately over each region (default: left half / right half of the
domain, matching the two_ice_grains_boundary IC convention -- small grain at
x=0, large grain at x=Lx) across every sol_*.dat snapshot, so the actual
mass transfer between grains is visible directly.

Usage
-----
  python track_grain_mass.py --dir /path/to/run
  python track_grain_mass.py --dir /path/to/run --split 3.0e-5   # custom x split
  python track_grain_mass.py --dir /path/to/run --save-csv mass.csv --save-fig mass.png
"""

import argparse
import glob
import os
import sys

import numpy as np
from igakit.io import PetIGA

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_permafrost_highres import _dense_uv


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dir", required=True, help="Run output directory (contains igasol.dat, sol_*.dat)")
    p.add_argument("--split", type=float, default=None,
                   help="x-coordinate splitting left/right regions (default: domain midpoint)")
    p.add_argument("--n-per-elem", type=int, default=3, help="Dense-sample points per element (default 3)")
    p.add_argument("--save-csv", default=None, help="Write time, left, right, total to this CSV path")
    p.add_argument("--save-fig", default=None, help="Save a plot to this path")
    p.add_argument("--no-show", action="store_true", help="Don't open an interactive plot window")
    args = p.parse_args()

    igasol = os.path.join(args.dir, "igasol.dat")
    if not os.path.isfile(igasol):
        sys.exit(f"Not found: {igasol}")
    nrb = PetIGA().read(igasol)
    u, v = _dense_uv(nrb, args.n_per_elem)

    Lx = float(nrb.points[..., 0].max())
    split = args.split if args.split is not None else Lx / 2.0

    sol_files = sorted(glob.glob(os.path.join(args.dir, "sol_*.dat")))
    if not sol_files:
        sys.exit(f"No sol_*.dat files found in {args.dir}")

    # Build the left/right mask once from the dense coordinate grid (same
    # grid is reused for every snapshot -- geometry doesn't move).
    C0, _ = nrb(u, v, fields=PetIGA().read_vec(sol_files[0], nrb))
    left_mask = C0[..., 0] < split
    right_mask = ~left_mask

    times, lefts, rights = [], [], []
    for f in sol_files:
        sol = PetIGA().read_vec(f, nrb)
        _, F = nrb(u, v, fields=sol)
        ice = F[..., 0]
        lefts.append(float(ice[left_mask].sum()))
        rights.append(float(ice[right_mask].sum()))
        times.append(len(times))  # placeholder; replaced below if outp.txt found

    lefts = np.array(lefts)
    rights = np.array(rights)
    totals = lefts + rights
    steps = np.arange(len(sol_files))

    print(f"{'step':>6} {'left':>12} {'right':>12} {'total':>12}")
    for s, l, r, t in zip(steps, lefts, rights, totals):
        print(f"{s:6d} {l:12.4f} {r:12.4f} {t:12.4f}")

    if len(lefts) > 1:
        dl = (lefts[-1] - lefts[0]) / lefts[0] * 100.0
        dr = (rights[-1] - rights[0]) / rights[0] * 100.0
        dt = (totals[-1] - totals[0]) / totals[0] * 100.0
        print(f"\nLeft region change:  {dl:+.2f}%")
        print(f"Right region change: {dr:+.2f}%")
        print(f"Total change:        {dt:+.4f}%  (should be ~0, conservation check)")

    if args.save_csv:
        with open(args.save_csv, "w") as fh:
            fh.write("step,left,right,total\n")
            for s, l, r, t in zip(steps, lefts, rights, totals):
                fh.write(f"{s},{l},{r},{t}\n")
        print(f"\nWrote {args.save_csv}")

    if args.save_fig or not args.no_show:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.plot(steps, lefts, label=f"left (x<{split:.2e})")
        ax.plot(steps, rights, label=f"right (x>{split:.2e})")
        ax.plot(steps, totals, label="total", linestyle="--", color="gray")
        ax.set_xlabel("step")
        ax.set_ylabel("integrated phi_ice")
        ax.legend()
        fig.tight_layout()
        if args.save_fig:
            fig.savefig(args.save_fig, dpi=140)
            print(f"Wrote {args.save_fig}")
        if not args.no_show:
            plt.show()


if __name__ == "__main__":
    main()
