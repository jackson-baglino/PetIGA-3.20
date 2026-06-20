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

    # SSA_evo.dat has one row per accepted TS step (column 2 = TIME [s]),
    # written by the same Monitor() callback that triggers sol_*.dat output
    # at -outp 1 (the project default) -- so row i's time matches sol_i.dat
    # exactly. Falls back to step index if SSA_evo.dat isn't found (e.g.
    # -outp > 1, breaking the 1:1 mapping) so the script still runs.
    ssa_path = os.path.join(args.dir, "SSA_evo.dat")
    times = None
    if os.path.isfile(ssa_path):
        ssa = np.loadtxt(ssa_path)
        if ssa.ndim == 1:
            ssa = ssa.reshape(1, -1)
        if len(ssa) == len(sol_files):
            times = ssa[:, 2]
        else:
            print(f"Warning: SSA_evo.dat has {len(ssa)} rows but {len(sol_files)} "
                  f"sol_*.dat files -- falling back to step index for the time axis "
                  f"(check -outp; the 1:1 mapping assumes -outp 1).", file=sys.stderr)

    # Build the left/right mask once from the dense coordinate grid (same
    # grid is reused for every snapshot -- geometry doesn't move).
    C0, _ = nrb(u, v, fields=PetIGA().read_vec(sol_files[0], nrb))
    left_mask = C0[..., 0] < split
    right_mask = ~left_mask

    lefts, rights = [], []
    for f in sol_files:
        sol = PetIGA().read_vec(f, nrb)
        _, F = nrb(u, v, fields=sol)
        ice = F[..., 0]
        lefts.append(float(ice[left_mask].sum()))
        rights.append(float(ice[right_mask].sum()))

    lefts = np.array(lefts)
    rights = np.array(rights)
    totals = lefts + rights
    steps = np.arange(len(sol_files))
    if times is None:
        times = steps.astype(float)

    print(f"{'step':>6} {'time [s]':>12} {'left':>12} {'right':>12} {'total':>12}")
    for s, ti, l, r, t in zip(steps, times, lefts, rights, totals):
        print(f"{s:6d} {ti:12.4e} {l:12.4f} {r:12.4f} {t:12.4f}")

    if len(lefts) > 1:
        dl = (lefts[-1] - lefts[0]) / lefts[0] * 100.0
        dr = (rights[-1] - rights[0]) / rights[0] * 100.0
        dt = (totals[-1] - totals[0]) / totals[0] * 100.0
        print(f"\nLeft region change:  {dl:+.2f}%")
        print(f"Right region change: {dr:+.2f}%")
        print(f"Total change:        {dt:+.4f}%  (should be ~0, conservation check)")

    if args.save_csv:
        with open(args.save_csv, "w") as fh:
            fh.write("step,time,left,right,total\n")
            for s, ti, l, r, t in zip(steps, times, lefts, rights, totals):
                fh.write(f"{s},{ti},{l},{r},{t}\n")
        print(f"\nWrote {args.save_csv}")

    if args.save_fig or not args.no_show:
        import matplotlib.pyplot as plt

        # Auto-pick a sensible time unit so the time-axis subplot isn't
        # plotted in seconds for a run spanning weeks/months.
        t_max = float(times.max())
        if t_max >= 86400:
            time_div, time_unit = 86400.0, "days"
        elif t_max >= 3600:
            time_div, time_unit = 3600.0, "hours"
        else:
            time_div, time_unit = 1.0, "s"
        times_scaled = times / time_div

        fig, (ax_step, ax_time) = plt.subplots(1, 2, figsize=(12, 5))
        for ax, x, xlabel in ((ax_step, steps, "step"),
                              (ax_time, times_scaled, f"simulation time [{time_unit}]")):
            ax.plot(x, lefts, label=f"left (x<{split:.2e})")
            ax.plot(x, rights, label=f"right (x>{split:.2e})")
            ax.plot(x, totals, label="total", linestyle="--", color="gray")
            ax.set_xlabel(xlabel)
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
