#!/usr/bin/env python3
"""Per-step extremes of phi_i -- a resolution and stability check.

The solver prints, every accepted step,

    BOUNDS: phi_ice [min, max]  phi_air [min, max]

and this reads those back. It is the cheapest way to see whether a run stayed
resolved: a healthy phase field sits inside [0, 1] to within rounding, and any
real excursion is a symptom -- an under-resolved interface, a timestep too long
for the interface speed, or a boundary term pushing the wrong way.

Two panels, because the two questions need different scales:

  top     phi_min and phi_max against 0 and 1. The at-a-glance view.
  bottom  the EXCURSION magnitude, max(-phi_min, 0) and max(phi_max - 1, 0), on
          a log axis, with the -phase_lo / -phase_hi guard band drawn in. This
          is the panel that matters: it separates -1e-23 (rounding, fine) from
          -1e-4 (something is wrong) which the linear panel cannot.

Usage:  python3 postprocess/plot_phi_bounds.py --dir <run> [--save fig.png]
"""
import argparse
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pplib import read_opts, opt_float, load_ssa, TIME    # noqa: E402

CSV_NAME = "phi_bounds.csv"
_RE = re.compile(r"BOUNDS:\s*phi_ice\s*\[\s*([-\d.eE+]+),\s*([-\d.eE+]+)\s*\]"
                 r"\s*phi_air\s*\[\s*([-\d.eE+]+),\s*([-\d.eE+]+)\s*\]")


def read_bounds(run_dir):
    """(ice_min, ice_max, air_min, air_max) per accepted step, from outp.txt."""
    log = os.path.join(run_dir, "outp.txt")
    if not os.path.exists(log):
        raise SystemExit(
            f"no outp.txt in {run_dir} -- the BOUNDS lines come from the solver's\n"
            "monitor, so the run needs -pf_monitor 1 and its stdout tee'd there.")
    vals = [tuple(float(g) for g in m.groups())
            for m in (_RE.search(l) for l in open(log, errors="replace")) if m]
    if not vals:
        raise SystemExit(
            f"no BOUNDS lines in {log}. Runs made before 2026-09-21 printed them\n"
            "at %.4f, which this parser still reads; a run with -pf_monitor 0\n"
            "has none at all.")
    return np.array(vals).T


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True)
    ap.add_argument("--save", default=None)
    args = ap.parse_args()

    imin, imax, amin, amax = read_bounds(args.dir)
    n = imin.size
    step = np.arange(n)

    # Prefer real time if SSA_evo.dat lines up step-for-step with the monitor.
    x, xlabel = step, "accepted step"
    try:
        ssa = load_ssa(os.path.join(args.dir, "SSA_evo.dat"))
        if ssa.shape[0] == n:
            x, xlabel = ssa[:, TIME] / 86400.0, "time [days]"
    except Exception:
        pass

    opts = read_opts(args.dir)
    lo = opt_float(opts, "-phase_lo", -0.05)
    hi = opt_float(opts, "-phase_hi", 1.05)

    under = np.maximum(-imin, 0.0)
    over = np.maximum(imax - 1.0, 0.0)
    worst = max(under.max(), over.max())

    with open(os.path.join(args.dir, CSV_NAME), "w") as fh:
        fh.write("step,x,phi_ice_min,phi_ice_max,phi_air_min,phi_air_max,"
                 "undershoot,overshoot\n")
        for k in range(n):
            fh.write("%d,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n"
                     % (k, x[k], imin[k], imax[k], amin[k], amax[k],
                        under[k], over[k]))

    print(f"\n  PHI BOUNDS  ({os.path.basename(args.dir.rstrip('/'))})")
    print(f"    steps recorded     {n}")
    print(f"    phi_i min          {imin.min():+.6e}   (guard {lo:+.3f})")
    print(f"    phi_i max          {imax.max():+.6e}   (guard {hi:+.3f})")
    print(f"    worst excursion    {worst:.3e}", end="")
    print("   <- rounding only" if worst < 1e-12 else
          "   <- REAL, inspect the run" if worst > 1e-6 else "   <- small")

    if args.save:
        fig, ax = plt.subplots(2, 1, figsize=(7.4, 6.2), sharex=True,
                               constrained_layout=True)
        ax[0].plot(x, imax, lw=1.2, color="#1f77b4", label=r"$\phi_i$ max")
        ax[0].plot(x, imin, lw=1.2, color="#d62728", label=r"$\phi_i$ min")
        ax[0].axhline(0.0, color="k", lw=0.8)
        ax[0].axhline(1.0, color="k", lw=0.8)
        ax[0].set_ylabel(r"$\phi_i$")
        ax[0].set_title("Phase-field extremes per step")
        ax[0].legend(fontsize=8)
        ax[0].grid(alpha=0.3)

        ax[1].axhspan(0, max(abs(lo), hi - 1.0), color="green", alpha=0.07)
        ax[1].axhline(abs(lo), color="green", ls="--", lw=1.2,
                      label=f"guard band  |phase_lo| = {abs(lo):.3g}")
        ax[1].plot(x, np.maximum(under, 1e-20), lw=1.2, color="#d62728",
                   label=r"undershoot  $-\min(\phi_i,0)$")
        ax[1].plot(x, np.maximum(over, 1e-20), lw=1.2, color="#1f77b4",
                   label=r"overshoot  $\max(\phi_i-1,0)$")
        ax[1].set_yscale("log")
        ax[1].set_ylim(1e-20, max(1e-3, 2 * max(worst, abs(lo))))
        ax[1].set_ylabel("excursion outside [0, 1]")
        ax[1].set_xlabel(xlabel)
        ax[1].legend(fontsize=8, loc="upper left")
        ax[1].grid(alpha=0.3, which="both")
        fig.savefig(args.save, dpi=150)
        print(f"    -> {args.save}")


if __name__ == "__main__":
    main()
