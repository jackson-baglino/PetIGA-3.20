#!/usr/bin/env python3
"""plot_ssa.py — ice-air surface area vs time, in physical units.

Reads SSA_evo.dat from a run folder. Where a packing file is recoverable the
curve is also normalised by the analytic initial surface area of the grains, so
the y-axis reads as "fraction of the starting surface remaining".

TURNING THE MONITOR COLUMN INTO A LENGTH
----------------------------------------
Column 0 is not a surface area. monitoring.c writes ``sub_interf / eps`` where
Integration() (src/assembly.c) accumulates ``S[1] = phi^2 * (1-phi)^2``. The
conversion is exact for a resolved equilibrium interface, and worth doing
rather than plotting the raw proxy:

The 1D equilibrium profile of this double well is logistic,
``phi(x) = 1/(1 + exp(-x/eps))``, so ``dphi/dx = phi(1-phi)/eps``. Substituting
``dx = eps*dphi / (phi(1-phi))`` into the integral across one interface,

    int phi^2 (1-phi)^2 dx = eps * int_0^1 phi(1-phi) dphi = eps/6

Every unit length of interface therefore contributes ``eps/6`` to sub_interf,
and since the file already divides by eps,

    interface length  L = 6 * column_0                     [m]   in 2D
    interface area    A = 6 * column_0                     [m^2] in 3D

ASSUMPTIONS, AND WHEN THIS OVERSTATES THE AREA. The factor 6 holds for an
interface at its equilibrium profile and adequately resolved by the mesh. A
profile smeared wider than equilibrium -- an under-resolved interface, or the
first few steps of a run whose initial condition was written with a different
eps -- inflates the integral, so early-time values read high. It is a
diagnostic of the diffuse field, not a substitute for measuring a contour.

Under -axisym the integrals carry a 2*pi*r weight, so column 0 is already a
true 3D area and the same factor applies.

Usage:
    python3 plot_ssa.py --dir <run> [--save <png>] [--time-unit h]
"""
from __future__ import annotations

import argparse
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")          # before pyplot: these run headless on HPC
import matplotlib.pyplot as plt

from pplib import (SSA, TIME, auto_time_unit, grain_radii, in_time_unit,
                   load_ssa, opt_float, read_opts)

# int phi^2(1-phi)^2 dx = eps/6 across one equilibrium interface; see module docstring.
INTERFACE_FACTOR = 6.0


def analytic_initial_area(run_dir):
    """Initial ice-air surface from the packing, or None.

    2D: sum of circumferences, 2*pi*r. 3D: sum of sphere areas, 4*pi*r^2.
    Grains in contact overlap slightly at their necks, so this is a mild
    OVERESTIMATE of the true free surface -- it is a reference scale, not a
    measurement.
    """
    radii = grain_radii(run_dir)
    if radii is None or len(radii) == 0:
        return None, 0
    dim = int(opt_float(read_opts(run_dir), "-dim", 2) or 2)
    if dim >= 3:
        return float(np.sum(4.0 * math.pi * radii ** 2)), len(radii)
    return float(np.sum(2.0 * math.pi * radii)), len(radii)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dir", default=".", help="run directory (default: cwd)")
    p.add_argument("--save", default=None,
                   help="output PNG (default: <dir>/ssa.png)")
    p.add_argument("--time-unit", choices=["s", "min", "h", "d"], default=None,
                   help="x-axis unit (default: chosen from the run length)")
    p.add_argument("--loglog", action="store_true",
                   help="log-log axes, for reading a coarsening exponent")
    p.add_argument("--title", default=None)
    a = p.parse_args(argv)

    data = load_ssa(a.dir)
    if data is None:
        print(f"❌ No usable SSA_evo.dat in {a.dir}", file=sys.stderr)
        return 1

    t = data[:, TIME]
    area = data[:, SSA] * INTERFACE_FACTOR
    if not np.any(area > 0):
        print(f"❌ Interface density is zero throughout {a.dir} — nothing to plot.",
              file=sys.stderr)
        return 1

    unit = a.time_unit or auto_time_unit(float(t.max()) if len(t) else 0.0)
    tt = in_time_unit(t, unit)

    dim = int(opt_float(read_opts(a.dir), "-dim", 2) or 2)
    unit_label = "m$^2$" if dim >= 3 else "m"

    a0, n_grains = analytic_initial_area(a.dir)

    fig, ax1 = plt.subplots(figsize=(10, 6))
    plot = ax1.loglog if a.loglog else ax1.plot
    plot(tt, area, "-", color="#1b4f72", lw=1.8)

    ax1.set_xlabel(f"Time [{unit}]", fontsize=14)
    ax1.set_ylabel(f"Ice-air surface  [{unit_label}]", fontsize=14)
    ax1.tick_params(axis="both", which="major", labelsize=11)
    if a.loglog:
        ax1.grid(True, which="both", alpha=0.3)

    # Second axis as a pure rescaling of the first, so the two can never
    # disagree about the shape of the curve.
    ref, ref_note = (a0, f"{n_grains} grains from -grains_file") if a0 else \
                    (area[0], "t = 0 value (no packing file found)")
    ax2 = ax1.twinx()
    ax2.set_ylim(np.array(ax1.get_ylim()) / ref)
    if a.loglog:
        ax2.set_yscale("log")
    ax2.set_ylabel(f"Fraction of initial surface\n({ref_note})", fontsize=12)
    ax2.tick_params(axis="both", which="major", labelsize=11)

    ax1.set_title(a.title or
                  f"Ice-air surface evolution\n"
                  f"{area[0]:.4g} → {area[-1]:.4g} {unit_label}  "
                  f"({(area[-1] / area[0] - 1) * 100:+.2f}%)   "
                  f"over {tt[-1]:.3g} {unit}",
                  fontsize=15)
    fig.tight_layout()

    out = a.save or os.path.join(a.dir, "ssa.png")
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"  surface: {area[0]:.6g} → {area[-1]:.6g} {unit_label}")
    if a0:
        print(f"  analytic initial surface from packing: {a0:.6g} {unit_label} "
              f"({n_grains} grain rows); measured/analytic = {area[0] / a0:.3f}")
    print(f"  wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
