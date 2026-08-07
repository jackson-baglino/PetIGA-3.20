#!/usr/bin/env python3
"""plot_porosity.py — porosity and interface density vs time, on shared axes.

Reads SSA_evo.dat from a run folder and plots porosity (left axis) against the
ice-air interface density (right axis), so the two halves of the sintering
story -- pore space closing, surface area falling -- are legible in one figure.

POROSITY IS EXACT HERE, NOT INFERRED FROM THE DOMAIN SIZE
---------------------------------------------------------
Integration() in src/assembly.c integrates both phase indicators over the same
quadrature:

    S[0] = phi        ->  column ICE   (tot_ice)
    S[2] = 1 - phi    ->  column AIR   (tot_air)

They sum to the domain measure identically, so

    porosity = tot_air / (tot_ice + tot_air)

needs no -Lx/-Ly and no -dim. It is also correct under -axisym, where both
integrals carry the same 2*pi*r weight and the ratio cancels it; computing
porosity as tot_ice / (Lx*Ly) instead -- which the previous version of this
script did -- is wrong for axisymmetric runs and wrong by a factor of the
domain area for everything else.

Older runs predate the tot_air column. For those the domain measure is taken
from the staged geometry .opts and porosity falls back to 1 - tot_ice/|Omega|,
which is exact for a Cartesian domain and NOT applied to axisymmetric runs.

Usage:
    python3 plot_porosity.py --dir <run> [--save <png>] [--time-unit h]
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")          # before pyplot: these run headless on HPC
import matplotlib.pyplot as plt

from pplib import (AIR, ICE, SSA, TIME, auto_time_unit, in_time_unit,
                   load_ssa, opt_float, read_opts)


def domain_measure(run_dir):
    """Domain area/volume [m^d] from the staged geometry opts, or None.

    Only used for legacy files that lack the tot_air column. Returns None for
    axisymmetric runs, where the integrals carry a 2*pi*r weight that Lx*Ly
    does not represent.
    """
    opts = read_opts(run_dir)
    if opts.get("-axisym", "0") not in ("0", ""):
        return None
    dim = int(opt_float(opts, "-dim", 2) or 2)
    extent = 1.0
    for key in ("-Lx", "-Ly", "-Lz")[:dim]:
        L = opt_float(opts, key)
        if L is None or L <= 0:
            return None
        extent *= L
    return extent


def porosity_series(data, run_dir):
    """(porosity, how) from an SSA_evo array. `how` names the method used."""
    if data.shape[1] > AIR:
        ice, air = data[:, ICE], data[:, AIR]
        total = ice + air
        if np.all(total > 0):
            return air / total, "tot_air / (tot_ice + tot_air)"

    vol = domain_measure(run_dir)
    if vol is None:
        return None, ("no tot_air column, and the domain measure is not "
                      "recoverable from the staged opts (axisymmetric, or "
                      "-Lx/-Ly missing)")
    return 1.0 - data[:, ICE] / vol, f"1 - tot_ice / {vol:.4e} m^d  (legacy)"


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dir", default=".", help="run directory (default: cwd)")
    p.add_argument("--save", default=None,
                   help="output PNG (default: <dir>/porosity.png)")
    p.add_argument("--time-unit", choices=["s", "min", "h", "d"], default=None,
                   help="x-axis unit (default: chosen from the run length)")
    p.add_argument("--title", default=None)
    a = p.parse_args(argv)

    data = load_ssa(a.dir)
    if data is None:
        print(f"❌ No usable SSA_evo.dat in {a.dir}", file=sys.stderr)
        return 1

    phi_pore, how = porosity_series(data, a.dir)
    if phi_pore is None:
        print(f"❌ Cannot compute porosity: {how}", file=sys.stderr)
        return 1

    t = data[:, TIME]
    unit = a.time_unit or auto_time_unit(float(t.max()) if len(t) else 0.0)
    tt = in_time_unit(t, unit)
    ssa = data[:, SSA]

    fig, ax1 = plt.subplots(figsize=(10, 6))

    ax1.plot(tt, phi_pore, "-", color="#c0392b", lw=1.8, label="porosity")
    ax1.set_xlabel(f"Time [{unit}]", fontsize=14)
    ax1.set_ylabel("Porosity  $\\phi_a$", fontsize=14, color="#c0392b")
    ax1.tick_params(axis="y", labelcolor="#c0392b")
    ax1.tick_params(axis="both", which="major", labelsize=11)

    ax2 = ax1.twinx()
    ax2.plot(tt, ssa, "-", color="#1b4f72", lw=1.8, label="interface density")
    ax2.set_ylabel("Ice-air interface density  $\\int\\phi^2\\phi_a^2\\,dV/\\epsilon$",
                   fontsize=13, color="#1b4f72")
    ax2.tick_params(axis="y", labelcolor="#1b4f72")
    ax2.tick_params(axis="both", which="major", labelsize=11)

    d0, d1 = phi_pore[0], phi_pore[-1]
    ax1.set_title(a.title or
                  f"Porosity and interface density\n"
                  f"$\\phi_a$: {d0:.4f} → {d1:.4f}  "
                  f"({(d1 - d0) * 100:+.2f} pp)   over {tt[-1]:.3g} {unit}",
                  fontsize=15)

    lines = ax1.get_lines() + ax2.get_lines()
    ax1.legend(lines, [ln.get_label() for ln in lines], loc="best", fontsize=11)
    fig.tight_layout()

    out = a.save or os.path.join(a.dir, "porosity.png")
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"  porosity: {d0:.6f} → {d1:.6f}   ({how})")
    print(f"  wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
