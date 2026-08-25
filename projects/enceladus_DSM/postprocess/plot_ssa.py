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
true 3D area in m^2 and the same factor applies -- the axis label and the
analytic reference both follow -axisym, not -dim.

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

    Planar 2D: sum of circumferences, 2*pi*r [m]. 3D AND AXISYMMETRIC: sum of
    sphere areas, 4*pi*r^2 [m^2].

    AXISYMMETRIC RUNS ARE 3D. Integration() in src/assembly.c carries a 2*pi*r
    weight under -axisym, so column 0 of SSA_evo.dat is a true 3D area even
    though -dim is 2. Branching on -dim alone divided an area in m^2 by a
    length in m -- for the Molaro pair that is 1.94e-07 over 1.09e-03, so the
    "fraction of starting surface" read ~1.8e-04 instead of ~1.
    Grains in contact overlap slightly at their necks, so this is a mild
    OVERESTIMATE of the true free surface -- it is a reference scale, not a
    measurement.
    """
    radii = grain_radii(run_dir)
    if radii is None or len(radii) == 0:
        return None, 0
    opts = read_opts(run_dir)
    dim = int(opt_float(opts, "-dim", 2) or 2)
    axisym = opts.get("-axisym", "0") not in ("0", "")
    if dim >= 3 or axisym:
        return float(np.sum(4.0 * math.pi * radii ** 2)), len(radii)
    return float(np.sum(2.0 * math.pi * radii)), len(radii)


def fit_power_law(t, y, t_lo, t_hi):
    """Least-squares fit of  y = C * t^n  over t_lo <= t <= t_hi.

    Done as OLS of log(y) on log(t), which is the convention every coarsening
    exponent in this repo is quoted in (see postprocess/fit_neck_growth.py and
    the Demmenie/Kuczynski comparisons). Returns None if the window holds
    fewer than three usable points.

    Returns a dict with the exponent, its standard error, R^2, the window
    actually used, and the LOCAL log-slope at each end of that window. The
    local slopes are the honest check on whether a single exponent means
    anything: if they disagree with n, the curve is not a power law and the
    fitted number is an artefact of the window, not a property of the run.
    """
    m = (t > 0) & (y > 0) & (t >= t_lo) & (t <= t_hi)
    if m.sum() < 3:
        return None
    x = np.log(t[m])
    z = np.log(y[m])

    n, logC = np.polyfit(x, z, 1)
    resid = z - (n * x + logC)
    dof = len(x) - 2
    ss_res = float(np.sum(resid ** 2))
    ss_tot = float(np.sum((z - z.mean()) ** 2))
    sxx = float(np.sum((x - x.mean()) ** 2))
    stderr = math.sqrt(ss_res / dof / sxx) if dof > 0 and sxx > 0 else float("nan")
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    # Local log-slope d(log y)/d(log t), smoothed over the ends of the window.
    if len(x) >= 8:
        loc = np.gradient(z, x)
        k = max(3, len(loc) // 20)
        slope_lo = float(np.median(loc[:k]))
        slope_hi = float(np.median(loc[-k:]))
    else:
        slope_lo = slope_hi = float("nan")

    return {"n": float(n), "C": float(np.exp(logC)), "stderr": stderr, "r2": r2,
            "t_lo": float(t[m].min()), "t_hi": float(t[m].max()),
            "npts": int(m.sum()), "slope_lo": slope_lo, "slope_hi": slope_hi}


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
    p.add_argument("--no-fit", action="store_true",
                   help="skip the power-law fit")
    p.add_argument("--fit-from", type=float, default=0.01, metavar="FRAC",
                   help="start of the fit window as a fraction of the run "
                        "length (default 0.01, i.e. drop the first 1%% as "
                        "initial-condition relaxation transient)")
    p.add_argument("--fit-to", type=float, default=1.0, metavar="FRAC",
                   help="end of the fit window as a fraction of the run length "
                        "(default 1.0)")
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

    _opts = read_opts(a.dir)
    dim = int(opt_float(_opts, "-dim", 2) or 2)
    _axisym = _opts.get("-axisym", "0") not in ("0", "")
    unit_label = "m$^2$" if (dim >= 3 or _axisym) else "m"

    a0, n_grains = analytic_initial_area(a.dir)

    # Power-law fit. Done in SECONDS so the exponent is independent of the
    # x-axis unit; only the prefactor C carries the unit.
    fit = None
    if not a.no_fit and len(t) > 3:
        span = float(t.max())
        fit = fit_power_law(t, area, a.fit_from * span, a.fit_to * span)

    fig, ax1 = plt.subplots(figsize=(10, 6))
    plot = ax1.loglog if a.loglog else ax1.plot
    plot(tt, area, "-", color="#1b4f72", lw=1.8, label="measured")

    if fit is not None:
        tf = np.linspace(fit["t_lo"], fit["t_hi"], 200)
        plot(in_time_unit(tf, unit), fit["C"] * tf ** fit["n"],
             "--", color="#c0392b", lw=1.6,
             label=(f"fit  $S \\propto t^{{{fit['n']:.4f}}}$   "
                    f"($R^2$ = {fit['r2']:.4f})"))
        ax1.legend(fontsize=11, loc="best")

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
    if fit is not None:
        lo_u = in_time_unit(np.array([fit["t_lo"], fit["t_hi"]]), unit)
        print(f"  power-law fit  S ~ C t^n   (OLS on log S vs log t)")
        print(f"    exponent n = {fit['n']:+.4f} ± {fit['stderr']:.4f}"
              f"   (R² = {fit['r2']:.4f})")
        print(f"    window     = {lo_u[0]:.4g} → {lo_u[1]:.4g} {unit}"
              f"  ({fit['npts']} points)")
        print(f"    local log-slope across that window: "
              f"{fit['slope_lo']:+.4f} → {fit['slope_hi']:+.4f}")
        drift = abs(fit["slope_hi"] - fit["slope_lo"])
        if not math.isnan(drift) and drift > 0.5 * max(abs(fit["n"]), 1e-12):
            print("    ⚠ the local slope moves by more than half the fitted "
                  "exponent across the window:")
            print("      this curve is not a clean power law, so n is a "
                  "property of the WINDOW as much as of the run.")
        print("    (a fitted exponent is a property of the curve AND the fit "
              "protocol — always quote n with its form and window)")
    if a0:
        print(f"  analytic initial surface from packing: {a0:.6g} {unit_label} "
              f"({n_grains} grain rows); measured/analytic = {area[0] / a0:.3f}")
    print(f"  wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
