#!/usr/bin/env python3
"""
eps_convergence.py — quantify what refining eps actually changes.

Written because three Molaro arms spanning eps = 0.35 to 0.86 um look identical
in ParaView. They are not: the neck GROWTH differs by 23%. Eyeballing a contour
cannot resolve this, and neither can comparing absolute quantities.

FOUR TRAPS THIS AVOIDS

1. Comparing ABSOLUTE quantities. Every arm starts from the same geometry, so
   total ice, SSA and neck width all agree to <0.1% at t=0 and stay dominated by
   that shared initial value. Compare the CHANGE (X(t) - X(0)) instead: on the
   Molaro runs, absolute SSA agrees to 0.10% while the SSA DROP differs by 3.8%.

2. Assuming the arms START the same. They frequently do not, and then you are
   measuring initial-condition error, not discretisation error. On the Molaro
   runs the t=0 neck is 32.85 / 34.30 / 35.13 um across the three eps -- a 7%
   spread -- because the diffuse interface rounds the sharp crease where the two
   spheres meet, over a scale ~eps. That rounding is inherent to the diffuse
   representation and no amount of geometry calibration removes it. So the
   report always prints the t=0 spread first; if it is comparable to the
   difference you are chasing, the comparison is not about resolution.

3. Normalising a time series by an instantaneous value that starts at zero.
   The relative error in DELTA-SSA(t) blows up to 65% at early times purely
   because the denominator is ~0. Normalise by the TOTAL change over the run.

4. Reporting only the finest-vs-coarsest gap. With three arms you can fit
   f(eps) = f0 + C*eps^p, which gives the observed order p AND an eps->0
   extrapolation to measure every arm against -- a far better reference than
   the finest arm, which has its own error.

USAGE
  # SSA / ice from SSA_evo.dat, auto-discovering arms under a results root
  ./postprocess/eps_convergence.py <results-root> --eps epsstrict=3.479e-7 \\
        epsmid=6.030e-7 epsloose=8.584e-7

  # add a neck-width series produced by neck_width.py
  ./postprocess/eps_convergence.py <results-root> --eps ... \\
        --neck-csv-dir <dir-with-{arm}.csv>
"""

from __future__ import annotations

import argparse
import glob
import os
import sys

import numpy as np
from pplib import load_ssa, SSA, ICE, TIME, MASS

try:
    from scipy.optimize import brentq
except ImportError:
    brentq = None


def power_fit(eps, f):
    """Fit f(eps) = f0 + C*eps^p through 3 points -> (p, f0). None if ill-posed.

    Requires the three values to be monotone in eps and to actually bend; a
    straight line or non-monotone data has no valid p and we say so rather than
    reporting a fitted number that means nothing.
    """
    if brentq is None or len(eps) != 3:
        return None
    e, y = np.asarray(eps, float), np.asarray(f, float)
    o = np.argsort(e)
    e, y = e[o], y[o]
    if not (np.all(np.diff(y) > 0) or np.all(np.diff(y) < 0)):
        return None

    def g(p):
        return (y[1] - y[0]) * (e[2] ** p - e[0] ** p) - \
               (y[2] - y[0]) * (e[1] ** p - e[0] ** p)

    try:
        lo, hi = 0.05, 8.0
        if g(lo) * g(hi) > 0:
            return None
        p = brentq(g, lo, hi)
        C = (y[1] - y[0]) / (e[1] ** p - e[0] ** p)
        return p, y[0] - C * e[0] ** p
    except Exception:
        return None


def report(name, arms, eps, series, unit="", scale=1.0):
    """arms: ordered finest->coarsest. series[arm] = (t, y)."""
    print(f"\n{'='*72}\n  {name}\n{'='*72}")
    y0 = {a: series[a][1][0] for a in arms}
    spread = (max(y0.values()) - min(y0.values())) / abs(np.mean(list(y0.values())))
    print(f"\n  t=0 values  (spread {spread*100:.2f}% -- if this is comparable to")
    print(f"              the differences below, you are measuring IC error)")
    for a in arms:
        print(f"    {a:>12}  eps={eps[a]*1e6:6.4f} um   {y0[a]*scale:12.5g} {unit}")

    print(f"\n  CHANGE over the run (the signal):")
    ch = {}
    for a in arms:
        t, y = series[a]
        ch[a] = y[-1] - y[0]
        rel = ch[a] / y[0] if y[0] else float("nan")
        print(f"    {a:>12}  {ch[a]*scale:12.5g} {unit}   ({rel*100:+7.3f}% of initial)")

    ref = arms[0]
    print(f"\n  vs finest ({ref}):")
    for a in arms[1:]:
        print(f"    {a:>12}  {abs(ch[a]/ch[ref]-1)*100:7.2f}%")

    fit = power_fit([eps[a] for a in arms], [ch[a] for a in arms])
    if fit:
        p, f0 = fit
        print(f"\n  Richardson f(eps) = f0 + C*eps^p:  order p = {p:.2f}")
        print(f"    extrapolated eps->0: {f0*scale:.5g} {unit}")
        for a in arms:
            print(f"    {a:>12}  err vs eps->0 = {abs(ch[a]/f0-1)*100:6.2f}%")
    else:
        print("\n  Richardson fit: not well posed (values not monotone in eps,")
        print("    or no curvature) -- usually means IC error dominates.")

    tmax = min(series[a][0][-1] for a in arms)
    tg = np.linspace(0, tmax, 300)
    gr = np.interp(tg, *series[ref])
    denom = abs(gr[-1] - gr[0]) or 1.0
    print(f"\n  time-resolved, normalised by the finest arm's total change:")
    for a in arms[1:]:
        g = np.interp(tg, *series[a])
        d = (g - gr) / denom
        print(f"    {a:>12}  L2 {np.sqrt(np.mean(d**2))*100:6.2f}%   "
              f"max {np.max(np.abs(d))*100:6.2f}%")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("results", help="root holding the run directories")
    ap.add_argument("--eps", nargs="+", required=True,
                    metavar="ARM=EPS", help="e.g. epsmid=6.03e-7")
    ap.add_argument("--neck-csv-dir", default=None,
                    help="directory with <arm>.csv from neck_width.py")
    args = ap.parse_args(argv)

    eps = {}
    for kv in args.eps:
        k, v = kv.split("=")
        eps[k] = float(v)
    arms = sorted(eps, key=lambda a: eps[a])          # finest first

    ssa, missing = {}, []
    for a in arms:
        hits = glob.glob(os.path.join(args.results, f"*{a}*", "*", "SSA_evo.dat")) \
            or glob.glob(os.path.join(args.results, f"*{a}*", "SSA_evo.dat"))
        if not hits:
            missing.append(a); continue
        a_ = load_ssa(os.path.dirname(hits[0]), min_cols=8)
        if a_ is not None:
            ssa[a] = dict(t=a_[:, TIME], ssa=a_[:, SSA],
                          ice=a_[:, ICE], mass=a_[:, MASS])
    if missing:
        print(f"  no SSA_evo.dat for: {', '.join(missing)}", file=sys.stderr)
    if len(ssa) == len(arms):
        report("SSA  (sub_interf/eps -- interfacial area proxy)", arms, eps,
               {a: (ssa[a]["t"], ssa[a]["ssa"]) for a in arms})
        icech = {a: abs(ssa[a]["ice"][-1] / ssa[a]["ice"][0] - 1) for a in arms}
        print(f"\n  ice-mass drift (conservation check, not a convergence metric):")
        for a in arms:
            print(f"    {a:>12}  {icech[a]*100:8.4f}%")
        if max(icech.values()) < 1e-9:
            print("    -> exactly conserved: closed system, so total ice cannot")
            print("       discriminate between arms. Use SSA or neck width.")

    if args.neck_csv_dir:
        neck = {}
        for a in arms:
            f = os.path.join(args.neck_csv_dir, f"{a}.csv")
            if os.path.isfile(f):
                d = np.genfromtxt(f, delimiter=",", names=True)
                neck[a] = (d["t_s"], d["neck_width_m"])
        if len(neck) == len(arms):
            report("NECK WIDTH", arms, eps, neck, unit="um", scale=1e6)
    return 0


if __name__ == "__main__":
    sys.exit(main())
