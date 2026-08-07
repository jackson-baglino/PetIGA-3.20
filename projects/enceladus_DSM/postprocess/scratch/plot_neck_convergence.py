#!/usr/bin/env python3
"""
plot_neck_convergence.py — neck width vs time across an eps sweep, plotted so
the convergence question is actually answerable.

WHY NOT JUST PLOT w(t) AND EYEBALL IT
The three Molaro arms are indistinguishable in ParaView, yet a naive
"neck growth" metric reports an 18% spread. Both readings are misleading, for
the same reason: at t = 0 the two spheres meet in a CREASE, and a C1 quadratic
B-spline basis cannot represent a kink. It rounds it over ~2 cells, so the
measured t = 0 neck is mesh-dependent -- 32.85 / 34.30 / 35.13 um at
eps = 0.348 / 0.603 / 0.858 um, a 6.9% spread against a sharp chord of 32.81 um.

That is an initial-condition artifact, not a difference in the dynamics:
  * it is localized AT the neck (RMS contour displacement over the whole grain
    is 0.04-0.06 um, while the max, at the neck, is 0.94-1.38 um = ~2.2 cells);
  * it decays as the crease rounds -- the spread falls from 6.9% to 0.66% by
    t = 7200 s, crossing 2% at t ~ 1290 s;
  * it is NOT a measurement error. Upsampling the chord scan 40x moves the t=0
    numbers by <0.25 um, so refining the measurement does not remove it.

So `growth = w(t_end) - w(0)` differences against a contaminated baseline and
inherits the whole 6.9%. This script therefore:
  1. plots w(t) directly against the Molaro data (absolute, no differencing);
  2. plots the inter-arm spread vs time, which shows the artifact decaying and
     makes the trustworthy time window explicit;
  3. reports growth from a REFERENCE TIME after the crease has rounded
     (--t-ref, default 1500 s), which is a fair convergence measure.

Colours are Okabe-Ito blue / vermillion / bluish-green: adjacent-pair OKLab
separation >= 18.7 normal and >= 11 under deuteranopia/protanopia. Blue vs green
under tritanopia is 5.4, so line style carries identity too -- never colour alone.

USAGE — point it at the run directories and it does everything:

  ./venv_enceladus/bin/python3 postprocess/plot_neck_convergence.py \\
      --runs ~/SimulationResults/enceladus_DSM/scratch/*/*/ \\
      --validation inputs/validation/molaro2019_fig11_T-20.csv \\
      --out neck_convergence.png

It measures each run with neck_width.py, reads eps from the geometry .opts the
run script staged into the output directory, and names each arm from that file.

The lower-level form is still available if you already have the CSVs:

  ... --csv epsstrict=a.csv epsmid=b.csv --eps epsstrict=3.479e-7 epsmid=6.03e-7
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Okabe-Ito; validated for CVD separation (see module docstring).
COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7"]
STYLES = ["-", "--", ":", "-."]

INK, MUTED, GRID = "#1a1a1a", "#5c5c5c", "#d8d8d8"


def _opt(text, key):
    import re
    m = re.search(rf"^{re.escape(key)}\s+(\S+)", text, re.MULTILINE)
    return m.group(1) if m else None


def discover(run_dirs, out_dir, axisym_override):
    """Measure each run and return ({arm: eps}, {arm: csv_path}).

    eps, the arm name and the axisymmetry flag all come from the geometry .opts
    that run_enceladus.sh stages into the output directory -- the same file the
    solver used -- so there is nothing to keep in sync by hand.
    """
    import subprocess
    here = Path(__file__).resolve().parent
    tmp = Path(out_dir) / "_neck_csv"
    tmp.mkdir(parents=True, exist_ok=True)

    # A geometry re-run leaves several timestamped directories per arm. Keep
    # only the NEWEST of each, so a fresh batch is compared against itself
    # rather than silently mixing runs from different builds.
    best = {}
    for rd in (Path(r) for r in run_dirs):
        if not (rd / "SSA_evo.dat").is_file():
            continue
        geo = [f for f in rd.glob("*.opts")
               if f.name != "solver.opts" and _opt(f.read_text(), "-eps")]
        if not geo:
            continue
        key = geo[0].stem
        mt = (rd / "SSA_evo.dat").stat().st_mtime
        if key not in best or mt > best[key][0]:
            best[key] = (mt, rd)
    if len(best) < len(list(run_dirs)):
        import datetime as _dt
        for k, (mt, rd) in sorted(best.items()):
            print(f"  using {rd.parent.name}/{rd.name}"
                  f"  ({_dt.datetime.fromtimestamp(mt):%Y-%m-%d %H:%M})",
                  file=sys.stderr)

    eps, csv = {}, {}
    for _mt, rd in sorted(best.values()):
        geo = [f for f in rd.glob("*.opts")
               if f.name != "solver.opts" and _opt(f.read_text(), "-eps")]
        txt = geo[0].read_text()
        arm = geo[0].stem
        eps[arm] = float(_opt(txt, "-eps"))
        ax = axisym_override
        if ax is None:
            ax = (_opt(txt, "-axisym") or "0") not in ("0", "false", "FALSE")
        out = tmp / f"{arm}__{int(_mt)}.csv"
        if not out.is_file():
            cmd = [sys.executable, str(here / "neck_width.py"), str(rd),
                   "--out", str(out)] + (["--axisym"] if ax else [])
            print(f"  measuring {arm} ...", file=sys.stderr)
            r = subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode != 0 or not out.is_file():
                print(f"  FAILED on {arm}:\n{r.stdout[-500:]}{r.stderr[-500:]}",
                      file=sys.stderr)
                eps.pop(arm, None)
                continue
        csv[arm] = out

    # Strip the prefix every arm shares -- the geometry names differ only in
    # their eps suffix, and the full names are too long to tabulate.
    if len(eps) > 1:
        names = list(eps)
        pre = os.path.commonprefix(names)
        pre = pre[:pre.rfind("_") + 1] if "_" in pre else ""
        if pre and all(len(n) > len(pre) for n in names):
            eps = {n[len(pre):]: v for n, v in eps.items()}
            csv = {n[len(pre):]: v for n, v in csv.items()}
    return eps, csv


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="+", default=None, metavar="RUN_DIR",
                    help="run directories; eps and the arm name are read from "
                         "the geometry .opts staged inside each one, and the "
                         "neck series is measured with neck_width.py")
    ap.add_argument("--csv", nargs="+", default=None, metavar="ARM=PATH",
                    help="alternative to --runs: neck_width.py output per arm")
    ap.add_argument("--eps", nargs="+", default=None, metavar="ARM=EPS")
    ap.add_argument("--axisym", action="store_true", default=None,
                    help="force axisymmetric measurement (auto-detected from "
                         "the staged geometry .opts when using --runs)")
    ap.add_argument("--validation", type=Path, default=None,
                    help="molaro2019_fig11_T-*.csv")
    ap.add_argument("--t-ref", type=float, default=1500.0,
                    help="reference time [s] for the growth metric; must be "
                         "after the initial crease has rounded")
    ap.add_argument("--out", type=Path, default=Path("neck_convergence.png"))
    args = ap.parse_args(argv)

    if args.runs:
        eps, csv = discover(args.runs, args.out.parent, args.axisym)
    elif args.csv and args.eps:
        eps = {k: float(v) for k, v in (s.split("=") for s in args.eps)}
        csv = {k: Path(v) for k, v in (s.split("=") for s in args.csv)}
    else:
        ap.error("give either --runs, or both --csv and --eps")
    if len(eps) < 2:
        ap.error(f"need at least 2 arms to compare, found {len(eps)}")
    arms = sorted(eps, key=lambda a: eps[a])           # finest first

    D = {}
    for a in arms:
        d = np.genfromtxt(csv[a], delimiter=",", names=True)
        D[a] = (d["t_s"], d["neck_width_m"] * 1e6)

    fig, (ax, ax2) = plt.subplots(
        2, 1, figsize=(7.2, 7.4), sharex=True,
        gridspec_kw=dict(height_ratios=[2.4, 1], hspace=0.12))

    # ---- panel 1: w(t), absolute -----------------------------------------
    if args.validation and args.validation.is_file():
        v = np.genfromtxt(args.validation, delimiter=",", comments="#")
        ax.errorbar(v[:, 0], v[:, 1], yerr=[v[:, 3], v[:, 2]],
                    fmt="o", ms=6, mfc="white", mec=INK, ecolor=MUTED,
                    elinewidth=1.2, capsize=3, zorder=5,
                    label="Molaro et al. 2019 (measured)")

    for i, a in enumerate(arms):
        t, w = D[a]
        ax.plot(t / 60.0, w, STYLES[i], color=COLORS[i], lw=2.0, zorder=4,
                label=f"$\\epsilon$ = {eps[a]*1e6:.3f} µm")

    ax.set_ylabel("neck width  [µm]", color=INK)
    ax.set_title("Neck width vs time — $\\epsilon$ convergence",
                 color=INK, fontsize=12, pad=10)
    ax.legend(frameon=False, loc="lower right", fontsize=9)
    ax.grid(True, color=GRID, lw=0.7)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)

    # ---- panel 2: inter-arm spread ---------------------------------------
    tmax = min(D[a][0][-1] for a in arms)
    tg = np.linspace(1.0, tmax, 400)
    W = np.array([np.interp(tg, *D[a]) for a in arms])
    spread = (W.max(axis=0) - W.min(axis=0)) / W[0] * 100.0

    ax2.plot(tg / 60.0, spread, "-", color=INK, lw=2.0)
    ax2.axhline(2.0, color=MUTED, lw=1.2, ls="--")
    ax2.text(tmax / 60.0, 2.15, "2% tolerance", ha="right", va="bottom",
             fontsize=8.5, color=MUTED)

    # Shade only up to where the spread actually falls below tolerance, and
    # put the two annotations on opposite sides of that line so they cannot
    # collide however the crossing time lands.
    below = tg[spread < 2.0]
    tc = below[0] if len(below) else tmax
    ax2.axvspan(0, tc / 60.0, color="#f2f2f2", zorder=0)
    ax2.axvline(tc / 60.0, color=MUTED, lw=1.0, ls=":")
    ytop = max(spread.max() * 1.15, 3)
    # When the crease resolves quickly the shaded band is too narrow to hold a
    # label; point at it from outside instead of writing over the axis.
    if tc / tmax > 0.08:
        ax2.text(tc / 120.0, ytop * 0.62, "initial crease\nunresolvable",
                 ha="center", va="center", fontsize=8.5, color=MUTED)
        ax2.annotate(f"spread < 2% from t = {tc:.0f} s",
                     xy=(tc / 60.0, 2.0), xytext=(tc / 60.0 + 8, ytop * 0.78),
                     fontsize=8.5, color=MUTED,
                     arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.9))
    else:
        ax2.annotate(f"initial crease unresolvable;\nspread < 2% from t = {tc:.0f} s",
                     xy=(tc / 60.0, 2.0),
                     xytext=(tmax / 60.0 * 0.16, ytop * 0.72),
                     fontsize=8.5, color=MUTED,
                     arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.9))

    ax2.set_xlabel("time  [min]", color=INK)
    ax2.set_ylabel("spread across\narms  [%]", color=INK)
    ax2.grid(True, color=GRID, lw=0.7)
    ax2.set_axisbelow(True)
    ax2.set_ylim(0, ytop)
    for s in ("top", "right"):
        ax2.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax2.spines[s].set_color(MUTED)

    for a in (ax, ax2):
        a.tick_params(colors=MUTED, labelcolor=INK)

    fig.savefig(args.out, dpi=160, bbox_inches="tight", facecolor="white")
    print(f"  wrote {args.out}")

    # ---- the table that goes with the figure -----------------------------
    print(f"\n{'arm':>12} {'eps[um]':>9} {'w(0)':>9} {'w(t_ref)':>10} "
          f"{'w(end)':>9} {'growth from t_ref':>18}")
    ref = arms[0]
    g = {}
    for a in arms:
        t, w = D[a]
        w0, wr, we = w[0], np.interp(args.t_ref, t, w), w[-1]
        g[a] = we - wr
        print(f"{a:>12} {eps[a]*1e6:>9.4f} {w0:>9.4f} {wr:>10.4f} {we:>9.4f} "
              f"{g[a]:>17.4f}")
    print(f"\n  vs finest ({ref}), using t_ref = {args.t_ref:g} s:")
    for a in arms[1:]:
        print(f"    {a:>12}  growth {abs(g[a]/g[ref]-1)*100:6.2f}%   "
              f"final width {abs(D[a][1][-1]/D[ref][1][-1]-1)*100:6.2f}%")
    print(f"\n  For comparison, growth measured from t=0 instead (contaminated):")
    for a in arms[1:]:
        g0 = D[a][1][-1] - D[a][1][0]
        gr0 = D[ref][1][-1] - D[ref][1][0]
        print(f"    {a:>12}  {abs(g0/gr0-1)*100:6.2f}%")
    return 0


if __name__ == "__main__":
    sys.exit(main())
