#!/usr/bin/env python3
"""plot_data.py -- the Thomas et al. (1994) neck data as digitised from
Molaro et al. (2019) Fig. 8(c), the d_fixed power-law fit over the window the
proposed run resolves, and that run's mesh resolution floor -- as relative
neck size x/a and as absolute neck width 2x [um].

Run from enceladus_DSM/ (after fit_power_law.py):
    python studies/historical_sintering/analysis/plot_data.py
"""

import csv
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, NullFormatter

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
FIG = HERE.parent / "figures"

INK, MUTED, GRID, BLUE = "#1a1a1a", "#5c5c5c", "#d8d8d8", "#0072B2"
FLOOR = 0.9 * 0.1811          # mesh floor: 0.9 x point 3
A_UM = 120.0                  # grain radius, Thomas: r ~ 120 um


def main():
    FIG.mkdir(exist_ok=True)
    lines = [l for l in open(DATA / "thomas1994_T-20_r120um.csv")
             if not l.startswith("#")]
    t, u = np.loadtxt(lines[1:], delimiter=",", unpack=True)
    fit = next(r for r in csv.DictReader(open(DATA / "thomas_powerlaw_fits.csv"))
               if r["window"] == "pts 3-11" and r["form"] == "d_fixed")
    a, C = float(fit["a"]), float(fit["C"])

    plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED,
                         "axes.labelcolor": INK, "xtick.color": MUTED,
                         "ytick.color": MUTED})
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 3.8), constrained_layout=True)
    # (a) relative neck x/a, as published; (b) absolute neck WIDTH 2x = 2a(x/a),
    # the quantity neck_width.py --axisym reports and our Molaro figures show.
    # a is Thomas's "r ~ 120 um", so (b) carries that radius uncertainty.
    panels = [
        (1.0, "relative neck size, x/a", "x/a", "{:.3f}",
         [0.1, 0.15, 0.2, 0.25, 0.3, 0.4], "(a)"),
        (2 * A_UM, "neck width, 2x [µm]", "2x", "{:.1f} µm",
         [25, 30, 40, 50, 60, 80, 100], "(b)"),
    ]
    tf = np.geomspace(t[2], t[-1], 50)
    plain = FuncFormatter(lambda v, _: f"{v:g}")
    for ax, (k, ylabel, sym, ffmt, yticks, tag) in zip(axes, panels):
        ax.plot(tf, k * C * tf**a, color=BLUE, lw=1.5, zorder=2)
        ax.text(tf[22], k * C * tf[22]**a * 1.10,
                f"{sym} = {k * C:.3g} t^{a:.2f}   (a ± {float(fit['a_ci']):.2f}, 95 %)",
                color=INK, fontsize=9, ha="right")
        ax.plot(t[2:], k * u[2:], "o", ms=6, mfc=BLUE, mec="white", mew=1.0,
                zorder=3)
        ax.plot(t[:2], k * u[:2], "o", ms=6, mfc="white", mec=BLUE, mew=1.2,
                zorder=3)
        ax.axhline(k * FLOOR, color=MUTED, lw=1.0, ls="--", zorder=1)
        ax.text(t.max(), k * FLOOR, f"mesh floor {sym} = {ffmt.format(k * FLOOR)}  ",
                va="bottom", ha="right", color=MUTED, fontsize=9)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("time as reported, t [h]")
        ax.set_ylabel(ylabel)
        for axis in (ax.xaxis, ax.yaxis):
            axis.set_major_formatter(plain)
            axis.set_minor_formatter(NullFormatter())
        ax.set_yticks(yticks)
        ax.set_ylim(k * 0.085, k * 0.53)
        ax.set_title(f"{tag} Thomas et al. (1994), −20 °C, a ≈ {A_UM:g} µm",
                     fontsize=10, color=INK, loc="left")
        ax.grid(True, which="major", color=GRID, lw=0.5)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"thomas_neck_data.{ext}", dpi=200)
    print(f"wrote {FIG}/thomas_neck_data.png/.pdf")


if __name__ == "__main__":
    main()
