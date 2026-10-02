#!/usr/bin/env python3
"""plot_data.py -- the Thomas et al. (1994) neck data as digitised from
Molaro et al. (2019) Fig. 8(c), the d_fixed power-law fit over the window the
proposed run resolves, and that run's mesh resolution floor.

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
    fig, ax = plt.subplots(figsize=(5.5, 3.8), constrained_layout=True)
    tf = np.geomspace(t[2], t[-1], 50)
    ax.plot(tf, C * tf**a, color=BLUE, lw=1.5, zorder=2)
    ax.text(tf[22], C * tf[22]**a * 1.10,
            f"x/a = {C:.3f} t^{a:.2f}   (a ± {float(fit['a_ci']):.2f}, 95 %)",
            color=INK, fontsize=9, ha="right")
    ax.plot(t[2:], u[2:], "o", ms=6, mfc=BLUE, mec="white", mew=1.0, zorder=3)
    ax.plot(t[:2], u[:2], "o", ms=6, mfc="white", mec=BLUE, mew=1.2, zorder=3)
    ax.axhline(FLOOR, color=MUTED, lw=1.0, ls="--", zorder=1)
    ax.text(t.max(), FLOOR, f"mesh floor x/a = {FLOOR:.3f}  ",
            va="bottom", ha="right", color=MUTED, fontsize=9)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("time as reported, t [h]")
    ax.set_ylabel("relative neck size, x/a")
    plain = FuncFormatter(lambda v, _: f"{v:g}")
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_formatter(plain)
        axis.set_minor_formatter(NullFormatter())
    ax.set_yticks([0.1, 0.15, 0.2, 0.25, 0.3, 0.4])
    ax.set_ylim(0.085, 0.53)
    ax.grid(True, which="major", color=GRID, lw=0.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"thomas_neck_data.{ext}", dpi=200)
    print(f"wrote {FIG}/thomas_neck_data.png/.pdf")


if __name__ == "__main__":
    main()
