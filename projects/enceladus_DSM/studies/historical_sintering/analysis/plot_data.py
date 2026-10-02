#!/usr/bin/env python3
"""plot_data.py -- the Kingery (1960) and Thomas et al. (1994) neck data as
digitised from Molaro et al. (2019) Fig. 8, with the mesh resolution floor
each proposed run is sized to.

Run from enceladus_DSM/:
    python studies/historical_sintering/analysis/plot_data.py
"""

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

PANELS = [
    # file, title, floor u (0.9 x first point), x-axis time unit
    ("kingery1960_T-17.8_r110um.csv",
     "Kingery (1960), T = -17.8 °C, a = 110 µm", 0.9 * 0.2485, "min", 60.0),
    ("thomas1994_T-20_r120um.csv",
     "Thomas et al. (1994), T = -20 °C, a ≈ 120 µm", 0.9 * 0.0992, "h", 1.0),
]


def main():
    FIG.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED,
                         "axes.labelcolor": INK, "xtick.color": MUTED,
                         "ytick.color": MUTED})
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.6), constrained_layout=True)
    for ax, (fname, title, floor, unit, mult), tag in zip(axes, PANELS, "ab"):
        lines = [l for l in open(DATA / fname) if not l.startswith("#")]
        t, u = np.loadtxt(lines[1:], delimiter=",", unpack=True)
        t = t * mult
        ax.plot(t, u, "o", ms=6, mfc=BLUE, mec="white", mew=1.0, color=BLUE,
                zorder=3)
        ax.axhline(floor, color=MUTED, lw=1.0, ls="--", zorder=2)
        ax.text(t.min(), floor, f"  mesh floor x/a = {floor:.3f}",
                va="bottom", ha="left", color=MUTED, fontsize=9)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(f"time as reported, t [{unit}]")
        ax.set_ylabel("relative neck size, x/a")
        ax.set_title(f"({tag}) {title}", fontsize=10, color=INK, loc="left")
        plain = FuncFormatter(lambda v, _: f"{v:g}")
        for axis in (ax.xaxis, ax.yaxis):
            axis.set_major_formatter(plain)
            axis.set_minor_formatter(NullFormatter())
        ax.set_yticks([0.1, 0.15, 0.2, 0.25, 0.3, 0.4] if u.min() < 0.2
                      else [0.2, 0.25, 0.3, 0.35])
        if unit == "min":
            ax.set_xticks([1.5, 2, 3, 4, 6, 8])
        ax.grid(True, which="major", color=GRID, lw=0.5)
        ax.set_ylim(floor * 0.85, u.max() * 1.2)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"historical_neck_data.{ext}", dpi=200)
    print(f"wrote {FIG}/historical_neck_data.png/.pdf")


if __name__ == "__main__":
    main()
