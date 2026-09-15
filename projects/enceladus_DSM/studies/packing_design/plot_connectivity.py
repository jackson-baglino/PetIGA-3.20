#!/usr/bin/env python3
"""Plot the connectivity sweep: why 2D cannot give both phases at once."""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "preprocess"))
import figstyle as fs                                      # noqa: E402

HERE = Path(__file__).parent
rows = list(csv.DictReader((HERE / "connectivity.csv").open()))
phi = np.array([float(r["porosity"]) for r in rows])
Zb = np.array([float(r["coordination_at_band"]) for r in rows])
frac = np.array([float(r["open_largest_frac"]) for r in rows])
sxy = np.array([int(r["solid_perc_x"]) + int(r["solid_perc_y"]) for r in rows])
oxy = np.array([int(r["open_perc_x"]) + int(r["open_perc_y"]) for r in rows])

fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.0, 4.2))

# (a) the two percolation counts against each other: the window that isn't.
# Both axes are integer-valued, so seeds at one porosity land on the same dot
# and a 2-vs-1 split would be invisible. Jitter in x by seed index.
seed_i = np.array([int(r["name"].split("seed")[1]) for r in rows])
jit = (seed_i - seed_i.mean()) * 0.004
a1.plot(phi + jit - 0.002, sxy, "o", color=fs.C[0], ms=8,
        markeredgecolor="white", markeredgewidth=1.0,
        label="solid matrix")
a1.plot(phi + jit + 0.002, oxy, "s", color=fs.C[1], ms=8,
        markeredgecolor="white", markeredgewidth=1.0,
        label="usable pore")
a1.axhline(2, color=fs.MUTED, lw=1.0, ls="--")
a1.text(phi.min() - 0.012, 2.08, "both directions — what an REV wants",
        fontsize=fs.FS_NOTE, color=fs.MUTED, va="bottom")
fs.style(a1, "porosity", "axes percolating (0, 1 or 2)",
         "(a)  no porosity connects both phases", logy=False)
a1.set_ylim(-0.3, 2.5)
a1.set_yticks([0, 1, 2])
a1.set_xticks(sorted(set(np.round(phi, 2))))
a1.legend(fontsize=fs.FS_LEG, frameon=False, loc="center left",
          title="directions percolating", title_fontsize=fs.FS_LEG)

# (b) coordination and usable-pore fraction move in opposite directions
a2.plot(phi, Zb, "o", color=fs.C[0], ms=7, markeredgecolor="white",
        markeredgewidth=1.0, label=r"$Z$ at the diffuse band")
a2.axhline(4.0, color=fs.C[0], lw=1.0, ls=":")
a2.text(phi.max(), 4.06, "2D isostatic", fontsize=fs.FS_NOTE,
        color=fs.C[0], ha="right", va="bottom")
a2b = a2.twinx()
a2b.plot(phi, frac, "s", color=fs.C[1], ms=7, markeredgecolor="white",
         markeredgewidth=1.0, label="largest usable-pore cluster")
a2b.set_ylabel("largest usable-pore cluster  [fraction of void]",
               fontsize=fs.FS_LABEL, color=fs.C[1])
a2b.tick_params(labelsize=fs.FS_TICK, labelcolor=fs.C[1])
a2b.set_ylim(0, 1)
a2b.spines[["top"]].set_visible(False)
fs.style(a2, "porosity", r"coordination number $Z$",
         "(b)  the trade is monotonic", logy=False)
a2.tick_params(labelcolor=fs.C[0])
a2.yaxis.label.set_color(fs.C[0])

fig.tight_layout()
fs.save(fig, HERE, "connectivity", dpi=200)
