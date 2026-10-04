#!/usr/bin/env python3
"""Why does k_yy overtake k_xx at high porosity? Ice chord lengths vs contacts.

    venv_enceladus/bin/python studies/rve_anisotropy/chord_anisotropy.py <campaign dir>

The 3a/3b runs (2026-10-04) gave seed-mean k_xx/k_yy of 1.02, 1.02, 0.91, 0.79,
0.64 for phi 0.275 ... 0.475, while the contact fabric F_xx/F_yy of the same
packings is > 1 throughout (1.04 ... 1.12). Contact ORIENTATION therefore does
not set the conduction anisotropy at high porosity. The hypothesis: gravity
deposition stacks grains into vertical columns, so the ice is CONTINUOUS for
longer along y even though more of its contacts lean sideways, and once the
solid thins (high phi) the long continuous paths are what conduct.

Test, on the initial packings (t = 0, sharp geometry, periodic raster):
  mean ice chord length along x and along y (lineal run lengths of solid),
  and its ratio L_x/L_y, against k_xx/k_yy at the opening sample and at 30 d
  (every production packing, seeds 1601-2005, -20 C) and against F_xx/F_yy.

Writes chord_anisotropy.csv and chord_anisotropy.png here.
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
PROJ = HERE.parents[1]
sys.path.insert(0, str(PROJ / "preprocess"))
from packing_lib import rasterize  # noqa: E402

N = 2000  # 1 um pixels on a 2 mm box: the smallest grains are ~10 um


def read_grains(f):
    rows = [l.split() for l in open(f) if l.strip() and not l.startswith("#")]
    Lx, Ly = map(float, rows[0][:2])
    a = np.array([[float(v) for v in r[:3]] for r in rows[1:]])
    return Lx, Ly, a[:, :2], a[:, 2]


def mean_chord(mask, axis):
    """Mean length (pixels) of solid runs along axis, periodic."""
    m = mask if axis == 1 else mask.T          # runs along each row
    lengths = []
    for row in m:
        if row.all():
            lengths.append(row.size); continue
        if not row.any():
            continue
        r = np.roll(row, -int(np.argmin(row)))  # start on a gap: no wrap split
        d = np.diff(np.concatenate(([0], r.astype(np.int8), [0])))
        lengths.extend(np.flatnonzero(d == -1) - np.flatnonzero(d == 1))
    return float(np.mean(lengths))


def main():
    camp = Path(sys.argv[1])
    fab = {int(r["base_seed"]): r for r in
           csv.DictReader(open(PROJ / "inputs/packings/keff_LR40/packings_summary.csv"))}
    rows = []
    for pk in sorted((PROJ / "inputs/packings/keff_LR40").glob("phi*")):
        seed = int(pk.name.split("seed")[1])
        phi = float(re.search(r"phi([\d.]+)", pk.name).group(1))
        Lx, Ly, c, r = read_grains(pk / "grains.dat")
        solid = rasterize(c, r, Lx, Ly, N, N, True, True)   # [y, x]
        cx, cy = mean_chord(solid, 1), mean_chord(solid, 0)
        runs = list(camp.glob(f"packing_2D_{pk.name}_*_T-20__*/k_eff.csv"))
        if not runs:
            continue
        k = np.genfromtxt(runs[0], delimiter=",", names=True)
        i0 = int(np.argmax(k["time"] >= 1.0))
        rows.append(dict(phi=phi, seed=seed, chord_x_um=cx * Lx / N * 1e6, chord_y_um=cy * Ly / N * 1e6,
                         chord_ratio=cx / cy, F_ratio=float(fab[seed]["F_ratio"]),
                         k_ratio_0=k["k_00"][i0] / k["k_11"][i0], k_ratio_30=k["k_00"][-1] / k["k_11"][-1]))
    with open(HERE / "chord_anisotropy.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    phis = sorted({r["phi"] for r in rows})
    print(f"{'phi':>6} {'L_x/L_y':>12} {'F_xx/F_yy':>10} {'kxx/kyy t0':>12} {'kxx/kyy 30d':>12}")
    for p in phis:
        g = [r for r in rows if r["phi"] == p]
        m = lambda k: (np.mean([x[k] for x in g]), np.std([x[k] for x in g], ddof=1))
        print(f"{p:6.3f} {m('chord_ratio')[0]:6.3f}±{m('chord_ratio')[1]:.3f} {m('F_ratio')[0]:10.3f} "
              f"{m('k_ratio_0')[0]:6.3f}±{m('k_ratio_0')[1]:.3f} {m('k_ratio_30')[0]:6.3f}±{m('k_ratio_30')[1]:.3f}")
    a = np.array([[r["chord_ratio"], r["F_ratio"], r["k_ratio_0"], r["k_ratio_30"]] for r in rows])
    lk0, lk30 = np.log(a[:, 2]), np.log(a[:, 3])
    for j, name in ((0, "chord L_x/L_y"), (1, "contact F_xx/F_yy")):
        x = np.log(a[:, j])
        print(f"  corr(log {name}, log k ratio): t0 {np.corrcoef(x, lk0)[0, 1]:+.2f}, "
              f"30 d {np.corrcoef(x, lk30)[0, 1]:+.2f}   (n = {len(a)})")

    cmap = plt.get_cmap("copper")
    col = {p: cmap(0.15 + 0.7 * i / (len(phis) - 1)) for i, p in enumerate(phis)}
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6), constrained_layout=True)
    for j, (xk, xl) in enumerate((("chord_ratio", r"ice chord ratio $L_x/L_y$ (t = 0)"),
                                  ("F_ratio", r"contact fabric $F_{xx}/F_{yy}$ (t = 0)"))):
        for p in phis:
            g = [r for r in rows if r["phi"] == p]
            ax[j].scatter([r[xk] for r in g], [r["k_ratio_30"] for r in g], s=36,
                          color=col[p], edgecolor="white", linewidth=0.8, label=f"φ = {p}")
        lim = [0.4, 1.6]
        ax[j].plot(lim, lim, color="#888888", lw=0.8, ls=":")
        ax[j].axhline(1, color="#bbbbbb", lw=0.6); ax[j].axvline(1, color="#bbbbbb", lw=0.6)
        ax[j].set(xlabel=xl, ylabel=r"$k_{xx}/k_{yy}$ at 30 d")
        ax[j].grid(True, alpha=0.25)
        for sp in ("top", "right"):
            ax[j].spines[sp].set_visible(False)
    ax[0].legend(frameon=False, fontsize=8)
    ax[0].set_title("(a) continuity of the ice", loc="left")
    ax[1].set_title("(b) orientation of contacts", loc="left")
    fig.savefig(HERE / "chord_anisotropy.png", dpi=170)
    print(f"wrote {HERE / 'chord_anisotropy.png'} and .csv")


if __name__ == "__main__":
    main()
