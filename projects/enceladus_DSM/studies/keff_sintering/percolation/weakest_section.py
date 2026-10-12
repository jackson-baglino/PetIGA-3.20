#!/usr/bin/env python3
"""weakest_section.py — how many ice bridges carry heat across the cell?

    venv_enceladus/bin/python studies/keff_sintering/percolation/weakest_section.py

For each production packing at -20 C, last snapshot (30 d): heat conducted in
x must cross every vertical cut (a column of the image); in y, every
horizontal cut. For each direction the script reports the weakest cut -- the
one with the least ice on it -- as the ice fraction on that cut and the number
of separate ice bridges crossing it, and the mean number of bridges over all
cuts. Few bridges at the weakest cut means the conductivity in that direction
is set by a handful of necks, so it is not an average over the domain.

Also k_xx / k_yy of each packing at the same instant (k_eff.csv).
Writes weakest_section.txt and weakest_section.csv beside this file.
"""
from __future__ import annotations
import csv, re
from collections import defaultdict
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
L = Path.home() / "SimulationResults/HPC_results/enceladus_DSM/keff_sintering_campaign"
RUN = re.compile(r"packing_2D_phi([0-9.]+)_Rave50um_LR40_seed(\d+)_L2mm_eps1000nm_perxy_T-20__snow_T-20_h1\.00_30d$")


def bridges(line) -> int:
    """Separate ice segments along a periodic line."""
    return int((np.diff(np.concatenate([line[-1:], line]).astype(int)) == 1).sum())


def main():
    rows = []
    for d in sorted(L.iterdir()):
        m = RUN.match(d.name)
        if not m or int(m.group(2)) < 1000:
            continue
        a = np.fromfile(sorted(d.glob("sol_*.dat"))[-1], dtype=">f8", offset=8)
        n = int(round((a.size / 3) ** 0.5))
        ice = a.reshape(n, n, 3)[: n - 2 : 2, : n - 2 : 2, 0] > 0.5          # [y, x], one period
        k = np.genfromtxt(d / "k_eff.csv", delimiter=",", names=True)
        col, row = ice.mean(axis=0), ice.mean(axis=1)
        ic, ir = int(col.argmin()), int(row.argmin())
        rows.append(dict(phi=float(m.group(1)), seed=int(m.group(2)),
                         x_min_ice=col.min(), x_bridges_min=bridges(ice[:, ic]),
                         x_bridges_mean=np.mean([bridges(ice[:, j]) for j in range(0, ice.shape[1], 10)]),
                         x_bridges_fewest=min(bridges(ice[:, j]) for j in range(0, ice.shape[1], 2)),
                         y_min_ice=row.min(), y_bridges_min=bridges(ice[ir, :]),
                         y_bridges_mean=np.mean([bridges(ice[j, :]) for j in range(0, ice.shape[0], 10)]),
                         y_bridges_fewest=min(bridges(ice[j, :]) for j in range(0, ice.shape[0], 2)),
                         kxx_over_kyy=k["k_00"][-1] / k["k_11"][-1]))
    with open(HERE / "weakest_section.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    by = defaultdict(list)
    for r in rows:
        by[r["phi"]].append(r)
    g = lambda v, key: np.array([r[key] for r in v])
    L_ = ["30 d, -20 C, five packings per porosity. Conduction in x crosses vertical cuts; in y, horizontal cuts.",
          "Mean over packings (smallest of the five in brackets).", "",
          "phi     x: least ice on a cut   fewest bridges on any cut   mean bridges | y: least ice   fewest bridges   mean bridges | k_xx/k_yy, each packing"]
    for p, v in sorted(by.items()):
        L_.append(f"{p:.3f}      {g(v,'x_min_ice').mean():.2f} ({g(v,'x_min_ice').min():.2f})            "
                  f"{g(v,'x_bridges_fewest').mean():4.1f} ({g(v,'x_bridges_fewest').min()})                {g(v,'x_bridges_mean').mean():4.1f}     |   "
                  f"{g(v,'y_min_ice').mean():.2f} ({g(v,'y_min_ice').min():.2f})     {g(v,'y_bridges_fewest').mean():4.1f} ({g(v,'y_bridges_fewest').min()})         "
                  f"{g(v,'y_bridges_mean').mean():4.1f}     | " + "  ".join(f"{q:.2f}" for q in g(v, 'kxx_over_kyy')))
    txt = "\n".join(L_)
    print(txt)
    (HERE / "weakest_section.txt").write_text(txt + "\n")


if __name__ == "__main__":
    main()
