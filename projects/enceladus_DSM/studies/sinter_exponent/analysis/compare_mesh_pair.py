#!/usr/bin/env python3
"""Neck width vs time: both mesh arms and the Molaro -20 C data, on one axes.

WHY THE SIMULATION CLOCKS ARE SHIFTED
-------------------------------------
Molaro et al. do not know when their grains were brought into contact. Their
t = 0 is the first frame they measured, at which the neck was ALREADY 32.81 um
wide -- so their clock zero corresponds to some t > 0 in a simulation that
starts from true tangency. Overlaying raw clocks would compare our first
moments against their mature neck.

Each simulation's clock is therefore shifted so that t = 0 falls where ITS neck
width equals their first measured width. That is the only defensible common
origin available, it fits nothing, and it is the convention Molaro et al.
themselves used in their Fig. 12.

    coarse (eps 0.87 um) reaches 32.81 um at  7.74 h  -> shifted by -7.74 h
    fine   (eps 0.24 um) reaches 32.81 um at 15.05 h  -> shifted by -15.05 h

THE DOTTED CURVES
-----------------
Best fit w = C*(t + t0)^a to each simulation, fitted on the arm's OWN raw clock
over the range being compared (w >= the anchor width), then re-rendered on the
shifted axis. Fitting on the raw clock and shifting afterwards keeps the fit a
description of the simulation's own trajectory rather than of the shift.

--axisym IS REQUIRED WHEN BUILDING THE INPUT CSVs
-------------------------------------------------
These runs are axisymmetric r-z: the grid's y is the radius and the axis is
y = 0, so the chord neck_width.py measures runs from the axis to the contour
and IS the neck radius. Without --axisym every width is half its true value.
Rebuild with:

    python postprocess/neck_width.py <run_dir> --axisym
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "postprocess"))
from fit_neck_growth import fit_d_free  # noqa: E402

COLORS = {"coarse": "#0072B2", "fine": "#009E73", "data": "#D55E00"}
INK, GRID = "#1a1a1a", "#d8d8d8"


def load_model(run_dir):
    d = np.loadtxt(Path(run_dir) / "neck_width.csv", delimiter=",", skiprows=1, ndmin=2)
    return d[:, 0], d[:, 1] * 1e6                       # s, um WIDTH


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch", type=Path, default=Path(
        "~/SimulationResults/HPC_results/enceladus_DSM/"
        "batch_2026-08-11__17.24.16_mesh_pair").expanduser())
    ap.add_argument("--data", type=Path,
                    default=Path("inputs/validation/molaro2019_fig11_T-20.csv"))
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--prefix", default="meshpair_compare")
    args = ap.parse_args()

    d = np.loadtxt(args.data, delimiter=",", comments="#", ndmin=2)
    td, wd, ep, em = d[:, 0] * 60.0, d[:, 1], d[:, 2], d[:, 3]
    w_a = float(wd[0])                                  # their first measurement

    arms = {}
    for arm, eps in (("coarse", 0.867), ("fine", 0.240)):
        hits = sorted(args.batch.glob(f"*_a1e-3_{arm}"))
        if not hits:
            sys.exit(f"no {arm} run under {args.batch}")
        t, w = load_model(hits[0])
        if w.max() < w_a:
            sys.exit(f"{arm} arm never reaches {w_a:.2f} um (max {w.max():.2f}) — "
                     f"was its CSV built with --axisym?")
        t_a = float(np.interp(w_a, w, t))
        k = (t > t_a) & (w >= w_a)
        f = fit_d_free(t[k], w[k])
        arms[arm] = {"t": t, "w": w, "t_a": t_a, "eps": eps, "k": k, "fit": f}
        print(f"{arm:<7} eps={eps:.3f} um  {w[0]:6.2f} -> {w[-1]:6.2f} um  "
              f"reaches {w_a:.2f} um at {t_a/3600:5.2f} h  "
              f"a={f['a']:.4f}+-{f['a_ci']:.4f}  R2={f['r2']:.4f}  n={k.sum()}")

    fd = fit_d_free(td[td > 0], wd[td > 0])
    print(f"{'Molaro':<7} {'':<14} {wd[0]:6.2f} -> {wd[-1]:6.2f} um  "
          f"{'':<22} a={fd['a']:.4f}+-{fd['a_ci']:.4f}  R2={fd['r2']:.4f}  n={(td>0).sum()}")

    args.out.mkdir(parents=True, exist_ok=True)
    with open(args.out / f"{args.prefix}_fits.csv", "w", newline="") as fh:
        wtr = csv.writer(fh, lineterminator="\n")
        wtr.writerow(["series", "eps_um", "a", "a_ci95", "C", "t0_s", "r2",
                      "n_pts", "t_anchor_s", "anchor_w_um"])
        for arm, A in arms.items():
            f = A["fit"]
            wtr.writerow([arm, A["eps"], f["a"], f["a_ci"], f["C"], f["t0"],
                          f["r2"], int(A["k"].sum()), A["t_a"], w_a])
        wtr.writerow(["Molaro 2019 (-20C)", "", fd["a"], fd["a_ci"], fd["C"],
                      fd["t0"], fd["r2"], int((td > 0).sum()), 0.0, w_a])

    # ---- figure -----------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.8, 5.4))

    for arm, A in arms.items():
        t, w, t_a, f = A["t"], A["w"], A["t_a"], A["fit"]
        k = A["k"]
        ax.plot(t[k] - t_a, w[k], "-", color=COLORS[arm], lw=2.0, zorder=4,
                label=f"{arm} mesh, $\\epsilon$ = {A['eps']:.2f} µm")
        # Fit lives on the raw clock; only its rendering is shifted.
        tg = np.geomspace(t[k].min(), t[k].max(), 200)
        ax.plot(tg - t_a, f["C"] * (tg + f["t0"]) ** f["a"], ":", color=COLORS[arm],
                lw=1.6, zorder=5,
                label=f"    fit  $a$ = {f['a']:.3f} $\\pm$ {f['a_ci']:.3f}")

    ax.errorbar(td[td > 0], wd[td > 0], yerr=[em[td > 0], ep[td > 0]], fmt="o",
                color=COLORS["data"], ms=6, lw=1.2, capsize=3, zorder=6,
                label=f"Molaro 2019 ($-20\\,^\\circ$C), $a$ = {fd['a']:.3f}")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(f"time since neck width = {w_a:.2f} µm  [s]")
    ax.set_ylabel("neck width [µm]")
    ax.set_title("Neck growth, clocks aligned on the experiment's first "
                 "measured neck", fontsize=11)
    ax.grid(alpha=0.25, lw=0.5, which="both", color=GRID)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(fontsize=8.5, frameon=False, loc="upper left")

    fig.tight_layout()
    p = args.out / f"{args.prefix}.png"
    fig.savefig(p, dpi=150)
    print(f"\nplot -> {p}")


if __name__ == "__main__":
    main()
