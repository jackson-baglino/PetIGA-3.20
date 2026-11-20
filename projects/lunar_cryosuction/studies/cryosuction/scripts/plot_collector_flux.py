#!/usr/bin/env python3
"""Water vapor flux out through the collector (the x = 0 face) of type 7 runs, against time.

(a) J_c(t) = J_vap_x0 from bflux.dat (outward, kg/s per m of depth in 2D), shown in
    µg m^-1 day^-1 after a running mean over --window days (the per-step series is noisy);
    > 0: water leaves the pore through the collector, < 0: the collector feeds the pore.
(b) harvested fraction: the time integral of J_c (trapezoid on the raw series) divided by the
    initial ice mass M0 = rho_ice * ice area at t = 0, from water_budget.csv of the run or of
    --m0-from (a run with the same initial layout) when the run has none.
Checked with figcheck before saving.

    venv_cryosuction/bin/python studies/cryosuction/scripts/plot_collector_flux.py \
        --runs RUN1 RUN2 ... --labels "C1 iso" ... --out flux.png [--window 1.0] [--m0-from RUN]
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "postprocess"))

DAY = 86400.0


def m0_of(run: Path):
    f = run / "water_budget.csv"
    if not f.exists():
        return None
    with f.open() as fh:
        return float(next(csv.DictReader(fh))["ice_kg_per_m"])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--runs", nargs="+", required=True)
    ap.add_argument("--labels", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--window", type=float, default=1.0, help="running-mean window [d]")
    ap.add_argument("--m0-from", default=None, help="run whose water_budget.csv gives M0 when a run has none")
    ap.add_argument("--legend", choices=("outside", "inside"), default="outside",
                    help="legend to the right of the plots, or inside the top panel (upper right)")
    ap.add_argument("--width", type=float, default=6.3)
    ap.add_argument("--height", type=float, default=4.7)
    a = ap.parse_args(argv)
    if len(a.labels) != len(a.runs):
        sys.exit("plot_collector_flux: one label per run")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9})
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(a.width, a.height), sharex=True)
    m0_ref = m0_of(Path(a.m0_from)) if a.m0_from else None
    for k, (run, lab) in enumerate(zip(a.runs, a.labels)):
        run = Path(run)
        d = np.loadtxt(run / "bflux.dat", comments="#")
        t, J = d[:, 0], d[:, 2]
        m0 = m0_of(run) or m0_ref
        if m0 is None:
            sys.exit(f"plot_collector_flux: no M0 for {run.name}; pass --m0-from")
        td = t / DAY
        Js = np.array([J[(td >= x - a.window / 2) & (td <= x + a.window / 2)].mean() for x in td])
        cum = np.concatenate([[0.0], np.cumsum(0.5 * (J[1:] + J[:-1]) * np.diff(t))])
        c = plt.cm.tab10(k % 10)
        ax1.plot(td, Js * 1e9 * DAY, color=c, lw=1.1, label=lab)
        ax2.plot(td, 100 * cum / m0, color=c, lw=1.1)
        print(f"{lab:32s} M0 {m0:.4e} kg/m, {td[-1]:6.1f} d, harvested {100 * cum[-1] / m0:+7.1f} %, "
              f"peak flux {Js.max() * 1e9 * DAY:8.2f} µg/m/d   {run.name}")
    ax1.axhline(0, color="0.5", lw=0.6)
    ax2.axhline(0, color="0.5", lw=0.6)
    ax1.set_ylabel("collector flux\n[µg m$^{-1}$ d$^{-1}$]")
    ax2.set_ylabel("harvested\n[% of initial ice]")
    ax2.set_xlabel("time [days]")
    if a.legend == "inside":
        ax1.legend(fontsize=8, frameon=False, loc="upper right")
    else:
        ax1.legend(fontsize=8, frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    fig.tight_layout()
    from figcheck import check_figure
    for problem in check_figure(fig):
        print(f"  figcheck: {problem}")
    fig.savefig(a.out, dpi=220)
    print(f"wrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
