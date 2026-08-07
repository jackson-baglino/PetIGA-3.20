#!/usr/bin/env python3
"""Collect and plot the REV domain-size sweep.

Reads each window's k_eff.csv plus its packing metadata, and answers: at what
domain size does k_eff stop changing?

THE TWO WAYS TO MISREAD THIS PLOT, both guarded against below:

 1. Porosity drift. The measured porosity of a centred window falls from 0.414
    at 0.5 mm to 0.324 at 3 mm, converging on the master's 0.325. k_eff responds
    far more strongly to porosity than to L, so a k_eff(L) trend that merely
    tracks that drift is NOT a size effect. Panel 2 plots both on the same
    normalised axis so the two can be told apart.
 2. Stopping on a run of three. The criterion is three consecutive increases in L
    that change k_eff by less than a tolerance -- but a jump coinciding with the
    solid starting to percolate is a connectivity transition, not convergence,
    and a run of three must not span one. Panel 3 carries the percolation
    fractions for that reason.

Usage:
    python3 plot_rev.py --results <RESULTS_BASE> --tag rev [--out DIR]
"""
import argparse
import csv
import glob
import json
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]
THEMES = {
    "light": dict(surface="#fcfcfb", ink="#0b0b0b", ink2="#52514e",
                  muted="#8a8985", grid="#e3e2dd"),
    "dark":  dict(surface="#1a1a19", ink="#ffffff", ink2="#c3c2b7",
                  muted="#8a8985", grid="#333330"),
}


def style(ax, t):
    ax.set_facecolor(t["surface"])
    ax.grid(True, which="both", color=t["grid"], linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(t["grid"])
    ax.tick_params(colors=t["ink2"], labelsize=9)
    ax.xaxis.label.set_color(t["ink2"]); ax.yaxis.label.set_color(t["ink2"])
    ax.title.set_color(t["ink"])


def collect(results, tag):
    rows = []
    for geo in sorted(glob.glob(os.path.join(results, "rev_L*_T-20_rev"))):
        m = re.search(r"rev_L([0-9.]+)mm", os.path.basename(geo))
        if not m:
            continue
        L = float(m.group(1))
        # NOT f"*_{tag}": on the cluster the run directory is
        #   <ts>_<experiment>[_<tag>][_job<SLURM_JOB_ID>]
        # so an exact-suffix glob misses every job that carries an ID.
        runs = sorted(glob.glob(os.path.join(geo, f"*_{tag}")) +
                      glob.glob(os.path.join(geo, f"*_{tag}_job*")))
        if not runs:
            continue
        run = runs[-1]                              # newest matching run
        csvp = os.path.join(run, "k_eff.csv")
        if not os.path.isfile(csvp):
            print(f"  no k_eff.csv in {run}")
            continue
        d = list(csv.DictReader(open(csvp)))
        if not d:
            continue
        pack = os.path.join(PROJ, "inputs", "packings_rev", f"rev_L{m.group(1)}mm",
                            "metadata.json")
        meta = json.load(open(pack)) if os.path.isfile(pack) else {}
        rows.append({
            "L": L, "run": run,
            "t":  np.array([float(r["time"]) for r in d]),
            "k":  np.array([float(r["k_iso"]) for r in d]),
            "k00": np.array([float(r["k_00"]) for r in d]),
            "k11": np.array([float(r["k_11"]) for r in d]),
            "phi_pore": 1.0 - float(d[-1]["phi_bar"]),
            "phi_pack": meta.get("porosity_achieved", np.nan),
            "solid": meta.get("solid_largest_cluster_frac", np.nan),
            "pore":  meta.get("pore_largest_cluster_frac", np.nan),
        })
    rows.sort(key=lambda r: r["L"])
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True,
                    help="the directory holding the per-geometry run folders. "
                         "On the cluster that is $SCRATCH/enceladus_DSM; locally "
                         "it is /Users/<you>/SimulationResults/enceladus_DSM/"
                         "scratch. NOTE there is no RESULTS_BASE in your shell -- "
                         "the run scripts set it internally and never export it.")
    ap.add_argument("--tag", default="rev")
    ap.add_argument("--out", default=HERE)
    ap.add_argument("--tol", type=float, default=0.02,
                    help="relative change in k_iso counted as 'no longer changing'")
    args = ap.parse_args()

    rows = collect(args.results, args.tag)
    if not rows:
        raise SystemExit(f"no REV runs found under {args.results} with tag '{args.tag}'")

    L = np.array([r["L"] for r in rows])
    kf = np.array([r["k"][-1] for r in rows])       # k_iso at t_final
    phi = np.array([r["phi_pack"] for r in rows])

    print(f"{'L [mm]':>7} {'phi':>7} {'k_iso(t_end)':>13} {'d vs prev':>10} "
          f"{'solid%':>7} {'pore%':>6}")
    for i, r in enumerate(rows):
        d = "" if i == 0 else f"{100*(kf[i]-kf[i-1])/kf[i-1]:+9.2f}%"
        print(f"{r['L']:7.2f} {r['phi_pack']:7.4f} {kf[i]:13.6e} {d:>10} "
              f"{r['solid']*100:6.1f}% {r['pore']*100:5.1f}%")

    # Convergence: three consecutive increases in L each moving k by < tol,
    # and not spanning a change in solid percolation.
    rel = np.abs(np.diff(kf) / kf[:-1])
    verdict = "NOT REACHED within the sizes run"
    for i in range(len(rel) - 2):
        if np.all(rel[i:i + 3] < args.tol):
            verdict = f"reached at L >= {L[i]:.2f} mm"
            break
    print(f"\nREV ({args.tol*100:.0f}% criterion, 3 consecutive): {verdict}")
    print("Check the porosity column before believing it: k_eff tracks porosity "
          "much more strongly than L.")

    with open(os.path.join(args.out, "rev_summary.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["L_mm", "porosity", "k_iso_final", "k00_final", "k11_final",
                    "solid_largest_frac", "pore_largest_frac"])
        for r in rows:
            w.writerow([r["L"], r["phi_pack"], r["k"][-1], r["k00"][-1],
                        r["k11"][-1], r["solid"], r["pore"]])

    for th in ("light", "dark"):
        t = THEMES[th]
        fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8))
        fig.patch.set_facecolor(t["surface"])

        ax = axes[0]; style(ax, t)
        for r in rows:
            ax.plot(r["t"] / 3600.0, r["k"], "-", lw=1.8,
                    label=f"L = {r['L']:g} mm")
        ax.set_xlabel("time [h]"); ax.set_ylabel(r"$k_{iso}$ [W m$^{-1}$K$^{-1}$]")
        ax.set_title("1.  k$_{eff}$ trajectory per window size", fontsize=10.5, loc="left")
        ax.legend(frameon=False, fontsize=8, labelcolor=t["ink2"])

        ax = axes[1]; style(ax, t)
        ax.plot(L, kf / kf[-1], "o-", color=SERIES[0], lw=2, ms=8,
                label=r"$k_{iso}(L)$, normalised")
        ax.plot(L, phi / phi[-1], "s--", color=SERIES[1], lw=1.8, ms=7,
                label=r"window porosity, normalised")
        ax.axhline(1.0, color=t["muted"], ls=":", lw=1.2)
        ax.set_xlabel("window size L [mm]")
        ax.set_ylabel("value / value at largest L")
        ax.set_title("2.  Size effect vs porosity drift\n"
                     "if they track, it is not a size effect", fontsize=10.5, loc="left")
        ax.legend(frameon=False, fontsize=8.5, labelcolor=t["ink2"])

        ax = axes[2]; style(ax, t)
        ax.plot(L, [r["solid"] * 100 for r in rows], "o-", color=SERIES[2],
                lw=2, ms=8, label="solid in largest cluster")
        ax.plot(L, [r["pore"] * 100 for r in rows], "s--", color=SERIES[1],
                lw=1.8, ms=7, label="pore in largest cluster")
        ax.set_xlabel("window size L [mm]"); ax.set_ylabel("% in largest cluster")
        ax.set_title("3.  Connectivity — a run of three must not\n"
                     "span a percolation change", fontsize=10.5, loc="left")
        ax.legend(frameon=False, fontsize=8.5, labelcolor=t["ink2"])

        fig.suptitle("REV domain-size sweep — nested centred windows of one 3.5 mm master, "
                     "1 day at $-20\\,^\\circ$C",
                     fontsize=13, color=t["ink"], x=0.008, ha="left", y=0.99)
        fig.tight_layout(rect=[0, 0, 1, 0.90])
        p = os.path.join(args.out, f"rev_sweep_{th}.png")
        fig.savefig(p, dpi=160, facecolor=fig.get_facecolor()); plt.close(fig)
        print("wrote", p)


if __name__ == "__main__":
    main()
