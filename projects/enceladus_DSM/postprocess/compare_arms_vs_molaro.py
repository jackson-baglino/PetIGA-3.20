#!/usr/bin/env python3
"""
compare_arms_vs_molaro.py — put several finished arms side by side against
Molaro et al. (2019), on the anchored clock, and score BOTH observables.

WHY THIS EXISTS SEPARATELY FROM plot_neck_vs_molaro.py
------------------------------------------------------
That script plots one run's neck. The question this campaign actually asks is a
trade-off between two observables that pull in opposite directions, and it can
only be seen with the arms on one pair of axes:

  * neck growth is driven by the INTERNAL curvature difference between the
    neck and the grain surfaces;
  * grain shrinkage is driven by the EXTERNAL wall undersaturation.

They are not independent, because both are fed by the same vapour. A wall that
is undersaturated enough to reproduce the observed grain recession competes with
the neck for that vapour and slows it; a wall saturated enough to let the neck
grow leaves the grains too fat. Scoring only the neck hides half the result.

THE WINDOW, WHICH IS EASY TO GET WRONG
--------------------------------------
Molaro's numbers are over 78 MINUTES: neck 32.81 -> 64.78 um, large grain
-2.93 % (least-squares slope over their nine points; their caption rounds it to
-3 %). Our runs are 120 min and start from a different neck, so BOTH quantities
must be read over [t*, t*+78 min], where t* is the model time at which the neck
first reaches their first measurement.

postprocess/run_batch_measure.sh anchors the neck that way but reports
`dR_large_pct` over the WHOLE run. Because R_large(t) is linear to R^2 = 1.0000
at fixed humidity, that overstates the shrinkage by exactly 120/78 = 1.54x, and
a humidity fitted against it lands 1.54x too saturated. This script reports the
anchored window, and prints the full-run number beside it so the two are never
confused again.

Usage
-----
  python postprocess/compare_arms_vs_molaro.py <run_dir> [<run_dir> ...] \\
      [--labels "a,b,c"] [--outdir DIR] [--data CSV]

Each <run_dir> needs neck_width.csv and grain_shrinkage.csv (produced by
postprocess/neck_width.py --axisym and postprocess/grain_shrinkage.py).
Writes summary.csv and comparison.png into --outdir (default: cwd).
"""

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_REPO = Path(__file__).resolve().parent.parent
DEFAULT_DATA = _REPO / "inputs/validation/molaro2019_fig11_T-20.csv"
ANCHOR = 32.81e-6        # their first measured neck WIDTH
WINDOW = 78 * 60.0       # their record length [s]

# Colour-blind-safe, and deliberately not a rainbow: the arms are an ordered
# family (how much tuning), so lightness carries that order.
COLORS = ["#B23A48", "#2E86AB", "#1B4965", "#E1A33B", "#5B8C5A", "#8E6C99"]


def read_opts(run: Path, key: str, default: str = "") -> str:
    for f in sorted(run.glob("*.opts")):
        for line in f.read_text().splitlines():
            p = line.split("#", 1)[0].split()
            if len(p) >= 2 and p[0] == key:
                return p[1]
    return default


def load_arm(run: Path) -> dict:
    nk = list(csv.DictReader((run / "neck_width.csv").open()))
    gs = list(csv.DictReader((run / "grain_shrinkage.csv").open()))
    t = np.array([float(r["t_s"]) for r in nk])
    w = np.array([float(r["neck_width_m"]) for r in nk])
    gt = np.array([float(r["t_s"]) for r in gs])
    Rl = np.array([float(r["R_large_m"]) for r in gs])
    Rs = np.array([float(r["R_small_m"]) for r in gs])

    if w.max() < ANCHOR:
        raise SystemExit(f"{run.name}: neck never reaches the {ANCHOR*1e6:.2f} um "
                         f"anchor (max {w.max()*1e6:.2f} um) — cannot be compared.")
    t_star = float(np.interp(ANCHOR, w, t))
    t_end = t_star + WINDOW
    if t[-1] < t_end:
        print(f"  ! {run.name}: run ends {(t_end-t[-1])/60:.1f} min before "
              f"t*+78 min; its 78-min numbers are EXTRAPOLATED.")

    R0, R1 = float(np.interp(t_star, gt, Rl)), float(np.interp(t_end, gt, Rl))
    S0, S1 = float(np.interp(t_star, gt, Rs)), float(np.interp(t_end, gt, Rs))
    return dict(
        run=run, name=run.name, t=t, w=w, gt=gt, Rl=Rl, Rs=Rs,
        t_star=t_star, t_end=t_end,
        w78=float(np.interp(t_end, t, w)),
        dRl=100.0 * (R1 / R0 - 1.0), dRs=100.0 * (S1 / S0 - 1.0),
        dRl_full=100.0 * (Rl[-1] / Rl[0] - 1.0),
        humidity=read_opts(run, "-humidity"), dif_vap=read_opts(run, "-dif_vap", "2.178e-05"),
        mob_scale=read_opts(run, "-mob_scale", "1"), alph_scale=read_opts(run, "-alph_scale", "1"),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dirs", type=Path, nargs="+")
    # '|' rather than ',' because arm labels routinely contain a comma
    # ("M_0 x5, alph_sub /100") and splitting on that silently miscounts.
    ap.add_argument("--labels", default=None,
                    help="'|'-separated (or comma, if none contain one), one per run")
    ap.add_argument("--data", type=Path, default=DEFAULT_DATA)
    ap.add_argument("--outdir", type=Path, default=Path("."))
    args = ap.parse_args()

    rows = [l.split(",") for l in args.data.read_text().splitlines()
            if l.strip() and not l.startswith("#")]
    m_t = np.array([float(r[0]) for r in rows]) * 60.0
    m_w = np.array([float(r[1]) for r in rows]) * 1e-6
    m_ep = np.array([float(r[2]) for r in rows]) * 1e-6
    m_em = np.array([float(r[3]) for r in rows]) * 1e-6
    m_lg = np.array([float(r[4]) for r in rows]) * 1e-6 / 2.0

    # Their large grain's least-squares trend is the shrinkage target: only it
    # clears the measurement noise (S/N 4.4; the small grain's is 1.1).
    c = np.polyfit(m_t, m_lg, 1)
    dR_target = 100.0 * c[0] * WINDOW / np.polyval(c, 0.0)

    arms = [load_arm(d) for d in args.run_dirs]
    if args.labels:
        sep = "|" if "|" in args.labels else ","
        labels = [s.strip() for s in args.labels.split(sep)]
    else:
        labels = [a["name"].split("__")[-1] for a in arms]
    if len(labels) != len(arms):
        raise SystemExit("--labels count does not match the number of runs")

    args.outdir.mkdir(parents=True, exist_ok=True)
    with (args.outdir / "summary.csv").open("w", newline="") as fh:
        wri = csv.writer(fh, lineterminator="\n")
        wri.writerow(["label", "run", "humidity", "dif_vap", "mob_scale", "alph_scale",
                      "t_star_s", "neck_w_at_78min_um", "neck_err_um",
                      "dR_large_78min_pct", "dR_small_78min_pct", "dR_large_fullrun_pct",
                      "neck_rms_um", "neck_chi"])
        for lab, a in zip(labels, arms):
            mod = np.interp(m_t + a["t_star"], a["t"], a["w"])
            resid = mod - m_w
            a["rms"] = float(np.sqrt(np.mean(resid ** 2))) * 1e6
            sig = np.where(resid > 0, m_ep, m_em)
            a["chi"] = float(np.sqrt(np.mean((resid / sig) ** 2)))
            a["mod"] = mod
            wri.writerow([lab, a["name"], a["humidity"], a["dif_vap"],
                          a["mob_scale"], a["alph_scale"],
                          f"{a['t_star']:.1f}", f"{a['w78']*1e6:.3f}",
                          f"{(a['w78']-m_w[-1])*1e6:+.3f}",
                          f"{a['dRl']:.4f}", f"{a['dRs']:.4f}", f"{a['dRl_full']:.4f}",
                          f"{a['rms']:.3f}", f"{a['chi']:.3f}"])

    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.7))
    fig.subplots_adjust(wspace=0.30, left=0.06, right=0.985, top=0.88, bottom=0.14)

    # --- A. neck width, anchored ------------------------------------------
    ax = axes[0]
    ax.errorbar(m_t / 60.0, m_w * 1e6, yerr=[m_em * 1e6, m_ep * 1e6], fmt="ko",
                ms=5, capsize=3, lw=1.2, zorder=5, label="Molaro et al. (2019)")
    for i, (lab, a) in enumerate(zip(labels, arms)):
        sel = (a["t"] >= a["t_star"]) & (a["t"] <= a["t_end"])
        ax.plot((a["t"][sel] - a["t_star"]) / 60.0, a["w"][sel] * 1e6,
                color=COLORS[i % len(COLORS)], lw=2.0, label=lab)
    ax.set_xlabel("t − t*  [min]")
    ax.set_ylabel("neck width  [µm]")
    ax.set_title("A. Neck growth, anchored at 32.81 µm", fontsize=11, loc="left")
    ax.legend(fontsize=8.5, frameon=False, loc="lower right")
    ax.grid(alpha=0.25)

    # --- B. large-grain recession, anchored --------------------------------
    ax = axes[1]
    ax.plot(m_t / 60.0, 100.0 * (m_lg / m_lg[0] - 1.0), "ko", ms=5,
            label="Molaro, large grain")
    ax.plot(m_t / 60.0, 100.0 * (np.polyval(c, m_t) / np.polyval(c, 0.0) - 1.0),
            "k--", lw=1.3, label=f"their LS fit ({dR_target:+.2f} % / 78 min)")
    for i, (lab, a) in enumerate(zip(labels, arms)):
        sel = (a["gt"] >= a["t_star"]) & (a["gt"] <= a["t_end"])
        R0 = float(np.interp(a["t_star"], a["gt"], a["Rl"]))
        ax.plot((a["gt"][sel] - a["t_star"]) / 60.0,
                100.0 * (a["Rl"][sel] / R0 - 1.0),
                color=COLORS[i % len(COLORS)], lw=2.0, label=lab)
    ax.axhline(0, color="0.6", lw=0.8)
    ax.set_xlabel("t − t*  [min]")
    ax.set_ylabel("large-grain radius change  [%]")
    ax.set_title("B. Grain recession, same window", fontsize=11, loc="left")
    ax.legend(fontsize=8.5, frameon=False, loc="lower left")
    ax.grid(alpha=0.25)

    # --- C. the trade-off --------------------------------------------------
    # The point of the whole campaign: no arm sits at the origin.
    ax = axes[2]
    for i, (lab, a) in enumerate(zip(labels, arms)):
        ax.scatter((a["w78"] - m_w[-1]) * 1e6, a["dRl"] - dR_target,
                   s=110, color=COLORS[i % len(COLORS)], zorder=4,
                   edgecolor="white", linewidth=1.2, label=lab)
    ax.scatter([0], [0], marker="*", s=320, color="k", zorder=5, label="Molaro")
    ax.axhline(0, color="0.6", lw=0.9)
    ax.axvline(0, color="0.6", lw=0.9)
    ax.set_xlabel("neck error at t*+78 min  [µm]")
    ax.set_ylabel("shrinkage error at t*+78 min  [% points]")
    ax.set_title("C. The trade-off: neck vs shrinkage", fontsize=11, loc="left")
    ax.legend(fontsize=8.5, frameon=False, loc="best")
    ax.grid(alpha=0.25)

    fig.suptitle("Molaro et al. (2019) T = −20 °C — three tuning options, "
                 "both observables on the anchored 78-minute window",
                 fontsize=12.5, y=0.975)
    out = args.outdir / "comparison.png"
    fig.savefig(out, dpi=200)
    print(f"wrote {args.outdir/'summary.csv'}")
    print(f"wrote {out}")

    w = max(len(l) for l in labels)
    print()
    print(f"{'arm':{w}s} {'t*[s]':>6s} {'neck@78':>8s} {'err':>7s} "
          f"{'dR@78':>7s} {'err':>7s} {'dR full':>8s} {'RMS':>6s} {'chi':>5s}")
    for lab, a in zip(labels, arms):
        print(f"{lab:{w}s} {a['t_star']:6.0f} {a['w78']*1e6:8.2f} "
              f"{(a['w78']-m_w[-1])*1e6:+7.2f} {a['dRl']:+7.2f} "
              f"{a['dRl']-dR_target:+7.2f} {a['dRl_full']:+8.2f} "
              f"{a['rms']:6.2f} {a['chi']:5.2f}")
    print(f"{'TARGET':{w}s} {'':6s} {m_w[-1]*1e6:8.2f} {0.0:+7.2f} {dR_target:+7.2f} {0.0:+7.2f}")


if __name__ == "__main__":
    main()
