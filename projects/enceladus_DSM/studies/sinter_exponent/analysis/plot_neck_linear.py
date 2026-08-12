#!/usr/bin/env python3
"""Neck WIDTH vs time on linear axes, with straight-line fits, model vs Molaro.

WHY A LINE, WHEN fit_neck_growth.py ALREADY FITS A POWER LAW
------------------------------------------------------------
Different question. The power-law fit asks which Kuczynski mechanism the curve
is consistent with; a straight line asks the much blunter "how many microns of
neck per hour, and does the rate hold steady over the measured range". The
second question is the one you need to compare a 79 h simulation against a
78 min experiment, and it is the one an experimentalist reads off a plot.

The line is NOT a competing physical model. r ~ t^0.2 (what these curves
actually do) is not a line, and this script says so rather than hiding it: the
fit window, R^2 and max residual are printed and drawn, and the residual panel
exists so the curvature is visible instead of buried in an R^2.

THE PARAMETERISATION, AND WHY THE REQUESTED ONE IS NOT FITTED AS WRITTEN
-----------------------------------------------------------------------
    w = A*(t - t0) + C

has three parameters but only two are identifiable: it expands to
A*t + (C - A*t0), so t0 and C trade off exactly and a fitter handed all three
returns whichever pair it stumbled into, with an infinite covariance ridge.
Both identifiable specialisations are reported instead:

    slope-intercept   w = A*t + C          (C = width extrapolated to t = 0)
    zero-crossing     w = A*(t - t0)       (t0 = time the line reaches w = 0)

They are the SAME LINE written two ways -- t0 = -C/A -- so the fit is done once
and both readings are printed. For these runs t0 is a large negative number
(the line, extended backwards, has the neck vanishing hundreds of hours before
the run began) which is a fact about the line, not about the ice: it is what
"the curve is concave, so a chord under-slopes the early growth" looks like
when you extrapolate it. Quote C, not t0.

THE TIME-SCALE FACTOR
---------------------
The model at alpha_c = 1e-3 is ~2 decades slower than the Molaro cryostage, so
on a shared linear time axis one of the two series is a vertical line at the
origin. The overlay therefore stretches the experiment's clock by

    S = A_data / A_model     (fitted over the SHARED WIDTH RANGE)

which is exactly the factor that makes the two fitted lines parallel on screen.
Everything left over after that -- curvature, offset, spread -- is shape
difference, which is the only thing the overlay is meant to show. S is printed
on the figure so it can never be read as a claim that the rates agree.

Usage
-----
    python studies/sinter_exponent/analysis/plot_neck_linear.py <run_dir> \
        --out studies/sinter_exponent/results/mesh_pair
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

# Okabe-Ito, matching fit_neck_growth.py.
COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7"]
INK, MUTED, GRID = "#1a1a1a", "#5c5c5c", "#d8d8d8"


def linfit(t, w):
    """w = A*t + C by OLS, with the 95% CI on A and the zero-crossing reading.

    Returns None below 3 points -- a two-parameter fit through two points has
    no residual and would report a meaningless R^2 of 1.
    """
    if len(t) < 3:
        return None
    A, C = np.polyfit(t, w, 1)
    pred = A * t + C
    resid = w - pred
    n = len(t)
    sxx = np.sum((t - t.mean()) ** 2)
    se = np.sqrt(np.sum(resid ** 2) / (n - 2) / sxx) if n > 2 and sxx > 0 else np.nan
    ss = np.sum((w - w.mean()) ** 2)
    return {
        "A": A, "A_ci95": 1.96 * se, "C": C,
        "t0": -C / A if A else np.nan,
        "r2": 1.0 - np.sum(resid ** 2) / ss if ss > 0 else np.nan,
        "max_resid": np.abs(resid).max(), "rms_resid": np.sqrt(np.mean(resid ** 2)),
        "n": n, "t_lo": t[0], "t_hi": t[-1], "w_lo": w[0], "w_hi": w[-1],
    }


def load_model(run_dir):
    d = np.loadtxt(Path(run_dir) / "neck_width.csv", delimiter=",", skiprows=1, ndmin=2)
    return d[:, 0], d[:, 1] * 1e6                      # s, um WIDTH


def load_data(path):
    d = np.loadtxt(path, delimiter=",", comments="#", ndmin=2)
    return d[:, 0] * 60.0, d[:, 1], d[:, 2], d[:, 3]   # s, um width, +err, -err


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--data", type=Path,
                    default=Path("inputs/validation/molaro2019_fig11_T-20.csv"))
    ap.add_argument("--prefix", default="meshpair_linear")
    ap.add_argument("--label", default="model coarse eps0.87um (tangent)")
    args = ap.parse_args()

    t, w = load_model(args.run_dir)
    td, wd, ep, em = load_data(args.data)

    # The shared window is a width range, not a time cut: the two clocks are
    # incommensurable but the necks are the same physical size, so width is the
    # only common coordinate the comparison has.
    #
    # It must be the INTERSECTION, clipped on BOTH series. Clipping only the
    # model leaves the data's chord spanning 32.8-64.8 um against the model's
    # 32.8-53.5; because both curves are concave, the longer chord sits at a
    # visibly shallower slope, and the ratio of the two comes out ~1.8x too
    # small purely from the mismatched spans.
    lo, hi = max(w.min(), wd.min()), min(w.max(), wd.max())
    shared, shared_d = (w >= lo) & (w <= hi), (wd >= lo) & (wd <= hi)

    fits = {
        "model (full run)": (t, w, linfit(t, w)),
        "model (shared width)": (t[shared], w[shared], linfit(t[shared], w[shared])),
        "Molaro (full range)": (td, wd, linfit(td, wd)),
        "Molaro (shared width)": (td[shared_d], wd[shared_d],
                                  linfit(td[shared_d], wd[shared_d])),
    }

    args.out.mkdir(parents=True, exist_ok=True)
    rows = []
    print(f"\n{'series':<24} {'A [um/h]':>10} {'+-95%':>9} {'C [um]':>9} "
          f"{'t0 [h]':>10} {'R2':>7} {'rms':>7} {'max':>7} {'window [um]':>14} {'n':>4}")
    print("-" * 112)
    for name, (tt, ww, f) in fits.items():
        if f is None:
            continue
        print(f"{name:<24} {f['A']*3600:>10.4f} {f['A_ci95']*3600:>9.4f} "
              f"{f['C']:>9.2f} {f['t0']/3600:>10.1f} {f['r2']:>7.4f} "
              f"{f['rms_resid']:>7.2f} {f['max_resid']:>7.2f} "
              f"{f['w_lo']:>6.1f}-{f['w_hi']:<7.1f} {f['n']:>4}")
        rows.append({"series": name, "A_um_per_s": f["A"], "A_um_per_h": f["A"] * 3600,
                     "A_ci95_um_per_h": f["A_ci95"] * 3600, "C_um": f["C"],
                     "t0_s": f["t0"], "t0_h": f["t0"] / 3600, "r2": f["r2"],
                     "rms_resid_um": f["rms_resid"], "max_resid_um": f["max_resid"],
                     "n_pts": f["n"], "t_lo_s": f["t_lo"], "t_hi_s": f["t_hi"],
                     "w_lo_um": f["w_lo"], "w_hi_um": f["w_hi"]})

    fm, fd = fits["model (shared width)"][2], fits["Molaro (shared width)"][2]

    # The overlay stretch is the ELAPSED-TIME ratio across [lo, hi], not the
    # ratio of fitted slopes. With only 9 experimental points the data's fit
    # window cannot land on the window edge -- it truncates at whatever sample
    # precedes it (46.8 um here, not 53.5) -- and a chord over that shorter,
    # steeper sub-range overstates the ratio. The elapsed-time measure
    # interpolates to the edges, uses no fit at all, and is what "how much
    # slower is the model at equal neck size" actually means. The slope ratio
    # is printed alongside as a consistency check, not used.
    trav = [float(np.interp(hi, ww, tt)) - float(np.interp(lo, ww, tt))
            for tt, ww in ((t, w), (td, wd))]
    S = trav[0] / trav[1]
    print(f"\ntime-scale factor over the shared width range {lo:.1f}-{hi:.1f} um:")
    print(f"    elapsed-time ratio   S = {S:.1f}x   (model {trav[0]/3600:.1f} h "
          f"vs data {trav[1]/60:.1f} min)   <- used for the overlay")
    print(f"    fitted-slope ratio       {fd['A']/fm['A']:.1f}x   (check; differs "
          f"because the data's window truncates at {fd['w_hi']:.1f} um)")

    csv_out = args.out / f"{args.prefix}_fits.csv"
    with open(csv_out, "w", newline="") as fh:
        wtr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), lineterminator="\n")
        wtr.writeheader()
        wtr.writerows(rows)
    print(f"csv  -> {csv_out}")

    # ---- figure -----------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax, axr) = plt.subplots(
        2, 1, figsize=(8.2, 6.6), sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.08})

    th, tdh = t / 3600.0, td * S / 3600.0             # both on the MODEL clock, hours

    ax.axvspan(th[shared].min(), th[shared].max(), color=COLORS[0], alpha=0.06, zorder=0,
               label=f"shared width range {lo:.1f}-{hi:.1f} $\\mu$m")
    ax.plot(th, w, "-", color=COLORS[0], lw=2.0, zorder=4, label=args.label)
    ax.errorbar(tdh, wd, yerr=[em, ep], fmt="o", color=COLORS[1], ms=5, lw=1.0,
                capsize=2.5, zorder=5,
                label=f"Molaro 2019 ($-20\\,^\\circ$C), clock $\\times${S:.0f}")

    # Each line is EVALUATED on the clock it was fitted on and only then moved
    # onto the shared axis. Evaluating the stretched data times against A_data
    # would inflate that line by the full factor S.
    for f, tt, stretch, c, lab in ((fm, t[shared], 1.0, COLORS[0], "model"),
                                   (fd, td[shared_d], S, COLORS[1], "Molaro")):
        tg = np.linspace(tt.min(), tt.max(), 2)                 # own clock, s
        ax.plot(tg * stretch / 3600.0, f["A"] * tg + f["C"], "--", color=c, lw=1.3,
                alpha=0.9, zorder=6,
                label=f"{lab} fit: $A$={f['A']*3600:.3f} $\\mu$m/h, "
                      f"$C$={f['C']:.1f} $\\mu$m, $R^2$={f['r2']:.3f}")

    ax.set_ylabel("neck width [$\\mu$m]")
    ax.set_title("Neck width vs time, linear fit  $w = A\\,t + C$", fontsize=11)
    ax.grid(alpha=0.25, lw=0.5, color=GRID)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(fontsize=8, frameon=False, loc="upper left")

    # The experiment's own clock, so the 78 min is never lost behind the stretch.
    axt = ax.secondary_xaxis("top", functions=(lambda x: x * 3600 / S / 60,
                                               lambda x: x * 60 * S / 3600))
    axt.set_xlabel("Molaro elapsed time [min]", fontsize=9, color=COLORS[1])
    axt.tick_params(colors=COLORS[1], labelsize=8)

    # Residuals: the curvature the R^2 rounds away.
    axr.axhline(0, color=MUTED, lw=0.8)
    axr.plot(th[shared], w[shared] - (fm["A"] * t[shared] + fm["C"]), "-",
             color=COLORS[0], lw=1.6)
    axr.plot(tdh[shared_d], wd[shared_d] - (fd["A"] * td[shared_d] + fd["C"]), "o",
             color=COLORS[1], ms=4)
    axr.plot(tdh[~shared_d], wd[~shared_d] - (fd["A"] * td[~shared_d] + fd["C"]), "o",
             color=COLORS[1], ms=4, alpha=0.25)      # extrapolated, out of window
    # Scale to the IN-WINDOW residuals. The faint out-of-window points run to
    # -33 um (the data curve leaving the model's range entirely) and would
    # otherwise flatten the +-4 um structure this panel exists to show.
    rin = np.concatenate([w[shared] - (fm["A"] * t[shared] + fm["C"]),
                          wd[shared_d] - (fd["A"] * td[shared_d] + fd["C"])])
    pad = max(1.0, 1.3 * np.abs(rin).max())
    axr.set_ylim(-pad, pad)
    axr.set_xlabel("model elapsed time [h]")
    axr.set_ylabel("resid.\n[$\\mu$m]", fontsize=9)
    axr.grid(alpha=0.25, lw=0.5, color=GRID)
    axr.spines[["top", "right"]].set_visible(False)
    axr.tick_params(labelsize=8)

    p = args.out / f"{args.prefix}.png"
    fig.savefig(p, dpi=150, bbox_inches="tight")
    print(f"plot -> {p}\n")


if __name__ == "__main__":
    main()
