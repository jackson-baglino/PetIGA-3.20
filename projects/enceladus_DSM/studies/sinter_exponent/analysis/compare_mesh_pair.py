#!/usr/bin/env python3
"""Coarse vs fine arm of the mesh-resolution pair, against the Molaro -20 C data.

THE QUESTION THE PAIR WAS BUILT TO ANSWER
-----------------------------------------
Same grains, same tangent contact, same alpha_c = 1e-3, same 79 h. Only eps and
the mesh it implies differ (0.87 um / 631x198 vs 0.24 um / 2281x713). If both
arms recover the same neck curve where they overlap, coarse meshes are adequate
and the campaign gets ~13x cheaper. If not, neck VISIBILITY is binding.

WHAT THIS SCRIPT DRAWS, AND WHY IT IS THREE PANELS AND NOT ONE
--------------------------------------------------------------
(a) Raw clocks, log-log. Everything, unshifted, with each arm's resolution
    floor sqrt(12*eps/R) drawn in its own colour. This is the panel that shows
    the two facts no amount of shifting can remove: the arms are far apart, and
    the fine arm's whole run stops BELOW the experiment's first measurement.

(b) The two arms re-zeroed at a neck width they BOTH reach. If the coarse arm
    were merely "ahead" -- starting with a bigger eps-scale artefact neck on the
    same underlying trajectory -- this shift would collapse them. It does not,
    and that is the convergence verdict.

(c) The requested comparison: clocks shifted so the neck equals the
    experiment's first measured width at t = 0 (Molaro et al.'s own Fig. 12
    convention). Only the coarse arm can appear here, because only the coarse
    arm ever reaches that width. The fine arm's extrapolated arrival time is
    annotated rather than plotted, because a curve drawn from an extrapolation
    reads as data.

WHY THE FINE ARM CANNOT BE ALIGNED TO THE EXPERIMENT
----------------------------------------------------
Molaro's series spans neck widths 32.81-64.78 um. The fine arm reaches 25.56 um
in 78.6 h. There is NO width common to the two, so any alignment between them is
extrapolation, not measurement -- roughly a factor of three in time beyond the
last computed point. Reported as a number with its assumption attached, never
as a plotted curve.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "postprocess"))
from fit_neck_growth import fit_d_free, fit_d_fixed, fit_kuczynski  # noqa: E402

COLORS = {"coarse": "#0072B2", "fine": "#009E73", "data": "#D55E00"}
INK, MUTED, GRID = "#1a1a1a", "#5c5c5c", "#d8d8d8"
DEMMENIE = (0.25, 0.34)          # alpha band, incl. their +-0.01
R0 = 8.675e-5                    # mean grain radius [m], both arms and the data


def load_model(run_dir):
    d = np.loadtxt(Path(run_dir) / "neck_width.csv", delimiter=",", skiprows=1, ndmin=2)
    return d[:, 0], d[:, 1] * 1e6                       # s, um WIDTH


def load_data(path):
    d = np.loadtxt(path, delimiter=",", comments="#", ndmin=2)
    return d[:, 0] * 60.0, d[:, 1], d[:, 2], d[:, 3]    # s, um width, +err, -err


def fits_on(t, w, label, window=None):
    """All three power-law forms on a (t, w) pair, optionally width-windowed.

    Fitted on the RADIUS (w/2), as fit_neck_growth.py does, so the exponents are
    directly comparable with that script's tables; the exponent itself is of
    course indifferent to the factor of two.
    """
    m = np.ones(len(t), bool) if window is None else (w >= window[0]) & (w <= window[1])
    m &= t > 0
    tt, rr = t[m], w[m] / 2.0
    out = []
    for f in (fit_d_free(tt, rr), fit_d_fixed(tt, rr), fit_kuczynski(tt, rr)):
        if f is not None:
            f = dict(f, series=label, n_pts=len(tt),
                     w_lo=w[m].min() if m.any() else np.nan,
                     w_hi=w[m].max() if m.any() else np.nan)
            out.append(f)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch", type=Path, default=Path(
        "~/SimulationResults/HPC_results/enceladus_DSM/"
        "batch_2026-08-11__17.24.16_mesh_pair").expanduser())
    ap.add_argument("--data", type=Path,
                    default=Path("inputs/validation/molaro2019_fig11_T-20.csv"))
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--prefix", default="meshpair_compare")
    ap.add_argument("--collapse-width", type=float, default=None,
                    help="neck width [um] at which panel (b) re-zeroes both "
                         "arms (default: midpoint of their shared range)")
    args = ap.parse_args()

    runs = {}
    for arm in ("coarse", "fine"):
        hits = sorted(args.batch.glob(f"*_a1e-3_{arm}"))
        if not hits:
            sys.exit(f"no {arm} run under {args.batch}")
        runs[arm] = load_model(hits[0])
    td, wd, ep, em = load_data(args.data)

    tc, wc = runs["coarse"]
    tf, wf = runs["fine"]

    # ---- windows ----------------------------------------------------------
    arm_lo, arm_hi = max(wc.min(), wf.min()), min(wc.max(), wf.max())
    w_col = args.collapse_width or float(np.sqrt(arm_lo * arm_hi))
    gap = wd.min() - wf.max()

    print(f"\ncoarse : {wc.min():6.2f} -> {wc.max():6.2f} um over {tc[-1]/3600:5.2f} h "
          f"({len(tc)} samples)")
    print(f"fine   : {wf.min():6.2f} -> {wf.max():6.2f} um over {tf[-1]/3600:5.2f} h "
          f"({len(tf)} samples)")
    print(f"Molaro : {wd.min():6.2f} -> {wd.max():6.2f} um over {td[-1]/60:5.2f} min "
          f"({len(td)} samples)")
    print(f"\narms share {arm_lo:.2f}-{arm_hi:.2f} um; collapse test at {w_col:.2f} um")
    print(f"fine arm vs experiment: NO shared width "
          f"(fine tops out {gap:.2f} um below the data's first point)")

    # ---- rate ratio at matched neck width ---------------------------------
    print(f"\n{'w [um]':>8} {'coarse [h]':>12} {'fine [h]':>10} {'fine/coarse':>12}")
    ratio_rows = []
    for W in np.linspace(arm_lo * 1.02, arm_hi, 7):
        a = float(np.interp(W, wc, tc)) / 3600.0
        b = float(np.interp(W, wf, tf)) / 3600.0
        print(f"{W:8.2f} {a:12.3f} {b:10.2f} {b/a:12.1f}")
        ratio_rows.append({"w_um": W, "t_coarse_h": a, "t_fine_h": b, "ratio": b / a})

    # ---- power-law fits ---------------------------------------------------
    rows = []
    rows += fits_on(tc, wc, "coarse (own range)")
    rows += fits_on(tf, wf, "fine (own range)")
    rows += fits_on(tc, wc, "coarse (shared arm window)", (arm_lo, arm_hi))
    rows += fits_on(tf, wf, "fine (shared arm window)", (arm_lo, arm_hi))
    rows += fits_on(td, wd, "Molaro (-20C)")

    print(f"\n{'series':<28} {'form':<8} {'a':>8} {'+-95%':>8} {'m':>7} {'R2':>8} "
          f"{'window [um]':>15} {'n':>4}")
    print("-" * 94)
    for r in rows:
        print(f"{r['series']:<28} {r['form']:<8} {r['a']:>8.4f} {r['a_ci']:>8.4f} "
              f"{r['m']:>7.2f} {r['r2']:>8.4f} "
              f"{r['w_lo']:>6.1f}-{r['w_hi']:<8.1f} {r['n_pts']:>4}")
    print(f"\nDemmenie et al. (2025): alpha = {DEMMENIE[0]:.2f}-{DEMMENIE[1]:.2f} (m = 3)")

    # ---- how far the fine arm is from the experiment ----------------------
    ff = next(r for r in rows if r["series"] == "fine (own range)" and r["form"] == "d_free")
    # r = C*(t+t0)^a  =>  t = (r/C)^(1/a) - t0
    t_reach = (wd.min() / 2.0 / ff["C"]) ** (1.0 / ff["a"]) - ff["t0"]
    print(f"\nfine arm extrapolated to the experiment's first width "
          f"({wd.min():.2f} um): t = {t_reach/3600:.0f} h "
          f"= {t_reach/tf[-1]:.1f}x the run's length. EXTRAPOLATION "
          f"(d_free, a = {ff['a']:.3f}) — not a measurement.")

    args.out.mkdir(parents=True, exist_ok=True)
    with open(args.out / f"{args.prefix}_fits.csv", "w", newline="") as fh:
        cols = ["series", "form", "a", "a_ci", "m", "C", "t0", "r2",
                "n_pts", "w_lo", "w_hi"]
        wtr = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore",
                             lineterminator="\n")
        wtr.writeheader()
        wtr.writerows(rows)
    with open(args.out / f"{args.prefix}_rate_ratio.csv", "w", newline="") as fh:
        wtr = csv.DictWriter(fh, fieldnames=list(ratio_rows[0].keys()),
                             lineterminator="\n")
        wtr.writeheader()
        wtr.writerows(ratio_rows)

    # ---- figure -----------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15.0, 5.0))

    # (a) raw clocks -------------------------------------------------------
    ax = axes[0]
    ax.plot(tc, wc, "-", color=COLORS["coarse"], lw=2.0, label="coarse, $\\epsilon$=0.87 µm")
    ax.plot(tf, wf, "-", color=COLORS["fine"], lw=2.0, label="fine, $\\epsilon$=0.24 µm")
    ax.errorbar(td[td > 0], wd[td > 0], yerr=[em[td > 0], ep[td > 0]], fmt="o",
                color=COLORS["data"], ms=5, lw=1.0, capsize=2.5,
                label="Molaro 2019 ($-20\\,^\\circ$C)")
    for arm, e in (("coarse", 8.6750e-07), ("fine", 2.4002e-07)):
        fl = 2 * R0 * np.sqrt(12 * e / R0) * 1e6           # floor as a WIDTH
        ax.axhline(fl, color=COLORS[arm], lw=0.9, ls=(0, (1, 3)), alpha=0.8)
        ax.annotate(f"{arm} floor {fl:.0f} µm", (tc[0], fl), fontsize=7,
                    color=COLORS[arm], va="bottom")
    ax.axhspan(wf.max(), wd.min(), color=INK, alpha=0.07, zorder=0)
    ax.annotate(f"no shared width\n({wf.max():.1f} → {wd.min():.1f} µm)",
                (tc[0] * 3, np.sqrt(wf.max() * wd.min())), fontsize=7.5,
                color=MUTED, va="center")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("time [s]"); ax.set_ylabel("neck width [µm]")
    ax.set_title("(a) raw clocks", fontsize=10)
    ax.legend(fontsize=7.5, frameon=False, loc="lower right")

    # (b) arms re-zeroed at a width both reach ------------------------------
    ax = axes[1]
    for arm, (t_, w_) in (("coarse", (tc, wc)), ("fine", (tf, wf))):
        t0_ = float(np.interp(w_col, w_, t_))
        k = t_ > t0_
        ax.plot(t_[k] - t0_, w_[k], "-", color=COLORS[arm], lw=2.0,
                label=f"{arm} (zeroed at {w_col:.1f} µm)")
    # Guide slopes are pinned to the collapse point and drawn only across the
    # span the arms actually occupy; y is clamped to the curves so the guides
    # cannot drag the axis open (a 1/3 slope over three decades reaches 300 um).
    lo_x, hi_x, hi_y = 3e2, 3e5, max(wc.max(), wf.max())
    tg = np.geomspace(lo_x, hi_x, 50)
    for a_, lab in ((1 / 3, "1/3  sublim.-cond."), (1 / 5, "1/5  bulk diff."),
                    (1 / 7, "1/7  surf. diff.")):
        yg = w_col * (tg / lo_x) ** a_
        ax.plot(tg, yg, color=MUTED, lw=0.7, ls=(0, (4, 3)), alpha=0.55, zorder=0)
        j = np.searchsorted(yg, hi_y * 0.98)
        j = min(j, len(tg) - 1)
        ax.annotate(lab, (tg[j], yg[j]), fontsize=7, color=MUTED,
                    ha="right", va="bottom")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lo_x, hi_x)
    ax.set_ylim(w_col * 0.92, hi_y * 1.12)
    ax.set_xlabel(f"time since neck = {w_col:.1f} µm [s]")
    ax.set_ylabel("neck width [µm]")
    ax.set_title("(b) arms re-zeroed at a shared neck — no collapse", fontsize=10)
    ax.legend(fontsize=7.5, frameon=False, loc="upper left")

    # (c) aligned on the experiment's first width ---------------------------
    ax = axes[2]
    w_a = float(wd.min())
    t_ac = float(np.interp(w_a, wc, tc))
    k = tc > t_ac
    ax.plot(tc[k] - t_ac, wc[k], "-", color=COLORS["coarse"], lw=2.0,
            label=f"coarse (shifted $-${t_ac/3600:.2f} h)")
    ax.errorbar(td[td > 0], wd[td > 0], yerr=[em[td > 0], ep[td > 0]], fmt="o",
                color=COLORS["data"], ms=5, lw=1.0, capsize=2.5,
                label="Molaro 2019 ($-20\\,^\\circ$C)")
    # No anchor marker here: the anchor is t = 0, which a log axis cannot show.
    # The axis label carries it instead.
    ax.axhline(w_a, color=INK, lw=0.8, ls=(0, (3, 3)), alpha=0.5, zorder=1)
    ax.annotate(f"fine arm absent: never reaches {w_a:.1f} µm\n"
                f"(extrapolates to ~{t_reach/3600:.0f} h, {t_reach/tf[-1]:.0f}× its run)",
                (0.03, 0.17), xycoords="axes fraction", fontsize=7.5, color=COLORS["fine"])
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(f"time since neck = {w_a:.1f} µm [s]")
    ax.set_ylabel("neck width [µm]")
    ax.set_title("(c) aligned on the experiment's first width", fontsize=10)
    ax.legend(fontsize=7.5, frameon=False, loc="lower right")

    for ax in axes:
        ax.grid(alpha=0.25, lw=0.5, which="both", color=GRID)
        ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    p = args.out / f"{args.prefix}.png"
    fig.savefig(p, dpi=150)
    print(f"\nplot -> {p}")
    print(f"csv  -> {args.out}/{args.prefix}_fits.csv, {args.prefix}_rate_ratio.csv\n")


if __name__ == "__main__":
    main()
