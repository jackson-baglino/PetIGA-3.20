"""
plot_neck_vs_molaro.py — model neck growth against Molaro et al. (2019), on a
common clock, with a power-law fit to each.

THE CLOCK PROBLEM
-----------------
Our run and their experiment do not share a t = 0. We start from a chosen
initial neck (r = 14 um in the current geometries); their record starts at
32.81 um WIDTH at an unknown time after the grains were placed in contact. So
raw times are not comparable and neither is any exponent fitted against them.

Both series are therefore shifted to the same anchor: t' = 0 is defined as the
moment each curve passes the SAME neck width (default 32.81 um, their first
measurement). For the model that instant is found by linear interpolation
between snapshots; for the experiment it is t = 0 by construction. After the
shift the two curves are directly comparable, which is the whole point.

THE FIT
-------
r = C*(t' + t0)^a with t0 free (fit_neck_growth.fit_d_free, Demmenie's
protocol). t0 MUST be free: at t' = 0 the neck is already 32.81 um, not zero,
so forcing t0 = 0 reports a spuriously small exponent -- the curve looks flat
because most of the growth happened before the window opens.

`a` is the exponent on the RADIUS. Widths are halved on read, so it is the same
number either way -- a is scale-invariant -- but the label says which.

Usage
-----
  python postprocess/plot_neck_vs_molaro.py <run_dir>
  python postprocess/plot_neck_vs_molaro.py <run_dir> --loglog
  python postprocess/plot_neck_vs_molaro.py <run_dir> --anchor-width 32.81e-6
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pplib  # noqa: E402
from fit_neck_growth import fit_d_free  # noqa: E402

# Validated categorical pair (OKLab/CVD checked: protan 19.7, deutan 24.1,
# normal 28.5). Colour is the ENTITY -- model vs experiment.
C_MODEL, C_DATA, C_INK, C_MUTED = "#3d74d9", "#d1495b", "#3f434a", "#8a8a8a"

DEFAULT_ANCHOR = 32.81e-6      # Molaro's first measured neck WIDTH [m]


def read_model(run_dir: Path):
    """(t_s, width_m) from neck_width.py's CSV."""
    csv_path = run_dir / "neck_width.csv"
    if not csv_path.is_file():
        sys.exit(f"{csv_path} not found — run:\n"
                 f"  python postprocess/neck_width.py {run_dir} --axisym")
    t, w = [], []
    with csv_path.open() as fh:
        for row in csv.DictReader(fh):
            t.append(float(row["t_s"]))
            w.append(float(row["neck_width_m"]))
    o = np.argsort(t)
    return np.asarray(t)[o], np.asarray(w)[o]


def read_experiment(path: Path):
    """(t_s, width_m) from the Fig. 11 validation CSV."""
    t, w = [], []
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            f = next(csv.reader([line]))
            t.append(float(f[0]) * 60.0)
            w.append(float(f[1]) * 1e-6)
    return np.asarray(t), np.asarray(w)


def anchor_time(t, w, target):
    """First time the curve reaches `target`, linearly interpolated."""
    for (t0, w0), (t1, w1) in zip(zip(t, w), zip(t[1:], w[1:])):
        if w0 <= target <= w1 and w1 > w0:
            return t0 + (target - w0) * (t1 - t0) / (w1 - w0)
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--data", type=Path, default=None,
                    help="validation CSV (default: the Molaro -20 C table)")
    ap.add_argument("--anchor-width", type=float, default=DEFAULT_ANCHOR,
                    help="neck WIDTH [m] at which both clocks are set to zero")
    ap.add_argument("--loglog", action="store_true",
                    help="log-log axes (default: linear)")
    ap.add_argument("--full-range", action="store_true",
                    help="show the whole model run, not just the window where "
                         "the experiment has data")
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    data_path = args.data or (Path(__file__).resolve().parent.parent
                              / "inputs/validation/molaro2019_fig11_T-20.csv")
    tm, wm = read_model(args.run_dir)
    td, wd = read_experiment(data_path)

    t_star = anchor_time(tm, wm, args.anchor_width)
    if t_star is None:
        print(f"  NOTE: the model neck never reaches the anchor "
              f"{args.anchor_width*1e6:.2f} um "
              f"(it spans {wm.min()*1e6:.2f}–{wm.max()*1e6:.2f} um).")
        print("  Falling back to anchoring on the model's FIRST point, so the")
        print("  curves start together but the comparison is not yet on their")
        print("  measured range. Run longer for a real comparison.")
        t_star, anchor = tm[0], wm[0]
    else:
        anchor = args.anchor_width
    td_star = anchor_time(td, wd, anchor) or td[0]

    tm_s, td_s = tm - t_star, td - td_star

    # Fit each series on its own post-anchor window, radius not width.
    fits = {}
    for key, ts, ws, col in (("model", tm_s, wm, C_MODEL),
                             ("Molaro 2019", td_s, wd, C_DATA)):
        m = ts >= 0
        fits[key] = fit_d_free(ts[m], ws[m] / 2.0) if m.sum() >= 4 else None

    print(f"  anchor width      {anchor*1e6:.2f} um")
    print(f"  model  t* = {t_star:9.1f} s   ({len(tm)} snapshots, "
          f"{wm.min()*1e6:.2f}–{wm.max()*1e6:.2f} um)")
    print(f"  Molaro t* = {td_star:9.1f} s")
    for k, f in fits.items():
        if f:
            print(f"  {k:12s} r ~ (t+t0)^a :  a = {f['a']:.4f} ± {f['a_ci']:.4f}"
                  f"   t0 = {f['t0']:8.1f} s   R2 = {f['r2']:.4f}")
        else:
            print(f"  {k:12s} fit skipped (fewer than 4 points past the anchor)")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.6, 5.0))
    tmax = max(tm_s.max(), td_s.max())

    for key, ts, ws, col, style in (
            ("model", tm_s, wm, C_MODEL, dict(ls="-", lw=2)),
            ("Molaro 2019", td_s, wd, C_DATA, dict(ls="none", marker="o",
                                                   mfc="none", mew=1.8, ms=8))):
        m = ts >= 0
        ax.plot(ts[m] / 60.0, ws[m] * 1e6, color=col, zorder=3,
                label=key, **style)
        f = fits[key]
        if f:
            tf = np.linspace(0, ts[m].max(), 300)
            # C is fitted on the RADIUS in metres: x2 for width, x1e6 for um.
            ax.plot(tf / 60.0, 2.0 * f["C"] * (tf + f["t0"]) ** f["a"] * 1e6,
                    color=col, ls="--", lw=1.4, alpha=0.85, zorder=2,
                    label=f"    fit:  a = {f['a']:.3f} ± {f['a_ci']:.3f}")

    # Default the x-range to the experiment's window: the model typically runs
    # hours past their 78 min, and letting that tail set the scale squashes the
    # part being compared into the first centimetre of the axis.
    if not args.full_range:
        span = td_s[td_s >= 0].max() if (td_s >= 0).any() else tmax
        ax.set_xlim(-0.03 * span / 60.0, 1.25 * span / 60.0)
        vis = wm[(tm_s >= 0) & (tm_s <= 1.25 * span)]
        top = max(wd.max(), vis.max() if len(vis) else 0.0)
        ax.set_ylim(0.9 * min(wd.min(), wm[tm_s >= 0].min() if (tm_s >= 0).any()
                              else wd.min()) * 1e6, 1.08 * top * 1e6)

    ax.axvline(0.0, color=C_MUTED, lw=0.8, ls=(0, (1, 3)), zorder=1)
    ax.annotate(f"clocks aligned at {anchor*1e6:.2f} $\\mu$m",
                (0.0, ax.get_ylim()[0]), textcoords="offset points",
                xytext=(6, 8), fontsize=8, color=C_INK)

    if args.loglog:
        ax.set_xscale("log")
        ax.set_yscale("log")
    ax.set_xlabel("time since both curves reached the anchor width  [min]")
    ax.set_ylabel("neck width  [$\\mu$m]")
    ax.set_title("Neck growth vs Molaro et al. (2019), $T = -20\\,^\\circ$C",
                 fontsize=11, color=C_INK)
    ax.grid(alpha=0.25, lw=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(colors=C_INK, labelsize=9)
    ax.legend(frameon=False, fontsize=9, loc="lower right")

    # The exponent is the headline: state it once, large, out of the way.
    if fits["model"] and fits["Molaro 2019"]:
        ax.text(0.03, 0.95,
                f"$a_{{\\rm model}}$ = {fits['model']['a']:.3f}"
                f"      $a_{{\\rm Molaro}}$ = {fits['Molaro 2019']['a']:.3f}",
                transform=ax.transAxes, fontsize=12, color=C_INK,
                va="top", ha="left")

    fig.tight_layout()
    out = args.out or (args.run_dir / "plots" /
                       ("neck_vs_molaro_loglog.png" if args.loglog
                        else "neck_vs_molaro.png"))
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=160)
    print(f"  wrote {out}")


if __name__ == "__main__":
    main()
