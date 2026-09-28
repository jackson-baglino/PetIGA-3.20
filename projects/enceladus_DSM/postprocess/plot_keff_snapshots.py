#!/usr/bin/env python3
"""plot_keff_snapshots.py — k_eff curve with microstructure snapshots inset.

    python3 plot_keff_snapshots.py --dir <run> [--times 1 10 20 30] [--steps ...]
        [--n-snapshots 4] [--width 7.2] [--absolute] [--iso-only]
        [--show-relaxation] [--title-time STR] [--title-ssa STR] [--no-title]
        [--xlabel-time STR] [--xlabel-ssa STR] [--ylabel STR]
        [--source vts|sol] [--save-dir DIR] [--formats pdf png]

Two page-width manuscript figures, each one k_eff curve with 3-4
microstructure snapshots inset along the bottom of the same axes, under the
curve. Each snapshot is lettered, and the same letter marks its instant on
the curve, joined to it by a leader line:

    keff_time_snapshots.{pdf,png}   k_eff / k_eff,0  vs  time [d]
    keff_ssa_snapshots.{pdf,png}    k_eff / k_eff,0  vs  SSA / SSA_0

The insets are ordered by where their markers fall on the x axis, so the
leaders never cross: left to right in time on the time figure, and RIGHT to
left in time on the SSA figure, where SSA falls as the packing sinters.

PANELS. They follow make_packing_movie.py, and they import its helpers so
the figure and the movie cannot drift apart: supersaturation
sigma = rho_v/rho_vs(T) - 1 on cmocean `balance`, re-centred so the pale
middle is sigma = 0, on an asinh scale; ice painted on top with cmocean
`ice`, transparent below phi = 0.5. The sigma range is taken over the pore
space of the SHOWN snapshots only, and shared by every panel.

NORMALIZATION. Each quantity is divided by its own value at the REFERENCE
sample: the first at t >= --baseline-days (default 1 d), as in plot_keff.py.
That is what the subscript 0 means -- k_xx by k_xx,0, k_iso by k_iso,0, SSA
by SSA_0 -- and the caption should say so. The reference time and values are
printed for it. --absolute plots k_eff in W m^-1 K^-1 on the time figure
instead; the SSA figure is always normalized on both axes.

ONLY MEASURED SAMPLES. Samples before the reference (the initial condition
relaxing) are not drawn and not annotated; that belongs in the manuscript
text. --show-relaxation draws them in grey, unlabelled. k_eff and SSA are
paired by step (plot_keff.load drops a k_eff sample with no SSA row rather
than borrowing a neighbour's), and the SSA figure plots the samples as
points, with no line joining them. Snapshots are only taken at steps that
have a k_eff sample, so every marker sits on a measured value.

SNAPSHOTS. --times (days) picks the eligible snapshot nearest each time;
--steps names them exactly. By default --n-snapshots are spaced evenly in
time from the reference to the end of the run.

OUTPUT. <dir>/plots/keff/snapshots/ (or --save-dir). PDF is the manuscript
file -- the curve is vector, the fields are embedded rasters -- and PNG is a
preview. Exits 0 without writing anything when the run has no k_eff CSV.
"""
from __future__ import annotations

import argparse
import glob
import os
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm, ListedColormap
from matplotlib.cm import ScalarMappable
from matplotlib.patches import ConnectionPatch
from matplotlib.ticker import MaxNLocator
import cmocean

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import pplib                                                # noqa: E402
from pplib import read_vts, step_times                      # noqa: E402
from plot_keff import load, DAY, C_XX, C_YY, C_ISO, C_RELAX  # noqa: E402
from make_neck_movie import (SIGMA_SCALE, ice_alpha_cmap,   # noqa: E402
                             centered_cmap, sigma_ticks)

WANT = ("IcePhase", "VaporDensity", "Temperature")
SOL_DOF = {"IcePhase": 0, "Temperature": 1, "VaporDensity": 2}
LETTERS = "abcdefgh"
INK, MUTED = "#1a1a1a", "#5c5c5c"
FS, FS_SMALL = 9, 8                 # page-width figure: 8-9 pt at print size

DEFAULTS = {
    "title_time": "Effective thermal conductivity during dry-snow metamorphism",
    "title_ssa": "Effective thermal conductivity vs specific surface area",
    "xlabel_time": "Time [d]",
    "xlabel_ssa": r"Normalized specific surface area  SSA$\,/\,$SSA$_0$",
    "ylabel_norm": r"Normalized conductivity  $k_\mathrm{eff}\,/\,k_{\mathrm{eff},0}$",
    "ylabel_abs": r"$k_\mathrm{eff}$  [W m$^{-1}$ K$^{-1}$]",
}


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
def snap_step(fn) -> int:
    return int(re.search(r"sol[V]?_(\d+)\.(?:vts|dat)$", str(fn)).group(1))


def make_reader(run: Path, source: str):
    """(files, reader) for the snapshots; see make_packing_movie --source."""
    if source == "sol":
        from igakit.io import PetIGA
        nrb = PetIGA().read(str(run / "igasol.dat"))
        P = nrb.points
        X, Y = P[..., 0].T, P[..., 1].T

        def read(fn, want=WANT):
            sol = PetIGA().read_vec(str(fn), nrb)
            return {k: sol[..., SOL_DOF[k]].T for k in want}, X, Y
        files = glob.glob(str(run / "sol_*.dat"))
    else:
        read = read_vts
        files = glob.glob(str(run / "vtkOut" / "solV_*.vts"))
    return sorted(files, key=snap_step), read


def pick_snapshots(files, tmap, d, times_d, steps, n, t_ref, t_min):
    """Indices into `files` for the panels, in time order, duplicates removed.

    Only snapshots whose step has a k_eff sample at t >= t_min are eligible,
    so every marker lands on a measured value, never a neighbour's.
    """
    fsteps = np.array([snap_step(f) for f in files])
    ftimes = np.array([tmap.get(int(s), np.nan) for s in fsteps])
    measured = {int(s) for s, t in zip(d["step"], d["t"]) if t >= t_min}
    cand = np.array([i for i, s in enumerate(fsteps) if int(s) in measured])
    if cand.size == 0:
        return [], fsteps, ftimes
    if steps:
        want = [int(cand[np.argmin(np.abs(fsteps[cand] - s))]) for s in steps]
    else:
        if not times_d:
            times_d = np.linspace(t_ref / DAY, d["t"][-1] / DAY, n)
        want = [int(cand[np.argmin(np.abs(ftimes[cand] / DAY - td))]) for td in times_d]
    return sorted(set(want), key=lambda i: fsteps[i]), fsteps, ftimes


# ---------------------------------------------------------------------------
# Drawing
# ---------------------------------------------------------------------------
def _field(ax, fl, X, Y, norm, vapcm, icecm):
    sig = SIGMA_SCALE * pplib.supersaturation(fl["VaporDensity"], fl["Temperature"])
    XX, YY = X * 1e3, Y * 1e3
    h = 0.5 * (XX[0, 1] - XX[0, 0])
    kw = dict(origin="lower", extent=(XX.min() - h, XX.max() + h,
                                      YY.min() - h, YY.max() + h),
              interpolation="antialiased", interpolation_stage="rgba")
    ax.imshow(sig, cmap=vapcm, norm=norm, **kw)
    ax.imshow(fl["IcePhase"], cmap=icecm, vmin=0.0, vmax=1.0, **kw)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_linewidth(0.6); sp.set_color(MUTED)
    return XX, YY


def _scalebar(ax, XX, YY):
    """A bar ~1/4 of the domain, rounded to 1/2/5 x 10^k mm, bottom left."""
    Lx = XX.max() - XX.min()
    raw = Lx / 4
    e = 10 ** np.floor(np.log10(raw))
    L = max(m * e for m in (1, 2, 5) if m * e <= raw)
    x0 = XX.min() + 0.06 * Lx
    y0 = YY.min() + 0.07 * (YY.max() - YY.min())
    for lw, c in ((3.4, "white"), (1.8, INK)):
        ax.plot([x0, x0 + L], [y0, y0], color=c, lw=lw, solid_capstyle="butt",
                zorder=5)
    txt = f"{L * 1e3:g} µm" if L < 1 else f"{L:g} mm"
    t = ax.text(x0 + L / 2, y0 + 0.025 * Lx, txt, ha="center", va="bottom",
                fontsize=FS_SMALL - 1, color=INK, zorder=5)
    t.set_bbox(dict(facecolor="white", alpha=0.75, lw=0, pad=0.8))


def _curve(ax, x, ys, shown, relax, keys, points):
    """k series. `shown` masks the samples drawn in colour; `relax` (a subset
    of the rest) is drawn grey and unlabelled when --show-relaxation is on.
    points=True draws each measured sample as a dot with no joining line."""
    lw = {"kxx": 1.0, "kyy": 1.0, "kiso": 1.8}
    ms = {"kxx": 1.6, "kyy": 1.6, "kiso": 2.4}
    col = {"kxx": C_XX, "kyy": C_YY, "kiso": C_ISO}
    lab = {"kxx": r"$k_{xx}$", "kyy": r"$k_{yy}$", "kiso": r"$k_\mathrm{iso}$"}
    for key in keys:
        y = ys[key]
        style = dict(ls="none", marker="o", ms=ms[key], mew=0) if points \
            else dict(ls="-", lw=lw[key])
        if relax.any():
            ax.plot(x[relax], y[relax], color=C_RELAX, zorder=1, **style)
        ax.plot(x[shown], y[shown], color=col[key], zorder=2, label=lab[key],
                **style)
    ax.tick_params(labelsize=FS_SMALL, width=0.6, length=3)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_linewidth(0.6)


def _colorbars(fig, cax_ice, cax_sig, norm, vapcm):
    """Two slim bars, each titled with its symbol so no rotated label is needed."""
    ice_lut = ListedColormap(cmocean.cm.ice(np.linspace(0.5, 1.0, 256)))
    cb = fig.colorbar(ScalarMappable(cmap=ice_lut, norm=plt.Normalize(0.5, 1.0)),
                      cax=cax_ice, ticks=[0.5, 0.75, 1.0], format="%g")
    cax_ice.set_title(r"$\phi_i$", fontsize=FS, pad=4)
    cb.ax.tick_params(labelsize=FS_SMALL - 1, width=0.5, length=2)
    cb.outline.set_linewidth(0.5)

    cb = fig.colorbar(ScalarMappable(cmap=vapcm, norm=norm), cax=cax_sig,
                      ticks=sigma_ticks(norm, min_gap=0.09), extend="both",
                      extendfrac=0.04)
    # Left-aligned: centred, it would run back over the ice bar's ticks.
    cax_sig.set_title(r"$\sigma$ [$\times10^{-4}$]", fontsize=FS_SMALL, pad=4,
                      loc="left")
    cb.ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _p: f"{v:.2g}"))
    cb.ax.tick_params(labelsize=FS_SMALL - 1, width=0.5, length=2)
    cb.outline.set_linewidth(0.5)


def build(kind, snaps, x, ys, shown, relax, keys, norm, vapcm, icecm, a):
    """One figure. `snaps` is [(fields, X, Y, t, row)] in time order.

    Laid out in INCHES: the insets must be square, sit in a band along the
    bottom of the curve's own axes, and line up with the colour bars beside
    it. The y limits are then chosen so the data fills the band ABOVE the
    insets, and the left spine and its ticks are cut to the data range, so
    the empty space the insets live in carries no misleading y scale.
    """
    n = len(snaps)
    W = a.width
    ml, mr = 0.62, 1.00            # y label | colour-bar strip
    axw = W - ml - mr
    padx, pady = 0.06, 0.06        # insets <-> axes frame
    gap = 0.16                     # between insets
    s_in = (axw - 2 * padx - (n - 1) * gap) / n
    t_band = 0.20                  # inset titles
    lead = 0.34                    # inset titles -> data band, for the leaders
    band = 1.55                    # data band height
    head = 0.14                    # headroom for the marker letters
    axh = pady + s_in + t_band + lead + band + head
    top = 0.34 if not a.no_title else 0.06
    bot = 0.46
    H = top + axh + bot
    fig = plt.figure(figsize=(W, H))
    F = lambda x0, y0, w, h: (x0 / W, y0 / H, w / W, h / H)

    ax = fig.add_axes(F(ml, bot, axw, axh))
    ax.patch.set_alpha(0.0)
    _curve(ax, x, ys, shown, relax, keys, points=(kind == "ssa"))

    # y: the drawn data fills [f0, f1] of the axes height.
    drawn = shown | relax
    yv = np.concatenate([ys[k][drawn] for k in keys])
    ymin, ymax = float(yv.min()), float(yv.max())
    f0 = (pady + s_in + t_band + lead) / axh
    f1 = 1.0 - head / axh
    span = (ymax - ymin) / (f1 - f0) if ymax > ymin else 1.0
    ax.set_ylim(ymin - f0 * span, ymin - f0 * span + span)
    ticks = [t for t in MaxNLocator(4, steps=[1, 2, 5, 10]).tick_values(ymin, ymax)
             if ymin - 1e-9 <= t <= ymax + 1e-9]
    ax.set_yticks(ticks)
    ax.spines["left"].set_bounds(min(ticks[0], ymin), max(ticks[-1], ymax))
    xv = x[drawn]
    xpad = 0.03 * (xv.max() - xv.min())
    if kind == "time":
        ax.set_xlim(0.0, xv.max() + xpad)
    else:
        ax.set_xlim(xv.min() - xpad, xv.max() + xpad)

    normalized = not (a.absolute and kind == "time")
    if normalized:
        ax.plot(ax.get_xlim(), [1.0, 1.0], color="#999999", lw=0.6, ls=":",
                zorder=0)
    ax.legend(fontsize=FS_SMALL, frameon=False, handlelength=1.6,
              markerscale=2.5 if kind == "ssa" else 1.0,
              loc="upper left" if kind == "time" else "upper right")
    if kind == "ssa":
        # SSA falls as the packing sinters, so time runs right to left.
        ax.annotate("", xy=(0.40, 0.985), xytext=(0.56, 0.985),
                    xycoords="axes fraction",
                    arrowprops=dict(arrowstyle="->", color=MUTED, lw=0.9))
        ax.text(0.57, 0.985, "time", transform=ax.transAxes, va="center",
                fontsize=FS_SMALL, color=MUTED)

    xlabel = (a.xlabel_time or DEFAULTS["xlabel_time"]) if kind == "time" \
        else (a.xlabel_ssa or DEFAULTS["xlabel_ssa"])
    ylabel = a.ylabel or (DEFAULTS["ylabel_norm"] if normalized
                          else DEFAULTS["ylabel_abs"])
    ax.set_xlabel(xlabel, fontsize=FS)
    ax.set_ylabel(ylabel, fontsize=FS)
    # Centre the y label on the data band, not on the whole axes.
    ax.yaxis.set_label_coords(-0.075, 0.5 * (f0 + f1))

    # Insets in the order their markers fall along x, so no leader crosses
    # another: time order on the time figure, reversed on the SSA figure.
    ymark = ys["kiso"]
    order = sorted(range(n), key=lambda i: x[snaps[i][4]])
    for slot, i in enumerate(order):
        fl, X, Y, t, row = snaps[i]
        axi = fig.add_axes(F(ml + padx + slot * (s_in + gap), bot + pady,
                             s_in, s_in))
        XX, YY = _field(axi, fl, X, Y, norm, vapcm, icecm)
        if slot == 0:
            _scalebar(axi, XX, YY)
        L = LETTERS[i]
        axi.set_title(f"({L})  t = {t / DAY:.1f} d", fontsize=FS, color=INK, pad=3)
        px, py = x[row], ymark[row]
        # The letter sits INSIDE the marker: beside it, it lands on k_xx or
        # k_yy, which run within a few percent of k_iso.
        ax.plot([px], [py], "o", ms=11, mfc="white", mec=INK, mew=1.0, zorder=6)
        ax.text(px, py, L, ha="center", va="center", fontsize=FS_SMALL - 1,
                color=INK, fontweight="bold", zorder=7)
        # From just above the inset's title to just below the marker's letter.
        con = ConnectionPatch(xyA=(0.5, 1.0), coordsA=axi.transAxes,
                              xyB=(px, py), coordsB=ax.transData,
                              color="#a0a0a0", lw=0.6, ls=(0, (3, 2)),
                              zorder=3, shrinkA=15, shrinkB=7)
        fig.add_artist(con)

    xc = ml + axw + 0.14
    cax_ice = fig.add_axes(F(xc, bot + pady, 0.08, s_in))
    cax_sig = fig.add_axes(F(xc + 0.08 + 0.36, bot + pady + 0.05, 0.08, s_in - 0.10))
    _colorbars(fig, cax_ice, cax_sig, norm, vapcm)

    if not a.no_title:
        title = (a.title_time or DEFAULTS["title_time"]) if kind == "time" \
            else (a.title_ssa or DEFAULTS["title_ssa"])
        fig.suptitle(title, x=(ml + 0.5 * axw) / W, y=1.0 - 0.06 / H,
                     ha="center", va="top", fontsize=FS + 1.5, color=INK)
    return fig


# ---------------------------------------------------------------------------
def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dir", default=".", help="run directory (default: cwd)")
    p.add_argument("--save-dir", default=None,
                   help="output directory (default: <dir>/plots/keff/snapshots)")
    p.add_argument("--times", type=float, nargs="+", default=None,
                   help="snapshot times [d]; the nearest snapshot to each is used")
    p.add_argument("--steps", type=int, nargs="+", default=None,
                   help="snapshot steps, overriding --times")
    p.add_argument("--n-snapshots", type=int, default=4,
                   help="default count, evenly spaced baseline..end (default 4)")
    p.add_argument("--baseline-days", type=float, default=1.0,
                   help="reference time [d] for the subscript-0 values; the first "
                        "sample at or after it (default 1, as in plot_keff.py)")
    p.add_argument("--source", choices=("vts", "sol"), default="vts",
                   help="vts: vtkOut/solV_*.vts (plenty at page width); "
                        "sol: full-resolution sol_*.dat via igakit")
    p.add_argument("--width", type=float, default=7.2,
                   help="figure width [in]; 7.2 is a full two-column page")
    p.add_argument("--absolute", action="store_true",
                   help="time figure: k_eff in W/m/K instead of k/k_0")
    p.add_argument("--iso-only", action="store_true",
                   help="plot k_iso only, not k_xx and k_yy")
    p.add_argument("--show-relaxation", action="store_true",
                   help="also draw the samples before the reference time, in "
                        "grey and unlabelled (default: not drawn)")
    p.add_argument("--sat-clip", type=float, default=0.1,
                   help="percentile clipped off each end of sigma (default 0.1)")
    p.add_argument("--title-time", default=None)
    p.add_argument("--title-ssa", default=None)
    p.add_argument("--no-title", action="store_true")
    p.add_argument("--xlabel-time", default=None)
    p.add_argument("--xlabel-ssa", default=None)
    p.add_argument("--ylabel", default=None, help="y label for both figures")
    p.add_argument("--formats", nargs="+", default=["pdf", "png"])
    p.add_argument("--dpi", type=int, default=400, help="PNG and raster dpi")
    a = p.parse_args(argv)

    run = Path(a.dir).resolve()
    d = load(run)
    if d is None:
        print(f"  no k_eff CSV (or no SSA_evo.dat / -Lx -Ly) in {run}; nothing to plot")
        return 0

    files, reader = make_reader(run, a.source)
    if not files:
        print(f"  no {'sol_*.dat' if a.source == 'sol' else 'vtkOut/solV_*.vts'} "
              f"in {run}; nothing to plot")
        return 0
    tmap = step_times(str(run))
    for s, t in zip(d["step"], d["t"]):          # k_eff CSV covers every step
        tmap.setdefault(int(s), float(t))

    early = d["t"] < a.baseline_days * DAY
    if early.all():
        print(f"  run ends before the {a.baseline_days:g} d reference time; "
              "nothing to plot", file=sys.stderr)
        return 1
    ib = int(np.argmax(~early))
    shown = ~early
    relax = early if a.show_relaxation else np.zeros_like(early)
    idx, fsteps, ftimes = pick_snapshots(
        files, tmap, d, a.times, a.steps, a.n_snapshots, d["t"][ib],
        t_min=0.0 if a.show_relaxation else d["t"][ib])
    if not idx:
        print("  no snapshot shares a step with a k_eff sample", file=sys.stderr)
        return 1
    if len(idx) > len(LETTERS):
        print(f"  at most {len(LETTERS)} snapshots", file=sys.stderr)
        return 1

    snaps, pore = [], []
    for i in idx:
        fl, X, Y = reader(files[i], want=WANT)
        t = ftimes[i] if np.isfinite(ftimes[i]) else tmap.get(int(fsteps[i]), 0.0)
        row = int(np.flatnonzero(d["step"] == fsteps[i])[0])
        snaps.append([fl, X, Y, float(d["t"][row]), row])
        s = SIGMA_SCALE * pplib.supersaturation(fl["VaporDensity"], fl["Temperature"])
        pore.append(s[fl["IcePhase"] < 0.5])
        print(f"  snapshot step {fsteps[i]}: t = {t / DAY:.2f} d")
    pore = np.concatenate(pore)
    lo, hi = np.percentile(pore, [a.sat_clip, 100.0 - a.sat_clip])
    norm = AsinhNorm(linear_width=max(max(abs(lo), abs(hi)) / 300.0, 1e-12),
                     vmin=lo, vmax=hi)
    vapcm = centered_cmap(cmocean.cm.balance, norm)
    icecm = ice_alpha_cmap()

    keys = ("kiso",) if a.iso_only else ("kxx", "kyy", "kiso")
    kn = {k: d[k] / d[k][ib] for k in ("kxx", "kyy", "kiso")}
    out = Path(a.save_dir) if a.save_dir else run / "plots" / "keff" / "snapshots"
    os.makedirs(out, exist_ok=True)

    plt.rcParams.update({"font.family": "sans-serif", "mathtext.fontset": "dejavusans",
                         "pdf.fonttype": 42, "svg.fonttype": "none"})
    figs = {
        "keff_time_snapshots": build(
            "time", snaps, d["t"] / DAY, d if a.absolute else kn, shown, relax,
            keys, norm, vapcm, icecm, a),
        "keff_ssa_snapshots": build(
            "ssa", snaps, d["ssa"] / d["ssa"][ib], kn, shown, relax,
            keys, norm, vapcm, icecm, a),
    }
    for stem, fig in figs.items():
        for fmt in a.formats:
            path = out / f"{stem}.{fmt}"
            fig.savefig(path, dpi=a.dpi, bbox_inches="tight", pad_inches=0.03)
            print(f"  wrote {path}")
        plt.close(fig)
    print(f"  reference (subscript 0): t_0 = {d['t'][ib] / DAY:.3f} d, step "
          f"{d['step'][ib]}: k_xx,0 = {d['kxx'][ib]:.4g}, k_yy,0 = {d['kyy'][ib]:.4g}, "
          f"k_iso,0 = {d['kiso'][ib]:.4g} W/m/K, SSA_0 = {d['ssa'][ib]:.4g} 1/m")
    print(f"  sigma x{SIGMA_SCALE:g} range {lo:+.3g} .. {hi:+.3g} over the shown pores")
    return 0


if __name__ == "__main__":
    sys.exit(main())
