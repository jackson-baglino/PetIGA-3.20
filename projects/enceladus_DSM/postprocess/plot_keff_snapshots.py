#!/usr/bin/env python3
"""plot_keff_snapshots.py — k_eff curve with microstructure snapshots inset.

    python3 plot_keff_snapshots.py --dir <run> [--times 1 10 20 30] [--steps ...]
        [--n-snapshots 4] [--width-mm 130] [--absolute] [--iso-only]
        [--title-time STR] [--title-ssa STR] [--no-title]
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

THE OPENING FRAME IS t = 0. The curve and the snapshots start at the same
frame the movies open on (pplib.opening_step): the first snapshot with
1 s <= t <= 1 h. From 2026-09-26 the solver writes one at t = 1 s
(-t_out_first), by which time the vapour field has relaxed onto the ice
geometry; the IC itself (step 0) and step 1 (t = dt) still carry a uniform
vapour field. Older runs fall back to their first log-spaced snapshot
(~71 s on the 2026-09-25 batch), then to step 1, then to the IC. Samples
before the opening frame are not drawn. Everything after it is, including
the fast first-day rise as the initial packing relaxes -- that is data, and
it is explained in the manuscript rather than annotated on the figure.

NORMALIZATION. Each quantity is divided by its own value at the opening
frame -- k_xx by k_xx,0, k_iso by k_iso,0, SSA by SSA_0 -- which is what
the subscript 0 means; the caption should say so. The reference time and
values are printed for it. --absolute plots k_eff in W m^-1 K^-1 on the time
figure instead; the SSA figure is always normalized on both axes.

ONLY MEASURED SAMPLES. k_eff and SSA are paired by step (plot_keff.load
drops a k_eff sample with no SSA row rather than borrowing a neighbour's),
and every curve runs from the first measured sample to the last -- nothing
is extrapolated past the range the simulation covered. Snapshots are only taken at steps that have a k_eff sample, so every marker
sits on a measured value.

Y AXIS. The axis runs the full height of the axes, down behind the insets,
so they read as drawn ON the plot. For k it starts at 0 whenever that
leaves the data at least DATA_MIN_IN of height; a curve too flat for that
gets a y range chosen to give it that height instead.

SNAPSHOTS. --times (days) picks the eligible snapshot nearest each time;
--steps names them exactly. By default --n-snapshots are spaced evenly in
time from the opening frame to the end of the run, so (a) is the opening
frame.

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
from pplib import read_vts, step_times, opening_step                      # noqa: E402
from plot_keff import load, DAY, C_XX, C_YY, C_ISO           # noqa: E402
from make_neck_movie import (SIGMA_SCALE, ice_alpha_cmap,   # noqa: E402
                             centered_cmap, sigma_ticks)

WANT = ("IcePhase", "VaporDensity", "Temperature")
SOL_DOF = {"IcePhase": 0, "Temperature": 1, "VaporDensity": 2}
LETTERS = "abcdefgh"
INK, MUTED = "#1a1a1a", "#5c5c5c"
# Type sizes are PRINT sizes: the figure is built at its final width
# (--width-mm) and saved without cropping, so nothing is rescaled on the page.
FS_TITLE, FS, FS_SMALL, FS_TINY = 11, 10, 9, 8
MM = 1.0 / 25.4
DATA_MIN_IN = 1.3                   # least height [in] the curve may occupy

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
                fontsize=FS_TINY, color=INK, zorder=5)
    t.set_bbox(dict(facecolor="white", alpha=0.75, lw=0, pad=0.8))


def _curve(ax, x, ys, keys):
    """k series, from the first measured sample to the last."""
    lw = {"kxx": 1.0, "kyy": 1.0, "kiso": 1.8}
    col = {"kxx": C_XX, "kyy": C_YY, "kiso": C_ISO}
    lab = {"kxx": r"$k_{xx}$", "kyy": r"$k_{yy}$", "kiso": r"$k_\mathrm{iso}$"}
    for key in keys:
        ax.plot(x, ys[key], "-", lw=lw[key], color=col[key], zorder=2, label=lab[key])
    ax.tick_params(labelsize=FS_SMALL, width=0.6, length=3, pad=2)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_linewidth(0.6)


def _colorbars(fig, cax_ice, cax_sig, norm, vapcm):
    """Two horizontal bars in a strip above the axes, each labelled on its
    left. Horizontal because at 130 mm a vertical pair beside the axes would
    cost the insets a quarter of their width."""
    def _label(cax, text):
        cax.text(-0.06, 0.5, text, transform=cax.transAxes, ha="right",
                 va="center", fontsize=FS_SMALL, color=INK)

    ice_lut = ListedColormap(cmocean.cm.ice(np.linspace(0.5, 1.0, 256)))
    cb = fig.colorbar(ScalarMappable(cmap=ice_lut, norm=plt.Normalize(0.5, 1.0)),
                      cax=cax_ice, orientation="horizontal",
                      ticks=[0.5, 0.75, 1.0], format="%g")
    _label(cax_ice, r"$\phi_i$")
    cb.ax.tick_params(labelsize=FS_TINY, width=0.5, length=2, pad=1.5)
    cb.outline.set_linewidth(0.5)

    # Thinned harder than the vertical bar did: horizontally each label is
    # as wide as three or four characters, not as tall as one.
    cb = fig.colorbar(ScalarMappable(cmap=vapcm, norm=norm), cax=cax_sig,
                      orientation="horizontal", extend="both", extendfrac=0.04,
                      ticks=sigma_ticks(norm, min_gap=0.16))
    _label(cax_sig, r"$\sigma$ [$\times10^{-4}$]")
    cb.ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _p: f"{v:.2g}"))
    cb.ax.tick_params(labelsize=FS_TINY, width=0.5, length=2, pad=1.5)
    cb.outline.set_linewidth(0.5)


def build(kind, snaps, x, ys, keys, norm, vapcm, icecm, a):
    """One figure. `snaps` is [(fields, X, Y, t, row)] in time order.

    Laid out in INCHES at the final print width: the insets must be square
    and sit in a band along the bottom of the curve's own axes, with the
    colour bars in a strip above. The data sits in the band ABOVE the insets; the y axis runs on down
    behind them, so the snapshots read as drawn on the plot.

    The axes height follows from the y range: with y starting at 0, the
    data's share of the height is fixed by (ymax - ymin) / ymax, and the
    inset band has to fit in the rest.
    """
    n = len(snaps)
    W = a.width_mm * MM
    ml, mr = 0.50, 0.08            # y label + ticks | right edge
    axw = W - ml - mr
    padx, pady = 0.05, 0.05        # insets <-> axes frame
    gap = 0.08                     # between insets
    s_in = (axw - 2 * padx - (n - 1) * gap) / n
    t_band = 0.21                  # inset titles
    lead = 0.30                    # inset titles -> data band, for the leaders
    head = 0.14                    # headroom for the marker letters
    below = pady + s_in + t_band + lead          # inset band, in inches

    yv = np.concatenate([ys[k] for k in keys])
    ymin, ymax = float(yv.min()), float(yv.max())
    # y from 0 if that leaves the data DATA_MIN_IN of height; otherwise a
    # lower floor that does. Either way ymin sits exactly at the band's top.
    band = below * (ymax - ymin) / ymin if ymin > 0 else 0.0
    if band >= DATA_MIN_IN:
        ybot = 0.0
    else:
        band = DATA_MIN_IN
        ybot = ymin - (ymax - ymin) * below / band
    axh = below + band + head
    ytop = ymax + (ymax - ybot) * head / (below + band)
    cb_h, cb_lab, cb_gap = 0.07, 0.17, 0.10     # bar | its tick labels | to axes
    strip = cb_h + cb_lab + cb_gap
    top = 0.30 if not a.no_title else 0.05
    bot = 0.40
    H = top + strip + axh + bot
    fig = plt.figure(figsize=(W, H))
    F = lambda x0, y0, w, h: (x0 / W, y0 / H, w / W, h / H)

    ax = fig.add_axes(F(ml, bot, axw, axh))
    ax.patch.set_alpha(0.0)
    _curve(ax, x, ys, keys)

    ax.set_ylim(ybot, ytop)
    ax.yaxis.set_major_locator(MaxNLocator(6, steps=[1, 2, 2.5, 5, 10]))
    xv = x
    xpad = 0.03 * (xv.max() - xv.min())
    if kind == "time":
        # A little room left of t = 0 so marker (a) is not cut by the spine.
        ax.set_xlim(-xpad, xv.max() + xpad)
    else:
        ax.set_xlim(xv.min() - xpad, xv.max() + xpad)

    normalized = not (a.absolute and kind == "time")
    if normalized:
        ax.plot(ax.get_xlim(), [1.0, 1.0], color="#999999", lw=0.6, ls=":",
                zorder=0)
    # One row, where neither curve nor leader runs: under the curve at the
    # left of the time figure (leader (a) stays below the band's floor), and
    # top right on the SSA figure, where the curve has already fallen away.
    f0 = below / axh
    if kind == "time":
        where = dict(loc="lower left", bbox_to_anchor=(0.05, f0 + 0.01))
    else:
        where = dict(loc="upper right", bbox_to_anchor=(1.0, 0.93))
    ax.legend(fontsize=FS_SMALL, frameon=False, handlelength=1.4, ncol=len(keys),
              columnspacing=1.0, handletextpad=0.5,
              **where)
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
    ax.set_xlabel(xlabel, fontsize=FS, labelpad=2)
    ax.set_ylabel(ylabel, fontsize=FS, labelpad=3)

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
        axi.set_title(f"({L}) {t / DAY:.1f} d", fontsize=FS_SMALL, color=INK, pad=2)
        px, py = x[row], ymark[row]
        # The letter sits INSIDE the marker: beside it, it lands on k_xx or
        # k_yy, which run within a few percent of k_iso.
        ax.plot([px], [py], "o", ms=10, mfc="white", mec=INK, mew=0.9, zorder=6)
        ax.text(px, py, L, ha="center", va="center", fontsize=FS_TINY,
                color=INK, fontweight="bold", zorder=7)
        # From just above the inset's title to just below the marker's letter.
        con = ConnectionPatch(xyA=(0.5, 1.0), coordsA=axi.transAxes,
                              xyB=(px, py), coordsB=ax.transData,
                              color="#a0a0a0", lw=0.6, ls=(0, (3, 2)),
                              zorder=3, shrinkA=13, shrinkB=6)
        fig.add_artist(con)

    # Colour bars: phi_i then sigma, left to right across the axes' width.
    y_cb = bot + axh + cb_gap + cb_lab
    lab_ice, lab_sig, sep = 0.22, 0.82, 0.30      # label widths, gap between
    tail = 0.14                   # sigma's end arrow and last tick label
    w_ice = 0.30 * (axw - lab_ice - lab_sig - sep - tail)
    w_sig = axw - lab_ice - lab_sig - sep - tail - w_ice
    cax_ice = fig.add_axes(F(ml + lab_ice, y_cb, w_ice, cb_h))
    cax_sig = fig.add_axes(F(ml + lab_ice + w_ice + sep + lab_sig, y_cb, w_sig, cb_h))
    _colorbars(fig, cax_ice, cax_sig, norm, vapcm)

    if not a.no_title:
        title = (a.title_time or DEFAULTS["title_time"]) if kind == "time" \
            else (a.title_ssa or DEFAULTS["title_ssa"])
        fig.suptitle(title, x=0.5, y=1.0 - 0.05 / H,
                     ha="center", va="top", fontsize=FS_TITLE, color=INK)
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
                   help="default count, evenly spaced opening frame..end (default 4)")
    p.add_argument("--source", choices=("vts", "sol"), default="vts",
                   help="vts: vtkOut/solV_*.vts (plenty at page width); "
                        "sol: full-resolution sol_*.dat via igakit")
    p.add_argument("--width-mm", type=float, default=130.0,
                   help="final printed width [mm] (default 130). The figure is "
                        "built at this size, so the fonts are print sizes")
    p.add_argument("--absolute", action="store_true",
                   help="time figure: k_eff in W/m/K instead of k/k_0")
    p.add_argument("--iso-only", action="store_true",
                   help="plot k_iso only, not k_xx and k_yy")
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

    # The opening frame: the movies' rule, over snapshots that have a k_eff
    # sample, falling back to step 1 and then the IC as pplib.drop_ic does.
    measured = set(int(v) for v in d["step"])
    snap_steps = [snap_step(f) for f in files if snap_step(f) in measured]
    if not snap_steps:
        print("  no snapshot shares a step with a k_eff sample", file=sys.stderr)
        return 1
    op = opening_step(snap_steps, [tmap.get(s_, np.inf) for s_ in snap_steps])
    if op is None:
        op = 1 if 1 in snap_steps else snap_steps[0]
        print(f"  note: no snapshot between 1 s and 1 h; opening on step {op}")
    ib = int(np.flatnonzero(d["step"] == op)[0])
    keep = d["t"] >= d["t"][ib]
    d = {k: (v[keep] if isinstance(v, np.ndarray) else v) for k, v in d.items()}
    ib = 0
    idx, fsteps, ftimes = pick_snapshots(
        files, tmap, d, a.times, a.steps, a.n_snapshots, d["t"][ib],
        t_min=d["t"][ib])
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
            "time", snaps, d["t"] / DAY, d if a.absolute else kn,
            keys, norm, vapcm, icecm, a),
        "keff_ssa_snapshots": build(
            "ssa", snaps, d["ssa"] / d["ssa"][ib], kn,
            keys, norm, vapcm, icecm, a),
    }
    for stem, fig in figs.items():
        for fmt in a.formats:
            path = out / f"{stem}.{fmt}"
            # No bbox_inches="tight": it would crop to the ink and the
            # saved width would no longer be --width-mm.
            fig.savefig(path, dpi=a.dpi)
            print(f"  wrote {path}")
        plt.close(fig)
    print(f"  reference (subscript 0) = opening frame: t_0 = {d['t'][ib]:.4g} s, step "
          f"{d['step'][ib]}: k_xx,0 = {d['kxx'][ib]:.4g}, k_yy,0 = {d['kyy'][ib]:.4g}, "
          f"k_iso,0 = {d['kiso'][ib]:.4g} W/m/K, SSA_0 = {d['ssa'][ib]:.4g} 1/m")
    print(f"  sigma x{SIGMA_SCALE:g} range {lo:+.3g} .. {hi:+.3g} over the shown pores")
    return 0


if __name__ == "__main__":
    sys.exit(main())
