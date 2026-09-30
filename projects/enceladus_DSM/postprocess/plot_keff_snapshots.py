#!/usr/bin/env python3
"""plot_keff_snapshots.py — k_eff curve with microstructure snapshots inset.

    python3 plot_keff_snapshots.py --dir <run> [--times 1 10 20 30] [--steps ...]
        [--n-snapshots 4] [--width-mm 170] [--absolute] [--iso-only]
        [--title] [--title-time STR] [--title-ssa STR] [--title-combined STR]
        [--xlabel-time STR] [--xlabel-ssa STR] [--ylabel STR]
        [--source vts|sol] [--save-dir DIR] [--formats pdf png]

Two page-width manuscript figures, each one k_eff curve with 3-4
microstructure snapshots inset along the bottom of the same axes, under the
curve. Each snapshot is numbered, and the same circled number marks its
instant on the curve, joined to it by a leader line:

    keff_time_snapshots.{pdf,png}   k_eff / k_eff,0  vs  time [d]
    keff_ssa_snapshots.{pdf,png}    k_eff / k_eff,0  vs  SSA / SSA_0
    keff_combined_snapshots.{pdf,png}
                                    three panels, (a) the snapshots,
                                    (b) k vs time, (c) k vs SSA; the
                                    curves' points numbered to match

INSTANTS ARE NUMBERED, PANELS ARE LETTERED. AGU asks for multi-part
figures to carry sequential lowercase panel labels, so letters belong to
the panels; the snapshot instants take circled numbers 1-4 so the two
never collide in a caption ("(b) ... at instants 1-4").

SIZE. AGU journals take figures 50-170 mm wide (105-170 mm for two
columns) and at most 228 mm tall; --width-mm defaults to the full 170 mm.
The fonts are print sizes (>= 8 pt, AGU's minimum) at that width.

The insets are ordered by where their markers fall on the x axis, so the
leaders never cross: left to right in time on the time figure, and RIGHT to
left in time on the SSA figure, where SSA falls as the packing sinters.

PANELS. They follow make_packing_movie.py, and they import its helpers so
the figure and the movie cannot drift apart: supersaturation
sigma = rho_v/rho_vs(T) - 1 on cmocean `balance`, re-centred so the pale
middle is sigma = 0, on an asinh scale; ice painted on top with cmocean
`ice`, transparent below phi = 0.5.

SIGMA RANGE. Symmetric about zero, so sigma = 0 sits at the middle of the
bar: [-v, +v] with v = min(|min sigma|, |max sigma|) over the pore space of
the SHOWN snapshots. v is the SMALLER extreme, so that side of the bar ends
exactly at the data (square cap) and the other side saturates (triangular
cap = values beyond the bar, drawn in its end colour). All panels share it.
This departs from the movies' asymmetric range on purpose: in a still figure
the reader must be able to read sigma = 0 off the middle of the bar.

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
is extrapolated past the range the simulation covered. Every marker sits
on a measured k_eff sample: the snapshot's own step when it has one, else
the nearest sample within --pair-tol (default 2 %) of its time -- k_eff is
sampled every N steps and snapshots are log-spaced, so after the first days
the two rarely share a step. Each pairing's offset is printed; a snapshot
title gives the image's own time.

Y AXIS. The axis runs the full height of the axes, down behind the insets,
so they read as drawn ON the plot. For k it starts at 0 whenever that
leaves the data at least DATA_MIN_IN of height; a curve too flat for that
gets a y range chosen to give it that height instead.

SNAPSHOTS. --times (days) picks the eligible snapshot nearest each time;
--steps names them exactly. By default --n-snapshots are spaced evenly in
time from the opening frame to the end of the run, so (a) is the opening
frame.

CLIP MAP (--clip-map, off by default). A diagnostic, not a manuscript
figure: the same snapshots, larger, in a 2 x 2 grid with the same colour
bars, plus an outline around the pore space where sigma is beyond the bar
(below -v in yellow; above +v in green, if that side clips too) and so is
drawn in the bar's end colour. Each panel's title gives that share of its
pore area. It outlines the region, not the sigma = -v level: small isolated
pores are often clipped whole, and their edge is the ice.
Written as sigma_out_of_range.{pdf,png} next to the other two.

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
from matplotlib.patches import ConnectionPatch, FancyArrowPatch
from matplotlib.ticker import MaxNLocator
import matplotlib.transforms as mtransforms
import cmocean

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import pplib                                                # noqa: E402
from pplib import read_vts, step_times, opening_step                      # noqa: E402
from plot_keff import (load, DAY, C_XX, C_YY, C_ISO,           # noqa: E402
                       LS_XX, LS_YY, LS_ISO)
from make_neck_movie import (SIGMA_SCALE, ice_alpha_cmap,   # noqa: E402
                             centered_cmap, sigma_ticks)

WANT = ("IcePhase", "VaporDensity", "Temperature")
SOL_DOF = {"IcePhase": 0, "Temperature": 1, "VaporDensity": 2}
LETTERS = "12345678"                # instant labels: circled numbers
PANELS = "abc"                      # panel labels of the combined figure
MAX_H_MM = 228.0                    # AGU's maximum figure height
INK, MUTED = "#1a1a1a", "#5c5c5c"
# Type sizes are PRINT sizes: the figure is built at its final width
# (--width-mm) and saved without cropping, so nothing is rescaled on the page.
FS_TITLE, FS, FS_SMALL, FS_TINY = 11, 10, 9, 8
MM = 1.0 / 25.4
DATA_MIN_IN = 1.1                   # least height [in] the curve may occupy

DEFAULTS = {
    "title_time": "Effective thermal conductivity during dry-snow metamorphism",
    "title_ssa": "Effective thermal conductivity vs specific surface area",
    "xlabel_time": "Time [d]",
    # Symbols only: these are manuscript figures, and the caption defines
    # the quantities. --xlabel-* / --ylabel give longer ones for slides.
    "xlabel_ssa": r"SSA$\,/\,$SSA$_0$",
    "ylabel_norm": r"$k_\mathrm{eff}\,/\,k_{\mathrm{eff},0}$",
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


def pair_rows(fsteps, ftimes, d, tol, t_min):
    """{file index: k_eff row} for every snapshot that has a k_eff sample on
    its own step, or -- with tol > 0 -- the nearest sample in time within
    tol * t of it. The marker goes on that measured sample, never on an
    interpolated value; only the image may be a step or two away from it.

    Why a tolerance exists: a run that samples k_eff every N steps (the
    campaign's kf5) and writes snapshots on a log-spaced schedule shares
    almost no steps between the two after the first days -- on the 3a
    shakedown, none after 3.4 d, where the nearest samples sit 2 steps away,
    0.5-1.9 % of t."""
    kt, ks = d["t"], d["step"].astype(int)
    out = {}
    for i, (s, t) in enumerate(zip(fsteps, ftimes)):
        exact = np.flatnonzero(ks == int(s))
        if exact.size:
            j = int(exact[0])
        elif tol > 0 and np.isfinite(t):
            j = int(np.argmin(np.abs(kt - t)))
            if abs(kt[j] - t) > tol * max(t, 0.0):
                continue
        else:
            continue
        if kt[j] >= t_min:
            out[i] = j
    return out


def pick_snapshots(files, tmap, d, times_d, steps, n, t_ref, t_min, tol=0.0):
    """Indices into `files` for the panels, in time order, duplicates removed,
    and the k_eff row each one's marker sits on (see pair_rows).
    """
    fsteps = np.array([snap_step(f) for f in files])
    ftimes = np.array([tmap.get(int(s), np.nan) for s in fsteps])
    rows = pair_rows(fsteps, ftimes, d, tol, t_min)
    cand = np.array(sorted(rows))
    if cand.size == 0:
        return [], fsteps, ftimes, rows
    if steps:
        want = [int(cand[np.argmin(np.abs(fsteps[cand] - s))]) for s in steps]
    else:
        if not times_d:
            times_d = np.linspace(t_ref / DAY, d["t"][-1] / DAY, n)
        want = [int(cand[np.argmin(np.abs(ftimes[cand] / DAY - td))]) for td in times_d]
    return sorted(set(want), key=lambda i: fsteps[i]), fsteps, ftimes, rows


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
    txt = rf"{L * 1e3:g} $\mu$m" if L < 1 else f"{L:g} mm"
    t = ax.text(x0 + L / 2, y0 + 0.025 * Lx, txt, ha="center", va="bottom",
                fontsize=FS_TINY, color=INK, zorder=5)
    t.set_bbox(dict(facecolor="white", alpha=0.75, lw=0, pad=0.8))


def _curve(ax, x, ys, keys):
    """k series, from the first measured sample to the last."""
    # Dashed / dotted / solid, as plot_keff.py: normalized, the three often
    # coincide, and the dash pattern keeps each visible.
    lw = {"kxx": 1.1, "kyy": 1.3, "kiso": 1.8}
    col = {"kxx": C_XX, "kyy": C_YY, "kiso": C_ISO}
    ls = {"kxx": LS_XX, "kyy": LS_YY, "kiso": LS_ISO}
    lab = {"kxx": r"$k_{xx}$", "kyy": r"$k_{yy}$", "kiso": r"$k_\mathrm{iso}$"}
    for key in keys:
        # k_iso UNDER the others: where they coincide, the dashes show on it.
        ax.plot(x, ys[key], ls=ls[key], lw=lw[key], color=col[key],
                zorder={"kiso": 2, "kyy": 3, "kxx": 4}[key],
                label=lab[key], dash_capstyle="round")
    ax.tick_params(labelsize=FS_SMALL, width=0.6, length=3, pad=2)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_linewidth(0.6)


MARK_MS = 10.0                  # circled-letter diameter [pt]


def _circled(ax, x, y, letter, transform=None, clip_on=True):
    """The circled letter: white disc, dark ring, bold letter. The ONE way a
    snapshot's letter is drawn, on the curves and over the snapshots alike."""
    kw = {} if transform is None else {"transform": transform}
    ax.plot([x], [y], "o", ms=MARK_MS, mfc="white", mec=INK, mew=0.9,
            zorder=6, clip_on=clip_on, **kw)
    ax.text(x, y, pplib.bold(letter), ha="center", va="center", fontsize=FS_TINY,
            color=INK, zorder=7, clip_on=clip_on, **kw)


def _mark(ax, x, y, letter):
    """A snapshot's letter on the curve."""
    _circled(ax, x, y, letter)


def _snap_title(ax, letter, label, gap_pt=3.0, lift_pt=2.0):
    """A snapshot's title: its circled letter, then `label`, centred as a
    group over the panel -- the same circled letter as on the curves."""
    fig = ax.figure
    txt = ax.annotate(label, (0.5, 1.0), xycoords="axes fraction",
                      xytext=(0.0, 0.0), textcoords="offset points",
                      ha="left", va="center", fontsize=FS_SMALL, color=INK,
                      annotation_clip=False)
    w_pt = txt.get_window_extent(fig.canvas.get_renderer()).width * 72.0 / fig.dpi
    x0 = -0.5 * (MARK_MS + gap_pt + w_pt)
    y = lift_pt + 0.5 * MARK_MS
    txt.xyann = (x0 + MARK_MS + gap_pt, y)
    _circled(ax, 0.5, 1.0, letter, clip_on=False,
             transform=mtransforms.offset_copy(ax.transAxes, fig=fig,
                                              x=x0 + 0.5 * MARK_MS, y=y,
                                              units="points"))


def _time_arrow(ax, x, y, s0=0.26, s1=0.52, off_pt=13.0, rad=0.08):
    """A curved arrow beside the curve, pointing the way the data move as
    time increases.

    Against SSA the curve runs right to left in time, so a horizontal arrow
    in a corner said "time goes left" but not along what. This one is laid
    along the curve itself: its ends are the curve's points at fractions
    s0 and s1 of its drawn ARC LENGTH (from the t = 0 end, in display space,
    so the placement does not depend on the axis units), shifted off_pt
    points to the curve's upper side, and bowed away from it. Samples are
    in time order, so s0 -> s1 is the direction of increasing time.

    Call it after the axis limits are final: the positions are computed in
    display space and stored in data coordinates.
    """
    fig = ax.figure
    P = ax.transData.transform(np.column_stack([x, y]))
    seg = np.hypot(*np.diff(P, axis=0).T)
    L = np.concatenate([[0.0], np.cumsum(seg)])
    if L[-1] <= 0:
        return
    i0 = int(np.searchsorted(L, s0 * L[-1]))
    i1 = int(np.searchsorted(L, s1 * L[-1]))
    if i1 <= i0:
        return
    A, B = P[i0], P[i1]
    d = (B - A) / np.hypot(*(B - A))
    nrm = np.array([-d[1], d[0]])
    if nrm[1] < 0:                    # the upper side of the curve
        nrm = -nrm
    off = off_pt * fig.dpi / 72.0
    inv = ax.transData.inverted()
    A2, B2 = inv.transform(A + nrm * off), inv.transform(B + nrm * off)
    # rad > 0 bows the arc toward +nrm, away from the curve (checked on the
    # rendered figure: the opposite sign bends it back across the data).
    ax.add_patch(FancyArrowPatch(A2, B2, connectionstyle=f"arc3,rad={rad}",
                                 arrowstyle="-|>", mutation_scale=9, lw=0.9,
                                 color=MUTED, zorder=5, shrinkA=0, shrinkB=0))
    # "time" just outside the arc's apex (the apex sits rad * chord off the
    # chord's midpoint).
    mid = 0.5 * (A + B) + nrm * (off + abs(rad) * np.hypot(*(B - A))
                                 + 6.0 * fig.dpi / 72.0)
    ax.text(*inv.transform(mid), "time", ha="center", va="center",
            fontsize=FS_SMALL, color=MUTED, zorder=5)


def _colorbars(fig, cax_ice, cax_sig, norm, vapcm, sig_extend):
    """Two horizontal bars in a strip above the axes, each labelled on its
    left. Horizontal because at page width a vertical pair beside the axes would
    cost the insets a quarter of their width."""
    def _label(cax, text):
        cax.text(-0.06, 0.5, text, transform=cax.transAxes, ha="right",
                 va="center", fontsize=FS_SMALL, color=INK)

    # The full [0, 1] ice map. The panels paint ice only where phi >= 0.5
    # (below that the vapour shows through), so the dark lower half of the
    # bar is not seen in them; it is the phase field's full range.
    cb = fig.colorbar(ScalarMappable(cmap=cmocean.cm.ice, norm=plt.Normalize(0.0, 1.0)),
                      cax=cax_ice, orientation="horizontal",
                      ticks=[0.0, 0.5, 1.0], format="%g")
    _label(cax_ice, r"$\phi_i$")
    cb.ax.tick_params(labelsize=FS_TINY, width=0.5, length=2, pad=1.5)
    cb.outline.set_linewidth(0.5)

    # Thinned harder than the vertical bar did: horizontally each label is
    # as wide as three or four characters, not as tall as one.
    cb = fig.colorbar(ScalarMappable(cmap=vapcm, norm=norm), cax=cax_sig,
                      orientation="horizontal", extend=sig_extend, extendfrac=0.04,
                      ticks=sigma_ticks(norm, min_gap=0.16))
    _label(cax_sig, r"$\sigma$ [$\times10^{-4}$]")
    # In mathtext, so a negative tick gets a true minus sign.
    cb.ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _p: f"${v:.2g}$"))
    cb.ax.tick_params(labelsize=FS_TINY, width=0.5, length=2, pad=1.5)
    cb.outline.set_linewidth(0.5)


def build(kind, snaps, x, ys, keys, norm, vapcm, icecm, sig_extend, a):
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
    padx, pady = 0.04, 0.04        # insets <-> axes frame
    gap = 0.06                     # between insets
    s_in = (axw - 2 * padx - (n - 1) * gap) / n
    t_band = 0.18                  # inset titles
    lead = 0.16                    # inset titles -> data band, for the leaders
    head = 0.09                    # headroom for the top marker
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
    cb_h, cb_lab, cb_gap = 0.07, 0.15, 0.04     # bar | its tick labels | to axes
    strip = cb_h + cb_lab + cb_gap
    top = 0.24 if not a.no_title else 0.05
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
    # under the curve at the left of the SSA figure, where it is highest --
    # the upper right is kept for the time arrow.
    f0 = below / axh
    if kind == "time":
        where = dict(loc="lower left", bbox_to_anchor=(0.05, f0 + 0.01))
    else:
        where = dict(loc="lower left", bbox_to_anchor=(0.02, f0 + 0.01))
    ax.legend(fontsize=FS_SMALL, frameon=False, handlelength=2.6, ncol=len(keys),
              columnspacing=1.0, handletextpad=0.5,
              **where)
    if kind == "ssa":
        _time_arrow(ax, x, ys["kiso"])

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
        _snap_title(axi, L, f"{t / DAY:.1f} d")
        px, py = x[row], ymark[row]
        # The letter sits INSIDE the marker: beside it, it lands on k_xx or
        # k_yy, which run within a few percent of k_iso.
        _mark(ax, px, py, L)
        # From just above the inset's title to just below the marker's letter.
        # Starts just ABOVE the title (circled letter + time, ~12 pt tall),
        # so it never strikes through it.
        con = ConnectionPatch(xyA=(0.5, 1.0), xyB=(px, py), coordsB=ax.transData,
                              coordsA=mtransforms.offset_copy(
                                  axi.transAxes, fig=fig, y=15.0, units="points"),
                              color="#a0a0a0", lw=0.6, ls=(0, (3, 2)),
                              zorder=3, shrinkA=0, shrinkB=6)
        fig.add_artist(con)

    # Colour bars: phi_i then sigma, left to right across the axes' width.
    y_cb = bot + axh + cb_gap + cb_lab
    lab_ice, lab_sig, sep = 0.22, 0.82, 0.30      # label widths, gap between
    tail = 0.14                   # sigma's end arrow and last tick label
    w_ice = 0.30 * (axw - lab_ice - lab_sig - sep - tail)
    w_sig = axw - lab_ice - lab_sig - sep - tail - w_ice
    cax_ice = fig.add_axes(F(ml + lab_ice, y_cb, w_ice, cb_h))
    cax_sig = fig.add_axes(F(ml + lab_ice + w_ice + sep + lab_sig, y_cb, w_sig, cb_h))
    _colorbars(fig, cax_ice, cax_sig, norm, vapcm, sig_extend)

    if not a.no_title:
        if kind == "time":
            title = a.title_time or DEFAULTS["title_time"]
        else:
            title = a.title_ssa or DEFAULTS["title_ssa"]
        fig.suptitle(title, x=0.5, y=1.0 - 0.03 / H,
                     ha="center", va="top", fontsize=FS_TITLE, color=INK)
    return fig


def _marked_panel(ax, kind, snaps, x, ys, keys, a):
    """One curve panel of the combined figure: k_eff / k_eff,0 against time
    or SSA, the snapshot instants marked with their letters. y spans the
    data only -- there are no insets to make room for."""
    ax.patch.set_alpha(0.0)
    _curve(ax, x, ys, keys)
    yv = np.concatenate([ys[k] for k in keys])
    ylo, yhi = float(yv.min()), float(yv.max())
    pad = 0.08 * (yhi - ylo)
    ax.set_ylim(ylo - pad, yhi + pad)
    ax.yaxis.set_major_locator(MaxNLocator(4, steps=[1, 2, 2.5, 5, 10]))
    xpad = 0.03 * (x.max() - x.min())
    ax.set_xlim(x.min() - xpad, x.max() + xpad)
    ax.plot(ax.get_xlim(), [1.0, 1.0], color="#999999", lw=0.6, ls=":", zorder=0)
    ax.set_ylabel(a.ylabel or DEFAULTS["ylabel_norm"], fontsize=FS, labelpad=3)
    ymark = ys["kiso"]
    for i, (_fl, _X, _Y, _t, row) in enumerate(snaps):
        _mark(ax, x[row], ymark[row], LETTERS[i])
    if kind == "time":
        # The curve rises left to right, so the lower right is empty.
        ax.legend(fontsize=FS_SMALL, frameon=False, handlelength=2.6,
                  ncol=len(keys), columnspacing=1.0, handletextpad=0.5,
                  loc="lower right")
        ax.set_xlabel(a.xlabel_time or DEFAULTS["xlabel_time"], fontsize=FS,
                      labelpad=2)
    else:
        _time_arrow(ax, x, ys["kiso"])
        ax.set_xlabel(a.xlabel_ssa or DEFAULTS["xlabel_ssa"], fontsize=FS,
                      labelpad=2)


def build_combined(snaps, xt, xs, kn, keys, norm, vapcm, icecm, sig_extend, a):
    """Both curves in one figure, read top to bottom as a story:

        colour bars
        (a)-(d) microstructure snapshots, in time order
        k_eff / k_eff,0  vs  time      lettered markers
        k_eff / k_eff,0  vs  SSA/SSA_0 lettered markers

    No leader lines: the letters tie each marked point to its snapshot.
    """
    n = len(snaps)
    W = a.width_mm * MM
    ml, mr = 0.50, 0.08
    axw = W - ml - mr
    gap = 0.06                              # between snapshots
    s_in = (axw - (n - 1) * gap) / n
    top = 0.24 if not a.no_title else 0.05
    cb_h, cb_lab, cb_gap = 0.07, 0.15, 0.06   # bar | tick labels | to titles
    t_band = 0.19                           # snapshot titles
    g_snap = 0.30                           # snapshots -> time panel (+ label)
    ph = 1.75                               # each curve panel
    g_x = 0.58                              # time x label -> SSA panel (+ label)
    bot = 0.40
    H = top + cb_h + cb_lab + cb_gap + t_band + s_in + g_snap + ph + g_x + ph + bot
    fig = plt.figure(figsize=(W, H))
    F = lambda x0, y0, w, h: (x0 / W, y0 / H, w / W, h / H)

    y_ssa = bot
    y_time = y_ssa + ph + g_x
    y_snap = y_time + ph + g_snap
    for i, (fl, X, Y, t, _row) in enumerate(snaps):
        axi = fig.add_axes(F(ml + i * (s_in + gap), y_snap, s_in, s_in))
        XX, YY = _field(axi, fl, X, Y, norm, vapcm, icecm)
        if i == 0:
            _scalebar(axi, XX, YY)
        _snap_title(axi, LETTERS[i], f"{t / DAY:.1f} d")

    _marked_panel(fig.add_axes(F(ml, y_time, axw, ph)), "time", snaps, xt, kn,
                  keys, a)
    _marked_panel(fig.add_axes(F(ml, y_ssa, axw, ph)), "ssa", snaps, xs, kn,
                  keys, a)

    # Panel labels, flush with the figure's left edge: (a) level with the
    # snapshot titles, (b) and (c) just ABOVE their axes, clear of the top
    # y tick label.
    for lab, y_top in zip(PANELS, (y_snap + s_in + 0.5 * t_band,
                                   y_time + ph + 0.14, y_ssa + ph + 0.14)):
        fig.text(0.02 / W, y_top / H, pplib.bold(f"({lab})"), ha="left",
                 va="center", fontsize=FS, color=INK)

    y_cb = y_snap + s_in + t_band + cb_gap + cb_lab
    lab_ice, lab_sig, sep, tail = 0.22, 0.82, 0.30, 0.14
    w_ice = 0.30 * (axw - lab_ice - lab_sig - sep - tail)
    w_sig = axw - lab_ice - lab_sig - sep - tail - w_ice
    cax_ice = fig.add_axes(F(ml + lab_ice, y_cb, w_ice, cb_h))
    cax_sig = fig.add_axes(F(ml + lab_ice + w_ice + sep + lab_sig, y_cb, w_sig, cb_h))
    _colorbars(fig, cax_ice, cax_sig, norm, vapcm, sig_extend)

    if not a.no_title:
        fig.suptitle(a.title_combined or DEFAULTS["title_time"], x=0.5,
                     y=1.0 - 0.03 / H, ha="center", va="top",
                     fontsize=FS_TITLE, color=INK)
    return fig


C_CLIP_LO, C_CLIP_HI = "#ffd21f", "#39d353"   # contours: below / above the bar


def build_clipmap(snaps, norm, vapcm, icecm, sig_extend, a):
    """Where sigma leaves the colour bar: 2 x 2 snapshots, same bars, with a
    contour at each clipped limit. For pointing at, not for the manuscript."""
    v = float(norm.vmax)
    n = len(snaps)
    ncol = 2 if n > 1 else 1
    nrow = int(np.ceil(n / ncol))
    W = a.width_mm * MM
    m, gap, t_band = 0.06, 0.10, 0.22
    pw = (W - 2 * m - (ncol - 1) * gap) / ncol
    cb_h, cb_lab, cb_gap = 0.07, 0.15, 0.10
    note = 0.22
    H = m + note + cb_h + cb_lab + cb_gap + nrow * (t_band + pw) \
        + (nrow - 1) * gap + m
    fig = plt.figure(figsize=(W, H))
    F = lambda x0, y0, w, h: (x0 / W, y0 / H, w / W, h / H)

    lo_clip = sig_extend in ("min", "both")
    hi_clip = sig_extend in ("max", "both")
    for i, (fl, X, Y, t, row) in enumerate(snaps):
        r, c = divmod(i, ncol)
        y0 = m + (nrow - 1 - r) * (pw + t_band + gap)
        ax = fig.add_axes(F(m + c * (pw + gap), y0, pw, pw))
        XX, YY = _field(ax, fl, X, Y, norm, vapcm, icecm)
        if i == 0:
            _scalebar(ax, XX, YY)
        sig = SIGMA_SCALE * pplib.supersaturation(fl["VaporDensity"], fl["Temperature"])
        pore = fl["IcePhase"] < 0.5
        # Outline the out-of-range REGION (pore AND beyond the bar), not the
        # sigma = -v level: most clipped pixels fill small isolated pores
        # whole, where sigma never crosses -v -- their edge is the ice. The
        # 0.5 contour of the indicator traces both kinds of edge.
        parts = []
        for on, mask, col, lab in (
                (lo_clip, pore & (sig < -v), C_CLIP_LO, "< $-v$"),
                (hi_clip, pore & (sig > v), C_CLIP_HI, "> $+v$")):
            if not on:
                continue
            if mask.any():
                ax.contour(XX, YY, mask.astype(float), levels=[0.5], colors=col,
                           linewidths=0.8)
            parts.append(f"{100 * mask.sum() / pore.sum():.1f}% {lab}")
        ax.set_title(f"{LETTERS[i]}: {t / DAY:.1f} d   " + ",  ".join(parts)
                     + " of pore", fontsize=FS_SMALL, color=INK, pad=3)

    # Same bars as the main figures, in the same strip arrangement.
    y_cb = H - m - cb_h
    lab_ice, lab_sig, sep, tail = 0.22, 0.82, 0.30, 0.14
    span = W - 2 * m
    w_ice = 0.30 * (span - lab_ice - lab_sig - sep - tail)
    w_sig = span - lab_ice - lab_sig - sep - tail - w_ice
    cax_ice = fig.add_axes(F(m + lab_ice, y_cb, w_ice, cb_h))
    cax_sig = fig.add_axes(F(m + lab_ice + w_ice + sep + lab_sig, y_cb, w_sig, cb_h))
    _colorbars(fig, cax_ice, cax_sig, norm, vapcm, sig_extend)

    txt = [f"$v$ = {v:.3g}" + r" $\times10^{-4}$ = min(|min $\sigma$|, |max $\sigma$|)"]
    if lo_clip:
        txt.append(r"yellow outline: pore where $\sigma < -v$")
    if hi_clip:
        txt.append(r"green outline: pore where $\sigma > +v$")
    fig.text(m / W, (H - m - cb_h - cb_lab - 0.05) / H, ";   ".join(txt),
             ha="left", va="top", fontsize=FS_TINY, color=INK)
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
    p.add_argument("--width-mm", type=float, default=170.0,
                   help="final printed width [mm] (default 170, the AGU "
                        "full-page maximum). The figure is "
                        "built at this size, so the fonts are print sizes")
    p.add_argument("--absolute", action="store_true",
                   help="time figure: k_eff in W/m/K instead of k/k_0")
    p.add_argument("--iso-only", action="store_true",
                   help="plot k_iso only, not k_xx and k_yy")
    p.add_argument("--title-time", default=None)
    p.add_argument("--title-ssa", default=None)
    p.add_argument("--title-combined", default=None)
    p.add_argument("--title", action="store_true",
                   help="draw the default titles (off by default: manuscript "
                        "figures carry no title). Any --title-* also turns "
                        "titles on")
    p.add_argument("--no-title", action="store_true",
                   help=argparse.SUPPRESS)       # the old flag; now the default
    p.add_argument("--xlabel-time", default=None)
    p.add_argument("--xlabel-ssa", default=None)
    p.add_argument("--ylabel", default=None, help="y label for both figures")
    p.add_argument("--clip-map", action="store_true",
                   help="also write sigma_out_of_range.{pdf,png}: the snapshots "
                        "with contours where sigma leaves the colour bar")
    p.add_argument("--pair-tol", type=float, default=0.02,
                   help="a snapshot with no k_eff sample on its own step pairs "
                        "with the nearest sample within this fraction of t "
                        "(default 0.02); 0 = same step only")
    p.add_argument("--formats", nargs="+", default=["pdf", "png"])
    p.add_argument("--dpi", type=int, default=400, help="PNG and raster dpi")
    a = p.parse_args(argv)
    a.no_title = not (a.title or a.title_time or a.title_ssa or a.title_combined)

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
    idx, fsteps, ftimes, rows = pick_snapshots(
        files, tmap, d, a.times, a.steps, a.n_snapshots, d["t"][ib],
        t_min=d["t"][ib], tol=a.pair_tol)
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
        row = rows[i]
        # The title gives the IMAGE's time; the marker sits on the sample.
        snaps.append([fl, X, Y, float(t), row])
        s = SIGMA_SCALE * pplib.supersaturation(fl["VaporDensity"], fl["Temperature"])
        pore.append(s[fl["IcePhase"] < 0.5])
        off = (d["t"][row] - t) / t if t > 0 else 0.0
        print(f"  snapshot step {fsteps[i]}: t = {t / DAY:.2f} d; marker on k_eff "
              f"step {int(d['step'][row])} ({off:+.2%} in t)")
    pore = np.concatenate(pore)
    smin, smax = float(pore.min()), float(pore.max())
    v = min(abs(smin), abs(smax))
    if v <= 0.0:                                 # one-signed field touching 0
        v = max(abs(smin), abs(smax))
    lo, hi = -v, v
    sig_extend = {(True, True): "both", (True, False): "min",
                  (False, True): "max", (False, False): "neither"}[
                      (smin < lo, smax > hi)]
    norm = AsinhNorm(linear_width=max(v / 300.0, 1e-12), vmin=lo, vmax=hi)
    vapcm = centered_cmap(cmocean.cm.balance, norm)
    icecm = ice_alpha_cmap()

    keys = ("kiso",) if a.iso_only else ("kxx", "kyy", "kiso")
    kn = {k: d[k] / d[k][ib] for k in ("kxx", "kyy", "kiso")}
    out = Path(a.save_dir) if a.save_dir else run / "plots" / "keff" / "snapshots"
    os.makedirs(out, exist_ok=True)

    plt.rcParams.update(pplib.MANUSCRIPT_RC)
    figs = {
        "keff_time_snapshots": build(
            "time", snaps, d["t"] / DAY, d if a.absolute else kn,
            keys, norm, vapcm, icecm, sig_extend, a),
        "keff_ssa_snapshots": build(
            "ssa", snaps, d["ssa"] / d["ssa"][ib], kn,
            keys, norm, vapcm, icecm, sig_extend, a),
    }
    # Both in one, three panels top to bottom: snapshots, k(t), k(SSA).
    # Always normalized -- the two curve panels share one y quantity.
    figs["keff_combined_snapshots"] = build_combined(
        snaps, d["t"] / DAY, d["ssa"] / d["ssa"][ib], kn, keys, norm, vapcm,
        icecm, sig_extend, a)
    if a.clip_map:
        figs["sigma_out_of_range"] = build_clipmap(snaps, norm, vapcm, icecm,
                                                   sig_extend, a)
    for stem, fig in figs.items():
        h_mm = fig.get_figheight() * 25.4
        if h_mm > MAX_H_MM:
            print(f"  WARNING: {stem} is {h_mm:.0f} mm tall, over AGU's "
                  f"{MAX_H_MM:.0f} mm limit", file=sys.stderr)
        for fmt in a.formats:
            path = out / f"{stem}.{fmt}"
            # No bbox_inches="tight": it would crop to the ink and the
            # saved width would no longer be --width-mm.
            # Transparent: the page (or whoever places the figure) supplies
            # the background. The field images themselves stay opaque.
            fig.savefig(path, dpi=a.dpi, transparent=True)
            print(f"  wrote {path}")
        plt.close(fig)
    print(f"  reference (subscript 0) = opening frame: t_0 = {d['t'][ib]:.4g} s, step "
          f"{d['step'][ib]}: k_xx,0 = {d['kxx'][ib]:.4g}, k_yy,0 = {d['kyy'][ib]:.4g}, "
          f"k_iso,0 = {d['kiso'][ib]:.4g} W/m/K, SSA_0 = {d['ssa'][ib]:.4g} 1/m")
    print(f"  sigma x{SIGMA_SCALE:g}: data {smin:+.3g} .. {smax:+.3g} over the shown "
          f"pores; bar +-{v:.3g} (extend={sig_extend})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
