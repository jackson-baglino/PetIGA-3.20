#!/usr/bin/env python3
"""plot_molaro_validation.py — manuscript figures for the Molaro (2019) comparison.

    python3 plot_molaro_validation.py --run-t20 <run> [--run-t5 <run>]
        [--steps S0 S1] [--crop-um Z0 Z1 RMAX] [--width-mm 170]
        [--save-dir DIR] [--formats pdf png]

Three page-width figures, built to match plot_keff_snapshots.py so the two
sets read as one: the same 170 mm width, print type sizes, ink and muted
greys, circled instant numbers, bold panel letters, cmocean `ice` over
cmocean `balance`, horizontal colour-bar strip, no titles, transparent
background. Its drawing helpers are IMPORTED from there, not copied, so
the sets cannot drift apart.

    molaro_neck_width.{pdf,png}       neck width vs time, model (lines) and
                                      Molaro et al. (2019) (points with their
                                      error bars), both temperatures
    molaro_microstructure.{pdf,png}   the -20 C meridional section at
                                      Molaro's t = 0 and at the end of their
                                      record, instants 1 and 2
    molaro_combined.{pdf,png}         (a) the sections, (b) the neck curves,
                                      instants 1-2 marked on the -20 C model
    molaro_grain_shrinkage.{pdf,png}  D / D_0 of (a) the large and (b) the
                                      small grain, both temperatures
    Figure3_molaro_validation.{pdf,png}             molaro_combined with the shrinkage
                                      under it as (c) and (d)

THE CLOCK. Molaro's record starts at an unknown time after contact, and our
runs start from a chosen r = 14 um neck. So each series' t = 0 is the moment
its neck first reaches Molaro's first measured width (32.81 um at -20 C,
32.51 um at -5 C) -- their Fig. 12 convention, and the one
plot_neck_vs_molaro.py and run_batch_measure.sh already use. Nothing before
t = 0 is drawn: every model curve opens on its sample nearest t*, and the
grain diameters are normalised by their values there (D_0).

THE INSTANTS. --instants-min (default 0 78) are anchored times; each picks
the nearest snapshot. 78 min is Molaro's last -20 C point, the end of the
window both the neck and the shrinkage are scored over. Only
snapshots with a neck_width.csv sample are eligible, so both markers sit on
measured values; the titles give each snapshot's own anchored time.

THE SECTIONS are mirrored across the symmetry axis, so the pair reads as two
grains rather than two half-discs, and cropped to the ice with a margin
(--crop-um overrides). Horizontal is the symmetry axis z, vertical r.

SIGMA RANGE. The k_eff figures' rule: symmetric about zero so sigma = 0 is
the bar's middle, sized to the SMALLER extreme. Here that is the neck's
+0.28 (x1e-4); the undersaturated far field (down to -27) is beyond the
bar and drawn in its end colour (triangular cap). The printout gives the
share of pore that clips. Shared by both snapshots.

MISSING -5 C RUN. Without --run-t5 the -5 C data are drawn with no model
curve; the figure is otherwise complete, and the note says so.

OUTPUT. studies/molaro_2019/manuscript/ (or --save-dir). PDF is the
manuscript file, PNG a preview. --copy-to also copies every file into the
manuscript's Figures/Figure<N>__<name>/ folder (Figure 2 for these).
"""
from __future__ import annotations

import argparse
import glob
import os
import shutil
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import cmocean

HERE = Path(__file__).parent
REPO = HERE.parent
sys.path.insert(0, str(HERE))
from pplib import step_times                                # noqa: E402
from plot_keff_snapshots import (INK, FS, FS_SMALL, MM, MAX_H_MM,  # noqa: E402
                                 LETTERS, PANELS, make_reader, snap_step,
                                 _field, _scalebar, _snap_title, _circled,
                                 _colorbars)
from make_neck_movie import (SIGMA_SCALE, ice_alpha_cmap,   # noqa: E402
                             centered_cmap, read_experiment)
from plot_neck_vs_molaro import read_model, anchor_time     # noqa: E402
import pplib                                                # noqa: E402

UM = 1e-6
WANT = ("IcePhase", "VaporDensity", "Temperature")
# One colour per TEMPERATURE, sampled from cmocean `balance` -- the map the
# sections' sigma already uses -- at 0.25 (blue, -20 C) and 0.75 (red,
# -5 C): cold blue, warm red, and the figure keeps one palette. 1/3 and 2/3
# were tried and are too pale (2.6:1 on white, protan dE 10.8); at
# 0.25/0.75 contrast is 3.9 / 4.5 : 1 and CVD dE >= 19.0.
# Line vs point is model vs experiment, so colour carries only the temperature.
C_T20, C_T5 = "#3888ba", "#bf573a"
SERIES = {
    "T-20": dict(label=r"$-20\,^{\circ}$C", color=C_T20, anchor_um=32.81, window_min=78.0,
                 data=REPO / "inputs/validation/molaro2019_fig11_T-20.csv"),
    "T-5":  dict(label=r"$-5\,^{\circ}$C",  color=C_T5, anchor_um=32.51, window_min=48.0,
                 data=REPO / "inputs/validation/molaro2019_fig11_T-5.csv"),
}
MAX_PX = 1800                       # raster columns per section, ~550 dpi


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
def load_model(run, anchor_um):
    """One run's anchored neck and grain curves. Times in minutes."""
    tm, wm = read_model(run)
    t_star = anchor_time(tm, wm, anchor_um * UM)
    if t_star is None:
        sys.exit(f"  {run}: the model neck ({wm.min()/UM:.2f}-{wm.max()/UM:.2f} um) "
                 f"never crosses the {anchor_um} um anchor")
    # Nothing before t = 0 is drawn. The curve opens on the sample nearest
    # t* -- a measured value, never an interpolated one -- and that is also
    # the sample instant 1 lands on.
    i0 = int(np.argmin(np.abs(tm - t_star)))
    tm, wm = tm[i0:], wm[i0:]
    return dict(t_s=tm, w_m=wm, t_star=t_star, t=(tm - t_star) / 60.0,
                w=wm / UM, grains=read_grains(run, tm[0], t_star))


def load_series(key, run, alt=None):
    """Anchored model(s) and experiment for one temperature. `alt` is a
    second run of the same temperature, drawn dashed (ALT_LS)."""
    s = dict(SERIES[key], key=key, run=run)
    td, wd, ep, em = read_experiment(s["data"])
    s["td"] = (td - td[0]) / 60.0
    s["wd"], s["ep"], s["em"] = wd / UM, ep / UM, em / UM
    s["model"] = load_model(run, s["anchor_um"]) if run is not None else None
    s["alt"] = load_model(alt, s["anchor_um"]) if alt is not None else None
    s["Dd"] = read_data_diameters(s["data"])
    return s


ALT_LS = (0, (4, 2.2))                   # the second run of a temperature


FIT_LS = (0, (3, 2))                     # least-squares line through the data


def legend_handles(series, data_fit=False):
    """Colour = temperature; line = model, dashed = its alternative wall;
    open circle = Molaro."""
    h = [Line2D([], [], color=s["color"], lw=1.8, label=s["label"]) for s in series]
    h.append(Line2D([], [], color=INK, lw=1.8, label="model"))
    alt = [s for s in series if s["alt"] is not None]
    if alt:
        h.append(Line2D([], [], color=INK, lw=1.8, ls=ALT_LS,
                        label=alt[0].get("alt_label", "model, alt. wall")))
    # "data", not the citation: the caption names Molaro et al. (2019).
    h.append(Line2D([], [], color=INK, ls="none", marker="o", ms=4.2, mfc="white",
                    mew=0.9, label="data"))
    if data_fit:
        h.append(Line2D([], [], color=INK, lw=1.0, ls=FIT_LS, label="fit"))
    return h


def _time_axis_from_zero(ax, xmax):
    """The time axis starts AT t = 0 -- nothing precedes Molaro's first
    measurement. The t = 0 points, their error bars and the instant markers
    are drawn unclipped, so they sit whole on top of the y axis rather than
    being cut by it (a padded axis read as time before t = 0)."""
    ax.set_xlim(0.0, 1.06 * xmax)


def read_grains(run, t_open, t_star):
    """Model grain diameters from grain_shrinkage.csv, from the opening sample
    on, normalised by their own values there. None if the CSV is missing."""
    f = run / "grain_shrinkage.csv"
    if not f.is_file():
        print(f"  {run.name}: no grain_shrinkage.csv -- run "
              f"postprocess/grain_shrinkage.py; no model shrinkage curve")
        return None
    g = np.genfromtxt(f, delimiter=",", names=True)
    keep = g["t_s"] >= t_open - 1e-9
    t = g["t_s"][keep]
    D_lg, D_sm = 2 * g["R_large_m"][keep], 2 * g["R_small_m"][keep]
    return dict(t=(t - t_star) / 60.0, large=D_lg / D_lg[0], small=D_sm / D_sm[0],
                D0=(D_lg[0] / UM, D_sm[0] / UM))


def read_data_diameters(path):
    """(t_min, D_large/D_large0, D_small/D_small0) from a Fig. 11 CSV."""
    rows = [ln.split(",") for ln in open(path)
            if ln.strip() and not ln.startswith("#")]
    a = np.array([[float(x) for x in r[:6]] for r in rows])
    return dict(t=a[:, 0] - a[0, 0], large=a[:, 4] / a[0, 4], small=a[:, 5] / a[0, 5])


def pick_instants(s, steps, instants_min):
    """[(file, step, row)] for instants 1 and 2, on neck_width.csv samples."""
    run, m = s["run"], s["model"]
    files, _ = make_reader(run, "sol" if glob.glob(str(run / "sol_*.dat")) else "vts")
    tmap = step_times(str(run))
    fsteps = np.array([snap_step(f) for f in files])
    ftimes = np.array([tmap.get(int(k), np.nan) for k in fsteps])
    # A snapshot is eligible only if the neck was measured at its time.
    row_of = {}
    for i, t in enumerate(ftimes):
        if np.isfinite(t):
            j = int(np.argmin(np.abs(m["t_s"] - t)))
            if abs(m["t_s"][j] - t) <= 1e-6 * max(1.0, abs(t)):
                row_of[i] = j
    if not row_of:
        sys.exit(f"  {run}: no snapshot has a neck_width.csv sample")
    cand = np.array(sorted(row_of))
    if steps:
        want = [int(cand[np.argmin(np.abs(fsteps[cand] - k))]) for k in steps]
    else:
        targets = [m["t_star"] + 60.0 * x for x in instants_min]
        want = [int(cand[np.argmin(np.abs(ftimes[cand] - t))]) for t in targets]
    return [(files[i], int(fsteps[i]), row_of[i]) for i in want]


def _mirror(fl, X, Y):
    """Reflect across the axis r = 0 (row 0), dropping the duplicate row."""
    m = lambda a: np.vstack([a[:0:-1], a])
    return {k: m(v) for k, v in fl.items()}, m(X), np.vstack([-Y[:0:-1], Y])


def read_section(run, fn, stride):
    """Fields mirrored across the axis (rows = r, cols = z), every `stride`."""
    _, reader = make_reader(run, "sol" if fn.endswith(".dat") else "vts")
    fl, X, Y = _mirror(*reader(fn, want=WANT))
    sub = lambda a: a[::stride, ::stride]
    return {k: sub(v) for k, v in fl.items()}, sub(X), sub(Y)


def crop(fl, X, Y, box_um):
    z0, z1, rmax = (b * UM for b in box_um)
    cols = (X[0] >= z0) & (X[0] <= z1)
    rows = np.abs(Y[:, 0]) <= rmax
    sub = lambda a: a[np.ix_(rows, cols)]
    return {k: sub(v) for k, v in fl.items()}, sub(X), sub(Y)


def auto_crop(fl, X, Y, margin=0.10):
    """The ice's bounding box, grown by `margin` of its axial length."""
    ice = fl["IcePhase"] >= 0.5
    zs, rs = X[ice], np.abs(Y[ice])
    pad = margin * (zs.max() - zs.min())
    return ((zs.min() - pad) / UM, (zs.max() + pad) / UM, (rs.max() + pad) / UM)


def load_sections(s, instants, box_um):
    run, fn0 = s["run"], instants[0][0]
    # Probe the grid once, for the crop and the stride.
    _, reader = make_reader(run, "sol" if fn0.endswith(".dat") else "vts")
    fl, X, Y = _mirror(*reader(fn0, want=("IcePhase",)))
    if box_um is None:
        box_um = auto_crop(fl, X, Y)
    frac = (box_um[1] - box_um[0]) * UM / (X.max() - X.min())
    stride = max(1, int(np.ceil(frac * X.shape[1] / MAX_PX)))
    out = []
    for fn, step, row in instants:
        out.append(crop(*read_section(run, fn, stride), box_um) + (step, row))
    return out, box_um, stride


# ---------------------------------------------------------------------------
# Drawing
# ---------------------------------------------------------------------------
def _neck_panel(ax, series, marks=()):
    """Neck width vs anchored time: model lines, Molaro's points."""
    ax.patch.set_alpha(0.0)
    xmax, ylo, yhi = 0.0, np.inf, -np.inf
    for s in series:
        c = s["color"]
        if s["model"] is not None:
            m = s["model"]
            ax.plot(m["t"], m["w"], "-", lw=1.8, color=c, zorder=2)
            vis = m["t"] <= SERIES["T-20"]["window_min"] * 1.08
            ylo, yhi = min(ylo, m["w"][vis].min()), max(yhi, m["w"][vis].max())
        if s["alt"] is not None:
            m2 = s["alt"]
            ax.plot(m2["t"], m2["w"], ls=ALT_LS, lw=1.5, color=c, zorder=2,
                    dash_capstyle="round")
        ax.errorbar(s["td"], s["wd"], yerr=[s["em"], s["ep"]], fmt="o", ms=4.2,
                    mfc="white", mec=c, mew=0.9, ecolor=c, elinewidth=0.7,
                    capsize=1.8, capthick=0.7, zorder=3, clip_on=False)
        xmax = max(xmax, float(s["td"].max()))
        ylo = min(ylo, float((s["wd"] - s["em"]).min()))
        yhi = max(yhi, float((s["wd"] + s["ep"]).max()))
    for x, y, lab in marks:
        _circled(ax, x, y, lab, clip_on=False)
    _time_axis_from_zero(ax, xmax)
    ypad = 0.08 * (yhi - ylo)
    ax.set_ylim(ylo - ypad, yhi + ypad)
    ax.xaxis.set_major_locator(MaxNLocator(8, steps=[1, 2, 2.5, 5, 10]))
    ax.yaxis.set_major_locator(MaxNLocator(5, steps=[1, 2, 2.5, 5, 10]))
    ax.tick_params(labelsize=FS_SMALL, width=0.6, length=3, pad=2)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_linewidth(0.6)
    ax.set_xlabel("Time [min]", fontsize=FS, labelpad=2)
    ax.set_ylabel(r"$w$  [$\mu$m]", fontsize=FS, labelpad=3)
    # Colour is the temperature, mark style the source: two short columns.
    h = legend_handles(series)
    ax.legend(handles=h, fontsize=FS_SMALL, frameon=False, handlelength=2.2,
              ncol=2, columnspacing=1.2, handletextpad=0.5, loc="lower right")


def _sections(fig, F, secs, x0, y0, s_w, s_h, gap, norm, vapcm, icecm, t_star):
    for i, (fl, X, Y, step, row, t) in enumerate(secs):
        axi = fig.add_axes(F(x0 + i * (s_w + gap), y0, s_w, s_h))
        XX, YY = _field(axi, fl, X, Y, norm, vapcm, icecm)
        axi.set_aspect("auto")
        if i == 0:
            _scalebar(axi, XX, YY)
        tm = (t - t_star) / 60.0
        _snap_title(axi, LETTERS[i], f"{round(tm) + 0:d} min")   # +0: no "-0"


def _strip(fig, F, x0, y_cb, axw, norm, vapcm, sig_extend):
    cb_h = 0.07
    lab_ice, lab_sig, sep, tail = 0.22, 0.82, 0.30, 0.14
    w_ice = 0.30 * (axw - lab_ice - lab_sig - sep - tail)
    w_sig = axw - lab_ice - lab_sig - sep - tail - w_ice
    cax_ice = fig.add_axes(F(x0 + lab_ice, y_cb, w_ice, cb_h))
    cax_sig = fig.add_axes(F(x0 + lab_ice + w_ice + sep + lab_sig, y_cb, w_sig, cb_h))
    _colorbars(fig, cax_ice, cax_sig, norm, vapcm, sig_extend)


def build_neck(series, a):
    W = a.width_mm * MM
    ml, mr, top, bot, ph = 0.50, 0.08, 0.08, 0.40, 2.30
    H = top + ph + bot
    fig = plt.figure(figsize=(W, H))
    F = lambda x0, y0, w, h: (x0 / W, y0 / H, w / W, h / H)
    _neck_panel(fig.add_axes(F(ml, bot, W - ml - mr, ph)), series)
    return fig


def _geom(W, secs, gap):
    ml, mr = 0.50, 0.08
    axw = W - ml - mr
    s_w = (axw - gap) / 2
    fl, X, Y = secs[0][:3]
    s_h = s_w * (Y.max() - Y.min()) / (X.max() - X.min())
    return ml, axw, s_w, s_h


def build_micro(secs, norm, vapcm, icecm, t_star, a):
    W = a.width_mm * MM
    gap = 0.12
    ml, axw, s_w, s_h = _geom(W, secs, gap)
    cb_h, cb_lab, cb_gap, t_band, top, bot = 0.07, 0.15, 0.06, 0.19, 0.05, 0.05
    H = top + cb_h + cb_lab + cb_gap + t_band + s_h + bot
    fig = plt.figure(figsize=(W, H))
    F = lambda x0, y0, w, h: (x0 / W, y0 / H, w / W, h / H)
    _sections(fig, F, secs, ml, bot, s_w, s_h, gap, norm, vapcm, icecm, t_star)
    _strip(fig, F, ml, bot + s_h + t_band + cb_gap + cb_lab, axw, norm, vapcm, a.sig_extend)
    return fig


def build_combined(secs, series, marks, norm, vapcm, icecm, t_star, a):
    """Read top to bottom: colour bars, (a) the sections at instants 1-2,
    (b) the neck curves with the same instants marked on the -20 C model."""
    W = a.width_mm * MM
    gap = 0.12
    ml, axw, s_w, s_h = _geom(W, secs, gap)
    cb_h, cb_lab, cb_gap, t_band, top = 0.07, 0.15, 0.06, 0.19, 0.05
    g_snap, ph, bot = 0.34, 2.30, 0.40
    H = top + cb_h + cb_lab + cb_gap + t_band + s_h + g_snap + ph + bot
    fig = plt.figure(figsize=(W, H))
    F = lambda x0, y0, w, h: (x0 / W, y0 / H, w / W, h / H)
    y_snap = bot + ph + g_snap
    _sections(fig, F, secs, ml, y_snap, s_w, s_h, gap, norm, vapcm, icecm, t_star)
    _neck_panel(fig.add_axes(F(ml, bot, axw, ph)), series, marks)
    for lab, y_top in zip(PANELS, (y_snap + s_h + 0.5 * t_band, bot + ph + 0.14)):
        fig.text(0.02 / W, y_top / H, pplib.bold(f"({lab})"), ha="left",
                 va="center", fontsize=FS, color=INK)
    _strip(fig, F, ml, y_snap + s_h + t_band + cb_gap + cb_lab, axw, norm, vapcm, a.sig_extend)
    return fig


def build_full(secs, series, marks, norm, vapcm, icecm, t_star, a):
    """The combined figure with the grain shrinkage under it: (a) sections,
    (b) neck curves, (c) large grain D / D_0, (d) small grain. One figure for
    the whole comparison -- the shrinkage is what shows the model loses mass
    at the measured rate even where the neck curves part from the data. The
    neck panel is lower than in molaro_combined to stay under the page."""
    W = a.width_mm * MM
    gap = 0.12
    ml, axw, s_w, s_h = _geom(W, secs, gap)
    cb_h, cb_lab, cb_gap, t_band, top = 0.07, 0.15, 0.06, 0.19, 0.05
    g_snap, ph, g_row, leg, ph2, bot = 0.34, 1.75, 0.62, 0.0, 1.55, 0.40
    pgap, ml2 = 0.78, ml + 0.17          # room for the 4-digit tick labels
    pw = (W - ml2 - 0.08 - pgap) / 2
    H = top + cb_h + cb_lab + cb_gap + t_band + s_h + g_snap + ph + g_row + leg + ph2 + bot
    fig = plt.figure(figsize=(W, H))
    F = lambda x0, y0, w, h: (x0 / W, y0 / H, w / W, h / H)
    y_neck = bot + ph2 + leg + g_row
    y_snap = y_neck + ph + g_snap
    _sections(fig, F, secs, ml, y_snap, s_w, s_h, gap, norm, vapcm, icecm, t_star)
    _neck_panel(fig.add_axes(F(ml, y_neck, axw, ph)), series, marks)
    xmax = max(float(s["Dd"]["t"].max()) for s in series)
    for i, (which, sym) in enumerate((("large", "l"), ("small", "s"))):
        ax = fig.add_axes(F(ml2 + i * (pw + pgap), bot, pw, ph2))
        _shrink_panel(ax, series, which, xmax)
        ax.set_ylabel(rf"$D_\mathrm{{{sym}}}\,/\,D_{{\mathrm{{{sym}}},0}}$",
                      fontsize=FS, labelpad=3)
        if i == 1:
            ax.legend(handles=[Line2D([], [], ls=FIT_LS, lw=1.0, color=INK, dash_capstyle="round",
                                      label="linear fit to data")], fontsize=FS_SMALL, frameon=False,
                      handlelength=2.2, loc="lower left", borderaxespad=0.3)
        fig.text((0.02 + i * (pw + pgap + (ml2 - 0.02) * 0)) / W + i * (ml2 - 0.72) / W, (bot + ph2 + 0.16) / H,
                 pplib.bold("(%s)" % "cd"[i]), ha="left", va="center",
                 fontsize=FS, color=INK)
    for lab, y_top in zip(PANELS, (y_snap + s_h + 0.5 * t_band, y_neck + ph + 0.14)):
        fig.text(0.02 / W, y_top / H, pplib.bold(f"({lab})"), ha="left",
                 va="center", fontsize=FS, color=INK)
    _strip(fig, F, ml, y_snap + s_h + t_band + cb_gap + cb_lab, axw, norm, vapcm, a.sig_extend)
    return fig


def _shrink_panel(ax, series, which, xmax):
    """D / D_0 of one grain: model lines, Molaro's points (no error bars --
    their table gives none for the diameters)."""
    ax.patch.set_alpha(0.0)
    vals = []
    for s in series:
        c, m = s["color"], s["model"]
        if m is not None and m["grains"] is not None:
            g = m["grains"]
            ax.plot(g["t"], g[which], "-", lw=1.8, color=c, zorder=2)
            vals.append(g[which][g["t"] <= xmax])
        if s["alt"] is not None and s["alt"]["grains"] is not None:
            g2 = s["alt"]["grains"]
            ax.plot(g2["t"], g2[which], ls=ALT_LS, lw=1.5, color=c, zorder=2,
                    dash_capstyle="round")
            vals.append(g2[which][g2["t"] <= xmax])
        d = s["Dd"]
        ax.plot(d["t"], d[which], "o", ms=4.2, mfc="white", mec=c, mew=0.9,
                ls="none", zorder=3, clip_on=False)
        # Least-squares line through the data (intercept free), over the span
        # the data cover only: the trend the wall humidity is scored against
        # (-2.93 %/78 min for the -20 C large grain).
        k, b = np.polyfit(d["t"], d[which], 1)
        tt = np.array([d["t"].min(), d["t"].max()])
        ax.plot(tt, k * tt + b, ls=FIT_LS, lw=1.0, color=c, zorder=2.5,
                dash_capstyle="round")
        print(f"  {s['label']} {which} grain, data fit: {100 * k * tt[1]:+.2f} % over "
              f"{tt[1]:.0f} min (slope {100 * k:+.4f} %/min)")
        vals.append(d[which])
    v = np.concatenate(vals)
    pad = 0.08 * (v.max() - v.min())
    ax.set_ylim(v.min() - pad, max(v.max(), 1.0) + pad)
    _time_axis_from_zero(ax, xmax)
    ax.xaxis.set_major_locator(MaxNLocator(5, steps=[1, 2, 2.5, 5, 10]))
    ax.yaxis.set_major_locator(MaxNLocator(5, steps=[1, 2, 2.5, 5, 10]))
    ax.tick_params(labelsize=FS_SMALL, width=0.6, length=3, pad=2)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_linewidth(0.6)
    ax.set_xlabel("Time [min]", fontsize=FS, labelpad=2)


def build_shrinkage(series, a):
    """(a) the large grain, (b) the small one, each D / D_0 against the
    anchored clock. Colour is the temperature, as on the neck figure."""
    W = a.width_mm * MM
    ml, mr, gap, top, bot, ph = 0.62, 0.08, 0.62, 0.42, 0.40, 2.10
    pw = (W - ml - mr - gap) / 2
    H = top + ph + bot
    fig = plt.figure(figsize=(W, H))
    F = lambda x0, y0, w, h: (x0 / W, y0 / H, w / W, h / H)
    xmax = max(float(s["Dd"]["t"].max()) for s in series)
    for i, (which, sym) in enumerate((("large", "l"), ("small", "s"))):
        ax = fig.add_axes(F(ml + i * (pw + gap), bot, pw, ph))
        _shrink_panel(ax, series, which, xmax)
        ax.set_ylabel(rf"$D_\mathrm{{{sym}}}\,/\,D_{{\mathrm{{{sym}}},0}}$",
                      fontsize=FS, labelpad=3)
        fig.text((ml + i * (pw + gap) - ml + 0.02) / W, (bot + ph + 0.10) / H,
                 pplib.bold(f"({PANELS[i]})"), ha="left", va="center",
                 fontsize=FS, color=INK)
    # One row above both panels: inside them every corner holds data.
    h = legend_handles(series, data_fit=True)
    fig.legend(handles=h, fontsize=FS_SMALL, frameon=False, handlelength=2.2,
               ncol=len(h), columnspacing=1.6, handletextpad=0.5,
               loc="center", bbox_to_anchor=(0.5, (bot + ph + top - 0.08) / H))
    return fig


# ---------------------------------------------------------------------------
def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run-t20", type=Path, required=True,
                   help="-20 C run directory (neck_width.csv + snapshots)")
    p.add_argument("--run-t5", type=Path, default=None,
                   help="-5 C run directory; omitted, its data are drawn alone")
    p.add_argument("--run-t5-alt", type=Path, default=None,
                   help="a second -5 C run, drawn dashed in the same colour "
                        "(the refit-wall arm)")
    p.add_argument("--alt-label", default=r"model, $-5\,^{\circ}$C wall refit",
                   help="legend text for the dashed run")
    p.add_argument("--steps", type=int, nargs=2, default=None,
                   help="snapshot steps for instants 1 and 2, overriding "
                        "--instants-min")
    p.add_argument("--instants-min", type=float, nargs=2, default=[0.0, 78.0],
                   help="anchored times [min] of instants 1 and 2; the nearest "
                        "measured snapshot is used (default 0 78: Molaro's "
                        "first and last -20 C points)")
    p.add_argument("--copy-to", type=Path, default=None,
                   help="also copy every figure here -- the manuscript's "
                        "Figures/Figure<N>__<name>/ folder")
    p.add_argument("--crop-um", type=float, nargs=3, default=None,
                   metavar=("Z0", "Z1", "RMAX"),
                   help="section window [um] (default: the ice + 10 %%)")
    p.add_argument("--width-mm", type=float, default=170.0)
    p.add_argument("--save-dir", type=Path, default=REPO / "studies/molaro_2019/manuscript")
    p.add_argument("--formats", nargs="+", default=["pdf", "png"])
    p.add_argument("--dpi", type=int, default=600,
                   help="PNG and raster dpi (default 600, the top of AGU's 300-600 ppi)")
    a = p.parse_args(argv)

    s20 = load_series("T-20", a.run_t20.resolve())
    s5 = load_series("T-5", a.run_t5.resolve() if a.run_t5 else None,
                     alt=a.run_t5_alt.resolve() if a.run_t5_alt else None)
    s5["alt_label"] = a.alt_label
    series = [s20, s5]
    for s in series:
        m = s["model"]
        if m is None:
            print(f"  {s['label']}: no run given -- Molaro's points only")
            continue
        wend = np.interp(m["t_star"] + s["window_min"] * 60.0, m["t_s"], m["w_m"]) / UM
        print(f"  {s['label']}: t* = {m['t_star']:.1f} s; w(t*+{s['window_min']:.0f} min)"
              f" = {wend:.2f} um vs Molaro {s['wd'][-1]:.2f} um "
              f"(growth {100*(wend-s['anchor_um'])/(s['wd'][-1]-s['anchor_um']):.0f} %)")

    inst = pick_instants(s20, a.steps, a.instants_min)
    secs, box, stride = load_sections(s20, inst, a.crop_um)
    m = s20["model"]
    secs = [sec + (float(m["t_s"][sec[4]]),) for sec in secs]
    marks = [(m["t"][row], m["w"][row], LETTERS[i])
             for i, (_f, _s, row) in enumerate(inst)]
    for i, (_f, step, row) in enumerate(inst):
        print(f"  instant {LETTERS[i]}: step {step}, t = {m['t_s'][row]:.1f} s "
              f"(t - t* = {m['t'][row]:+.2f} min), w = {m['w'][row]:.2f} um")
    print(f"  section window z {box[0]:.1f}-{box[1]:.1f} um, |r| <= {box[2]:.1f} um, "
          f"stride {stride}")

    pore = np.concatenate([SIGMA_SCALE * pplib.supersaturation(
        fl["VaporDensity"], fl["Temperature"])[fl["IcePhase"] < 0.5]
        for fl, *_ in secs])
    smin, smax = float(pore.min()), float(pore.max())
    # plot_keff_snapshots' rule: symmetric about 0, sized to the SMALLER
    # extreme, so sigma = 0 is the bar's middle; the other end saturates.
    v = min(abs(smin), abs(smax))
    if v <= 0.0:
        v = max(abs(smin), abs(smax))
    sig_extend = {(True, True): "both", (True, False): "min",
                  (False, True): "max", (False, False): "neither"}[(smin < -v, smax > v)]
    norm = AsinhNorm(linear_width=max(v / 300.0, 1e-12), vmin=-v, vmax=v)
    vapcm = centered_cmap(cmocean.cm.balance, norm)
    icecm = ice_alpha_cmap()
    clipped = 100.0 * np.mean((pore < -v) | (pore > v))
    print(f"  sigma x{SIGMA_SCALE:g} over the shown pore: {smin:+.3g} .. {smax:+.3g}; "
          f"bar +-{v:.3g} (extend={sig_extend}, {clipped:.0f} % of pore beyond it)")

    a.sig_extend = sig_extend
    plt.rcParams.update(pplib.MANUSCRIPT_RC)
    figs = {
        "molaro_neck_width": build_neck(series, a),
        "molaro_microstructure": build_micro(secs, norm, vapcm, icecm, m["t_star"], a),
        "molaro_combined": build_combined(secs, series, marks, norm, vapcm, icecm,
                                          m["t_star"], a),
        "molaro_grain_shrinkage": build_shrinkage(series, a),
        "Figure3_molaro_validation": build_full(secs, series, marks, norm, vapcm, icecm, m["t_star"], a),
    }
    os.makedirs(a.save_dir, exist_ok=True)
    for stem, fig in figs.items():
        h_mm = fig.get_figheight() * 25.4
        if h_mm > MAX_H_MM:
            print(f"  WARNING: {stem} is {h_mm:.0f} mm tall, over AGU's "
                  f"{MAX_H_MM:.0f} mm limit", file=sys.stderr)
        for fmt in a.formats:
            path = a.save_dir / f"{stem}.{fmt}"
            fig.savefig(path, dpi=a.dpi, transparent=True)
            print(f"  wrote {path}")
            if a.copy_to is not None:
                a.copy_to.mkdir(parents=True, exist_ok=True)
                shutil.copy2(path, a.copy_to / path.name)
        plt.close(fig)
    if a.copy_to is not None:
        print(f"  copied to {a.copy_to}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
