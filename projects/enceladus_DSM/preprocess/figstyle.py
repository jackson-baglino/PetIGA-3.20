"""figstyle.py — the figure machinery shared by the constraint figure sets.

Extracted so that two families of plots (Libbrecht sigma0(T), and constant
alpha_c) cannot drift apart in type scale, legend grouping, off-scale marking
or width policy. Nothing in here knows about ice; it knows about slides.

THE POLICY IT ENFORCES
----------------------
  * Width is a hard ceiling (SLIDE_W_IN): these go on PowerPoint slides, and an
    image wider than the slide is useless however good it looks. `save()`
    asserts it rather than trusting a figsize. Height is free.
  * No type below 10 pt, at either size.
  * Legends live UNDER the axes, never over a curve, and are laid out as
    titled columns -- a flat strip of entries is a lookup problem.
  * A panel narrower than COMPACT_BELOW_IN is drawn in COMPACT mode: short
    title, short y-label, no prose annotations, no in-frame legend. That is a
    flag the panel functions consult, not a second set of drawing code, so a
    small panel and its full-size twin cannot disagree.
  * A curve clipped by the display window is marked on the frame with the
    extreme it reaches, so clipping never reads as missing data.
"""

from __future__ import annotations

import math

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Okabe-Ito, the repo standard (plot_alpha_kinetics.py, fit_neck_growth.py).
C = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"]
INK, MUTED, GRID = "#1a1a1a", "#5c5c5c", "#d8d8d8"

BOUND_COLOR = {"B-HEAT": C[0], "B-VAPOR": C[5], "B-KINETIC": C[1], "B-CURV": C[2]}
# Short, because these ride in a legend under the axes. The formulas they used
# to carry live in docs/tex/constraints_iguanatex.txt.
BOUND_LABEL = {
    "B-HEAT":    "B-HEAT  (K&P 43a/b)",
    "B-VAPOR":   "B-VAPOR  (K&P 43c)",
    "B-KINETIC": "B-KINETIC  (K&P 45)",
    "B-CURV":    "B-CURV  (geometric)",
}
GRP_BIND = "Which bound binds"
GRP_REF = "Reference"

# Type scale. Sized for a projector, not for a screen an arm's length away.
FS_TITLE, FS_LABEL, FS_TICK, FS_LEG, FS_NOTE = 15, 13, 11.5, 10, 10
FS_C_TITLE, FS_C_LABEL, FS_C_TICK = 12, 10.5, 10

SLIDE_W_IN = 10.0
ASPECT = 0.58              # 10 x 5 of plot, plus the legend strip underneath
COMPACT_ASPECT = 0.80      # a small panel needs proportionally more height
COMPACT_BELOW_IN = 6.0

_COMPACT = [False]
# Separate from _COMPACT: only a legend SHARED across panels needs several of
# them to agree on one label for their reference lines. A panel's own legend
# has the room, and the value, so it keeps the specific wording.
_SHARED = [False]


def compact():
    return _COMPACT[0]


def set_mode(compact=False, shared=False):
    _COMPACT[0], _SHARED[0] = compact, shared


def lbl(full, short):
    """Legend label: the panel's own wording, or the wording a shared legend
    needs."""
    return short if _SHARED[0] else full


def sig_tex(v):
    """1e-2 -> 10^{-2}, 2.7e-3 -> 2.7\\times10^{-3}."""
    e = int(math.floor(math.log10(v)))
    m = v / 10.0 ** e
    return f"10^{{{e}}}" if abs(m - 1.0) < 1e-9 else f"{m:.1f}\\times10^{{{e}}}"


def mark_offscale(a, x, y, color, fmt="{:.1e}", row=0):
    """Flag where a curve leaves the frame, and say how far off scale it goes.

    Clipping to the usable window is the point of the tight limits, but a
    curve that simply vanishes at the frame reads as missing data rather
    than as a number too large to run. One marker on the edge plus the
    extreme value fixes that without restoring sixty decades of axis.
    """
    lo, hi = a.get_ylim()
    for out, edge, mk, va, extreme, arrow in (
            (y > hi, hi, "^", "top", np.nanmax, "↑"),
            (y < lo, lo, "v", "bottom", np.nanmin, "↓")):
        if not out.any() or out.all():
            continue
        # Mark the run that actually contains the extreme, at the end of that
        # run touching the frame -- a curve can leave and re-enter (the sigma0
        # kink does exactly that), and a marker parked in the middle of some
        # other excursion points at the wrong place.
        k = int(np.nanargmax(y) if va == "top" else np.nanargmin(y))
        lo_i = hi_i = k
        while lo_i > 0 and out[lo_i - 1]:
            lo_i -= 1
        while hi_i < len(out) - 1 and out[hi_i + 1]:
            hi_i += 1
        j = hi_i if hi_i < len(out) - 1 else (lo_i if lo_i > 0 else k)
        # One row lower in the master: the region names sit on the same edge
        # and there is a third of the width to keep them apart in.
        pad = 8.0 + 13.0 * (row + (1 if compact() else 0))
        a.plot([x[j]], [edge], mk, ms=9, color=color, zorder=7,
               clip_on=compact())
        a.annotate(f"{arrow} {fmt.format(extreme(y))}", (x[j], edge),
                   xytext=(0, -pad if va == "top" else pad),
                   textcoords="offset points",
                   fontsize=FS_C_TICK if compact() else FS_NOTE, color=color,
                   ha="center", va=va, zorder=7,
                   bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none",
                             alpha=0.85))



def style(a, xlabel, ylabel, title, logy=True, title_size=FS_TITLE, pad=10,
           short=None, short_y=None):
    lab = FS_C_LABEL if compact() else FS_LABEL
    a.set_xlabel(xlabel, fontsize=lab)
    a.set_ylabel(short_y if (compact() and short_y) else ylabel, fontsize=lab)
    if compact():
        title, title_size, pad = (short or title), FS_C_TITLE, 6
    if title:
        a.set_title(title, fontsize=title_size, color=INK, loc="left", pad=pad)
    if logy:
        a.set_yscale("log")
    a.grid(alpha=0.25, lw=0.5, color=GRID, which="both")
    a.spines[["top", "right"]].set_visible(False)
    a.tick_params(labelsize=FS_C_TICK if compact() else FS_TICK)



def blank():
    """An invisible handle, so a group heading can occupy a legend row."""
    return Line2D([], [], linestyle="none", marker="none")


def entries(*axes):
    """label -> handle, unioned across axes. Panels that draw the same line
    with the same label collapse to one entry, which is what lets the master
    grid carry a single legend for all six."""
    have = {}
    for a in axes:
        for h, lbl in zip(*a.get_legend_handles_labels()):
            have.setdefault(lbl, h)
    return have


def rows_stacked(have, groups):
    """Rows a single-column stacked legend will occupy: a heading plus its
    entries, per group that actually has something plotted."""
    n = 0
    for _t, keys in groups:
        keys = [k for k in keys if k in have]
        if keys:
            n += 1 + len(keys)
    return n


def grouped(target, have, groups, size, ncol=None, **kw):
    """A legend laid out as titled columns.

    A flat strip of entries is a lookup problem: the reader scans the whole row
    to find one curve. matplotlib fills a multi-column legend column-major, so
    padding every group to the same length puts each group in its own column
    under its own bold heading -- and the headings carry the provenance the
    short entry labels leave out.

    `groups` is [(heading, [label, ...]), ...]; labels that were never plotted
    are skipped, so one group list can serve several panels.
    """
    cols = [(t, [k for k in keys if k in have]) for t, keys in groups]
    cols = [c for c in cols if c[1]]
    # ncol=1 stacks the groups one after another instead of side by side --
    # what a 3 in wide legend needs, and padding there would only add blanks.
    stack = ncol == 1
    rows = 0 if stack else max(len(keys) for _t, keys in cols)
    handles, labels, headings = [], [], set()
    for title, keys in cols:
        handles.append(blank()); labels.append(title); headings.add(title)
        for k in keys:
            handles.append(have[k]); labels.append(k)
        if not stack:
            handles += [blank()] * (rows - len(keys))
            labels += [""] * (rows - len(keys))
    leg = target.legend(handles, labels, ncol=1 if stack else len(cols),
                        fontsize=size,
                        frameon=False, handlelength=3.0 if stack else 2.0,
                        handletextpad=0.6,
                        columnspacing=1.2, labelspacing=0.45,
                        borderaxespad=0.0, alignment="left", **kw)
    for t in leg.get_texts():
        if t.get_text() in headings:
            t.set_fontweight("bold")
    return leg



def legend(a, groups, y=-0.165):
    """Per-panel legend under the axes.

    A compact panel has no room for one, so it records what its legend WOULD
    have said on the axes instead. The caller reads that back to build the
    panel's own legend -- inline underneath it, or as a separate image -- so
    the two can never list different things than the panel actually drew.
    """
    a._lk_groups = groups
    if compact():
        return None
    return grouped(a, entries(a), groups, FS_LEG, loc="upper center",
                    bbox_to_anchor=(0.5, y))



def shade_binding(a, T, binding, y_frac=0.975):
    """Shade contiguous T-runs sharing a binding bound; name each region once.

    y_frac puts the region names on whichever edge the off-scale markers are
    not using -- they are on the bottom in panel 5 and the top in panel 6, and
    at 3 in wide there is no room to share an edge.
    """
    y, va = y_frac, ("top" if y_frac > 0.5 else "bottom")
    seen = set()
    i = 0
    while i < len(T):
        j = i
        while j + 1 < len(T) and binding[j + 1] == binding[i]:
            j += 1
        name = binding[i]
        a.axvspan(T[i], T[j], color=BOUND_COLOR.get(name, MUTED), alpha=0.12,
                  lw=0, zorder=0,
                  label=BOUND_LABEL[name] if name not in seen else None)
        seen.add(name)
        # Name the region inside it, so the shading is never colour-alone.
        a.annotate(name, (0.5 * (T[i] + T[j]), y),
                   xycoords=("data", "axes fraction"),
                   fontsize=FS_C_TICK if compact() else FS_NOTE + 1,
                   color=BOUND_COLOR.get(name, MUTED), ha="center", va=va,
                   fontweight="bold", zorder=7,
                   bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none",
                             alpha=0.85))
        i = j + 1




# =========================================================================
# The driver: panels, their legends, the master grid, the legend strip.
# =========================================================================

def leg_in(rows, size=FS_C_TICK):
    """Inches a stacked legend of `rows` rows needs. A row is about 1.45 line
    heights; the constant is the box padding."""
    return rows * size * 1.45 / 72.0 + 0.18


def save(fig, out, name, dpi, max_w=SLIDE_W_IN):
    """Write the figure, asserting the width ceiling.

    The ceiling is the whole point of the layout, so it is checked rather than
    left to whoever last edited a figsize.
    """
    w, h = fig.get_size_inches()
    assert w <= max_w + 1e-6, f"{name} is {w:.2f} in wide, > {max_w}"
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / f"{name}.png", dpi=dpi)
    plt.close(fig)
    print(f"plot -> {out / (name + '.png')}  ({w:.2f} x {h:.2f} in @ {dpi} dpi)")


def emit(out, panels, T, ctx, *, width, panel_width, dpi, panel_legend,
         master_groups, suptitle, ncols=3):
    """Write every figure for one family of panels.

    `panels` is [(name, fn(ax, T, ctx)), ...]. Each fn draws itself and calls
    legend() with its groups; legend() records them on the axes even when
    compact mode suppresses drawing, which is what lets this build a panel's
    own key without ever being able to list something the panel did not draw.
    """
    small = panel_width < COMPACT_BELOW_IN
    figsize = (panel_width, panel_width * (COMPACT_ASPECT if small else ASPECT))

    set_mode(compact=small)
    try:
        for name, fn in panels:
            if not (small and panel_legend == "inline"):
                fig, a = plt.subplots(figsize=figsize, layout="constrained")
                fn(a, T, ctx)
                save(fig, out, name, dpi, width)
            else:
                # Draw once on a throwaway axes purely to learn how many legend
                # rows this panel needs, so the figure can be sized for them
                # rather than stealing the height from the plot.
                probe, pa = plt.subplots(figsize=figsize)
                fn(pa, T, ctx)
                rows = rows_stacked(entries(pa), getattr(pa, "_lk_groups", []))
                plt.close(probe)
                lh = leg_in(rows)
                fig = plt.figure(figsize=(figsize[0], figsize[1] + lh),
                                 layout="constrained")
                gs = fig.add_gridspec(2, 1, height_ratios=[figsize[1], lh])
                a, al = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
                al.axis("off")
                fn(a, T, ctx)
                grouped(al, entries(a), getattr(a, "_lk_groups", []),
                        FS_C_TICK, ncol=1, loc="center")
                save(fig, out, name, dpi, width)

            if small and panel_legend == "separate":
                have, groups = entries(a), getattr(a, "_lk_groups", [])
                figL = plt.figure(
                    figsize=(figsize[0], leg_in(rows_stacked(have, groups))),
                    layout="constrained")
                grouped(figL, have, groups, FS_C_TICK, ncol=1, loc="center")
                save(figL, out, f"{name}_legend", dpi, width)
    finally:
        set_mode()

    # --- the master: every panel inside the same width ---------------------
    nrows = -(-len(panels) // ncols)
    set_mode(compact=True, shared=True)
    try:
        fig, axes = plt.subplots(nrows, ncols, layout="constrained",
                                 figsize=(width, width * 0.33 * nrows))
        for ax, (_n, fn) in zip(np.atleast_1d(axes).ravel(), panels):
            fn(ax, T, ctx)
        for ax in np.atleast_1d(axes).ravel()[len(panels):]:
            ax.axis("off")
        fig.suptitle(suptitle, fontsize=FS_C_LABEL, color=INK)
        have = entries(*np.atleast_1d(axes).ravel()[:len(panels)])
        grouped(fig, have, master_groups(ctx), FS_C_TICK,
                loc="outside lower center")
        save(fig, out, "fig0_master", dpi, width)

        # The same legend on its own: the small panels carry none in 'none'
        # mode, and a slide may lay several out under one key.
        figL = plt.figure(figsize=(width, 1.15), layout="constrained")
        grouped(figL, have, master_groups(ctx), FS_C_TICK, loc="center")
        save(figL, out, "fig7_legend", dpi, width)
    finally:
        set_mode()
