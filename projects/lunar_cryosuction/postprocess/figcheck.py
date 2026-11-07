"""Deterministic layout check for a matplotlib figure: overlaps and font size.

The user's standing rule: no figure is shown with overlapping elements, and no
font is smaller than 7 pt at the size the figure is finally shown. This module
measures both after drawing, so a figure script can refuse to pass a bad
layout instead of relying on someone looking.

    from figcheck import check_figure
    problems = check_figure(fig, min_pt=7.0, placed_width_in=None)
    # [] when clean; otherwise one human-readable string per problem

What is checked:
  1. every pair of visible text artists in the figure (titles, axis labels,
     tick labels, legend entries, annotations, colorbar labels): no overlap;
  2. axes (including colorbars) against each other: no overlap;
  3. free text inside an axes (ax.texts) and the legend box against the data
     drawn in that axes (lines, patches, images, collections): no overlap;
  4. the smallest font, scaled by placed_width_in / figure width when the
     figure is placed at another width: at least min_pt.

Mark an intentional overlap with artist.set_gid("overlap_ok") (for example a
label written inside its own shape). Tick labels are compared with other
texts but not with data, because they live outside the axes.
"""
from __future__ import annotations

import matplotlib
import matplotlib.collections
from matplotlib.text import Text


def _tick_labels(fig):
    """(all tick-label Texts, the ones actually drawn). matplotlib keeps tick
    objects for ticks outside the view limits and a hidden label2 on each; only
    the ticks returned by Axis._update_ticks() are drawn."""
    every, drawn = set(), set()
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            for tk in axis.majorTicks + axis.minorTicks:
                every.update((id(tk.label1), id(tk.label2)))
            for tk in axis._update_ticks():
                for lab in (tk.label1, tk.label2):
                    if lab.get_visible():
                        drawn.add(id(lab))
    return every, drawn


def _visible_texts(fig):
    every, drawn = _tick_labels(fig)
    out = []
    for t in fig.findobj(Text):
        if not t.get_visible() or not t.get_text().strip():
            continue
        if t.get_gid() == "overlap_ok":
            continue
        if id(t) in every and id(t) not in drawn:
            continue
        out.append(t)
    return out


def _tbox(t, rnd):
    bp = t.get_bbox_patch()
    bb = bp.get_window_extent(rnd) if bp is not None else Text.get_window_extent(t, rnd)
    return bb


def _shrink(bb, px):
    from matplotlib.transforms import Bbox
    return Bbox.from_extents(bb.x0 + px, bb.y0 + px, bb.x1 - px, bb.y1 - px)


def _data_paths(ax):
    paths = []
    for ln in ax.lines:
        if ln.get_visible() and ln.get_gid() != "overlap_ok":
            paths.append((f"line {ln.get_label()}", ln.get_path().transformed(ln.get_transform()), False))
    for pa in ax.patches:
        if pa.get_visible() and pa.get_gid() != "overlap_ok":
            paths.append((f"patch {type(pa).__name__}",
                          pa.get_path().transformed(pa.get_transform()), True))
    return paths


def check_figure(fig, min_pt: float = 7.0, placed_width_in: float | None = None,
                 pad_px: float = 0.5) -> list[str]:
    fig.canvas.draw()
    rnd = fig.canvas.get_renderer()
    bad: list[str] = []

    texts = _visible_texts(fig)
    boxes = [(t, _shrink(_tbox(t, rnd), pad_px)) for t in texts]
    for i, (t1, b1) in enumerate(boxes):
        for t2, b2 in boxes[i + 1:]:
            if b1.overlaps(b2):
                bad.append(f"text/text {t1.get_text()[:30]!r} x {t2.get_text()[:30]!r}")

    axes = [ax for ax in fig.axes if ax.get_visible()]
    for i, a1 in enumerate(axes):
        for a2 in axes[i + 1:]:
            if _shrink(a1.get_window_extent(rnd), pad_px).overlaps(_shrink(a2.get_window_extent(rnd), pad_px)):
                bad.append(f"axes/axes {a1.get_label() or a1.get_title()!r} x {a2.get_label() or a2.get_title()!r}")

    for ax in axes:
        free = [(t, _shrink(_tbox(t, rnd), pad_px)) for t in ax.texts
                if t.get_visible() and t.get_text().strip() and t.get_gid() != "overlap_ok"]
        leg = ax.get_legend()
        if leg is not None and leg.get_visible():
            free.append((leg, _shrink(leg.get_window_extent(rnd), pad_px)))
        paths = _data_paths(ax)
        for art, bb in free:
            name = art.get_text()[:30] if isinstance(art, Text) else "legend"
            for pname, pth, filled in paths:
                if pth.intersects_bbox(bb, filled=filled):
                    bad.append(f"{name!r} x {pname}")
                    break
            for im in ax.images:
                if im.get_visible() and im.get_window_extent(rnd).overlaps(bb):
                    bad.append(f"{name!r} x image")
            for col in ax.collections:
                if not col.get_visible() or col.get_gid() == "overlap_ok":
                    continue
                try:
                    if col.get_window_extent(rnd).overlaps(bb) and isinstance(
                            col, matplotlib.collections.QuadMesh):
                        bad.append(f"{name!r} x {type(col).__name__}")
                except Exception:
                    continue

    if texts:
        scale = 1.0 if placed_width_in is None else placed_width_in / fig.get_size_inches()[0]
        small = min(t.get_fontsize() for t in texts) * scale
        if small < min_pt:
            bad.append(f"smallest printed font {small:.1f} pt < {min_pt} pt (scale {scale:.2f})")
    return bad
