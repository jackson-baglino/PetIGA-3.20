#!/usr/bin/env python3
"""Manuscript figure panels for an axisymmetric grain-pair sintering run.

Produces the pieces of a multi-panel figure as SEPARATE vector files, so the
figure itself is assembled in Inkscape:

    xsec_<tag>.pdf/.svg          full-domain meridional cross-section:
                                 vapour density in the air, opaque ice body
    ice3d_<tag>.pdf/.svg         shaded 3-D rendering of the ice body (the
                                 phi = 0.5 surface revolved about the
                                 symmetry axis), transparent background
    ice3d_vapour_<tag>.pdf/.svg  the same 3-D body above the symmetry axis,
                                 the vapour cross-section below it

    cbar_rhov_{light,dark}       vapour bar, absolute rho_v
    cbar_supersat_{light,dark}   the same bar labelled as supersaturation
    cbar_ice_{light,dark}        ice-phase (phi) bar
    scalebar_{light,dark}        standalone scale bar (see below)

EVERY PANEL IS DRAWN AT ONE COMMON SCALE (--um-per-inch), so panels drop into
Inkscape side by side without rescaling and a single scale bar serves all of
them. That is also why the scale bar ships as its own file rather than being
burnt into each panel: the 3-D panels have transparent backgrounds and the
right ink for them depends on the page.

WHY EVERY PANEL SHARES ONE PROJECTION. The 3-D renderer views the body along
+z, exactly perpendicular to the symmetry axis. For a surface of revolution
that choice has two properties nothing else has: the silhouette IS the
phi = 0.5 profile, and the surface cannot occlude itself (at fixed x the
section is a single circle and the view direction has no x-component, so
every sight line meets the surface exactly twice, once front, once back).
The first is what lets the 3-D body and the cross-section register exactly
when stacked across the axis; the second is what makes the render correct
without a z-buffer, which is what lets it stay vector.

WHY THE SHADING IS CONTOURED, NOT TESSELLATED. Tessellating the revolved
surface into quads would put ~10^5 polygons into the PDF and make the file
unusable in Inkscape. Instead the Lambert+Blinn intensity is evaluated on a
dense (profile x theta) grid and drawn as filled BANDS of that intensity in
projected coordinates. Visually that is a finely posterised render;
structurally it is a few dozen paths.

WHY THE BANDS ARE NESTED REGIONS, NOT contourf. Two reasons, both learned the
hard way. (1) Adjacent contourf polygons share an edge, and in vector output
that seam shows as a hairline; the usual cure, set_edgecolor("face"), strokes
contourf's *compound* path and draws a visible line across the figure joining
disjoint pieces of the same band. (2) contourf keeps every vertex, which put
11 MB into a single SVG. Here band i is drawn as the whole region {Z >= l_i},
Douglas-Peucker simplified, painted lightest-first: the regions are nested by
construction, so each one covers the last and there is no seam to hide and
nothing to stroke.

Usage
-----
  python postprocess/manuscript_grainpair_figure.py <run_dir> \
      --steps 62 236 --tags t0 t78
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.collections import LineCollection, PathCollection
from matplotlib.path import Path
from contourpy import contour_generator, FillType

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pplib import read_opts, opt_float, rho_vs, step_times  # noqa: E402

UM = 1e-6
BIG = 1e30

# Text stays live text in the SVG (Inkscape can edit it) and is embedded as a
# TrueType subset in the PDF rather than converted to outlines.
matplotlib.rcParams.update({
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
})

# ---------------------------------------------------------------------------
# palette
# ---------------------------------------------------------------------------
# Sequential = one hue, light -> dark. Blue ramp steps 100..700; low rho_v (the
# undersaturated chamber wall) is the light end, high rho_v (the vapour collar
# hugging the ice) the dark end.
BLUE_RAMP = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
             "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281",
             "#0d366b"]
CMAP_RHOV = LinearSegmentedColormap.from_list("rhov_blue", BLUE_RAMP)

# Ice ramp: phi = 0 (air) at the dark end, phi = 1 (ice) at the near-white end.
# The 3-D shading reads its tones off this same ramp, so the phi colourbar and
# the rendered body belong to one visual family.
ICE_RAMP = ["#33465a", "#546a80", "#7b90a4", "#a8bccd", "#d5e3ee", "#f7fbff"]
CMAP_ICE = LinearSegmentedColormap.from_list("ice", ICE_RAMP)

ICE_FILL = "#fbfdff"     # flat ice colour in the 2-D cross-section
ICE_EDGE = "#16324f"     # phi = 0.5 stroke
AXIS_INK = "#16324f"

INK_LIGHT = "#0b0b0b"    # text/rules on a light background
INK_DARK = "#f4f3ee"     # text/rules on a dark background

N_BAND_RHOV = 40         # vapour bands
N_BAND_SHADE = 48        # shading bands on the 3-D body


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------
def sample(run_dir, step, nu=1801, nv=901):
    """Evaluate (phi, rho_v) of one snapshot on a uniform grid over the whole
    meridional half-domain.

    Returns (x, y, phi, rhov); x and y are 1-D in metres, the fields are
    (ny, nx) -- nrb() hands back (nu, nv, ndof), which is the transpose of
    what contouring wants.
    """
    from igakit.io import PetIGA

    nrb = PetIGA().read(os.path.join(run_dir, "igasol.dat"))
    opts = read_opts(run_dir)
    u = np.linspace(0.0, opt_float(opts, "-Lx"), nu)
    v = np.linspace(0.0, opt_float(opts, "-Ly"), nv)
    sol = PetIGA().read_vec(os.path.join(run_dir, f"sol_{step:05d}.dat"), nrb)
    _, F = nrb(u, v, fields=sol)
    return u, v, F[..., 0].T, F[..., 2].T


def ice_profile(x, y, phi, simplify_um=0.02):
    """The phi = 0.5 meridional profile, left to right, as (M, 2) in microns.

    The ice sits on the symmetry axis, so the level set is an OPEN curve with
    both ends on y = 0. Only the longest line is kept: at these eps the level
    set is a single curve and anything else is a stray fragment.
    """
    gen = contour_generator(x=x / UM, y=y / UM, z=phi, name="serial")
    segs = [np.asarray(s) for s in gen.lines(0.5) if len(np.asarray(s)) > 2]
    if not segs:
        raise RuntimeError("no phi = 0.5 contour found")
    p = max(segs, key=lambda s: np.hypot(*np.diff(s, axis=0).T).sum())
    if p[0, 0] > p[-1, 0]:
        p = p[::-1]
    p = p.copy()
    p[:, 1] = np.maximum(p[:, 1], 0.0)
    return rdp(p, simplify_um)


def rdp(pts, tol):
    """Ramer-Douglas-Peucker, iterative (a 10k-vertex ring would blow the
    recursion limit). This is what keeps the vector files small enough to
    edit."""
    if tol <= 0 or len(pts) < 4:
        return pts
    keep = np.zeros(len(pts), dtype=bool)
    keep[0] = keep[-1] = True
    stack = [(0, len(pts) - 1)]
    while stack:
        i, j = stack.pop()
        if j <= i + 1:
            continue
        seg = pts[j] - pts[i]
        norm = float(np.hypot(*seg))
        rel = pts[i + 1:j] - pts[i]
        if norm == 0.0:
            d = np.hypot(rel[:, 0], rel[:, 1])
        else:
            d = np.abs(rel[:, 0] * seg[1] - rel[:, 1] * seg[0]) / norm
        k = int(np.argmax(d))
        if d[k] > tol:
            k += i + 1
            keep[k] = True
            stack += [(i, k), (k, j)]
    return pts[keep]


# ---------------------------------------------------------------------------
# filled bands as nested regions
# ---------------------------------------------------------------------------
def filled_bands(ax, X, Y, Z, levels, colors, tol, zorder=1, label=None):
    """Paint {Z >= levels[i]} in colors[i], lightest level first.

    The regions nest, so painting in order builds the ramp with no shared
    edges: no seams, no strokes, and one PathCollection per band (which is
    also one tidy group per band in Inkscape).

    X, Y may be curvilinear -- that is what lets the same routine draw the
    projected 3-D shading and the flat cross-section.
    """
    gen = contour_generator(x=X, y=Y, z=Z, name="serial",
                            fill_type=FillType.OuterOffset)
    artists = []
    for i, (lv, col) in enumerate(zip(levels, colors)):
        pts_list, off_list = gen.filled(lv, BIG)
        paths = []
        for pts, offs in zip(pts_list, off_list):
            verts, codes = [], []
            for a, b in zip(offs[:-1], offs[1:]):
                ring = rdp(np.asarray(pts[a:b]), tol)
                if len(ring) < 3:
                    continue
                verts.append(ring)
                codes.append(np.r_[Path.MOVETO,
                                   np.full(len(ring) - 2, Path.LINETO),
                                   Path.CLOSEPOLY])
            if verts:
                paths.append(Path(np.vstack(verts), np.concatenate(codes)))
        if not paths:
            continue
        pc = PathCollection(paths, facecolors=[col], edgecolors="none",
                            linewidths=0.0, antialiaseds=True,
                            zorder=zorder + i * 1e-3)
        if label:
            pc.set_label(f"{label}_{i:02d}")
        ax.add_collection(pc)
        artists.append(pc)
    return artists


def band_levels(vmin, vmax, n, cmap):
    """n band thresholds spanning [vmin, vmax] and the colour for each.

    The first threshold sits slightly below vmin so the lightest band covers
    every unmasked cell -- otherwise the region at exactly vmin is empty and
    the figure shows holes where the field is at its floor.
    """
    lv = np.linspace(vmin, vmax, n + 1)[:-1]
    lv[0] -= 1e-9 * max(abs(vmax - vmin), 1.0)
    return lv, cmap(np.linspace(0.0, 1.0, n))


# ---------------------------------------------------------------------------
# 3-D rendering of the revolved ice body
# ---------------------------------------------------------------------------
def _upper_hull(p):
    """Upper convex hull of the profile (Andrew monotone chain), returned as
    r_hull(x) sampled at the profile's own x. Used only as an ambient-occlusion
    proxy: where the profile dips below its own hull it is inside a crevice."""
    def turn(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    h = []
    for q in p:
        while len(h) >= 2 and turn(h[-2], h[-1], q) >= 0:
            h.pop()
        h.append(q)
    h = np.array(h)
    return np.interp(p[:, 0], h[:, 0], h[:, 1])


def shade_revolved(ax, profile, n_theta=721, light=(-0.42, 0.66, 0.62),
                   ambient=0.24, diffuse=0.58, specular=0.14, shine=30.0,
                   ao_min=0.44, ao_pow=0.70, n_bands=N_BAND_SHADE,
                   tol=0.12, zorder=3):
    """Draw the ice body: the profile revolved about y = 0, Lambert + Blinn
    shaded with a crevice-darkening term, as filled bands of the intensity.

    The camera looks along -z orthographically, so screen (X, Y) = (x, y) in
    micrometres and the silhouette is the profile itself.

    AMBIENT OCCLUSION IS NOT COSMETIC HERE. Without it the neck renders
    *brighter* than the grains: its fillet sweeps the normal through the
    mirror direction, so the specular term spikes right where the surface is
    most enclosed and should be darkest. The proxy used is the profile's
    depth below its own upper convex hull, which is zero on the grains and
    large exactly in the neck.
    """
    xp, rp = profile[:, 0], profile[:, 1]
    dx, dr = np.gradient(xp), np.gradient(rp)

    th = np.linspace(0.0, 2.0 * np.pi, n_theta)
    c, s = np.cos(th)[None, :], np.sin(th)[None, :]

    # Outward normal of a surface of revolution about x:
    #   n  proportional to  (-r', x' cos(theta), x' sin(theta))
    nx = np.broadcast_to(-dr[:, None], (len(xp), n_theta))
    ny = dx[:, None] * c
    nz = dx[:, None] * s
    nn = np.sqrt(nx ** 2 + ny ** 2 + nz ** 2)
    nn[nn == 0.0] = 1.0
    nx, ny, nz = nx / nn, ny / nn, nz / nn

    view = np.array([0.0, 0.0, 1.0])
    L = np.asarray(light, float)
    L /= np.linalg.norm(L)
    H = L + view
    H /= np.linalg.norm(H)

    hull = _upper_hull(np.column_stack([xp, rp]))
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = np.clip(np.where(hull > 0, rp / hull, 1.0), 0.0, 1.0)
    ao = (ao_min + (1.0 - ao_min) * frac ** ao_pow)[:, None]

    lam = np.clip(nx * L[0] + ny * L[1] + nz * L[2], 0.0, None)
    spc = np.clip(nx * H[0] + ny * H[1] + nz * H[2], 0.0, None) ** shine
    inten = np.clip(ao * (ambient + diffuse * lam + specular * spc), 0.0, 1.0)

    # Screen coordinates, back half masked away. At this camera the surface
    # cannot occlude itself, so back-face culling is the whole of the
    # hidden-surface problem.
    X = np.broadcast_to(xp[:, None], (len(xp), n_theta))
    Y = rp[:, None] * c
    Z = np.ma.masked_where(nz <= 1e-9, inten)

    lv, cols = band_levels(0.0, 1.0, n_bands, CMAP_ICE)
    art = filled_bands(ax, X, Y, Z, lv, cols, tol, zorder=zorder,
                       label="shade")

    # Exact silhouette, on top of the banded fill.
    sil = np.vstack([np.column_stack([xp, rp]),
                     np.column_stack([xp[::-1], -rp[::-1]])])
    art.append(ax.plot(sil[:, 0], sil[:, 1], color=ICE_EDGE, lw=0.7,
                       solid_joinstyle="round", zorder=zorder + 1)[0])
    return art


# ---------------------------------------------------------------------------
# 2-D vapour cross-section
# ---------------------------------------------------------------------------
def vapour_section(ax, x, y, phi, rhov, vmin, vmax, stride=2, tol=0.25,
                   n_bands=N_BAND_RHOV, zorder=1):
    """Filled bands of rho_v in the air phase, mirrored about the axis.

    Masking at phi > 0.5 leaves the band edges ragged at grid resolution; the
    opaque ice polygon drawn on top is what makes the boundary crisp, so the
    mask only has to be approximately right -- which is why the field can be
    decimated (`stride`) without any visible cost.
    """
    xu, yu = x[::stride] / UM, y[::stride] / UM
    ph, rh = phi[::stride, ::stride], rhov[::stride, ::stride]

    yfull = np.concatenate([-yu[:0:-1], yu])
    ph = np.vstack([ph[:0:-1, :], ph])
    rh = np.vstack([rh[:0:-1, :], rh])

    X, Y = np.meshgrid(xu, yfull)
    Z = np.ma.masked_where(ph > 0.5, rh)
    lv, cols = band_levels(vmin, vmax, n_bands, CMAP_RHOV)
    return filled_bands(ax, X, Y, Z, lv, cols, tol, zorder=zorder,
                        label="rhov")


def draw_ice_body(ax, profile, fill=ICE_FILL, lw=0.7, zorder=4):
    """Opaque ice cross-section from the phi = 0.5 profile, mirrored."""
    xp, rp = profile[:, 0], profile[:, 1]
    poly = np.vstack([np.column_stack([xp, rp]),
                      np.column_stack([xp[::-1], -rp[::-1]])])
    return ax.fill(poly[:, 0], poly[:, 1], facecolor=fill, edgecolor=ICE_EDGE,
                   lw=lw, joinstyle="round", zorder=zorder)


def draw_axis_line(ax, x0, x1, ink=AXIS_INK, lw=0.5, zorder=6):
    lc = LineCollection([np.array([[x0, 0.0], [x1, 0.0]])], colors=ink,
                        linewidths=lw, linestyles=[(0, (6, 3))], zorder=zorder)
    ax.add_collection(lc)
    return [lc]


# ---------------------------------------------------------------------------
# figure plumbing
# ---------------------------------------------------------------------------
def _blank_axes(fig, xlim, ylim):
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal")
    ax.set_autoscale_on(False)
    ax.axis("off")
    return ax


def _panel_figure(xlim, ylim, um_per_inch):
    """A figure whose data window maps to inches at the common scale."""
    w = (xlim[1] - xlim[0]) / um_per_inch
    h = (ylim[1] - ylim[0]) / um_per_inch
    fig = plt.figure(figsize=(w, h))
    return fig, _blank_axes(fig, xlim, ylim)


def _clip(artists, ax, x0, y0, w, h):
    rect = plt.Rectangle((x0, y0), w, h, transform=ax.transData)
    for a in artists:
        a.set_clip_path(rect)


def _save(fig, out, stem, formats=("pdf", "svg")):
    paths = []
    for fmt in formats:
        p = os.path.join(out, f"{stem}.{fmt}")
        fig.savefig(p, format=fmt, transparent=True, pad_inches=0)
        paths.append(p)
    return paths


# ---------------------------------------------------------------------------
# panels
# ---------------------------------------------------------------------------
def panel_xsec(out, tag, x, y, phi, rhov, profile, vmin, vmax, ums, stride):
    """Full-domain meridional cross-section, mirrored about the axis."""
    Lx, Ly = x[-1] / UM, y[-1] / UM
    fig, ax = _panel_figure((0, Lx), (-Ly, Ly), ums)
    vapour_section(ax, x, y, phi, rhov, vmin, vmax, stride=stride)
    draw_ice_body(ax, profile)
    draw_axis_line(ax, 0, Lx)
    return _save(fig, out, f"xsec_{tag}"), fig


def panel_ice3d(out, tag, profile, crop, ums):
    """The revolved ice body on its own, transparent background."""
    (x0, x1), ymax = crop
    fig, ax = _panel_figure((x0, x1), (-ymax, ymax), ums)
    shade_revolved(ax, profile)
    return _save(fig, out, f"ice3d_{tag}"), fig


def panel_ice3d_vapour(out, tag, x, y, phi, rhov, profile, vmin, vmax, ums,
                       stride):
    """3-D body above the symmetry axis, meridional section below it, in the
    same frame as panel_xsec.

    The vapour field is painted across the WHOLE frame and the 3-D body sits
    opaquely in the upper half of it. That is legitimate rather than a cheat:
    the plane being coloured is the meridional plane through the axis, the
    same plane the section below is cut from, and the body's silhouette is
    exactly the phi = 0.5 profile -- so the ice covers precisely the part of
    that plane the ice occupies, and the two halves register with no seam.
    """
    Lx, Ly = x[-1] / UM, y[-1] / UM
    fig, ax = _panel_figure((0, Lx), (-Ly, Ly), ums)

    vapour_section(ax, x, y, phi, rhov, vmin, vmax, stride=stride)
    _clip(draw_ice_body(ax, profile), ax, 0.0, -Ly, Lx, Ly)
    _clip(shade_revolved(ax, profile), ax, 0.0, 0.0, Lx, Ly)
    draw_axis_line(ax, 0, Lx, lw=0.6)
    return _save(fig, out, f"ice3d_vapour_{tag}"), fig


# ---------------------------------------------------------------------------
# colourbars and scale bars
# ---------------------------------------------------------------------------
def colourbar(out, stem, cmap, vmin, vmax, ticks, ticklabels, label, ink,
              n_bands=256, horizontal=True):
    """A standalone vector colourbar: band rectangles, so the file carries
    paths rather than an embedded raster."""
    if horizontal:
        fig = plt.figure(figsize=(2.6, 0.64))
        ax = fig.add_axes([0.04, 0.45, 0.92, 0.33])
    else:
        fig = plt.figure(figsize=(1.0, 2.7))
        ax = fig.add_axes([0.08, 0.05, 0.28, 0.90])

    edges = np.linspace(vmin, vmax, n_bands + 1)
    cols = cmap(np.linspace(0.0, 1.0, n_bands))
    for i in range(n_bands):
        lo, hi = edges[i], edges[i + 1]
        xy, w, h = ((lo, 0.0), hi - lo, 1.0) if horizontal \
            else ((0.0, lo), 1.0, hi - lo)
        # A hairline stroke in the band's own colour closes the sub-pixel gaps
        # a pure fill would leave between 256 adjacent rectangles.
        ax.add_patch(plt.Rectangle(xy, w, h, facecolor=cols[i],
                                   edgecolor=cols[i], lw=0.4))

    if horizontal:
        ax.set_xlim(vmin, vmax)
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.set_xticks(ticks)
        ax.set_xticklabels(ticklabels)
        ax.tick_params(axis="x", colors=ink, labelsize=7, length=2.5,
                       width=0.6, pad=2)
        ax.set_xlabel(label, color=ink, fontsize=7.5, labelpad=2)
    else:
        ax.set_ylim(vmin, vmax)
        ax.set_xlim(0, 1)
        ax.set_xticks([])
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
        ax.set_yticks(ticks)
        ax.set_yticklabels(ticklabels)
        ax.tick_params(axis="y", colors=ink, labelsize=7, length=2.5,
                       width=0.6, pad=2)
        ax.set_ylabel(label, color=ink, fontsize=7.5, labelpad=4)

    for sp in ax.spines.values():
        sp.set_color(ink)
        sp.set_linewidth(0.6)

    paths = _save(fig, out, stem)
    plt.close(fig)
    return paths


def write_colourbars(out, vmin, vmax, rvs, orientations=("h", "v")):
    made = []
    for orient in orientations:
        horiz = orient == "h"
        sfx = "" if horiz else "_vert"
        for mode, ink in (("light", INK_LIGHT), ("dark", INK_DARK)):
            t = np.array([8.465e-4, 8.470e-4, 8.475e-4, 8.480e-4])
            made += colourbar(
                out, f"cbar_rhov_{mode}{sfx}", CMAP_RHOV, vmin, vmax,
                t, [f"{v * 1e4:.3f}" for v in t],
                r"$\rho_v$  ($10^{-4}$ kg m$^{-3}$)", ink, horizontal=horiz)

            # The same bar, relabelled: supersaturation w.r.t. flat ice is the
            # quantity the growth law actually responds to.
            sig = np.array([-2.5, -2.0, -1.5, -1.0, -0.5])
            made += colourbar(
                out, f"cbar_supersat_{mode}{sfx}", CMAP_RHOV, vmin, vmax,
                rvs * (1.0 + sig * 1e-3), [f"{s:.1f}" for s in sig],
                r"$(\rho_v-\rho_{vs})/\rho_{vs}$  ($\times 10^{-3}$)", ink,
                horizontal=horiz)

            made += colourbar(
                out, f"cbar_ice_{mode}{sfx}", CMAP_ICE, 0.0, 1.0,
                [0.0, 0.5, 1.0], ["0", "0.5", "1"],
                r"Ice phase  $\phi$", ink, horizontal=horiz)
    return made


def write_scalebars(out, um_per_inch, length_um=100.0):
    """Scale bars at the panels' own scale, so one of these is correct for
    every panel as long as the panels are not rescaled in Inkscape."""
    made = []
    for mode, ink in (("light", INK_LIGHT), ("dark", INK_DARK)):
        w = length_um / um_per_inch
        fig = plt.figure(figsize=(w, 0.30))
        ax = fig.add_axes([0, 0, 1, 1])
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        ax.add_line(plt.Line2D([0, 1], [0.30, 0.30], color=ink, lw=1.8,
                               solid_capstyle="butt"))
        ax.text(0.5, 0.42, f"{length_um:g} " + "µm", ha="center",
                va="bottom", color=ink, fontsize=8)
        made += _save(fig, out, f"scalebar_{mode}")
        plt.close(fig)
    return made


README = """\
Manuscript figure pieces -- {run}

Generated by postprocess/manuscript_grainpair_figure.py on the snapshots
listed below. Everything is vector (.pdf and .svg, identical content); assemble
the figure in Inkscape.

SNAPSHOTS
{snaps}
The two times are the ones Molaro et al. (2019) measured at, transferred to the
model clock by their own overlay convention: t* is where the simulated neck
first reaches their first measured width (2 r = 32.81 um), and the second panel
is t* + 78 min, the span their last point covers. Anchoring on neck size rather
than on the solver's t = 0 is what makes the comparison fair -- the run starts
from a pre-necked r = 14 um pair, which is not the state the experiment opens
in.

SCALE
Every panel is drawn at {ums:g} um per inch, so panels drop in side by side with
no rescaling and one scale bar is correct for all of them. 1 um = {ptum:.4f} pt.
Do not rescale a panel in Inkscape unless you rescale all of them and the
scale bar together.

PANELS  (one pair of files per time tag)
  xsec_<tag>          full domain, meridional section mirrored about the axis:
                      rho_v in the air, opaque ice body, phi = 0.5 outline
  ice3d_<tag>         the ice body alone -- the phi = 0.5 surface revolved
                      about the symmetry axis, shaded, TRANSPARENT background
  ice3d_vapour_<tag>  the two combined in the xsec frame: 3-D body above the
                      axis, cut section below, vapour field behind both

KEYS  (light = dark ink for a light page; dark = light ink for a dark page;
       _vert = vertical bar)
  cbar_rhov_*         vapour density, 10^-4 kg m^-3
  cbar_supersat_*     the SAME bar relabelled as (rho_v - rho_vs)/rho_vs;
                      use one or the other, not both
  cbar_ice_*          ice phase phi, over the ice tone ramp the 3-D body and
                      the cut face are drawn from
  scalebar_*          {sblen:g} um at the panels' own scale

The rho_v colour scale is SHARED across every time shown, so panels are
directly comparable: [{vmin:.6e}, {vmax:.6e}] kg m^-3, i.e.
rho_v/rho_vs from {smin:.6f} to {smax:.6f} at T = {T:g} C.

preview_*.png are checking renders, not figure material.
"""


def write_readme(out, run, snaps, ums, vmin, vmax, rvs, T, sblen=100.0):
    txt = README.format(
        run=run, snaps="".join(f"  {s}\n" for s in snaps), ums=ums,
        ptum=72.0 / ums, sblen=sblen, vmin=vmin, vmax=vmax,
        smin=vmin / rvs, smax=vmax / rvs, T=T)
    p = os.path.join(out, "README.txt")
    with open(p, "w") as fh:
        fh.write(txt)
    return [p]


# ---------------------------------------------------------------------------
def _neck_radius(prof):
    """Minimum profile radius between the two grain apices, in micrometres."""
    r = prof[:, 1]
    half = len(r) // 2
    i = int(np.argmax(r[:half]))
    j = half + int(np.argmax(r[half:]))
    return float(r[i:j + 1].min()) if j > i else float(r.min())


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--steps", type=int, nargs="+", required=True)
    ap.add_argument("--tags", nargs="+", default=None,
                    help="one label per step (default: the step numbers)")
    ap.add_argument("--out", default=None,
                    help="default: <run_dir>/plots/manuscript")
    ap.add_argument("--crop", type=float, nargs=3,
                    metavar=("X0", "X1", "YMAX"), default=(35.0, 415.0, 118.0),
                    help="micrometre window for the 3-D panels")
    ap.add_argument("--um-per-inch", type=float, default=125.0,
                    help="common scale for every panel and the scale bars")
    ap.add_argument("--nu", type=int, default=1801)
    ap.add_argument("--nv", type=int, default=901)
    ap.add_argument("--field-stride", type=int, default=2,
                    help="decimation of the vapour field before banding")
    ap.add_argument("--no-preview", action="store_true")
    args = ap.parse_args()

    run = os.path.abspath(args.run_dir)
    out = args.out or os.path.join(run, "plots", "manuscript")
    os.makedirs(out, exist_ok=True)
    tags = args.tags or [f"step{s:05d}" for s in args.steps]
    if len(tags) != len(args.steps):
        sys.exit("--tags must have one entry per --steps")

    T = opt_float(read_opts(run), "-temp", -20.0)
    rvs = float(rho_vs(T))
    tmap = step_times(run)
    snaps = []

    # One pass first, so the colour scale is shared across times: panels at
    # different times are only comparable if it is.
    data = {}
    for step, tag in zip(args.steps, tags):
        x, y, phi, rhov = sample(run, step, args.nu, args.nv)
        prof = ice_profile(x, y, phi)
        data[tag] = (x, y, phi, rhov, prof)
        t = tmap.get(step, float("nan"))
        line = (f"{tag}: step {step}  t = {t:9.1f} s ({t / 60:6.2f} min)  "
                f"neck radius = {_neck_radius(prof):5.2f} um  "
                f"R_large = {prof[:, 1].max():6.2f} um")
        snaps.append(line)
        print("  " + line)

    vmin = min(float(np.min(d[3][d[2] <= 0.5])) for d in data.values())
    vmax = max(float(np.max(d[3][d[2] <= 0.5])) for d in data.values())
    print(f"  shared rho_v scale [{vmin:.6e}, {vmax:.6e}] kg/m^3 = "
          f"S [{vmin / rvs:.6f}, {vmax / rvs:.6f}]")

    crop = ((args.crop[0], args.crop[1]), args.crop[2])
    ums, stride = args.um_per_inch, args.field_stride
    written = []

    for tag, (x, y, phi, rhov, prof) in data.items():
        made = []
        p, f1 = panel_xsec(out, tag, x, y, phi, rhov, prof, vmin, vmax, ums,
                           stride)
        made.append((p, f1, f"xsec_{tag}"))
        p, f2 = panel_ice3d(out, tag, prof, crop, ums)
        made.append((p, f2, f"ice3d_{tag}"))
        p, f3 = panel_ice3d_vapour(out, tag, x, y, phi, rhov, prof, vmin,
                                   vmax, ums, stride)
        made.append((p, f3, f"ice3d_vapour_{tag}"))
        for p, fig, stem in made:
            written += p
            if not args.no_preview:
                fig.savefig(os.path.join(out, f"preview_{stem}.png"), dpi=220,
                            transparent=False, facecolor="#ffffff")
            plt.close(fig)

    written += write_colourbars(out, vmin, vmax, rvs)
    written += write_scalebars(out, ums)
    written += write_readme(out, run, snaps, ums, vmin, vmax, rvs, T)

    print(f"\n  {len(written)} files written to {out}")


if __name__ == "__main__":
    main()
