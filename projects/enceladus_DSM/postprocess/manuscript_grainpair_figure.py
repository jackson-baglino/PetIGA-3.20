#!/usr/bin/env python3
"""Manuscript figure panels for an axisymmetric grain-pair sintering run.

Produces the pieces of a multi-panel figure as SEPARATE vector files, so the
figure itself is assembled in Inkscape:

    xsec_<tag>.pdf/.svg          full-domain meridional cross-section:
                                 supersaturation in the air, opaque ice body
    ice3d_<tag>.pdf/.svg         shaded 3-D rendering of the ice body (the
                                 phi = 0.5 surface revolved about the
                                 symmetry axis), transparent background
    ice3d_vapour_<tag>.pdf/.svg  the same 3-D body above the symmetry axis,
                                 the vapour cross-section below it

    cbar_supersat_{light,dark}   supersaturation bar (cmocean balance,
                                 white anchored at sigma = 0)
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

import cmocean
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
# The vapour field is drawn as SUPERSATURATION sigma = rho_v/rho_vs - 1 on
# cmocean's `balance`: a diverging map needs a white point that means
# something, and sigma = 0 -- equilibrium with flat ice -- is the only value
# here that does.
#
# ONLY THE BLUE ARM IS EVER USED, and that is a result, not an oversight: at
# this wall humidity the entire domain is undersaturated (sigma runs about
# -2.85e-3 at the chamber wall to -0.47e-3 at the ice), so nothing in the run
# reaches white, let alone red. sigma = 0 still anchors the scale, and the
# colourbar runs up to it so the reader can see how far short the field falls.
CMAP_SAT = cmocean.cm.balance

# Ice ramp: phi = 0 (air) at the dark end, phi = 1 (ice) at the near-white end.
# The 3-D shading reads its tones off this same ramp, so the phi colourbar and
# the rendered body belong to one visual family.
ICE_RAMP = ["#33465a", "#546a80", "#7b90a4", "#a8bccd", "#d5e3ee", "#f7fbff"]
CMAP_ICE = LinearSegmentedColormap.from_list("ice", ICE_RAMP)

# The cut face is a MID tone off the ice ramp, not a near-white one. With the
# phi = 0.5 outline gone the boundary has to be carried by fill contrast
# alone, and the ice is always adjacent to the brightest end of the vapour
# scale (sigma is highest right at the surface) -- so the fill has to be
# clearly darker than near-white. It does not matter that the far field is
# darker still: the ice is never adjacent to that.
ICE_FILL = "#9db1c3"     # = CMAP_ICE(0.55), so the flat cut face and the
                         # shaded 3-D body read as the same material
ICE_EDGE = "#2c3a48"     # silhouette of the 3-D body only; the 2-D cut face
                         # is drawn unstroked
AXIS_INK = "#3a4650"

INK_LIGHT = "#0b0b0b"    # text/rules on a light background
INK_DARK = "#f4f3ee"     # text/rules on a dark background

N_BAND_SAT = 72          # vapour bands
N_BAND_SHADE = 120       # shading bands on the 3-D body -- this is what sets
                         # how smooth the rendered ice reads, so it is high


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------
def sample(run_dir, step, nu=3001, nv=1501):
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


def ice_profile(x, y, phi, simplify_um=0.004):
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


def band_levels(vmin, vmax, n, cmap, pos=None):
    """n band thresholds spanning [vmin, vmax] and the colour for each.

    `pos` maps a data value to its position in the colormap. It exists so a
    diverging map can be anchored on a value rather than on the data range:
    for supersaturation the white point has to land on sigma = 0 wherever the
    data happens to stop, which a plain linear vmin..vmax stretch cannot do.

    The first threshold sits slightly below vmin so the first band covers
    every unmasked cell -- otherwise the region at exactly vmin is empty and
    the figure shows holes where the field is at its floor.
    """
    lv = np.linspace(vmin, vmax, n + 1)[:-1]
    mid = 0.5 * (lv + np.r_[lv[1:], vmax])
    lv = lv.copy()
    lv[0] -= 1e-9 * max(abs(vmax - vmin), 1.0)
    f = pos or (lambda v: (v - vmin) / (vmax - vmin))
    return lv, cmap(np.clip(f(mid), 0.0, 1.0))


# balance's extreme blue is nearly black. The Dirichlet wall sits at exactly
# sigma_min and wraps the whole domain boundary, so running the scale to the
# very end floods the corners with near-black. Starting the arm a little way
# in keeps sigma = 0 anchored at white and still prints.
SAT_DARK_TRIM = 0.09


def sat_pos(sigma, sigma_min):
    """Colormap position for supersaturation on a diverging map whose white
    point is sigma = 0. sigma_min (the most undersaturated value in the run)
    fixes the blue end; sigma = 0 always lands on the midpoint."""
    t = np.clip(np.asarray(sigma, float) / sigma_min, 0.0, 1.0)
    return SAT_DARK_TRIM + (0.5 - SAT_DARK_TRIM) * (1.0 - t)


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


def shade_revolved(ax, profile, n_theta=1441, light=(-0.44, 0.62, 0.65),
                   ambient=0.22, diffuse=0.52, specular=0.13, shine=34.0,
                   rim=0.20, rim_pow=3.0, ao_min=0.44, ao_pow=0.70,
                   n_bands=N_BAND_SHADE, tol=0.06, zorder=3, outline=False):
    """Draw the ice body: the profile revolved about y = 0, Lambert + Blinn
    shaded with a crevice-darkening term, as filled bands of the intensity.

    The camera looks along -z orthographically, so screen (X, Y) = (x, y) in
    micrometres and the silhouette is the profile itself.

    THE CAMERA IS DELIBERATELY STRAIGHT ON, not tilted: Molaro et al.'s
    micrographs are side-on views of the pair, so this is the orientation the
    figure has to be compared against. All of the depth cue therefore has to
    come out of the lighting, which is why there are four terms and not one.

    AMBIENT OCCLUSION IS NOT COSMETIC HERE. Without it the neck renders
    *brighter* than the grains: its fillet sweeps the normal through the
    mirror direction, so the specular term spikes right where the surface is
    most enclosed and should be darkest. The proxy used is the profile's
    depth below its own upper convex hull, which is zero on the grains and
    large exactly in the neck.

    THE RIM TERM is what makes a straight-on view read as a solid rather than
    a flat disc. It brightens where the surface turns away from the camera
    (n . view -> 0, i.e. the silhouette), gated by the diffuse term so only
    the lit limb picks it up -- an unmodulated rim glows all the way round and
    reads as a halo instead of curvature.
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
    rimt = rim * (1.0 - np.clip(nz, 0.0, 1.0)) ** rim_pow * (0.3 + 0.7 * lam)
    inten = np.clip(ao * (ambient + diffuse * lam + specular * spc + rimt),
                    0.0, 1.0)

    # Screen coordinates, back half masked away. At this camera the surface
    # cannot occlude itself, so back-face culling is the whole of the
    # hidden-surface problem.
    X = np.broadcast_to(xp[:, None], (len(xp), n_theta))
    Y = rp[:, None] * c
    Z = np.ma.masked_where(nz <= 1e-9, inten)

    lv, cols = band_levels(0.0, 1.0, n_bands, CMAP_ICE)
    art = filled_bands(ax, X, Y, Z, lv, cols, tol, zorder=zorder,
                       label="shade")

    # Optional silhouette, off by default: the body is meant to read as a lit
    # solid, and a stroke around it reads as a drawn contour. Kept available
    # only as an escape hatch for a background that swallows the limb --
    # checked unstroked on white, cream and near-black, where the edge holds
    # because the highlight sits inboard of the limb, not on it.
    if outline:
        sil = np.vstack([np.column_stack([xp, rp]),
                         np.column_stack([xp[::-1], -rp[::-1]])])
        art.append(ax.plot(sil[:, 0], sil[:, 1], color=ICE_EDGE, lw=0.55,
                           alpha=0.85, solid_joinstyle="round",
                           zorder=zorder + 1)[0])
    return art


# ---------------------------------------------------------------------------
# 2-D vapour cross-section
# ---------------------------------------------------------------------------
def vapour_section(ax, x, y, phi, sigma, smin, stride=2, tol=0.18,
                   n_bands=N_BAND_SAT, zorder=1):
    """Filled bands of supersaturation in the air phase, mirrored about the
    axis. White (sigma = 0) is anchored by `sat_pos`, not by the data range.

    THE FIELD IS DRAWN THROUGH THE ICE, not clipped at phi = 0.5, and the
    opaque ice polygon is what hides it. contourpy drops any cell with a
    masked corner, so a mask at phi = 0.5 pulls the bands up to one cell SHORT
    of the boundary and leaves a ragged sliver of bare page around the ice --
    which the unstroked cut face has nothing to cover. Letting the bands run
    under the ice instead guarantees the overlap. rho_v inside the ice is the
    local equilibrium value, in range and never seen, so nothing is distorted;
    the mask at phi > 0.99 is only a guard.

    Because the boundary comes from the ice polygon and not from the mask, the
    field itself can be decimated (`stride`) at no visible cost.
    """
    xu, yu = x[::stride] / UM, y[::stride] / UM
    ph, sg = phi[::stride, ::stride], sigma[::stride, ::stride]

    yfull = np.concatenate([-yu[:0:-1], yu])
    ph = np.vstack([ph[:0:-1, :], ph])
    sg = np.vstack([sg[:0:-1, :], sg])

    X, Y = np.meshgrid(xu, yfull)
    Z = np.ma.masked_where(ph > 0.99, sg)
    lv, cols = band_levels(smin, 0.0, n_bands, CMAP_SAT,
                           pos=lambda v: sat_pos(v, smin))
    return filled_bands(ax, X, Y, Z, lv, cols, tol, zorder=zorder,
                        label="sigma")


def draw_ice_body(ax, profile, fill=ICE_FILL, zorder=4):
    """Opaque ice cross-section from the phi = 0.5 profile, mirrored.

    Unstroked on purpose: an outline on a flat cut face reads as an annotation
    drawn over the figure rather than as the edge of a solid. The boundary is
    carried by the fill's contrast against the vapour field instead, which is
    why ICE_FILL is warm and the field is not.
    """
    xp, rp = profile[:, 0], profile[:, 1]
    poly = np.vstack([np.column_stack([xp, rp]),
                      np.column_stack([xp[::-1], -rp[::-1]])])
    return ax.fill(poly[:, 0], poly[:, 1], facecolor=fill, edgecolor="none",
                   lw=0.0, joinstyle="round", zorder=zorder)


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


# Everything this script is allowed to delete on a rerun. Narrow on purpose:
# a rerun with different --tags leaves panels from the previous tags behind,
# and stale panels in a figure folder are worse than no panels -- you cannot
# tell by looking which colour scale or which snapshot they came from. Only
# these patterns are swept, so anything else in the folder is left alone.
_OWNED = ("xsec_*.pdf", "xsec_*.svg", "ice3d_*.pdf", "ice3d_*.svg",
          "cbar_*.pdf", "cbar_*.svg", "scalebar_*.pdf", "scalebar_*.svg",
          "preview_*.png", "README.txt")


def clean_output(out):
    """Remove this script's own products from `out`. Returns the count."""
    import glob

    n = 0
    for pat in _OWNED:
        for p in glob.glob(os.path.join(out, pat)):
            try:
                os.remove(p)
                n += 1
            except OSError:
                pass
    return n


# ---------------------------------------------------------------------------
# panels
# ---------------------------------------------------------------------------
def panel_xsec(out, tag, x, y, phi, sigma, profile, smin, ums, stride):
    """Full-domain meridional cross-section, mirrored about the axis."""
    Lx, Ly = x[-1] / UM, y[-1] / UM
    fig, ax = _panel_figure((0, Lx), (-Ly, Ly), ums)
    vapour_section(ax, x, y, phi, sigma, smin, stride=stride)
    draw_ice_body(ax, profile)
    draw_axis_line(ax, 0, Lx)
    return _save(fig, out, f"xsec_{tag}"), fig


def panel_ice3d(out, tag, profile, crop, ums, n_bands=N_BAND_SHADE,
                outline=False):
    """The revolved ice body on its own, transparent background."""
    (x0, x1), ymax = crop
    fig, ax = _panel_figure((x0, x1), (-ymax, ymax), ums)
    shade_revolved(ax, profile, n_bands=n_bands, outline=outline)
    return _save(fig, out, f"ice3d_{tag}"), fig


def panel_ice3d_vapour(out, tag, x, y, phi, sigma, profile, smin, ums,
                       stride, n_bands=N_BAND_SHADE, outline=False):
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

    vapour_section(ax, x, y, phi, sigma, smin, stride=stride)
    _clip(draw_ice_body(ax, profile), ax, 0.0, -Ly, Lx, Ly)
    _clip(shade_revolved(ax, profile, n_bands=n_bands, outline=outline),
          ax, 0.0, 0.0, Lx, Ly)
    draw_axis_line(ax, 0, Lx, lw=0.6)
    return _save(fig, out, f"ice3d_vapour_{tag}"), fig


# ---------------------------------------------------------------------------
# colourbars and scale bars
# ---------------------------------------------------------------------------
def colourbar(out, stem, cmap, vmin, vmax, ticks, ticklabels, label, ink,
              n_bands=256, horizontal=True, pos=None):
    """A standalone vector colourbar: band rectangles, so the file carries
    paths rather than an embedded raster.

    `pos` is the same value -> colormap-position map the panels use, so an
    anchored diverging scale is keyed by exactly the colours it is drawn in.
    """
    if horizontal:
        fig = plt.figure(figsize=(3.0, 0.64))
        ax = fig.add_axes([0.03, 0.45, 0.94, 0.33])
    else:
        fig = plt.figure(figsize=(1.05, 2.7))
        ax = fig.add_axes([0.07, 0.05, 0.26, 0.90])

    edges = np.linspace(vmin, vmax, n_bands + 1)
    mid = 0.5 * (edges[:-1] + edges[1:])
    f = pos or (lambda v: (v - vmin) / (vmax - vmin))
    cols = cmap(np.clip(f(mid), 0.0, 1.0))
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


def write_colourbars(out, smin, orientations=("h", "v")):
    """One bar for the vapour field and one for the ice phase, each in light
    and dark ink and in both orientations.

    The vapour bar runs to sigma = 0 even though the data stops short of it:
    clipped at the data maximum the bar would show blue fading to pale blue
    and read as an ordinary sequential scale, hiding the one thing the
    diverging map is there to say -- where equilibrium is, and that the run
    never gets there.
    """
    made = []
    tk = np.arange(-2.5e-3, 1e-9, 0.5e-3)
    for orient in orientations:
        horiz = orient == "h"
        sfx = "" if horiz else "_vert"
        for mode, ink in (("light", INK_LIGHT), ("dark", INK_DARK)):
            made += colourbar(
                out, f"cbar_supersat_{mode}{sfx}", CMAP_SAT, smin, 0.0,
                tk, [f"{v * 1e3:.1f}" for v in tk],
                r"$\rho_v/\rho_{vs}-1$   ($\times 10^{-3}$)", ink,
                horizontal=horiz, pos=lambda v: sat_pos(v, smin))

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
  cbar_supersat_*     supersaturation sigma = rho_v/rho_vs - 1, x 10^-3
  cbar_ice_*          ice phase phi, over the ice tone ramp the 3-D body is
                      shaded from
  scalebar_*          {sblen:g} um at the panels' own scale

VAPOUR SCALE
sigma = rho_v/rho_vs - 1 on cmocean `balance`, white anchored at sigma = 0
(equilibrium with flat ice at T = {T:g} C). Shared across every time shown, so
the panels are directly comparable.

Over both snapshots sigma runs [{smin:+.3f}, {smax:+.3f}] x 10^-3: the whole
domain is UNDERSATURATED, so only the blue arm of the map appears and nothing
reaches white. That is the result, not a plotting choice -- the pair is net
sublimating into the chamber wall over this window (R_large falls from 100.9
to 97.7 um). The colourbar runs all the way to sigma = 0 so the reader can see
how far short of equilibrium the field stays; the panels themselves stop at
{smax:+.3f} x 10^-3, the value hugging the ice.

NOTHING IS STROKED. Neither the 2-D cut face nor the 3-D body carries an
outline; every edge in these panels is fill contrast. That is what sets the
ice colour: sigma is highest right AT the ice, so the ice is always adjacent
to the palest colour on the scale, and the cut face therefore has to be a MID
blue-grey rather than a near-white one. It is a mid tone of the ramp the 3-D
body is shaded from, so the flat face and the solid read as one material.

`ice3d_*` has a transparent background and was checked unstroked on white,
cream and near-black: the body is a mid blue-grey overall and its highlight
sits inboard of the limb, so the edge holds on all three. If some other ground
does swallow it, rerun with --outline-3d rather than dropping the rim term --
the rim is what keeps a straight-on sphere from reading as a flat disc.

The blue arm stops a little short of balance's darkest end, so the Dirichlet
wall prints rather than going to near-black. The colourbar uses the same
mapping, so it keys the panels exactly.

preview_*.png are checking renders, not figure material.
"""


def write_readme(out, run, snaps, ums, smin, smax, T, sblen=100.0):
    txt = README.format(
        run=run, snaps="".join(f"  {s}\n" for s in snaps), ums=ums,
        ptum=72.0 / ums, sblen=sblen, smin=smin * 1e3, smax=smax * 1e3, T=T)
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
    ap.add_argument("--nu", type=int, default=3001,
                    help="sampling grid across Lx (sets profile smoothness)")
    ap.add_argument("--nv", type=int, default=1501)
    ap.add_argument("--field-stride", type=int, default=2,
                    help="decimation of the vapour field before banding")
    ap.add_argument("--shade-bands", type=int, default=N_BAND_SHADE,
                    help="intensity bands on the 3-D body; lower for smaller "
                         "files, higher for a smoother surface")
    ap.add_argument("--outline-3d", action="store_true",
                    help="stroke the 3-D body's silhouette (off by default)")
    ap.add_argument("--keep-existing", action="store_true",
                    help="do not sweep this script's earlier output from the "
                         "target folder first")
    ap.add_argument("--no-preview", action="store_true")
    args = ap.parse_args()

    run = os.path.abspath(args.run_dir)
    out = args.out or os.path.join(run, "plots", "manuscript")
    os.makedirs(out, exist_ok=True)
    if not args.keep_existing:
        n = clean_output(out)
        if n:
            print(f"  swept {n} earlier file(s) from {out}")
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
        data[tag] = (x, y, phi, rhov / rvs - 1.0, prof)
        t = tmap.get(step, float("nan"))
        line = (f"{tag}: step {step}  t = {t:9.1f} s ({t / 60:6.2f} min)  "
                f"neck radius = {_neck_radius(prof):5.2f} um  "
                f"R_large = {prof[:, 1].max():6.2f} um")
        snaps.append(line)
        print("  " + line)

    smin = min(float(np.min(d[3][d[2] <= 0.5])) for d in data.values())
    smax = max(float(np.max(d[3][d[2] <= 0.5])) for d in data.values())
    print(f"  shared supersaturation scale: sigma = rho_v/rho_vs - 1 in "
          f"[{smin * 1e3:+.3f}, {smax * 1e3:+.3f}] x 1e-3 "
          f"(white anchored at sigma = 0; the field never reaches it)")

    crop = ((args.crop[0], args.crop[1]), args.crop[2])
    ums, stride = args.um_per_inch, args.field_stride
    written = []

    for tag, (x, y, phi, sigma, prof) in data.items():
        made = []
        p, f1 = panel_xsec(out, tag, x, y, phi, sigma, prof, smin, ums,
                           stride)
        made.append((p, f1, f"xsec_{tag}"))
        p, f2 = panel_ice3d(out, tag, prof, crop, ums,
                            n_bands=args.shade_bands,
                            outline=args.outline_3d)
        made.append((p, f2, f"ice3d_{tag}"))
        p, f3 = panel_ice3d_vapour(out, tag, x, y, phi, sigma, prof, smin,
                                   ums, stride, n_bands=args.shade_bands,
                                   outline=args.outline_3d)
        made.append((p, f3, f"ice3d_vapour_{tag}"))
        for p, fig, stem in made:
            written += p
            if not args.no_preview:
                fig.savefig(os.path.join(out, f"preview_{stem}.png"), dpi=220,
                            transparent=False, facecolor="#ffffff")
            plt.close(fig)

    written += write_colourbars(out, smin)
    written += write_scalebars(out, ums)
    written += write_readme(out, run, snaps, ums, smin, smax, T)

    print(f"\n  {len(written)} files written to {out}")


if __name__ == "__main__":
    main()
