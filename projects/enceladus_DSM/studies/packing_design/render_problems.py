#!/usr/bin/env python3
"""Render what the solver actually sees, to make the packing problems visible.

Three figures, one per problem:

  problem1_throats.png   a near-tangent pair is SOLID to the solver
  problem2_fragments.png the usable pore is in hundreds of disconnected pockets
  problem3_trade.png     opening the pore breaks the solid, at every porosity

phi is built the way the solver builds it -- the ADDITIVE Molaro convention
(-ic_grain_union 0), phi = clamp(sum_k [0.5 - 0.5 tanh(0.5 (r_k - R_k)/eps)]).
That matters: two grains each contribute ~0.5 at a near-tangency, so their
tails SUM to ~1 and the gap fills in. At the midpoint of a gap g the sum is
1 - tanh(g/4eps), so phi only falls below 0.01 once g > 10.6 eps.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "preprocess"))
import packing_lib as pl                                    # noqa: E402
import figstyle as fs                                       # noqa: E402

HERE = Path(__file__).parent
EPS = 45e-9                      # the eps these packings were sized for
BAND = 9.2 * EPS                 # nominal 1%-99% width of one interface
RAVE = 2.5e-6
RASTER = 2048            # must match connectivity.csv: percolation is
                         # decided at pixel scale, so a figure at another
                         # raster can contradict the table


def load(d: Path):
    m = json.loads((d / "metadata.json").read_text())
    rows = [l.split() for l in (d / "grains.dat").read_text().splitlines()
            if l and not l.startswith("#")]
    raw = np.array([[float(v) for v in r] for r in rows[1:]])
    c, r = raw[:, :2], raw[:, 2]
    k = ((c[:, 0] >= 0) & (c[:, 0] < m["Lx"])
         & (c[:, 1] >= 0) & (c[:, 1] < m["Ly"]))
    return m, c[k].copy(), r[k].copy()


def phi_field(cen, rad, Lx, Ly, eps, x0, y0, w, h, n):
    """Additive phi on an n x n window at (x0, y0) of size (w, h)."""
    xs = np.linspace(x0, x0 + w, n)
    ys = np.linspace(y0, y0 + h, n)
    X, Y = np.meshgrid(xs, ys)
    phi = np.zeros_like(X)
    offs = [(dx * Lx, dy * Ly) for dx in (-1, 0, 1) for dy in (-1, 0, 1)]
    reach = 10.0 * eps
    for (cx, cy), R in zip(cen, rad):
        for ox, oy in offs:
            px, py = cx + ox, cy + oy
            if (px + R + reach < x0 or px - R - reach > x0 + w
                    or py + R + reach < y0 or py - R - reach > y0 + h):
                continue
            r = np.hypot(X - px, Y - py)
            phi += 0.5 - 0.5 * np.tanh(0.5 * (r - R) / eps)
    return np.clip(phi, 0.0, 1.0)


# ---------------------------------------------------------------- problem 1
def _mark(ax, pa, pb, mid):
    """Ring the throat under discussion, so the reader knows which gap."""
    t = ((pa + pb) * 0.5 - mid) * 1e6
    rr = 0.34 * np.hypot(*(pb - pa)) * 1e6
    ax.add_patch(plt.Circle(tuple(t), rr, fill=False, ec="k", lw=1.8,
                            ls=(0, (4, 3)), zorder=6))
    for pt in (pa, pb):
        ax.plot(*((pt - mid) * 1e6), "+", color="k", ms=9, mew=1.6, zorder=6)


def problem1(d: Path):
    m, cen, rad = load(d)
    Lx, Ly = m["Lx"], m["Ly"]
    bonds = pl.delaunay_bonds(cen, Lx, Ly, True, True)
    i, j = bonds[:, 0], bonds[:, 1]
    dd = cen[j] + bonds[:, 2:4] * np.array([Lx, Ly]) - cen[i]
    gap = np.hypot(dd[:, 0], dd[:, 1]) - (rad[i] + rad[j])

    # one ambiguous pair (a gap the band swallows) and one honest channel
    amb = np.flatnonzero((gap > 0.25 * BAND) & (gap < 0.6 * BAND))
    opn = np.flatnonzero((gap > 3.0 * BAND) & (gap < 6.0 * BAND))
    picks = [(amb[len(amb) // 2], "gap = %.2f x band\nSOLID to the solver"),
             (opn[len(opn) // 2], "gap = %.2f x band\na real channel")]

    fig, axes = plt.subplots(2, 2, figsize=(9.0, 8.4))
    for row, (k, tmpl) in enumerate(picks):
        a, b = int(i[k]), int(j[k])
        pa, pb = cen[a], cen[b] + bonds[k, 2:4] * np.array([Lx, Ly])
        mid = 0.5 * (pa + pb)
        half = 1.15 * (np.hypot(*(pb - pa)) * 0.5 + max(rad[a], rad[b]))
        x0, y0 = mid[0] - half, mid[1] - half
        w = h = 2 * half
        n = 420
        phi = phi_field(cen, rad, Lx, Ly, EPS, x0, y0, w, h, n)
        ext = [(x0 - mid[0]) * 1e6, (x0 + w - mid[0]) * 1e6,
               (y0 - mid[1]) * 1e6, (y0 + h - mid[1]) * 1e6]

        # left: what you drew
        ax = axes[row][0]
        for (cx, cy), R in zip(cen, rad):
            for ox in (-Lx, 0, Lx):
                for oy in (-Ly, 0, Ly):
                    ax.add_patch(plt.Circle(((cx + ox - mid[0]) * 1e6,
                                             (cy + oy - mid[1]) * 1e6),
                                            R * 1e6, facecolor=fs.C[0],
                                            alpha=0.30, edgecolor=fs.C[0],
                                            lw=1.2))
        ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3])
        ax.set_aspect("equal")
        _mark(ax, pa, pb, mid)
        fs.style(ax, "x  [um]", "y  [um]",
                 "what you drew" if row == 0 else "", logy=False)

        # right: what the solver sees
        ax = axes[row][1]
        cmap = ListedColormap([fs.C[1], "#dddddd", fs.C[0]])
        cls = np.digitize(phi, [0.01, 0.99])       # 0 pore, 1 band, 2 solid
        ax.imshow(cls, origin="lower", extent=ext, cmap=cmap, vmin=0, vmax=2,
                  interpolation="nearest")
        ax.contour(np.linspace(ext[0], ext[1], n), np.linspace(ext[2], ext[3], n),
                   phi, levels=[0.5], colors="k", linewidths=1.0)
        ax.set_aspect("equal")
        _mark(ax, pa, pb, mid)
        fs.style(ax, "x  [um]", "",
                 "what the solver sees" if row == 0 else "", logy=False)
        ax.text(0.5, 0.02, tmpl % (gap[k] / BAND), transform=ax.transAxes,
                ha="center", va="bottom", fontsize=fs.FS_NOTE, color=fs.INK,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="none",
                          alpha=0.85))

    handles = [plt.Rectangle((0, 0), 1, 1, fc=fs.C[0], label=r"solid  $\phi>0.99$"),
               plt.Rectangle((0, 0), 1, 1, fc="#dddddd", label=r"diffuse band"),
               plt.Rectangle((0, 0), 1, 1, fc=fs.C[1], label=r"open pore  $\phi<0.01$")]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               fontsize=fs.FS_LEG, bbox_to_anchor=(0.5, -0.01))
    fig.suptitle("Problem 1: a gap narrower than the diffuse band is solid",
                 fontsize=fs.FS_TITLE, x=0.02, ha="left")
    fig.tight_layout(rect=[0, 0.045, 1, 0.96])
    fs.save(fig, HERE, "problem1_throats", dpi=170)


# ---------------------------------------------------------------- problem 2
def problem2(d: Path):
    m, cen, rad = load(d)
    Lx, Ly = m["Lx"], m["Ly"]
    R = RASTER
    solid = pl.rasterize(cen, rad, Lx, Ly, R, R, True, True)
    op = pl.open_pore(solid, BAND, Lx / R, True, True)
    lab, n = ndimage.label(op, structure=np.ones((3, 3), dtype=bool))
    sizes = np.bincount(lab.ravel())[1:]
    big = int(np.argmax(sizes)) + 1

    rng = np.random.default_rng(3)
    cols = rng.uniform(0.25, 0.95, size=(n + 1, 3))
    cols[0] = [0.13, 0.13, 0.13]                   # solid + band = dark
    cols[big] = [0.85, 0.30, 0.05]                 # the largest pocket
    img = ListedColormap(cols)(lab / max(n, 1) * 0 + lab)  # index directly
    ext = [0, Lx * 1e6, 0, Ly * 1e6]

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.0, 5.4))
    a1.imshow(~solid, origin="lower", extent=ext, cmap="gray",
              interpolation="nearest")
    a1.set_aspect("equal")
    fs.style(a1, "x  [um]", "y  [um]",
             "(a)  the pore as drawn — looks connected", logy=False)
    a1.text(0.02, 0.98, "white = pore\nblack = ice", transform=a1.transAxes,
            va="top", ha="left", fontsize=fs.FS_NOTE, color=fs.INK,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="none", alpha=0.9))

    a2.imshow(ListedColormap(cols)(lab), origin="lower", extent=ext,
              interpolation="nearest")
    a2.set_aspect("equal")
    fs.style(a2, "x  [um]", "",
             "(b)  the pore the solver can use", logy=False)
    a2.text(0.02, 0.98,
            f"{n} disconnected pockets\nlargest (orange) holds "
            f"{sizes.max()/op.sum():.0%} of it\nnone reaches across the cell",
            transform=a2.transAxes, va="top", ha="left", fontsize=fs.FS_NOTE,
            color=fs.INK,
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="none", alpha=0.9))
    fig.suptitle("Problem 2: vapour cannot cross the domain",
                 fontsize=fs.FS_TITLE, x=0.02, ha="left")
    fig.tight_layout(rect=[0, 0.02, 1, 0.94])
    fs.save(fig, HERE, "problem2_fragments", dpi=170)


# ---------------------------------------------------------------- problem 3
def _largest(mask, diagonal):
    """Boolean mask of the single largest connected cluster."""
    st = np.ones((3, 3), dtype=bool) if diagonal else None
    lab, n = ndimage.label(mask, structure=st)
    if n == 0:
        return np.zeros_like(mask)
    sizes = np.bincount(lab.ravel())[1:]
    return lab == (int(np.argmax(sizes)) + 1)


def problem3(dirs):
    """Colour ONLY the largest cluster of each phase, so 'spans the cell' is
    something you can see rather than something you take on trust."""
    fig, axes = plt.subplots(1, len(dirs), figsize=(10.0, 4.6))
    for ax, d in zip(axes, dirs):
        m, cen, rad = load(Path(d))
        Lx, Ly = m["Lx"], m["Ly"]
        R = RASTER
        solid = pl.rasterize(cen, rad, Lx, Ly, R, R, True, True)
        op = pl.open_pore(solid, BAND, Lx / R, True, True)
        sx, sy, _, _ = pl.percolates(solid)
        ox, oy, _, _ = pl.percolates(op, diagonal=True)

        rgb = np.full(solid.shape + (3,), 0.97)
        rgb[solid] = [0.78, 0.85, 0.92]                 # ice, not in the backbone
        rgb[op] = [0.97, 0.88, 0.80]                    # usable pore, not backbone
        rgb[_largest(solid, False)] = [0.11, 0.35, 0.62]   # the ice backbone
        rgb[_largest(op, True)] = [0.84, 0.33, 0.04]       # the biggest pore pocket
        ax.imshow(rgb, origin="lower", extent=[0, Lx * 1e6, 0, Ly * 1e6],
                  interpolation="nearest")
        ax.set_aspect("equal")

        def span(x, y):
            return "x and y" if x and y else ("y only" if y else
                                              ("x only" if x else "NEITHER"))
        phi_v = 1.0 - solid.mean()
        fs.style(ax, "x  [um]", "y  [um]" if ax is axes[0] else "",
                 f"porosity {phi_v:.2f}", logy=False)
        ax.text(0.5, -0.30,
                f"ice backbone spans:  {span(sx, sy)}\n"
                f"pore pocket spans:   {span(ox, oy)}",
                transform=ax.transAxes, va="top", ha="center",
                fontsize=fs.FS_NOTE, color=fs.INK, family="monospace")

    handles = [
        plt.Rectangle((0, 0), 1, 1, fc=[0.11, 0.35, 0.62],
                      label="largest ice cluster"),
        plt.Rectangle((0, 0), 1, 1, fc=[0.78, 0.85, 0.92], label="other ice"),
        plt.Rectangle((0, 0), 1, 1, fc=[0.84, 0.33, 0.04],
                      label="largest usable-pore cluster"),
        plt.Rectangle((0, 0), 1, 1, fc=[0.97, 0.88, 0.80], label="other pore"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
               fontsize=fs.FS_LEG, bbox_to_anchor=(0.5, 0.005))
    fig.suptitle("Problem 3: no porosity gives both a spanning ice network "
                 "and a spanning pore", fontsize=fs.FS_TITLE, x=0.02, ha="left")
    fig.tight_layout(rect=[0, 0.21, 1, 0.93])
    fs.save(fig, HERE, "problem3_trade", dpi=170)


if __name__ == "__main__":
    ref = ROOT / "projects/enceladus_DSM/inputs/packings/periodic_xy/phi0.325_Rave2.5um_seed1"
    ref = Path(str(ref).replace("/projects/enceladus_DSM/projects/", "/projects/"))
    problem1(ref)
    problem2(ref)
    packs = HERE / "packings"
    problem3([packs / "phi0.30_seed1", packs / "phi0.45_seed1", packs / "phi0.50_seed2"])
