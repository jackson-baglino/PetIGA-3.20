#!/usr/bin/env python3
"""
build_geometry_sediment_grain.py — prototype non-square domain with igakit.

Builds a single-patch NURBS surface for a rectangular domain [0,Lx]x[0,Ly]
with a semicircular "bite" removed from the bottom edge, representing the
footprint of a sediment grain. The boundary curves are:

    left    : line   (0,0)    -> (0,Ly)
    right   : line   (Lx,0)   -> (Lx,Ly)
    top     : line   (0,Ly)   -> (Lx,Ly)
    bottom  : line   (0,0)    -> (Lx/2-R,0)
              + arc   (Lx/2-R,0) -> (Lx/2+R,0)   (concave, center (Lx/2,0))
              + line  (Lx/2+R,0) -> (Lx,0)

The bottom curve is assembled with cad.join() and the four boundary curves
are assembled into a 2D surface patch with cad.coons(). The result is
refined to match the solver's discretization (p=2, C=1) at a target
element count, plotted for a visual sanity check, and written out as a
PetIGA-readable geometry file.

This is a standalone prototype: it does not modify permafrost2.c or any
solver input files. Run it and inspect the PNG before doing anything else.
"""

import numpy as np
import matplotlib.pyplot as plt

from igakit import cad
from igakit.io import PetIGA

# ---------- domain parameters (same order of magnitude as
# inputs/geometry/2D_two_ice_grains_boundary.opts) ----------
Lx = 4.0e-5     # domain width  [m]
Ly = 4.0e-5     # domain height [m]
R_sed = 1.0e-5  # sediment-grain bite radius [m]

# target discretization (match solver.opts: p=2, C=1)
p_target = 2
Nx = 80   # elements in x
Ny = 80   # elements in y


def build_surface():
    xm = 0.5 * Lx

    # --- bottom edge: line + concave arc + line ---
    left_seg = cad.line((0.0, 0.0), (xm - R_sed, 0.0))
    arc_seg = cad.circle(radius=R_sed, center=(xm, 0.0), angle=(np.pi, 0.0))
    right_seg = cad.line((xm + R_sed, 0.0), (Lx, 0.0))

    bottom = cad.join(cad.join(left_seg, arc_seg, axis=0), right_seg, axis=0)

    # --- other three edges: straight lines ---
    left = cad.line((0.0, 0.0), (0.0, Ly))
    right = cad.line((Lx, 0.0), (Lx, Ly))
    top = cad.line((0.0, Ly), (Lx, Ly))

    surf = cad.coons(((left, right), (bottom, top)))
    return surf


def refine_surface(surf):
    # elevate to the target polynomial degree (lines/coons start at p=1)
    for axis in (0, 1):
        p_cur = surf.degree[axis]
        if p_cur < p_target:
            surf = surf.elevate(axis, p_target - p_cur)

    # insert interior knots for the requested element counts
    for axis, n_elem in ((0, Nx), (1, Ny)):
        breaks = surf.breaks(axis)
        u0, u1 = breaks[0], breaks[-1]
        targets = np.linspace(u0, u1, n_elem + 1)[1:-1]
        new_knots = np.setdiff1d(targets, breaks)
        if new_knots.size:
            surf = surf.refine(axis, new_knots)

    return surf


def plot_surface(surf, fname):
    nu, nv = 200, 200
    u = np.linspace(surf.knots[0][0], surf.knots[0][-1], nu)
    v = np.linspace(surf.knots[1][0], surf.knots[1][-1], nv)
    pts = surf(u, v)  # shape (nu, nv, 3)

    fig, ax = plt.subplots(figsize=(6, 6))

    # physical-space boundary
    ax.plot(pts[:, 0, 0], pts[:, 0, 1], 'k-', lw=1.5)
    ax.plot(pts[:, -1, 0], pts[:, -1, 1], 'k-', lw=1.5)
    ax.plot(pts[0, :, 0], pts[0, :, 1], 'k-', lw=1.5)
    ax.plot(pts[-1, :, 0], pts[-1, :, 1], 'k-', lw=1.5)

    # element grid lines (at knot breaks)
    for ub in surf.breaks(0):
        line = surf([ub], v)[0]
        ax.plot(line[:, 0], line[:, 1], 'b-', lw=0.3)
    for vb in surf.breaks(1):
        line = surf(u, [vb])[:, 0]
        ax.plot(line[:, 0], line[:, 1], 'b-', lw=0.3)

    # control points
    ctrl = surf.points
    ax.plot(ctrl[..., 0].ravel(), ctrl[..., 1].ravel(), 'r.', ms=3)

    ax.set_aspect('equal')
    ax.set_title('Rectangle with sediment-grain bite (control mesh in red)')
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    fig.tight_layout()
    fig.savefig(fname, dpi=150)
    print(f"wrote {fname}")


def main():
    surf = build_surface()
    surf = refine_surface(surf)

    print("degree:", surf.degree)
    print("shape (control points):", surf.shape)
    print("breaks axis0:", surf.breaks(0).size - 1, "elements")
    print("breaks axis1:", surf.breaks(1).size - 1, "elements")

    plot_surface(surf, "preprocess/sediment_grain_geometry.png")

    out = "inputs/geometry/sediment_grain_test.dat"
    PetIGA().write(out, surf)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
