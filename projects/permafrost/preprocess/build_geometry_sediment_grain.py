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

# target element counts (the resulting geometry's own degree/continuity --
# whatever cad.coons/cad.join produce, typically p=2 from the circular arc --
# is read directly by IGARead and overrides solver.opts' -p/-C for this run)
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


def _arclength_knots(curve, breaks, n_elem, n_sample=2000):
    """Interior knot values giving ~equal physical arc-length per element,
    found by sampling `curve` densely and inverting cumulative arc length."""
    u0, u1 = breaks[0], breaks[-1]
    s = np.linspace(u0, u1, n_sample)
    pts = curve(s)
    seglen = np.linalg.norm(np.diff(pts[:, :2], axis=0), axis=1)
    cum = np.concatenate([[0.0], np.cumsum(seglen)])
    total = cum[-1]
    targets_len = np.linspace(0.0, total, n_elem + 1)[1:-1]
    targets_u = np.interp(targets_len, cum, s)
    return np.setdiff1d(targets_u, breaks)


def refine_surface(surf):
    # bottom edge (v=0) is line/arc/line -> refine axis 0 by arc length
    # along it so elements are roughly uniform in physical space across
    # the line/arc/line transitions.
    bottom_edge = surf.boundary(1, 0)  # curve at v = v_min, varies in u
    breaks0 = surf.breaks(0)
    new_u = _arclength_knots(bottom_edge, breaks0, Nx)
    if new_u.size:
        surf = surf.refine(0, new_u)

    # left edge (u=0) is a straight line -> arc-length == uniform parameter,
    # but compute it the same way for consistency.
    left_edge = surf.boundary(0, 0)  # curve at u = u_min, varies in v
    breaks1 = surf.breaks(1)
    new_v = _arclength_knots(left_edge, breaks1, Ny)
    if new_v.size:
        surf = surf.refine(1, new_v)

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


def write_vtk(surf, fname, n_per_elem=4):
    """Write the physical mesh as a legacy VTK structured grid, sampling
    n_per_elem points per element (in each direction) so curved element
    edges (e.g. the arc) are visible, not just straight lines between
    element corners."""
    breaks0 = surf.breaks(0)
    breaks1 = surf.breaks(1)

    u = np.unique(np.concatenate([
        np.linspace(breaks0[i], breaks0[i + 1], n_per_elem, endpoint=False)
        for i in range(len(breaks0) - 1)
    ] + [[breaks0[-1]]]))
    v = np.unique(np.concatenate([
        np.linspace(breaks1[i], breaks1[i + 1], n_per_elem, endpoint=False)
        for i in range(len(breaks1) - 1)
    ] + [[breaks1[-1]]]))

    pts = surf(u, v)  # shape (nu, nv, 3)
    nu, nv = pts.shape[0], pts.shape[1]

    with open(fname, "w") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("sediment_grain_geometry\n")
        f.write("ASCII\n")
        f.write("DATASET STRUCTURED_GRID\n")
        f.write(f"DIMENSIONS {nu} {nv} 1\n")
        f.write(f"POINTS {nu * nv} float\n")
        # VTK structured-grid point order: x fastest, then y, then z
        for j in range(nv):
            for i in range(nu):
                x, y, _ = pts[i, j]
                f.write(f"{x:.8e} {y:.8e} 0.0\n")

    print(f"wrote {fname}")


def main():
    surf = build_surface()
    surf = refine_surface(surf)

    print("degree:", surf.degree)
    print("shape (control points):", surf.shape)
    print("breaks axis0:", surf.breaks(0).size - 1, "elements")
    print("breaks axis1:", surf.breaks(1).size - 1, "elements")

    plot_surface(surf, "preprocess/sediment_grain_geometry.png")
    write_vtk(surf, "preprocess/sediment_grain_geometry.vtk")

    out = "inputs/geometry/sediment_grain_test.dat"
    PetIGA().write(out, surf)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
