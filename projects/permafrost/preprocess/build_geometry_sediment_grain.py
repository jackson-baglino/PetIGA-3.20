#!/usr/bin/env python3
"""
build_geometry_sediment_grain.py — prototype non-square domain with igakit.

Builds a single-patch NURBS surface for a rectangular domain [0,Lx]x[0,Ly]
whose bottom edge has a smooth "bump" raised into the domain, representing
the footprint of a sediment grain (a smoothed-out semicircular bite).

Construction (avoids the Coons-patch folding seen with a sharp
line-arc-line bite + cad.coons -- see git history for that attempt):

    bottom(x) = a C-infinity bump function, zero (with all derivatives
                vanishing) outside |x - Lx/2| < R_sed, rising smoothly to
                height R_sed at the center -- NO cusps/reflex corners.
    top(x)    = Ly  (flat)

Both curves are parametrized directly by u = x/Lx (same knot vector,
piecewise-linear, uniformly spaced control points), so
cad.ruled(bottom, top) produces perfectly VERTICAL v-isolines at every x
-- the Jacobian is guaranteed positive everywhere since bottom(x) < Ly.
The v-direction (2 control points from ruled()) is then knot-refined to
get Ny elements.

This is a standalone prototype: it does not modify permafrost2.c or any
solver input files (beyond the already-committed -geom_file IGARead path).
Run it and inspect the PNG/VTK before doing anything else.
"""

import numpy as np
import matplotlib.pyplot as plt

from igakit import cad
from igakit.nurbs import NURBS
from igakit.io import PetIGA

# ---------- domain parameters (same order of magnitude as
# inputs/geometry/2D_two_ice_grains_boundary.opts) ----------
Lx = 4.0e-5     # domain width  [m]
Ly = 4.0e-5     # domain height [m]
R_sed = 1.0e-5  # sediment-grain bump half-width / height [m]

# target element counts -- the resulting geometry is degree (1,1) with
# C0 interior knots, matching solver.opts' -p 1 -C 0; -geom_file/IGARead
# reads this directly and overrides -p/-C/-Nx/-Ny for this run.
Nx = 80   # elements in x
Ny = 80   # elements in y


def _bump(x, center, R, height):
    """C-infinity bump: height*exp(1 - 1/(1-t^2)) for |t|<1 (t=(x-c)/R),
    0 outside -- touches 0 with all derivatives at t=+-1 (no cusp).

    NOTE: duplicated in src/initial_conditions.c as SedimentBump() so the
    -ic_type two_ice_grains_boundary IC can place grains at the correct
    physical (x,y) on this geometry. Keep R_sed here in sync with
    -geom_bump_R in inputs/geometry/2D_sediment_grain_test.opts."""
    t = (x - center) / R
    out = np.zeros_like(x)
    mask = np.abs(t) < 1.0
    out[mask] = height * np.exp(1.0 - 1.0 / (1.0 - t[mask] ** 2))
    return out


def build_surface():
    xm = 0.5 * Lx

    # bottom curve: graph y = bump(x), piecewise-linear B-spline through
    # Nx+1 uniformly-spaced points -- since x_i are uniform and the knot
    # vector is uniform/clamped, x(u) = u*Lx exactly.
    x = np.linspace(0.0, Lx, Nx + 1)
    y = _bump(x, xm, R_sed, R_sed)
    ctrl_bottom = np.zeros((Nx + 1, 3))
    ctrl_bottom[:, 0] = x
    ctrl_bottom[:, 1] = y
    interior = np.linspace(0.0, 1.0, Nx + 1)[1:-1]
    U = np.concatenate([[0.0, 0.0], interior, [1.0, 1.0]])
    bottom = NURBS([U], control=ctrl_bottom)

    # top curve: flat line, same u <-> x parametrization
    top = cad.line((0.0, Ly), (Lx, Ly))

    surf = cad.ruled(bottom, top)  # v: bottom -> top, vertical isolines
    return surf


def refine_surface(surf):
    # v-direction currently has 1 element (2 control points, degree 1);
    # insert interior knots for Ny elements.
    new_v = np.linspace(0.0, 1.0, Ny + 1)[1:-1]
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
