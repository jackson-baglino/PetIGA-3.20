#!/usr/bin/env python3
"""
build_geometry_multi_grain.py — non-square domain with multiple sediment
"bumps" along the bottom edge, generalizing build_geometry_sediment_grain.py
from one bump to N.

Builds a single-patch NURBS surface for a rectangular domain [0,Lx]x[0,Ly]
whose bottom edge is the SUM of several C-infinity bump functions (one per
sediment grain), each rising smoothly from 0 to its own height/half-width --
no cusps/reflex corners, so cad.ruled(bottom, top) still produces a
positive-Jacobian surface as long as bottom(x) < Ly everywhere.

    bottom(x) = sum_k bump(x; center_k, R_k)
    top(x)    = Ly  (flat)

Both curves are parametrized directly by u = x/Lx (same knot vector,
piecewise-linear, uniformly spaced control points), so x(u) = u*Lx exactly --
same construction as build_geometry_sediment_grain.py.

NOTE: the SEDIMENT_GRAINS list below must match -sed_grain_x/-sed_grain_R in
inputs/geometry/2D_multi_grain_test.opts, and _bump() here must match
SedimentBump() in src/initial_conditions.c (summed by SedimentBumpField()).
"""

import numpy as np
import matplotlib.pyplot as plt

from igakit import cad
from igakit.nurbs import NURBS
from igakit.io import PetIGA

# ---------- domain parameters (same as build_geometry_sediment_grain.py) ----------
Lx = 4.0e-5     # domain width  [m]
Ly = 4.0e-5     # domain height [m]

# Sediment grains: (center_x [m], R [m]) -- each contributes a bump of
# half-width == height R, centered at center_x. Keep supports
# (center_x +/- R) within [0,Lx] and non-overlapping for distinct,
# well-separated grains. Must match -sed_grain_x/-sed_grain_R in
# inputs/geometry/2D_multi_grain_test.opts.
SEDIMENT_GRAINS = [
    (0.6e-5, 0.6e-5),   # support [0.0e-5, 1.2e-5] -- touches x=0
    (2.0e-5, 0.6e-5),   # support [1.4e-5, 2.6e-5]
    (3.4e-5, 0.6e-5),   # support [2.8e-5, 4.0e-5] -- touches x=Lx
]

# target element counts -- the resulting geometry is degree (1,1) with
# C0 interior knots, matching solver.opts' -p 1 -C 0; -geom_file/IGARead
# reads this directly and overrides -p/-C/-Nx/-Ny for this run.
# 160x160 (vs. 80x80) doubles eps/h to ~1.9 -- the previous 80x80 mesh had
# eps (4.6648e-07) barely smaller than h (5.0e-7), giving a jagged,
# under-resolved ice-air interface.
Nx = 160   # elements in x
Ny = 160   # elements in y


def _bump(x, center, R, height):
    """C-infinity bump: height*exp(1 - 1/(1-t^2)) for |t|<1 (t=(x-c)/R),
    0 outside -- touches 0 with all derivatives at t=+-1 (no cusp).

    NOTE: duplicated in src/initial_conditions.c as SedimentBump()."""
    t = (x - center) / R
    out = np.zeros_like(x)
    mask = np.abs(t) < 1.0
    out[mask] = height * np.exp(1.0 - 1.0 / (1.0 - t[mask] ** 2))
    return out


def _bump_field(x):
    """Sum of all SEDIMENT_GRAINS bumps -- must match SedimentBumpField()."""
    y = np.zeros_like(x)
    for center, R in SEDIMENT_GRAINS:
        y = y + _bump(x, center, R, R)
    return y


def build_surface():
    # bottom curve: graph y = bump_field(x), piecewise-linear B-spline
    # through Nx+1 uniformly-spaced points -- since x_i are uniform and the
    # knot vector is uniform/clamped, x(u) = u*Lx exactly.
    x = np.linspace(0.0, Lx, Nx + 1)
    y = _bump_field(x)

    if np.any(y >= Ly):
        raise ValueError(f"bump field reaches/exceeds Ly={Ly:.3e} "
                          f"(max={y.max():.3e}) -- reduce grain heights/Ly")

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
    ax.set_title('Rectangle with multiple sediment-grain bumps (control mesh in red)')
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    fig.tight_layout()
    fig.savefig(fname, dpi=150)
    print(f"wrote {fname}")


def write_vtk(surf, fname, n_per_elem=4):
    """Write the physical mesh as a legacy VTK structured grid, sampling
    n_per_elem points per element (in each direction) so curved element
    edges are visible, not just straight lines between element corners."""
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
        f.write("multi_grain_geometry\n")
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

    plot_surface(surf, "preprocess/multi_grain_geometry.png")
    write_vtk(surf, "preprocess/multi_grain_geometry.vtk")

    out = "inputs/geometry/multi_grain_test.dat"
    PetIGA().write(out, surf)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
