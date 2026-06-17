#!/usr/bin/env python3
"""
build_geometry_multi_grain.py — non-square domain with multiple sediment
"bumps" along the bottom edge, generalizing build_geometry_sediment_grain.py
from one bump to N.

Usage
-----
Run from the project root (PetIGA-3.20/projects/permafrost/):

    python3 preprocess/build_geometry_multi_grain.py [options]

Options:
    --Nx INT    Number of elements in x (default: 240)
    --Ny INT    Number of elements in y (default: 240)
    --P  INT    B-spline degree; continuity is C^{P-1} (default: 2)
    --out PATH  Output PetIGA binary mesh file
                  (default: inputs/geometry/multi_grain_test.dat)
    --plot PATH Output control-mesh PNG plot
                  (default: preprocess/multi_grain_geometry.png)
    --vtk PATH  Output VTK structured grid for visualization
                  (default: preprocess/multi_grain_geometry.vtk)

Examples:
    # Default 240x240, P=2/C1 mesh (the standard production geometry):
    python3 preprocess/build_geometry_multi_grain.py

    # Finer mesh for convergence study:
    python3 preprocess/build_geometry_multi_grain.py --Nx 360 --Ny 360

    # P=1/C0 comparison mesh written to a separate file:
    python3 preprocess/build_geometry_multi_grain.py --P 1 \\
        --out inputs/geometry/multi_grain_test_p1.dat \\
        --plot preprocess/multi_grain_geometry_p1.png \\
        --vtk  preprocess/multi_grain_geometry_p1.vtk

After regenerating the .dat file, update two things in the matching
inputs/geometry/*.opts file:
    1. # DOF_GRID: X X  →  Nx+P in each direction (printed by this script
                            as "shape (control points): (X, X)")
    2. -eps VALUE        →  recompute via preprocess/comp_eps.py if needed

Builds a single-patch NURBS surface for a rectangular domain [0,Lx]x[0,Ly]
whose bottom edge is the SUM of several C-infinity bump functions (one per
sediment grain), each rising smoothly from 0 to its own height/half-width --
no cusps/reflex corners, so the surface has a positive Jacobian as long as
bottom(x) < Ly everywhere.

    bottom(x) = sum_k bump(x; center_k, R_k)
    top(x)    = Ly  (flat)

Geometry degree/continuity: (P, P) with C^{P-1} (single interior knots).
Both axes use an open-uniform knot vector whose Greville abscissae g_i
satisfy sum_i N_i(t) g_i = t exactly (B-splines reproduce linear functions
at their Greville abscissae) -- so:
  - x(u) = Lx * u exactly (control-point x-coords = Lx * g_i),
  - the v-direction is an exact linear interpolation between the bottom and
    top curves: y(u,v) = bottom_y(u) + v*(Ly - bottom_y(u)).
  - bottom_y(u) = bump_field(Lx*u) is sampled at the Greville abscissae,
    giving a high-order (not exact) approximation of the analytic bump
    curve for P>1.

src/initial_conditions.c's FormInitialMultiGrains2D() computes the matching
physical (x,y) for each DOF via GrevilleAbscissae() (reads the actual knot
vector from the IGA axis), so the IC stays consistent with this geometry for
any (P, C^{P-1}) choice made here.

NOTE: the SEDIMENT_GRAINS list below must match -sed_grain_x/-sed_grain_R in
inputs/geometry/2D_multi_grain_test.opts, and _bump() here must match
SedimentBump() in src/initial_conditions.c (summed by SedimentBumpField()).
"""

import argparse

import numpy as np
import matplotlib.pyplot as plt

from igakit.nurbs import NURBS
from igakit.io import PetIGA

# ---------- domain parameters (same as build_geometry_sediment_grain.py) ----------
Lx = 1.0e-4     # domain width  [m]  (2.5× wider than old 4.0e-5 — channel geometry)
Ly = 4.0e-5     # domain height [m]

# Bottom sediment grains: (center_x [m], R [m], height [m])
#   center_x : bump centre position along x
#   R        : half-width (bump support is [center_x-R, center_x+R])
#   height   : peak height of the bump (independent of R)
# Keep supports within [0,Lx] and non-overlapping.
# Must match -sed_grain_x / -sed_grain_R / -sed_grain_h in the .opts file.
#
# 12 bumps tiling [0, Lx] (2R = 0.8e-5 m, same density as proven 5-bump run),
# uniform h=0.2e-5 (same height as proven 5-bump run; max slope ≈0.82, below 45°).
SEDIMENT_GRAINS = [
    (0.4e-5, 0.4e-5, 0.2e-5),    # support [0.0e-5, 0.8e-5]
    (1.2e-5, 0.4e-5, 0.2e-5),    # support [0.8e-5, 1.6e-5]
    (2.0e-5, 0.4e-5, 0.2e-5),    # support [1.6e-5, 2.4e-5]
    (2.8e-5, 0.4e-5, 0.2e-5),    # support [2.4e-5, 3.2e-5]
    (3.6e-5, 0.4e-5, 0.2e-5),    # support [3.2e-5, 4.0e-5]
    (4.4e-5, 0.4e-5, 0.2e-5),    # support [4.0e-5, 4.8e-5]
    (5.2e-5, 0.4e-5, 0.2e-5),    # support [4.8e-5, 5.6e-5]
    (6.0e-5, 0.4e-5, 0.2e-5),    # support [5.6e-5, 6.4e-5]
    (6.8e-5, 0.4e-5, 0.2e-5),    # support [6.4e-5, 7.2e-5]
    (7.6e-5, 0.4e-5, 0.2e-5),    # support [7.2e-5, 8.0e-5]
    (8.4e-5, 0.4e-5, 0.2e-5),    # support [8.0e-5, 8.8e-5]
    (9.2e-5, 0.4e-5, 0.2e-5),    # support [8.8e-5, 1.0e-4]
]

# Ceiling (top-wall) grains — empty: flat top, bumps on bottom only.
# Set to a list of (center, R, height) tuples to add ceiling bumps.
TOP_GRAINS = []

# target element counts. eps is fixed by physics (preprocess/comp_eps.py);
# Nx/Ny set the number of elements -- override via --Nx/--Ny on the CLI.
Nx = 600   # elements in x (Lx/Nx = 1.667e-7 m, same h as old 240x240 on 4.0e-5)
Ny = 240   # elements in y

# basis-function degree; geometry is (P,P) with C^{P-1} (single interior
# knots, maximal smoothness). P=2 gives quadratic, C1 basis functions --
# smoother ice-air interfaces than the previous P=1/C0 mesh at the same
# element count. Override via --P (e.g. --P 1 for a P=1/C0 comparison mesh).
P = 2


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
    for center, R, height in SEDIMENT_GRAINS:
        y = y + _bump(x, center, R, height)
    return y


def _top_bump_field(x):
    """Sum of all TOP_GRAINS bumps (downward displacement from Ly) -- must match TopBumpField()."""
    y = np.zeros_like(x)
    for center, R, height in TOP_GRAINS:
        y = y + _bump(x, center, R, height)
    return y


def open_uniform_knots(N, p):
    """Open-uniform knot vector for N elements, degree p, C^{p-1}
    (single interior knots): p+1-fold end knots, N-1 single interior knots."""
    interior = np.linspace(0.0, 1.0, N + 1)[1:-1]
    return np.concatenate([np.zeros(p + 1), interior, np.ones(p + 1)])


def greville_abscissae(U, p):
    """g_i = mean(U[i+1 .. i+p]) for i = 0 .. (len(U)-p-2) -- matches
    GrevilleAbscissae() in src/initial_conditions.c."""
    n = len(U) - p - 1
    return np.array([np.mean(U[i + 1:i + p + 1]) for i in range(n)])


def build_surface():
    Ux = open_uniform_knots(Nx, P)
    Uy = open_uniform_knots(Ny, P)
    gx = greville_abscissae(Ux, P)   # (Nx+P,) parametric x-DOF positions
    gy = greville_abscissae(Uy, P)   # (Ny+P,) parametric y-DOF positions

    x        = Lx * gx
    bottom_y = _bump_field(x)            # floor rises from 0
    top_y    = Ly - _top_bump_field(x)  # ceiling drops from Ly

    if np.any(bottom_y >= Ly):
        raise ValueError(f"floor bump reaches/exceeds Ly={Ly:.3e} "
                          f"(max={bottom_y.max():.3e}) -- reduce SEDIMENT_GRAINS heights")
    if np.any(top_y <= 0):
        raise ValueError(f"ceiling bump drops to/below 0 "
                          f"(min top_y={top_y.min():.3e}) -- reduce TOP_GRAINS heights")
    if np.any(top_y <= bottom_y):
        raise ValueError(f"floor and ceiling cross "
                          f"(min gap={( top_y - bottom_y).min():.3e}) -- bumps too large")

    nx, ny = len(gx), len(gy)
    ctrl = np.zeros((nx, ny, 3))
    for i in range(nx):
        ctrl[i, :, 0] = x[i]
        ctrl[i, :, 1] = bottom_y[i] + gy * (top_y[i] - bottom_y[i])

    surf = NURBS([Ux, Uy], control=ctrl)
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
    global P, Nx, Ny

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--Nx", type=int, default=Nx,
                         help=f"elements in x, default {Nx}")
    parser.add_argument("--Ny", type=int, default=Ny,
                         help=f"elements in y, default {Ny}")
    parser.add_argument("--P", type=int, default=P,
                         help=f"basis-function degree (C^{{P-1}} continuity), default {P}")
    parser.add_argument("--out", default="inputs/geometry/multi_grain_test.dat",
                         help="output PetIGA mesh file")
    parser.add_argument("--plot", default="preprocess/multi_grain_geometry.png",
                         help="output control-mesh plot")
    parser.add_argument("--vtk", default="preprocess/multi_grain_geometry.vtk",
                         help="output VTK structured grid")
    args = parser.parse_args()
    Nx = args.Nx
    Ny = args.Ny
    P = args.P

    surf = build_surface()

    print("degree:", surf.degree)
    print("shape (control points):", surf.shape)
    print("breaks axis0:", surf.breaks(0).size - 1, "elements")
    print("breaks axis1:", surf.breaks(1).size - 1, "elements")

    plot_surface(surf, args.plot)
    write_vtk(surf, args.vtk)

    PetIGA().write(args.out, surf)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
