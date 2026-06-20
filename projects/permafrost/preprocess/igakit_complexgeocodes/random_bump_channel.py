#!/usr/bin/env python3
"""
random_bump_channel.py — rectangular domain with randomly-sized/spaced bumps
along the bottom edge, to look like the rough wall of a fracture in
rock/sediment rather than a row of identical, evenly-spaced grains.

This is a first, deliberately small step: same flat-top-rectangle /
additive-bump / Greville-abscissae construction as
../build_geometry_multi_grain.py (the project's production geometry
builder), just with the bump list drawn at random instead of hand-picked.
No boundary-layer y-refinement, no ceiling bumps, no wiring into
initial_conditions.c/opts files yet -- that comes later, once the basic
shape looks right.

Bump primitive (same as build_geometry_multi_grain.py's _bump / src/
initial_conditions.c's SedimentBump): a C-infinity bump
    height * exp(1 - 1/(1 - t^2)),   t = (x - center) / R,   |t| < 1
so the floor is smooth (no cusps) everywhere, and the surface Jacobian
stays positive as long as bottom(x) < Ly everywhere.

Bump placement: rejection-sampled along [0, Lx] so that no two bump
supports [center-R, center+R] overlap by more than `gap_factor` of their
combined half-widths -- keeps the random asperities visually distinct
(like a rough fracture wall) instead of frequently merging into wider
plateaus.

Usage:
    python3 random_bump_channel.py [--seed INT] [--n-bumps INT] ...
    (see --help for the full list; defaults produce a small, quick mesh)
"""

import argparse

import numpy as np
import matplotlib.pyplot as plt

from igakit.nurbs import NURBS
from igakit.io import PetIGA, VTK


def _bump(x, center, R, height):
    t = (x - center) / R
    out = np.zeros_like(x)
    mask = np.abs(t) < 1.0
    out[mask] = height * np.exp(1.0 - 1.0 / (1.0 - t[mask] ** 2))
    return out


def _bump_field(x, bumps):
    y = np.zeros_like(x)
    for center, R, height in bumps:
        y = y + _bump(x, center, R, height)
    return y


def random_bumps(rng, Lx, n_bumps, R_range, h_range, gap_factor, margin, max_tries=2000):
    """Rejection-sample n_bumps non-overlapping (center, R, height) triples
    within [margin, Lx-margin]. Raises if it can't place them all -- that
    means the domain is too crowded for the requested count/sizes."""
    bumps = []
    for _ in range(n_bumps):
        for _try in range(max_tries):
            R = rng.uniform(*R_range)
            height = rng.uniform(*h_range)
            center = rng.uniform(margin + R, Lx - margin - R)
            ok = True
            for c2, R2, _h2 in bumps:
                if abs(center - c2) < gap_factor * (R + R2):
                    ok = False
                    break
            if ok:
                bumps.append((center, R, height))
                break
        else:
            raise RuntimeError(
                f"could not place bump {len(bumps) + 1}/{n_bumps} without overlap -- "
                "reduce --n-bumps, shrink --R-range, or shrink --gap-factor"
            )
    bumps.sort(key=lambda b: b[0])
    return bumps


def open_uniform_knots(N, p):
    """Open-uniform knot vector for N elements, degree p, C^{p-1}."""
    interior = np.linspace(0.0, 1.0, N + 1)[1:-1]
    return np.concatenate([np.zeros(p + 1), interior, np.ones(p + 1)])


def greville_abscissae(U, p):
    n = len(U) - p - 1
    return np.array([np.mean(U[i + 1:i + p + 1]) for i in range(n)])


def build_surface(Lx, Ly, Nx, Ny, P, bumps):
    Ux = open_uniform_knots(Nx, P)
    Uy = open_uniform_knots(Ny, P)
    gx = greville_abscissae(Ux, P)
    gy = greville_abscissae(Uy, P)

    x = Lx * gx
    bottom_y = _bump_field(x, bumps)
    top_y = np.full_like(x, Ly)

    if np.any(bottom_y >= Ly):
        raise ValueError(
            f"floor bump reaches/exceeds Ly={Ly:.3e} (max={bottom_y.max():.3e}) "
            "-- shrink --h-range"
        )

    nx, ny = len(gx), len(gy)
    ctrl = np.zeros((nx, ny, 3))
    for i in range(nx):
        ctrl[i, :, 0] = x[i]
        ctrl[i, :, 1] = bottom_y[i] + gy * (top_y[i] - bottom_y[i])

    return NURBS([Ux, Uy], control=ctrl)


def plot_surface(surf, bumps, fname):
    nu, nv = 200, 80
    u = np.linspace(surf.knots[0][0], surf.knots[0][-1], nu)
    v = np.linspace(surf.knots[1][0], surf.knots[1][-1], nv)
    pts = surf(u, v)

    fig, ax = plt.subplots(figsize=(8, 3.5))
    ax.plot(pts[:, 0, 0], pts[:, 0, 1], 'k-', lw=1.5)
    ax.plot(pts[:, -1, 0], pts[:, -1, 1], 'k-', lw=1.5)
    ax.plot(pts[0, :, 0], pts[0, :, 1], 'k-', lw=1.5)
    ax.plot(pts[-1, :, 0], pts[-1, :, 1], 'k-', lw=1.5)

    for ub in surf.breaks(0):
        line = surf([ub], v)[0]
        ax.plot(line[:, 0], line[:, 1], 'b-', lw=0.3)
    for vb in surf.breaks(1):
        line = surf(u, [vb])[:, 0]
        ax.plot(line[:, 0], line[:, 1], 'b-', lw=0.3)

    ctrl = surf.points
    ax.plot(ctrl[..., 0].ravel(), ctrl[..., 1].ravel(), 'r.', ms=2)

    for center, R, height in bumps:
        ax.plot(center, height, 'g+', ms=8, mew=1.5)

    ax.set_aspect('equal')
    ax.set_title(f'Random fracture-wall channel ({len(bumps)} bumps)')
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    fig.tight_layout()
    fig.savefig(fname, dpi=150)
    print(f"wrote {fname}")


def write_vtk_legacy(surf, fname, n_per_elem=4):
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
    pts = surf(u, v)
    nu, nv = pts.shape[0], pts.shape[1]
    with open(fname, "w") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("random_bump_channel\n")
        f.write("ASCII\n")
        f.write("DATASET STRUCTURED_GRID\n")
        f.write(f"DIMENSIONS {nu} {nv} 1\n")
        f.write(f"POINTS {nu * nv} float\n")
        for j in range(nv):
            for i in range(nu):
                x, y, _ = pts[i, j]
                f.write(f"{x:.8e} {y:.8e} 0.0\n")
    print(f"wrote {fname}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--n-bumps", type=int, default=8)
    p.add_argument("--Lx", type=float, default=1.0e-4)
    p.add_argument("--Ly", type=float, default=4.0e-5)
    p.add_argument("--Nx", type=int, default=300)
    p.add_argument("--Ny", type=int, default=150)
    p.add_argument("--P", type=int, default=2)
    p.add_argument("--R-range", type=float, nargs=2, default=(0.15e-5, 0.6e-5),
                    help="min/max bump half-width [m]")
    p.add_argument("--h-range", type=float, nargs=2, default=(0.1e-5, 0.35e-5),
                    help="min/max bump height [m]")
    p.add_argument("--gap-factor", type=float, default=1.2,
                    help="min center spacing, as a multiple of (R_i+R_j)")
    p.add_argument("--margin", type=float, default=0.3e-5,
                    help="keep-out zone at the left/right domain edges [m]")
    p.add_argument("--out", default="random_bump_channel.dat")
    p.add_argument("--vtk", default="random_bump_channel.vtk")
    p.add_argument("--plot", default="random_bump_channel.png")
    args = p.parse_args()

    rng = np.random.default_rng(args.seed)
    bumps = random_bumps(rng, args.Lx, args.n_bumps, tuple(args.R_range),
                          tuple(args.h_range), args.gap_factor, args.margin)
    print("bumps (center, R, height):")
    for b in bumps:
        print(f"  {b[0]:.4e}  {b[1]:.4e}  {b[2]:.4e}")

    surf = build_surface(args.Lx, args.Ly, args.Nx, args.Ny, args.P, bumps)
    print("degree:", surf.degree)
    print("shape (control points):", surf.shape[:2])

    plot_surface(surf, bumps, args.plot)
    write_vtk_legacy(surf, args.vtk)
    PetIGA().write(args.out, surf)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
