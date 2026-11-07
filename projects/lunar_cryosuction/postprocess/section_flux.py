#!/usr/bin/env python3
"""Vapor and heat flux through any line in the domain (simulation plan tool
T3; measurement M9), for every snapshot.

    J = integral over the line of  -D_v(T) phi_aef grad(rho_v) . n  ds   [kg/(m s)]
    Q = integral over the line of  -k(phi)       grad(T)     . n  ds     [W/m]

with the solver's own laws (pplib.dif_vap, thermal_cond, vapor_air_fraction,
the same as BoundaryFluxMonitor in src/monitoring.c): UNSCALED D_v and k, and
phi_aef = max(1 - phi, air_lim). n is the unit normal of the line pointing to
the RIGHT of the direction of travel from P0 to P1. Positive J is flow
across the line in the +n direction:

    line from (x, 0) to (x, Ly)  ->  n = +x_hat  (J > 0: vapor moving deeper)
    line from (0, y) to (Lx, y)  ->  n = -y_hat

The field values come from the evaluated spline (pplib.iter_snapshots), and
the gradients from second-order differences on that grid, mapped to physical
space (pplib.grid_gradient), so a mapped mesh is handled. The line is sampled
at --n points (default: about 4 per grid spacing) with the midpoint rule; parts
of the line outside the domain are dropped and reported.

CHECKS BUILT IN

On a domain face the result can be compared with bflux.dat, which the solver
computes at its quadrature points: --check-face x0 integrates along that face
(a hair inside it) and prints both numbers. On the stage-0c channel they agree
to about 1e-10 relative (postprocess/tests/test_new_tools.py).

Usage:
    python3 postprocess/section_flux.py --dir <run> --x 2e-4          # vertical line at x
    python3 postprocess/section_flux.py --dir <run> --line X0 Y0 X1 Y1
    python3 postprocess/section_flux.py --dir <run> --check-face x0
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pplib  # noqa: E402


def line_flux(snap, p0, p1, n_pts=None, grads=None, sampler=None):
    """(J, Q, length_inside, length_total) through the segment p0 -> p1."""
    X, Y = snap.X, snap.Y
    p0 = np.asarray(p0, dtype=float)
    p1 = np.asarray(p1, dtype=float)
    d = p1 - p0
    L = float(np.hypot(*d))
    if L == 0:
        raise ValueError("the line has zero length")
    if n_pts is None:
        h = min(np.min(np.hypot(np.diff(X, axis=1), np.diff(Y, axis=1))),
                np.min(np.hypot(np.diff(X, axis=0), np.diff(Y, axis=0))))
        n_pts = int(min(max(50, 4 * L / max(h, 1e-300)), 200000))
    s = (np.arange(n_pts) + 0.5) / n_pts          # midpoint rule
    px, py = p0[0] + s * d[0], p0[1] + s * d[1]
    nx, ny = d[1] / L, -d[0] / L                   # right-hand normal
    if grads is None:
        grads = (pplib.grid_gradient(snap.rhov, X, Y), pplib.grid_gradient(snap.T, X, Y))
    (rx, ry), (tx, ty) = grads
    S = sampler or pplib.GridSampler(X, Y)
    phi = S(snap.phi, px, py)
    T = S(snap.T, px, py)
    grv = S(rx, px, py) * nx + S(ry, px, py) * ny
    gT = S(tx, px, py) * nx + S(ty, px, py) * ny
    inside = np.isfinite(phi) & np.isfinite(grv) & np.isfinite(gT)
    ds = L / n_pts
    j = -pplib.dif_vap(T[inside]) * pplib.vapor_air_fraction(phi[inside]) * grv[inside]
    q = -pplib.thermal_cond(np.clip(phi[inside], 0.0, 1.0)) * gT[inside]
    return float(j.sum() * ds), float(q.sum() * ds), float(inside.sum() * ds), L


def face_line(snap, face, inset_frac=1e-6):
    """End points of a straight domain face, moved inside by inset_frac of the
    domain size so the sampler sees the face's own cells."""
    X, Y = snap.X, snap.Y
    x0, x1, y0, y1 = X.min(), X.max(), Y.min(), Y.max()
    dx, dy = inset_frac * (x1 - x0), inset_frac * (y1 - y0)
    # Orientation chosen so the right-hand normal points OUT of the domain,
    # the bflux.dat sign convention.
    return {"x0": ((x0 + dx, y1), (x0 + dx, y0)),
            "x1": ((x1 - dx, y0), (x1 - dx, y1)),
            "y0": ((x0, y0 + dy), (x1, y0 + dy)),
            "y1": ((x1, y1 - dy), (x0, y1 - dy))}[face]


def run(run_dirs, lines, source="auto", every=1, n_pts=None):
    """Rows (step, t, line index, J, Q, inside, total) over all snapshots."""
    rows = []
    sampler = None
    for snap in pplib.iter_snapshots(run_dirs, source=source, every=every):
        if sampler is None:
            sampler = pplib.GridSampler(snap.X, snap.Y)
        grads = (pplib.grid_gradient(snap.rhov, snap.X, snap.Y), pplib.grid_gradient(snap.T, snap.X, snap.Y))
        for k, ln in enumerate(lines):
            p0, p1 = ln(snap) if callable(ln) else ln
            J, Q, li, lt = line_flux(snap, p0, p1, n_pts, grads, sampler)
            rows.append((snap.step, snap.t, k, J, Q, li, lt))
    return rows


def bflux_at(run_dir, face, steps):
    """J and Q of a face from bflux.dat at the given steps (NaN if absent)."""
    path = os.path.join(run_dir, "bflux.dat")
    if not os.path.isfile(path):
        return {}
    names = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("# t step"):
                names = line[1:].split()
                break
    a = np.genfromtxt(path, comments="#", ndmin=2)
    if f"J_vap_{face}" not in names:
        return {}
    cj, cq = names.index(f"J_vap_{face}"), names.index(f"Q_heat_{face}")
    by = {int(r[1]): (r[cj], r[cq]) for r in a}
    return {s: by.get(s, (np.nan, np.nan)) for s in steps}


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dir", action="append", required=True)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--x", type=float, action="append", help="vertical line across the pore at this x [m]")
    g.add_argument("--line", type=float, nargs=4, action="append", metavar=("X0", "Y0", "X1", "Y1"))
    g.add_argument("--check-face", choices=("x0", "x1", "y0", "y1"))
    p.add_argument("--n", type=int, default=None, help="sample points on each line")
    p.add_argument("--every", type=int, default=1)
    p.add_argument("--source", choices=("auto", "nurbs", "vts"), default="auto")
    p.add_argument("--out", default=None)
    a = p.parse_args(argv)

    if a.check_face:
        lines = [lambda s, f=a.check_face: face_line(s, f)]
        labels = [f"face {a.check_face} (outward normal)"]
    elif a.x:
        lines = [lambda s, x=x: ((x, s.Y.min()), (x, s.Y.max())) for x in a.x]
        labels = [f"x = {x:.4e} m (normal +x)" for x in a.x]
    else:
        lines = [((l[0], l[1]), (l[2], l[3])) for l in a.line]
        labels = [f"({l[0]:.3e},{l[1]:.3e}) -> ({l[2]:.3e},{l[3]:.3e})" for l in a.line]

    rows = run(a.dir, lines, a.source, a.every, a.n)
    out = a.out or a.dir[0]
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, "section_flux.csv")
    with open(path, "w") as fh:
        fh.write("# " + " | ".join(f"line {k}: {lab}" for k, lab in enumerate(labels)) + "\n")
        fh.write("step,t_s,line,J_vap_kg_per_m_s,Q_heat_W_per_m,length_inside_m,length_m\n")
        for r in rows:
            fh.write(f"{r[0]},{r[1]:.10e},{r[2]},{r[3]:.10e},{r[4]:.10e},{r[5]:.6e},{r[6]:.6e}\n")
    print(f"section_flux: {len({r[0] for r in rows})} snapshots, {len(lines)} line(s); wrote {path}")
    short = [r for r in rows if r[5] < 0.999 * r[6]]
    if short:
        print(f"  note: {len(short)} line samples run partly outside the domain "
              f"(inside {min(r[5] / r[6] for r in short):.3f} of the length at worst)")
    if a.check_face:
        ref = bflux_at(a.dir[0], a.check_face, [r[0] for r in rows])
        print(f"  {'step':>6} {'t [s]':>11} {'J here':>12} {'J bflux':>12} {'rel':>9} "
              f"{'Q here':>12} {'Q bflux':>12}")
        for r in rows:
            jb, qb = ref.get(r[0], (np.nan, np.nan))
            rel = abs(r[3] - jb) / abs(jb) if jb else np.nan
            print(f"  {r[0]:6d} {r[1]:11.4e} {r[3]:12.4e} {jb:12.4e} {rel:9.2e} {r[4]:12.4e} {qb:12.4e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
