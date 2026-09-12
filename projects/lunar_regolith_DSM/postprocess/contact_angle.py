#!/usr/bin/env python3
"""Measure the equilibrium contact angle where ice meets a regolith wall.

WHAT IS BEING TESTED

The solver is given three surface energies and derives cos(theta) from Young's
equation, cos(theta) = (gamma_as - gamma_is)/gamma_ia. It then enforces that
angle only LOCALLY, through a wall free-energy term whose natural boundary
condition is dphi/dn = cos(theta)*phi(1-phi)/eps. Nothing in the solver makes
the MACROSCOPIC shape of the ice come out right; that has to emerge from the
relaxation. This script measures the macroscopic angle and compares it to
Young's, which is the actual validation.

HOW

1. Evaluate the NURBS basis directly rather than reading the .vts control net.
   The .vts files carry p=2 B-spline CONTROL COEFFICIENTS, not field values;
   reading them as values costs ~1e-9 m of interface position that drifts in
   phase as the interface crosses the grid (wedge_gt_velocity.py:298-333).

2. Extract the phi=0.5 contour row by row (pplib.contour_points).

3. DISCARD contour points within --exclude-eps decay lengths of either wall.
   The diffuse profile bends within a few eps of the contact line; the
   macroscopic angle lives in the arc, not at the contact line. Fitting through
   the bent region is the single easiest way to get a wrong answer here.

4. Least-squares circle through the remaining arc (pplib.circle_radius), then
   extrapolate that circle to the wall.

5. At the intersection P, take the radial unit vector r, and decide which way
   is "into the ice" by SAMPLING phi a short distance either side of P along r
   -- m points toward increasing phi, by definition. Then

       cos(theta) = m . n        (n = outward wall normal)

   which is precisely the quantity the boundary condition enforces, measured
   from the macroscopic shape instead of from the BC. Sampling phi rather than
   assuming concave/convex keeps this correct at theta = 0 and 180, where a
   sign convention based on the arc's curvature degenerates.

A confined bridge gives four independent estimates (two menisci x two walls);
their spread is the error bar. Also reported, as a cross-check on the circle
fit only, is the parallel-plate result cos(theta) = d_perp/R_arc.

Usage:
    python3 postprocess/contact_angle.py --dir <run_dir> [--save fig.png]
"""
import argparse
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pplib import (                                    # noqa: E402
    auto_time_unit,
    contour_points,
    in_time_unit,
    opt_float,
    read_opts,
    step_times,
)

CSV_NAME = "contact_angle.csv"


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------
def read_snapshots(run_dir, nu, nv):
    """Yield (step, x1d, y1d, phi[ix, iy]) on the exact NURBS field.

    Assumes a rectangular patch: x depends only on the u index and y only on
    the v index. That holds for the plain Cartesian channel this test uses, and
    is checked rather than assumed -- a curvilinear mesh (a wedge, a bumpy
    floor) would silently give wrong coordinates otherwise.
    """
    from igakit.io import PetIGA
    io = PetIGA()
    nrb = io.read(os.path.join(run_dir, "igasol.dat"))
    u = np.linspace(0.0, 1.0, nu)
    v = np.linspace(0.0, 1.0, nv)
    files = sorted(f for f in os.listdir(run_dir)
                   if re.fullmatch(r"sol_\d+\.dat", f))
    if not files:
        raise SystemExit(f"no sol_*.dat in {run_dir}")

    checked = False
    for f in files:
        C, F = nrb(u, v, fields=io.read_vec(os.path.join(run_dir, f), nrb))
        if not checked:
            dx = np.ptp(C[:, :, 0], axis=1).max()
            dy = np.ptp(C[:, :, 1], axis=0).max()
            span = max(np.ptp(C[:, :, 0]), np.ptp(C[:, :, 1]))
            if dx > 1e-9 * span or dy > 1e-9 * span:
                raise SystemExit(
                    "the patch is not a rectangle (x varies along v by "
                    f"{dx:.3e} m, y varies along u by {dy:.3e} m). This script "
                    "assumes the flat channel geometry; a curvilinear mesh "
                    "needs the contour mapped through the geometry instead.")
            checked = True
        yield (int(re.search(r"sol_(\d+)", f).group(1)),
               C[:, 0, 0].copy(), C[0, :, 1].copy(), F[:, :, 0].copy())


def bilinear(x1d, y1d, phi, px, py):
    """phi at (px, py) by bilinear interpolation on the separable grid."""
    i = np.clip(np.searchsorted(x1d, px) - 1, 0, x1d.size - 2)
    j = np.clip(np.searchsorted(y1d, py) - 1, 0, y1d.size - 2)
    tx = (px - x1d[i]) / (x1d[i + 1] - x1d[i])
    ty = (py - y1d[j]) / (y1d[j + 1] - y1d[j])
    return ((1 - tx) * (1 - ty) * phi[i, j] + tx * (1 - ty) * phi[i + 1, j]
            + (1 - tx) * ty * phi[i, j + 1] + tx * ty * phi[i + 1, j + 1])


# ---------------------------------------------------------------------------
# The measurement
# ---------------------------------------------------------------------------
def fit_implicit_circle(px, py):
    """Fit A(x^2+y^2) + Bx + Cy + D = 0 to the contour points.

    Returns (A, B, C, D) normalised to unit length, or None.

    Why not the affine fit in pplib.circle_radius (2ax + 2by + c = x^2 + y^2):
    that one cannot represent a STRAIGHT interface at all. It is the a -> inf
    limit, so at theta = 90 -- where the meniscus in a channel is exactly flat
    and R is infinite -- it returns numerical noise. Measured on a synthetic
    90-degree bridge it gave 78 +/- 25 degrees.

    The implicit form has no such limit: A -> 0 IS the line, and the unit-norm
    constraint keeps the solution well posed there. The normal direction comes
    straight from the gradient,
        grad F = (2Ax + B, 2Ay + C),
    which stays finite and correct whether the arc is curved or flat.

    Coordinates are centred and scaled before the fit, purely for conditioning,
    and the coefficients are mapped back afterwards.
    """
    if px.size < 5:
        return None
    x0, y0 = px.mean(), py.mean()
    sc = max(np.ptp(px), np.ptp(py))
    if not np.isfinite(sc) or sc <= 0.0:
        return None
    u, v = (px - x0) / sc, (py - y0) / sc

    M = np.column_stack((u * u + v * v, u, v, np.ones(u.size)))
    try:
        _, _, Vt = np.linalg.svd(M, full_matrices=False)
    except np.linalg.LinAlgError:
        return None
    a, b, c, d = Vt[-1]                       # smallest singular value

    # Undo (x - x0)/sc: A(u^2+v^2) + Bu + Cv + D  ->  coefficients in x, y.
    A = a / (sc * sc)
    B = b / sc - 2.0 * a * x0 / (sc * sc)
    C = c / sc - 2.0 * a * y0 / (sc * sc)
    D = d - b * x0 / sc - c * y0 / sc + a * (x0 * x0 + y0 * y0) / (sc * sc)
    n = np.hypot(np.hypot(A, B), np.hypot(C, D))
    if not np.isfinite(n) or n == 0.0:
        return None
    return A / n, B / n, C / n, D / n


def circle_from_implicit(coef):
    """(xc, yc, R) from (A, B, C, D), or (nan, nan, inf) for a straight arc."""
    A, B, C, D = coef
    if abs(A) < 1e-12:
        return float("nan"), float("nan"), float("inf")
    xc, yc = -B / (2.0 * A), -C / (2.0 * A)
    disc = xc * xc + yc * yc - D / A
    return xc, yc, (float(np.sqrt(disc)) if disc > 0.0 else float("nan"))


def angle_at_wall(coef, y_wall, n_out, arc_x, x1d, y1d, phi, eps):
    """Contact angle where the fitted arc meets the wall y = y_wall.

    Returns (theta_deg, x_contact) or None when the arc does not reach the wall
    -- which happens legitimately for a drop that has pulled away from it.
    """
    A, B, C, D = coef
    k = A * y_wall * y_wall + C * y_wall + D          # F(x, y_wall) = A x^2 + B x + k

    if abs(A) < 1e-12:                                 # straight arc
        if abs(B) < 1e-30:
            return None
        roots = np.array([-k / B])
    else:
        disc = B * B - 4.0 * A * k
        if disc < 0.0:
            return None
        sq = np.sqrt(disc)
        roots = np.array([(-B - sq) / (2.0 * A), (-B + sq) / (2.0 * A)])

    # Take the intersection the measured arc actually runs toward.
    px = float(roots[np.argmin(np.abs(roots - np.mean(arc_x)))])

    g = np.array([2.0 * A * px + B, 2.0 * A * y_wall + C])
    ng = np.hypot(*g)
    if not np.isfinite(ng) or ng == 0.0:
        return None
    r = g / ng                                         # unit normal to the arc

    # Which way is into the ice? phi increases that way, by definition. Probing
    # the field rather than assuming concave/convex is what keeps this correct
    # at theta = 0 and 180, where a curvature-based sign convention degenerates.
    d = 2.0 * eps
    probe = np.array([px, y_wall]) + np.outer([+1.0, -1.0], r) * d
    probe[:, 0] = np.clip(probe[:, 0], x1d[0], x1d[-1])
    probe[:, 1] = np.clip(probe[:, 1], y1d[0], y1d[-1])
    vals = bilinear(x1d, y1d, phi, probe[:, 0], probe[:, 1])
    m = r if vals[0] >= vals[1] else -r

    cos_t = float(np.clip(np.dot(m, n_out), -1.0, 1.0))
    return float(np.degrees(np.arccos(cos_t))), px


def measure(x1d, y1d, phi, eps, walls, exclude_eps, level=0.5):
    """All contact-angle estimates for one snapshot.

    `walls` is a list of (y_wall, outward_normal) pairs.
    Returns a list of dicts, one per (meniscus, wall) pair.
    """
    ny, nx = y1d.size, x1d.size
    Y = np.broadcast_to(y1d[:, None], (ny, nx))
    xl, yl, xr, yr = contour_points(x1d, Y, phi.T, level=level)
    if xl.size == 0:
        return []

    out = []
    for name, ax, ay in (("left", xl, yl), ("right", xr, yr)):
        for y_wall, n_out in walls:
            keep = np.abs(ay - y_wall) > exclude_eps * eps
            # Stay on this side of the channel: for a bridge between two walls
            # the far wall's bent region is excluded by its own filter, but the
            # arc must still be long enough to pin a circle.
            for yw2, _ in walls:
                if yw2 != y_wall:
                    keep &= np.abs(ay - yw2) > exclude_eps * eps
            if keep.sum() < 8:
                continue
            coef = fit_implicit_circle(ax[keep], ay[keep])
            if coef is None:
                continue
            res = angle_at_wall(coef, y_wall, n_out, ax[keep], x1d, y1d, phi, eps)
            if res is None:
                continue
            theta, px = res
            xc, yc, R = circle_from_implicit(coef)
            # Parallel-plate cross-check: cos(theta) = d_perp/R_arc. Validates
            # the circle fit, not the physics -- it assumes the very shape the
            # fit produced.
            if np.isfinite(R) and np.isfinite(yc) and R > 0.0:
                cos_pp = np.clip(abs(yc - y_wall) / R, -1.0, 1.0)
                theta_pp = float(np.degrees(np.arccos(cos_pp)))
                if theta > 90.0:
                    theta_pp = 180.0 - theta_pp
            else:
                theta_pp = 90.0        # straight arc: the flat-meniscus limit
            out.append(dict(meniscus=name, y_wall=y_wall, theta=theta,
                            theta_plate=float(theta_pp), R_arc=float(R),
                            xc=float(xc), yc=float(yc), x_contact=px,
                            n_points=int(keep.sum())))
    return out


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True, help="run directory")
    ap.add_argument("--save", default=None, help="figure path (png)")
    ap.add_argument("--nu", type=int, default=601, help="NURBS samples along x")
    ap.add_argument("--nv", type=int, default=301, help="NURBS samples along y")
    ap.add_argument("--exclude-eps", type=float, default=5.0,
                    help="discard contour points within this many eps of a wall")
    ap.add_argument("--level", type=float, default=0.5)
    args = ap.parse_args()

    run = args.dir
    opts = read_opts(run)
    eps = opt_float(opts, "-eps")
    Ly = opt_float(opts, "-Ly")
    if eps is None or Ly is None:
        raise SystemExit("could not read -eps and -Ly from the staged .opts")

    faces = str(opts.get("-wall_faces", "") or "")
    walls = []
    if "y0" in faces:
        walls.append((0.0, np.array([0.0, -1.0])))
    if "y1" in faces:
        walls.append((Ly, np.array([0.0, +1.0])))
    if not walls:
        raise SystemExit(
            "-wall_faces names no y face in this run, so there is no regolith "
            "wall to measure a contact angle against.")

    # What the solver was told, for the comparison column.
    g_ia = opt_float(opts, "-gamma_ia")
    g_is = opt_float(opts, "-gamma_is", 0.0)
    g_as = opt_float(opts, "-gamma_as", 0.0)
    theta_cli = opt_float(opts, "-contact_angle_deg")
    if theta_cli is not None:
        theta_young = theta_cli
        source = "-contact_angle_deg (Young bypassed)"
    else:
        if g_ia is None:
            g_ia = opt_float(opts, "-Sigma_i", 0.109)
        theta_young = float(np.degrees(np.arccos(
            np.clip((g_as - g_is) / g_ia, -1.0, 1.0))))
        source = "Young from gamma_is/gamma_as"

    times = step_times(run)
    rows = []
    for step, x1d, y1d, phi in read_snapshots(run, args.nu, args.nv):
        for r in measure(x1d, y1d, phi, eps, walls, args.exclude_eps, args.level):
            r.update(step=step, t=times.get(step, np.nan))
            rows.append(r)
    if not rows:
        raise SystemExit("no contour could be measured in any snapshot")

    steps = sorted({r["step"] for r in rows})
    csv = os.path.join(run, CSV_NAME)
    with open(csv, "w") as fh:
        fh.write("# gamma_ia,gamma_is,gamma_as = %.6g,%.6g,%.6g   "
                 "theta_young = %.4f deg   (%s)\n"
                 % (g_ia, g_is, g_as, theta_young, source))
        fh.write("step,time_s,meniscus,y_wall,theta_deg,theta_plate_deg,"
                 "R_arc_m,x_contact_m,n_points,theta_young_deg,error_deg\n")
        for r in rows:
            fh.write("%d,%.6e,%s,%.6e,%.4f,%.4f,%.6e,%.6e,%d,%.4f,%.4f\n"
                     % (r["step"], r["t"], r["meniscus"], r["y_wall"],
                        r["theta"], r["theta_plate"], r["R_arc"],
                        r["x_contact"], r["n_points"], theta_young,
                        r["theta"] - theta_young))

    last = [r for r in rows if r["step"] == steps[-1]]
    th = np.array([r["theta"] for r in last])
    print(f"\n  CONTACT ANGLE  ({os.path.basename(run.rstrip('/'))})")
    print(f"    gamma_ia = {g_ia:.6g}  gamma_is = {g_is:.6g}  "
          f"gamma_as = {g_as:.6g}   J/m^2")
    print(f"    theta_Young    = {theta_young:8.3f} deg   [{source}]")
    print(f"    theta_measured = {th.mean():8.3f} +/- {th.std():.3f} deg "
          f"({th.size} estimates, final snapshot)")
    print(f"    error          = {th.mean() - theta_young:+8.3f} deg")
    print(f"    -> {csv}")

    if args.save:
        unit = auto_time_unit(max(r["t"] for r in rows if np.isfinite(r["t"]))
                              if any(np.isfinite(r["t"]) for r in rows) else 1.0)
        fig, ax = plt.subplots(figsize=(7.0, 4.2), constrained_layout=True)
        for name in ("left", "right"):
            for y_wall, _ in walls:
                sel = [r for r in rows
                       if r["meniscus"] == name and r["y_wall"] == y_wall]
                if not sel:
                    continue
                t = in_time_unit([r["t"] for r in sel], unit)
                ax.plot(t, [r["theta"] for r in sel], marker="o", ms=3, lw=1.2,
                        label=f"{name} meniscus, wall y={y_wall:.3g}")
        ax.axhline(theta_young, color="k", ls="--", lw=1.4,
                   label=f"Young  {theta_young:.1f}°")
        ax.set_xlabel(f"time [{unit}]")
        ax.set_ylabel("contact angle [deg]")
        ax.set_title("Measured contact angle vs Young's equation")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
        fig.savefig(args.save, dpi=150)
        print(f"    -> {args.save}")


if __name__ == "__main__":
    main()
