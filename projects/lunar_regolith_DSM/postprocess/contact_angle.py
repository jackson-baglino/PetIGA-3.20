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
    """Yield (step, X, Y, phi) on the exact NURBS field, all shaped (nu, nv).

    X and Y are the PHYSICAL coordinates of each parametric sample, so this
    works for a curvilinear patch (the wedge) as well as a rectangle. The one
    structural assumption, checked below, is that the patch is RULED in the
    sense that x depends on u alone -- true for every geometry
    build_geometry_*.py produces, since they warp y between two wall curves and
    leave x untouched. It is what lets the phi = 0.5 contour be found by
    scanning along u within each v row, and what makes the inverse map in
    sample_phi() exact.
    """
    from igakit.io import PetIGA
    io = PetIGA()
    nrb = io.read(os.path.join(run_dir, "igasol.dat"))

    def prange(k, deg):
        return float(k[deg]), float(k[-deg - 1])

    (u0, u1) = prange(nrb.knots[0], nrb.degree[0])
    (v0, v1) = prange(nrb.knots[1], nrb.degree[1])
    u = np.linspace(u0, u1, nu)
    v = np.linspace(v0, v1, nv)

    files = sorted(f for f in os.listdir(run_dir)
                   if re.fullmatch(r"sol_\d+\.dat", f))
    if not files:
        raise SystemExit(f"no sol_*.dat in {run_dir}")

    checked = False
    for f in files:
        C, F = nrb(u, v, fields=io.read_vec(os.path.join(run_dir, f), nrb))
        if not checked:
            span = max(np.ptp(C[:, :, 0]), np.ptp(C[:, :, 1]))
            dx = np.ptp(C[:, :, 0], axis=1).max()
            if dx > 1e-9 * span:
                raise SystemExit(
                    f"x varies by {dx:.3e} m along a column, so the patch is not "
                    "ruled in x and the row-scan contour extraction does not "
                    "apply. A genuinely 2D-warped mesh needs the contour taken "
                    "in parametric space and mapped through the geometry.")
            checked = True
        yield (int(re.search(r"sol_(\d+)", f).group(1)),
               C[:, :, 0].copy(), C[:, :, 1].copy(), F[:, :, 0].copy())


def wall_y(px, w):
    """y of the wall curve at x = px.  w = (y0, slope)."""
    return w[0] + w[1] * px


def wall_normal(w, side):
    """Outward unit normal of an affine wall.  side 0 = bottom, 1 = top."""
    s = w[1]
    n = np.array([s, -1.0]) if side == 0 else np.array([-s, 1.0])
    return n / np.hypot(*n)


def wall_distance(px, py, w):
    """Perpendicular distance from a point to the wall line."""
    return abs(py - wall_y(px, w)) / np.hypot(1.0, w[1])


def sample_phi(X, Y, phi, px, py, wb, wt):
    """phi at physical (px, py), by inverting the ruled map exactly.

    x depends only on u, so u comes from interpolating the column positions;
    the patch is ruled between the two wall curves, so
    v = (y - y_bot(x)) / (y_top(x) - y_bot(x)) is exact. Then bilinear in index
    space. For flat walls this reduces to the separable rectangular case.
    """
    xs = X[:, 0]
    iu = np.interp(px, xs, np.arange(xs.size))
    yb, yt = wall_y(px, wb), wall_y(px, wt)
    iv = (py - yb) / (yt - yb) * (Y.shape[1] - 1)
    iu = np.clip(iu, 0, X.shape[0] - 1.001)
    iv = np.clip(iv, 0, Y.shape[1] - 1.001)
    i0, j0 = int(iu), int(iv)
    tu, tv = iu - i0, iv - j0
    return ((1 - tu) * (1 - tv) * phi[i0, j0] + tu * (1 - tv) * phi[i0 + 1, j0]
            + (1 - tu) * tv * phi[i0, j0 + 1] + tu * tv * phi[i0 + 1, j0 + 1])


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


def angle_at_wall(coef, w, side, arc_x, X, Y, phi, eps, wb, wt, ice_dx):
    """Contact angle where the fitted arc meets the wall line y = y0 + s*x.

    Generalises the flat-wall case: substituting the wall line into
    A(x^2+y^2) + Bx + Cy + D = 0 gives a quadratic in x,

        A(1+s^2) x^2 + (2 A s y0 + B + C s) x + (A y0^2 + C y0 + D) = 0,

    which reduces to the old A x^2 + B x + (A y_w^2 + C y_w + D) when s = 0.

    Returns (theta_deg, x_contact, y_contact) or None if the arc misses the wall.
    """
    A, B, C, D = coef
    y0, sl = w
    a2 = A * (1.0 + sl * sl)
    a1 = 2.0 * A * sl * y0 + B + C * sl
    a0 = A * y0 * y0 + C * y0 + D

    if abs(a2) < 1e-12:                     # straight arc
        if abs(a1) < 1e-30:
            return None
        roots = np.array([-a0 / a1])
    else:
        disc = a1 * a1 - 4.0 * a2 * a0
        if disc < 0.0:
            return None
        sq = np.sqrt(disc)
        roots = np.array([(-a1 - sq) / (2.0 * a2), (-a1 + sq) / (2.0 * a2)])

    px = float(roots[np.argmin(np.abs(roots - np.mean(arc_x)))])
    py = wall_y(px, w)

    g = np.array([2.0 * A * px + B, 2.0 * A * py + C])
    ng = np.hypot(*g)
    if not np.isfinite(ng) or ng == 0.0:
        return None
    r = g / ng                               # unit normal to the arc

    # Orient toward the ice. Probing phi is exact here thanks to the ruled
    # inverse map, and unlike a curvature-based convention it stays correct at
    # theta = 0 and 180. Fall back on the scan direction if both probes land
    # outside the patch.
    d = 2.0 * eps
    best, m = None, None
    for sgn in (+1.0, -1.0):
        q = np.array([px, py]) + sgn * r * d
        if wall_distance(q[0], q[1], wb) < 0 or not np.isfinite(q).all():
            continue
        val = sample_phi(X, Y, phi, q[0], q[1], wb, wt)
        if best is None or val > best:
            best, m = val, sgn * r
    if m is None:
        m = r if np.dot(r, [ice_dx, 0.0]) > 0 else -r

    n_out = wall_normal(w, side)
    cos_t = float(np.clip(np.dot(m, n_out), -1.0, 1.0))
    return float(np.degrees(np.arccos(cos_t))), px, py


# ---------------------------------------------------------------------------
# The measurement
# ---------------------------------------------------------------------------
def measure(X, Y, phi, eps, walls, bounds, exclude_eps, level=0.5):
    """All contact-angle estimates for one snapshot.

    `walls` is a list of (w, side) with w = (y0, slope) and side 0 = bottom.
    Contour points are taken row by row along u, in physical coordinates, and
    the wall exclusion uses PERPENDICULAR distance to the wall line -- on a
    wedge a vertical offset would under-cut the near wall and over-cut the far
    one.
    """
    # The ruled inverse map needs BOTH bounding curves, even when only one of
    # them carries a contact angle.
    wb, wt = bounds
    xs = X[:, 0]
    xl, yl, xr, yr = contour_points(xs, Y.T, phi.T, level=level)
    if xl.size == 0:
        return []

    # Ice lies to the RIGHT of the left meniscus and to the left of the right
    # one; used only as a fallback if both phi probes fall outside the patch.
    ice_left, ice_right = +1.0, -1.0
    out = []
    for name, ax, ay, idx in (("left", xl, yl, ice_left),
                              ("right", xr, yr, ice_right)):
        for w, side in walls:
            keep = np.ones(ax.shape, dtype=bool)
            for w2 in (wb, wt):
                keep &= wall_distance(ax, ay, w2) > exclude_eps * eps
            if keep.sum() < 8:
                continue
            coef = fit_implicit_circle(ax[keep], ay[keep])
            if coef is None:
                continue
            res = angle_at_wall(coef, w, side, ax[keep], X, Y, phi, eps,
                                wb, wt, idx)
            if res is None:
                continue
            theta, px, py = res
            xc, yc, R = circle_from_implicit(coef)
            if np.isfinite(R) and np.isfinite(yc) and R > 0.0:
                cos_pp = np.clip(wall_distance(xc, yc, w) / R, -1.0, 1.0)
                theta_pp = float(np.degrees(np.arccos(cos_pp)))
                if theta > 90.0:
                    theta_pp = 180.0 - theta_pp
            else:
                theta_pp = 90.0
            out.append(dict(meniscus=name, y_wall=wall_y(px, w), theta=theta,
                            theta_plate=float(theta_pp), R_arc=float(R),
                            xc=float(xc), yc=float(yc), x_contact=px,
                            n_points=int(keep.sum())))
    return out


# ---------------------------------------------------------------------------
def extrapolate(t, th):
    """Equilibrium angle from theta(t) = theta_inf + A*exp(-t/tau).

    Relaxation here is vapour-diffusion limited and slow: the pilot run
    measured tau = 14.4 days, so reaching the true equilibrium takes months of
    simulated time. It is also an extremely clean single exponential -- residual
    RMS 0.003 degrees over the last 40% of that run -- which makes the
    three-parameter fit a legitimate way to read the equilibrium off a run that
    stopped short, and reporting it guards the sweep against an under-estimated
    t_final.

    theta_inf is FREE. That is the whole point: it is compared to Young's
    prediction afterwards, never assumed.

    Returns (theta_inf, sigma, tau_seconds, rms, n) or None.
    """
    try:
        from scipy.optimize import curve_fit
    except ImportError:
        return None
    t = np.asarray(t, float)
    th = np.asarray(th, float)
    ok = np.isfinite(t) & np.isfinite(th)
    t, th = t[ok], th[ok]
    if t.size < 8 or t.max() <= 0:
        return None

    # Drop the early transient: the grain first shrinks under Gibbs-Thomson
    # until the sealed box's vapour saturates, and theta briefly moves the wrong
    # way. Fitting through that biases tau badly.
    sel = t > 0.4 * t.max()
    if sel.sum() < 6:
        return None

    def f(x, a, A, tau):
        return a + A * np.exp(-x / tau)

    span = max(th[sel].max() - th[sel].min(), 1e-6)
    try:
        p, cov = curve_fit(f, t[sel], th[sel],
                           p0=[th[sel][-1], span, 0.3 * t.max()], maxfev=40000)
    except Exception:
        return None
    if not np.all(np.isfinite(cov)) or p[2] <= 0:
        return None
    rms = float(np.sqrt(np.mean((th[sel] - f(t[sel], *p)) ** 2)))
    return float(p[0]), float(np.sqrt(cov[0, 0])), float(p[2]), rms, int(sel.sum())


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

    # Wall curves y = y0 + slope*x. The defaults (0, 0, Ly, 0) reproduce a flat
    # channel exactly, which is what initial_conditions.c documents; a wedge
    # sets the slopes and the same code handles it.
    wb = (opt_float(opts, "-wall_bot_y0", 0.0) or 0.0,
          opt_float(opts, "-wall_bot_slope", 0.0) or 0.0)
    wt = (opt_float(opts, "-wall_top_y0", None) if
          opt_float(opts, "-wall_top_y0", None) is not None else Ly,
          opt_float(opts, "-wall_top_slope", 0.0) or 0.0)
    bounds = (wb, wt)

    faces = str(opts.get("-wall_faces", "") or "")
    walls = []
    if "y0" in faces:
        walls.append((wb, 0))
    if "y1" in faces:
        walls.append((wt, 1))
    if not walls:
        raise SystemExit(
            "-wall_faces names no y face in this run, so there is no regolith "
            "wall to measure a contact angle against.")

    # What the solver was told, for the comparison column.
    # Runs written before 2026-09-18 name these -sigma_*; the solver has used
    # -gamma_* since. Read both, newest first, or every archived run silently
    # reports gamma_is = gamma_as = 0 -> theta_Young = 90 regardless of what it
    # was actually given.
    def _energy(new_key, old_key, default=None):
        v = opt_float(opts, new_key)
        if v is None:
            v = opt_float(opts, old_key)
        return default if v is None else v

    g_ia = _energy("-gamma_ia", "-sigma_ia")
    g_is = _energy("-gamma_is", "-sigma_is", 0.0)
    g_as = _energy("-gamma_as", "-sigma_as", 0.0)
    theta_cli = opt_float(opts, "-contact_angle_deg")
    if theta_cli is not None:
        theta_young = theta_cli
        source = "-contact_angle_deg (Young bypassed)"
    else:
        if g_ia is None:
            g_ia = _energy("-gamma_ia_bulk", "-sigma_ia_bulk", 0.109)
        theta_young = float(np.degrees(np.arccos(
            np.clip((g_as - g_is) / g_ia, -1.0, 1.0))))
        source = "Young from gamma_is/gamma_as"

    times = step_times(run)
    rows = []
    for step, X, Y, phi in read_snapshots(run, args.nu, args.nv):
        for r in measure(X, Y, phi, eps, walls, bounds, args.exclude_eps, args.level):
            r.update(step=step, t=times.get(step, np.nan))
            rows.append(r)
    if not rows:
        raise SystemExit("no contour could be measured in any snapshot")

    steps = sorted({r["step"] for r in rows})
    # Snapshot-averaged trajectory, for the extrapolation.
    traj_t, traj_th = [], []
    for st in steps:
        sel = [r for r in rows if r["step"] == st]
        traj_t.append(sel[0]["t"])
        traj_th.append(float(np.mean([r["theta"] for r in sel])))
    ext = extrapolate(traj_t, traj_th)
    csv = os.path.join(run, CSV_NAME)
    with open(csv, "w") as fh:
        fh.write("# gamma_ia,gamma_is,gamma_as = %.6g,%.6g,%.6g   "
                 "theta_young = %.4f deg   (%s)\n"
                 % (g_ia, g_is, g_as, theta_young, source))
        if ext is not None:
            fh.write("# theta_inf = %.4f +/- %.4f deg, tau = %.6e s, "
                     "fit_rms = %.5f, n = %d, error_inf = %+.4f deg\n"
                     % (ext[0], ext[1], ext[2], ext[3], ext[4],
                        ext[0] - theta_young))
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
    if ext is not None:
        a_inf, sig, tau, rms, n = ext
        print(f"    theta_inf      = {a_inf:8.3f} +/- {sig:.3f} deg   "
              f"[extrapolated, tau = {tau/86400.0:.2f} d, "
              f"fit RMS {rms:.4f} over {n} pts]")
        print(f"    error(inf)     = {a_inf - theta_young:+8.3f} deg")
        if abs(th.mean() - a_inf) > 1.0:
            print(f"    NOTE: the run stopped {abs(th.mean()-a_inf):.1f} deg short "
                  f"of its own extrapolated equilibrium; t_final is too small to\n"
                  f"          read the angle directly. Trust theta_inf, or extend to "
                  f"~{-tau*np.log(0.5/max(abs(th.mean()-a_inf),1e-9))/86400:.0f} days.")
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
        if ext is not None:
            a_inf, sig, tau, _, _ = ext
            tt = np.linspace(0.0, max(traj_t), 300)
            A = np.mean([th_ - a_inf for th_, t_ in zip(traj_th, traj_t)
                         if t_ > 0.4 * max(traj_t)]
                        ) / np.mean([np.exp(-t_ / tau) for t_ in traj_t
                                     if t_ > 0.4 * max(traj_t)])
            ax.plot(in_time_unit(tt, unit), a_inf + A * np.exp(-tt / tau),
                    color="0.45", ls=":", lw=1.4,
                    label=f"fit  $\\theta_\\infty$={a_inf:.2f}°, $\\tau$={tau/86400:.1f} d")
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
