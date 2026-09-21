#!/usr/bin/env python3
"""Gate the contact-angle MEASUREMENT against fields whose angle is known exactly.

The solver-side gates (verify_wall_bc.sh) prove the wall term is assembled
correctly. They say nothing about whether contact_angle.py reads the right
angle back out of a field -- and a measurement bug would look exactly like a
physics failure in the sweep. So, following the pattern of
enceladus_DSM/studies/molaro_2019/verification/verify_grain_shrinkage.py and
scripts/paraview_macros/verify_curvature.py: synthesise a phase field whose
answer is known in closed form, push it through the PRODUCTION measurement
function (imported, not reimplemented), and check what comes back.

Construction. A meniscus at equilibrium in a channel of height H is a circular
arc meeting both walls at theta. Put its centre on the midplane at (xc, H/2);
then cos(theta) = (H/2)/R fixes R = H/(2|cos theta|). The field is the model's
own equilibrium profile in the true distance to that circle,

    phi = 1/2 (1 + tanh(d / (2 eps))),    d = +-(|P - C| - R)

with the sign chosen so the ice lies to the LEFT of the arc. theta = 90 is the
degenerate R -> infinity case and is built as a flat interface.

Writes verify_contact_angle_measure.csv next to this script; exits non-zero on
failure.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "postprocess"))
from contact_angle import measure                       # noqa: E402

TOL_DEG = 1.0
H = 1.0e-4                 # channel height
EPS = 4.0e-7               # decay length; H/eps = 250, a well-resolved interface
NX, NY = 1400, 700


def arc(theta_deg, x_contact, side):
    """Signed-distance function of one meniscus arc, positive on the ice side.

    The arc meets BOTH walls at theta_deg and crosses y = 0 and y = H at
    x = x_contact. Its centre sits on the midplane, so cos(theta) = (H/2)/R
    fixes R = H/(2|cos theta|); whether the centre lies on the vapour side
    (wetting, concave) or the ice side (non-wetting, convex) is what carries
    the sign of cos(theta).

    `side` is +1 for a meniscus whose ice lies to the LEFT, -1 to the RIGHT.

    Returned distance is EXACT: the magnitude is the true distance to the
    circle, and the sign comes from which side of the arc the point is on,
    x vs x_s(y). Using (x_s - x) directly would overestimate the distance by
    1/cos(local slope) and bias the recovered angle.
    """
    c = np.cos(np.radians(theta_deg))

    if abs(c) < 1e-12:                       # theta = 90: R -> infinity, flat
        def f(X, Y):
            return side * (x_contact - X)
        return f

    R = H / (2.0 * abs(c))
    half = np.sqrt(max(R * R - 0.25 * H * H, 0.0))
    # c > 0 (wetting): centre on the vapour side, arc concave, s = -side.
    # c < 0: centre on the ice side, arc convex, s = +side.
    s = -side if c > 0.0 else side
    xc = x_contact - s * half

    def f(X, Y):
        x_s = xc + s * np.sqrt(np.maximum(R * R - (Y - 0.5 * H) ** 2, 0.0))
        return np.sign(side * (x_s - X)) * np.abs(np.hypot(X - xc, Y - 0.5 * H) - R)
    return f


def synth(theta_deg, slope=0.0):
    """(x1d, y1d, phi) for a BRIDGE spanning the channel, both menisci at theta.

    Two crossings per row is what pplib.contour_points needs, and it is also
    the geometry the real validation runs use -- a grain confined between two
    walls, relaxing into a bridge.
    """
    Lx = 3.0 * H
    x1d = np.linspace(0.0, Lx, NX)
    # Ruled between two wall curves, exactly as build_geometry_*.py does:
    # y_bot = -slope*x, y_top = H + slope*x, so slope=0 is the flat channel and
    # slope>0 opens like a wedge. The bridge is built in the SHEARED frame, so
    # its menisci meet each wall at theta_deg whatever the slope -- which is
    # what makes this a test of the sloped-wall path rather than a new geometry.
    X = np.repeat(x1d[:, None], NY, axis=1)
    yb = -slope * x1d
    yt = H + slope * x1d
    v = np.linspace(0.0, 1.0, NY)
    Y = yb[:, None] + (yt - yb)[:, None] * v[None, :]

    # Work in the UNSHEARED frame: map (x, y) back to a flat channel, build the
    # bridge there, and the shear carries it onto the sloped walls with the
    # angles preserved (the shear is y-only and the walls move with it).
    Yflat = (Y - yb[:, None]) / (yt - yb)[:, None] * H
    w = 0.7 * H
    d_right = arc(theta_deg, 0.5 * Lx + w, side=+1)(X, Yflat)
    d_left = arc(theta_deg, 0.5 * Lx - w, side=-1)(X, Yflat)
    d = np.minimum(d_right, d_left)          # ice = inside both menisci

    gap = d.max()
    assert gap > 20.0 * EPS, (
        f"theta={theta_deg}: the two menisci are only {gap/EPS:.1f} eps apart; "
        "they would overlap and the field would not be two clean interfaces")
    return X, Y, 0.5 * (1.0 + np.tanh(d / (2.0 * EPS)))


def build_walls(slope):
    wb = (0.0, -slope)
    wt = (H, slope)
    return [(wb, 0), (wt, 1)], (wb, wt)


def main():
    rows, fail = [], 0

    print("\n  CONTACT-ANGLE MEASUREMENT GATE")
    print("  synthetic arcs, H/eps = %.0f, grid %dx%d, tol %.1f deg" %
          (H / EPS, NX, NY, TOL_DEG))
    print("  slope 0 is a flat channel; slope 0.25 is the wedge taper.\n")
    print("   slope   theta_true   theta_meas   err     spread   n   result")

    # Flat walls only. A sloped-wall case was attempted here and WITHDRAWN: it
    # built the bridge in an unsheared frame and mapped it onto the tapered
    # walls, but that map scales y by an x-dependent factor and so does not
    # preserve angles -- the synthetic never had theta at the walls, and the
    # gate was measuring the distortion rather than the code. It reported
    # errors of up to 25 deg on a measurement that is in fact correct.
    #
    # Constructing a correct wedge synthetic is doable -- a circle centred on
    # the wedge axis at distance c from the apex meets both walls at
    # cos(theta) = c*sin(alpha)/rho, with c = 0 giving the apex-centred
    # 90-degree arc -- but it is a real piece of geometry and is not attempted
    # here.
    #
    # Until then the sloped-wall path is validated against the SOLVER instead,
    # which is the stronger test anyway: batch 2026-09-19 prescribes 60 and 120
    # degrees on a wedge mesh, and this measurement reads back 60.039 and
    # 119.953 from the sealed runs. Solver and measurement are independent, so
    # agreement to 0.05 deg is hard to get by accident.
    cases = [(0.0, th) for th in (15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165)]

    for slope, th in cases:
        X, Y, phi = synth(th, slope)
        walls, bounds = build_walls(slope)
        res = measure(X, Y, phi, EPS, walls, bounds, exclude_eps=5.0)
        if not res:
            print(f"   {slope:5.2f}   {th:8.1f}   {'--':>10}   no contour      FAIL")
            rows.append((slope, th, float("nan"), float("nan"), 0)); fail = 1
            continue
        meas = np.array([r["theta"] for r in res])
        err = meas.mean() - th
        ok = abs(err) < TOL_DEG
        fail |= (not ok)
        print("   %5.2f   %8.1f   %10.3f   %+6.3f  %6.3f  %2d   %s"
              % (slope, th, meas.mean(), err, meas.std(), meas.size,
                 "PASS" if ok else "FAIL"))
        rows.append((slope, th, meas.mean(), err, meas.size))

    csv = os.path.join(HERE, "verify_contact_angle_measure.csv")
    with open(csv, "w") as fh:
        fh.write("wall_slope,theta_true_deg,theta_measured_deg,error_deg,"
                 "n_estimates,tolerance_deg,result\n")
        for sl, th, m, e, n in rows:
            fh.write("%.2f,%.1f,%.4f,%.4f,%d,%.1f,%s\n"
                     % (sl, th, m, e, n, TOL_DEG,
                        "PASS" if abs(e) < TOL_DEG else "FAIL"))
    print(f"\n   -> {csv}")
    print("   " + ("ALL PASSED" if not fail else "FAILURES"))
    return int(fail)


if __name__ == "__main__":
    sys.exit(main())
