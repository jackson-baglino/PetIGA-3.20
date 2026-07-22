#!/usr/bin/env python3
"""Solve for the grain-centre separation that yields a target t=0 neck width.

WHY THIS EXISTS
---------------
The two-grain IC does not place a sharp lens. initial_conditions.c
(multi_grains branch) builds

    phi(x,y) = clamp( sum_k [ 0.5 - 0.5*tanh( 0.5*(r_k - R_k)/eps ) ], 0, 1 )

i.e. the grains' diffuse fields ADD. In the overlap region both terms
contribute, so the measured phi=0.5 waist is systematically WIDER than the
sharp-sphere chord 2*x0. For the -20 C geometry the documented offset is
x0 = 13.496 um (sharp) -> 16.405 um (measured radius) -> 32.81 um width:
nearly 3 um of radius, ~22%. Sizing a geometry from the sharp chord alone
therefore lands the initial neck well off target, and the offset depends on
eps, so it must be re-solved whenever eps changes.

This script replicates the C IC and postprocess/neck_width.py's measurement
convention exactly, then bisects on the centre separation d until the
measured width hits the target. No solver run required.

MEASUREMENT CONVENTION (mirrors neck_width.py)
  planar : chord = span between the outermost phi=0.5 crossings of a column.
  axisym : grains sit on the axis (cy = 0); a column's chord runs from the
           axis to the contour and IS the neck radius, so width = 2*chord.
  The neck is the interior minimum of w(x) strictly between the two grain
  centres.

VERIFY BEFORE TRUSTING: --check reproduces the committed -20 C geometry
(2D_molaro_axisym.opts) and compares against its documented 32.81 um.

Usage:
    python calibrate_neck_geometry.py --check
    python calibrate_neck_geometry.py --R0 76.5e-6 --R1 98.5e-6 \
        --neck 32.51e-6 --eps 3.9e-7 --axisym
"""

from __future__ import annotations

import argparse
import math

import numpy as np

PAD = 20.0e-6          # domain padding on every side [m] (matches the
                       # committed molaro geometries: Lx = d+R0+R1+2*PAD,
                       # Ly = R1+PAD)


def phi_field(x, y, cx, R, eps):
    """The C IC's ice field: additive tanh grains, clamped. cy = 0 (axisym)."""
    phi = np.zeros((len(y), len(x)))
    X, Y = np.meshgrid(x, y)
    for cxk, Rk in zip(cx, R):
        r = np.sqrt((X - cxk) ** 2 + Y ** 2)
        phi += 0.5 - 0.5 * np.tanh(0.5 * (r - Rk) / eps)
    return np.clip(phi, 0.0, 1.0)


def chord_width(col, y, level=0.5):
    """Span between outermost level-crossings — verbatim neck_width.py logic."""
    above = col >= level
    if not above.any():
        return 0.0
    idx = np.flatnonzero(above)
    lo_i, hi_i = idx[0], idx[-1]
    y_lo = y[lo_i]
    if lo_i > 0:
        f = (level - col[lo_i - 1]) / (col[lo_i] - col[lo_i - 1])
        y_lo = y[lo_i - 1] + f * (y[lo_i] - y[lo_i - 1])
    y_hi = y[hi_i]
    if hi_i < len(y) - 1:
        f = (col[hi_i] - level) / (col[hi_i] - col[hi_i + 1])
        y_hi = y[hi_i] + f * (y[hi_i + 1] - y[hi_i])
    return y_hi - y_lo


def measure_neck(d, R0, R1, eps, axisym=True, nx=3000, ny=1200):
    """Measured t=0 neck width for centre separation d."""
    cx = [PAD + R0, PAD + R0 + d]
    Lx = d + R0 + R1 + 2 * PAD
    Ly = R1 + PAD
    x = np.linspace(0.0, Lx, nx)
    y = np.linspace(0.0, Ly, ny)
    phi = phi_field(x, y, cx, [R0, R1], eps)

    w = np.array([chord_width(phi[:, j], y) for j in range(len(x))])
    if axisym:
        w = 2.0 * w

    # interior minimum strictly between the grain centres
    j0 = int(np.searchsorted(x, cx[0]))
    j1 = int(np.searchsorted(x, cx[1]))
    seg = w[j0 + 1:j1]
    jm = int(np.argmin(seg)) + j0 + 1
    return w[jm], x[jm], cx, Lx, Ly


def sharp_chord(d, R0, R1):
    """Sharp two-sphere lens chord half-width.

    With the ADDITIVE IC this is NOT the measured neck (the diffuse skin adds
    to it, eps-dependently). With -ic_grain_union it IS the measured neck, for
    every eps -- which is the whole point of that flag.
    """
    a = (d * d + R0 * R0 - R1 * R1) / (2 * d)
    v = R0 * R0 - a * a
    return math.sqrt(v) if v > 0 else 0.0


def solve_d_union(R0, R1, target_width):
    """Separation d whose SHARP chord equals target_width.

    This is the -ic_grain_union sizing: phi=0.5 sits on the sharp union surface
    at any eps, so the sharp chord IS the initial neck and no eps-dependent
    calibration is needed. One d serves an entire eps series.
    """
    lo, hi = 0.5 * (R0 + R1), R0 + R1
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if 2.0 * sharp_chord(mid, R0, R1) > target_width:
            lo = mid          # chord too wide -> less overlap -> larger d
        else:
            hi = mid
    return 0.5 * (lo + hi)


def solve_d(R0, R1, eps, target, axisym=True, **kw):
    """Bisect on d for measured neck width == target."""
    lo, hi = 0.80 * (R0 + R1), R0 + R1        # more overlap -> wider neck
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        w, *_ = measure_neck(mid, R0, R1, eps, axisym, **kw)
        if w > target:
            lo = mid          # too wide -> less overlap -> larger d
        else:
            hi = mid
        if hi - lo < 1e-12:
            break
    return 0.5 * (lo + hi)


def emit(name, R0, R1, eps, target, axisym):
    d = solve_d(R0, R1, eps, target, axisym)
    w, xn, cx, Lx, Ly = measure_neck(d, R0, R1, eps, axisym)
    x0 = sharp_chord(d, R0, R1)
    Nx = math.ceil(Lx * math.sqrt(2) / eps)
    Ny = math.ceil(Ly * math.sqrt(2) / eps)
    print(f"\n=== {name}")
    print(f"  R0/R1            = {R0*1e6:.2f} / {R1*1e6:.2f} um")
    print(f"  eps              = {eps:.4e} m")
    print(f"  target neck      = {target*1e6:.2f} um")
    print(f"  -> separation d  = {d:.6e} m")
    print(f"  measured neck    = {w*1e6:.3f} um   (err {(w-target)*1e9:+.1f} nm)")
    print(f"  sharp chord 2*x0 = {2*x0*1e6:.3f} um  "
          f"(x0 = {x0*1e6:.3f} um; diffuse adds {(w-2*x0)*1e6:+.3f} um)")
    print(f"  neck plane x     = {xn*1e6:.3f} um")
    print(f"  -ice_grain_cx {cx[0]:.5e},{cx[1]:.5e}")
    print(f"  -Lx {Lx:.5e}   -Ly {Ly:.5e}")
    print(f"  -Nx {Nx}   -Ny {Ny}")
    return dict(d=d, w=w, cx=cx, Lx=Lx, Ly=Ly, Nx=Nx, Ny=Ny)


def check():
    """Reproduce the committed -20 C axisym geometry as a self-test."""
    R0, R1, eps = 7.25e-5, 1.01e-4, 4.64e-7
    cx = [9.25e-5, 2.6383e-4]
    d = cx[1] - cx[0]
    w, xn, _, Lx, Ly = measure_neck(d, R0, R1, eps, axisym=True)
    x0 = sharp_chord(d, R0, R1)
    print("=== SELF-CHECK vs committed 2D_molaro_axisym.opts")
    print(f"  d = {d:.5e} m (from its -ice_grain_cx)")
    print(f"  sharp x0        = {x0*1e6:.3f} um   [opts documents 13.496]")
    print(f"  measured radius = {w/2*1e6:.3f} um  [opts documents 16.405]")
    print(f"  measured width  = {w*1e6:.3f} um    [opts documents 32.81]")
    print(f"  Lx = {Lx:.5e} [opts 3.848e-4]   Ly = {Ly:.5e} [opts 1.21e-4]")
    ok = abs(w * 1e6 - 32.81) < 0.15 and abs(x0 * 1e6 - 13.496) < 0.01
    print(f"  => replication {'OK' if ok else 'MISMATCH — do not trust output'}")
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true",
                    help="self-test against the committed -20 C geometry")
    ap.add_argument("--R0", type=float, help="smaller grain radius [m]")
    ap.add_argument("--R1", type=float, help="larger grain radius [m]")
    ap.add_argument("--neck", type=float, help="target t=0 neck WIDTH [m]")
    ap.add_argument("--eps", type=float, help="interface parameter [m]")
    ap.add_argument("--axisym", action="store_true")
    ap.add_argument("--name", default="geometry")
    args = ap.parse_args()

    if args.check:
        check()
        return
    for req in ("R0", "R1", "neck", "eps"):
        if getattr(args, req) is None:
            ap.error(f"--{req} is required (or use --check)")
    emit(args.name, args.R0, args.R1, args.eps, args.neck, args.axisym)


if __name__ == "__main__":
    main()
