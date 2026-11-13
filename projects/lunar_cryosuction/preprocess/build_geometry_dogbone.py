#!/usr/bin/env python3
"""build_geometry_dogbone.py: a dog-bone pore (plan type 6) with a flat-faced
ice lens in one wide section.

THE PORE. Wide sections of width W at both ends, a straight neck of width wc
in the middle, symmetric about the pore axis y = W/2:

    x:  0 ......... wide ......... | shoulder | neck (straight) | shoulder | ......... wide ......... Lx

The half-width is

    h(x) = W/2 - (W - wc)/2 * b(x),   b(x) = psi((x - x_n0 + Ls)/Ls) * psi((x_n1 + Ls - x)/Ls)

with the neck on [x_n0, x_n1] and shoulders of length Ls on either side.
psi is the C-infinity step

    psi(t) = f(t) / (f(t) + f(1 - t)),   f(t) = exp(-1/t) for t > 0, else 0,

which is exactly 0 for t <= 0 and exactly 1 for t >= 1. So the neck is
exactly straight on [x_n0, x_n1] and each wide section is exactly flat
outside the shoulders. The steepest wall slope is 2 (W - wc)/(2 Ls), at the
middle of each shoulder.

THE ICE. A lens of length --lens-len (measured on the pore axis, apex to
apex) that fills the pore (a plug). Its faces start at --lens-theta, the
contact angle measured through the ice, which should be the run's wall-BC
angle (user, 2026-11-09: a lens first relaxes to the BC angle, so start it
there):
  theta = 90   exactly flat faces: -ice_flat_x/-ice_flat_R/-ice_flat_height
               (IceFlatField, src/initial_conditions.c), ice on |x - x_l| < L/2
               below a height well above the top wall. Until 2026-11-12 this
               case used two far, crossed discs (R = 7.5 mm), whose faces sat
               (W/2)^2/(2R) = 0.375 um out of flat with the axis ahead; the
               user ruled that seed unacceptable (6B-1 is void).
  theta < 90   concave faces, arcs of radius R = (W/2)/cos(theta) centered on
               the axis outside the ice (-bridge_concaveL/R 1, FREEZE ledger
               18270f6f). The ice reaches sag = R - sqrt(R^2 - (W/2)^2) past
               each apex along the walls.
  theta > 90   convex faces, R = (W/2)/|cos(theta)|, centered inside the ice;
               the wall contacts sit sag inside the apexes.

WHY THE OPTS FILE GIVES FLAT WALLS. The solver lays the initial condition out
by mapping each node through WallBottom(x)/WallTop(x) (src/initial_conditions.c),
and those know only an affine wall plus C-infinity bumps. They cannot express
this profile, and adding a wall option to the solver is a change to the
frozen source. It is not needed: the initial ice field is 0 everywhere except
in the lens, and the lens is placed where the walls are exactly flat
(b(x) = 0 exactly, asserted below with a margin). The temperature IC depends
only on x. So the initial condition is exact with the flat-wall default
(walls y = 0 and y = W). Everything after t = 0 (wall boundary condition,
contact angle, transport) uses the mesh in the .dat file, which carries the
true dog-bone walls. Postprocessing that reads wall positions from the opts
file must not be used inside the shoulders or neck; read them from the mesh
or from the "# DOGBONE:" line this script writes.

Usage (from the project root):
    venv_cryosuction/bin/python preprocess/build_geometry_dogbone.py --wc 4e-5
    venv_cryosuction/bin/python preprocess/build_geometry_dogbone.py --wc 8e-5 --tag wc80 \\
        --out inputs/geometry/meshes/dogbone_wc80.dat --opts inputs/geometry/cryosuction/x.opts

Every number used is printed.
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "preprocess"))
sys.path.insert(0, str(ROOT / "postprocess"))

T0_C = -20.0
EPS = 8.5840e-7          # production eps (issue 1 closed 2026-11-07)
P, C = 2, 1
LX, W, NECK_LEN, SHOULDER = 1.1e-3, 1.5e-4, 3.0e-4, 1.0e-4
LENS_LEN = 1.5e-4
MIN_MARGIN_EPS = 20.0    # lens to shoulder and lens to pore end, in eps


def psi(t):
    """C-infinity step: 0 for t <= 0, 1 for t >= 1."""
    t = np.asarray(t, float)

    def f(s):
        out = np.zeros_like(s)
        m = s > 0
        out[m] = np.exp(-1.0 / s[m])
        return out
    a, b = f(t), f(1.0 - t)
    return a / (a + b)


def neck_weight(x, x_n0, x_n1, ls):
    """b(x): 1 on the neck, 0 on the wide sections, psi on the shoulders."""
    x = np.asarray(x, float)
    return psi((x - x_n0 + ls) / ls) * psi((x_n1 + ls - x) / ls)


def walls(x, w, wc, x_n0, x_n1, ls):
    h = 0.5 * w - 0.5 * (w - wc) * neck_weight(x, x_n0, x_n1, ls)
    return 0.5 * w - h, 0.5 * w + h


def build_surface(lx, w, wc, x_n0, x_n1, ls, nx, ny):
    from build_geometry_multi_grain import greville_abscissae, open_uniform_knots
    from igakit.nurbs import NURBS
    ux, uy = open_uniform_knots(nx, P, C), open_uniform_knots(ny, P, C)
    gx, gy = greville_abscissae(ux, P), greville_abscissae(uy, P)
    x = lx * gx
    yb, yt = walls(x, w, wc, x_n0, x_n1, ls)
    ctrl = np.zeros((len(gx), len(gy), 3))
    for i in range(len(gx)):
        ctrl[i, :, 0] = x[i]
        ctrl[i, :, 1] = yb[i] + gy * (yt[i] - yb[i])
    return NURBS([ux, uy], control=ctrl)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--Lx", type=float, default=LX, help="pore length [m]")
    ap.add_argument("--W", type=float, default=W, help="width of the wide sections [m]")
    ap.add_argument("--wc", type=float, required=True, help="neck width [m]; = W for no constriction")
    ap.add_argument("--neck-len", type=float, default=NECK_LEN, help="length of the straight neck [m]")
    ap.add_argument("--shoulder", type=float, default=SHOULDER, help="length of each shoulder [m]")
    ap.add_argument("--neck-x", type=float, default=None, help="neck center [m] (default Lx/2)")
    ap.add_argument("--lens-len", type=float, default=LENS_LEN, help="lens length [m]")
    ap.add_argument("--lens-theta", type=float, default=90.0,
                    help="initial face angle through the ice [deg] (use the wall-BC angle)")
    ap.add_argument("--lens-x", type=float, default=None,
                    help="lens center [m] (default: middle of the flat wide section at x = Lx)")
    ap.add_argument("--eps", type=float, default=EPS, help="interface parameter [m]")
    ap.add_argument("--out", required=True, help="mesh .dat path")
    ap.add_argument("--opts", required=True, help="geometry .opts path")
    ap.add_argument("--png", default=None, help="preview PDF (default preprocess/<opts name>.pdf)")
    ap.add_argument("--header", default="", help="extra comment lines for the top of the .opts")
    a = ap.parse_args(argv)

    lx, w, wc, ls, eps = a.Lx, a.W, a.wc, a.shoulder, a.eps
    assert 0 < wc <= w, f"need 0 < wc <= W (wc {wc:.3e}, W {w:.3e})"
    xc_n = 0.5 * lx if a.neck_x is None else a.neck_x
    x_n0, x_n1 = xc_n - 0.5 * a.neck_len, xc_n + 0.5 * a.neck_len
    flat_lo, flat_hi = x_n0 - ls, x_n1 + ls          # wide sections: [0, flat_lo] and [flat_hi, Lx]
    assert flat_lo > 0 and flat_hi < lx, "neck and shoulders do not fit in the pore"
    x_l = 0.5 * (flat_hi + lx) if a.lens_x is None else a.lens_x
    x_a, x_b = x_l - 0.5 * a.lens_len, x_l + 0.5 * a.lens_len   # apexes on the axis
    th = a.lens_theta
    assert 0.0 < th < 180.0, "--lens-theta must be in (0, 180) deg"
    ct = math.cos(math.radians(th))
    flat = abs(ct) < 1e-12
    if flat:                                    # exactly flat faces (IceFlatField); no arc
        R, concave = math.inf, False
    else:
        R, concave = 0.5 * w / abs(ct), ct > 0
    sag = 0.0 if flat else R - math.sqrt(R * R - (0.5 * w) ** 2)
    # extent of the ice along x, at the walls for concave faces, on the axis otherwise
    x_a, x_b, ax_a, ax_b = (x_a - sag, x_b + sag, x_a, x_b) if concave else (x_a, x_b, x_a, x_b)
    if x_a > flat_hi:
        m_shoulder, m_end = x_a - flat_hi, lx - x_b
    elif x_b < flat_lo:
        m_shoulder, m_end = flat_lo - x_b, x_a
    else:
        raise SystemExit(f"lens {x_a:.3e}-{x_b:.3e} m is not inside a flat wide section")
    lim = MIN_MARGIN_EPS * eps
    assert m_shoulder >= lim and m_end >= lim, \
        f"lens margins {m_shoulder:.3e} (shoulder), {m_end:.3e} (end) m; need >= {lim:.3e}"
    xs = np.linspace(x_a, x_b, 2001)
    assert np.all(neck_weight(xs, x_n0, x_n1, ls) == 0.0), "walls are not exactly flat at the lens"

    nx = math.ceil(lx * math.sqrt(2) / eps)
    ny = math.ceil(w * math.sqrt(2) / eps)
    slope = 2.0 * 0.5 * (w - wc) / ls
    print(f"dog-bone pore  Lx {lx*1e6:.1f} um, W {w*1e6:.1f} um, neck wc {wc*1e6:.1f} um "
          f"on x {x_n0*1e6:.1f}-{x_n1*1e6:.1f} um, shoulders {ls*1e6:.1f} um")
    print(f"  wide sections flat on x 0-{flat_lo*1e6:.1f} and {flat_hi*1e6:.1f}-{lx*1e6:.1f} um")
    print(f"  steepest wall slope {slope:.3f} ({math.degrees(math.atan(slope)):.1f} deg), mid-shoulder")
    hw = 0.5 * w        # sharp lens area: rectangle -/+ 2 circular segments beyond the axis chord
    seg = 0.0 if flat else 2 * hw * R - hw * math.sqrt(R * R - hw * hw) - R * R * math.asin(hw / R)
    area = a.lens_len * w + (2 * seg if concave else -2 * seg)
    print(f"  lens apexes {ax_a*1e6:.1f}-{ax_b*1e6:.1f} um, ice x {x_a*1e6:.1f}-{x_b*1e6:.1f} um, "
          f"faces at {th:g} deg (R {R:.3e} m, {'concave' if concave else 'convex/flat'})")
    print(f"  lens area {area:.4e} m^2; margins "
          f"{m_shoulder*1e6:.1f} um to the shoulder, {m_end*1e6:.1f} um to the end "
          f"({m_shoulder/eps:.0f} and {m_end/eps:.0f} eps)")
    print(f"  eps {eps:.4e} m, Nx {nx}, Ny {ny}, DOF {(nx+P)*(ny+P)*3/1e6:.2f} M; element "
          f"{lx/nx*1e6:.3f} x {w/ny*1e6:.3f} um wide, {lx/nx*1e6:.3f} x {wc/ny*1e6:.3f} um in the neck "
          f"(aspect {(lx/nx)/(wc/ny):.2f})")
    if wc / w > 0.999:
        print("  wc = W: no constriction (a flat pore), kept on the same mesh path for uniformity")

    surf = build_surface(lx, w, wc, x_n0, x_n1, ls, nx, ny)
    # The mesh wall is the B-spline through control points on the profile; how far is it from the profile?
    u = np.linspace(0, 1, 20001)
    pb = surf(u, [0.0])[:, 0]
    pt = surf(u, [1.0])[:, 0]
    eb = np.abs(pb[:, 1] - walls(pb[:, 0], w, wc, x_n0, x_n1, ls)[0]).max()
    et = np.abs(pt[:, 1] - walls(pt[:, 0], w, wc, x_n0, x_n1, ls)[1]).max()
    print(f"  mesh wall vs profile: max deviation {max(eb, et):.3e} m ({max(eb, et)/eps:.2e} eps)")
    assert max(eb, et) < 0.01 * eps, "mesh wall departs from the profile by more than 0.01 eps"

    from igakit.io import PetIGA
    dat, opts = Path(a.out).resolve(), Path(a.opts).resolve()
    dat.parent.mkdir(parents=True, exist_ok=True)
    PetIGA().write(str(dat), surf)
    print(f"wrote {dat}")
    if concave:      # disc centers outside the ice; apex = cxL + R = ax_a, cxR - R = ax_b
        bridge = (f"-bridge_cxL {ax_a - R:.6e}\n-bridge_rL  {R:.6e}\n-bridge_cxR {ax_b + R:.6e}\n"
                  f"-bridge_rR  {R:.6e}\n-bridge_cy  {0.5*w:.6e}\n-bridge_concaveL 1\n-bridge_concaveR 1\n")
        face = f"concave faces at {th:g} deg (R = {R:.3e} m, the ice reaches {sag*1e6:.3f} um past each apex at the walls)"
    elif flat:       # ice on [ax_a, ax_b], full height (the top wall is at y = w)
        bridge = (f"-ice_flat_x {x_l:.6e}\n-ice_flat_R {0.5*a.lens_len:.6e}\n-ice_flat_height {2*w:.6e}\n")
        face = f"exactly flat faces (90 deg; -ice_flat_*, ice on x = {ax_a*1e6:.1f}-{ax_b*1e6:.1f} um at every y)"
    else:            # disc centers inside the ice; apex = cxL - R = ax_a, cxR + R = ax_b
        bridge = (f"-bridge_cxL {ax_a + R:.6e}\n-bridge_rL  {R:.6e}\n-bridge_cxR {ax_b - R:.6e}\n"
                  f"-bridge_rR  {R:.6e}\n-bridge_cy  {0.5*w:.6e}\n")
        face = (f"convex faces at {th:g} deg (R = {R:.3e} m, wall contacts {sag*1e6:.3f} um inside each apex)")
    try:
        rel = dat.relative_to(ROOT)
    except ValueError:
        rel = dat
    hdr = "".join(f"# {ln}\n" for ln in a.header.splitlines())
    opts.write_text(f"""# =============================================================================
# geometry {opts.name}   (preprocess/build_geometry_dogbone.py)
{hdr}#
# Dog-bone pore: wide sections W = {w*1e6:.1f} um, straight neck wc = {wc*1e6:.1f} um on
# x = {x_n0*1e6:.1f}-{x_n1*1e6:.1f} um, C-infinity shoulders {ls*1e6:.1f} um long (steepest wall
# slope {slope:.3f}); walls flat on x = 0-{flat_lo*1e6:.1f} and {flat_hi*1e6:.1f}-{lx*1e6:.1f} um.
# The walls are in the mesh ({rel.name}); y0 and y1 are wall faces, x0 and x1 are
# the pore ends.
#
# The -Ly / default flat walls below are ONLY for the initial-condition layout.
# The lens sits where the walls are exactly flat (y = 0 and y = W), so the IC
# is exact there; elsewhere the IC is ice-free. See the builder's docstring.
# Postprocessing must take wall positions from the mesh or from this line:
# DOGBONE: W {w:.6e} wc {wc:.6e} x_n0 {x_n0:.6e} x_n1 {x_n1:.6e} shoulder {ls:.6e}
#
# Lens: apexes on the axis x = {ax_a*1e6:.1f}-{ax_b*1e6:.1f} um, ice x = {x_a*1e6:.1f}-{x_b*1e6:.1f} um,
# {face}.
# eps {eps:.4e} m (production value, issue 1), h = eps/sqrt(2): Nx {nx}, Ny {ny}.
# =============================================================================
# DOF_GRID: {nx+P} {ny+P}
-geom_file {rel}
-p {P}
-C {C}
-ic_type multi_grains
-dim 2
-Lx {lx:.6e}
-Ly {w:.6e}
-Lz 0
{bridge}-delt_t 1.0e-4
-eps {eps:.4e}
-eps_valid_temp {T0_C:g}
-periodic 0
-wall_faces y0,y1
""")
    print(f"wrote {opts}")
    from plot_geometry import render
    png = Path(a.png).resolve() if a.png else ROOT / "preprocess" / (opts.stem + ".pdf")
    bad = render(opts, png, f"dog-bone: W {w*1e6:.0f} µm, neck {wc*1e6:.0f} µm x {a.neck_len*1e6:.0f} µm, "
                            f"wide sections {flat_lo*1e6:.0f} µm, lens {ax_a*1e6:.0f}-{ax_b*1e6:.0f} µm at {th:g}° (blue)")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
