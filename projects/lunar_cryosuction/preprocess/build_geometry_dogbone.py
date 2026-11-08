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

THE ICE. A lens of length --lens-len that fills the pore (a plug), with flat
faces at t = 0 (90 degrees), built as in plan types 3 and 4: the
intersection of two far, crossed discs (-bridge_*). It relaxes to the run's
Young angle within the first days.

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
LX, W, NECK_LEN, SHOULDER = 1.2e-3, 1.5e-4, 3.0e-4, 1.0e-4
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


def preview(surf, lx, w, wc, x_n0, x_n1, ls, x_a, x_b, eps, fname, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from figcheck import check_figure

    x = np.linspace(0, lx, 3000)
    yb, yt = walls(x, w, wc, x_n0, x_n1, ls)
    s = 1e6
    fig, ax = plt.subplots(2, 1, figsize=(7.5, 4.6), constrained_layout=True,
                           gridspec_kw=dict(height_ratios=[1.0, 1.25]))
    a = ax[0]
    a.fill_between(x * s, 0, yb * s, color="0.6", lw=0)
    a.fill_between(x * s, yt * s, w * s, color="0.6", lw=0)
    a.plot(x * s, yb * s, "k-", lw=0.8)
    a.plot(x * s, yt * s, "k-", lw=0.8)
    xl = np.linspace(x_a, x_b, 50)
    a.fill_between(xl * s, 0, w * s, color="#66b3ff", lw=0)
    for xx in (x_n0, x_n1):
        a.axvline(xx * s, color="0.3", lw=0.6, ls=":")
    a.set_aspect("equal")
    a.set_xlim(0, lx * s)
    a.set_ylim(0, w * s)
    a.set_xlabel("x [µm]")
    a.set_ylabel("y [µm]")
    a.set_title(title, fontsize=8)

    # mesh close-up of one shoulder: every 4th element line
    b = ax[1]
    u = np.linspace(0, 1, 1500)
    xs0, xs1 = x_n1 - 0.2 * ls, x_n1 + 1.2 * ls
    for ub in surf.breaks(0)[::4]:
        if xs0 <= ub * lx <= xs1:
            pts = surf([ub], u)[0]
            b.plot(pts[:, 0] * s, pts[:, 1] * s, color="#1f6fd0", lw=0.3)
    for vb in surf.breaks(1)[::4]:
        pts = surf(u, [vb])[:, 0]
        m = (pts[:, 0] >= xs0) & (pts[:, 0] <= xs1)
        b.plot(pts[m, 0] * s, pts[m, 1] * s, color="#1f6fd0", lw=0.3)
    b.plot(x * s, yb * s, "k-", lw=0.8)
    b.plot(x * s, yt * s, "k-", lw=0.8)
    b.set_xlim(xs0 * s, xs1 * s)
    b.set_ylim(0, w * s)
    b.set_aspect("equal")
    b.set_xlabel("x [µm]")
    b.set_ylabel("y [µm]")
    b.set_title("mesh at the cold-side shoulder (every 4th element line)", fontsize=8)
    bad = check_figure(fig)
    fig.savefig(fname, dpi=200)
    plt.close(fig)
    print(f"wrote {fname}  figcheck: {bad or 'OK'}")
    return bad


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--Lx", type=float, default=LX, help="pore length [m]")
    ap.add_argument("--W", type=float, default=W, help="width of the wide sections [m]")
    ap.add_argument("--wc", type=float, required=True, help="neck width [m]; = W for no constriction")
    ap.add_argument("--neck-len", type=float, default=NECK_LEN, help="length of the straight neck [m]")
    ap.add_argument("--shoulder", type=float, default=SHOULDER, help="length of each shoulder [m]")
    ap.add_argument("--neck-x", type=float, default=None, help="neck center [m] (default Lx/2)")
    ap.add_argument("--lens-len", type=float, default=LENS_LEN, help="lens length [m]")
    ap.add_argument("--lens-x", type=float, default=None,
                    help="lens center [m] (default: middle of the flat wide section at x = Lx)")
    ap.add_argument("--eps", type=float, default=EPS, help="interface parameter [m]")
    ap.add_argument("--out", required=True, help="mesh .dat path")
    ap.add_argument("--opts", required=True, help="geometry .opts path")
    ap.add_argument("--png", default=None, help="preview figure (default: next to the .dat)")
    ap.add_argument("--header", default="", help="extra comment lines for the top of the .opts")
    a = ap.parse_args(argv)

    lx, w, wc, ls, eps = a.Lx, a.W, a.wc, a.shoulder, a.eps
    assert 0 < wc <= w, f"need 0 < wc <= W (wc {wc:.3e}, W {w:.3e})"
    xc_n = 0.5 * lx if a.neck_x is None else a.neck_x
    x_n0, x_n1 = xc_n - 0.5 * a.neck_len, xc_n + 0.5 * a.neck_len
    flat_lo, flat_hi = x_n0 - ls, x_n1 + ls          # wide sections: [0, flat_lo] and [flat_hi, Lx]
    assert flat_lo > 0 and flat_hi < lx, "neck and shoulders do not fit in the pore"
    x_l = 0.5 * (flat_hi + lx) if a.lens_x is None else a.lens_x
    x_a, x_b = x_l - 0.5 * a.lens_len, x_l + 0.5 * a.lens_len
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
    print(f"  lens x {x_a*1e6:.1f}-{x_b*1e6:.1f} um (area {a.lens_len*w:.4e} m^2); margins "
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
    png = Path(a.png).resolve() if a.png else dat.with_suffix(".png")
    bad = preview(surf, lx, w, wc, x_n0, x_n1, ls, x_a, x_b, eps, png,
                  f"dog-bone: W {w*1e6:.0f} µm, neck {wc*1e6:.0f} µm x {a.neck_len*1e6:.0f} µm, "
                  f"lens {x_a*1e6:.0f}-{x_b*1e6:.0f} µm (blue)")

    R = 50.0 * max(w, 1.25e-4)                  # far discs: face sag (W/2)^2/(2R)
    sag = (0.5 * w) ** 2 / (2 * R)
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
# Lens: x = {x_a*1e6:.1f}-{x_b*1e6:.1f} um, area {a.lens_len*w:.4e} m^2, flat faces (90 deg at
# t = 0) from two far, crossed discs (R = {R:.3e} m, face sag {sag*1e6:.3f} um).
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
-bridge_cxL {x_a + R:.6e}
-bridge_rL  {R:.6e}
-bridge_cxR {x_b - R:.6e}
-bridge_rR  {R:.6e}
-bridge_cy  {0.5*w:.6e}
-delt_t 1.0e-4
-eps {eps:.4e}
-eps_valid_temp {T0_C:g}
-periodic 0
-wall_faces y0,y1
""")
    print(f"wrote {opts}")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
