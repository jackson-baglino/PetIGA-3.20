#!/usr/bin/env python3
"""build_geometry_two_throat.py — a pore channel with TWO pore throats of
different aperture, each holding a frozen capillary bridge, for the icy-lunar
regolith study (studies/icy_regolith/implicit_pore_domain/, Effort 1).

THE QUESTION. The wedge run showed that a single ice deposit migrates along a
taper (-18.75 um/yr toward the narrow end). This asks the network version:
given two ice bridges plugging two throats of DIFFERENT aperture, does ice
redistribute from one to the other, and which one wins?

PREDICTION: ice moves from the WIDE throat to the NARROW one.

Each bridge is bounded by two menisci that meet both walls at 90 degrees, the
model's natural (Neumann) contact angle. Because the wall bump is convex INTO
the channel, those menisci come out CONCAVE -- the waisted hourglass shape of a
real capillary bridge -- so kappa < 0. The narrower throat holds its ice at
MORE negative kappa, hence at LOWER equilibrium vapour density by Gibbs-Thomson,
which makes it the vapour SINK. Ice therefore collects in the tighter throat,
the same sign as capillary condensation filling small pores first.

NOTE THE CONTRAST WITH THE WEDGE, which is easy to get backwards. In the wedge
the walls are STRAIGHT, the 90-degree meniscus is an apex-centred arc, and a
plug's outer meniscus is CONVEX. Here the walls are curved and their curvature
dominates, flipping the meniscus sign. A throat is still a stable POSITION for
a plug (that is the wedge result applied locally), but between two throats it
is the aperture that decides, and the tighter one wins. Do not carry the
straight-wall intuition over; this script prints the measured kappa of every
meniscus so the direction is never assumed.

GEOMETRY. Flat walls at y=0 and y=Ly with one C-infinity bump per wall per
throat, mirrored top/bottom so the channel is symmetric about y=Ly/2. Both
throats use the SAME bump half-width and differ only in height, so aperture is
the only variable. Gap at a throat = Ly - 2*bump_height.

THE 90-DEGREE MENISCUS. For a channel symmetric about its axis, an arc centred
ON the axis at (xc, Ly/2) meets a wall perpendicularly exactly when the centre
lies on that wall's normal at the contact point:

    xc = xp + y_b'(xp) * (y_b(xp) - Ly/2)

so choosing the wall contact xp fixes both xc and rho = |P - C|. By symmetry the
same arc meets the opposite wall at 90 degrees too. Unlike the wedge this is NOT
an apex-centred circle -- the radius has to be solved against the actual wall
slope, which is why these numbers must come from this script and never be
hand-written into an .opts.

Usage (from the project root):
    python3 preprocess/build_geometry_two_throat.py --orientation narrow-left
    python3 preprocess/build_geometry_two_throat.py --orientation narrow-right
"""

import argparse
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "preprocess"))
from build_geometry_multi_grain import generate_random_bumps  # noqa: E402

# --------------------------------------------------------------------------- #
# Domain and throat parameters
# --------------------------------------------------------------------------- #
Ly = 1.6e-4              # bounding-box height [m]
Lx = 8.5e-4              # domain length [m]
X_THROAT = (2.5e-4, 6.0e-4)   # throat centres; separation 3.5e-4 m

# Wall roughness: the rest of the channel is real pore space, not flat. Bumps are
# drawn from a lunar-regolith GSD (docs/material_parameters.md 2.3) exactly as in
# build_geometry_regolith_pore.py, then the two throats are imposed on top.
#
# The two walls are MIRRORED about y = Ly/2. That is a deliberate simplification,
# not an oversight: it puts the channel axis on a straight line, which is what
# makes the 90-degree meniscus construction exact (the arc centre can be placed
# on the axis and the two wall contacts are then symmetric). With independently
# rough walls the axis wanders, the arc centre is the intersection of the two
# wall TANGENTS, and the contacts must be co-solved -- doable, but it trades an
# exactly-90-degree initial condition for an approximate one, and every wedge run
# showed that the contact angle is exactly what must not be approximated.
N_ROUGH = 20             # roughness bumps per wall
ROUGH_R = (2.2e-5, 4.0e-5)
ROUGH_HR = (0.14, 0.24)  # height/half-width: gentle, to leave slope budget for the throats
ROUGH_OV = (0.20, 0.32)
SEED_ROUGH = 5

# Each throat is one extra wide bump whose height is SOLVED so the gap hits the
# target exactly, on top of whatever roughness is already there.
THROAT_R = 1.6e-4

# Same half-width for both throats so APERTURE is the only difference. R is set
# by the mesh-quality budget: a C-infinity bump's max slope is 2.17*h/R, and the
# taller (narrow-throat) bump must stay under ~1.05 (46 deg), the slope the
# validated regolith geometries run at.
GAP_NARROW = 3.0e-5      # 3x contrast
GAP_WIDE = 9.0e-5
SLOPE_BUDGET = 1.10

# Bridge sizing. Both bridges are given the SAME axial span, i.e. the same ice
# extent measured along the channel axis, and the wall contact point needed to
# achieve it is solved per throat. This is what makes aperture the single
# variable: with the span matched, the narrow throat necessarily holds LESS ice
# (its gap is smaller), so it is the smaller bridge by mass without any second
# knob being turned.
#
# The span must also comfortably exceed the diffuse interface, or the bridge has
# no flat core and simply dissolves instead of exchanging: the 1%-99% band is
# 9.2*eps = 7.9e-6 m, and a fixed contact fraction of 0.5 gave the WIDE throat
# only 1.9e-5 (2.4 bands) while handing the narrow throat 6.8e-5. Matching the
# span fixes both at ~8.6 bands.
BRIDGE_SPAN = 6.8e-5
MIN_SPAN_BANDS = 4.0        # assert span >= this many 9.2*eps interface widths

# eps: the comp_eps.py (Kaempfer & Plapp) value for T=-20 C, alpha_c=1.341e-2 --
# the same kinetic-bound value every other -20 C geometry in this study uses.
T0_C, ALPHA_C, EPS = -20.0, 1.341e-2, 8.5840e-7
P, C = 2, 1

D0 = 1.0166e-9           # capillary length at -20 C (matches -d0_sub0)
LAT, RV, TK = 2.83e6, 461.5, 253.15   # for the thermal-drive comparison


def bump(x, c, R, h):
    """C-infinity bump -- must match SedimentBump() in src/initial_conditions.c."""
    x = np.atleast_1d(np.asarray(x, float))
    t = (x - c) / R
    out = np.zeros_like(x)
    m = np.abs(t) < 1.0
    out[m] = h * np.exp(1.0 - 1.0 / (1.0 - t[m] ** 2))
    return out


def wall_bot(x, bumps):
    """Bottom wall y, i.e. the summed bump field rising from y=0."""
    x = np.atleast_1d(np.asarray(x, float))
    y = np.zeros_like(x)
    for c, R, h in bumps:
        y = y + bump(x, c, R, h)
    return y


def dwall_bot(x, bumps, e=1e-10):
    return (wall_bot(x + e, bumps) - wall_bot(x - e, bumps)) / (2.0 * e)


def meniscus(xp, bumps, cy):
    """The 90-degree meniscus arc whose bottom-wall contact is at xp.

    Ninety degrees means the arc's TANGENT is perpendicular to the wall's, i.e.
    the arc's RADIUS is parallel to the wall tangent, i.e. the centre lies along
    the wall's TANGENT LINE at the contact:

        xc = xp + (cy - y_wall(xp)) / y_wall'(xp)

    Putting the centre on the wall NORMAL instead is the classic sign error: it
    makes the arc TANGENT to the wall, a 0-degree contact. For a straight wall
    this returns the wedge apex, which is the cross-check that the sign is right.

    Returns (xc, rho). The centre sits on the axis, so by mirror symmetry the
    same arc meets the opposite wall at 90 degrees too."""
    yb = float(wall_bot(xp, bumps)[0])
    dp = float(dwall_bot(xp, bumps)[0])
    if abs(dp) < 1e-9:
        raise AssertionError(
            f"wall is flat at x={xp:.3e} (slope {dp:.2e}); the 90-degree meniscus "
            "there is a straight line of infinite radius -- move the contact off "
            "the throat crest")
    xc = xp + (cy - yb) / dp
    rho = math.hypot(xp - xc, yb - cy)
    return xc, rho


def contact_angle_deg(xc, rho, xp, bumps, cy, top=False):
    """Measured angle between the arc and the wall at the contact -- the gate."""
    yb = float(wall_bot(xp, bumps)[0])
    y = (Ly - yb) if top else yb
    d = float(dwall_bot(xp, bumps)[0]) * (-1.0 if top else 1.0)
    rx, ry = xp - xc, y - cy
    L = math.hypot(rx, ry)
    tx, ty = -ry / L, rx / L                      # arc tangent
    wx, wy = 1.0 / math.hypot(1.0, d), d / math.hypot(1.0, d)
    return math.degrees(math.acos(min(1.0, abs(tx * wx + ty * wy))))


def build_walls(gap_narrow, gap_wide, orientation):
    """Rough regolith wall + two throat bumps whose heights are SOLVED so each
    throat's gap hits its target exactly on top of the roughness."""
    rough = generate_random_bumps(np.random.default_rng(SEED_ROUGH), Lx, N_ROUGH,
                                  ROUGH_R, ROUGH_HR, ROUGH_OV)
    gaps = ((gap_narrow, gap_wide) if orientation == "narrow-left"
            else (gap_wide, gap_narrow))
    bumps = list(rough)
    for xt, g in zip(X_THROAT, gaps):
        want = 0.5 * (Ly - g)                       # wall height needed at xt
        have = float(wall_bot(xt, bumps)[0])        # what the roughness already gives
        assert want > have, \
            (f"throat at x={xt:.3e} wants wall height {want:.3e} m but the "
             f"roughness is already {have:.3e} m -- lower ROUGH_HR or widen the gap")
        bumps.append((xt, THROAT_R, want - have))
    return bumps, gaps


def preview(bumps, bridges, fname, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x = np.linspace(0, Lx, 4000)
    yb, yt = wall_bot(x, bumps), Ly - wall_bot(x, bumps)
    fig, ax = plt.subplots(figsize=(13, 13 * Ly / Lx + 1.6))
    ax.fill_between(x, 0, yb, color="0.55", zorder=1)
    ax.fill_between(x, yt, Ly, color="0.55", zorder=1)
    ax.plot(x, yb, "k-", lw=1.2, zorder=3)
    ax.plot(x, yt, "k-", lw=1.2, zorder=3)

    gx, gy = np.meshgrid(np.linspace(0, Lx, 2200), np.linspace(0, Ly, 420))
    wb = np.zeros_like(gx)
    for c, R, h in bumps:
        wb += bump(gx.ravel(), c, R, h).reshape(gx.shape)
    inside = (gy >= wb) & (gy <= Ly - wb)
    ice = np.zeros_like(gx, dtype=bool)
    for b in bridges:
        dL = np.hypot(gx - b["cxL"], gy - b["cy"])
        dR = np.hypot(gx - b["cxR"], gy - b["cy"])
        ice |= (dL <= b["rL"]) & (dR <= b["rR"])
    ice &= inside
    ax.contourf(gx, gy, ice.astype(float), levels=[0.5, 1.5],
                colors=["#66b3ff"], zorder=2)
    ax.contour(gx, gy, ice.astype(float), levels=[0.5],
               colors=["#1f6fd0"], linewidths=1.0, zorder=3)
    ax.set_xlim(0, Lx); ax.set_ylim(0, Ly); ax.set_aspect("equal")
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")
    fig.tight_layout(); fig.savefig(fname, dpi=140); plt.close(fig)
    print(f"wrote {fname}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--orientation", choices=("narrow-left", "narrow-right"),
                    default="narrow-left",
                    help="which throat is the narrow one (default: narrow-left)")
    ap.add_argument("--gap-narrow", type=float, default=GAP_NARROW)
    ap.add_argument("--gap-wide", type=float, default=GAP_WIDE)
    ap.add_argument("--contact-frac", type=float, default=0.55,
                    help="meniscus wall contact, as a fraction of THROAT_R from the "
                         "throat centre; larger = more ice and gentler menisci")
    ap.add_argument("--tag", default=None, help="override the output name suffix")
    args = ap.parse_args()

    tag = args.tag if args.tag else args.orientation.replace("-", "_")
    dat = ROOT / f"inputs/geometry/two_throat_{tag}.dat"
    opts = ROOT / f"inputs/geometry/2D_two_throat_{tag}.opts"
    png = ROOT / f"preprocess/two_throat_{tag}.png"

    cy = 0.5 * Ly
    bumps, gaps = build_walls(args.gap_narrow, args.gap_wide, args.orientation)

    xs = np.linspace(0, Lx, 30000)
    yb = wall_bot(xs, bumps)
    gap_prof = Ly - 2.0 * yb
    slope = float(np.max(np.abs(np.gradient(yb, xs))))
    assert slope <= SLOPE_BUDGET, \
        (f"wall slope {math.degrees(math.atan(slope)):.0f} deg exceeds the "
         f"{math.degrees(math.atan(SLOPE_BUDGET)):.0f} deg budget -- widen THROAT_R "
         "or soften ROUGH_HR")
    assert gap_prof.min() > 0.0, "walls touch"

    # ---- Plugs: one per throat, contact at +-contact_frac*THROAT_R ----
    a_off = args.contact_frac * THROAT_R
    bridges, report = [], []
    for i, xt in enumerate(X_THROAT):
        xcL, rL = meniscus(xt - a_off, bumps, cy)
        xcR, rR = meniscus(xt + a_off, bumps, cy)
        assert abs(xcR - xcL) < rL + rR, \
            (f"throat {i}: meniscus discs do not overlap -- plug would be empty; "
             "raise --contact-frac")
        # THE GATE: 90 degrees at all four wall contacts of this plug.
        for xc, rho, xp, lab in ((xcL, rL, xt - a_off, "left"),
                                 (xcR, rR, xt + a_off, "right")):
            for top in (False, True):
                ang = contact_angle_deg(xc, rho, xp, bumps, cy, top)
                assert abs(ang - 90.0) < 1e-6, \
                    (f"throat {i} {lab} meniscus x {'top' if top else 'bottom'} wall: "
                     f"contact angle {ang:.4f} deg, not 90 -- the meniscus "
                     "construction is wrong")
        bridges.append({"cxL": xcL, "rL": rL, "cxR": xcR, "rR": rR, "cy": cy})

        # curvature seen by the ice: convex (ice inside both discs)
        kappa = 0.5 * (1.0 / rL + 1.0 / rR)
        # ice area by quadrature on the actual clipped shape
        gx, gy = np.meshgrid(np.linspace(xt - 2 * THROAT_R, xt + 2 * THROAT_R, 1600),
                            np.linspace(0.0, Ly, 800))
        wb = wall_bot(gx.ravel(), bumps).reshape(gx.shape)
        m = ((np.hypot(gx - xcL, gy - cy) <= rL) & (np.hypot(gx - xcR, gy - cy) <= rR)
             & (gy >= wb) & (gy <= Ly - wb))
        cell = (gx[0, 1] - gx[0, 0]) * (gy[1, 0] - gy[0, 0])
        xin = gx[m]
        span = (float(xin.max()) - float(xin.min())) if xin.size else 0.0
        assert span / (9.2 * EPS) >= MIN_SPAN_BANDS, \
            (f"throat {i}: plug spans {span:.2e} m = {span/(9.2*EPS):.1f} interface "
             f"widths (< {MIN_SPAN_BANDS}); it has no flat core and will dissolve "
             "rather than exchange -- raise --contact-frac")
        assert EPS / min(rL, rR) < 0.05, \
            f"throat {i}: eps/rho = {EPS/min(rL,rR):.1%} (meniscus under-resolved)"
        report.append({"x": xt, "gap": gaps[i], "rL": rL, "rR": rR, "kappa": kappa,
                       "area": float(m.sum()) * cell, "span": span})

    Nx = math.ceil(Lx * math.sqrt(2) / EPS)
    Ny = math.ceil(Ly * math.sqrt(2) / EPS)

    print(f"two-throat pore channel  Lx={Lx:.4e} Ly={Ly:.4e}  Nx={Nx} Ny={Ny}  "
          f"nodes {(Nx+P)*(Ny+P)/1e6:.2f}M")
    print(f"  orientation: {args.orientation}   {N_ROUGH} roughness bumps/wall "
          f"(mirrored) + 2 throats;  wall max slope "
          f"{math.degrees(math.atan(slope)):.0f} deg (budget "
          f"{math.degrees(math.atan(SLOPE_BUDGET)):.0f})")
    print(f"  gap profile: min {gap_prof.min():.3e}  max {gap_prof.max():.3e} m")
    for i, r in enumerate(report):
        print(f"  throat {i} at x={r['x']:.3e}: gap {r['gap']:.3e} m, meniscus rho "
              f"{r['rL']:.3e}/{r['rR']:.3e}  (CONVEX, contact 90.000000 deg verified)")
        print(f"      d0*kappa {D0*r['kappa']:+.3e}   ice area {r['area']:.4e} m^2   "
              f"plug spans {r['span']:.3e} m ({r['span']/(9.2*EPS):.1f} interface widths)")
    src = 0 if report[0]["kappa"] > report[1]["kappa"] else 1
    drive = D0 * abs(report[0]["kappa"] - report[1]["kappa"])
    print(f"\n  DRIVE |d0*dkappa| = {drive:.3e}")
    print(f"  Convex menisci -> HIGHER curvature = HIGHER equilibrium vapour density,")
    print(f"  so throat {src} (gap {report[src]['gap']:.1e}, "
          f"{'NARROW' if report[src]['gap'] < report[1-src]['gap'] else 'WIDE'}) is the "
          f"SOURCE and throat {1-src} is the SINK.")
    sep = abs(X_THROAT[1] - X_THROAT[0])
    print(f"\n  Thermal comparison over the {sep:.1e} m throat separation:")
    for g in (3.0, 0.3):
        dT = g * sep
        print(f"    dT/dx = {g:>4} K/m -> drho/rho = {(LAT/(RV*TK*TK))*dT:.2e} "
              f"({(LAT/(RV*TK*TK))*dT/drive:.1f}x the geometric drive)")

    # ---- Mesh ----
    def bstr(b):
        return ";".join(f"{c:.6e},{R:.6e},{h:.6e}" for c, R, h in b)
    cmd = [sys.executable,
           str(ROOT / "preprocess/build_geometry_multi_grain.py"),
           "--bumps", bstr(bumps), "--top-bumps", bstr(bumps),
           "--Lx", f"{Lx}", "--Ly", f"{Ly}", "--Nx", f"{Nx}", "--Ny", f"{Ny}",
           "--P", str(P), "--C", str(C), "--out", str(dat),
           "--plot", str(png.with_name(png.stem + "_mesh.png")), "--vtk", "/dev/null"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    print(r.stdout[-200:])
    if r.returncode != 0:
        print("MESH BUILD FAILED:\n", r.stderr[-2000:]); sys.exit(1)

    preview(bumps, bridges, png,
            f"two-throat pore channel ({args.orientation}) — gaps "
            f"{gaps[0]:.1e} / {gaps[1]:.1e} m, throat plugs at 90 deg; "
            f"source = throat {src} ({'narrow' if report[src]['gap'] < report[1-src]['gap'] else 'wide'})")

    def arr(v):
        return ",".join(f"{z:.6e}" for z in v)
    rel = dat.relative_to(ROOT)
    with open(opts, "w") as f:
        f.write(f"""# =============================================================================
# geometry/{opts.name} — pore channel with TWO throats of different aperture,
# each plugged by a frozen capillary bridge.
# Generated by preprocess/build_geometry_two_throat.py (deterministic).
# Study: studies/icy_regolith/implicit_pore_domain/ (Effort 1).
#
# THE QUESTION: given two ice bridges in throats of different aperture, does ice
# redistribute between them, and which throat wins? The network version of the
# wedge result (a single deposit migrated -18.75 um/yr toward the narrow end).
#
# PAIR WITH A CLOSED (Neumann-vapour) EXPERIMENT — tgrad_T-20_h1.00_90d_G0 for the
# isothermal case, tgrad_T-20_h1.00_90d_G0.3 for the gradient. Do NOT add
# -flag_BC_rhovfix here: each bridge plugs its throat, so the domain is three
# pockets, and the MIDDLE one already connects the two bridges' inner menisci.
# That is the transport path, it needs no reservoir, and keeping the system
# closed makes total ice conservation a working correctness gate.
#
# GEOMETRY. A full rough pore channel: {N_ROUGH} C-infinity roughness bumps per wall
# drawn from a lunar-regolith GSD, MIRRORED top/bottom, with two wide throat
# bumps (half-width {THROAT_R:.1e} m) whose heights are solved so each throat's gap
# hits its target exactly on top of that roughness. Aperture is the only
# deliberate difference between the two throats:
#   throat 0  x={report[0]['x']:.3e}  gap {report[0]['gap']:.3e} m
#   throat 1  x={report[1]['x']:.3e}  gap {report[1]['gap']:.3e} m
#   gap ranges {gap_prof.min():.3e}..{gap_prof.max():.3e} m; wall max slope {math.degrees(math.atan(slope)):.0f} deg
#
# The two walls are mirrored about y=Ly/2 on purpose: it puts the channel axis on
# a straight line, which is what lets the 90-degree arc centre sit on the axis and
# makes the contact exact. Independently rough walls would need the arc centre at
# the intersection of the two wall tangents with both contacts co-solved -- and
# would trade an exactly-90-degree IC for an approximate one, which every wedge
# run showed is the thing not to approximate.
#
# ICE. Each plug is the INTERSECTION of two axis-centred meniscus discs
# (-bridge_cxL/-bridge_rL/-bridge_cxR/-bridge_rR/-bridge_cy), so both menisci are
# CONVEX and the plug is BARREL-shaped -- wider at the axis than at its wall
# contacts. That is what a 90-degree contact angle forces in a converging-
# diverging channel; the hourglass of a wetting capillary bridge needs theta < 90,
# which this model's Neumann wall BC cannot produce. Contact is exactly 90 deg at
# all four contacts per plug, so there is no relaxation transient to confuse with
# transport:
#   throat 0  d0*kappa {D0*report[0]['kappa']:+.3e}   ice area {report[0]['area']:.4e} m^2
#   throat 1  d0*kappa {D0*report[1]['kappa']:+.3e}   ice area {report[1]['area']:.4e} m^2
#
# PREDICTION: throat {src} (the {'NARROW' if report[src]['gap'] < report[1-src]['gap'] else 'WIDE'} one) is the SOURCE, throat {1-src} the SINK.
# Both menisci are CONVEX, so higher curvature means HIGHER equilibrium vapour
# density: the tighter throat holds its ice at higher kappa and loses it. This is
# Ostwald ripening with aperture playing the role of grain size -- NOT the
# capillary-condensation sign, which would need theta < 90 and an hourglass
# meniscus that a Neumann wall BC cannot produce.
# Drive |d0*dkappa| = {drive:.3e}, vs 1.7e-06 for the wedge run that moved
# measurably at -18.75 um/yr.
#
# 90 DEGREES REQUIRES THE ARC CENTRE ON THE WALL'S TANGENT, not its normal --
# the normal gives TANGENCY, a 0-degree contact, and that error silently
# produced a plausible-looking hourglass shape with the direction reversed. The
# builder now asserts the measured contact angle at all four contacts per plug.
#
# eps={EPS:.4e} is the comp_eps.py (Kaempfer&Plapp) value for T={T0_C:g}C,
# alpha_c={ALPHA_C:g}. RECOMPUTE eps for any other temperature.
# =============================================================================
# DOF_GRID: {Nx+P} {Ny+P}
-geom_file {rel}
-p {P}
-C {C}
-ic_type multi_grains
-dim 2
-Lx {Lx:.6e}
-Ly {Ly:.6e}
-Lz 0
-sed_grain_x {arr([b[0] for b in bumps])}
-sed_grain_R {arr([b[1] for b in bumps])}
-sed_grain_h {arr([b[2] for b in bumps])}
-top_grain_x {arr([b[0] for b in bumps])}
-top_grain_R {arr([b[1] for b in bumps])}
-top_grain_h {arr([b[2] for b in bumps])}
-bridge_cxL {arr([b['cxL'] for b in bridges])}
-bridge_rL  {arr([b['rL'] for b in bridges])}
-bridge_cxR {arr([b['cxR'] for b in bridges])}
-bridge_rR  {arr([b['rR'] for b in bridges])}
-bridge_cy  {arr([b['cy'] for b in bridges])}
-delt_t 1.0e-4
-eps {EPS:.4e}
-eps_valid_temp {T0_C:g}   # C: temperature eps/mesh were sized for; solver ABORTS if -temp differs (override: -eps_temp_override 1)
-periodic 0
""")
    print(f"\nwrote {opts}\nwrote {dat}")


if __name__ == "__main__":
    main()
