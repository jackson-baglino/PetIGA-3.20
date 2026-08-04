#!/usr/bin/env python3
"""build_geometry_wedge.py — a tapered (wedge) pore channel with perfectly flat
walls, and a single ice band bridging it, for the icy-lunar-regolith study
(studies/icy_regolith/implicit_pore_domain/, Effort 1).

THE QUESTION. Does pore taper ALONE redistribute ice? Every other geometry in
this study answers a question about rough, disordered walls under a thermal
gradient. This one strips all of that away: two straight walls, no gradient, no
second grain. If the ice moves, the geometry moved it.

THE PREDICTION: it migrates toward the NARROW end. The walls are two rays from
a virtual apex off the left edge of the domain, and a circular arc centred on
that apex is perpendicular to every ray from it -- so at the model's natural
90-degree (Neumann) contact angle, both menisci of a bridging band are
apex-centred arcs whose radius of curvature is simply their distance from the
apex. The inner meniscus has the ice on the OUTSIDE of its circle, so it is
concave (kappa = -1/r1); the outer one is convex (kappa = +1/r2). The concave
face has the lower equilibrium vapour density, so vapour flows inward and the
band walks toward the apex.

The energy accounting agrees, independently of any flux argument. At 90 degrees
Young's equation gives gamma_wall-ice = gamma_wall-vapour, so sliding the band
along the wall costs no wall energy and only the two menisci count:

    interfacial length  L = 2*alpha*(r1 + r2)
    ice area            A = alpha*(r2^2 - r1^2) = alpha*(r2-r1)*(r2+r1)

Substituting S = r1 + r2 gives L = 2*alpha*S with r2 - r1 = A/(alpha*S), so at
fixed area L falls monotonically as S shrinks -- i.e. toward the apex.
Curvature RISES toward the apex while interfacial AREA falls, and it is area
that carries the energy. Same sign as two familiar cases: sintering necks grow
rather than pinch off, and capillary condensation fills the smallest pores
first.

This is the analytic clean-room version of the 'throat_bridge' hypothesis in
the study README. The run is the TEST of that prediction, not a demonstration
of it -- a drift toward the wide end would be a real finding, or a sign error.

TWO ICE SHAPES (--ice-shape), and the trade-off between them:

  lens (default) — an oversized DISC clipped by the two walls: a grain larger
      than the channel, cut off top and bottom. TOTALLY CONVEX; both menisci
      are arcs of the same circle. This is what an ice grain grown in a pore
      actually looks like. Uses the ordinary -ice_grain_* circle IC; the walls
      clip it because the domain simply ends there.

  band — the annulus about the wedge apex. The EQUILIBRIUM shape: an
      apex-centred arc is perpendicular to every ray from the apex, hence to
      both walls, so this sits at exactly 90 degrees with no relaxation
      transient. Its inner meniscus is concave. Uses -wedge_band_r1/r2.

CONVEXITY AND 90 DEGREES ARE MUTUALLY EXCLUSIVE HERE. A clipped disc meets the
wall at arccos(h/R), where h is the PERPENDICULAR distance from the centre to
the wall. Driving that to 90 degrees needs R/h -> infinity, and the lens is
2*sqrt(R^2 - h^2) long, so a convex shape at 90 degrees is infinitely long. In
a domain of finite length the lens is therefore always well below 90 degrees
(this script prints the angle it achieved).

CONSEQUENCE FOR THE LENS RUN: the shape will NOT stay convex. The wall BC is
natural Neumann, so the equilibrium contact angle is 90 degrees, and the ice
relaxes toward the band shape -- meaning the INNER (left) meniscus goes concave
early in the run. That is expected, and is itself a useful first check that the
run is behaving. The sustained migration rate quoted below is the one set by
the relaxed shape, not by the initial lens.

Unlike the bump walls, a straight wall is represented EXACTLY by the mesh:
B-splines reproduce linear functions at their Greville abscissae. The affine
wall baseline is carried by -wall_{bot,top}_{y0,slope}, mirrored by
--bot-y0/--bot-slope/--top-y0/--top-slope in build_geometry_multi_grain.py.

Usage (from the project root):
    python3 preprocess/build_geometry_wedge.py
    python3 preprocess/build_geometry_wedge.py --w-throat 3.0e-5 --Lx 2.5e-4 \\
        --tag steep
"""

import argparse
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent

# --------------------------------------------------------------------------- #
# Domain defaults — the 14-degree wedge (see the driving-force note below)
# --------------------------------------------------------------------------- #
Ly = 2.0e-4          # bounding-box height [m]; the mouth exactly fills it
Lx = 3.0e-4          # domain length [m]
W_THROAT = 5.0e-5    # channel width at x=0   (narrow end) [m]
W_MOUTH = 2.0e-4     # channel width at x=Lx  (wide end)   [m]

ICE_X = 0.5          # ice centre as a fraction of Lx
BAND_DR = 1.0e-4     # band radial extent [m] (--ice-shape band)
# Lens radius as a multiple of the PERPENDICULAR wall distance h (--ice-shape
# lens). Larger = contact angle closer to 90 deg but a longer lens with less
# room to migrate; the disc must overhang the walls (>1) to be clipped at all.
# 1.5 gives ~48 deg and leaves ~0.4 channel-widths of travel in the default
# domain. The script prints both, so tune against those, not against this number.
LENS_R_OVER_H = 1.5

# eps is the comp_eps.py (Kaempfer & Plapp) value for T=-20 C, alpha_c=1.341e-2.
# The binding constraint is the temperature-dependent kinetic bound, so it is
# the same value every other -20 C geometry in this study uses. Verified for
# this domain:
#   python3 preprocess/comp_eps.py --Lx 3.0e-4 --Ly 2.0e-4 --Rave 6.25e-5 \
#           --T0 -20 --alpha 1.341e-2      ->  eps 8.5840e-07, Nx 495, Ny 330
T0_C = -20.0
ALPHA_C = 1.341e-2
EPS = 8.5840e-7
P, C = 2, 1

# Driving-force reference: the two-grain ripening pair (50 um vs 150 um) that
# showed clear change in ~30 days. d0 is the capillary length the -20 C
# experiments set via -d0_sub0.
D0 = 1.0166e-9
RIPENING_REF = D0 * (1.0 / 5.0e-5 - 1.0 / 1.5e-4)

# Keep the band clear of the Dirichlet-T end walls, and leave it somewhere to go.
EDGE_MARGIN = 2.0e-5


def preview(y_bot, y_top, ice_mask, lx, ly, fname, title):
    """Render the wedge + the ice, clipped to the channel exactly as the solver
    clips it (the domain simply ends at the walls), so the shape on screen is
    the shape that gets solved."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x = np.linspace(0, lx, 2000)
    fig, ax = plt.subplots(figsize=(11, 11 * ly / lx + 1.4))
    ax.fill_between(x, 0, y_bot(x), color="0.55", zorder=1)
    ax.fill_between(x, y_top(x), ly, color="0.55", zorder=1)
    ax.plot(x, y_bot(x), "k-", lw=1.2, zorder=3)
    ax.plot(x, y_top(x), "k-", lw=1.2, zorder=3)

    gx, gy = np.meshgrid(np.linspace(0, lx, 1200), np.linspace(0, ly, 700))
    inside = ice_mask(gx, gy) & (gy >= y_bot(gx)) & (gy <= y_top(gx))
    ax.contourf(gx, gy, inside.astype(float), levels=[0.5, 1.5],
                colors=["#66b3ff"], zorder=2)
    ax.contour(gx, gy, inside.astype(float), levels=[0.5],
               colors=["#1f6fd0"], linewidths=1.0, zorder=3)

    ax.set_xlim(0, lx)
    ax.set_ylim(0, ly)
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    fig.tight_layout()
    fig.savefig(fname, dpi=130)
    plt.close(fig)
    print(f"wrote {fname}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--Lx", type=float, default=Lx, help="domain length [m]")
    ap.add_argument("--Ly", type=float, default=Ly, help="bounding-box height [m]")
    ap.add_argument("--w-throat", type=float, default=W_THROAT,
                    help="channel width at the narrow (left) end [m]")
    ap.add_argument("--w-mouth", type=float, default=W_MOUTH,
                    help="channel width at the wide (right) end [m]")
    ap.add_argument("--ice-shape", choices=("lens", "band"), default="lens",
                    help="lens = oversized disc clipped by the walls, totally convex "
                         "(default); band = apex-centred annulus, the 90-degree "
                         "equilibrium shape")
    ap.add_argument("--ice-x", type=float, default=ICE_X,
                    help="ice centre as a fraction of Lx (default 0.5)")
    ap.add_argument("--lens-R-over-h", type=float, default=LENS_R_OVER_H,
                    help="lens radius / perpendicular wall distance; >1 so the disc "
                         "overhangs the walls. Larger = contact angle nearer 90 deg "
                         "but a longer lens with less room to migrate")
    ap.add_argument("--band-dr", type=float, default=BAND_DR,
                    help="band radial extent [m] (--ice-shape band only)")
    ap.add_argument("--tag", default=None,
                    help="suffix for output basenames, e.g. --tag steep -> "
                         "2D_wedge_band_steep.opts / wedge_steep.dat")
    ap.add_argument("--out", default=None, help="override mesh .dat path")
    ap.add_argument("--opts", default=None, help="override .opts path")
    args = ap.parse_args()

    lx, ly = args.Lx, args.Ly
    w_t, w_m = args.w_throat, args.w_mouth
    tag = f"_{args.tag}" if args.tag else ""
    dat = Path(args.out) if args.out else ROOT / f"inputs/geometry/wedge{tag}.dat"
    opts = Path(args.opts) if args.opts else ROOT / f"inputs/geometry/2D_wedge_band{tag}.opts"
    png = ROOT / f"preprocess/wedge{tag}.png"

    assert w_m > w_t > 0, f"need w_mouth > w_throat > 0 (got {w_t:.3e}, {w_m:.3e})"
    assert w_m <= ly, f"mouth {w_m:.3e} exceeds Ly={ly:.3e}"

    # ---- Walls: straight lines, symmetric about y = Ly/2 ----
    # width(x) = w_t + 2*m*x, so the channel opens at rate 2*m.
    m = (w_m - w_t) / (2.0 * lx)
    bot_y0 = 0.5 * (ly - w_t)
    top_y0 = 0.5 * (ly + w_t)
    y_bot = lambda x: bot_y0 - m * x
    y_top = lambda x: top_y0 + m * x
    alpha = math.atan(m)                      # half-angle

    # ---- Virtual apex: where the two walls would meet (outside the domain) ----
    apex_x = -w_t / (2.0 * m)
    apex_y = 0.5 * ly
    assert apex_x < 0.0, "apex must lie outside the domain (left of x=0)"

    xs = np.linspace(0.0, lx, 4000)
    assert (y_top(xs) > y_bot(xs)).all(), "walls cross"
    assert (y_bot(xs) >= -1e-15).all() and (y_top(xs) <= ly + 1e-15).all(), \
        "walls leave the [0, Ly] bounding box"

    # ---- Ice ----
    ice_cx = args.ice_x * lx
    ice_cy = 0.5 * ly                         # the wedge axis
    r_c = ice_cx - apex_x                     # centre distance from the apex
    w_at_ice = w_t + 2.0 * m * ice_cx
    # PERPENDICULAR distance from the axis to a wall -- this, not the vertical
    # half-width, is what sets the contact angle of a clipped disc.
    h_perp = 0.5 * w_at_ice * math.cos(alpha)

    if args.ice_shape == "lens":
        R = args.lens_R_over_h * h_perp
        assert R > h_perp, \
            (f"--lens-R-over-h must exceed 1 so the disc overhangs the walls and is "
             f"clipped (got {args.lens_R_over_h})")
        contact_deg = math.degrees(math.acos(h_perp / R))
        ice_x_lo, ice_x_hi = ice_cx - R, ice_cx + R
        half_chord = math.sqrt(R * R - h_perp * h_perp)
        # Clipped-disc area: disc minus the two segments cut off by the walls.
        area = 2.0 * (h_perp * half_chord + R * R * math.asin(h_perp / R))
        # Both menisci are arcs of the SAME circle, so their curvatures are
        # equal and there is no instantaneous driving force. The sustained rate
        # is the one the relaxed (90-degree) shape sets, so quote the band's.
        r1, r2 = r_c - 0.5 * w_at_ice, r_c + 0.5 * w_at_ice
        dkappa = 1.0 / r1 - 1.0 / r2
        assert EPS / R < 0.05, f"eps/R = {EPS/R:.1%} (interface under-resolved)"
    else:
        r1 = r_c - 0.5 * args.band_dr
        r2 = r_c + 0.5 * args.band_dr
        assert r1 > 0.0, f"band inner radius {r1:.3e} <= 0 — band_dr too large"
        R = None
        contact_deg = 90.0
        ice_x_lo, ice_x_hi = r1 + apex_x, r2 + apex_x
        area = alpha * (r2 ** 2 - r1 ** 2)
        dkappa = 1.0 / r1 - 1.0 / r2
        assert EPS / r1 < 0.05, f"eps/r1 = {EPS/r1:.1%} (interface under-resolved)"

    assert ice_x_lo > EDGE_MARGIN, \
        (f"ice reaches x={ice_x_lo:.3e}, inside the {EDGE_MARGIN:.1e} m keep-out at "
         "the Dirichlet-T left wall — shrink it or move it right")
    assert ice_x_hi < lx - EDGE_MARGIN, \
        (f"ice reaches x={ice_x_hi:.3e}, too close to the right wall")

    drive = D0 * dkappa

    Nx = math.ceil(lx * math.sqrt(2) / EPS)
    Ny = math.ceil(ly * math.sqrt(2) / EPS)

    print(f"wedge pore channel  Lx={lx:.4e} Ly={ly:.4e}  Nx={Nx} Ny={Ny}  "
          f"nodes {(Nx+P)*(Ny+P)/1e6:.2f}M")
    print(f"  walls: flat, width {w_t:.3e} -> {w_m:.3e}, half-angle "
          f"{math.degrees(alpha):.2f} deg (slope {m:+.4f})")
    print(f"         y_bot(x) = {bot_y0:.4e} {-m:+.4f}*x")
    print(f"         y_top(x) = {top_y0:.4e} {m:+.4f}*x")
    print(f"  apex:  ({apex_x:.4e}, {apex_y:.4e})  [virtual, outside the domain]")
    if args.ice_shape == "lens":
        print(f"  ice:   'lens' — disc R={R:.4e} centred ({ice_cx:.4e}, {ice_cy:.4e}), "
              f"clipped by both walls")
        print(f"         R/h = {args.lens_R_over_h:.2f} (h_perp={h_perp:.4e}) -> "
              f"contact angle {contact_deg:.1f} deg, TOTALLY CONVEX")
    else:
        print(f"  ice:   'band' — annulus r1={r1:.4e} r2={r2:.4e} about the apex, "
              f"contact angle 90 deg (equilibrium)")
    print(f"         spans x {ice_x_lo:.4e}..{ice_x_hi:.4e}, local channel width "
          f"{w_at_ice:.3e}, area {area:.3e} m^2")
    print(f"         room to migrate inward: {ice_x_lo - EDGE_MARGIN:.3e} m "
          f"({(ice_x_lo - EDGE_MARGIN)/w_at_ice:.2f} channel widths)")
    print(f"  DRIVE: d_kappa = {dkappa:.4g} 1/m -> d0*d_kappa = {drive:.3e} "
          f"({drive/RIPENING_REF:.2f}x the 50/150 um ripening pair)")
    if args.ice_shape == "lens":
        print("         (both lens menisci share one circle, so d_kappa is 0 at t=0; "
              "the figure above is the")
        print("          SUSTAINED rate, set by the relaxed 90-degree shape the lens "
              "decays into)")
    if drive < 0.05 * RIPENING_REF:
        print("  WARNING: driving force is <5% of the ripening reference — expect "
              "no measurable motion. Steepen the wedge (force scales as 4*alpha^2/w).")
    if args.ice_shape == "lens" and contact_deg < 60.0:
        print(f"  NOTE: contact angle {contact_deg:.0f} deg is well below the 90 deg "
              "equilibrium, so expect a visible")
        print("        early relaxation (the inner meniscus will go CONCAVE) before "
              "any migration. Raise")
        print("        --lens-R-over-h for a blunter start, at the cost of travel room.")

    # ---- Build the mesh via the shared igakit template ----
    cmd = [str(ROOT / "venv_pf311/bin/python3"),
           str(ROOT / "preprocess/build_geometry_multi_grain.py"),
           "--bumps", "", "--top-bumps", "",
           "--bot-y0", f"{bot_y0}", "--bot-slope", f"{-m}",
           "--top-y0", f"{top_y0}", "--top-slope", f"{m}",
           "--Lx", f"{lx}", "--Ly", f"{ly}", "--Nx", f"{Nx}", "--Ny", f"{Ny}",
           "--P", str(P), "--C", str(C), "--out", str(dat),
           "--plot", str(png.with_name(png.stem + "_mesh.png")),
           "--vtk", "/dev/null"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    print(r.stdout[-300:])
    if r.returncode != 0:
        print("MESH BUILD FAILED:\n", r.stderr[-2000:])
        sys.exit(1)

    if args.ice_shape == "lens":
        mask = lambda gx, gy: np.hypot(gx - ice_cx, gy - ice_cy) <= R
        shape_txt = f"convex lens, R/h={args.lens_R_over_h:.2f}, contact {contact_deg:.0f} deg"
    else:
        mask = lambda gx, gy: ((np.hypot(gx - apex_x, gy - apex_y) >= r1) &
                               (np.hypot(gx - apex_x, gy - apex_y) <= r2))
        shape_txt = "apex-centred band, contact 90 deg"
    preview(y_bot, y_top, mask, lx, ly, png,
            f"wedge pore — half-angle {math.degrees(alpha):.1f} deg, {shape_txt}, "
            f"at x={ice_cx:.2e} m (sustained drive {drive/RIPENING_REF:.2f}x ripening)")

    rel_dat = dat.relative_to(ROOT)
    if args.ice_shape == "lens":
        ice_note = f"""# ICE — 'lens': a disc of radius {R:.4e} m centred on the wedge axis at
# ({ice_cx:.4e}, {ice_cy:.4e}), LARGER than the channel so both walls cut it.
# The result is TOTALLY CONVEX: both menisci are arcs of the same circle, which
# is what an ice grain grown inside a pore actually looks like. No special IC is
# needed — it is an ordinary -ice_grain_* circle, and the walls clip it because
# the domain simply ends there.
#
# CONTACT ANGLE {contact_deg:.1f} deg, = arccos(h/R) with h={h_perp:.4e} m the
# perpendicular distance from the axis to a wall. This is NOT 90 deg, and it
# cannot be: driving arccos(h/R) to 90 deg needs R/h -> infinity, and the lens is
# 2*sqrt(R^2-h^2) long, so a convex shape at 90 deg is infinitely long.
# Convexity and 90-degree contact are mutually exclusive in a converging channel.
#
# CONSEQUENCE: the shape will NOT stay convex. The wall BC is natural Neumann,
# so equilibrium is 90 deg, and the ice relaxes toward the apex-centred band —
# the inner (left) meniscus goes CONCAVE early in the run. Expected, and a
# useful first check that the run is behaving. Raise --lens-R-over-h for a
# blunter start, at the cost of travel room.
#"""
        ice_block = (f"-ice_grain_cx {ice_cx:.6e}\n"
                     f"-ice_grain_cy {ice_cy:.6e}\n"
                     f"-ice_grain_R  {R:.6e}\n"
                     f"-ice_grain_ax {R:.6e}\n"
                     f"-ice_grain_ay {R:.6e}")
    else:
        ice_note = f"""# ICE — 'band': the annulus r1 <= |X - apex| <= r2 about the VIRTUAL apex at
# ({apex_x:.4e}, {apex_y:.4e}), where the two walls would meet. An apex-centred
# arc is perpendicular to every ray from the apex, hence to both walls, so this
# is the shape that already sits at the natural 90-degree contact angle, with no
# relaxation transient. Its inner meniscus is concave; it is not convex.
#"""
        ice_block = (f"-wedge_band_r1 {r1:.6e}\n"
                     f"-wedge_band_r2 {r2:.6e}")

    with open(opts, "w") as f:
        f.write(f"""# =============================================================================
# geometry/{opts.name} — tapered (wedge) pore channel with a bridging ice band.
# Generated by preprocess/build_geometry_wedge.py (deterministic).
# Study: studies/icy_regolith/implicit_pore_domain/ (Effort 1).
#
# THE QUESTION: does pore taper ALONE move ice? Two perfectly flat walls, no
# thermal gradient, no second grain. PAIR WITH A ZERO-GRADIENT EXPERIMENT
# (Tgrad_T-20_G0_90d / Tgrad_T-20_G0_10s) — a gradient would swamp the effect.
#
# GEOMETRY. Straight walls, symmetric about y = Ly/2, opening left to right:
#   width {w_t:.3e} m at x=0  ->  {w_m:.3e} m at x=Lx
#   half-angle {math.degrees(alpha):.2f} deg
#   y_bot(x) = {bot_y0:.6e} {-m:+.6f}*x
#   y_top(x) = {top_y0:.6e} {m:+.6f}*x
# The walls are carried by -wall_{{bot,top}}_{{y0,slope}}, NOT by the bump lists:
# bumps have compact support and cannot express a ramp. A straight wall is
# represented EXACTLY by the mesh (B-splines reproduce linear functions at
# their Greville abscissae), unlike the bump geometries.
#
{ice_note}
# PREDICTION: the ice migrates toward the NARROW (left) end. Once relaxed to
# the 90-degree contact angle its inner meniscus is concave (kappa = -1/r1) and
# its outer convex (+1/r2), so the inner face has the lower equilibrium vapour
# density and vapour flows inward. Sustained driving force
# d0*(1/r1 - 1/r2) = {drive:.3e}, i.e. {drive/RIPENING_REF:.2f}x the 50/150 um ripening
# pair that showed clear change in ~30 days. Expect modest but measurable
# displacement over 90 days, not a dramatic one. If it is too slow, steepen the
# wedge — the force scales as 4*alpha^2/w, NOT with run time.
#
# The counterintuitive part is that the narrow end has the LARGER |kappa|. True,
# but kappa carries a sign and the sign is what sets transport direction. The
# energy accounting agrees with no flux argument at all: at 90 degrees Young
# gives gamma_wall-ice = gamma_wall-vap, so sliding costs no wall energy and only
# the menisci count — L = 2*alpha*(r1+r2) at fixed A = alpha*(r2^2-r1^2), which
# with S = r1+r2 is L = 2*alpha*S and falls monotonically toward the apex.
#
# eps={EPS:.4e} is the comp_eps.py (Kaempfer&Plapp) value for T={T0_C:g}C,
# alpha_c={ALPHA_C:g} — the same kinetic-bound value every other {T0_C:g}C geometry
# here uses. RECOMPUTE eps for any other temperature.
# =============================================================================
# DOF_GRID: {Nx+P} {Ny+P}
-geom_file {rel_dat}
-p {P}
-C {C}
-ic_type multi_grains
-dim 2
-Lx {lx:.6e}
-Ly {ly:.6e}
-Lz 0
-wall_bot_y0 {bot_y0:.6e}
-wall_bot_slope {-m:.6e}
-wall_top_y0 {top_y0:.6e}
-wall_top_slope {m:.6e}
# The apex is emitted for BOTH ice shapes: postprocess/track_wedge_band.py
# reads it to compute the centroid's distance from the apex, which is the
# number the whole run is about. With -wedge_band_r2 <= r1 the band IC is off.
-wedge_apex_x {apex_x:.6e}
-wedge_apex_y {apex_y:.6e}
{ice_block}
-delt_t 1.0e-4
-eps {EPS:.4e}
-eps_valid_temp {T0_C:g}   # C: temperature eps/mesh were sized for; solver ABORTS if -temp differs (override: -eps_temp_override 1)
-periodic 0
""")
    print(f"\nwrote {opts}\nwrote {dat}")
    print(f"eps={EPS:.4e} is the comp_eps.py value for T={T0_C:g}C "
          f"(eps/r1={EPS/r1*100:.1f}%). Valid ONLY at T={T0_C:g}C — pair with a "
          f"{T0_C:g}C ZERO-GRADIENT experiment.")


if __name__ == "__main__":
    main()
