#!/usr/bin/env python3
"""build_geometry_ripening.py — two-sided bumpy domain with semicircular ice
caps over the bumps, a very large grain at one end, and a squared ice block at
the other. For an Ostwald-ripening study.

Derived from 2D_random_bumpy_floor_molaroscale.opts:
  * DOMAIN 75% as long (Lx *= 0.75), same height.
  * BUMPS on BOTH edges (floor + ceiling). CENTRAL bumps are amplified into a
    CHANNEL PINCH (height x2, half-width widened in step so the floor slope /
    mesh shear stays ~42deg, matching the un-pinched run). Different seeds
    top/bottom -> similar but not identical.
  * ICE CAPS as NON-TOUCHING semicircles ON each boundary, over the FLANK bumps
    only -- each is checked against the local bump so none is buried, which
    leaves the tall central pinch bare. One "winner" cap per boundary.
  * A VERY LARGE ice grain centred at (Lx, Ly/2) -- a semicircle at the RIGHT
    edge.
  * A squared ICE BLOCK (flat top, no curvature; -ice_flat) at the LEFT edge,
    spanning the FULL domain height.

Loose eps = 8.584e-7. The generator ASSERTS every constraint before writing:
floor+ceiling clearance, no buried cap, eps/R < 5%, and NO overlaps at all.

Writes the .dat mesh, the .opts, and a preview PNG. Deterministic (seeded).
"""

import math
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "preprocess"))
from build_geometry_multi_grain import generate_random_bumps  # noqa: E402

# ---------------------------------------------------------------------------
Lx = 0.75 * 1.3467e-3
Ly = 2.6933e-4
EPS = 8.5840e-7
P, C = 2, 1

N_BUMPS = 18
BUMP_R = (3.3667e-5, 4.7133e-5)
BUMP_HR = (0.30, 0.56)
BUMP_OV = (0.20, 0.32)
SEED_BOT_BUMP, SEED_TOP_BUMP = 0, 7

# Central amplification -> a channel constriction in the middle. Each bump's
# height is scaled by up to CENTRAL_H_MULT and its half-width by up to
# CENTRAL_R_MULT, tapering to 1x at the flanks via a Gaussian of width
# CENTRAL_SIG. Doubling the height alone steepens the floor to ~59deg (skew);
# widening the central bumps in step keeps the slope near the flank value
# (~45deg) while still narrowing the channel -- a broad pinch, not a spike.
CENTRAL_H_MULT = 2.0            # central bump HEIGHT doubled (the pinch)
CENTRAL_R_MULT = 2.8            # widen central bumps in step so the floor slope
                                # (mesh shear) stays ~42deg -- the same as the
                                # clean current run, not the 59deg of height-only
CENTRAL_SIG_FRAC = 0.14         # Gaussian sigma as a fraction of Lx (flanks near
                                # normal height so caps are not buried there)
MIN_CHANNEL_GAP_FRAC = 0.30     # assert the pinch leaves >= this fraction of Ly

R_CAP = (2.6e-5, 4.0e-5)        # cap radii (checked against the LOCAL bump so
R_WIN = (5.0e-5, 6.5e-5)        # none is buried; central pinch is left bare)
SEED_BOT_CAP, SEED_TOP_CAP = 1, 4
MIN_GAP = 1.8e-5                # min edge-to-edge gap between caps (non-touching)
BURY_MARGIN = 0.6e-5           # cap must clear the local bump by this much

# End features
R_BIG = 1.0e-4                  # very large grain, semicircle at (Lx, Ly/2)
# Squared ice block, LEFT, full height. Left face pushed PAST x=0 so the block
# is flush with the domain edge (no air gap): span [XF-RF, XF+RF] = [-3e-5, 1.3e-4].
BLOCK_XF, BLOCK_RF = 5.0e-5, 8.0e-5


def place_caps(rng, y, x_lo, x_hi, local_h, tag):
    """NON-TOUCHING caps across [x_lo,x_hi], MIN_GAP apart, sizes varied. A cap
    is placed only where it clears the LOCAL bump (R > local_h(x)+BURY_MARGIN),
    so the tall central pinch is left bare instead of burying caps there.
    Promotes one interior cap to a winner where it fits. Returns [x,y,R,tag]."""
    caps = []
    cursor = x_lo
    while True:
        R = float(rng.uniform(*R_CAP))
        xc = cursor + R
        if xc + R > x_hi:
            break
        if R > local_h(xc) + BURY_MARGIN:          # clears the local bump -> place
            caps.append([xc, y, R, None])
            cursor = xc + R + float(rng.uniform(MIN_GAP, MIN_GAP + 2.2e-5))
        else:                                       # buried (tall bump) -> skip past
            cursor = xc + R

    def slack(j):  # largest radius cap j could take without touching a neighbour
        s = 1e9
        if j > 0:
            s = min(s, caps[j][0] - (caps[j - 1][0] + caps[j - 1][2]))
        if j < len(caps) - 1:
            s = min(s, (caps[j + 1][0] - caps[j + 1][2]) - caps[j][0])
        return s

    # winner: interior cap with the most room, if a winner radius fits (a bigger
    # R only clears burial further, so no extra check needed)
    cand = [(slack(j), j) for j in range(1, len(caps) - 1) if slack(j) >= R_WIN[0]]
    if cand:
        _, j = max(cand)
        caps[j][2] = float(rng.uniform(R_WIN[0], min(R_WIN[1], slack(j))))
    return caps


def overlaps(a, b):
    (xa, ya, Ra, _), (xb, yb, Rb, _) = a, b
    return math.hypot(xa - xb, ya - yb) < (Ra + Rb)


def amplify_central(bumps):
    """Scale each bump's height (x CENTRAL_H_MULT) and half-width (x
    CENTRAL_R_MULT) near the domain centre, tapering to 1x at the flanks, to
    form the channel constriction. Widening in step with the height keeps the
    floor slope (mesh skew) from blowing up."""
    sig = CENTRAL_SIG_FRAC * Lx
    out = []
    for cx, R, h in bumps:
        env = math.exp(-((cx - 0.5 * Lx) / sig) ** 2)
        out.append((cx, R * (1 + (CENTRAL_R_MULT - 1) * env),
                    h * (1 + (CENTRAL_H_MULT - 1) * env)))
    return out


def bump_field(bumps, xq):
    """Summed bump height at scalar/array xq (matches SedimentBumpField)."""
    xq = np.atleast_1d(np.asarray(xq, float))
    f = np.zeros_like(xq)
    for cx, R, h in bumps:
        t = (xq - cx) / R
        m = np.abs(t) < 1
        f[m] += h * np.exp(1 - 1.0 / (1 - t[m] ** 2))
    return f


def main():
    bot_b = amplify_central(generate_random_bumps(
        np.random.default_rng(SEED_BOT_BUMP), Lx, N_BUMPS, BUMP_R, BUMP_HR, BUMP_OV))
    top_b = amplify_central(generate_random_bumps(
        np.random.default_rng(SEED_TOP_BUMP), Lx, N_BUMPS, BUMP_R, BUMP_HR, BUMP_OV))

    # Local bump-height lookups so caps skip the tall central pinch.
    floor_h = lambda x: float(bump_field(bot_b, x)[0])
    ceil_h = lambda x: float(bump_field(top_b, x)[0])

    x_lo = BLOCK_XF + BLOCK_RF + 4.0e-5
    x_hi = Lx - R_BIG - 4.0e-5
    bot_c = place_caps(np.random.default_rng(SEED_BOT_CAP), 0.0, x_lo, x_hi, floor_h, "bot")
    top_c = place_caps(np.random.default_rng(SEED_TOP_CAP), Ly, x_lo, x_hi, ceil_h, "top")

    big = [Lx, 0.5 * Ly, R_BIG, "big"]          # right-edge semicircle
    BLOCK_H = Ly + 1.0e-5                        # block spans the full height
    caps = bot_c + top_c + [big]

    # ---- Constraints + mesh-skew report ----
    xs = np.linspace(0, Lx, 8000)
    fb, ft = bump_field(bot_b, xs), bump_field(top_b, xs)
    gap = Ly - (fb + ft)
    min_gap, x_pinch = gap.min(), xs[np.argmin(gap)]
    max_slope = float(np.max(np.abs(np.gradient(fb, xs))))
    Rmin = min(c[2] for c in caps)
    assert min_gap > MIN_CHANNEL_GAP_FRAC * Ly, \
        f"channel pinch too tight: {min_gap/Ly:.0%} Ly < {MIN_CHANNEL_GAP_FRAC:.0%}"
    assert EPS / Rmin < 0.05, f"eps/R = {EPS/Rmin:.1%}"
    for c in bot_c:
        assert c[2] > floor_h(c[0]) + 0.0, "a bottom cap is buried"
    for c in top_c:
        assert c[2] > ceil_h(c[0]) + 0.0, "a top cap is buried"
    bad = [(i, k) for i in range(len(caps)) for k in range(i + 1, len(caps))
           if overlaps(caps[i], caps[k])]
    assert not bad, f"unintended overlaps: {bad}"

    Nx = math.ceil(Lx * math.sqrt(2) / EPS)
    Ny = math.ceil(Ly * math.sqrt(2) / EPS)
    print(f"Lx={Lx:.6e} Ly={Ly:.6e}  Nx={Nx} Ny={Ny}  nodes {(Nx+P)*(Ny+P)/1e6:.2f}M")
    print(f"CHANNEL PINCH: min gap {min_gap:.3e} ({min_gap/Ly*100:.0f}% Ly) at "
          f"x={x_pinch/Lx:.2f}Lx;  max floor slope {max_slope:.2f} "
          f"({math.degrees(math.atan(max_slope)):.0f} deg mesh shear)")
    print(f"caps: {len(bot_c)} bottom + {len(top_c)} top (non-touching, off the pinch)")
    print(f"  + very large grain (R={R_BIG:.2e}, right semicircle) + full-height "
          f"flat block (span [{BLOCK_XF-BLOCK_RF:.2e},{BLOCK_XF+BLOCK_RF:.2e}], flush left)")
    allR = [c[2] for c in caps]
    print(f"  R {min(allR):.3e}..{max(allR):.3e} ({max(allR)/min(allR):.1f}x); "
          f"eps/R_smallest {EPS/Rmin*100:.1f}%")

    # ---- Mesh ----
    def bstr(b): return ";".join(f"{cx:.6e},{R:.6e},{h:.6e}" for cx, R, h in b)
    dat = ROOT / "inputs/geometry/ripening_two_sided.dat"
    cmd = [str(ROOT / "venv_pf311/bin/python3"),
           str(ROOT / "preprocess/build_geometry_multi_grain.py"),
           "--bumps", bstr(bot_b), "--top-bumps", bstr(top_b),
           "--Lx", f"{Lx}", "--Ly", f"{Ly}", "--Nx", f"{Nx}", "--Ny", f"{Ny}",
           "--P", str(P), "--C", str(C), "--out", str(dat),
           "--plot", str(ROOT / "preprocess/ripening_two_sided.png")]
    r = subprocess.run(cmd, capture_output=True, text=True)
    print(r.stdout[-400:])
    if r.returncode != 0:
        print("MESH BUILD FAILED:\n", r.stderr[-2000:]); sys.exit(1)

    # ---- Opts ----
    def arr(v): return ",".join(f"{z:.6e}" for z in v)
    gx = [c[0] for c in caps]; gy = [c[1] for c in caps]; gR = [c[2] for c in caps]
    opts = ROOT / "inputs/geometry/2D_ripening_two_sided.opts"
    with open(opts, "w") as f:
        f.write(f"""# =============================================================================
# geometry/2D_ripening_two_sided.opts — two-sided bumpy domain, semicircular ice
# caps over the bumps, a very large grain (right), a squared ice block (left).
# Generated by preprocess/build_geometry_ripening.py (deterministic; re-run to
# regenerate). Ostwald-ripening study.
#
# Domain 75% of the molaroscale length; bumps on BOTH edges. Central bumps are
# amplified (x{CENTRAL_H_MULT:g} height, x{CENTRAL_R_MULT:g} half-width) into a
# CHANNEL PINCH: min gap {min_gap/Ly*100:.0f}% Ly at x={x_pinch/Lx:.2f}Lx, max
# floor slope {math.degrees(math.atan(max_slope)):.0f} deg (mesh shear). Ice caps
# sit ON each boundary as NON-TOUCHING semicircles over the flank bumps (the
# pinch is left bare so caps are not buried). A very large grain (R={R_BIG:.2e})
# is a semicircle at the right edge (Lx, Ly/2); a full-height flat ice block
# (-ice_flat, no curvature) is flush with the left edge. Loose eps = {EPS:.4e},
# eps/R_smallest = {EPS/Rmin*100:.1f}%.
# Pair with an experiment that sets kinetics, e.g. 30day_T-20_h1.00_arrh.opts.
# =============================================================================
# DOF_GRID: {Nx+P} {Ny+P}
-geom_file inputs/geometry/ripening_two_sided.dat
-p {P}
-C {C}
-ic_type multi_grains
-dim 2
-Lx {Lx:.6e}
-Ly {Ly:.6e}
-Lz 0
-sed_grain_x {arr([b[0] for b in bot_b])}
-sed_grain_R {arr([b[1] for b in bot_b])}
-sed_grain_h {arr([b[2] for b in bot_b])}
-top_grain_x {arr([b[0] for b in top_b])}
-top_grain_R {arr([b[1] for b in top_b])}
-top_grain_h {arr([b[2] for b in top_b])}
-ice_grain_cx {arr(gx)}
-ice_grain_cy {arr(gy)}
-ice_grain_R  {arr(gR)}
-ice_grain_ax {arr(gR)}
-ice_grain_ay {arr(gR)}
-ice_flat_x {BLOCK_XF:.6e}
-ice_flat_R {BLOCK_RF:.6e}
-ice_flat_height {BLOCK_H:.6e}
-delt_t 1.0e-4
-eps {EPS:.4e}
-periodic 0
""")
    print(f"\nwrote {opts}\nwrote {dat}")


if __name__ == "__main__":
    main()
