#!/usr/bin/env python3
"""build_geometry_ripening.py — two-sided bumpy domain with semicircular ice
caps over the bumps, a very large grain at one end, and a squared ice block at
the other. For an Ostwald-ripening study.

Derived from 2D_random_bumpy_floor_molaroscale.opts:
  * DOMAIN 75% as long (Lx *= 0.75), same height.
  * BUMPS on BOTH edges (floor + ceiling), magnitude x2, same bump width,
    count ~0.75x. Different seeds top/bottom -> similar but not identical.
  * ICE CAPS as semicircles ON each boundary (cy = 0 / Ly), sitting over the
    bumps. R > tallest bump so none is buried. ALL NON-TOUCHING (>= MIN_GAP
    apart). Wide size spread with one large "winner" cap per boundary.
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
BUMP_HR = (0.30, 0.56)          # x2 magnitude
BUMP_OV = (0.20, 0.32)
SEED_BOT_BUMP, SEED_TOP_BUMP = 0, 7

R_CAP = (2.6e-5, 4.0e-5)        # cap radii (all > tallest bump so not buried)
R_WIN = (5.0e-5, 6.5e-5)        # one "winner" cap per boundary
SEED_BOT_CAP, SEED_TOP_CAP = 1, 4
MIN_GAP = 1.8e-5                # min edge-to-edge gap between caps (non-touching)

# End features
R_BIG = 1.0e-4                  # very large grain, semicircle at (Lx, Ly/2)
BLOCK_XF, BLOCK_RF = 7.0e-5, 6.0e-5     # squared ice block, left; full height


def place_caps(rng, y, x_lo, x_hi, tag):
    """A handful of NON-TOUCHING caps across [x_lo,x_hi], MIN_GAP apart, sizes
    varied; promote one interior cap to a winner where it fits without touching
    a neighbour. Returns [x, y, R, tag]."""
    caps = []
    while True:
        R = float(rng.uniform(*R_CAP))
        if not caps:
            xc = x_lo + R
        else:
            xp, _, Rp, _ = caps[-1]
            xc = xp + Rp + R + float(rng.uniform(MIN_GAP, MIN_GAP + 2.2e-5))
        if xc + R > x_hi:
            break
        caps.append([xc, y, R, None])

    def slack(j):  # largest radius cap j could take without touching a neighbour
        s = 1e9
        if j > 0:
            s = min(s, caps[j][0] - (caps[j - 1][0] + caps[j - 1][2]))
        if j < len(caps) - 1:
            s = min(s, (caps[j + 1][0] - caps[j + 1][2]) - caps[j][0])
        return s

    # winner: the interior cap with the most room, if a winner radius fits
    cand = [(slack(j), j) for j in range(1, len(caps) - 1) if slack(j) >= R_WIN[0]]
    if cand:
        _, j = max(cand)
        caps[j][2] = float(rng.uniform(R_WIN[0], min(R_WIN[1], slack(j))))
    return caps


def overlaps(a, b):
    (xa, ya, Ra, _), (xb, yb, Rb, _) = a, b
    return math.hypot(xa - xb, ya - yb) < (Ra + Rb)


def main():
    bot_b = generate_random_bumps(np.random.default_rng(SEED_BOT_BUMP), Lx,
                                  N_BUMPS, BUMP_R, BUMP_HR, BUMP_OV)
    top_b = generate_random_bumps(np.random.default_rng(SEED_TOP_BUMP), Lx,
                                  N_BUMPS, BUMP_R, BUMP_HR, BUMP_OV)
    maxH = max(max(h for *_, h in bot_b), max(h for *_, h in top_b))

    # Caps live in the middle span; ends reserved for the full-height block
    # (left) and the very large semicircle grain (right, centred on x=Lx).
    x_lo = BLOCK_XF + BLOCK_RF + 4.0e-5
    x_hi = Lx - R_BIG - 4.0e-5
    bot_c = place_caps(np.random.default_rng(SEED_BOT_CAP), 0.0, x_lo, x_hi, "bot")
    top_c = place_caps(np.random.default_rng(SEED_TOP_CAP), Ly, x_lo, x_hi, "top")

    # Very large grain centred exactly at (Lx, Ly/2): a semicircle at the right
    # edge (half in-domain).
    big = [Lx, 0.5 * Ly, R_BIG, "big"]
    BLOCK_H = Ly + 1.0e-5          # block spans the ENTIRE domain height

    caps = bot_c + top_c + [big]

    # ---- Constraints ----
    xs = np.linspace(0, Lx, 6000)
    def field(bumps):
        f = np.zeros_like(xs)
        for cx, R, h in bumps:
            t = (xs - cx) / R
            m = np.abs(t) < 1
            f[m] += h * np.exp(1 - 1.0 / (1 - t[m] ** 2))
        return f
    clearance = Ly - (field(bot_b) + field(top_b)).max()
    Rmin = min(c[2] for c in caps)
    assert clearance > 0, f"floor+ceiling collide ({clearance:.2e})"
    assert min(c[2] for c in bot_c) + 0.0 > maxH, "bottom cap buried"
    assert min(c[2] for c in top_c) > maxH, "top cap buried"
    assert EPS / Rmin < 0.05, f"eps/R = {EPS/Rmin:.1%}"
    # no UNINTENDED overlaps (tagged *-pair partners may overlap)
    bad = []
    for i in range(len(caps)):
        for k in range(i + 1, len(caps)):
            if overlaps(caps[i], caps[k]):
                ti, tk = caps[i][3], caps[k][3]
                if not (ti and tk and ti == tk):
                    bad.append((i, k, caps[i][3], caps[k][3]))
    assert not bad, f"unintended overlaps: {bad}"

    Nx = math.ceil(Lx * math.sqrt(2) / EPS)
    Ny = math.ceil(Ly * math.sqrt(2) / EPS)
    print(f"Lx={Lx:.6e} Ly={Ly:.6e}  Nx={Nx} Ny={Ny}  nodes {(Nx+P)*(Ny+P)/1e6:.2f}M")
    print(f"bumps {len(bot_b)}+{len(top_b)}, maxH {maxH:.3e} ({maxH/Ly*100:.1f}% Ly), "
          f"floor+ceiling gap {clearance/Ly*100:.0f}% Ly")
    print(f"caps: {len(bot_c)} bottom + {len(top_c)} top (all non-touching)")
    print(f"  + 1 very large grain (R={R_BIG:.2e}, right) + 1 flat block "
          f"(xf={BLOCK_XF:.2e}, Rf={BLOCK_RF:.2e}, H={BLOCK_H:.2e}, left)")
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
# Domain 75% of the molaroscale length; bumps on BOTH edges at x2 magnitude
# (seeds top!=bottom). Ice caps sit ON each boundary (cy=0 / Ly) as semicircles
# over the bumps, R > tallest bump ({maxH:.2e}) so none is buried; MOSTLY
# ISOLATED with one touching pair per boundary. A very large grain (R={R_BIG:.2e})
# sits mid-height at the right; a flat-topped squared ice block (-ice_flat,
# no curvature) at the left. Loose eps = {EPS:.4e}, eps/R_smallest = {EPS/Rmin*100:.1f}%.
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
