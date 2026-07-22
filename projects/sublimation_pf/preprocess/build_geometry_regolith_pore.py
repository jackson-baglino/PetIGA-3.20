#!/usr/bin/env python3
"""build_geometry_regolith_pore.py — a two-sided pore channel whose walls are
lined with lunar-regolith grains, for the icy-lunar-regolith study
(studies/icy_regolith/implicit_pore_domain/, Effort 1).

The regolith is modelled *implicitly*: it is not a simulated field but the
deformed top and bottom boundaries of the (single-patch) domain. Ice, vapour,
and temperature evolve in the pore space between the walls. Interior regolith
grains are not representable in a single tensor-product patch — that is Effort 2
(explicit sediment phase).

Derived from build_geometry_ripening.py, with two deliberate changes:

  1. WALLS FROM A REGOLITH GRAIN-SIZE DISTRIBUTION. Wall bumps are drawn with
     half-widths bracketing the lunar-regolith median grain radius (~50 µm;
     docs/material_parameters.md §2.3), rather than the ripening study's tuned
     sizes. A central throat pinch is kept (the defining pore-throat feature)
     and its tightness is a CLI knob.

  2. ICE PLACEMENT IS A SWAPPABLE STRATEGY (--ice-placement). How ice adheres to
     regolith as a function of pore geometry is the open question this study
     probes, so the initial ice configuration is not hard-coded. Strategies:
       flank_caps   — ice caps on the flank grains, throat left bare (baseline).
       throat_bridge— a single ice grain bridging the central throat (tests
                      whether ice preferentially persists/adheres at the pinch),
                      plus non-buried flank caps.
       pore_lining  — many small ice grains lining both walls across the pore
                      bodies (distributed adhesion).
     Add more by registering a function in ICE_STRATEGIES.

Writes inputs/geometry/regolith_pore.dat and inputs/geometry/2D_regolith_pore.opts
(so scripts/Studio/run_permafrost.sh finds them by the usual convention), plus a
preview PNG. Deterministic (seeded). Every geometric constraint is asserted
before anything is written.

Usage (from the project root):
    python3 preprocess/build_geometry_regolith_pore.py --ice-placement flank_caps
    python3 preprocess/build_geometry_regolith_pore.py --ice-placement throat_bridge \\
        --tag throat --seed 3

eps is loose here (reuses the validated ~8.58e-7 at this scale). ALWAYS
recompute it for the actual run temperature/geometry with preprocess/comp_eps.py
before a production run — see the note printed at the end.
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
# Regolith / domain parameters
# --------------------------------------------------------------------------- #
# Lunar regolith: median grain radius ~50 µm, range ~20–100 µm
# (docs/material_parameters.md §2.3). Wall bumps represent grains protruding
# into the pore; half-width R brackets the median, height is tied to R.
R_MED = 5.0e-5                       # median regolith grain radius [m]
BUMP_R = (3.5e-5, 7.0e-5)            # wall-grain half-width range [m] (brackets R_MED)
BUMP_HR = (0.35, 0.60)              # bump height / half-width (bounds curvature ~h/R^2)
BUMP_OV = (0.20, 0.32)             # support overlap fraction (full wall coverage)

Ly = 2.6933e-4                       # channel height [m] (~5 median grains tall)
N_BUMPS = 12                         # wall grains per side
# Domain width so N_BUMPS median grains tile it with typical overlap:
# (n-1)*step ~= Lx, step = 2*R_MED*(1-avg_overlap).
Lx = (N_BUMPS - 1) * (2 * R_MED * (1 - 0.26))

# eps is the comp_eps.py (Kaempfer & Plapp) value FOR T=-20 C, alpha_c=1.341e-2
# — binding constraint is the T-dependent kinetic bound, so it is grain-size
# independent here (all strategies pass eps/R_ave < 5%). Verified: reproduces
# the validated 2D_ripening_two_sided reference run exactly, and pairs with the
# 30day_T-20_h1.00_arrh experiment's beta_sub0=5.9216e5 / d0_sub0=1.0166e-9.
# RECOMPUTE for any other run temperature:
#   python3 preprocess/comp_eps.py --Lx {Lx} --Ly {Ly} --Rave <R_smallest_ice> \
#           --T0 <degC> --alpha <alpha_c>
T0_C = -20.0                         # temperature eps is valid for [deg C]
ALPHA_C = 1.341e-2                   # attachment coefficient used for eps/kinetics
EPS = 8.5840e-7                      # comp_eps.py value at T0_C (NOT arbitrary)
P, C = 2, 1

SEED_BOT_BUMP, SEED_TOP_BUMP = 0, 7  # different seeds -> similar-but-distinct walls

# Central throat pinch (same mechanism as build_geometry_ripening.py):
# amplify central wall grains, tapering to 1x at the flanks.
CENTRAL_H_MULT = 2.0
CENTRAL_R_MULT = 2.8
CENTRAL_SIG_FRAC = 0.14
MIN_CHANNEL_GAP_FRAC = 0.25          # assert the throat leaves >= this fraction of Ly

# Ice-cap sizing (used by flank_caps / throat_bridge flanks)
R_CAP = (2.6e-5, 4.0e-5)
SEED_BOT_CAP, SEED_TOP_CAP = 1, 4
MIN_GAP = 1.8e-5                     # min edge-to-edge gap between caps
BURY_MARGIN = 0.6e-5               # a cap must clear its local bump by this much
EDGE_MARGIN = 4.0e-5               # keep ice this far from the left/right edges


# --------------------------------------------------------------------------- #
# Shared geometry helpers (match SedimentBumpField / SedimentBump in
# src/initial_conditions.c and _bump_field in build_geometry_multi_grain.py)
# --------------------------------------------------------------------------- #
def amplify_central(bumps):
    """Scale central bump height (x CENTRAL_H_MULT) and half-width
    (x CENTRAL_R_MULT), tapering to 1x at the flanks, forming the throat."""
    sig = CENTRAL_SIG_FRAC * Lx
    out = []
    for cx, R, h in bumps:
        env = math.exp(-((cx - 0.5 * Lx) / sig) ** 2)
        out.append((cx, R * (1 + (CENTRAL_R_MULT - 1) * env),
                    h * (1 + (CENTRAL_H_MULT - 1) * env)))
    return out


def bump_field(bumps, xq):
    """Summed C-infinity bump height at xq (matches SedimentBumpField)."""
    xq = np.atleast_1d(np.asarray(xq, float))
    f = np.zeros_like(xq)
    for cx, R, h in bumps:
        t = (xq - cx) / R
        m = np.abs(t) < 1
        f[m] += h * np.exp(1 - 1.0 / (1 - t[m] ** 2))
    return f


def overlaps(a, b):
    (xa, ya, Ra), (xb, yb, Rb) = a[:3], b[:3]
    return math.hypot(xa - xb, ya - yb) < (Ra + Rb)


def preview_ice(bot_b, top_b, ice, fname, title):
    """Render the pore geometry + ice placement: regolith walls (filled solid),
    the open pore channel, and each ice grain (clipped to the pore). Far more
    useful for the adhesion study than the bare control mesh."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle

    x = np.linspace(0, Lx, 4000)
    floor = bump_field(bot_b, x)
    ceil = Ly - bump_field(top_b, x)

    fig, ax = plt.subplots(figsize=(12, 12 * Ly / Lx + 1))
    # regolith solids (below floor, above ceiling)
    ax.fill_between(x, 0, floor, color="0.55", zorder=1)
    ax.fill_between(x, ceil, Ly, color="0.55", zorder=1)
    ax.plot(x, floor, "k-", lw=1, zorder=3)
    ax.plot(x, ceil, "k-", lw=1, zorder=3)
    # ice grains (clipped to the domain box so wall-seated semicircles read right)
    for cx, cy, R in ice:
        ax.add_patch(Circle((cx, cy), R, facecolor="#66b3ff", edgecolor="#1f6fd0",
                            lw=1.0, alpha=0.9, zorder=2))
    ax.set_xlim(0, Lx)
    ax.set_ylim(0, Ly)
    ax.set_aspect("equal")
    ax.set_title(title)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    fig.tight_layout()
    fig.savefig(fname, dpi=130)
    plt.close(fig)
    print(f"wrote {fname}")


# --------------------------------------------------------------------------- #
# Ice-placement strategies. Each returns a list of [cx, cy, R] ice grains
# (cy=0 -> sits on the floor wall, cy=Ly -> ceiling wall, 0<cy<Ly -> interior).
# floor_h(x) / ceil_h(x) give the local wall protrusion so caps are not buried.
# --------------------------------------------------------------------------- #
def _walk_caps(rng, y, x_lo, x_hi, local_h):
    """Non-touching caps across [x_lo, x_hi], each placed only where it clears
    its local wall bump (so tall throat grains are skipped, not buried)."""
    caps, cursor = [], x_lo
    while True:
        R = float(rng.uniform(*R_CAP))
        xc = cursor + R
        if xc + R > x_hi:
            break
        if R > local_h(xc) + BURY_MARGIN:
            caps.append([xc, y, R])
            cursor = xc + R + float(rng.uniform(MIN_GAP, MIN_GAP + 2.2e-5))
        else:
            cursor = xc + R
    return caps


def ice_flank_caps(rng, bot_b, top_b, floor_h, ceil_h, x_lo, x_hi):
    """Baseline: ice caps on flank grains of both walls; throat left bare."""
    caps = _walk_caps(rng.spawn(1)[0], 0.0, x_lo, x_hi, floor_h)
    caps += _walk_caps(rng.spawn(1)[0], Ly, x_lo, x_hi, ceil_h)
    return caps


def ice_throat_bridge(rng, bot_b, top_b, floor_h, ceil_h, x_lo, x_hi):
    """A single ice grain bridging the central throat, plus flank caps. Tests
    whether ice preferentially persists at the constriction."""
    xs = np.linspace(x_lo, x_hi, 4000)
    gap = Ly - (bump_field(bot_b, xs) + bump_field(top_b, xs))
    x_pinch = float(xs[int(np.argmin(gap))])
    gap_min = float(gap.min())
    # Bridge grain centred in the throat, radius ~40% of the open gap so it
    # sits in the pore without immediately overlapping either wall.
    bridge = [x_pinch, 0.5 * Ly, 0.40 * gap_min]
    caps = ice_flank_caps(rng, bot_b, top_b, floor_h, ceil_h, x_lo, x_hi)
    # drop any flank cap that would overlap the bridge
    caps = [c for c in caps if not overlaps(c, bridge)]
    return [bridge] + caps


def _wall_troughs(bumps, x_lo, x_hi, n=3000):
    """x-positions of local minima of the wall protrusion — the valleys
    (reentrant pore corners) between adjacent regolith grains."""
    x = np.linspace(x_lo, x_hi, n)
    h = bump_field(bumps, x)
    return [float(x[i]) for i in range(1, n - 1)
            if h[i] <= h[i - 1] and h[i] < h[i + 1]]


def ice_pore_lining(rng, bot_b, top_b, floor_h, ceil_h, x_lo, x_hi):
    """Small ice grains seated in the wall troughs of both walls — ice
    adhering in the reentrant corners between regolith grains (distributed
    adhesion). Each grain is placed TANGENT to the local wall surface
    (cy = wall_height ± R), so it is never buried."""
    R = 1.8e-5
    ice = []
    for xt in _wall_troughs(bot_b, x_lo, x_hi):
        g = [xt, floor_h(xt) + R, R]
        if not any(overlaps(g, o) for o in ice):
            ice.append(g)
    for xt in _wall_troughs(top_b, x_lo, x_hi):
        g = [xt, Ly - ceil_h(xt) - R, R]
        if not any(overlaps(g, o) for o in ice):
            ice.append(g)
    return ice


ICE_STRATEGIES = {
    "flank_caps": ice_flank_caps,
    "throat_bridge": ice_throat_bridge,
    "pore_lining": ice_pore_lining,
}


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ice-placement", choices=sorted(ICE_STRATEGIES), default="flank_caps",
                    help="initial ice configuration (default: flank_caps)")
    ap.add_argument("--seed", type=int, default=None,
                    help="override the ice-placement RNG seed (walls stay fixed)")
    ap.add_argument("--throat-gap-frac", type=float, default=None,
                    help="target minimum throat gap as a fraction of Ly "
                         "(default: whatever the amplified walls give, asserted "
                         f">= {MIN_CHANNEL_GAP_FRAC})")
    ap.add_argument("--tag", default=None,
                    help="suffix for the output basenames, e.g. --tag throat -> "
                         "2D_regolith_pore_throat.opts / regolith_pore_throat.dat")
    ap.add_argument("--out", default=None, help="override mesh .dat path")
    ap.add_argument("--opts", default=None, help="override .opts path")
    args = ap.parse_args()

    tag = f"_{args.tag}" if args.tag else ""
    # The mesh depends only on the walls (fixed seeds/GSD), NOT on the ice
    # placement, so every strategy shares ONE mesh .dat; only the .opts (ice
    # grains) is per-strategy.
    dat = Path(args.out) if args.out else ROOT / "inputs/geometry/regolith_pore.dat"
    opts = Path(args.opts) if args.opts else ROOT / f"inputs/geometry/2D_regolith_pore{tag}.opts"
    png = ROOT / f"preprocess/regolith_pore{tag}.png"

    # ---- Walls from the regolith GSD, with a central throat ----
    bot_b = amplify_central(generate_random_bumps(
        np.random.default_rng(SEED_BOT_BUMP), Lx, N_BUMPS, BUMP_R, BUMP_HR, BUMP_OV))
    top_b = amplify_central(generate_random_bumps(
        np.random.default_rng(SEED_TOP_BUMP), Lx, N_BUMPS, BUMP_R, BUMP_HR, BUMP_OV))
    floor_h = lambda x: float(bump_field(bot_b, x)[0])
    ceil_h = lambda x: float(bump_field(top_b, x)[0])

    # ---- Throat report + constraints ----
    xs = np.linspace(0, Lx, 8000)
    gap = Ly - (bump_field(bot_b, xs) + bump_field(top_b, xs))
    min_gap, x_pinch = float(gap.min()), float(xs[int(np.argmin(gap))])
    max_slope = float(np.max(np.abs(np.gradient(bump_field(bot_b, xs), xs))))
    assert min_gap > MIN_CHANNEL_GAP_FRAC * Ly, \
        f"throat too tight: {min_gap/Ly:.0%} Ly < {MIN_CHANNEL_GAP_FRAC:.0%}"
    if args.throat_gap_frac is not None:
        assert abs(min_gap / Ly - args.throat_gap_frac) < 0.08, \
            (f"throat gap {min_gap/Ly:.0%} Ly not near requested "
             f"{args.throat_gap_frac:.0%}; adjust CENTRAL_*_MULT")

    # ---- Ice placement ----
    ice_seed = args.seed if args.seed is not None else (SEED_BOT_CAP * 100 + SEED_TOP_CAP)
    x_lo, x_hi = EDGE_MARGIN, Lx - EDGE_MARGIN
    ice = ICE_STRATEGIES[args.ice_placement](
        np.random.default_rng(ice_seed), bot_b, top_b, floor_h, ceil_h, x_lo, x_hi)
    assert ice, f"ice-placement '{args.ice_placement}' produced no grains"

    # no unintended ice-ice overlaps
    bad = [(i, k) for i in range(len(ice)) for k in range(i + 1, len(ice))
           if overlaps(ice[i], ice[k])]
    assert not bad, f"overlapping ice grains: {bad}"
    Rmin = min(c[2] for c in ice)
    assert EPS / Rmin < 0.06, f"eps/R_smallest = {EPS/Rmin:.1%} (interface under-resolved)"

    Nx = math.ceil(Lx * math.sqrt(2) / EPS)
    Ny = math.ceil(Ly * math.sqrt(2) / EPS)

    print(f"regolith pore channel  Lx={Lx:.4e} Ly={Ly:.4e}  Nx={Nx} Ny={Ny}  "
          f"nodes {(Nx+P)*(Ny+P)/1e6:.2f}M")
    print(f"  walls: {N_BUMPS}+{N_BUMPS} regolith grains (R_med={R_MED:.1e}); "
          f"throat min gap {min_gap/Ly*100:.0f}% Ly at x={x_pinch/Lx:.2f}Lx, "
          f"floor slope {math.degrees(math.atan(max_slope)):.0f} deg")
    print(f"  ice: '{args.ice_placement}' -> {len(ice)} grains, "
          f"R {Rmin:.2e}..{max(c[2] for c in ice):.2e}, eps/R_min {EPS/Rmin*100:.1f}%")

    # ---- Build the mesh via the shared igakit template ----
    def bstr(b):
        return ";".join(f"{cx:.6e},{R:.6e},{h:.6e}" for cx, R, h in b)
    cmd = [str(ROOT / "venv_pf311/bin/python3"),
           str(ROOT / "preprocess/build_geometry_multi_grain.py"),
           "--bumps", bstr(bot_b), "--top-bumps", bstr(top_b),
           "--Lx", f"{Lx}", "--Ly", f"{Ly}", "--Nx", f"{Nx}", "--Ny", f"{Ny}",
           "--P", str(P), "--C", str(C), "--out", str(dat),
           "--plot", str(png.with_name(png.stem + "_mesh.png")),
           "--vtk", "/dev/null"]  # skip the large structured-grid VTK dump
    r = subprocess.run(cmd, capture_output=True, text=True)
    print(r.stdout[-300:])
    if r.returncode != 0:
        print("MESH BUILD FAILED:\n", r.stderr[-2000:])
        sys.exit(1)

    preview_ice(bot_b, top_b, ice, png,
                f"regolith pore channel — ice placement: {args.ice_placement} "
                f"({len(ice)} grains)")

    # ---- Write the .opts ----
    def arr(v):
        return ",".join(f"{z:.6e}" for z in v)
    gx = [c[0] for c in ice]
    gy = [c[1] for c in ice]
    gR = [c[2] for c in ice]
    rel_dat = dat.relative_to(ROOT)
    with open(opts, "w") as f:
        f.write(f"""# =============================================================================
# geometry/{opts.name} — regolith pore channel (implicit regolith).
# Generated by preprocess/build_geometry_regolith_pore.py (deterministic).
# Study: studies/icy_regolith/implicit_pore_domain/ (Effort 1).
#
# Two-sided channel whose top/bottom walls are {N_BUMPS}+{N_BUMPS} lunar-regolith
# grains (R_med={R_MED:.1e} m), with a central throat pinch (min gap
# {min_gap/Ly*100:.0f}% Ly at x={x_pinch/Lx:.2f}Lx). Ice placement strategy:
# '{args.ice_placement}' -> {len(ice)} ice grains in the pore space. The regolith
# is the deformed boundary, NOT a simulated field (that is Effort 2).
#
# eps={EPS:.4e} is the comp_eps.py (Kaempfer&Plapp) value for T={T0_C:g}C,
# alpha_c={ALPHA_C:g} (kinetic-bound limited, so grain-size independent here;
# eps/R_smallest={EPS/Rmin*100:.1f}% < 5%). Reproduces the validated
# 2D_ripening_two_sided reference exactly. PAIR ONLY with a T={T0_C:g}C
# experiment (e.g. 30day_T-20_h1.00_arrh, which sets the matching
# beta_sub0/d0_sub0). RECOMPUTE eps for any other run temperature.
# =============================================================================
# DOF_GRID: {Nx+P} {Ny+P}
-geom_file {rel_dat}
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
-delt_t 1.0e-4
-eps {EPS:.4e}
-periodic 0
""")
    print(f"\nwrote {opts}\nwrote {dat}")
    print(f"eps={EPS:.4e} is the comp_eps.py value for T={T0_C:g}C "
          f"(eps/R_smallest={EPS/Rmin*100:.1f}%). Valid ONLY at T={T0_C:g}C — "
          f"pair with a {T0_C:g}C experiment; RECOMPUTE eps for other temperatures.")


if __name__ == "__main__":
    main()
