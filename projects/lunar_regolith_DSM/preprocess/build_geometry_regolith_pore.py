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
       wall_divots  — grains adhering in smooth DIVOTS carved into the walls,
                      each seat sized to its grain (see below).
     Add more by registering a function in ICE_STRATEGIES.

A strategy may also reshape the walls: it returns (ice, bot_extra, top_extra)
and the extra bump triples are carved in before the mesh is built. wall_divots
uses this to cut each grain a seat — a bump with NEGATIVE height, which the
solver's SedimentBumpField() already sums correctly, so no C change is needed.
Strategies listed in CARVES_WALLS get their own mesh .dat instead of sharing
the pore-only mesh, since their walls genuinely differ.

Writes inputs/geometry/meshes/regolith_pore.dat and inputs/geometry/2D_regolith_pore.opts
(so scripts/Studio/run_lunar.sh finds them by the usual convention), plus a
preview PNG. Deterministic (seeded). Every geometric constraint is asserted
before anything is written — including the wall-slope budget that keeps the
wall-conforming mesh from becoming over-skewed.

Usage (from the project root):
    python3 preprocess/build_geometry_regolith_pore.py --ice-placement flank_caps
    python3 preprocess/build_geometry_regolith_pore.py --ice-placement throat_bridge \\
        --tag throat --seed 3
    python3 preprocess/build_geometry_regolith_pore.py --ice-placement wall_divots \\
        --tag divots

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
from build_geometry_multi_grain import MAX_SED_GRAINS, generate_random_bumps  # noqa: E402

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
# the validated ripening_2D_L1.01mm_eps0.86um_2sided reference run exactly, and pairs with the
# base_T-20_h1.00_30d_a1.34e-2 experiment's beta_sub0=5.9216e5 / d0_sub0=1.0166e-9.
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
MAX_CHANNEL_GAP_FRAC = 0.60          # ...and no more, or it is not a "throat" any more

# Ice-cap sizing (used by flank_caps / throat_bridge flanks)
R_CAP = (2.6e-5, 4.0e-5)
SEED_BOT_CAP, SEED_TOP_CAP = 1, 4
MIN_GAP = 1.8e-5                     # min edge-to-edge gap between caps
BURY_MARGIN = 0.6e-5               # a cap must clear its local bump by this much
EDGE_MARGIN = 4.0e-5               # keep ice this far from the left/right edges

# --------------------------------------------------------------------------- #
# wall_divots parameters — ice grains seated in smooth depressions carved into
# the wall (see ice_wall_divots()).
# --------------------------------------------------------------------------- #
# A divot is just another C-infinity bump with NEGATIVE height, summed into the
# same -sed_grain_*/-top_grain_* field the solver reads (SedimentBump() in
# src/initial_conditions.c has no sign restriction), so the wall stays
# C-infinity and no solver change is needed.
#
# The binding constraint is MESH QUALITY, not the divot shape in isolation. The
# mesh maps y(u,v) = y_bot(x) + v*(y_top(x)-y_bot(x)), so element skew near a
# wall tracks that wall's slope. The validated stiff-wall geometry peaks at
# 46 deg, and a divot's own slope ADDS to whatever the wall already has there.
# Two things keep the total in budget:
#   1. softer regolith relief (BUMP_HR_SOFT) frees ~6 deg of headroom;
#   2. divots are only seated at wall EXTREMA, where the local wall slope is 0.
# Seating them in troughs alone (the pore_lining rule) is NOT enough: measured
# 55-59 deg, because trough bottoms are not flat enough and the divot flank
# lands on the neighbouring grain's rise.
BUMP_HR_SOFT = (0.22, 0.36)          # softened height/half-width for wall_divots
DIVOT_H_MULT = 2.8                   # central amplification that restores the throat
R_ICE_DIVOT = (1.8e-5, 3.0e-5)       # seated ice grain radius range [m]
DIVOT_DEPTH_F = 0.5                  # divot depth / R_ice
DIVOT_HW_F = 1.75                    # divot half-width / R_ice (seat ~1.75 grain
                                     # diameters wide; wider seats crowd each
                                     # other out before 6 fit across the wall)
DIVOT_SEP = 6.0e-6                   # min gap between adjacent divot supports [m]
SLOPE_BUDGET = 1.10                  # max |dy/dx| of the carved wall (~47.7 deg)
DIVOT_N_TARGET = 6                   # per wall; fewer is OK if the budget binds
DIVOT_OPP_CLEAR_F = 1.0              # grain must clear the opposing wall by this x R


# --------------------------------------------------------------------------- #
# Shared geometry helpers (match SedimentBumpField / SedimentBump in
# src/initial_conditions.c and _bump_field in build_geometry_multi_grain.py)
# --------------------------------------------------------------------------- #
def amplify_central(bumps, h_mult=CENTRAL_H_MULT):
    """Scale central bump height (x h_mult) and half-width
    (x CENTRAL_R_MULT), tapering to 1x at the flanks, forming the throat."""
    sig = CENTRAL_SIG_FRAC * Lx
    out = []
    for cx, R, h in bumps:
        env = math.exp(-((cx - 0.5 * Lx) / sig) ** 2)
        out.append((cx, R * (1 + (CENTRAL_R_MULT - 1) * env),
                    h * (1 + (h_mult - 1) * env)))
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


def bump_field_slope(bumps, x0=0.0, x1=None, n=40000):
    """max |d/dx| of the summed bump field over [x0, x1] — the mesh-quality
    proxy. Element skew near a wall tracks that wall's slope (the mesh
    interpolates linearly in v between floor and ceiling), so this is what
    bounds how aggressively a divot may be carved.

    The window matters when SITING divots: a whole-domain max would let the
    first divot that reaches the budget veto every later one anywhere else on
    the wall, no matter how flat that spot is. Divot supports never overlap, so
    checking each one over its own support is equivalent — and main() still
    asserts the whole-domain max on the finished wall."""
    x1 = Lx if x1 is None else x1
    span = max(x1 - x0, 1e-12)
    xq = np.linspace(x0, x1, max(int(n * span / Lx), 64))
    return float(np.max(np.abs(np.gradient(bump_field(bumps, xq), xq))))


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
# Ice-placement strategies. Each returns (ice, bot_extra, top_extra):
#   ice       list of [cx, cy, R] ice grains (cy=0 -> sits on the floor wall,
#             cy=Ly -> ceiling wall, 0<cy<Ly -> interior)
#   bot_extra extra (cx, R, h) wall bumps the strategy needs CARVED INTO the
#   top_extra bottom/top wall. Empty for every strategy that only places ice
#             into the pore; wall_divots uses them (with h<0) to cut the seats
#             its grains sit in. main() appends them to the wall bump lists
#             BEFORE the mesh is built, so mesh and .opts always agree.
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
    return caps, [], []


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
    caps, _, _ = ice_flank_caps(rng, bot_b, top_b, floor_h, ceil_h, x_lo, x_hi)
    # drop any flank cap that would overlap the bridge
    caps = [c for c in caps if not overlaps(c, bridge)]
    return [bridge] + caps, [], []


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
    return ice, [], []


def _wall_flat_sites(bumps, x_lo, x_hi, n=4000):
    """x-positions of ALL local extrema of the wall protrusion — trough bottoms
    AND grain crests. These are the only places a divot can be carved without
    its own slope stacking on top of an existing wall slope (see the
    SLOPE_BUDGET note in the parameter block).

    Crests are legitimate seats physically as well as numerically: a dimple on
    top of a regolith grain is as plausible an adhesion site as one in the
    valley between two grains, and near the throat — where the wall is a single
    wide amplified hump — crests are the ONLY candidate sites."""
    x = np.linspace(x_lo, x_hi, n)
    d = np.gradient(bump_field(bumps, x), x)
    return [float(x[i]) for i in range(1, n - 1) if d[i - 1] * d[i + 1] <= 0.0]


def _carve_divots(wall, x_lo, x_hi):
    """Carve DIVOT_N_TARGET seats into one wall, spread across its whole length.

    Sites are anchored to evenly spaced targets rather than taken greedily
    left-to-right: a greedy walk takes the widest seat it can at every chance
    and runs out of room before the right-hand end, leaving the COLD quarter of
    the domain bare — precisely where a warm->cold gradient is expected to pile
    ice up. Each target claims the nearest usable wall extremum instead, so
    both ends of the channel are seeded.

    At the chosen site take the LARGEST ice radius that still fits: the divot
    support must stay inside [x_lo, x_hi] (so no divot touches the Dirichlet-T
    end walls), must not overlap an already-placed divot's support, and the
    carved wall must stay within SLOPE_BUDGET. Testing the budget on the
    *carved* field, not on the divot alone, is what keeps the sum bounded.

    Returns (divots, seats): divots as (cx, R_div, -depth) bump triples,
    seats as [cx, R_ice] for the grains that will sit in them."""
    sites = _wall_flat_sites(wall, x_lo, x_hi)
    step = (x_hi - x_lo) / DIVOT_N_TARGET
    targets = [x_lo + (k + 0.5) * step for k in range(DIVOT_N_TARGET)]

    divots, seats = [], []          # kept sorted by x via `placed`
    placed = []                     # (cx, half-width) of accepted divots
    for x_t in targets:
        for x_d in sorted(sites, key=lambda s: abs(s - x_t)):
            hit = False
            for R_ice in np.linspace(R_ICE_DIVOT[1], R_ICE_DIVOT[0], 13):
                hw = DIVOT_HW_F * float(R_ice)
                if x_d - hw < x_lo or x_d + hw > x_hi:
                    continue
                if any(abs(x_d - px) < hw + phw + DIVOT_SEP for px, phw in placed):
                    continue
                trial = divots + [(x_d, hw, -DIVOT_DEPTH_F * float(R_ice))]
                if bump_field_slope(wall + trial,
                                    x_d - hw, x_d + hw) <= SLOPE_BUDGET:
                    divots = trial
                    placed.append((x_d, hw))
                    seats.append([x_d, float(R_ice)])
                    hit = True
                    break
            if hit:
                break
    order = sorted(range(len(seats)), key=lambda i: seats[i][0])
    return [divots[i] for i in order], [seats[i] for i in order]


def ice_wall_divots(rng, bot_b, top_b, floor_h, ceil_h, x_lo, x_hi):
    """Ice grains adhering to the regolith walls, each seated in a smooth divot
    of its own size — the configuration for the dT/dx = -50 K/m aggregation
    study (studies/icy_regolith/implicit_pore_domain/).

    Deliberately NOT a big-ice-block-vs-big-grain setup: every grain is a
    comparable, wall-adhering deposit, so where ice ends up is an outcome of
    the gradient and the pore geometry rather than of the initial condition.

    Each grain's centre goes at the BOTTOM of its own divot, so the wall clips
    the disc on both sides and the exposed ice is a little less than a
    semicircle nestled in the seat. Floor and ceiling walls use different
    seeds, so the seats are naturally staggered in x — grains face open
    channel, not each other."""
    bot_dv, bot_seats = _carve_divots(bot_b, x_lo, x_hi)
    top_dv, top_seats = _carve_divots(top_b, x_lo, x_hi)

    bot_carved, top_carved = bot_b + bot_dv, top_b + top_dv
    ice = [[x, float(bump_field(bot_carved, x)[0]), R] for x, R in bot_seats]
    ice += [[x, Ly - float(bump_field(top_carved, x)[0]), R] for x, R in top_seats]
    return ice, bot_dv, top_dv


ICE_STRATEGIES = {
    "flank_caps": ice_flank_caps,
    "throat_bridge": ice_throat_bridge,
    "pore_lining": ice_pore_lining,
    "wall_divots": ice_wall_divots,
}

# Per-strategy wall parameters, resolved BEFORE the walls are generated.
#
# The first three strategies differ only in where ice is placed, so they share
# one mesh (inputs/geometry/meshes/regolith_pore.dat) and MUST keep the original stiff
# wall values — re-cutting those walls would silently invalidate every run
# already done against them. wall_divots spends slope budget on the divots, so
# it softens the regolith relief to buy that budget back and raises the central
# amplification to keep a real throat despite the softer grains.
_STIFF_WALL = (BUMP_HR, CENTRAL_H_MULT)
WALL_PARAMS = {
    "flank_caps": _STIFF_WALL,
    "throat_bridge": _STIFF_WALL,
    "pore_lining": _STIFF_WALL,
    "wall_divots": (BUMP_HR_SOFT, DIVOT_H_MULT),
}

# Strategies that reshape the walls need their own mesh, not the shared one.
CARVES_WALLS = {"wall_divots"}


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
    carves = args.ice_placement in CARVES_WALLS
    # The mesh depends only on the WALLS, not on where ice is placed, so the
    # strategies that only fill the pore all share ONE mesh .dat and differ
    # only in their .opts. A strategy that carves the walls (wall_divots) has a
    # different mesh by construction and gets its own .dat — overwriting the
    # shared one would silently change the other three geometries.
    if args.out:
        dat = Path(args.out)
    elif carves:
        dat = ROOT / f"inputs/geometry/regolith_pore{tag or '_' + args.ice_placement}.dat"
    else:
        dat = ROOT / "inputs/geometry/meshes/regolith_pore.dat"
    opts = Path(args.opts) if args.opts else ROOT / f"inputs/geometry/regolithpore_2D_L814um_eps0.86um{tag}.opts"
    png = ROOT / f"preprocess/regolith_pore{tag}.png"

    # ---- Walls from the regolith GSD, with a central throat ----
    bump_hr, h_mult = WALL_PARAMS[args.ice_placement]
    bot_b = amplify_central(generate_random_bumps(
        np.random.default_rng(SEED_BOT_BUMP), Lx, N_BUMPS, BUMP_R, bump_hr, BUMP_OV), h_mult)
    top_b = amplify_central(generate_random_bumps(
        np.random.default_rng(SEED_TOP_BUMP), Lx, N_BUMPS, BUMP_R, bump_hr, BUMP_OV), h_mult)
    floor_h = lambda x: float(bump_field(bot_b, x)[0])
    ceil_h = lambda x: float(bump_field(top_b, x)[0])

    # ---- Ice placement (may carve divots into the walls) ----
    ice_seed = args.seed if args.seed is not None else (SEED_BOT_CAP * 100 + SEED_TOP_CAP)
    x_lo, x_hi = EDGE_MARGIN, Lx - EDGE_MARGIN
    ice, bot_dv, top_dv = ICE_STRATEGIES[args.ice_placement](
        np.random.default_rng(ice_seed), bot_b, top_b, floor_h, ceil_h, x_lo, x_hi)
    assert ice, f"ice-placement '{args.ice_placement}' produced no grains"
    bot_b, top_b = bot_b + bot_dv, top_b + top_dv   # walls the mesh will be cut from

    # A wall's bump list goes straight into -sed_grain_*/-top_grain_*, which
    # PETSc's option parser SILENTLY TRUNCATES past MAX_SED_GRAINS.
    for side, b in (("bottom", bot_b), ("top", top_b)):
        assert len(b) <= MAX_SED_GRAINS, \
            (f"{side} wall has {len(b)} bumps > MAX_SED_GRAINS={MAX_SED_GRAINS} "
             "(include/NASA_types.h) — PETSc would truncate them without error")

    # ---- Throat + mesh-quality constraints (on the FINAL, carved walls) ----
    xs = np.linspace(0, Lx, 8000)
    floor_y, ceil_drop = bump_field(bot_b, xs), bump_field(top_b, xs)
    gap = Ly - (floor_y + ceil_drop)
    min_gap, x_pinch = float(gap.min()), float(xs[int(np.argmin(gap))])
    bot_slope, top_slope = bump_field_slope(bot_b), bump_field_slope(top_b)
    assert min_gap > MIN_CHANNEL_GAP_FRAC * Ly, \
        f"throat too tight: {min_gap/Ly:.0%} Ly < {MIN_CHANNEL_GAP_FRAC:.0%}"
    assert min_gap < MAX_CHANNEL_GAP_FRAC * Ly, \
        (f"throat too open: {min_gap/Ly:.0%} Ly > {MAX_CHANNEL_GAP_FRAC:.0%} — "
         "not a pore throat any more; raise the central height multiplier")
    # Element skew near a wall tracks that wall's slope (the mesh interpolates
    # linearly in v between floor and ceiling), so this is the mesh-quality gate.
    for side, s in (("floor", bot_slope), ("ceiling", top_slope)):
        assert s <= SLOPE_BUDGET, \
            (f"{side} max slope {math.degrees(math.atan(s)):.0f} deg exceeds the "
             f"{math.degrees(math.atan(SLOPE_BUDGET)):.0f} deg budget — mesh would be "
             "over-skewed; soften BUMP_HR or shrink DIVOT_DEPTH_F")
    assert float((Ly - ceil_drop - floor_y).min()) > 0.0, "floor and ceiling cross"
    if args.throat_gap_frac is not None:
        assert abs(min_gap / Ly - args.throat_gap_frac) < 0.08, \
            (f"throat gap {min_gap/Ly:.0%} Ly not near requested "
             f"{args.throat_gap_frac:.0%}; adjust CENTRAL_*_MULT")

    # no unintended ice-ice overlaps
    bad = [(i, k) for i in range(len(ice)) for k in range(i + 1, len(ice))
           if overlaps(ice[i], ice[k])]
    assert not bad, f"overlapping ice grains: {bad}"
    Rmin = min(c[2] for c in ice)
    assert EPS / Rmin < 0.06, f"eps/R_smallest = {EPS/Rmin:.1%} (interface under-resolved)"

    # A WALL-SEATED grain must leave the pore open across from it. Grains a
    # strategy deliberately suspends in mid-channel (throat_bridge) are exempt:
    # spanning the gap is the whole point there, so the strategy owns that call.
    for cx, cy, R in ice:
        y_bot = float(bump_field(bot_b, cx)[0])
        y_top = Ly - float(bump_field(top_b, cx)[0])
        on_bot, on_top = abs(cy - y_bot) <= R, abs(cy - y_top) <= R
        if not (on_bot or on_top):
            continue                       # interior grain, not wall-seated
        clear = (y_top - (cy + R)) if on_bot else ((cy - R) - y_bot)
        assert clear >= DIVOT_OPP_CLEAR_F * R, \
            (f"wall-seated ice grain at x={cx:.3e} clears the opposing wall by only "
             f"{clear/R:.2f} R (< {DIVOT_OPP_CLEAR_F} R) — it would bridge the channel")

    Nx = math.ceil(Lx * math.sqrt(2) / EPS)
    Ny = math.ceil(Ly * math.sqrt(2) / EPS)

    print(f"regolith pore channel  Lx={Lx:.4e} Ly={Ly:.4e}  Nx={Nx} Ny={Ny}  "
          f"nodes {(Nx+P)*(Ny+P)/1e6:.2f}M")
    print(f"  walls: {N_BUMPS}+{N_BUMPS} regolith grains (R_med={R_MED:.1e}, "
          f"h/R {bump_hr[0]:.2f}-{bump_hr[1]:.2f}); "
          f"throat min gap {min_gap/Ly*100:.0f}% Ly at x={x_pinch/Lx:.2f}Lx, "
          f"slope floor/ceiling {math.degrees(math.atan(bot_slope)):.0f}/"
          f"{math.degrees(math.atan(top_slope)):.0f} deg "
          f"(budget {math.degrees(math.atan(SLOPE_BUDGET)):.0f})")
    if bot_dv or top_dv:
        print(f"  divots: {len(bot_dv)} floor + {len(top_dv)} ceiling, "
              f"depth {DIVOT_DEPTH_F:g}*R_ice, half-width {DIVOT_HW_F:g}*R_ice; "
              f"bumps/wall {len(bot_b)}/{len(top_b)} of {MAX_SED_GRAINS}")
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
    divot_note = ""
    if bot_dv or top_dv:
        divot_note = f"""
#
# WALL DIVOTS. {len(bot_dv)} floor + {len(top_dv)} ceiling smooth depressions are carved
# into the walls, one per ice grain, each depth {DIVOT_DEPTH_F:g}*R_ice and half-width
# {DIVOT_HW_F:g}*R_ice, so every grain sits nestled in a seat of its own size rather
# than perched on a bare wall. A divot is just another C-infinity bump with
# NEGATIVE height in the -sed_grain_*/-top_grain_* lists below, so the wall
# stays C-infinity and the solver needs no special handling.
#
# Divots cost mesh slope budget, so this geometry uses SOFTER regolith relief
# (h/R {bump_hr[0]:.2f}-{bump_hr[1]:.2f} vs {BUMP_HR[0]:.2f}-{BUMP_HR[1]:.2f} for the other strategies) and a higher
# central amplification to keep a real throat. Carved slopes are {math.degrees(math.atan(bot_slope)):.0f}/{math.degrees(math.atan(top_slope)):.0f} deg
# (floor/ceiling), within the {math.degrees(math.atan(SLOPE_BUDGET)):.0f} deg budget.
#
# Because the walls differ, this geometry has its OWN mesh — it does not share
# inputs/geometry/meshes/regolith_pore.dat with flank_caps/throat_bridge/pore_lining."""
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
# is the deformed boundary, NOT a simulated field (that is Effort 2).{divot_note}
#
# eps={EPS:.4e} is the comp_eps.py (Kaempfer&Plapp) value for T={T0_C:g}C,
# alpha_c={ALPHA_C:g} (kinetic-bound limited, so grain-size independent here;
# eps/R_smallest={EPS/Rmin*100:.1f}% < 5%). Reproduces the validated
# ripening_2D_L1.01mm_eps0.86um_2sided reference exactly. PAIR ONLY with a T={T0_C:g}C
# experiment (e.g. base_T-20_h1.00_30d_a1.34e-2, which sets the matching
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
-eps_valid_temp {T0_C:g}   # C: temperature eps/mesh were sized for; solver ABORTS if -temp differs (override: -eps_temp_override 1)
-periodic 0
""")
    print(f"\nwrote {opts}\nwrote {dat}")
    print(f"eps={EPS:.4e} is the comp_eps.py value for T={T0_C:g}C "
          f"(eps/R_smallest={EPS/Rmin*100:.1f}%). Valid ONLY at T={T0_C:g}C — "
          f"pair with a {T0_C:g}C experiment; RECOMPUTE eps for other temperatures.")


if __name__ == "__main__":
    main()
