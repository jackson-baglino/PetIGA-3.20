#!/usr/bin/env python3
"""generate_packing_gravity.py — 2D grain packings that look gravity-deposited,
at a prescribed mean grain size, size distribution and porosity.

    python3 preprocess/generate_packing_gravity.py \\
        --Lx 2.0e-3 --porosity 0.325 --mean-r 100e-6 --sigma-ln 0.5 \\
        --seed 4 --out inputs/packings/grav_phi0.325_seed4

WHAT "GRAVITY" MEANS HERE, AND WHAT IT DOES NOT
===============================================
There is no DEM, no contact law, no time integration. Grains are placed one at
a time by a purely geometric drop-and-roll rule: released above the bed, lowered
until they touch something, then allowed to roll downhill around that contact
until they find a second support or run out of rolling budget.

That is enough to reproduce the structural signature that matters -- every
grain rests on something below it, so load paths run to the floor and the solid
matrix is connected BY CONSTRUCTION rather than by luck. It is not a settling
simulation and should not be reported as one.

THE POROSITY KNOB IS PHYSICAL
-----------------------------
`--roll-tol` is how far (in radians) a grain may roll around its first contact
before it locks:

    roll-tol -> 0        lock on first touch. Ballistic deposition, the loosest
                         packing the rule can make (porosity ~ 0.45-0.50).
    roll-tol -> pi/2     roll to a stable two-point seat every time. Dense
                         settling (porosity ~ 0.16-0.20).

The range brackets the 0.25-0.40 the snow-thermal study uses from both sides,
so porosity is dialled by bisecting a real deposition parameter rather than by
distorting the geometry afterwards. Radii are never rescaled to hit a target --
that would break contact tangency, which is the whole point of these packings.

MATRIX FIRST, THEN FILLERS
--------------------------
Deposition is bisected to land somewhat LOOSER than the target
(`--fill-margin`), and the remaining porosity is removed by inserting grains
into the largest voids, biggest hole first.

This single pass serves two goals at once. It brings porosity down to target,
and because it always attacks the largest inscribed empty circle, it drives the
pore heterogeneity down monotonically at the same time. Filler grains are not
required to rest on anything -- some of them float.

FLOATING GRAINS ARE INTENTIONAL. A 2D section through a 3D packing shows
grains with no visible support all the time; more to the point, ice sinters, so
a grain that is connected out of plane is mechanically real even where this
slice cannot show it. The matrix is built connected first precisely so that the
floaters are a garnish on a percolating skeleton rather than the structure
itself. `n_matrix_grains` / `n_filler_grains` in metadata.json make the split
visible, and the preview draws fillers in a different colour.

HOMOGENEITY IS A GATE, NOT A HOPE
---------------------------------
For k_eff by homogenization the domain has to be a usable REV, so a packing
that is well connected on one side and sparse on the other is worse than
useless -- it produces a number that looks fine and means nothing. The
generator therefore refuses to emit a packing that fails

    max_void_radius_per_mean_r <= --max-void-ratio
    local_density_cv           <= --max-density-cv
    solid percolates in x and y

and re-seeds up to `--max-tries` instead. Rejecting is deliberate: the retired
select_packing_seeds.py ranked candidates instead, and see packing_lib's
accept_reasons() for why ranking on PORE connectivity actively selects for the
clumped configurations this gate exists to exclude.

Those limits are not arbitrary numbers, they are calibrated against the 25
accepted jammed packings in inputs/scratch/packings/:

  --max-void-ratio   defaults to a porosity-dependent envelope. BOTH this and
  --max-density-cv   the CV grow with porosity in perfectly good packings
                     (void 1.04 -> 1.58 and CV 0.19 -> 0.42 going from phi 0.25
                     to 0.40 in that set), so a FLAT gate is silently a
                     porosity gate. Both defaults started flat and both were
                     wrong: 1.30 and 0.40 reject the reference standard above
                     phi ~ 0.33, and phi 0.400 failed in all three periodicity
                     modes while hitting its porosity target to four decimals.
  --max-asymmetry    half-domain density difference, in x and in y. The
                     reference packings run to 28% with a median of 11%.

The asymmetry gate exists because local_density_cv does NOT see the failure it
is named for. On an 8x8 grid the CV measures fine-scale clumping, so a packing
can sit well inside the CV gate while one half of the domain carries 20% more
solid than the other -- exactly "well connected on one side, sparse on the
other".

Represent macro-pores by running a HIGHER-porosity REV, not by letting one REV
go patchy.

THE DOMAIN IS A WINDOW CUT FROM A LARGER BED
--------------------------------------------
Deposition runs in a strip taller than the domain; the output window is cropped
from its interior, discarding the floor contact layer and the loose top
surface, which are both boundary artifacts. Grains are kept when ANY part of
them falls inside -- centres span [-R, L+R] -- so the grains do not "see" the
boundary. See packing_lib.write_outputs for why dropping edge-straddling grains
biases k_eff low.

PERIODICITY
-----------
    --periodic x     (default) the deposition strip wraps sideways. Gravity
                     acts in y, so the x-wrap costs nothing and removes the
                     side-wall ordering a closed box would impose.
    --periodic xy    additionally wrap in y, healing the seam by relaxing
                     overlaps under the minimum image. Only grains near the
                     seam move; the reported max displacement says how much.
    --periodic y     wraps along gravity but not across it. Provided for BC
                     comparisons; x or xy is the more natural REV.
    --periodic none  plain window, no edge images.

The flag controls the OUTPUT cell -- which edge images are written, how the
raster and the metrics wrap, and how the window is cropped. Deposition itself
always wraps in x, because the bed being modelled is wide and the domain is a
window cut from it; a strip with free sides is a generator artifact rather than
a microstructure anyone wants (see the note in deposit_at_porosity).

Because the deposition is geometric there is no compaction under load, so
density has no systematic drift with depth and the y-wrap is defensible. That
is checked rather than assumed: across seeds the bottom-vs-top density
difference changes SIGN and does not shrink when the crop margin is widened
from 3 to 10 mean radii, which is what configurational variance looks like and
not what a depth gradient looks like. The magnitude is not small, so it is
gated (--max-asymmetry) rather than waved through.

OUTPUT is grains.dat + metadata.json + preview.png, the same contract the
solver already reads (-ic_type multi_grains_file -grains_file ...).

Output is grouped by periodicity, with a periodic_<mode>/ level inserted above
the leaf name:

    --out inputs/packings/phi0.325_seed4 --periodic x
      -> inputs/packings/periodic_x/phi0.325_seed4/

The same parameters under different wrapping are different microstructures and
belong to different runs, but they naturally get the same descriptive name, so
without the split the second would silently overwrite the first. Pass
--no-periodic-subdir to write straight to --out.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import packing_lib as pl                                          # noqa: E402


# =========================================================================
# Drop-and-roll deposition
# =========================================================================

def _rest_height(x, r, cen, rad, Lx, px):
    """Highest resting position for a grain of radius r released at x.

    Returns (y, support_index); support_index is -1 for the floor. A falling
    grain stops at the FIRST thing it meets, which is the maximum over all
    candidate contact heights.
    """
    y, k = r, -1                                   # the floor
    if len(rad):
        dx = cen[:, 0] - x
        if px:
            dx -= Lx * np.round(dx / Lx)
        reach = rad + r
        m = np.abs(dx) < reach
        if np.any(m):
            cand = cen[m, 1] + np.sqrt(reach[m] ** 2 - dx[m] ** 2)
            j = int(np.argmax(cand))
            if cand[j] > y:
                y, k = float(cand[j]), int(np.nonzero(m)[0][j])
    return y, k


def _clearance(x, y, r, cen, rad, Lx, px, skip=-1):
    """Smallest surface gap between this grain and the settled ones.

    Negative means overlap. `skip` excludes the current support.
    """
    if not len(rad):
        return np.inf
    d = np.stack([cen[:, 0] - x, cen[:, 1] - y], axis=1)
    if px:
        d[:, 0] -= Lx * np.round(d[:, 0] / Lx)
    gap = np.hypot(d[:, 0], d[:, 1]) - (rad + r)
    if 0 <= skip < len(gap):
        gap[skip] = np.inf
    return float(gap.min())


def _settle(x, r, cen, rad, Lx, px, roll_tol, n_steps=24, max_hops=4):
    """Place one grain by dropping, then rolling downhill around its contact.

    Rolling stops at whichever comes first: a second contact (a stable seat),
    the rolling budget `roll_tol`, the floor, or rolling past the support's
    equator -- in which case the grain has fallen off and is re-dropped from
    where it left, up to `max_hops` times.

    Returns (x, y).
    """
    for _ in range(max_hops):
        y, k = _rest_height(x, r, cen, rad, Lx, px)
        if k < 0 or roll_tol <= 0.0:
            return x, y                            # on the floor, or no rolling

        cx, cy = cen[k]
        dx = x - cx
        if px:
            dx -= Lx * np.round(dx / Lx)
        R = rad[k] + r
        theta = math.atan2(dx, y - cy)             # 0 = directly above support
        if abs(theta) < 1e-9:                      # balanced on the apex
            theta = 1e-3 * (1.0 if (hash((round(x, 12), round(r, 15))) & 1) else -1.0)
        sgn = 1.0 if theta > 0 else -1.0

        budget = min(roll_tol, math.pi / 2 - abs(theta))
        if budget <= 0.0:
            return x, y
        step = budget / n_steps

        rolled = 0.0
        bx, by = x, y
        fell_off = False
        while rolled < budget - 1e-12:
            t = theta + sgn * (rolled + step)
            nx, ny = cx + R * math.sin(t), cy + R * math.cos(t)
            if ny < r:                             # would go through the floor
                break
            if _clearance(nx, ny, r, cen, rad, Lx, px, skip=k) < 0.0:
                break                              # second contact: stable seat
            bx, by, rolled = nx, ny, rolled + step
            if abs(t) >= math.pi / 2 - 1e-9:
                fell_off = True
                break

        if not fell_off:
            return bx, max(by, r)
        x = bx                                     # rolled off; drop again
    return x, max(by, r)


def surface_profile(cen, rad, Lx, px, nbins=256):
    """Height of the bed surface sampled on a coarse x-grid.

    The top of a deposited bed is rough on the scale of a grain, so its
    "height" is a profile, not a number.
    """
    if not len(rad):
        return np.zeros(nbins)
    xs = (np.arange(nbins) + 0.5) * (Lx / nbins)
    h = np.zeros(nbins)
    for (cx, cy), r in zip(cen, rad):
        dx = xs - cx
        if px:
            dx -= Lx * np.round(dx / Lx)
        m = np.abs(dx) < r
        if np.any(m):
            h[m] = np.maximum(h[m], cy + np.sqrt(r * r - dx[m] ** 2))
    return h


def deposit(rng, Lx, height, mean_r, sigma_ln, clip_frac, roll_tol, px,
            max_grains=200000, check_every=20):
    """Fill an Lx-wide strip until the bed covers `height` EVERYWHERE.

    The stop test is the MINIMUM of the surface profile, not the highest grain.
    Stopping on the highest grain looks equivalent and is not: the surface is
    rough on the scale of a grain, so a single tall grain ends deposition while
    the mean surface is still several radii lower. The crop window then reaches
    up into the loose, half-filled surface layer, and every subsequent step
    inherits the damage -- the filler pass finds all its "largest voids" along
    the top edge and piles the fillers there in a band, which is visible in the
    preview and is not a packing anyone wants.

    Radii are drawn one at a time from the same log-normal that sample_radii
    uses, so the realized size distribution is the requested one -- no
    grow-to-touch step overwrites it.

    Returns (centres, radii).
    """
    mu = math.log(mean_r) - 0.5 * sigma_ln ** 2
    lo, hi = (1.0 - clip_frac) * mean_r, (1.0 + clip_frac) * mean_r
    if lo <= 0.0:
        lo = 1e-3 * mean_r

    cen = np.empty((0, 2)); rad = np.empty(0)
    cx, cr = [], []
    while len(cr) < max_grains:
        r = float(np.clip(rng.lognormal(mu, sigma_ln), lo, hi))
        x = float(rng.uniform(0.0, Lx))
        x, y = _settle(x, r, cen, rad, Lx, px, roll_tol)
        if px:
            x = float(np.mod(x, Lx))
        cx.append([x, y]); cr.append(r)
        cen = np.asarray(cx); rad = np.asarray(cr)

        # The profile is O(N) to build, so amortise it.
        if len(cr) % check_every == 0:
            if surface_profile(cen, rad, Lx, px).min() >= height:
                break
    return cen, rad


# =========================================================================
# Window + fillers
# =========================================================================

def crop_window(cen, rad, Lx, Ly, y0, py=False):
    """Grains from the window [0,Lx] x [y0, y0+Ly], shifted so y0 -> 0.

    OPEN y (py=False): keep every grain with ANY part inside. These packings
    represent a window cut from a much larger bed, and the grains must not see
    the edge -- a grain centred just outside still puts solid, and a conduction
    path, inside the domain.

    PERIODIC y (py=True): keep only grains whose CENTRE is in [0, Ly). The
    caller is about to identify y=0 with y=Ly, and under that identification a
    grain poking out of the top IS the grain re-entering at the bottom -- it is
    its own image. Keeping the overlap set instead would fold the top and
    bottom straddlers onto each other and count that boundary layer twice,
    which showed up as porosity landing 0.03-0.09 BELOW target no matter how
    the rolling budget was bisected.
    """
    y = cen[:, 1] - y0
    if py:
        keep = (y >= 0.0) & (y < Ly)
    else:
        keep = (y + rad > 0.0) & (y - rad < Ly)
    return np.stack([cen[keep, 0], y[keep]], axis=1), rad[keep]


def fill_voids(rng, cen, rad, Lx, Ly, px, py, target_porosity, mean_r,
               sigma_ln, clip_frac, raster, context=None, max_fill=4000):
    """Insert grains into the largest voids until porosity reaches the target.

    Biggest hole first, so porosity falls to target while pore heterogeneity
    falls with it -- the two goals are served by the same ordering.

    EACH GRAIN IS CAPPED BY THE REMAINING AREA DEFICIT, not just by the hole it
    is filling. Sizing a filler to its void alone overshoots badly: a single
    100 um grain is 6% of a 1 mm domain, so two of them took a test packing
    from phi 0.356 straight past a 0.325 target to 0.270. Porosity is what the
    caller asked for, so it is the hard constraint; the last grain is shrunk to
    land on it and the loop stops rather than overshooting.

    Returns (centres, radii, n_inserted).
    """
    mu = math.log(mean_r) - 0.5 * sigma_ln ** 2
    lo, hi = (1.0 - clip_frac) * mean_r, (1.0 + clip_frac) * mean_r
    if lo <= 0.0:
        lo = 1e-3 * mean_r
    r_floor = 0.15 * mean_r          # below this a "grain" is mesh noise
    domain = Lx * Ly

    added = 0
    for _ in range(max_fill):
        solid = pl.rasterize(cen, rad, Lx, Ly, raster, raster, px, py)
        porosity = 1.0 - float(solid.mean())
        if porosity <= target_porosity:
            break

        # Area still to be removed, and the largest radius that fits in it.
        deficit = (porosity - target_porosity) * domain
        r_deficit = math.sqrt(deficit / math.pi)

        # `context` is the surrounding bed the window crop discarded. Without
        # it an open edge reads as empty beyond the boundary and every filler
        # ends up in a band along that edge.
        # Same convention as grade(): the parent bed, x always wrapped.
        dist, dx, dy = pl.void_map(cen, rad, Lx, Ly, raster, raster, True, py,
                                   context=None if py else context)
        big = float(dist.max())
        if big < r_floor or r_deficit < r_floor:
            break                     # no room, or nothing meaningful left

        j, i = np.unravel_index(int(np.argmax(dist)), dist.shape)
        x, y = (i + 0.5) * dx, (j + 0.5) * dy

        r = float(np.clip(rng.lognormal(mu, sigma_ln), lo, hi))
        r = min(r, big, r_deficit)
        if r < r_floor:
            break

        cen = np.vstack([cen, [x, y]])
        rad = np.append(rad, r)
        added += 1
    return cen, rad, added


# =========================================================================
# Grading
# =========================================================================

def auto_void_ratio(target_porosity):
    """Default max_void_radius_per_mean_r gate, scaled with target porosity.

    Same trap as auto_density_cv: the largest inscribed pore circle grows with
    the porosity, so a flat gate silently becomes a porosity gate. Worst case
    per group over the 25 accepted jammed packings in inputs/scratch/packings/:

        phi 0.250 -> 0.97    phi 0.325 -> 1.07    phi 0.400 -> 1.66
        phi 0.287 -> 1.06    phi 0.362 -> 1.45

    which fits void_max ~ 4.7*phi - 0.29. The default is that envelope plus a
    little headroom. A flat 1.30 -- the value this started with -- rejects the
    reference standard for every phi >= 0.33, and did exactly that: phi 0.400
    failed in all three periodicity modes while hitting its porosity target to
    four decimals.

    These are the reference packings RE-MEASURED with packing_lib.periodic_edt,
    not the numbers in their metadata.json. Those were computed with a plain
    distance transform that cannot see solid across a periodic face, and are
    off by 0.86x to 1.05x -- in both directions, so it does not cancel. Gating
    a corrected measurement against an uncorrected calibration would be
    comparing two different quantities.
    """
    return max(0.9, 4.7 * target_porosity - 0.19)


def auto_density_cv(target_porosity):
    """Default local_density_cv gate, scaled with the target porosity.

    A FIXED CV gate is wrong, because the CV of a perfectly good packing rises
    with its porosity -- more void means more variance in the coarse-cell solid
    fraction. Measured on the 25 accepted jammed packings in
    inputs/scratch/packings/, the worst CV in each porosity group runs

        phi 0.250 -> 0.19    phi 0.325 -> 0.32    phi 0.400 -> 0.42
        phi 0.287 -> 0.24    phi 0.362 -> 0.33

    which is close to linear: cv_worst ~ 1.5*phi - 0.19. The default here is
    that envelope plus ~0.05 of headroom, so the gate means "no worse than the
    packings already in use at this porosity" rather than an arbitrary number.
    A flat 0.40 gate, for comparison, would have rejected phi0.400_seed6
    (cv 0.416), which is a packing the study accepted.
    """
    return max(0.15, 1.5 * target_porosity - 0.14)


def coarse_grid_for(Lx, mean_r):
    """Coarse-cell count for local_density_cv, scaled to the grain size.

    A FIXED 8x8 grid makes the CV domain-size dependent, which quietly turns
    the homogeneity gate into a domain-size gate. On the 2 mm / 100 um
    reference packing an 8x8 grid gives 250 um cells, ~2.5 mean radii across,
    holding a few grains each -- a real density measurement. The same 8x8 on a
    1 mm domain holds well under one grain per cell, so the CV measures
    Poisson noise and lands near 0.45 for perfectly good packings.

    Fixing the cell size in mean radii instead keeps the number comparable
    across domains, and reproduces ncoarse = 8 at the reference geometry.
    """
    return max(4, int(round(Lx / (2.5 * mean_r))))


def grade(cen, rad, Lx, Ly, px, py, mean_r, raster, context=None):
    """Every metric that goes into metadata.json and the acceptance gate.

    POROSITY is measured on the OUTPUT CELL, with its actual periodicity --
    that is the porosity the solver will see, and it is what the caller asked
    for.

    HOMOGENEITY IS MEASURED ON THE PARENT BED, which is x-periodic because
    that is how deposition ran, and y-periodic only if the output is. Judging
    it on the output cell instead confuses a window-cutting artifact with a
    microstructural defect: under an open axis `rasterize` CLIPS every grain
    that overhangs the wall, which severs the solid there (so percolation
    fails) and leaves nothing beyond it for the distance transform (so the
    void against the wall reads enormous). With --periodic y that rejected 2
    of 3 seeds for voids up to 3.4 mean radii and for x-percolation, on
    packings whose interiors were fine.

    A half-pore sitting against an open wall is a property of where the window
    was cut, not of the bed. The bed is the thing that is homogeneous or not,
    and the bed is what has to be a valid REV.
    """
    solid = pl.rasterize(cen, rad, Lx, Ly, raster, raster, px, py)
    porosity = 1.0 - float(solid.mean())

    # The parent bed: x always wraps (deposition did), y follows the output.
    mpx, mpy = True, py
    parent = pl.rasterize(cen, rad, Lx, Ly, raster, raster, mpx, mpy)
    ncoarse = coarse_grid_for(Lx, mean_r)
    cv, _ = pl.uniformity(parent, Lx, Ly, mpx, mpy, ncoarse=ncoarse)
    dist, _, _ = pl.void_map(cen, rad, Lx, Ly, raster, raster, mpx, mpy,
                             context=None if mpy else context)
    max_void = float(dist.max())
    sx, sy, sfrac, sn = pl.percolates(parent)
    px_, py_, pfrac, pn = pl.percolates(~parent)
    coord, ncont, throat = pl.descriptors(cen, rad, Lx, Ly, px, py)

    # HALF-DOMAIN ASYMMETRY -- the largest-scale homogeneity mode, and the one
    # that matters for an REV: "well connected on one side, sparse on the
    # other". local_density_cv does NOT catch it. On an 8x8 grid the CV
    # measures fine-scale clumping, and a packing can sit comfortably inside
    # the CV gate while one half carries 20% more solid than the other.
    #
    # Measured in y as well as x. In y it doubles as the check that the
    # geometric deposition has no compaction-with-depth (which would make the
    # --periodic xy seam indefensible): across seeds the sign flips and a
    # larger crop margin does not shrink it, so what is left is configurational
    # variance rather than a depth drift.
    hy, hx = parent.shape[0] // 2, parent.shape[1] // 2
    bot, top = float(parent[:hy].mean()), float(parent[hy:].mean())
    lef, rig = float(parent[:, :hx].mean()), float(parent[:, hx:].mean())
    asym_y = abs(bot - top) / max(bot + top, 1e-30) * 2.0
    asym_x = abs(lef - rig) / max(lef + rig, 1e-30) * 2.0

    return {
        "porosity_raster": porosity,
        "local_density_cv": cv,
        "max_void_radius_m": max_void,
        "max_void_radius_per_mean_r": max_void / mean_r,
        "solid_percolates_x": sx, "solid_percolates_y": sy,
        "solid_largest_cluster_frac": sfrac, "solid_n_clusters": sn,
        "pore_percolates_x": px_, "pore_percolates_y": py_,
        "pore_largest_cluster_frac": pfrac, "pore_n_clusters": pn,
        "coordination_number": coord, "n_contacts": ncont,
        "throat_gap_m": throat,
        "density_cv_ncoarse": ncoarse,
        "solid_frac_bottom_half": bot, "solid_frac_top_half": top,
        "solid_frac_left_half": lef, "solid_frac_right_half": rig,
        "half_domain_asymmetry_y": asym_y,
        "half_domain_asymmetry_x": asym_x,
        "half_domain_asymmetry": max(asym_x, asym_y),
    }


def deposit_at_porosity(seed, a, px, py, want, verbose=True):
    """Deposit a skeleton whose porosity is ~`want`, bisecting the roll budget.

    THE BISECTION IS PER SEED, and has to be. At a fixed rolling budget,
    different seeds land at genuinely different porosities -- on a 1 mm domain
    (~25 grains) one budget produced 0.356, 0.275 and 0.225 across three seeds.
    Bisecting once on seed A and reusing that budget for seeds B and C misses
    the target by more than the fill pass can repair, and the fill pass can only
    ever REMOVE porosity, so a skeleton that comes out too dense is unrecoverable.

    Porosity falls monotonically as the budget grows, so a plain bisection on
    [0, pi/2] works. Deposition is re-run from the same seed each time, so this
    is deterministic.

    Returns (centres, radii, roll_tol, skeleton_porosity, strip), where
    `strip` is the full deposited bed in window coordinates -- the filler pass
    needs the material outside the window (see packing_lib.void_map).
    """
    margin = a.crop_margin * a.mean_r
    strip_h = a.Ly + 2.0 * margin

    def attempt(rt):
        rng = np.random.default_rng(seed)
        # DEPOSITION ALWAYS WRAPS IN X, whatever the output periodicity. The
        # bed being modelled is wide and the domain is a window cut from it, so
        # a strip with free sides is a generator artifact, not a microstructure
        # anyone wants: grains at x=0 and x=Lx have no neighbours beyond the
        # edge, which opens large voids there AND makes the min-surface stop
        # criterion over-deposit to fill the sparse edge columns, driving the
        # interior too dense. With --periodic y that showed up as porosity
        # landing at 0.22-0.30 against a 0.325 target with voids up to 2.8
        # mean radii, and two of three seeds failing outright.
        #
        # Output x-periodicity (`px`) is a separate question -- it controls
        # edge images, rasterization and the crop -- and is applied below.
        cen, rad = deposit(rng, a.Lx, strip_h, a.mean_r, a.sigma_ln,
                           a.radius_clip_frac, rt, True)
        # Keep the WHOLE strip alongside the window. The grains the crop
        # discards are the bed that physically continues past the domain
        # edge, and the filler pass needs them so an open boundary does not
        # read as empty space (see packing_lib.void_map).
        strip = (cen.copy() - np.array([0.0, margin]), rad.copy())
        cen, rad = crop_window(cen, rad, a.Lx, a.Ly, margin, py)
        solid = pl.rasterize(cen, rad, a.Lx, a.Ly, a.fill_raster,
                             a.fill_raster, px, py)
        return cen, rad, 1.0 - float(solid.mean()), strip

    if a.roll_tol is not None:
        cen, rad, phi, strip = attempt(a.roll_tol)
        return cen, rad, a.roll_tol, phi, strip

    lo, hi = 0.0, math.pi / 2
    c_lo, r_lo, p_lo, s_lo = attempt(lo)
    c_hi, r_hi, p_hi, s_hi = attempt(hi)
    if verbose:
        print(f"    budget bracket: phi({lo:.2f})={p_lo:.4f}  "
              f"phi({hi:.2f})={p_hi:.4f}  want ~{want:.4f}", file=sys.stderr)

    if want >= p_lo:                 # even ballistic is denser than we want
        return c_lo, r_lo, lo, p_lo, s_lo
    if want <= p_hi:                 # even full rolling is looser than we want
        return c_hi, r_hi, hi, p_hi, s_hi

    best = (c_hi, r_hi, hi, p_hi, s_hi)
    for _ in range(a.bisect_iters):
        mid = 0.5 * (lo + hi)
        cen, rad, phi, strip = attempt(mid)
        if abs(phi - want) < abs(best[3] - want):
            best = (cen, rad, mid, phi, strip)
        if phi > want:
            lo = mid
        else:
            hi = mid
        if abs(phi - want) < 0.004:
            return cen, rad, mid, phi, strip
    return best


def build_once(seed, a, px, py, verbose=True):
    """One full attempt: deposit to porosity, fill, seam-heal, grade."""
    rng = np.random.default_rng(seed + 7919)      # fillers get their own stream
    want = min(a.porosity + a.fill_margin, 0.52)

    cen, rad, roll_tol, skeleton_porosity, strip = deposit_at_porosity(
        seed, a, px, py, want, verbose)
    n_matrix = len(rad)

    # SEAM HEAL BEFORE FILLING, not after. Identifying y=0 with y=Ly leaves
    # grains near the seam overlapping, and a rasterized overlap is counted
    # once -- so the geometry reads looser than it really is. Filling against
    # that inflated porosity and only then relaxing the overlaps apart adds
    # solid the fill pass never budgeted for, and porosity lands ~0.02 BELOW
    # target however the rolling budget is bisected. Healing first means the
    # fill pass measures the final geometry.
    seam_shift = 0.0
    if py:
        before = cen.copy()
        # Centres are already inside [0,Ly) (crop_window with py=True), so the
        # wrap is a no-op; what matters is the relax, which resolves the
        # overlaps created by y=0 and y=Ly now being the same line.
        cen[:, 1] = np.mod(cen[:, 1], a.Ly)
        if px:
            cen[:, 0] = np.mod(cen[:, 0], a.Lx)
        # Relax with x WRAPPED, not with the output px. On an open axis
        # pl.relax clips centres to keep every grain wholly inside, which here
        # would bodily yank each x-overhanging grain inward and sever the
        # solid across that boundary -- under --periodic y that made
        # x-percolation fail on 23 of 24 attempts, against 0-1 of 24 for the
        # other three modes. This step is healing the y seam of an x-periodic
        # bed; x must wrap for it.
        pl.relax(cen, rad, a.Lx, a.Ly, True, py,
                 tol_frac=a.tol_frac, max_iter=a.max_iter)
        d = pl.min_image(cen - before, a.Lx, a.Ly, True, py)
        seam_shift = float(np.hypot(d[:, 0], d[:, 1]).max())

    cen, rad, n_fill = fill_voids(rng, cen, rad, a.Lx, a.Ly, px, py,
                                  a.porosity, a.mean_r, a.sigma_ln,
                                  a.radius_clip_frac, a.fill_raster,
                                  # Context is for an OPEN axis only. Where y
                                  # wraps, the cell's own images already ARE
                                  # the surroundings and adding the strip is
                                  # phantom solid that hides real voids (1.45
                                  # reported against a true 2.27). x never
                                  # needs it: deposition wrapped in x, so the
                                  # images are the real neighbours.
                                  context=None if py else strip)

    meta = grade(cen, rad, a.Lx, a.Ly, px, py, a.mean_r, a.raster,
                 context=strip)
    meta.update({
        "n_grains": len(rad),
        "n_matrix_grains": n_matrix,
        "n_filler_grains": n_fill,
        "filler_fraction": n_fill / max(len(rad), 1),
        "roll_tol_rad": roll_tol,
        "skeleton_porosity": skeleton_porosity,
        "seam_max_shift_m": seam_shift,
        "mean_r_m_realized": float(rad.mean()),
        "min_r_m": float(rad.min()), "max_r_m": float(rad.max()),
    })
    if verbose:
        print(f"    roll={roll_tol:.3f}  skeleton phi={skeleton_porosity:.4f}"
              f"  +{n_fill} fillers -> phi={meta['porosity_raster']:.4f}"
              f"  void={meta['max_void_radius_per_mean_r']:.2f}"
              f"  cv={meta['local_density_cv']:.3f}"
              f"  asym={meta['half_domain_asymmetry'] * 100:.0f}%", file=sys.stderr)
    return cen, rad, meta


# =========================================================================

def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--Lx", type=float, required=True, help="domain width [m]")
    p.add_argument("--Ly", type=float, default=None,
                   help="domain height [m] (default: square)")
    p.add_argument("--porosity", type=float, required=True, help="target porosity")
    p.add_argument("--mean-r", dest="mean_r", type=float, default=100e-6,
                   help="mean grain radius [m] (default 100 um)")
    p.add_argument("--sigma-ln", dest="sigma_ln", type=float, default=0.5,
                   help="log-normal sigma of the radius distribution")
    p.add_argument("--radius-clip-frac", dest="radius_clip_frac", type=float,
                   default=1.0, help="clip radii to (1 +/- f)*mean_r")
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--out", required=True,
                   help="output directory. A periodic_<mode>/ level is inserted "
                        "above the leaf name (see --no-periodic-subdir)")
    p.add_argument("--no-periodic-subdir", dest="periodic_subdir",
                   action="store_false",
                   help="write straight to --out instead of grouping by "
                        "periodicity")

    p.add_argument("--periodic", choices=["x", "y", "xy", "none"], default="x",
                   help="which axes wrap (default: x)")
    p.add_argument("--roll-tol", dest="roll_tol", type=float, default=None,
                   help="rolling budget [rad]; default: bisected to hit the "
                        "target porosity")
    p.add_argument("--fill-margin", dest="fill_margin", type=float, default=0.03,
                   help="how far ABOVE target the deposited skeleton should "
                        "land, leaving room for the filler pass (default 0.03)")
    p.add_argument("--crop-margin", dest="crop_margin", type=float, default=3.0,
                   help="strip margin above and below the window, in mean "
                        "radii, discarded as floor/free-surface artifact")

    p.add_argument("--max-void-ratio", dest="max_void_ratio", type=float,
                   default=None,
                   help="reject if the largest inscribed pore circle exceeds "
                        "this many mean grain radii. Default: scaled with the "
                        "target porosity (see auto_void_ratio)")
    p.add_argument("--max-density-cv", dest="max_density_cv", type=float,
                   default=None,
                   help="reject above this local-density CV. Default: scaled "
                        "with the target porosity (see auto_density_cv)")
    p.add_argument("--max-asymmetry", dest="max_asymmetry", type=float,
                   default=0.22,
                   help="reject when one half of the domain carries this much "
                        "more solid than the other, in either axis "
                        "(default 0.22)")
    p.add_argument("--porosity-tol", dest="porosity_tol", type=float,
                   default=0.01,
                   help="reject if achieved porosity misses the target by more "
                        "than this (default 0.01). The target is the point of "
                        "the exercise, so it is gated like the rest.")
    p.add_argument("--no-percolation-gate", dest="percolation_gate",
                   action="store_false",
                   help="do not require the solid to percolate in both axes")
    p.add_argument("--max-tries", dest="max_tries", type=int, default=24,
                   help="re-seed this many times before giving up (default 24). "
                        "Passing the gate is a per-seed lottery and the odds "
                        "worsen with porosity: at phi 0.40 a loose bed forms "
                        "vertical chains easily but a horizontal path across "
                        "the full width often fails, and one sweep needed 14 "
                        "attempts. An attempt costs ~0.15 s at 2 mm.")

    p.add_argument("--raster", type=int, default=1024,
                   help="raster for the reported metrics")
    p.add_argument("--fill-raster", dest="fill_raster", type=int, default=512,
                   help="raster inside the fill loop (speed)")
    p.add_argument("--bisect-iters", dest="bisect_iters", type=int, default=12)
    p.add_argument("--tol-frac", dest="tol_frac", type=float, default=1e-4)
    p.add_argument("--max-iter", dest="max_iter", type=int, default=4000)
    p.add_argument("--no-preview", dest="preview", action="store_false")
    a = p.parse_args(argv)

    if a.Ly is None:
        a.Ly = a.Lx
    px = a.periodic in ("x", "xy")
    py = a.periodic in ("y", "xy")

    # Group output by periodicity. The same packing parameters under different
    # wrapping are different microstructures and belong in different runs, but
    # they naturally get the same descriptive name -- so without this the
    # second one silently overwrites the first. Insert the level ABOVE the leaf
    # so the leaf keeps meaning "this packing":
    #   --out inputs/packings/phi0.325_seed4  --periodic x
    #     -> inputs/packings/periodic_x/phi0.325_seed4
    if a.periodic_subdir:
        parent, leaf = os.path.split(os.path.normpath(a.out))
        a.out = os.path.join(parent, f"periodic_{a.periodic}", leaf)
    if a.max_void_ratio is None:
        a.max_void_ratio = auto_void_ratio(a.porosity)
    if a.max_density_cv is None:
        a.max_density_cv = auto_density_cv(a.porosity)
    print(f"  gates for target phi={a.porosity:.3f}:  "
          f"max_void_ratio {a.max_void_ratio:.2f}   "
          f"local_density_cv {a.max_density_cv:.3f}   "
          f"asymmetry {a.max_asymmetry * 100:.0f}%", file=sys.stderr)

    print(f"Gravity-style packing: {a.Lx * 1e3:.3f} x {a.Ly * 1e3:.3f} mm, "
          f"target phi={a.porosity}, mean r={a.mean_r * 1e6:.1f} um, "
          f"sigma_ln={a.sigma_ln}, periodic={a.periodic}", file=sys.stderr)

    best = None
    for attempt in range(a.max_tries):
        seed = a.seed + 1000 * attempt
        print(f"  attempt {attempt + 1}/{a.max_tries} (seed {seed})", file=sys.stderr)
        cen, rad, meta = build_once(seed, a, px, py)
        meta["seed"] = seed
        meta["attempt"] = attempt + 1
        meta["porosity_achieved"] = meta["porosity_raster"]
        bad = pl.accept_reasons(meta, a.max_void_ratio, a.max_density_cv,
                                a.percolation_gate)
        # The porosity the caller asked for is the point of the whole exercise,
        # so it is part of the gate. Without this the generator happily emitted
        # a 0.2249 packing against a 0.3250 target.
        d_phi = abs(meta["porosity_achieved"] - a.porosity)
        if d_phi > a.porosity_tol:
            bad.append(f"porosity {meta['porosity_achieved']:.4f} misses target "
                       f"{a.porosity:.4f} by {d_phi:.4f} > {a.porosity_tol}")
        if meta["half_domain_asymmetry"] > a.max_asymmetry:
            bad.append(f"half-domain asymmetry "
                       f"{meta['half_domain_asymmetry'] * 100:.1f}% "
                       f"(x {meta['half_domain_asymmetry_x'] * 100:.1f}, "
                       f"y {meta['half_domain_asymmetry_y'] * 100:.1f}) "
                       f"> {a.max_asymmetry * 100:.0f}%")
        if not bad:
            best = (cen, rad, meta, bad, (0.0, 0.0))
            print("    ACCEPTED", file=sys.stderr)
            break
        print("    rejected: " + "; ".join(bad), file=sys.stderr)
        # Keep the least-bad attempt to report on failure: porosity error
        # first, since that is the hard requirement, then void size.
        rank = (round(d_phi / max(a.porosity_tol, 1e-9), 3),
                meta["max_void_radius_per_mean_r"])
        if best is None or rank < best[4]:
            best = (cen, rad, meta, bad, rank)

    cen, rad, meta, bad, _ = best
    if bad:
        print(f"\n❌ No packing passed the gate in {a.max_tries} tries.",
              file=sys.stderr)
        print(f"   Best attempt (seed {meta['seed']}) still fails:", file=sys.stderr)
        for b in bad:
            print(f"     - {b}", file=sys.stderr)
        print("   Nothing was written. Loosen --max-void-ratio / "
              "--max-density-cv if these limits are wrong for this porosity, "
              "or raise --max-tries.", file=sys.stderr)
        return 1

    meta.update({
        "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "generator": "preprocess/generate_packing_gravity.py",
        "algorithm": "geometric drop-and-roll deposition, interior window "
                     "crop, largest-void filling",
        "Lx": a.Lx, "Ly": a.Ly,
        "porosity_target": a.porosity,
        "porosity_achieved": meta["porosity_raster"],
        "mean_r_m_requested": a.mean_r,
        "sigma_ln": a.sigma_ln,
        "radius_clip_frac": a.radius_clip_frac,
        "fill_margin": a.fill_margin,
        "crop_margin_mean_r": a.crop_margin,
        "gate_max_void_ratio": a.max_void_ratio,
        "gate_max_density_cv": a.max_density_cv,
        "gate_max_asymmetry": a.max_asymmetry,
        "gate_porosity_tol": a.porosity_tol,
        "raster_nx": a.raster, "raster_ny": a.raster,
    })

    dat = pl.write_outputs(a.out, cen, rad, a.Lx, a.Ly, meta, px, py, a.preview)

    print(f"\n✅ {a.out}", file=sys.stderr)
    print(f"   porosity      {meta['porosity_achieved']:.4f} "
          f"(target {a.porosity:.4f})", file=sys.stderr)
    print(f"   grains        {meta['n_grains']}  "
          f"({meta['n_matrix_grains']} matrix + {meta['n_filler_grains']} filler)",
          file=sys.stderr)
    print(f"   mean radius   {meta['mean_r_m_realized'] * 1e6:.2f} um "
          f"(requested {a.mean_r * 1e6:.2f})", file=sys.stderr)
    print(f"   largest void  {meta['max_void_radius_per_mean_r']:.2f} mean radii "
          f"(gate {a.max_void_ratio:.2f})", file=sys.stderr)
    print(f"   density cv    {meta['local_density_cv']:.3f} "
          f"(gate {a.max_density_cv:.3f})", file=sys.stderr)
    print(f"   coordination  {meta['coordination_number']:.2f}", file=sys.stderr)
    print(f"   half-domain   {meta['half_domain_asymmetry'] * 100:.1f}% asymmetry "
          f"(x {meta['half_domain_asymmetry_x'] * 100:.1f}, "
          f"y {meta['half_domain_asymmetry_y'] * 100:.1f}; gate "
          f"{a.max_asymmetry * 100:.0f})", file=sys.stderr)
    print(f"   rows written  {meta['n_grain_rows_emitted']} "
          f"({meta['n_edge_images']} edge images)  -> {dat}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
