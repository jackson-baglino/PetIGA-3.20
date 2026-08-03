#!/usr/bin/env python3
"""
generate_packing.py — periodic, jammed, polydisperse 2D disk packings for the
snow-thermal study.

Why not the old generateCircularGrains.py
-----------------------------------------
That script did random sequential addition (RSA) with a "grow the new disk
until it touches its nearest neighbour" step. Two problems for this study:

1. RSA saturates near solid fraction 0.547. The targets here run to 0.80
   (porosity 0.20), which is why it needed MAX_ATTEMPTS = 2e6 and still had to
   distort radii to get there. The grow-to-touch step overwrites the log-normal
   radius it just sampled, so the realized size distribution is not the one
   requested -- and grain size is a control variable in this study.
2. It set ALLOW_CROSS_BOUNDARY without wrapping, so boundary grains were
   clipped with no matching partial grain on the opposite face. Both the solver
   (-periodic 1) and effective_thermal_cond (periodic cell problem,
   setup.c:45,51,58) treat the domain as periodic, so that microstructure was
   discontinuous across every face it touched.

What this does instead
----------------------
Collective rearrangement (Jodrey-Tory family) under minimum-image periodicity,
with the radii grown gradually rather than the centres re-drawn:

  1. Sample N log-normal radii; N is chosen so their total area hits the target
     solid fraction exactly.
  2. Place all centres uniformly at random, overlaps allowed.
  3. Scale all radii by s < 1, then relax: repeatedly push overlapping pairs
     apart along their centre line (minimum image) until no pair overlaps.
  4. Raise s and re-relax. On failure, halve the step and retry.
  5. Stop when s reaches 1 (target porosity met) or the step underflows
     (jammed short of target -- reported, not silently accepted).

Contacts are then natural: at jamming, neighbours are touching because the
configuration cannot be compressed further, not because a radius was forced.
The realized radius distribution is the sampled one, uniformly scaled.

The contact gap, and why it is not optional in 2D
-------------------------------------------------
At EXACT tangency a jammed 2D packing has a DISCONNECTED pore space. This is
not a rasterization artifact: jamming means the contact network spans the
domain in both directions, and in 2D a spanning contact network cuts the
complement into isolated cells. Measured on a 2 mm, 433-grain, porosity-0.30
packing: the pore fell into 1412 clusters with the largest holding 11% of the
void, and dilating or eroding the solid did not change that.

The consequence is that a perfectly tangent 2D packing cannot transport vapour
across the domain at all, which defeats the purpose of the simulation. (This
obstruction is specific to 2D -- in 3D both phases percolate simultaneously.)

--contact-gap jams the grains against a radius r + gap/2 while reporting the
true r. A POSITIVE gap leaves a uniform surface-to-surface separation at every
contact; ZERO gives exact tangency; a NEGATIVE gap makes neighbours overlap.
Porosity is computed from the true radii and is unaffected either way.

WHY THE DEFAULT IS ZERO -- EXACT TANGENCY (changed 2026-08-01)
--------------------------------------------------------------
The initial condition gives each grain its own diffuse field

    phi_k = 0.5 - 0.5*tanh( 0.5*(r_k - R_k)/eps )

which crosses 0.5 at exactly r_k = R_k, the sharp grain surface, FOR ANY eps.
"In contact" therefore has an eps-independent meaning: the grains are tangent,
centre distance = R_1 + R_2, so that at the contact point phi_1 = phi_2 = 0.5
and the additive sum is 1. Nothing needs calibrating and the radii never depend
on eps -- the geometry is fixed by the sharp surfaces, and only the width of the
diffuse skin around them changes.

Do NOT use a negative gap here. Overlapping grains are what the UNION form
(-ic_grain_union 1) was built for; with the additive sum used for packings, an
overlap has to be calibrated against a target neck (the Molaro axisymmetric
geometries do this only because they reproduce a specific experiment's initial
neck) and that calibration is eps-dependent, which reintroduces exactly the
eps sensitivity contact was meant to remove.

WHY NOT A POSITIVE GAP
----------------------
A positive gap bought vapour percolation at t=0 and it is not worth the price:

  * It does not last. Sintering closes the gap within the first part of any run,
    so a positive gap only defers the disconnection -- it does not avoid it.
  * It corrupts the observable this study exists to measure. k_eff is largely a
    function of GRAIN CONNECTIVITY, and the rise of k_eff through a run IS the
    growth of that connectivity. Starting from a pack whose grains do not touch
    means the early k_eff rise is partly a numerical gap closing rather than
    real neck growth, which skews the very trend being reported.
  * With a 4 um gap and a diffuse interface 9.2-18.4 um wide, the phase field
    bridges the gap anyway -- so conduction at t=0 was carried by the interface
    representation, not by the microstructure. Measured: k_eff differed by
    31-39% between eps = 1 um and 2 um on an identical packing
    (studies/snow_thermal/verification/).
  * For Enceladus specifically, the vapour budget is small: at tiger-stripe
    pressures and temperatures there is little vapour to transport, so paying
    for pore percolation with a disconnected solid is the wrong trade.

The way to prolong vapour transport is therefore NOT to hold the grains apart,
but to SELECT configurations that retain large pores -- see
preprocess/select_packing_seeds.py, which generates many seeds and ranks them on
pore retention. This is a real limitation of 2D and belongs in the write-up: in
3D both phases percolate simultaneously and no such choice is forced.

Measured on 2 mm packings, mean radius 45 um, porosity 0.20 to 0.40 (see the
table below): the pore reconnects between 2 and 4 um of gap, and 4 um gives a
single cluster holding 100% of the void at EVERY porosity. The transition is
set by grain-scale contact geometry, not by the global density, which is why
it barely moves with porosity. Hence the 4 um default.

  gap    porosity 0.25 .. 0.40
  0 um   pore blocked  (largest cluster 10-50% of void), solid percolates
  2 um   pore blocked  (15-66%),                         solid broken up
  4 um   pore PERCOLATES (100%, one cluster),            solid broken up
  6 um   as 4 um

Pore and solid percolation are mutually exclusive throughout. That is the
correct starting point for this study rather than a defect: an unsintered
packing has no solid conduction path, and k_eff should RISE from near the air
value as vapour deposition builds necks. A 4 um gap is also narrower than the
5-95% diffuse band (about 6*eps = 5.2 um at eps = 0.86 um), so neighbouring
grains' interfaces already overlap and necks form immediately.

LOWEST ACHIEVABLE POROSITY. Jamming against r + gap/2 costs solid fraction:
with 2D random close packing near 0.84, the target must satisfy
(1-phi)*((r+gap/2)/r)^2 <~ 0.84. At r = 45 um and gap = 4 um that is
phi >~ 0.23, and porosity 0.20 indeed fails to jam (it stops at s = 0.975 and
overshoots to 0.24). To go below ~0.23, raise --mean_r_m: the penalty is
(1 + gap/2r)^2, so r = 90 um admits phi >~ 0.20.

Checks (all failures are reported in metadata.json and on stderr):
  - achieved area porosity vs target
  - zero overlaps under minimum image
  - pore phase percolates (else vapour cannot diffuse -- see above)
  - solid phase percolates (informational: a gapped packing has no
    solid-contact conduction path at t=0, which is the physically correct
    starting point -- k_eff should rise as necks form)

Usage
-----
  ./venv_enceladus/bin/python3 preprocess/generate_packing.py \
      --Lx 2e-3 --Ly 2e-3 --porosity 0.30 --seed 1 \
      --out inputs/packings/phi0.30_seed1
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
from datetime import datetime, timezone

import numpy as np
from scipy.spatial import cKDTree, Delaunay
from scipy import ndimage


# =========================================================================
# Radius sampling
# =========================================================================

def sample_radii(rng, Lx, Ly, porosity, mean_r, sigma_ln, clip_frac):
    """Draw log-normal radii until their total area hits the solid target.

    mu is set so E[r] == mean_r exactly: for a log-normal,
    E[r] = exp(mu + sigma^2/2), hence mu = ln(mean_r) - sigma^2/2.
    Radii are clipped to (1 +/- clip_frac)*mean_r so a single freak draw
    cannot dominate the domain.
    """
    target_area = (1.0 - porosity) * Lx * Ly
    mu = math.log(mean_r) - 0.5 * sigma_ln ** 2
    lo, hi = (1.0 - clip_frac) * mean_r, (1.0 + clip_frac) * mean_r
    if lo <= 0.0:
        lo = 1e-3 * mean_r

    radii, area = [], 0.0
    while area < target_area:
        r = float(np.clip(rng.lognormal(mu, sigma_ln), lo, hi))
        radii.append(r)
        area += math.pi * r * r

    radii = np.array(radii)
    # Drop the last grain if doing so lands closer to the target area.
    if len(radii) > 1:
        over = area - target_area
        last = math.pi * radii[-1] ** 2
        if over > 0.5 * last:
            radii = radii[:-1]
    return radii


# =========================================================================
# Periodic overlap relaxation
# =========================================================================

# Set once from --non-periodic before any packing work. A module-level cell
# rather than a threaded argument because _worst_overlap is called from four
# places and the alternative is plumbing a boolean through all of them.
PERIODIC = [True]


def _worst_overlap(centres, radii, Lx, Ly):
    """Return (max_overlap, pairs, deltas) for all overlapping pairs.

    Overlap is (r_i + r_j) - d under the minimum image. cKDTree's boxsize
    argument does the periodic wrapping for us, so no 3x3 replication is
    needed; it requires coordinates inside [0, L).
    """
    # boxsize=None drops the minimum image, so grains no longer interact with
    # their images across the wall -- required when the cell is a closed box
    # with zero-flux walls rather than a periodic tile.
    tree = (cKDTree(centres, boxsize=[Lx, Ly]) if PERIODIC[0]
            else cKDTree(centres))
    rmax = float(radii.max())
    pairs = np.array(list(tree.query_pairs(2.0 * rmax)), dtype=np.int64)
    if pairs.size == 0:
        return 0.0, pairs, None, None

    i, j = pairs[:, 0], pairs[:, 1]
    d = centres[j] - centres[i]
    d[:, 0] -= Lx * np.round(d[:, 0] / Lx)      # minimum image
    d[:, 1] -= Ly * np.round(d[:, 1] / Ly)
    dist = np.hypot(d[:, 0], d[:, 1])
    overlap = (radii[i] + radii[j]) - dist

    keep = overlap > 0.0
    if not np.any(keep):
        return 0.0, pairs[:0], None, None
    return float(overlap[keep].max()), pairs[keep], d[keep], dist[keep]


def relax(centres, radii, Lx, Ly, tol_frac, max_iter, damping=0.5):
    """Push overlapping pairs apart until no overlap exceeds tol.

    Each overlapping pair is separated along its centre line, splitting the
    correction between the two disks. Displacements from all pairs touching a
    given disk are accumulated before moving anything, so the sweep is
    Jacobi-style (order-independent) rather than Gauss-Seidel.
    """
    tol = tol_frac * float(radii.min())
    for it in range(max_iter):
        worst, pairs, delta, dist = _worst_overlap(centres, radii, Lx, Ly)
        if worst <= tol:
            return True, it, worst

        i, j = pairs[:, 0], pairs[:, 1]
        overlap = (radii[i] + radii[j]) - dist
        # Unit vector i->j; if two centres coincide, pick a random direction.
        with np.errstate(invalid="ignore", divide="ignore"):
            u = delta / dist[:, None]
        bad = ~np.isfinite(u).all(axis=1)
        if np.any(bad):
            ang = np.random.default_rng(it).uniform(0, 2 * np.pi, bad.sum())
            u[bad] = np.stack([np.cos(ang), np.sin(ang)], axis=1)

        shift = damping * overlap[:, None] * u
        disp = np.zeros_like(centres)
        np.add.at(disp, i, -0.5 * shift)
        np.add.at(disp, j, +0.5 * shift)

        centres += disp
        if PERIODIC[0]:
            centres[:, 0] = np.mod(centres[:, 0], Lx)   # stay inside [0, L)
            centres[:, 1] = np.mod(centres[:, 1], Ly)
        else:
            # Keep every grain WHOLLY inside the box, so nothing is clipped by
            # the wall. This is what makes a zero-flux run legitimate: the
            # multi_grains IC only applies minimum-image wrapping under
            # -periodic 1, so a grain straddling the boundary would otherwise
            # be cut, leaving a flat artificial face from t = 0.
            centres[:, 0] = np.clip(centres[:, 0], radii, Lx - radii)
            centres[:, 1] = np.clip(centres[:, 1], radii, Ly - radii)

    worst, _, _, _ = _worst_overlap(centres, radii, Lx, Ly)
    return worst <= tol, max_iter, worst


def jam(centres, radii, Lx, Ly, tol_frac, max_iter, s_start, verbose=True):
    """Grow radii from s_start to 1.0, relaxing overlaps at each step.

    Returns the largest scale factor that could be relaxed overlap-free.
    """
    s_ok = 0.0
    best = centres.copy()
    s, step = s_start, 0.5 * (1.0 - s_start)

    # Seed: relax at the (easy) starting scale.
    ok, _, _ = relax(centres, radii * s, Lx, Ly, tol_frac, max_iter)
    if not ok:
        raise RuntimeError(f"cannot relax even at s={s:.3f}; lower --s-start")
    s_ok, best = s, centres.copy()

    while step > 1e-4:
        s_try = min(1.0, s_ok + step)
        trial = best.copy()
        ok, iters, worst = relax(trial, radii * s_try, Lx, Ly, tol_frac, max_iter)
        if ok:
            s_ok, best = s_try, trial
            if verbose:
                print(f"    s = {s_ok:.4f} ok ({iters} sweeps)", file=sys.stderr)
            if s_ok >= 1.0:
                break
        else:
            step *= 0.5
            if verbose:
                print(f"    s = {s_try:.4f} failed (worst overlap "
                      f"{worst:.2e}); step -> {step:.4f}", file=sys.stderr)

    return best, s_ok


# =========================================================================
# Rasterization and percolation
# =========================================================================

def rasterize(centres, radii, Lx, Ly, nx, ny):
    """Boolean solid mask on an nx-by-ny grid.

    Honours PERIODIC. Under PERIODIC[0] = False the mask is a plain box: grains
    whose centres lie OUTSIDE [0,L] contribute only the part that falls inside,
    and nothing wraps.

    That distinction is not cosmetic. A window cut from a larger pack keeps every
    grain overlapping it, so centres span [-R, L+R]; rasterising those with the
    periodic branch wraps each overhanging grain onto the OPPOSITE face, adding
    solid that is not there and removing it from where it is. Porosity,
    percolation and coordination are all corrupted by it.
    """
    dx, dy = Lx / nx, Ly / ny
    solid = np.zeros((ny, nx), dtype=bool)
    xs = (np.arange(nx) + 0.5) * dx
    ys = (np.arange(ny) + 0.5) * dy

    for (cx, cy), r in zip(centres, radii):
        # Only touch the bounding box of this grain.
        i0, i1 = int(np.floor((cx - r) / dx)) - 1, int(np.ceil((cx + r) / dx)) + 1
        j0, j1 = int(np.floor((cy - r) / dy)) - 1, int(np.ceil((cy + r) / dy)) + 1

        if PERIODIC[0]:
            ii = np.arange(i0, i1 + 1)
            jj = np.arange(j0, j1 + 1)
            ddx = xs[np.mod(ii, nx)] - cx
            ddy = ys[np.mod(jj, ny)] - cy
            ddx -= Lx * np.round(ddx / Lx)      # minimum image
            ddy -= Ly * np.round(ddy / Ly)
            inside = (ddx[None, :] ** 2 + ddy[:, None] ** 2) <= r * r
            if inside.any():
                solid[np.ix_(np.mod(jj, ny), np.mod(ii, nx))] |= inside
        else:
            ii = np.arange(max(i0, 0), min(i1, nx - 1) + 1)
            jj = np.arange(max(j0, 0), min(j1, ny - 1) + 1)
            if ii.size == 0 or jj.size == 0:
                continue                        # grain lies wholly outside
            ddx = xs[ii] - cx
            ddy = ys[jj] - cy
            inside = (ddx[None, :] ** 2 + ddy[:, None] ** 2) <= r * r
            if inside.any():
                solid[np.ix_(jj, ii)] |= inside
    return solid


def percolates(mask):
    """Does `mask` connect to its own periodic translate in x and in y?

    Tiles the mask 2x2 and labels it WITHOUT periodicity. If a pixel and the
    copy one period to the right carry the same label, the cluster joins its
    own image, i.e. it wraps -- which is what percolation means on a torus.

    Returns (perc_x, perc_y, largest_cluster_fraction, n_clusters).
    """
    ny, nx = mask.shape
    if not mask.any():
        return False, False, 0.0, 0
    tiled = np.tile(mask, (2, 2))
    lab, _ = ndimage.label(tiled)          # 4-connectivity
    centre = lab[:ny, :nx]
    sel = centre > 0
    perc_x = bool(np.any(centre[sel] == lab[:ny, nx:2 * nx][sel]))
    perc_y = bool(np.any(centre[sel] == lab[ny:2 * ny, :nx][sel]))

    single, n = ndimage.label(mask)
    sizes = np.bincount(single.ravel())[1:]
    frac = float(sizes.max()) / float(mask.sum()) if n else 0.0
    return perc_x, perc_y, frac, int(n)


# =========================================================================
# Structural descriptors
# =========================================================================

def descriptors(centres, radii, Lx, Ly, gap=0.0, contact_tol_frac=0.02):
    """Coordination number and throat (gap) statistics.

    Contacts: pairs whose surface gap is within contact_tol_frac * r_min of
    zero. Throats: gaps along Delaunay edges of the periodically replicated
    centres, restricted to edges with at least one endpoint in the base tile.
    """
    tol = contact_tol_frac * float(radii.min())

    # boxsize=None drops the minimum image, so grains no longer interact with
    # their images across the wall -- required when the cell is a closed box
    # with zero-flux walls rather than a periodic tile.
    tree = (cKDTree(centres, boxsize=[Lx, Ly]) if PERIODIC[0]
            else cKDTree(centres))
    pairs = np.array(list(tree.query_pairs(2.0 * radii.max() + gap + tol)), dtype=np.int64)
    n_contact = 0
    if pairs.size:
        i, j = pairs[:, 0], pairs[:, 1]
        d = centres[j] - centres[i]
        d[:, 0] -= Lx * np.round(d[:, 0] / Lx)
        d[:, 1] -= Ly * np.round(d[:, 1] / Ly)
        sep = np.hypot(d[:, 0], d[:, 1]) - (radii[i] + radii[j])
        n_contact = int(np.sum(np.abs(sep - gap) <= tol))
    coord = 2.0 * n_contact / len(radii)

    # Throats via Delaunay on a 3x3 replication (so wrap-around edges exist).
    n = len(centres)
    offs = np.array([[a * Lx, b * Ly] for a in (-1, 0, 1) for b in (-1, 0, 1)])
    rep = (centres[None, :, :] + offs[:, None, :]).reshape(-1, 2)
    rep_r = np.tile(radii, len(offs))
    rep_id = np.tile(np.arange(n), len(offs))
    base = np.full(len(rep), False)
    base[4 * n:5 * n] = True                      # the (0,0) offset block

    tri = Delaunay(rep)
    edges = set()
    for s in tri.simplices:
        for a in range(3):
            for b in range(a + 1, 3):
                p, q = int(s[a]), int(s[b])
                if base[p] or base[q]:
                    edges.add((min(p, q), max(p, q)))
    if edges:
        e = np.array(sorted(edges))
        gaps = (np.linalg.norm(rep[e[:, 0]] - rep[e[:, 1]], axis=1)
                - (rep_r[e[:, 0]] + rep_r[e[:, 1]]))
        gaps = gaps[gaps > 0.0]                   # true openings only
    else:
        gaps = np.array([])

    stats = {}
    if gaps.size:
        for p in (5, 25, 50, 75, 95):
            stats[f"p{p}"] = float(np.percentile(gaps, p))
        stats["mean"] = float(gaps.mean())
        stats["n"] = int(gaps.size)
    return coord, n_contact, stats, rep_id


# =========================================================================
# Output
# =========================================================================

def write_outputs(out, centres, radii, Lx, Ly, meta, preview=True):
    os.makedirs(out, exist_ok=True)

    # EXPLICIT EDGE IMAGES.
    #
    # The pack is jammed under minimum-image periodicity, so a grain near a wall
    # genuinely continues on the opposite face. The multi_grains IC only
    # reconstructs that under -periodic 1; under -periodic 0 (zero-flux walls,
    # which is what the phase-field, temperature and vapour equations want) it
    # would simply CUT the grain and leave a flat artificial face.
    #
    # So the wrapping is baked into the file instead of left to the solver: any
    # grain reaching past an edge is emitted a second time at its periodic image,
    # giving a centre list that spans [-R, Lx+R] x [-R, Ly+R]. A grain then needs
    # only to be PARTLY inside the domain, which is the point -- these packings
    # represent a window cut from a much larger pack, and deleting the grains
    # that straddle the edge would sever the conduction paths that leave the
    # domain and bias k_eff low. Connectivity across the boundary is exactly what
    # k_eff is most sensitive to.
    #
    # Because the images are the true periodic partners, the microstructure also
    # stays periodic across every face, which removes the seam mismatch the
    # corrector's periodic BCs would otherwise see. Not required -- Calonne et
    # al. (2011) apply periodic BCs to plainly non-periodic tomographic images --
    # but it costs nothing here and removes one source of window bias.
    emitted = []
    for (x, y), r in zip(centres, radii):
        for ox in (-Lx, 0.0, Lx):
            for oy in (-Ly, 0.0, Ly):
                cx, cy = x + ox, y + oy
                # keep any copy whose disk overlaps the domain at all
                if (-r <= cx <= Lx + r) and (-r <= cy <= Ly + r):
                    emitted.append((cx, cy, r))
    meta["n_grain_rows_emitted"] = len(emitted)
    meta["n_edge_images"] = len(emitted) - len(radii)

    dat = os.path.join(out, "grains.dat")
    with open(dat, "w") as f:
        f.write(f"# periodic 2D packing, generated {meta['generated_utc']}\n")
        f.write(f"# porosity_achieved = {meta['porosity_achieved']:.6f}"
                f"  (target {meta['porosity_target']:.6f})\n")
        f.write(f"# {len(radii)} distinct grains; {len(emitted)} rows including "
                f"{len(emitted)-len(radii)} periodic edge images.\n")
        f.write("# Centres span [-R, L+R]: a grain need only be PARTLY inside.\n")
        f.write("# Do NOT drop the out-of-box rows -- they carry the conduction\n")
        f.write("# paths that cross the domain edge, and k_eff depends on them.\n")
        f.write("# header row is 'Lx Ly'; grain rows are 'x y r', all in metres\n")
        f.write(f"{Lx:.12e} {Ly:.12e}\n")
        for (x, y, r) in emitted:
            f.write(f"{x:.12e} {y:.12e} {r:.12e}\n")

    with open(dat, "rb") as f:
        meta["grains_dat_sha256"] = hashlib.sha256(f.read()).hexdigest()
    with open(os.path.join(out, "metadata.json"), "w") as f:
        json.dump(meta, f, indent=2)

    if preview:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.patches import Circle
            fig, ax = plt.subplots(figsize=(6.4, 6.4 * Ly / Lx))
            # Draw the margin OUTSIDE the domain too, so the grains that hang
            # past an edge are visible -- they are written to grains.dat and the
            # solver sees them, and they carry the conduction paths that leave
            # the domain. A preview cropped to [0,L] hides exactly the part of
            # the configuration that k_eff is most sensitive to.
            pad = float(radii.max())
            for (x, y, r), is_img in ((e, out_of_box) for e, out_of_box in
                                      ((e, not (0 <= e[0] <= Lx and 0 <= e[1] <= Ly))
                                       for e in emitted)):
                ax.add_patch(Circle((x, y), r,
                                    facecolor="#c9d9e6" if is_img else "#7fb3d5",
                                    edgecolor="#1b4f72",
                                    lw=0.3, alpha=0.55 if is_img else 1.0,
                                    zorder=2 if is_img else 3))
            ax.add_patch(plt.Rectangle((0, 0), Lx, Ly, fill=False,
                                       edgecolor="#c0392b", lw=1.8, zorder=4))
            ax.set_xlim(-pad, Lx + pad); ax.set_ylim(-pad, Ly + pad)
            ax.set_aspect("equal")
            ax.set_xlabel("x [m]", fontsize=8); ax.set_ylabel("y [m]", fontsize=8)
            ax.tick_params(labelsize=7)
            ax.set_title(
                f"phi={meta['porosity_achieved']:.4f}  N={len(radii)}  "
                f"seed={meta['seed']}  L={Lx*1e3:.2f} mm\n"
                f"{len(emitted)-len(radii)} edge images (pale)   "
                f"solid in largest cluster {meta['solid_largest_cluster_frac']*100:.0f}%   "
                f"pore {meta['pore_largest_cluster_frac']*100:.0f}%",
                fontsize=8.5)
            fig.savefig(os.path.join(out, "preview.png"), dpi=140,
                        bbox_inches="tight")
            plt.close(fig)
        except Exception as e:                       # preview is optional
            print(f"  (preview skipped: {e})", file=sys.stderr)


# =========================================================================

def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Periodic jammed polydisperse 2D disk packing.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--Lx", type=float, required=True, help="domain width [m]")
    ap.add_argument("--Ly", type=float, default=None,
                    help="domain height [m] (default: square)")
    ap.add_argument("--porosity", type=float, required=True,
                    help="target AREA porosity (void fraction)")
    ap.add_argument("--mean_r_m", type=float, default=45e-6,
                    help="mean grain radius [m]")
    ap.add_argument("--sigma_ln", type=float, default=0.5,
                    help="log-normal shape parameter of the radii")
    ap.add_argument("--radius_clip_frac", type=float, default=1.0,
                    help="clip radii to (1 +/- f)*mean_r_m")
    ap.add_argument("--contact-gap", type=float, default=0.0,
                    help="surface-to-surface separation held at every contact "
                         "[m]. 0 (the default) is EXACT TANGENCY: the grains "
                         "touch where each one's phi crosses 0.5, which is "
                         "eps-independent. Positive holds them apart; negative "
                         "overlaps them and then needs an eps-dependent neck "
                         "calibration -- see the module docstring")
    ap.add_argument("--non-periodic", action="store_true",
                    help="build a CLOSED-BOX packing: no minimum image, and "
                         "every grain kept wholly inside the walls so nothing "
                         "is clipped. Use with -periodic 0 runs, where the "
                         "multi_grains IC does not wrap grains and a straddling "
                         "grain would be cut")
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--raster", type=int, default=1024,
                    help="grid for the porosity/percolation checks")
    ap.add_argument("--max-iter", type=int, default=4000,
                    help="relaxation sweeps allowed per growth step")
    ap.add_argument("--tol-frac", type=float, default=1e-4,
                    help="overlap tolerance as a fraction of the smallest radius")
    ap.add_argument("--s-start", type=float, default=0.60,
                    help="initial radius scale factor")
    ap.add_argument("--no-preview", action="store_true")
    args = ap.parse_args(argv)

    Lx = args.Lx
    Ly = args.Ly if args.Ly is not None else args.Lx
    PERIODIC[0] = not args.non_periodic
    rng = np.random.default_rng(args.seed)

    print(f"[{args.out}] target porosity {args.porosity:.3f} on "
          f"{Lx*1e3:.2f} x {Ly*1e3:.2f} mm, seed {args.seed}", file=sys.stderr)

    radii = sample_radii(rng, Lx, Ly, args.porosity, args.mean_r_m,
                         args.sigma_ln, args.radius_clip_frac)
    n = len(radii)
    print(f"  sampled {n} grains", file=sys.stderr)

    # Jam against an inflated radius so a uniform surface gap survives at every
    # contact; porosity is always computed from the TRUE radii.
    gap = args.contact_gap
    centres = np.column_stack([rng.uniform(0, Lx, n), rng.uniform(0, Ly, n)])
    centres, s_ok = jam(centres, radii + 0.5 * gap, Lx, Ly, args.tol_frac,
                        args.max_iter, args.s_start)
    radii = radii * s_ok
    radii_eff = radii + 0.5 * gap

    # ---- checks ----------------------------------------------------------
    worst, _, _, _ = _worst_overlap(centres, radii_eff, Lx, Ly)
    tol = args.tol_frac * float(radii_eff.min())
    overlap_ok = worst <= tol

    ny = int(round(args.raster * Ly / Lx))
    solid = rasterize(centres, radii, Lx, Ly, args.raster, ny)
    phi_raster = 1.0 - solid.mean()
    phi_exact = 1.0 - float(np.pi * np.sum(radii ** 2)) / (Lx * Ly)

    sx, sy, s_frac, s_n = percolates(solid)
    px, py, p_frac, p_n = percolates(~solid)

    coord, n_contact, throats, _ = descriptors(centres, radii, Lx, Ly, gap)

    problems = []
    if not overlap_ok:
        problems.append(f"overlaps remain (worst {worst:.3e} m)")
    if abs(phi_exact - args.porosity) > 0.01:
        problems.append(f"porosity {phi_exact:.4f} misses target "
                        f"{args.porosity:.4f} by more than 0.01 "
                        f"(jammed at s={s_ok:.4f})")
    # Pore non-percolation is a FAILURE only when a positive gap was requested,
    # because then reconnecting the pore was the whole point of paying for it.
    # At gap <= 0 the grains are meant to be in contact, and in 2D a spanning
    # contact network necessarily cuts the void into cells -- so this is an
    # accepted, recorded consequence of the geometry, not a defect. Keeping it
    # in `problems` there would make generate_study_opts.py skip every packing.
    warnings = []
    if not (px and py):
        msg = (f"pore space does not percolate (x={px}, y={py}; largest cluster "
               f"holds {p_frac*100:.1f}% of the void in {p_n} clusters)")
        if gap > 0.0:
            problems.append(
                msg + " -- vapour cannot diffuse across the domain, and a "
                f"positive --contact-gap ({gap*1e6:.2f} um) was requested "
                "precisely to prevent this. Raise it further, or set it to <= 0 "
                "and accept contact (see the module docstring).")
        else:
            warnings.append(
                msg + f" -- EXPECTED at --contact-gap {gap*1e6:.2f} um. In 2D a "
                "spanning contact network cuts the void into cells; this is the "
                "accepted cost of having the grains actually touch. Rank seeds "
                "on pore retention with preprocess/select_packing_seeds.py.")

    meta = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "generator": "preprocess/generate_packing.py",
        "algorithm": "periodic collective rearrangement (Jodrey-Tory family), "
                     "radii grown to jamming",
        "periodic": bool(PERIODIC[0]),
        "Lx": Lx, "Ly": Ly,
        "seed": args.seed,
        "n_grains": int(n),
        "porosity_target": args.porosity,
        "porosity_achieved": phi_exact,
        "porosity_raster": float(phi_raster),
        "raster_nx": args.raster, "raster_ny": ny,
        "scale_factor_reached": float(s_ok),
        "contact_gap_m": float(gap),
        "mean_r_m_requested": args.mean_r_m,
        "mean_r_m_realized": float(radii.mean()),
        "min_r_m": float(radii.min()), "max_r_m": float(radii.max()),
        "sigma_ln": args.sigma_ln,
        "radius_clip_frac": args.radius_clip_frac,
        "max_overlap_m": float(worst),
        "overlap_free": bool(overlap_ok),
        "solid_percolates_x": sx, "solid_percolates_y": sy,
        "solid_largest_cluster_frac": s_frac, "solid_n_clusters": s_n,
        "pore_percolates_x": px, "pore_percolates_y": py,
        "pore_largest_cluster_frac": p_frac, "pore_n_clusters": p_n,
        "coordination_number": float(coord),
        "n_contacts": int(n_contact),
        "throat_gap_m": throats,
        "problems": problems,
        "warnings": warnings,
    }

    write_outputs(args.out, centres, radii, Lx, Ly, meta,
                  preview=not args.no_preview)

    print(f"  N={n}  phi={phi_exact:.4f} (raster {phi_raster:.4f})  "
          f"s={s_ok:.4f}  Z={coord:.2f}", file=sys.stderr)
    print(f"  gap={gap*1e6:.2f} um | pore perc x={px} y={py} "
          f"(largest {p_frac*100:.1f}%, {p_n} clusters) | solid perc "
          f"x={sx} y={sy} (largest {s_frac*100:.1f}%)", file=sys.stderr)
    for w in warnings:
        print(f"  note: {w}", file=sys.stderr)
    if problems:
        for p in problems:
            print(f"  PROBLEM: {p}", file=sys.stderr)
        return 1
    print("  OK", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
