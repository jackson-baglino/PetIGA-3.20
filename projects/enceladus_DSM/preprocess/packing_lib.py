#!/usr/bin/env python3
"""packing_lib.py — shared machinery for building and grading 2D disk packings.

Everything here except the periodicity plumbing is lifted unchanged from
preprocess/scratch/generate_packing.py. That script's *algorithm* (collective
rearrangement to jamming) was retired in favour of gravity-style deposition,
but its measurement and file I/O were not: how you grade a packing and how you
write it out are independent of how you built it, and these had already been
debugged against real k_eff runs.

WHAT LIVES HERE

    sample_radii    log-normal radii, count set by the solid-area target
    deposit-side helpers:
    relax           minimum-image overlap removal (used to heal a periodic seam)
    grading:
    rasterize       solid mask
    percolates      does a phase wrap around the torus
    uniformity      local_density_cv and max_void_radius -- the homogeneity pair
    descriptors     coordination number, contact count, throat-gap percentiles
    output:
    write_outputs   grains.dat + metadata.json + preview.png

PERIODICITY IS PER-AXIS AND EXPLICIT

The original carried a module-level ``PERIODIC = [True]`` cell that every
function read. Gravity deposition needs x and y to differ -- the strip wraps
sideways while gravity gives the vertical direction a floor and a free surface
-- which a single flag cannot express, and a hidden global that silently
changes what `rasterize` means is a bad way to carry it besides. Every function
that cares now takes ``px`` and ``py`` booleans.

That distinction is not cosmetic, and the original comment is worth repeating:
a window cut from a larger pack keeps every grain overlapping it, so centres
span [-R, L+R]. Rasterising those through the periodic branch wraps each
overhanging grain onto the OPPOSITE face -- adding solid where there is none
and removing it from where it is. Porosity, percolation and coordination are
all corrupted by getting this wrong.

THE GRAINS.DAT CONTRACT is fixed by ReadGrainsFromFile in
src/initial_conditions.c: an optional two-field ``Lx Ly`` header, then ``x y r``
rows in metres, with the row width fixed by the first data row. Do not change
it here without changing the reader.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import sys

import numpy as np
from scipy import ndimage
from scipy.spatial import Delaunay, cKDTree


# =========================================================================
# Radii
# =========================================================================

def sample_radii(rng, Lx, Ly, porosity, mean_r, sigma_ln, clip_frac=1.0):
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
# Minimum image
# =========================================================================

def min_image(d, Lx, Ly, px, py):
    """Fold separation vectors into the minimum image, per axis.

    `d` is modified in place and returned.
    """
    if px:
        d[:, 0] -= Lx * np.round(d[:, 0] / Lx)
    if py:
        d[:, 1] -= Ly * np.round(d[:, 1] / Ly)
    return d


def _offsets(Lx, Ly, px, py):
    """Image offsets to replicate by: 3 per periodic axis, 1 per open axis."""
    xs = (-Lx, 0.0, Lx) if px else (0.0,)
    ys = (-Ly, 0.0, Ly) if py else (0.0,)
    return np.array([[ox, oy] for ox in xs for oy in ys])


def _candidate_pairs(centres, radii, Lx, Ly, px, py, reach):
    """Unique (i, j) base-index pairs whose centres are within `reach`.

    Built by querying the base points against a replicated point set, so pairs
    that are only close across a periodic face are found. Returns an (M, 2)
    int array; self-image pairs (a grain wider than half the box meeting its
    own image) are dropped rather than silently double-counted.
    """
    n = len(centres)
    offs = _offsets(Lx, Ly, px, py)
    rep = (centres[None, :, :] + offs[:, None, :]).reshape(-1, 2)
    rep_id = np.tile(np.arange(n), len(offs))

    near = cKDTree(centres).query_ball_tree(cKDTree(rep), reach)
    seen = set()
    for i, lst in enumerate(near):
        for j in lst:
            k = int(rep_id[j])
            if k != i:
                seen.add((i, k) if i < k else (k, i))
    if not seen:
        return np.empty((0, 2), dtype=np.int64)
    return np.array(sorted(seen), dtype=np.int64)


# =========================================================================
# Overlap relaxation
# =========================================================================

def worst_overlap(centres, radii, Lx, Ly, px, py):
    """(max_overlap, pairs, deltas, dists) over all overlapping pairs.

    Overlap is (r_i + r_j) - d under the minimum image on the periodic axes.
    """
    if len(centres) < 2:
        return 0.0, np.empty((0, 2), dtype=np.int64), None, None

    pairs = _candidate_pairs(centres, radii, Lx, Ly, px, py,
                             2.0 * float(radii.max()))
    if pairs.size == 0:
        return 0.0, pairs, None, None

    i, j = pairs[:, 0], pairs[:, 1]
    d = min_image(centres[j] - centres[i], Lx, Ly, px, py)
    dist = np.hypot(d[:, 0], d[:, 1])
    overlap = (radii[i] + radii[j]) - dist

    keep = overlap > 0.0
    if not np.any(keep):
        return 0.0, pairs[:0], None, None
    return float(overlap[keep].max()), pairs[keep], d[keep], dist[keep]


def relax(centres, radii, Lx, Ly, px, py, tol_frac=1e-4, max_iter=4000,
          damping=0.5):
    """Push overlapping pairs apart until no overlap exceeds tol.

    Each overlapping pair is separated along its centre line, splitting the
    correction between the two disks. Displacements from all pairs touching a
    given disk are accumulated before moving anything, so the sweep is
    Jacobi-style (order-independent) rather than Gauss-Seidel.

    Centres are wrapped on periodic axes and clipped on open ones -- on an open
    axis a grain is kept WHOLLY inside, since the multi_grains IC only applies
    minimum-image wrapping under -periodic 1 and would otherwise cut a
    boundary-straddling grain, leaving a flat artificial face from t = 0.

    Returns (converged, iterations, worst_remaining_overlap).
    """
    tol = tol_frac * float(radii.min())
    for it in range(max_iter):
        worst, pairs, delta, dist = worst_overlap(centres, radii, Lx, Ly, px, py)
        if worst <= tol:
            return True, it, worst

        i, j = pairs[:, 0], pairs[:, 1]
        overlap = (radii[i] + radii[j]) - dist
        # Unit vector i->j; if two centres coincide, pick a random direction.
        with np.errstate(invalid="ignore", divide="ignore"):
            u = delta / dist[:, None]
        bad = ~np.isfinite(u).all(axis=1)
        if np.any(bad):
            ang = np.random.default_rng(it).uniform(0, 2 * np.pi, int(bad.sum()))
            u[bad] = np.stack([np.cos(ang), np.sin(ang)], axis=1)

        shift = damping * overlap[:, None] * u
        disp = np.zeros_like(centres)
        np.add.at(disp, i, -0.5 * shift)
        np.add.at(disp, j, +0.5 * shift)
        centres += disp

        if px:
            centres[:, 0] = np.mod(centres[:, 0], Lx)
        else:
            centres[:, 0] = np.clip(centres[:, 0], radii, Lx - radii)
        if py:
            centres[:, 1] = np.mod(centres[:, 1], Ly)
        else:
            centres[:, 1] = np.clip(centres[:, 1], radii, Ly - radii)

    worst, _, _, _ = worst_overlap(centres, radii, Lx, Ly, px, py)
    return worst <= tol, max_iter, worst


# =========================================================================
# Rasterization and percolation
# =========================================================================

def rasterize(centres, radii, Lx, Ly, nx, ny, px, py):
    """Boolean solid mask on an nx-by-ny grid.

    On an OPEN axis the mask is a plain box: a grain whose centre lies outside
    [0, L] contributes only the part that falls inside, and nothing wraps. See
    the module docstring for why feeding window-cut centres through the
    periodic branch corrupts every derived quantity.
    """
    dx, dy = Lx / nx, Ly / ny
    solid = np.zeros((ny, nx), dtype=bool)
    xs = (np.arange(nx) + 0.5) * dx
    ys = (np.arange(ny) + 0.5) * dy

    for (cx, cy), r in zip(centres, radii):
        i0, i1 = int(np.floor((cx - r) / dx)) - 1, int(np.ceil((cx + r) / dx)) + 1
        j0, j1 = int(np.floor((cy - r) / dy)) - 1, int(np.ceil((cy + r) / dy)) + 1

        if px:
            ii = np.arange(i0, i1 + 1)
            ddx = xs[np.mod(ii, nx)] - cx
            ddx -= Lx * np.round(ddx / Lx)
            ci = np.mod(ii, nx)
        else:
            ii = np.arange(max(i0, 0), min(i1, nx - 1) + 1)
            if ii.size == 0:
                continue
            ddx = xs[ii] - cx
            ci = ii

        if py:
            jj = np.arange(j0, j1 + 1)
            ddy = ys[np.mod(jj, ny)] - cy
            ddy -= Ly * np.round(ddy / Ly)
            cj = np.mod(jj, ny)
        else:
            jj = np.arange(max(j0, 0), min(j1, ny - 1) + 1)
            if jj.size == 0:
                continue
            ddy = ys[jj] - cy
            cj = jj

        inside = (ddx[None, :] ** 2 + ddy[:, None] ** 2) <= r * r
        if inside.any():
            solid[np.ix_(cj, ci)] |= inside
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


def uniformity(solid, Lx, Ly, ncoarse=8):
    """How evenly the solid is spread, and how big the worst void is.

    Returns (local_density_cv, max_void_radius_m).

    local_density_cv : coefficient of variation of the solid fraction over an
        ncoarse x ncoarse coarse grid. 0 would be perfectly uniform; clumping
        raises it because some cells fill up while others empty out.
    max_void_radius_m : radius of the largest circle that fits entirely in the
        pore space, from a Euclidean distance transform of the void. This is the
        direct measure of "overly large voids" -- a packing can have an
        acceptable mean porosity and still be unusable because all of it is in
        one hole.

    These two ARE the homogeneity criterion. Do not substitute pore-cluster
    connectivity for them; see the note in accept_reasons().
    """
    ny, nx = solid.shape
    ky, kx = max(1, ny // ncoarse), max(1, nx // ncoarse)
    cells = []
    for j in range(0, ny - ky + 1, ky):
        for i in range(0, nx - kx + 1, kx):
            cells.append(solid[j:j + ky, i:i + kx].mean())
    cells = np.asarray(cells, dtype=float)
    cv = float(cells.std() / cells.mean()) if cells.mean() > 0 else float("inf")

    dpix = ndimage.distance_transform_edt(~solid)
    return cv, float(dpix.max() * (Lx / nx))


def void_map(solid, Lx, Ly):
    """(distance_to_solid [m], dx, dy) for the pore space.

    The filler pass uses this to find the biggest hole: the global maximum is
    the largest inscribed empty circle, and its location is where a grain
    should go.
    """
    ny, nx = solid.shape
    dpix = ndimage.distance_transform_edt(~solid)
    return dpix * (Lx / nx), Lx / nx, Ly / ny


# =========================================================================
# Structural descriptors
# =========================================================================

def descriptors(centres, radii, Lx, Ly, px, py, gap=0.0, contact_tol_frac=0.02):
    """Coordination number, contact count, and throat (gap) statistics.

    Contacts: pairs whose surface gap is within contact_tol_frac * r_min of
    `gap`. Throats: gaps along Delaunay edges of the replicated centres,
    restricted to edges with at least one endpoint in the base tile.
    """
    tol = contact_tol_frac * float(radii.min())
    n = len(centres)
    if n < 2:
        return 0.0, 0, {}

    pairs = _candidate_pairs(centres, radii, Lx, Ly, px, py,
                             2.0 * float(radii.max()) + gap + tol)
    n_contact = 0
    if pairs.size:
        i, j = pairs[:, 0], pairs[:, 1]
        d = min_image(centres[j] - centres[i], Lx, Ly, px, py)
        sep = np.hypot(d[:, 0], d[:, 1]) - (radii[i] + radii[j])
        n_contact = int(np.sum(np.abs(sep - gap) <= tol))
    coord = 2.0 * n_contact / n

    offs = _offsets(Lx, Ly, px, py)
    rep = (centres[None, :, :] + offs[:, None, :]).reshape(-1, 2)
    rep_r = np.tile(radii, len(offs))
    base = np.zeros(len(rep), dtype=bool)
    zero = int(np.argmin(np.abs(offs).sum(axis=1)))       # the (0,0) block
    base[zero * n:(zero + 1) * n] = True

    stats = {}
    try:
        tri = Delaunay(rep)
    except Exception:
        return coord, n_contact, stats

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
        if gaps.size:
            for p in (5, 25, 50, 75, 95):
                stats[f"p{p}"] = float(np.percentile(gaps, p))
            stats["mean"] = float(gaps.mean())
            stats["n"] = int(gaps.size)
    return coord, n_contact, stats


# =========================================================================
# Acceptance
# =========================================================================

def accept_reasons(meta, max_void_ratio, max_density_cv, require_percolation=True):
    """List of reasons this packing should be REJECTED. Empty means accept.

    RANK ON UNIFORMITY, NEVER ON PORE CONNECTIVITY. This is carried over from
    the retired select_packing_seeds.py, which learned it the hard way, and it
    is easy to get backwards:

    Maximising the fraction of void held in one connected cluster selects for
    exactly the configurations that concentrate all the empty space in a single
    hole -- so the winners are clumped grains beside a large void, the opposite
    of a representative microstructure. One packing was picked under that rule
    BECAUSE it had a void running the full height of the domain (pore 69%, the
    highest in its group).

    So the gate is:
      max_void_radius_per_mean_r  the largest circle fitting in the pore, in
                                  mean grain radii -- "one overly large void"
      local_density_cv            clumping over an 8x8 grid
      solid percolation           the matrix must actually connect

    Pore connectivity is still measured and reported. It is not optimised for,
    and in 2D it cannot be had together with contacting grains anyway.
    """
    bad = []
    if meta["max_void_radius_per_mean_r"] > max_void_ratio:
        bad.append(f"largest void {meta['max_void_radius_per_mean_r']:.2f} "
                   f"mean-radii > {max_void_ratio:.2f}")
    if meta["local_density_cv"] > max_density_cv:
        bad.append(f"local_density_cv {meta['local_density_cv']:.3f} "
                   f"> {max_density_cv:.3f}")
    if require_percolation and not (meta["solid_percolates_x"]
                                    and meta["solid_percolates_y"]):
        bad.append(f"solid does not percolate "
                   f"(x={meta['solid_percolates_x']}, "
                   f"y={meta['solid_percolates_y']})")
    return bad


# =========================================================================
# Output
# =========================================================================

def write_outputs(out, centres, radii, Lx, Ly, meta, px, py, preview=True):
    """Write grains.dat, metadata.json and preview.png into `out`.

    EXPLICIT EDGE IMAGES. On a periodic axis a grain near the wall genuinely
    continues on the opposite face. The multi_grains IC only reconstructs that
    under -periodic 1; under -periodic 0 (zero-flux walls, which is what the
    phase-field, temperature and vapour equations want) it would simply CUT the
    grain and leave a flat artificial face.

    So the wrapping is baked into the file instead of left to the solver: any
    grain reaching past a periodic edge is emitted a second time at its image,
    giving centres that span [-R, L+R]. A grain then needs only to be PARTLY
    inside the domain, which is the point -- these packings represent a window
    cut from a much larger pack, and deleting the grains that straddle the edge
    would sever the conduction paths leaving the domain and bias k_eff low.
    Connectivity across the boundary is what k_eff is most sensitive to.

    Because the images are the true periodic partners, the microstructure also
    stays periodic across those faces, removing the seam mismatch the
    corrector's periodic BCs would otherwise see.
    """
    os.makedirs(out, exist_ok=True)

    xoffs = (-Lx, 0.0, Lx) if px else (0.0,)
    yoffs = (-Ly, 0.0, Ly) if py else (0.0,)
    emitted = []
    for (x, y), r in zip(centres, radii):
        for ox in xoffs:
            for oy in yoffs:
                cx, cy = x + ox, y + oy
                if (-r <= cx <= Lx + r) and (-r <= cy <= Ly + r):
                    emitted.append((cx, cy, r))

    meta["n_grain_rows_emitted"] = len(emitted)
    meta["n_edge_images"] = len(emitted) - len(radii)
    meta["periodic_x"] = bool(px)
    meta["periodic_y"] = bool(py)

    axes = "".join(a for a, f in (("x", px), ("y", py)) if f) or "none"
    dat = os.path.join(out, "grains.dat")
    with open(dat, "w") as f:
        f.write(f"# 2D packing, generated {meta['generated_utc']}\n")
        f.write(f"# generator = {meta.get('generator', 'packing_lib')}\n")
        f.write(f"# periodic axes = {axes}\n")
        f.write(f"# porosity_achieved = {meta['porosity_achieved']:.6f}"
                f"  (target {meta['porosity_target']:.6f})\n")
        f.write(f"# {len(radii)} distinct grains; {len(emitted)} rows including "
                f"{len(emitted) - len(radii)} periodic edge images.\n")
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
        _preview(out, emitted, radii, Lx, Ly, meta)
    return dat


def _preview(out, emitted, radii, Lx, Ly, meta):
    """Preview PNG. Optional -- a failure here must not lose the packing."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.patches import Circle

        fig, ax = plt.subplots(figsize=(6.4, 6.4 * Ly / Lx))
        # Draw the margin OUTSIDE the domain too, so grains hanging past an
        # edge are visible -- they are in grains.dat, the solver sees them, and
        # they carry the conduction paths that leave the domain. A preview
        # cropped to [0,L] hides exactly what k_eff is most sensitive to.
        pad = float(np.max(radii))
        n_filler = meta.get("n_filler_grains", 0)
        n_matrix = len(radii) - n_filler
        for k, (x, y, r) in enumerate(emitted):
            img = not (0 <= x <= Lx and 0 <= y <= Ly)
            filler = k < len(radii) and k >= n_matrix
            face = "#c9d9e6" if img else ("#e8b84b" if filler else "#7fb3d5")
            ax.add_patch(Circle((x, y), r, facecolor=face, edgecolor="#1b4f72",
                                lw=0.3, alpha=0.55 if img else 1.0,
                                zorder=2 if img else 3))
        ax.add_patch(plt.Rectangle((0, 0), Lx, Ly, fill=False,
                                   edgecolor="#c0392b", lw=1.8, zorder=4))
        ax.set_xlim(-pad, Lx + pad); ax.set_ylim(-pad, Ly + pad)
        ax.set_aspect("equal")
        ax.set_xlabel("x [m]", fontsize=8); ax.set_ylabel("y [m]", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.set_title(
            f"phi={meta['porosity_achieved']:.4f}  N={len(radii)}"
            f"  ({n_matrix} matrix + {n_filler} filler)  seed={meta['seed']}"
            f"  L={Lx * 1e3:.2f} mm\n"
            f"{len(emitted) - len(radii)} edge images (pale)   "
            f"max void {meta['max_void_radius_per_mean_r']:.2f} r_mean   "
            f"cv {meta['local_density_cv']:.3f}   "
            f"Z {meta['coordination_number']:.2f}",
            fontsize=8.5)
        fig.savefig(os.path.join(out, "preview.png"), dpi=140, bbox_inches="tight")
        plt.close(fig)
    except Exception as e:
        print(f"  (preview skipped: {e})", file=sys.stderr)
