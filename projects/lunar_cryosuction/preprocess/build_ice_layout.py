#!/usr/bin/env python3
"""build_ice_layout.py: a rough-walled pore (plan type 7) and a mixed initial
ice layout, written as a mesh, a geometry .opts, a node ice-field file for
-ic_type phi_file, and a preview PDF.

THE PORE. "walls": {"flat": true} gives straight walls (types 1, 2, 5).
Otherwise a channel of length Lx whose walls are chains of overlapping
C-infinity bumps (grain-like roughness), optionally with throats (a large
bump on both walls; "throats" takes a list, each sized after the ones
before it). Bumps are placed by a left-to-right walk so they always overlap
(no flat stretches). The baseline is set so the mean gap is --mean-gap (spec
"mean_gap"). Each wall is the sum of its bumps:

    y_bot(x) = sum_k b_k(x),   y_top(x) = H - sum_k t_k(x),
    b(x) = h exp(1 - 1/(1 - s^2)), s = (x - c)/R, |s| < 1.

"walls": {"profile": "cos" | "sinhyp", "W": ..., "wt": ...} gives a smooth
converging-diverging cell (plan type 5): gap g(x) about a throat at x_t
("x_throat", default Lx/2), walls symmetric about y = W/2 (y_bot = (W - g)/2,
y_top = (W + g)/2), narrowest (g = wt) at the throat and widest (g = W) half a
period away (at x = 0 and x = Lx for the default x_t). With xi = x - x_t,

    cos     g = wt + (W - wt) (1 - cos(2 pi xi / Lx)) / 2
    sinhyp  g = sqrt(wt^2 + (c Lx/pi)^2 sin^2(pi xi / Lx)),  c = pi sqrt(W^2 - wt^2) / Lx

sinhyp is the hyperbolic throat g = sqrt(wt^2 + c^2 xi^2) with xi replaced by
(Lx/pi) sin(pi xi/Lx): equal to it to O(xi^2/Lx^2) near the throat, but
Lx-periodic and C-infinity with zero slope at the seam. Both are periodic, so
they serve "periodic": true (one cell of an infinite chain of pores).

The mesh maps the unit square onto the pore: control points at the Greville
abscissae, x = Lx g_x, y = y_bot(x) + g_y (y_top(x) - y_bot(x)).

PERIODIC ("periodic": true, walls flat or a profile). The x knot vector is the
periodic uniform one PetIGA builds itself (petigaaxis.c, IGAAxisInitUniform):
U = (-p, ..., N + p)/N, N + p control columns, of which the last p repeat the
first p shifted by Lx (x + Lx, same y), so the mapping, the walls and their
slope are continuous across the seam. The solver (-periodic_x 1) wraps the
last p basis functions onto the first p: it has N x nodes, and the phi file
holds those N columns only. The ice is drawn on a raster three periods long
([-Lx, 2Lx), an integer number of cells per period) and folded: a feature may
cross the seam (a plug with x0 < 0 or x1 > Lx). Rounding and distances are
taken on the periodic (tiled) field.

THE ICE (a JSON layout file; see preprocess/layouts/). Features, each drawn
as a sharp region and then combined:

  plugs     ice lenses that span the pore. Built by the equal-curvature
            (morphological) rule: gas is the union of discs of radius r whose
            centers keep a distance r cos(theta) from the walls (a disc whose
            center is r cos(theta) from a flat wall meets it at theta through
            the ice), allowed only outside the plug, x <= x0 - r and
            x >= x1 + r; the ice is the pore between, in [x0 - r, x1 + r].
            Faces are therefore arcs of radius r meeting the rough walls at
            theta. r = g/(2 cos theta), with g the mean gap over the plug
            ("r" in the spec overrides). theta >= 90 is not handled.
            A plug with x1 >= Lx ends at the dead end.
  caps      droplets on a wall: a disc of radius R whose center sits R cos(theta)
            behind the local wall along its normal, cut by the pore.
  fill      (optional, at most one) equal-curvature fill of the whole pore: one
            face radius r for every body, so ice spans each stretch whose gap
            is below 2 r cos(theta), as plugs (above) with that r.
            {"r": ...} or {"ice_frac": ...} (r found by bisection), with
            optional "theta", "x0", "x1" (ice kept only in [x0, x1]).
  grains    free ice discs (or ellipses ax, ay) at (x, y); y may be given as
            "axis" (the local mid-gap).

theta for plugs and caps is the layout's "theta_ice_init" (the user's rule,
2026-11-09: droplets start at liquid water's angle, lenses at the wall-BC
angle; a feature may override with its own "theta").

ROUNDING. The sharp union is then opened and closed with a disc of radius
"round_eps" x eps (default 3): no ice-air corner sharper than that radius and
no spike or notch narrower than twice it survives. A first type-7 test without
this step (2026-11-09) drove phi past the 1.01 band within 4 s at corners
where lens faces meet rough walls, and the run aborted.

THE FIELD. phi = 1/2 - 1/2 tanh(d / (2 eps)) at each control point, with d
the signed distance to the sharp ice (positive in air), computed on a raster
of spacing eps/4. Walls are not interfaces: distances are measured to ice and
to air only. This is the profile every shape IC in src/initial_conditions.c
uses (PhaseSlope = 1/(2 eps)).

Checks printed: ice fraction, number of ice bodies, the thinnest body (max
half-thickness) against 3 eps, and the narrowest air gap between two bodies
against 6 eps (closer bodies merge as soon as the run starts).

    venv_cryosuction/bin/python preprocess/build_ice_layout.py preprocess/layouts/type7_mixed_A.json
        [--eps 8.584e-7] [--name NAME] [--no-pdf]

Writes inputs/geometry/meshes/<name>.dat, inputs/geometry/meshes/<name>_phi.bin,
inputs/geometry/cryosuction/<name>.opts and preprocess/<name>.pdf.
"""
from __future__ import annotations

import argparse
import json
import math
import struct
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "preprocess"))
sys.path.insert(0, str(ROOT / "postprocess"))
from build_geometry_multi_grain import _bump, greville_abscissae, open_uniform_knots  # noqa: E402

P, C = 2, 1
EPS = 8.584e-7
T0_C = -20.0
MAGIC = b"PHIIC001"


# ---------------------------------------------------------------- walls ----
def bump_walk(rng, lx, n, R_range, h_range, ov_range, span_tol=0.03, tries=2000):
    """n overlapping bumps covering [0, Lx] (centers at 0 and Lx)."""
    for _ in range(tries):
        R = rng.uniform(*R_range, size=n)
        h = R * rng.uniform(*h_range, size=n)
        step = (R[:-1] + R[1:]) * (1.0 - rng.uniform(*ov_range, size=n - 1))
        if abs(step.sum() - lx) <= span_tol * lx:
            c = np.concatenate([[0.0], np.cumsum(step * lx / step.sum())])
            return np.column_stack([c, R, h])
    raise RuntimeError("bump walk did not fit Lx; change n or the ranges")


def bump_sum(x, bumps):
    y = np.zeros_like(np.asarray(x, float))
    for c, R, h in bumps:
        y += _bump(np.asarray(x, float), c, R, h)
    return y


def gap_profile(x, lx, W, wt, kind, xt=None):
    """gap of the converging-diverging cell: wt at the throat xt (default Lx/2),
    W half a period away."""
    xi = np.asarray(x, float) - (0.5 * lx if xt is None else xt)
    if kind == "cos":
        return wt + (W - wt) * 0.5 * (1.0 - np.cos(2.0 * np.pi * xi / lx))
    if kind == "sinhyp":
        a2 = W * W - wt * wt                       # (c Lx / pi)^2
        return np.sqrt(wt * wt + a2 * np.sin(np.pi * xi / lx) ** 2)
    raise ValueError(f"unknown wall profile {kind!r} (cos or sinhyp)")


class Pore:
    def __init__(self, spec):
        self.lx = spec["Lx"]
        self.periodic = bool(spec.get("periodic", False))
        w = spec["walls"]
        self.profile = w.get("profile")
        if self.profile:                           # converging-diverging cell (type 5)
            self.W, self.wt = float(w["W"]), float(w["wt"])
            self.xt = float(w.get("x_throat", 0.5 * self.lx))
            assert 0 < self.wt < self.W, "profile needs 0 < wt < W"
            gap_profile(0.0, self.lx, self.W, self.wt, self.profile)   # name check
            self.bot = self.top = np.zeros((0, 3))
            self.H = self.W
            self.throats, self.throat = [], None
            return
        if w.get("flat"):                          # straight walls (types 1, 2, 5)
            self.bot = self.top = np.zeros((0, 3))
        else:
            rb, rt = np.random.default_rng(w["seed_bottom"]), np.random.default_rng(w["seed_top"])
            args = (w["n_bumps"], w["R_range"], w["height_ratio_range"], w["overlap_range"])
            self.bot = bump_walk(rb, self.lx, *args)
            self.top = bump_walk(rt, self.lx, *args)
        xs = np.linspace(0, self.lx, 40001)
        self.H = spec["mean_gap"] + bump_sum(xs, self.bot).mean() + bump_sum(xs, self.top).mean()
        # "throat" (one) or "throats" (a list, added in order; each later throat is
        # sized against the walls with the earlier throats already in place)
        self.throats = []
        for th in spec.get("throats") or ([spec["throat"]] if spec.get("throat") else []):
            # throat height chosen so the NARROWEST gap within the throat bump is th["gap"]
            xw = np.linspace(th["x"] - th["R"], th["x"] + th["R"], 4001)
            yb0, yt0 = self.walls(xw)
            g0 = yt0 - yb0
            assert g0.min() > th["gap"], f"throat gap {th['gap']:.3e} is wider than the pore there"
            lo, hi = 0.0, 0.5 * g0.max()
            for _ in range(60):
                hth = 0.5 * (lo + hi)
                gmin = (g0 - 2 * _bump(xw, th["x"], th["R"], hth)).min()
                lo, hi = (hth, hi) if gmin > th["gap"] else (lo, hth)
            self.throats.append((th["x"], th["R"], 0.5 * (lo + hi)))
        self.throat = self.throats[0] if self.throats else None
        if self.periodic:
            assert not (len(self.bot) or self.throats), \
                "periodic pores take flat walls or a wall profile (bumps are not periodic)"

    def walls(self, x):
        x = np.asarray(x, float)
        if self.profile:
            g = gap_profile(x, self.lx, self.W, self.wt, self.profile, self.xt)
            return 0.5 * (self.W - g), 0.5 * (self.W + g)
        yb, yt = bump_sum(x, self.bot), self.H - bump_sum(x, self.top)
        for t_ in self.throats:
            yb = yb + _bump(x, *t_)
            yt = yt - _bump(x, *t_)
        return yb, yt


def periodic_uniform_knots(N, p):
    """Periodic uniform knots on [0, 1] (continuity p - 1): the vector PetIGA's
    IGAAxisInitUniform builds for a periodic axis, U[k] = (k - p)/N, k = 0..N + 2p.
    N + p basis functions; the last p are the first p shifted by one period."""
    return (np.arange(N + 2 * p + 1, dtype=float) - p) / N


def build_surface(pore, nx, ny):
    from igakit.nurbs import NURBS
    uy = open_uniform_knots(ny, P, C)
    gy = greville_abscissae(uy, P)
    if pore.periodic:
        assert C == P - 1
        ux = periodic_uniform_knots(nx, P)
        gx = greville_abscissae(ux, P)                  # nx + P values, first is -1/(2 nx)
        x = pore.lx * gx[:nx]
        yb, yt = pore.walls(x)
        # the last P columns are the first P shifted by exactly one period
        x = np.concatenate([x, x[:P] + pore.lx])
        yb, yt = np.concatenate([yb, yb[:P]]), np.concatenate([yt, yt[:P]])
    else:
        ux = open_uniform_knots(nx, P, C)
        gx = greville_abscissae(ux, P)
        x = pore.lx * gx
        yb, yt = pore.walls(x)
    ctrl = np.zeros((len(gx), len(gy), 3))
    ctrl[:, :, 0] = x[:, None]
    ctrl[:, :, 1] = yb[:, None] + gy[None, :] * (yt - yb)[:, None]
    return NURBS([ux, uy], control=ctrl)


# ------------------------------------------------------------------ ice ----
class Raster:
    """Pixel grid of the pore. Open pores: x in [0, Lx]. Periodic pores: three
    periods, x in [-Lx, 2Lx), with M = Lx/dx cells per period exactly (dx is
    rounded so), so a field is periodic when its three thirds are equal;
    middle() is the period [0, Lx) and periodize() re-tiles it."""

    def __init__(self, pore, dx):
        self.periodic = pore.periodic
        if self.periodic:
            self.M = int(round(pore.lx / dx))
            dx = pore.lx / self.M
            self.x = -pore.lx + dx * np.arange(3 * self.M)
        else:
            self.x = np.arange(0.0, pore.lx + 0.5 * dx, dx)
            self.M = len(self.x)
        self.dx = dx
        self.x0 = float(self.x[0])
        yb, yt = pore.walls(self.x)
        self.yb, self.yt = yb, yt
        self.y = np.arange(0.0, pore.H + 0.5 * dx, dx)
        self.X, self.Y = np.meshgrid(self.x, self.y, indexing="ij")
        self.inside = (self.Y >= yb[:, None]) & (self.Y <= yt[:, None])
        # distance to the walls: the pore ends are not walls (extend them by edge copies);
        # the solid beyond the raster's top and bottom rows is (pad with solid rows)
        pad = int(np.ceil(2e-4 / dx))
        m = np.pad(self.inside, ((pad, pad), (0, 0)), mode="edge")
        m = np.pad(m, ((0, 0), (1, 1)), constant_values=False)
        self.D = ndi.distance_transform_edt(m, sampling=dx)[pad:-pad, 1:-1]

    def idx(self, x):
        return int(round((x - self.x0) / self.dx))

    def middle(self, a):
        return a[self.M:2 * self.M] if self.periodic else a

    def periodize(self, a):
        return np.tile(self.middle(a), (3, 1)) if self.periodic else a

    def fold(self, a):
        """periodic image of a feature drawn on the 3-period raster: OR of the thirds."""
        if not self.periodic:
            return a
        M = self.M
        return np.tile(a[:M] | a[M:2 * M] | a[2 * M:], (3, 1))


def plug(r, rs: Raster, x0, x1, theta):
    """equal-curvature plug between x0 and x1 (sharp mask)."""
    ct = math.cos(math.radians(theta))
    assert ct > 0, "plugs need theta < 90"
    g = rs.yt - rs.yb
    m = (rs.x >= x0) & (rs.x <= min(x1, rs.x[-1]))
    rad = r if r else 0.5 * g[m].mean() / ct
    # A face at x0 exists only if a gas disc can be centered at (x0 - r, mid-gap) with its
    # center r cos(theta) from the walls (same for x1); shrink r (5 % steps) until both
    # faces form. Both faces keep the same r (one curvature, one pressure).
    def fits(xf, side):
        if xf >= rs.x[-1] - rs.dx:
            return True                          # the plug ends at the dead end
        xc = xf - side * rad
        if xc < rs.x0 or xc > rs.x[-1]:
            return True                          # face center beyond the pore end
        i = rs.idx(xc)
        return rs.D[i].max() >= rad * ct
    if not r:
        while not (fits(x0, +1) and fits(x1, -1)):
            rad *= 0.95
    centers = rs.inside & (rs.D >= rad * ct) & ((rs.X <= x0 - rad) | (rs.X >= x1 + rad))
    gas = ndi.distance_transform_edt(~centers, sampling=rs.dx) <= rad
    ice = rs.inside & ~gas & (rs.X >= x0 - rad) & (rs.X <= x1 + rad)
    return ice, rad


def fill(rs: Raster, r, theta, x0, x1):
    """equal-curvature fill (sharp mask): one ice-gas face radius r for the whole
    pore. In a channel a face of radius r meeting the walls at theta stands where
    the local gap is g_c = 2 r cos(theta) (flat-wall estimate; a sloped wall shifts
    it, theta -/+ the wall angle). Ice therefore spans every stretch with gap < g_c:
    each stretch becomes a plug with faces of radius r (plug()); stretches closer
    than about 2 r merge, since no gas body of that curvature fits between them.
    Stretches are looked for in [x0 + r, x1] so the faces stay at x >= x0.
    (A global disc opening was tried first, 2026-11-12: discs of radius r centered
    in a wide body reach through a short throat, leaving the throat gas.)"""
    ct = math.cos(math.radians(theta))
    assert ct > 0, "fill needs theta < 90"
    g = rs.yt - rs.yb
    narrow = (g < 2.0 * r * ct) & (rs.x >= x0 + r) & (rs.x <= x1)
    lab, n = ndi.label(narrow)
    ivs = [(rs.x[sl[0]][0], rs.x[sl[0]][-1]) for sl in ndi.find_objects(lab)]
    merged = []
    for a, b in ivs:
        if merged and a - merged[-1][1] < 2.0 * r:
            merged[-1] = (merged[-1][0], b)
        else:
            merged.append((a, b))
    ice = np.zeros_like(rs.inside)
    for a, b in merged:
        ice |= plug(r, rs, a, b, theta)[0]
    return ice


def fill_for_fraction(rs: Raster, target, theta, x0, x1):
    """r that gives the ice fraction target (of the whole pore) by bisection on r;
    the ice fraction grows with r."""
    g = rs.yt - rs.yb
    ct = math.cos(math.radians(theta))
    lo, hi = 0.5 * float(g.min()) / ct, 0.5 * float(g.max()) / ct
    area = rs.inside.sum()
    for _ in range(30):
        r = 0.5 * (lo + hi)
        f = fill(rs, r, theta, x0, x1).sum() / area
        lo, hi = (r, hi) if f < target else (lo, r)
    return 0.5 * (lo + hi)


def wall_normal_point(rs: Raster, xc, wall):
    """wall point at xc and the unit normal pointing into the pore."""
    i = rs.idx(xc)
    yw = (rs.yb if wall == "bottom" else rs.yt)
    s = np.gradient(yw, rs.dx)[i]
    n = np.array([-s, 1.0]) / math.hypot(s, 1.0)
    if wall == "top":
        n = -n
    return np.array([xc, yw[i]]), n


def cap(rs: Raster, xc, wall, R, theta):
    p, n = wall_normal_point(rs, xc, wall)
    c = p - n * R * math.cos(math.radians(theta))
    return rs.inside & (np.hypot(rs.X - c[0], rs.Y - c[1]) <= R)


def grain(rs: Raster, xc, yc, ax, ay):
    if yc == "axis":
        i = rs.idx(xc)
        yc = 0.5 * (rs.yb[i] + rs.yt[i])
    return rs.inside & (((rs.X - xc) / ax) ** 2 + ((rs.Y - yc) / ay) ** 2 <= 1.0), yc


def _erode(a, dx, rho):
    return ndi.distance_transform_edt(a, sampling=dx) > rho


def _dilate(a, dx, rho):
    return ndi.distance_transform_edt(~a, sampling=dx) <= rho


def round_features(ice, inside, dx, rho):
    """Opening then closing of (ice OR solid) by a disc of radius rho, cut back to
    the pore. Removes ice spikes and air notches narrower than 2 rho and rounds
    every ice-air corner to a radius of at least rho. The solid counts as ice, so
    the ice-wall contacts are not eroded (only ice-air corners are rounded)."""
    a = ice | ~inside
    a = _dilate(_erode(a, dx, rho), dx, rho)      # opening: rounds convex ice corners
    a = _erode(_dilate(a, dx, rho), dx, rho)      # closing: rounds concave ice corners
    return a & inside


def label_bodies(ice, rs=None):
    """ndi.label, with bodies that cross the seam of a periodic raster merged.
    Returns labels on the raster's layout (tiled for a periodic one) and n."""
    if rs is None or not rs.periodic:
        return ndi.label(ice)
    mid = rs.middle(ice)
    lab, n = ndi.label(mid)
    parent = list(range(n + 1))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for a, b in zip(lab[0], lab[-1]):               # x = 0 column against x = Lx - dx
        if a and b:
            parent[find(a)] = find(b)
    roots = sorted({find(k) for k in range(1, n + 1)})
    new = np.zeros(n + 1, int)
    for k in range(1, n + 1):
        new[k] = roots.index(find(k)) + 1
    return np.tile(new[lab], (3, 1)), len(roots)


def bodies_report(ice, inside, dx, eps, rs=None):
    lab, n = label_bodies(ice, rs)
    half = ndi.distance_transform_edt(ice, sampling=dx)
    mid = rs.middle if rs is not None else (lambda a: a)
    thick = np.array(ndi.maximum(mid(half), mid(lab), index=np.arange(1, n + 1))) if n else np.array([])
    # narrowest air gap between two different bodies (through the pore)
    min_gap = np.inf
    if n > 1:
        for k in range(1, n + 1):
            others = mid(ice & (lab != k))
            if not others.any():
                continue
            dk = mid(ndi.distance_transform_edt(lab != k, sampling=dx))
            min_gap = min(min_gap, float(dk[others].min()))
    return n, thick, min_gap


def node_field(surf, rs: Raster, ice, eps):
    """(x, y, phi) at every control point, x index fastest."""
    d_ice = ndi.distance_transform_edt(~ice, sampling=rs.dx)
    air = rs.inside & ~ice
    d_air = ndi.distance_transform_edt(~air, sampling=rs.dx)
    sd = np.where(ice, -d_air, d_ice)
    # Control points on a wall can land in a raster cell just outside the pore. Give
    # every outside cell the value of its nearest pore cell, so a wall point inside a
    # lens is not read as air (a first build did that: isolated phi = 0.6 dips on the
    # wall drove phi past the 1.01 band within 4 s, 2026-11-09).
    _, (ii, jj) = ndi.distance_transform_edt(~rs.inside, return_indices=True)
    sd = sd[ii, jj]
    cp = surf.control[:, :, :2]                       # (mx, my, 2)
    fi = (cp[..., 0] - rs.x0) / rs.dx
    fj = cp[..., 1] / rs.dx
    d = ndi.map_coordinates(sd, [fi.ravel(), fj.ravel()], order=1, mode="nearest").reshape(fi.shape)
    phi = 0.5 - 0.5 * np.tanh(d / (2.0 * eps))
    rec = np.stack([cp[..., 0], cp[..., 1], phi], axis=-1)   # (mx, my, 3)
    return rec.transpose(1, 0, 2)                              # (my, mx, 3): x fastest when raveled


def write_phi(path, rec):
    my, mx, _ = rec.shape
    with open(path, "wb") as f:
        f.write(MAGIC)
        f.write(struct.pack("<qq", mx, my))
        f.write(np.ascontiguousarray(rec, dtype="<f8").tobytes())


def read_phi(path):
    with open(path, "rb") as f:
        assert f.read(8) == MAGIC
        mx, my = struct.unpack("<qq", f.read(16))
        return np.frombuffer(f.read(), dtype="<f8").reshape(my, mx, 3)


def throat_txt(pore):
    if len(pore.throats) == 1:
        return ", throat at x = %.1f um" % (pore.throat[0] * 1e6)
    return "".join(", throat at x = %.1f um" % (t_[0] * 1e6) for t_ in pore.throats)


# ----------------------------------------------------------------- main ----
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("layout", help="layout JSON (preprocess/layouts/)")
    ap.add_argument("--eps", type=float, default=None, help="interface parameter [m] (default: layout, else 8.584e-7)")
    ap.add_argument("--name", default=None, help="output name (default: layout 'name')")
    ap.add_argument("--outdir", default=None, help="write all files here instead of the project folders")
    ap.add_argument("--no-pdf", action="store_true")
    a = ap.parse_args(argv)

    spec = json.loads(Path(a.layout).read_text())
    eps = a.eps or spec.get("eps", EPS)
    name = a.name or spec["name"]
    th0 = spec["theta_ice_init"]
    pore = Pore(spec)
    lx, H = pore.lx, pore.H
    xs = np.linspace(0, lx, 40001)
    yb, yt = pore.walls(xs)
    gap = yt - yb
    assert gap.min() > 10 * eps, f"pore pinches to {gap.min():.3e} m"
    nx = math.ceil(lx * math.sqrt(2) / eps)
    ny = math.ceil(gap.max() * math.sqrt(2) / eps)
    per = pore.periodic
    if per:                                       # value and slope across the seam
        hd = 1e-3 * eps
        y0, y1 = pore.walls(np.array([-hd, 0.0, hd]))[0], pore.walls(lx + np.array([-hd, 0.0, hd]))[0]
        print(f"  periodic in x: seam wall jump {abs(y1[1] - y0[1]):.2e} m, slope jump "
              f"{abs((y1[2] - y1[0]) - (y0[2] - y0[0])) / (2 * hd):.2e}")
    slope = max(np.abs(np.gradient(yb, xs)).max(), np.abs(np.gradient(yt, xs)).max())
    print(f"pore {name}: Lx {lx*1e6:.1f} um, H {H*1e6:.1f} um, gap {gap.min()*1e6:.1f}-{gap.max()*1e6:.1f} um "
          f"(mean {gap.mean()*1e6:.1f}), steepest wall slope {slope:.2f}, "
          f"{len(pore.bot)}+{len(pore.top)} wall bumps" + throat_txt(pore))
    nxd = nx if per else nx + P                   # x nodes the solver has
    print(f"  eps {eps:.4e} m, Nx {nx}, Ny {ny}, DOF {nxd*(ny+P)*3/1e6:.2f} M; element "
          f"{lx/nx*1e6:.3f} um long, {gap.max()/ny*1e6:.3f}-{gap.min()/ny*1e6:.3f} um across")

    rs = Raster(pore, eps / 4.0)
    ice = np.zeros_like(rs.inside)
    parts = []
    for f in spec.get("plugs", []):
        th = f.get("theta", th0)
        m, rad = plug(f.get("r"), rs, f["x0"], f["x1"], th)
        ice |= m
        parts.append(f"plug {f['x0']*1e6:.0f}-{f['x1']*1e6:.0f} um, theta {th:g}, faces r {rad*1e6:.1f} um")
    if spec.get("fill"):
        f = spec["fill"]
        th = f.get("theta", th0)
        fx0, fx1 = f.get("x0", 0.0), f.get("x1", lx)
        rad = f.get("r") or fill_for_fraction(rs, f["ice_frac"], th, fx0, fx1)
        ice |= fill(rs, rad, th, fx0, fx1)
        parts.append(f"fill {fx0*1e6:.0f}-{fx1*1e6:.0f} um, theta {th:g}, one face radius r {rad*1e6:.1f} um"
                     + (f" (set for {f['ice_frac']:.0%} ice before the other features)" if "r" not in f else ""))
    for f in spec.get("caps", []):
        th = f.get("theta", th0)
        ice |= cap(rs, f["x"], f["wall"], f["R"], th)
        parts.append(f"cap x {f['x']*1e6:.0f} um on the {f['wall']} wall, R {f['R']*1e6:.0f} um, theta {th:g}")
    for f in spec.get("grains", []):
        ax_, ay_ = f.get("ax", f.get("R")), f.get("ay", f.get("R"))
        m, yc = grain(rs, f["x"], f["y"], ax_, ay_)
        ice |= m
        parts.append(f"grain ({f['x']*1e6:.0f}, {yc*1e6:.0f}) um, {ax_*1e6:.0f} x {ay_*1e6:.0f} um")
    for p_ in parts:
        print("  " + p_)
    ice = rs.fold(ice)                           # periodic: features wrap at the seam

    rho = spec.get("round_eps", 3.0) * eps
    before = ice.sum()
    ice = rs.periodize(round_features(ice, rs.inside, rs.dx, rho))
    print(f"  rounded ice-air corners to radius >= {rho*1e6:.2f} um ({rho/eps:g} eps): area change "
          f"{(ice.sum() - before) / max(before, 1):+.2%}")
    lab, nb = label_bodies(ice, rs)
    half = ndi.distance_transform_edt(ice, sampling=rs.dx)
    th_b = np.array(ndi.maximum(half, lab, index=np.arange(1, nb + 1)))
    drop = np.flatnonzero(th_b < 3 * eps) + 1
    if len(drop):
        area = [float((lab == k).sum() * rs.dx ** 2) for k in drop]
        ice &= ~np.isin(lab, drop)
        print(f"  dropped {len(drop)} slivers thinner than 3 eps (wall crevices next to plug faces), "
              f"total area {sum(area):.3e} m^2")
    frac = rs.middle(ice).sum() / rs.middle(rs.inside).sum()
    n, thick, min_gap = bodies_report(ice, rs.inside, rs.dx, eps, rs)
    nthin = int((thick < 3 * eps).sum())
    print(f"  ice {frac:.1%} of the pore area ({rs.middle(ice).sum()*rs.dx**2:.4e} m^2), {n} bodies; thinnest "
          f"half-thickness {thick.min()*1e6:.2f} um ({thick.min()/eps:.1f} eps); {nthin} below 3 eps")
    if n > 1:
        print(f"  narrowest air gap between bodies {min_gap*1e6:.2f} um ({min_gap/eps:.1f} eps)")
    ok = nthin == 0 and (n < 2 or min_gap >= 6 * eps)
    if not ok:
        print("  WARNING: a body is below 3 eps or two bodies are closer than 6 eps (see above)")

    outdir = Path(a.outdir) if a.outdir else None
    mesh_dir = outdir or ROOT / "inputs/geometry/meshes"
    opts_dir = outdir or ROOT / "inputs/geometry/cryosuction"
    pdf_dir = outdir or ROOT / "preprocess"
    for d_ in (mesh_dir, opts_dir, pdf_dir):
        d_.mkdir(parents=True, exist_ok=True)
    dat, phi_path = mesh_dir / f"{name}.dat", mesh_dir / f"{name}_phi.bin"
    opts = opts_dir / f"{name}.opts"

    surf = build_surface(pore, nx, ny)
    rec = node_field(surf, rs, ice, eps)
    if per:
        rec = rec[:, :nx]                         # the solver's N periodic x nodes
    from igakit.io import PetIGA
    PetIGA().write(str(dat), surf)
    write_phi(phi_path, rec)
    sharp_area = float(rs.middle(ice).sum() * rs.dx ** 2)

    def rel(p):
        try:
            return p.resolve().relative_to(ROOT)
        except ValueError:
            return p.resolve()
    opts.write_text(f"""# =============================================================================
# geometry {opts.name}   (preprocess/build_ice_layout.py {rel(Path(a.layout))})
#
# """ + (f"Converging-diverging cell, wall profile {pore.profile} (throat {pore.wt*1e6:.1f} um at x = "
           f"{pore.xt*1e6:.1f} um, widest {pore.W*1e6:.1f} um)" if pore.profile else "Rough-walled pore") + f""": Lx {lx*1e6:.1f} um, gap {gap.min()*1e6:.1f}-{gap.max()*1e6:.1f} um (mean {gap.mean()*1e6:.1f}),
# {len(pore.bot)}+{len(pore.top)} wall bumps{throat_txt(pore)}.
# Walls are in the mesh; wall faces {spec.get("wall_faces", "y0,y1")}; """ + (
        "periodic in x (-periodic_x 1): x = 0 and x = Lx are one seam."
        if per else "x0 and x1 are the pore ends.") + f"""
# Ice from {rel(phi_path)} (-ic_type phi_file), initial angle {th0:g} deg:
""" + "".join(f"#   {p_}\n" for p_ in parts) + f"""# ice {frac:.1%} of the pore, sharp area {sharp_area:.4e} m^2, {n} bodies.
# eps {eps:.4e} m, h = eps/sqrt(2): Nx {nx}, Ny {ny}.
# =============================================================================
# DOF_GRID: {nxd} {ny+P}
-geom_file {rel(dat)}
-p {P}
-C {C}
-ic_type phi_file
-phi_file {rel(phi_path)}
-dim 2
-Lx {lx:.6e}
-Ly {H:.6e}
-Lz 0
-delt_t 1.0e-4
-eps {eps:.4e}
-eps_valid_temp {T0_C:g}
-periodic 0
""" + ("-periodic_x 1\n" if per else "") + f"""-wall_faces {spec.get("wall_faces", "y0,y1")}
""" + (f"-dtmax {spec['dtmax']:.4e}                     # 2 tau_sub at this eps (layout 'dtmax')\n"
       if "dtmax" in spec else ""))
    print(f"wrote {dat}\nwrote {phi_path}  ({rec.shape[1]} x {rec.shape[0]} control points)\nwrote {opts}")
    if not a.no_pdf:
        render(pore, rs, rs.middle(ice), rec, name, pdf_dir / f"{name}.pdf", nx, ny, eps, frac, n)
    return 0 if ok else 2


def render(pore, rs, ice, rec, name, out, nx, ny, eps, frac, n):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from figcheck import check_figure
    s = 1e6
    fig, axs = plt.subplots(2, 1, figsize=(7.5, 3.6), constrained_layout=True)
    a, b = axs
    for ax in axs:
        ax.fill_between(rs.x * s, 0, rs.yb * s, color="0.6", lw=0)
        ax.fill_between(rs.x * s, rs.yt * s, pore.H * s, color="0.6", lw=0)
        ax.plot(rs.x * s, rs.yb * s, "k-", lw=0.5)
        ax.plot(rs.x * s, rs.yt * s, "k-", lw=0.5)
        ax.set_xlim(0, pore.lx * s)
        ax.set_ylim(0, pore.H * s)
        ax.set_aspect("equal")
        ax.set_ylabel("y [µm]")
    st = max(1, int(round(1e-6 / rs.dx)))
    X, Y = rs.middle(rs.X), rs.middle(rs.Y)
    a.contourf(X[::st, ::st] * s, Y[::st, ::st] * s, ice[::st, ::st].astype(float),
               levels=[0.5, 1.5], colors=["#66b3ff"])
    a.set_title(f"{name}: sharp ice, {frac:.0%} of the pore, {n} bodies", fontsize=8)
    b.contourf(rec[..., 0] * s, rec[..., 1] * s, rec[..., 2], levels=np.linspace(0, 1, 11), cmap="Blues")
    b.set_title(f"phi at the {rec.shape[1]} x {rec.shape[0]} control points (what the solver reads); "
                f"Nx {nx}, Ny {ny}, eps {eps*1e6:.3f} µm", fontsize=7.5)
    b.set_xlabel("x [µm]   (face x0 at x = 0, face x1 at x = Lx)")
    bad = check_figure(fig)
    fig.savefig(out)
    plt.close(fig)
    print(f"wrote {out}  figcheck: {bad or 'OK'}")


if __name__ == "__main__":
    raise SystemExit(main())
