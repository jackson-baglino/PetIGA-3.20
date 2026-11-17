#!/usr/bin/env python3
"""vapor_movie.py -- movie of the ice and the vapor saturation of a 2D run, optionally with tracer particles
that follow the vapor flux. Pure matplotlib + cmocean + scipy (no ParaView), one frame per vtkOut snapshot.

WHAT IT DRAWS
  * Air (IcePhase < 0.5, Gouraud shading): supersaturation sigma = rho_v / rho_vs^I(T) - 1 on cmocean 'balance', symmetric
    about 0, so white is saturated, blue is undersaturated (sublimating), red is supersaturated (depositing).
    rho_vs^I is pplib.rho_vs, which mirrors the solver's RhoVS_I, with -rho_air read from the staged opts.
  * Ice (IcePhase >= 0.5): filled contours of IcePhase on cmocean 'ice'; black line at IcePhase = 0.5.
  * Regolith: anything outside the mapped mesh, drawn as a flat tan background.
  * --tracers: gold particles carried along the model's vapor flux direction, -phi_a grad(rho_v) (the vapor
    diffusion term of the residual, src/assembly.c; its constant factor xi_v*D_v drops out).

THE TRACERS ARE A VISUALIZATION, NOT PART OF THE MODEL
  The solver has no particles. For each snapshot the flux is interpolated to a regular grid and particles are
  advected through it (RK2, --subframes steps per snapshot). Their SPEED IS RESCALED FOR VISIBILITY: a particle
  moves --dmax um per step at the reference flux (the median over snapshots of the 90th percentile of |flux| in
  the air), and no faster. With the physical D_v the vapor drift speed D_v |grad ln rho_v| is ~0.1-0.3 mm/s in the
  stage-0b runs, so vapor crosses a pore in seconds, far less than the hours between snapshots, so the absolute speed means nothing; direction and relative speed (fast in a throat, stalled behind
  ice) are what the particles show. Particles respawn at random air locations when they enter ice, leave the
  mesh, or reach the domain edge moving outward (an open vapor face such as a reservoir), plus a 1 % random
  refresh per step so they do not bunch up.

  Flux numerics: rho_v and IcePhase are interpolated linearly (Delaunay triangulation of the mesh nodes) to a
  --grid um grid; points outside the mesh take the nearest inside value, so the one-cell Gaussian smoothing
  applied before the gradient does not mix fake values in at the walls (without this fill, spurious wall
  gradients set the speed scale and the real flow barely moves).

LIMITS
  2D only. Reads vtkOut/solV_*.vts through pf.pvd (snapshot times); run plot_fields.py first. The .vts files hold
  B-spline control-point values, not point values (see README), which is fine at this resolution.
  Needs ffmpeg on PATH for the mp4 (H.264, yuv420p, plays in PowerPoint and QuickTime).

Usage:
    python3 postprocess/vapor_movie.py --dir <run>                 # -> <run>/plots/vapor_movie.mp4
    python3 postprocess/vapor_movie.py --dir <run> --tracers       # -> <run>/plots/vapor_movie_tracers.mp4
Options: --save, --fps, --subframes, --n-particles, --dmax, --grid, --sigma-lim, --every, --seed, --dpi.
Also writes <movie>_first.png (use as the PowerPoint poster frame) and <movie>_mid.png.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile

import matplotlib
import numpy as np

matplotlib.use("Agg")
import cmocean  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.path import Path  # noqa: E402
from scipy.interpolate import LinearNDInterpolator, RegularGridInterpolator  # noqa: E402
from scipy.ndimage import distance_transform_edt, gaussian_filter, label  # noqa: E402
from scipy.spatial import Delaunay  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pplib import RHO_AIR_DEFAULT, auto_time_unit, in_time_unit, opt_float, read_opts, read_vts, supersaturation  # noqa: E402

REGOLITH = "#b8ab95"
TRACER = "#f2b134"
TAIL = 5  # positions kept per tracer tail


def snapshots(run):
    """[(t [s], path)] from pf.pvd, which plot_fields.py writes next to vtkOut/."""
    pvd = os.path.join(run, "pf.pvd")
    if not os.path.isfile(pvd):
        sys.exit(f"vapor_movie: {pvd} not found; run plot_fields.py first")
    txt = open(pvd).read()
    steps = [(float(t), os.path.join(run, f)) for t, f in re.findall(r'timestep="([^"]+)"[^>]*file="([^"]+)"', txt)]
    if not steps:
        sys.exit(f"vapor_movie: no snapshots listed in {pvd}")
    return steps


def ice_area(X, Y, p):
    """Integral of phi over the mesh: shoelace area of each quad times the corner-mean phi [m^2 if X, Y in m]."""
    x0, x1, x2, x3 = X[:-1, :-1], X[:-1, 1:], X[1:, 1:], X[1:, :-1]
    y0, y1, y2, y3 = Y[:-1, :-1], Y[:-1, 1:], Y[1:, 1:], Y[1:, :-1]
    A = 0.5 * np.abs((x0 * y1 - x1 * y0) + (x1 * y2 - x2 * y1) + (x2 * y3 - x3 * y2) + (x3 * y0 - x0 * y3))
    return float((A * 0.25 * (p[:-1, :-1] + p[:-1, 1:] + p[1:, 1:] + p[1:, :-1])).sum())


class FluxGrid:
    """Regular grid over the mesh, with the domain mask and the flux-direction field of one snapshot."""

    def __init__(self, X, Y, dx):
        self.dx = dx
        poly = np.vstack([np.c_[X[0], Y[0]], np.c_[X[:, -1], Y[:, -1]],
                          np.c_[X[-1, ::-1], Y[-1, ::-1]], np.c_[X[::-1, 0], Y[::-1, 0]]])
        self.xg = np.arange(X.min(), X.max() + dx / 2, dx)
        self.yg = np.arange(Y.min(), Y.max() + dx / 2, dx)
        self.XG, self.YG = np.meshgrid(self.xg, self.yg)
        self.inside = Path(poly).contains_points(np.c_[self.XG.ravel(), self.YG.ravel()]).reshape(self.XG.shape)
        _, (self.fi, self.fj) = distance_transform_edt(~self.inside, return_indices=True)
        # distance to the mesh edge [um]; pad with outside cells so the array border counts as an edge
        dist = distance_transform_edt(np.pad(self.inside, 1))[1:-1, 1:-1] * dx
        self.dist = dist
        gy, gx = np.gradient(dist, dx, dx)                        # points inward
        self.nin = (gx, gy)
        self.tri = Delaunay(np.c_[X.ravel(), Y.ravel()])

    def field(self, rhov, phi, sig=None):
        rv = LinearNDInterpolator(self.tri, rhov.ravel())(self.XG, self.YG)
        ph = LinearNDInterpolator(self.tri, np.clip(phi, 0, 1).ravel())(self.XG, self.YG)
        rv = np.where(np.isnan(rv), np.nanmean(rv), rv)           # convex-hull gaps only
        ph = np.where(self.inside & ~np.isnan(ph), ph, 1.0)
        rvf = gaussian_filter(rv[self.fi, self.fj], 1.0)
        gy, gx = np.gradient(rvf, self.dx, self.dx)
        air = self.inside & (ph < 0.5)
        ux = np.where(air, -(1 - ph) * gx, 0.0)
        uy = np.where(air, -(1 - ph) * gy, 0.0)
        if sig is None:
            return ux.astype(np.float32), uy.astype(np.float32), air
        sg = LinearNDInterpolator(self.tri, np.abs(sig).ravel())(self.XG, self.YG)
        sg = np.where(air & ~np.isnan(sg), sg, 0.0)
        return ux.astype(np.float32), uy.astype(np.float32), air, sg.astype(np.float32)

    def interp(self, f):
        return RegularGridInterpolator((self.yg, self.xg), f, bounds_error=False, fill_value=0.0)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dir", required=True, help="run folder (with pf.pvd and vtkOut/)")
    ap.add_argument("--save", help="output .mp4 (default <run>/plots/vapor_movie[_tracers].mp4)")
    ap.add_argument("--tracers", action="store_true", help="add tracer particles along the vapor flux")
    ap.add_argument("--fps", type=float, default=None, help="frames per second (default 10, or 10*subframes with tracers)")
    ap.add_argument("--subframes", type=int, default=3, help="tracer steps per snapshot (default 3)")
    ap.add_argument("--n-particles", type=int, default=None,
                    help="most tracers on screen (default: --density times the air area in view)")
    ap.add_argument("--density", type=float, default=2.0, help="tracers per 1000 um^2 of air in view (default 2)")
    ap.add_argument("--dmax", type=float, default=None,
                    help="tracer step at the reference flux [um] (default: from --speed)")
    ap.add_argument("--speed", type=float, default=0.05,
                    help="tracer speed at the reference flux, as a fraction of the view width per second of video "
                         "(default 0.05: 20 s to cross the view)")
    ap.add_argument("--uniform-tracers", action="store_true",
                    help="old behavior: tracers spread uniformly over the air and a fixed number of them; default "
                         "places them in proportion to |sigma| and scales their number with the frame's mean |sigma|")
    ap.add_argument("--zoom-ice", action="store_true",
                    help="crop to where ice is during the run (phi >= 0.5 in any snapshot), plus --margin; each "
                         "separate ice cluster gets its own panel, all at the same scale and height")
    ap.add_argument("--margin", type=float, default=0.6,
                    help="zoom margin as a fraction of the ice extent (default 0.6; at least 20 um)")
    ap.add_argument("--grid", type=float, default=1.0, help="flux grid spacing [um] (default 1)")
    ap.add_argument("--sigma-lim", type=float, default=None,
                    help="color limit for |sigma| (default: 99.5th percentile of |sigma| in the air over the run)")
    ap.add_argument("--every", type=int, default=1, help="use every Nth snapshot")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--dpi", type=int, default=160,
                    help="frame resolution (default 160: a 10 in wide frame is 1600 px; use ~320 for full-resolution vtkOut)")
    a = ap.parse_args()

    run = os.path.abspath(a.dir)
    if shutil.which("ffmpeg") is None:
        sys.exit("vapor_movie: ffmpeg not on PATH")
    steps = snapshots(run)[:: a.every]
    rho_air = float(opt_float(read_opts(run), "-rho_air", RHO_AIR_DEFAULT) or RHO_AIR_DEFAULT)
    sub = a.subframes if a.tracers else 1
    fps = a.fps or (10.0 * sub)
    out = a.save or os.path.join(run, "plots", "vapor_movie_tracers.mp4" if a.tracers else "vapor_movie.mp4")
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    rng = np.random.default_rng(a.seed)

    f0, X, Y = read_vts(steps[0][1], want=("IcePhase",))
    X, Y = X * 1e6, Y * 1e6
    if np.ptp(Y) == 0 or "IcePhase" not in f0:
        sys.exit("vapor_movie: needs a 2D run with IcePhase")
    unit = auto_time_unit(steps[-1][0])
    print(f"vapor_movie: {len(steps)} snapshots, mesh {X.shape}, rho_air {rho_air}, tracers {a.tracers}")

    # pass 1: color limit and (with tracers) the flux fields and their reference magnitude
    samp, flux = [], []
    grid = FluxGrid(X, Y, a.grid) if a.tracers else None
    union = np.zeros(X.shape, bool)
    for _, fn in steps:
        f, _, _ = read_vts(fn, want=("IcePhase", "VaporDensity", "Temperature"))
        sg = supersaturation(f["VaporDensity"], f["Temperature"], rho_air)
        s = np.abs(sg)[f["IcePhase"] < 0.5]
        samp.append(rng.choice(s, min(s.size, 50000), replace=False) if s.size else s)
        union |= f["IcePhase"] >= 0.5
        if a.tracers:
            flux.append(grid.field(f["VaporDensity"], f["IcePhase"], None if a.uniform_tracers else sg))
    lim = a.sigma_lim or float(np.percentile(np.concatenate(samp), 99.5)) or 1e-12
    print(f"  sigma color limit +/-{lim:.4g}")
    # view window: whole mesh, or the ice's extent over the run plus a margin
    views = [(X.min(), X.max())]            # x ranges of the panels; all panels share y range vy0..vy1
    vy0, vy1 = Y.min(), Y.max()
    if a.zoom_ice and union.any():
        lab, nl = label(union, structure=np.ones((3, 3)))
        boxes = []
        for j in range(1, nl + 1):
            m = lab == j
            x0, x1, y0, y1 = X[m].min(), X[m].max(), Y[m].min(), Y[m].max()
            mg = max(20.0, a.margin * max(x1 - x0, y1 - y0))
            boxes.append([max(X.min(), x0 - mg), min(X.max(), x1 + mg), max(Y.min(), y0 - mg), min(Y.max(), y1 + mg)])
        boxes.sort()
        merged = [boxes[0]]
        for b in boxes[1:]:                 # overlapping x ranges become one panel
            if b[0] <= merged[-1][1]:
                merged[-1] = [merged[-1][0], max(merged[-1][1], b[1]), min(merged[-1][2], b[2]), max(merged[-1][3], b[3])]
            else:
                merged.append(b)
        views = [(b[0], b[1]) for b in merged]
        vy0, vy1 = min(b[2] for b in merged), max(b[3] for b in merged)
        print(f"  zoom: {len(views)} panel(s), x " + ", ".join(f"{u:.0f}-{v:.0f}" for u, v in views)
              + f" um, y {vy0:.0f}-{vy1:.0f} um (ice over the run plus margin)")
    vx0, vx1 = views[0][0], views[-1][1]
    vwidth = sum(v - u for u, v in views)

    def inview(q):
        ok = (q[:, 1] >= vy0) & (q[:, 1] <= vy1)
        return ok & np.any([(q[:, 0] >= u) & (q[:, 0] <= v) for u, v in views], axis=0)
    uref = 1.0
    if a.tracers:
        p90 = [np.percentile(np.hypot(ux, uy)[air], 90) if air.any() else 0.0 for ux, uy, air, *_ in flux]
        uref = float(np.median(p90)) or 1.0
        if a.dmax is None:
            a.dmax = a.speed * vwidth / fps
        gin = inview(np.c_[grid.XG.ravel(), grid.YG.ravel()]).reshape(grid.XG.shape)
        air_area = float((flux[0][2] & gin).sum()) * a.grid ** 2
        nmax = a.n_particles or max(20, int(round(a.density * air_area / 1000.0)))
        if a.uniform_tracers:
            ntrace = [nmax] * len(flux)
        else:
            msig = np.array([float(fl[3][fl[2] & gin].mean()) if (fl[2] & gin).any() else 0.0 for fl in flux])
            ntrace = [max(nmax // 10, int(round(nmax * m / (msig.max() or 1.0)))) for m in msig]
        print(f"  tracer reference |phi_a grad rho_v| = {uref:.4g} kg m^-3 um^-1 ({a.dmax:.3g} um per step, "
              f"{a.dmax * fps:.3g} um per s of video); up to {nmax} tracers "
              f"({'uniform' if a.uniform_tracers else 'placed by |sigma|, count scaled by mean |sigma|'})")

    def seed(n, air_i, wsig=None):
        if wsig is not None:   # pick grid cells with probability ~ |sigma| (air in view only), jitter inside the cell
            w = np.clip(np.nan_to_num(np.where(gin, wsig, 0.0).ravel().astype(float), nan=0.0, posinf=0.0), 0.0, None)
            if np.isfinite(w.sum()) and w.sum() > 0:
                k = rng.choice(w.size, n, p=w / w.sum())
                return np.c_[grid.XG.ravel()[k] + rng.uniform(-0.5, 0.5, n) * a.grid,
                             grid.YG.ravel()[k] + rng.uniform(-0.5, 0.5, n) * a.grid]
        pts = []
        for tries in range(200):   # bounded: once ice fills (nearly) all air in view there is no spot to find
            if len(pts) >= n:
                break
            c = np.c_[rng.uniform(vx0, vx1, 4 * n + 16), rng.uniform(vy0, vy1, 4 * n + 16)]
            pts.extend(c[(air_i(c[:, ::-1]) > (0.99 if tries < 100 else 0.5)) & inview(c)][: n - len(pts)])
        if len(pts) < n:           # no air left: park the rest in view; they sit under the ice and respawn
            c = np.c_[rng.uniform(vx0, vx1, 64 * n + 64), rng.uniform(vy0, vy1, 64 * n + 64)]
            c = c[inview(c)]
            pts.extend(c[rng.integers(0, len(c), n - len(pts))] if len(c) else
                       np.repeat([[0.5 * (vx0 + vx1), 0.5 * (vy0 + vy1)]], n - len(pts), axis=0))
        return np.array(pts)[:n]

    nin_x, nin_y = (grid.interp(grid.nin[0]), grid.interp(grid.nin[1])) if a.tracers else (None, None)
    dist_i = grid.interp(grid.dist) if a.tracers else None
    P, hist = None, []

    figw = 10.0
    figh = min(8.0, max(2.4, figw * (vy1 - vy0) / vwidth + 1.1))
    tmp = tempfile.mkdtemp(dir=os.path.dirname(os.path.abspath(out)))
    k, A0 = 0, None
    try:
        for i, (t, fn) in enumerate(steps):
            f, _, _ = read_vts(fn, want=("IcePhase", "VaporDensity", "Temperature"))
            p = np.clip(f["IcePhase"], 0, 1)
            sig = supersaturation(f["VaporDensity"], f["Temperature"], rho_air)
            ia = ice_area(X, Y, p)
            A0 = A0 or ia
            if a.tracers:
                ux, uy, air = flux[i][:3]
                wsig = None if a.uniform_tracers else flux[i][3]
                fx, fy, fa = grid.interp(ux / uref), grid.interp(uy / uref), grid.interp(air.astype(np.float32))
                if P is None:
                    P = seed(nmax, fa, wsig)
                    hist = [P.copy()]
            for _ in range(sub):
                if a.tracers:
                    def vel(q, fx=fx, fy=fy):
                        v = np.c_[fx(q[:, ::-1]), fy(q[:, ::-1])]
                        s = np.hypot(v[:, 0], v[:, 1])
                        return v * np.minimum(1.0, 1.0 / np.maximum(s, 1e-30))[:, None] * a.dmax

                    v0 = vel(P)
                    P = P + vel(P + 0.5 * v0)
                    q = P[:, ::-1]
                    outward = (dist_i(q) < 2.5 * a.grid) & ((v0[:, 0] * nin_x(q) + v0[:, 1] * nin_y(q)) < 0)
                    bad = (fa(q) < 0.5) | outward | ~inview(P) | (rng.random(len(P)) < 0.01)
                    if bad.any():
                        P[bad] = seed(int(bad.sum()), fa, wsig)
                        for h in hist:
                            h[bad] = P[bad]
                    hist = (hist + [P.copy()])[-TAIL:]

                fig, axs = plt.subplots(1, len(views), figsize=(figw, figh), dpi=a.dpi, squeeze=False,
                                        gridspec_kw={"width_ratios": [v - u for u, v in views], "wspace": 0.08})
                axs = axs[0]
                for j, (ax, (u, v)) in enumerate(zip(axs, views)):
                    ax.set_facecolor(REGOLITH)
                    # air: sigma with Gouraud (bilinear) shading, smooth even on a strided viewing copy; the ice
                    # (phi >= 0.5) is drawn on top as filled contours, so sigma under the ice is never seen
                    pm = ax.pcolormesh(X, Y, np.clip(sig, -lim, lim), cmap=cmocean.cm.balance, vmin=-lim, vmax=lim,
                                       shading="gouraud", rasterized=True)
                    if (p >= 0.5).any():
                        ax.contourf(X, Y, p, levels=np.linspace(0.5, 1.0, 11), cmap=cmocean.cm.ice, vmin=0.0, vmax=1.0)
                    ax.contour(X, Y, p, levels=[0.5], colors="k", linewidths=0.6)
                    if a.tracers:
                        nt = ntrace[i]   # draw the first nt tracers; the rest move unseen, so the count can change smoothly
                        vis = fa(P[:nt, ::-1]) >= 0.5   # tracers parked under ice (no air left in view) are hidden
                        ax.add_collection(LineCollection(np.stack(hist, axis=1)[:nt][vis], colors=TRACER,
                                                         linewidths=0.8, alpha=0.9))
                        ax.scatter(P[:nt][vis, 0], P[:nt][vis, 1], s=2.5, c=TRACER, edgecolors="none", zorder=5)
                    ax.set_aspect("equal")
                    ax.set_xlim(u, v)
                    ax.set_ylim(vy0 - 2, vy1 + 2)
                    ax.set_xlabel(r"$x$ [$\mu$m]")
                    if j == 0:
                        ax.set_ylabel(r"$y$ [$\mu$m]")
                    else:
                        ax.tick_params(labelleft=False)
                cb = fig.colorbar(pm, ax=list(axs), pad=0.01, fraction=0.035)
                cb.set_label(r"$\sigma=\rho_v/\rho_{vs}^I(T)-1$", fontsize=10)
                cb.formatter.set_powerlimits((-2, 2))
                tv = float(in_time_unit(t, unit))
                axs[0].set_title(rf"$t = {tv:.2f}$ {unit},   ice area {100 * ia / A0:.0f} % of initial", fontsize=13,
                                 loc="left")
                fig.savefig(os.path.join(tmp, f"f{k:05d}.png"), facecolor="white", bbox_inches="tight", pad_inches=0.08)
                plt.close(fig)
                k += 1
        subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(fps),
                        "-i", os.path.join(tmp, "f%05d.png"),
                        "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2:color=white", "-c:v", "libx264", "-pix_fmt", "yuv420p",
                        "-crf", "18", "-movflags", "+faststart", out], check=True)
        stem = os.path.splitext(out)[0]
        shutil.copy(os.path.join(tmp, "f00000.png"), stem + "_first.png")
        shutil.copy(os.path.join(tmp, f"f{k // 2:05d}.png"), stem + "_mid.png")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print(f"  {k} frames at {fps:g} fps; ice area {100 * ia / A0:.1f} % of initial at the end")
    print(f"  wrote {out}")


if __name__ == "__main__":
    main()
