#!/usr/bin/env python3
"""vapor_movie.py -- movie of the ice and the vapor saturation of a 2D run, optionally with tracer particles
that follow the vapor flux. Pure matplotlib + cmocean + scipy (no ParaView), one frame per vtkOut snapshot.

WHAT IT DRAWS
  * Air (IcePhase < 0.5): supersaturation sigma = rho_v / rho_vs^I(T) - 1 on cmocean 'balance', symmetric
    about 0, so white is saturated, blue is undersaturated (sublimating), red is supersaturated (depositing).
    rho_vs^I is pplib.rho_vs, which mirrors the solver's RhoVS_I, with -rho_air read from the staged opts.
  * Ice (IcePhase >= 0.5): IcePhase on cmocean 'ice'; black line at IcePhase = 0.5.
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
Options: --save, --fps, --subframes, --n-particles, --dmax, --grid, --sigma-lim, --every, --seed.
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
from scipy.ndimage import distance_transform_edt, gaussian_filter  # noqa: E402
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

    def field(self, rhov, phi):
        rv = LinearNDInterpolator(self.tri, rhov.ravel())(self.XG, self.YG)
        ph = LinearNDInterpolator(self.tri, np.clip(phi, 0, 1).ravel())(self.XG, self.YG)
        rv = np.where(np.isnan(rv), np.nanmean(rv), rv)           # convex-hull gaps only
        ph = np.where(self.inside & ~np.isnan(ph), ph, 1.0)
        rvf = gaussian_filter(rv[self.fi, self.fj], 1.0)
        gy, gx = np.gradient(rvf, self.dx, self.dx)
        air = self.inside & (ph < 0.5)
        ux = np.where(air, -(1 - ph) * gx, 0.0)
        uy = np.where(air, -(1 - ph) * gy, 0.0)
        return ux.astype(np.float32), uy.astype(np.float32), air

    def interp(self, f):
        return RegularGridInterpolator((self.yg, self.xg), f, bounds_error=False, fill_value=0.0)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dir", required=True, help="run folder (with pf.pvd and vtkOut/)")
    ap.add_argument("--save", help="output .mp4 (default <run>/plots/vapor_movie[_tracers].mp4)")
    ap.add_argument("--tracers", action="store_true", help="add tracer particles along the vapor flux")
    ap.add_argument("--fps", type=float, default=None, help="frames per second (default 10, or 10*subframes with tracers)")
    ap.add_argument("--subframes", type=int, default=3, help="tracer steps per snapshot (default 3)")
    ap.add_argument("--n-particles", type=int, default=700)
    ap.add_argument("--dmax", type=float, default=8.0, help="tracer step at the reference flux [um] (default 8)")
    ap.add_argument("--grid", type=float, default=1.0, help="flux grid spacing [um] (default 1)")
    ap.add_argument("--sigma-lim", type=float, default=None,
                    help="color limit for |sigma| (default: 99.5th percentile of |sigma| in the air over the run)")
    ap.add_argument("--every", type=int, default=1, help="use every Nth snapshot")
    ap.add_argument("--seed", type=int, default=0)
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
    for _, fn in steps:
        f, _, _ = read_vts(fn, want=("IcePhase", "VaporDensity", "Temperature"))
        s = np.abs(supersaturation(f["VaporDensity"], f["Temperature"], rho_air))[f["IcePhase"] < 0.5]
        samp.append(rng.choice(s, min(s.size, 50000), replace=False) if s.size else s)
        if a.tracers:
            flux.append(grid.field(f["VaporDensity"], f["IcePhase"]))
    lim = a.sigma_lim or float(np.percentile(np.concatenate(samp), 99.5)) or 1e-12
    print(f"  sigma color limit +/-{lim:.4g}")
    uref = 1.0
    if a.tracers:
        p90 = [np.percentile(np.hypot(ux, uy)[air], 90) if air.any() else 0.0 for ux, uy, air in flux]
        uref = float(np.median(p90)) or 1.0
        print(f"  tracer reference |phi_a grad rho_v| = {uref:.4g} kg m^-3 um^-1 ({a.dmax} um per step)")

    def seed(n, air_i):
        pts = []
        while len(pts) < n:
            c = np.c_[rng.uniform(X.min(), X.max(), 4 * n + 16), rng.uniform(Y.min(), Y.max(), 4 * n + 16)]
            pts.extend(c[air_i(c[:, ::-1]) > 0.99][: n - len(pts)])
        return np.array(pts)

    nin_x, nin_y = (grid.interp(grid.nin[0]), grid.interp(grid.nin[1])) if a.tracers else (None, None)
    dist_i = grid.interp(grid.dist) if a.tracers else None
    P, hist = None, []

    figw = 10.0
    figh = min(8.0, max(2.4, figw * np.ptp(Y) / np.ptp(X) + 1.1))
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
                ux, uy, air = flux[i]
                fx, fy, fa = grid.interp(ux / uref), grid.interp(uy / uref), grid.interp(air.astype(np.float32))
                if P is None:
                    P = seed(a.n_particles, fa)
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
                    bad = (fa(q) < 0.5) | outward | (rng.random(len(P)) < 0.01)
                    if bad.any():
                        P[bad] = seed(int(bad.sum()), fa)
                        for h in hist:
                            h[bad] = P[bad]
                    hist = (hist + [P.copy()])[-TAIL:]

                fig, ax = plt.subplots(figsize=(figw, figh), dpi=160)
                ax.set_facecolor(REGOLITH)
                pm = ax.pcolormesh(X, Y, np.ma.masked_where(p >= 0.5, sig), cmap=cmocean.cm.balance,
                                   vmin=-lim, vmax=lim, shading="nearest", rasterized=True)
                ax.pcolormesh(X, Y, np.ma.masked_where(p < 0.5, p), cmap=cmocean.cm.ice, vmin=0.0, vmax=1.0,
                              shading="nearest", rasterized=True)
                ax.contour(X, Y, p, levels=[0.5], colors="k", linewidths=0.6)
                if a.tracers:
                    ax.add_collection(LineCollection(np.stack(hist, axis=1), colors=TRACER, linewidths=0.8, alpha=0.9))
                    ax.scatter(P[:, 0], P[:, 1], s=2.5, c=TRACER, edgecolors="none", zorder=5)
                cb = fig.colorbar(pm, ax=ax, pad=0.01, fraction=0.035)
                cb.set_label(r"$\sigma=\rho_v/\rho_{vs}^I(T)-1$", fontsize=10)
                cb.formatter.set_powerlimits((-2, 2))
                tv = float(in_time_unit(t, unit))
                ax.set_title(rf"$t = {tv:.2f}$ {unit},   ice area {100 * ia / A0:.0f} % of initial", fontsize=13, loc="left")
                ax.set_aspect("equal")
                ax.set_xlim(X.min(), X.max())
                ax.set_ylim(Y.min() - 2, Y.max() + 2)
                ax.set_xlabel(r"$x$ [$\mu$m]")
                ax.set_ylabel(r"$y$ [$\mu$m]")
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
