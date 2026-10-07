#!/usr/bin/env python3
"""Movie of one velocity-study run with its interface velocity underneath.

Top panel: the domain at each snapshot. Ice is drawn with the cmocean "ice"
map; the vapour is coloured by supersaturation sigma = (rho_v - rho_vs)/rho_vs
on the diverging cmocean "balance" map, symmetric about zero with one range for
the whole movie. These are the conventions of postprocess/make_movie.py.

Bottom panel: the interface velocity against time, exactly as
compare_velocity.py plots it (simulation solid, theory dashed; the mean
meniscus velocity on the channel, the two centreline velocities on the wedge),
with a marker that follows the frame.

One frame per stored snapshot. The t = 0 snapshot is skipped, so the movie
opens on the first solved step.

Usage:
    python3 movie_with_velocity.py --dir <run folder> [--out movie.mp4] [--fps 8]
"""
import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cmocean

import compare_velocity as CV
import make_theory_figures as T

sys.path.insert(0, os.path.join(T.HERE, "..", "..", "..", "postprocess"))
from pplib import supersaturation, step_times             # noqa: E402

DAY = T.DAY
NM_DAY = T.NM_DAY


def snapshots(run_dir, n_long=900):
    """Yield (step, X, Y, phi, sigma) on the true NURBS geometry."""
    from igakit.io import PetIGA
    io = PetIGA()
    nrb = io.read(os.path.join(run_dir, "igasol.dat"))
    rng = lambda i: (float(nrb.knots[i][nrb.degree[i]]),
                     float(nrb.knots[i][-nrb.degree[i] - 1]))
    u = np.linspace(*rng(0), n_long)
    v = np.linspace(*rng(1), max(120, n_long // 3))
    for f in sorted(x for x in os.listdir(run_dir) if re.fullmatch(r"sol_\d+\.dat", x)):
        step = int(f[4:9])
        if step == 0:
            continue
        C, F = nrb(u, v, fields=io.read_vec(os.path.join(run_dir, f), nrb))
        yield (step, C[..., 0] * 1e6, C[..., 1] * 1e6, F[..., 0],
               supersaturation(F[..., 2], F[..., 1]))


def load_run(run_dir):
    """This run's velocity series, via compare_velocity's loader."""
    run_dir = os.path.abspath(run_dir.rstrip("/"))
    root = os.path.dirname(os.path.dirname(run_dir))
    for r in CV.load(root):
        if os.path.abspath(r["dir"]) == run_dir:
            return r
    raise SystemExit("no velocity CSV for %s -- run the per-run post-processing first" % run_dir)


def describe(r):
    geom = "channel" if r["geom"] == "channel" else "wedge"
    return r"%s,  $\theta = %d^\circ$,  $\sigma_\infty = %s$,  %s" % (
        geom, r["theta"],
        "0" if r["sigma"] == 0 else r"%+d\times10^{-5}" % round(r["sigma"] * 1e5),
        CV.AC_LABEL[r["ac"]])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True, help="run folder")
    ap.add_argument("--out", default=None, help="output .mp4 (default: <run>/plots/)")
    ap.add_argument("--fps", type=int, default=8)
    ap.add_argument("--still", type=float, default=None,
                    help="also save the frame nearest this time [days] as a PNG")
    args = ap.parse_args()

    r = load_run(args.dir)
    out = args.out or os.path.join(args.dir, "plots", "movie_with_velocity.mp4")
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)

    times = step_times(args.dir)
    frames = [(times[s] / DAY, X, Y, phi, sig)
              for s, X, Y, phi, sig in snapshots(args.dir) if s in times]
    if len(frames) < 3:
        raise SystemExit("need at least 3 snapshots")
    X, Y = frames[0][1], frames[0][2]
    t_end = max(r["t"].max() / DAY, frames[-1][0])      # a meniscus may vanish early

    # one supersaturation range for the whole movie, from the vapour only
    vap = np.concatenate([np.abs(f[4][f[3] < 0.05]).ravel() for f in frames
                          if (f[3] < 0.05).any()])
    smax = float(np.percentile(vap, 99.5)) or 1e-6

    keys = list(r["v"])
    colors = {"mean": T.BLUE, "inner": T.BLUE, "outer": T.ORANGE}
    aspect = (X.max() - X.min()) / (Y.max() - Y.min())
    width = 11.0
    h_top = width * 0.86 / aspect
    fig = plt.figure(figsize=(width, h_top + 4.3), dpi=150)
    gs = fig.add_gridspec(2, 2, height_ratios=[h_top, 3.0], width_ratios=[1, 0.018],
                          left=0.115, right=0.9, top=0.9, bottom=0.115, hspace=0.42, wspace=0.03)
    ax = fig.add_subplot(gs[0, 0])
    cax = fig.add_subplot(gs[0, 1])
    av = fig.add_subplot(gs[1, 0])

    # --- velocity panel (static part) ---------------------------------------
    td = r["t"] / DAY
    m = td >= 1.0
    allv = np.concatenate([r["v"][k][m] for k in keys]) * NM_DAY
    pad = 0.1 * max(allv.max() - allv.min(), 1e-9)
    for k in keys:
        av.plot(td[m], r["v_th"][k][m] * NM_DAY, color=colors[k], lw=1.1, ls=(0, (4, 3)))
        av.plot(td[m], r["v"][k][m] * NM_DAY, color=colors[k], lw=2.2,
                label={"mean": "simulation", "inner": "inner meniscus",
                       "outer": "outer meniscus"}[k])
    av.plot([], [], color=T.MUTED, lw=1.1, ls=(0, (4, 3)), label="theory")
    av.axhline(0, color=T.MUTED, lw=0.8, zorder=1)
    av.set_xlim(0, t_end)
    av.set_ylim(allv.min() - pad, allv.max() + pad)
    av.set_xlabel("time [days]")
    av.set_ylabel("velocity [nm/day]")
    av.legend(fontsize=11, ncols=len(keys) + 1, loc="lower left", bbox_to_anchor=(0, 1.0))
    cursor = av.axvline(0, color=T.INK, lw=1.0)
    dots = [av.plot([], [], "o", ms=9, mfc=colors[k], mec="white", mew=2.0, zorder=5)[0]
            for k in keys]

    # --- field panel ---------------------------------------------------------
    ax.set_aspect("equal")
    ax.grid(False)
    ax.set_xlabel("x [µm]")
    ax.set_ylabel("y [µm]")
    fig.suptitle(describe(r), fontsize=14, x=0.115, ha="left", color=T.MUTED)
    stamp = ax.text(0.99, 1.03, "", transform=ax.transAxes, ha="right", va="bottom",
                    fontsize=13, color=T.INK)

    tmp = tempfile.mkdtemp(prefix="velmovie_")
    still_i = None if args.still is None else int(np.argmin([abs(f[0] - args.still) for f in frames]))
    try:
        for i, (t, X, Y, phi, sig) in enumerate(frames):
            for art in list(ax.collections):
                art.remove()
            im = ax.pcolormesh(X, Y, np.ma.masked_where(phi >= 0.5, sig), cmap=cmocean.cm.balance,
                               vmin=-smax, vmax=smax, shading="gouraud", rasterized=True)
            ax.pcolormesh(X, Y, np.ma.masked_where(phi < 0.5, phi), cmap=cmocean.cm.ice,
                          vmin=-0.6, vmax=1.25, shading="gouraud", rasterized=True)
            ax.contour(X, Y, phi, levels=[0.5], colors=[T.INK], linewidths=1.0)
            if i == 0:
                cb = fig.colorbar(im, cax=cax)
                cb.set_label(r"supersaturation $\sigma$", fontsize=12)
                cb.formatter.set_powerlimits((0, 0))
                cb.outline.set_visible(False)
            stamp.set_text("day %.0f" % t)
            cursor.set_xdata([t, t])
            for k, d in zip(keys, dots):
                if td[m][0] <= t <= td[m][-1]:
                    d.set_data([t], [np.interp(t, td[m], r["v"][k][m] * NM_DAY)])
                else:
                    d.set_data([], [])
            fig.savefig(os.path.join(tmp, "f_%05d.png" % i), facecolor="white")
            if i == still_i:
                fig.savefig(os.path.splitext(out)[0] + "_day%03d.png" % round(t), facecolor="white")
        subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(args.fps),
                        "-i", os.path.join(tmp, "f_%05d.png"),
                        "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2:color=white",
                        "-c:v", "libx264", "-pix_fmt", "yuv420p", out], check=True)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print("wrote %s  (%d frames, %.1f s at %d fps)" % (out, len(frames), len(frames) / args.fps, args.fps))


if __name__ == "__main__":
    main()
