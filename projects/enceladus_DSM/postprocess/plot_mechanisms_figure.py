#!/usr/bin/env python3
"""plot_mechanisms_figure.py — the sintering-mechanisms schematic (Figure 1).

    python3 plot_mechanisms_figure.py --run <grain-pair run> [--step N]
        [--width-mm 85] [--save-dir DIR] [--copy-to DIR]

The phi = 0.5 outline of a simulated grain pair late in sintering, mirrored
across the symmetry axis, with the two mass-transport paths drawn on it:

    surface diffusion   an arrow laid ALONG the interface -- tangent to it,
                        a small constant distance outside it -- running
                        toward the neck, labelled just above
    vapor transport     an arrow that leaves the interface along its outward
                        normal, crosses the pore, and ends pointing into the
                        neck

The outline is the solver's own interface, not a drawing: the last snapshot
of --run (the largest neck simulated), read at full resolution from
sol_*.dat, contoured at phi = 0.5. Both arrows are built from that contour --
the surface-diffusion arc is the contour itself, offset outward; the vapour
arrow starts on it and ends at the measured neck -- so they sit exactly on
the interface at any size.

Style follows plot_keff_snapshots.py / plot_molaro_validation.py: Computer
Modern throughout (pplib.MANUSCRIPT_RC), ice in cmocean `ice`, ink and muted
greys, colours from the cool family plot_keff.py uses. One AGU column
(85 mm) by default. Transparent background; PDF for the manuscript, PNG
preview, 600 dpi.
"""
from __future__ import annotations

import argparse
import glob
import os
import shutil
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Polygon
from matplotlib.path import Path as MPath
from contourpy import contour_generator
import cmocean

HERE = Path(__file__).parent
REPO = HERE.parent
sys.path.insert(0, str(HERE))
import pplib                                                   # noqa: E402
from plot_keff import C_ISO, C_YY                              # noqa: E402
from plot_keff_snapshots import make_reader, snap_step, INK, FS_SMALL, MM  # noqa: E402

UM = 1e-6
C_SURF = C_ISO           # surface diffusion: deep indigo
C_VAP = C_YY             # vapour transport: steel blue
C_ICE = cmocean.cm.ice(0.93)
C_AXIS = "#b0b0b0"


def last_snapshot(run: Path, step):
    files, reader = make_reader(run, "sol" if glob.glob(str(run / "sol_*.dat")) else "vts")
    if step is not None:
        files = [f for f in files if snap_step(f) == step]
        if not files:
            sys.exit(f"  no snapshot for step {step} in {run}")
    return files[-1], reader


def outline(run: Path, step):
    """Mirrored phi = 0.5 loops [um] and the snapshot's step."""
    fn, reader = last_snapshot(run, step)
    fl, X, Y = reader(fn, want=("IcePhase",))
    phi = np.vstack([fl["IcePhase"][:0:-1], fl["IcePhase"]])
    x = X[0] / UM
    y = np.concatenate([-Y[:0:-1, 0], Y[:, 0]]) / UM
    gen = contour_generator(x, y, phi)
    loops = [ln for ln in gen.lines(0.5) if len(ln) > 50]
    return loops, snap_step(fn)


def grains(run: Path):
    """Grain centres and radii [um] from the staged opts, small grain first."""
    o = pplib.read_opts(str(run))
    cx = [float(v) / UM for v in o["-ice_grain_cx"].split(",")]
    R = [float(v) / UM for v in o["-ice_grain_R"].split(",")]
    return sorted(zip(R, cx))                     # [(R_small, cx), (R_large, cx)]


def arc_on_contour(loop, centre, th0, th1, offset, upper=True):
    """The part of `loop` between polar angles th0 -> th1 [deg] about
    `centre`, in that order, pushed `offset` um outward along the radius."""
    c = np.asarray(centre)
    d = loop - c
    th = np.degrees(np.arctan2(d[:, 1], d[:, 0]))
    lo, hi = min(th0, th1), max(th0, th1)
    m = (th >= lo) & (th <= hi) & ((loop[:, 1] > 0) if upper else (loop[:, 1] < 0))
    pts, t = loop[m], th[m]
    order = np.argsort(t)
    if th0 > th1:
        order = order[::-1]
    pts = pts[order]
    r = np.linalg.norm(pts - c, axis=1, keepdims=True)
    return c + (pts - c) * (1.0 + offset / r)


def point_on_contour(loop, centre, th):
    """The loop point nearest polar angle th [deg] about centre."""
    d = loop - np.asarray(centre)
    a = np.degrees(np.arctan2(d[:, 1], d[:, 0]))
    return loop[int(np.argmin(np.abs((a - th + 180) % 360 - 180)))]


def build(run: Path, a):
    loops, step = outline(run, a.step)
    (Rs, cs), (Rl, cl) = grains(run)
    tmap = pplib.step_times(str(run))
    n = np.genfromtxt(run / "neck_width.csv", delimiter=",", names=True)
    t_snap = tmap.get(step, float("nan"))
    i = int(np.argmin(np.abs(n["t_s"] - t_snap)))
    w = n["neck_width_m"][i] / UM
    x_neck = n["x_neck_m"][i] / UM
    print(f"  snapshot step {step}, t = {t_snap:.0f} s: neck width {w:.2f} um at "
          f"x = {x_neck:.1f} um; grains R = {Rs:.1f} / {Rl:.1f} um")

    big = max(loops, key=len)
    allp = np.vstack(loops)
    pad = 0.06 * (allp[:, 0].max() - allp[:, 0].min())
    x0, x1 = allp[:, 0].min() - pad, allp[:, 0].max() + pad
    y0, y1 = allp[:, 1].min() - pad, allp[:, 1].max() + pad

    W = a.width_mm * MM
    H = W * (y1 - y0) / (x1 - x0)
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_aspect("equal")
    ax.axis("off"); ax.patch.set_alpha(0.0)

    # Symmetry axis: the mirrored half is the same section, reflected.
    ax.plot([x0, x1], [0, 0], color=C_AXIS, lw=0.5, ls=(0, (6, 2, 1, 2)), zorder=0)
    for ln in loops:
        ax.add_patch(Polygon(ln, closed=True, fc=C_ICE, ec=INK, lw=0.9, zorder=1))

    # --- surface diffusion: along the upper surface of the small grain ----
    sd = arc_on_contour(big, (cs, 0.0), a.sd_from, a.sd_to, a.offset_um)
    ax.add_patch(FancyArrowPatch(path=MPath(sd), arrowstyle="-|>", mutation_scale=8,
                                 lw=1.1, color=C_SURF, zorder=3,
                                 shrinkA=0, shrinkB=0, capstyle="round"))
    k = len(sd) // 2
    tang = sd[min(k + 3, len(sd) - 1)] - sd[max(k - 3, 0)]
    ang = np.degrees(np.arctan2(tang[1], tang[0]))
    if ang > 90: ang -= 180
    if ang < -90: ang += 180
    nrm = (sd[k] - np.array([cs, 0.0])); nrm /= np.linalg.norm(nrm)
    lab = sd[k] + nrm * a.label_gap_um
    ax.text(*lab, "surface diffusion", rotation=ang, rotation_mode="anchor",
            ha="center", va="bottom", fontsize=FS_SMALL, color=INK, zorder=4)

    # --- vapour transport: leaves the small grain's lower surface, ends in the
    # lower neck. Cubic Bezier: out along the normal, back in toward the neck.
    p0 = point_on_contour(big, (cs, 0.0), a.vap_from)
    n0 = (p0 - np.array([cs, 0.0])); n0 /= np.linalg.norm(n0)
    p0 = p0 + n0 * a.offset_um
    p3 = np.array([x_neck, -0.5 * w - a.offset_um])
    L = a.vap_reach_um
    p1 = p0 + n0 * L
    p2 = p3 + np.array([0.0, -L])
    ax.add_patch(FancyArrowPatch(path=MPath([p0, p1, p2, p3],
                                            [MPath.MOVETO, MPath.CURVE4,
                                             MPath.CURVE4, MPath.CURVE4]),
                                 arrowstyle="-|>", mutation_scale=8, lw=1.1,
                                 color=C_VAP, zorder=3, shrinkA=0, shrinkB=0,
                                 capstyle="round"))
    # Label under the arc's lowest point, clear of both grains.
    tt = np.linspace(0, 1, 200)[:, None]
    bez = ((1 - tt) ** 3 * p0 + 3 * (1 - tt) ** 2 * tt * p1
           + 3 * (1 - tt) * tt ** 2 * p2 + tt ** 3 * p3)
    lo = bez[int(np.argmin(bez[:, 1]))]
    ax.text(lo[0], lo[1] - a.label_gap_um, "vapor transport", ha="center",
            va="top", fontsize=FS_SMALL, color=INK, zorder=4)
    return fig


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run", type=Path, required=True,
                   help="axisymmetric grain-pair run (sol_*.dat or vtkOut/, "
                        "neck_width.csv)")
    p.add_argument("--step", type=int, default=None,
                   help="snapshot step (default: the last, largest neck)")
    p.add_argument("--width-mm", type=float, default=85.0,
                   help="printed width [mm] (default 85, one AGU column)")
    p.add_argument("--sd-from", type=float, default=128.0,
                   help="surface-diffusion arc start, polar angle about the "
                        "small grain's centre [deg]")
    p.add_argument("--sd-to", type=float, default=36.0, help="... and end [deg]")
    p.add_argument("--vap-from", type=float, default=-125.0,
                   help="vapour arrow start, polar angle on the small grain [deg]")
    p.add_argument("--vap-reach-um", type=float, default=60.0,
                   help="how far the vapour arrow bows into the pore [um]")
    p.add_argument("--offset-um", type=float, default=4.0,
                   help="gap between the arrows and the interface [um]")
    p.add_argument("--label-gap-um", type=float, default=4.0)
    p.add_argument("--save-dir", type=Path, default=REPO / "studies/molaro_2019/manuscript")
    p.add_argument("--copy-to", type=Path, default=None,
                   help="also copy the figure here (the manuscript Figures folder)")
    p.add_argument("--formats", nargs="+", default=["pdf", "png"])
    p.add_argument("--dpi", type=int, default=600)
    a = p.parse_args(argv)

    plt.rcParams.update(pplib.MANUSCRIPT_RC)
    fig = build(a.run.resolve(), a)
    os.makedirs(a.save_dir, exist_ok=True)
    for fmt in a.formats:
        out = a.save_dir / f"figure1_mechanisms.{fmt}"
        fig.savefig(out, dpi=a.dpi, transparent=True)
        print(f"  wrote {out}")
        if a.copy_to is not None:
            a.copy_to.mkdir(parents=True, exist_ok=True)
            shutil.copy2(out, a.copy_to / out.name)
    print(f"  {fig.get_figwidth() * 25.4:.0f} x {fig.get_figheight() * 25.4:.0f} mm")
    return 0


if __name__ == "__main__":
    sys.exit(main())
