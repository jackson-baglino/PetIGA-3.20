#!/usr/bin/env python3
"""plot_mechanisms_figure.py — the sintering-mechanisms schematic (Figure 1).

    python3 plot_mechanisms_figure.py --run <grain-pair run> [--step N]
        [--width-mm 85] [--save-dir DIR] [--copy-to DIR]

A SCHEMATIC of the two dominant sintering mechanisms, for explaining what an
experiment would see -- not a simulation result, so nothing is filled or
colour-mapped. It shows:

    the sintered pair   the phi = 0.5 outline of a simulated pair late in
                        sintering (the last snapshot of --run: the largest
                        neck simulated), mirrored across the symmetry axis
    the initial pair    two dashed circles, each a least-squares fit to its
                        grain's still-circular surface (points more than
                        --fit-exclude-deg from the neck direction)
    surface diffusion   an arrow laid along the interface, a constant
                        --offset-um outside it, running into the neck; its
                        head is drawn on the line's own last segment, so it
                        ends the line exactly
    vapor transport     a cubic spline that leaves the small grain along its
                        outward normal, crosses the pore and arrives in the
                        neck heading straight in

Labels are horizontal. The outline is the solver's interface only because it
is the right SHAPE; the figure is a diagram of mechanisms.

Style follows plot_keff_snapshots.py / plot_molaro_validation.py: Computer
Modern throughout (pplib.MANUSCRIPT_RC), ink and muted greys, arrow colours
from the cool family plot_keff.py uses. One AGU column (85 mm) by default.
Transparent background; PDF for the manuscript, PNG preview, 600 dpi.
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
from matplotlib.patches import Circle, FancyArrowPatch, Polygon
from contourpy import contour_generator
from scipy.interpolate import CubicSpline

HERE = Path(__file__).parent
REPO = HERE.parent
sys.path.insert(0, str(HERE))
import pplib                                                   # noqa: E402
from plot_keff import C_ISO, C_YY                              # noqa: E402
from plot_keff_snapshots import make_reader, snap_step, INK, FS_SMALL, MM  # noqa: E402

UM = 1e-6
C_SURF = C_ISO           # surface diffusion: deep indigo
C_VAP = C_YY             # vapour transport: steel blue
C_INIT = "#8a8a8a"       # initial-condition circles, dashed
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


def fit_circle(pts):
    """Algebraic least-squares circle (Kasa): (xc, yc, r)."""
    x, y = pts[:, 0], pts[:, 1]
    A = np.column_stack([x, y, np.ones_like(x)])
    b = x ** 2 + y ** 2
    c, *_ = np.linalg.lstsq(A, b, rcond=None)
    xc, yc = c[0] / 2, c[1] / 2
    return xc, yc, float(np.sqrt(c[2] + xc ** 2 + yc ** 2))


def arrow(ax, pts, color, lw=1.1, head=8):
    """A line through `pts` with its head on the last segment, so the head
    ends the line exactly (a FancyArrowPatch on a many-vertex path puts the
    head on a vanishing final segment and leaves the stroke running on)."""
    pts = np.asarray(pts)
    seg = np.hypot(*np.diff(pts, axis=0).T)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    # the head's straight run: the last ~3 % of the arc, at least 2 points
    k = max(1, int(np.searchsorted(s, s[-1] * 0.97)) - 1)
    ax.plot(pts[:k + 1, 0], pts[:k + 1, 1], color=color, lw=lw,
            solid_capstyle="round", zorder=3)
    ax.add_patch(FancyArrowPatch(pts[k], pts[-1], arrowstyle="-|>",
                                 mutation_scale=head, lw=lw, color=color,
                                 shrinkA=0, shrinkB=0, zorder=3))


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
          f"x = {x_neck:.1f} um")

    big = max(loops, key=len)
    # Initial condition: a circle fitted to each grain's circular surface,
    # everything within fit_exclude_deg of the neck direction left out.
    circles = []
    for cx, side in ((cs, -1), (cl, +1)):
        d = big - np.array([cx, 0.0])
        th = np.degrees(np.arctan2(d[:, 1], d[:, 0] * -side))   # 0 = away from neck
        m = (np.abs(th) < 180.0 - a.fit_exclude_deg) & (np.sign(big[:, 0] - x_neck) == side)
        xc, yc, r = fit_circle(big[m])
        circles.append((xc, yc, r))
        print(f"  initial circle: centre ({xc:.1f}, {yc:.2f}) um, r = {r:.2f} um "
              f"(opts R = {Rs if side < 0 else Rl:.1f})")

    allp = np.vstack(loops + [np.array([[xc - r, yc - r], [xc + r, yc + r]])
                              for xc, yc, r in circles])
    pad = 0.07 * (allp[:, 0].max() - allp[:, 0].min())
    x0, x1 = allp[:, 0].min() - pad, allp[:, 0].max() + pad
    y0, y1 = allp[:, 1].min() - pad, allp[:, 1].max() + pad

    W = a.width_mm * MM
    H = W * (y1 - y0) / (x1 - x0)
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_aspect("equal")
    ax.axis("off"); ax.patch.set_alpha(0.0)

    ax.plot([x0, x1], [0, 0], color=C_AXIS, lw=0.5, ls=(0, (6, 2, 1, 2)), zorder=0)
    for xc, yc, r in circles:
        ax.add_patch(Circle((xc, yc), r, fill=False, ec=C_INIT, lw=0.8,
                            ls=(0, (3, 2.5)), zorder=1))
    for ln in loops:
        ax.add_patch(Polygon(ln, closed=True, fill=False, ec=INK, lw=1.0, zorder=2))

    # --- surface diffusion: along the upper surface of the small grain ----
    sd = arc_on_contour(big, (cs, 0.0), a.sd_from, a.sd_to, a.offset_um)
    arrow(ax, sd, C_SURF)
    top = sd[int(np.argmax(sd[:, 1]))]
    ax.text(top[0], top[1] + a.label_gap_um, "surface diffusion", ha="center",
            va="bottom", fontsize=FS_SMALL, color=INK, zorder=4)

    # --- vapour transport: a spline from the small grain's lower surface into
    # the lower neck: leaves along the surface normal, arrives heading +y.
    p0 = point_on_contour(big, (cs, 0.0), a.vap_from)
    n0 = p0 - np.array([cs, 0.0]); n0 /= np.linalg.norm(n0)
    p0 = p0 + n0 * a.offset_um
    p3 = np.array([x_neck, -0.5 * w - a.offset_um])
    mid = np.array([0.5 * (p0[0] + p3[0]), min(p0[1], p3[1]) - a.vap_reach_um])
    P = np.array([p0, mid, p3])
    u = np.concatenate([[0.0], np.cumsum(np.hypot(*np.diff(P, axis=0).T))])
    L = u[-1]
    bc = ((1, n0 * L), (1, np.array([0.0, 1.0]) * L))          # d/du end slopes
    spl = [CubicSpline(u / L, P[:, k], bc_type=((1, bc[0][1][k]), (1, bc[1][1][k])))
           for k in (0, 1)]
    tt = np.linspace(0, 1, 400)
    vap = np.column_stack([spl[0](tt), spl[1](tt)])
    arrow(ax, vap, C_VAP)
    lo = vap[int(np.argmin(vap[:, 1]))]
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
    p.add_argument("--vap-from", type=float, default=-105.0,
                   help="vapour arrow start, polar angle on the small grain [deg]")
    p.add_argument("--vap-reach-um", type=float, default=14.0,
                   help="how far below its ends the vapour arrow dips [um]")
    p.add_argument("--fit-exclude-deg", type=float, default=70.0,
                   help="initial-circle fit ignores points within this angle of "
                        "the neck direction [deg]")
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
