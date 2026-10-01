#!/usr/bin/env python3
"""Show how a packing is built: gravity deposition, crop, seam heal, fill.

    venv_enceladus/bin/python postprocess/make_deposition_movie.py \\
        inputs/packings/keff_LR40/phi0.325_Rave50um_LR40_seed1701 [--out <dir>]
        [--fps 30] [--storyboard-only]

Writes <out>/deposition.mp4 (1920x1080) and <out>/deposition_storyboard.png
(+ .pdf), default <out> = the packing directory.

NOT AN ILLUSTRATION -- A REPLAY. The packing's own build is re-run from its
metadata.json (accepted seed, porosity-salted stream, rolling budget) with the
generator's route trace switched on, and the final grains are checked against
its grains.dat before anything is drawn. Every fall and every roll in the movie
is the one that produced the packing the solver runs on. The trace does not
change the geometry (generate_packing_gravity._settle).

WHAT IT SHOWS, in order
  1. Deposition in a strip taller than the domain. Each grain is dropped at a
     random x, falls until it touches the bed, then rolls downhill around that
     contact until it finds a second support or uses up its rolling budget
     (the porosity knob). The strip wraps in x, so a grain leaving one side
     re-enters the other; its image is drawn faintly.
  2. The domain is cut from the interior, away from the floor layer and the
     loose top surface.
  3. y is made periodic: bottom and top edges are identified and overlaps
     across that seam are relaxed apart.
  4. Fillers go into the largest voids until the porosity target is met.
  5. The finished cell tiled 3x3: it is periodic both ways.

The first grains are shown slowly with their routes; deposition then speeds up
(several grains in flight at once), which is a display choice only -- the
routes are computed strictly one grain at a time.
"""
from __future__ import annotations

import argparse
import json
import math
import shutil
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle, Rectangle

HERE = Path(__file__).resolve().parent
PROJ = HERE.parent
sys.path.insert(0, str(PROJ / "preprocess"))
sys.path.insert(0, str(HERE))
import generate_packing_gravity as gpg                      # noqa: E402
import packing_lib as pl                                    # noqa: E402
from make_packing_movie import assemble                     # noqa: E402

# The packing previews' palette (packing_lib._preview), so slide and movie match.
C_MATRIX, C_FILLER, C_IMAGE, C_EDGE = "#7fb3d5", "#e8b84b", "#c9d9e6", "#1b4f72"
C_WINDOW, C_ROUTE, C_MOVED, C_CUT = "#c0392b", "#1b4f72", "#8e44ad", "#d5d8dc"
INK, MUTED = "#1a1a1a", "#555555"


# ----------------------------------------------------------------------------
# Replay the build
# ----------------------------------------------------------------------------

def replay(pdir: Path):
    m = json.loads((pdir / "metadata.json").read_text())
    Lx, Ly = m["Lx"], m["Ly"]
    R, phi = m["mean_r_m_requested"], m["porosity_target"]
    if not (m.get("periodic_x") and m.get("periodic_y")):
        sys.exit("only --periodic xy packings are supported")
    margin = m["crop_margin_mean_r"] * R
    strip_h = Ly + 2.0 * margin

    trace = []
    rng = gpg._stream(m["seed"], phi)
    cen, rad = gpg.deposit(rng, Lx, strip_h, R, m["sigma_ln"], m["radius_clip_frac"],
                           m["roll_tol_rad"], True, trace=trace)
    strip = (cen.copy(), rad.copy())

    y = cen[:, 1] - margin
    keep = (y >= 0.0) & (y < Ly)                     # crop_window(py=True)
    wc = np.stack([np.mod(cen[keep, 0], Lx), y[keep]], axis=1)
    wr = rad[keep].copy()
    before = wc.copy()
    pl.relax(wc, wr, Lx, Ly, True, True, tol_frac=1e-4, max_iter=4000)
    healed = wc.copy()
    fc, fr, n_fill = gpg.fill_voids(gpg._stream(m["seed"], phi, salt=7919), wc, wr,
                                    Lx, Ly, True, True, phi, R, m["sigma_ln"],
                                    m["radius_clip_frac"], 512, context=None)

    # --- check against grains.dat --------------------------------------------
    rows = [l.split() for l in (pdir / "grains.dat").read_text().splitlines()
            if l.strip() and not l.startswith("#")]
    g = np.array([[float(v) for v in r] for r in rows[1:]])
    gin = g[(g[:, 0] >= 0) & (g[:, 0] < Lx) & (g[:, 1] >= 0) & (g[:, 1] < Ly)]
    if len(gin) != len(fr):
        sys.exit(f"replay mismatch: {len(fr)} grains vs {len(gin)} in grains.dat")
    a = np.lexsort((fc[:, 1], fc[:, 0])); b = np.lexsort((gin[:, 1], gin[:, 0]))
    err = max(np.abs(fc[a] - gin[b, :2]).max(), np.abs(fr[a] - gin[b, 2]).max())
    if err > 1e-3 * R:
        sys.exit(f"replay mismatch: max |difference| {err:.3e} m > 1e-3 R")
    print(f"  replay matches grains.dat: {len(fr)} grains, max |diff| {err:.1e} m")

    return dict(m=m, Lx=Lx, Ly=Ly, R=R, margin=margin, strip_h=strip_h, trace=trace,
                strip=strip, keep=keep, before=before, healed=healed, rad_w=rad[keep],
                final_c=fc, final_r=fr, n_fill=n_fill, n_matrix=int(keep.sum()))


# ----------------------------------------------------------------------------
# Drawing helpers
# ----------------------------------------------------------------------------

def circles(ax, cen, rad, face, alpha=1.0, lw=0.35, z=2, edge=C_EDGE):
    if len(rad) == 0:
        return
    pc = PatchCollection([Circle((x, y), r) for (x, y), r in zip(cen, rad)],
                         facecolor=face, edgecolor=edge, linewidth=lw, alpha=alpha,
                         zorder=z)
    ax.add_collection(pc)


def with_x_images(cen, rad, Lx):
    """Grains plus their x-periodic images where they overhang an edge."""
    c, r = [cen], [rad]
    for s in (-Lx, Lx):
        m = (cen[:, 0] + s + rad > 0) & (cen[:, 0] + s - rad < Lx)
        c.append(cen[m] + [s, 0.0]); r.append(rad[m])
    return np.concatenate(c), np.concatenate(r), len(rad)


def along(pts, f):
    """Point a fraction f of the way along a polyline, by arclength."""
    pts = np.asarray(pts)
    if len(pts) == 1:
        return pts[0]
    seg = np.hypot(*np.diff(pts, axis=0).T)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    t = f * s[-1]
    i = min(int(np.searchsorted(s, t, side="right")) - 1, len(seg) - 1)
    w = 0.0 if seg[i] == 0 else (t - s[i]) / seg[i]
    return pts[i] + w * (pts[i + 1] - pts[i])


MOVED_TOL = 0.1          # seam heal: highlight grains displaced by > 0.1 R


def roll_length(rec):
    p = np.asarray(rec[1])
    return float(np.hypot(*np.diff(p, axis=0).T).sum()) if len(p) > 1 else 0.0


def pick_showcase(trace, R, Lx, n_show=5, frac=0.35):
    """First run of n_show consecutive grains past `frac` of the build in which
    every grain rolls at least 1 R and stays 2 R clear of the x edges -- the
    ones the slow-motion segment follows. The earliest grains land on a bare
    floor and barely roll, so they would show the drop but not the roll; a
    route crossing the x wrap would draw as a jump across the frame."""
    def ok(rec):
        x = np.asarray(rec[1])[:, 0]
        return (roll_length(rec) >= R and x.min() > 2 * R and x.max() < Lx - 2 * R)
    start = int(frac * len(trace))
    for g in range(start, len(trace) - n_show):
        if all(ok(trace[g + i]) for i in range(n_show)):
            return g
    return start


def route(rec, top):
    """Full display route: release above the strip, then the traced path."""
    r, path = rec
    x0, y0 = path[0]
    return np.array([(x0, top + r)] + list(path))


# ----------------------------------------------------------------------------
# Frame layout
# ----------------------------------------------------------------------------

def new_frame(d, title, body, stats=""):
    fig = plt.figure(figsize=(16, 9), dpi=120)
    ax = fig.add_axes([0.04, 0.05, 0.50, 0.90])
    mm = 1e3
    ax.set_xlim(-0.12 * d["Lx"] * mm, 1.12 * d["Lx"] * mm)
    ax.set_ylim(-0.05 * d["strip_h"] * mm, 1.12 * d["strip_h"] * mm)
    ax.set_aspect("equal")
    ax.set_xlabel("x  [mm]", fontsize=12)
    ax.set_ylabel("height  [mm]", fontsize=12)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    tx = fig.add_axes([0.58, 0.05, 0.40, 0.90]); tx.axis("off")
    tx.text(0.0, 0.92, title, fontsize=26, fontweight="bold", color=INK, va="top")
    tx.text(0.0, 0.80, body, fontsize=15, color=INK, va="top", linespacing=1.5, wrap=True)
    if stats:
        tx.text(0.0, 0.10, stats, fontsize=13, color=MUTED, va="bottom", family="monospace")
    return fig, ax


def mm(a):
    return np.asarray(a) * 1e3


# ----------------------------------------------------------------------------
# The movie
# ----------------------------------------------------------------------------

def frames(d, out: Path, fps: int):
    fdir = out / "deposition_frames"
    if fdir.exists():
        shutil.rmtree(fdir)
    fdir.mkdir(parents=True)
    Lx, Ly, H, R = d["Lx"], d["Ly"], d["strip_h"], d["R"]
    trace = d["trace"]
    n = len(trace)
    top = 1.06 * H
    routes = [route(rec, top) for rec in trace]
    rads = np.array([rec[0] for rec in trace])
    rest = np.array([rt[-1] for rt in routes])
    rest[:, 0] = np.mod(rest[:, 0], Lx)
    phi = d["m"]["porosity_target"]
    rollbud = math.degrees(d["m"]["roll_tol_rad"])
    k = [0]

    def save(fig):
        fig.savefig(fdir / f"frame_{k[0]:05d}.png", dpi=120)
        plt.close(fig)
        k[0] += 1

    def floor(ax):
        ax.add_patch(Rectangle((mm(-0.12 * Lx), mm(-0.05 * H)), mm(1.24 * Lx), mm(0.05 * H),
                               facecolor="#e5e7e9", edgecolor="none", zorder=0))
        for xe in (0.0, Lx):
            ax.axvline(mm(xe), color=MUTED, lw=0.8, ls=(0, (4, 4)), zorder=1)

    # Timeline. Moderate pace for the first grains (the floor layer), fast
    # through most of the bed, and slow motion for a few consecutive grains
    # partway up that each roll >= 1 R, with their routes drawn.
    g0 = pick_showcase(trace, R, Lx)
    slow = set(range(g0, g0 + 5))
    dur_slow, dur_fast = int(2.0 * fps), max(4, fps // 5)
    start, dur, t = [], [], 0.0
    for g in range(n):
        start.append(int(round(t)))
        if g in slow:
            dur.append(dur_slow)
            t += dur_slow + fps // 2
        elif g < 20:
            dur.append(int(0.6 * fps))
            t += 0.25 * fps
        else:
            dur.append(dur_fast)
            t += 0.3 if g < g0 else 0.25
        if g == g0 - 1:
            t = max(t, max(s_ + d_ for s_, d_ in zip(start, dur)) + fps // 2)
    start, dur = np.array(start), np.array(dur)
    t_end = int((start + dur).max()) + fps // 2

    body = ("Grains are dropped one at a time at a random x.\n"
            "Each falls until it touches the bed, then ROLLS\n"
            "downhill around that contact until it finds a\n"
            "second support or uses up its rolling budget.\n\n"
            f"Rolling budget {rollbud:.0f}° sets the porosity\n"
            "(no roll = loosest, full roll = densest).\n\n"
            "The strip wraps in x: a grain leaving one side\n"
            "re-enters on the other.")
    for f in range(t_end):
        landed = (start + dur) <= f
        flying = (start <= f) & ~landed
        nl = int(landed.sum())
        stats = f"grains placed  {nl:4d} / {n}\nporosity target {phi:.3f}"
        if any(g in slow for g in np.nonzero(flying)[0]):
            stats = "SLOW MOTION: fall -> contact -> roll to a seat\n" + stats
        fig, ax = new_frame(d, "1  Gravity deposition", body, stats)
        floor(ax)
        c, r, _ = with_x_images(rest[landed], rads[landed], Lx)
        base = int(landed.sum())
        circles(ax, mm(c[:base]), mm(r[:base]), C_MATRIX)
        circles(ax, mm(c[base:]), mm(r[base:]), C_IMAGE, alpha=0.55)
        for g in np.nonzero(flying)[0]:
            fr = (f - start[g]) / dur[g]
            # slow grains: 60% of the time falling, 40% rolling
            rt = routes[g]
            if g in slow and len(rt) > 2:
                fall = np.hypot(*(rt[1] - rt[0]))
                total = fall + np.hypot(*np.diff(rt[1:], axis=0).T).sum()
                cut = fall / max(total, 1e-30)
                ff = (fr / 0.6) * cut if fr < 0.6 else cut + (fr - 0.6) / 0.4 * (1 - cut)
            else:
                ff = fr
            p = along(rt, min(ff, 1.0))
            pc = np.array([[np.mod(p[0], Lx), p[1]]])
            cc, rr, _ = with_x_images(pc, np.array([rads[g]]), Lx)
            circles(ax, mm(cc[:1]), mm(rr[:1]), C_MATRIX, z=4, lw=0.8)
            circles(ax, mm(cc[1:]), mm(rr[1:]), C_IMAGE, alpha=0.55, z=4)
            if g in slow:
                seen = rt[:1 + int(np.searchsorted(
                    np.concatenate([[0], np.cumsum(np.hypot(*np.diff(rt, axis=0).T))])
                    / max(np.hypot(*np.diff(rt, axis=0).T).sum(), 1e-30), ff))]
                trail = np.vstack([seen, p])
                ax.plot(mm(np.mod(trail[:, 0], Lx)), mm(trail[:, 1]), "-", color=C_ROUTE,
                        lw=1.2, alpha=0.8, zorder=3)
        save(fig)

    # 2: crop
    y0, y1 = d["margin"], d["margin"] + Ly
    keep = d["keep"]
    body2 = ("The simulation domain is cut from the INTERIOR\n"
             "of the bed, discarding the flat layer on the\n"
             "floor and the loose top surface -- both are\n"
             "boundary artifacts, not snow.")
    for f in range(int(2.0 * fps)):
        a = min(1.0, f / fps)
        fig, ax = new_frame(d, "2  Cut the domain", body2)
        floor(ax)
        c, r, base = with_x_images(rest, rads, Lx)
        inside = np.concatenate([keep] + [keep[(rest[:, 0] + s + rads > 0) & (rest[:, 0] + s - rads < Lx)]
                                          for s in (-Lx, Lx)])
        circles(ax, mm(c[inside]), mm(r[inside]), C_MATRIX)
        out_face = C_CUT if a > 0.5 else C_MATRIX
        circles(ax, mm(c[~inside]), mm(r[~inside]), out_face, alpha=1.0 - 0.6 * a)
        ax.add_patch(Rectangle((0, mm(y0)), mm(Lx), mm(Ly), fill=False, edgecolor=C_WINDOW,
                               lw=2.5 * a + 0.01, zorder=6))
        save(fig)

    # 3: seam heal
    Lh = Ly
    before, healed, rw = d["before"], d["healed"], d["rad_w"]
    disp = pl.min_image(healed - before, Lx, Ly, True, True)
    moved = np.hypot(*disp.T) > MOVED_TOL * R
    body3 = ("y is made periodic by identifying the bottom\n"
             "edge with the top. Grains overlapping across\n"
             "that seam are pushed apart; the push spreads\n"
             "into the bed (purple: moved > 0.1 R).\n\n"
             "A seam that ends up contact-poor is rejected\n"
             "and the packing is rebuilt.")
    for f in range(int(2.5 * fps)):
        w = min(1.0, max(0.0, (f - 0.5 * fps) / fps))
        cur = before + w * pl.min_image(healed - before, Lx, Ly, True, True)
        fig, ax = new_frame(d, "3  Close the y seam", body3,
                            f"grains moved > 0.1 R  {int(moved.sum())}\n"
                            f"largest move  {np.hypot(*disp.T).max() / R:.2f} R")
        shift = np.array([0.0, d["margin"]])
        allc, allr = [], []
        for sy in (-Lh, 0.0, Lh):
            for sx in (-Lx, 0.0, Lx):
                allc.append(cur + [sx, sy] + shift); allr.append(rw)
        allc, allr = np.concatenate(allc), np.concatenate(allr)
        centre = np.zeros(len(allr), bool); centre[4 * len(rw):5 * len(rw)] = True
        vis = (allc[:, 1] + allr > y0 - 3 * R) & (allc[:, 1] - allr < y1 + 3 * R) & \
              (allc[:, 0] + allr > -0.12 * Lx) & (allc[:, 0] - allr < 1.12 * Lx)
        circles(ax, mm(allc[vis & ~centre]), mm(allr[vis & ~centre]), C_IMAGE, alpha=0.5)
        mv = np.zeros(len(allr), bool); mv[4 * len(rw):5 * len(rw)] = moved
        circles(ax, mm(allc[centre & ~mv]), mm(allr[centre & ~mv]), C_MATRIX)
        circles(ax, mm(allc[centre & mv]), mm(allr[centre & mv]),
                C_MOVED if w > 0 else C_MATRIX, z=3)
        ax.add_patch(Rectangle((0, mm(y0)), mm(Lx), mm(Ly), fill=False, edgecolor=C_WINDOW,
                               lw=2.5, zorder=6))
        for yy in (y0, y1):
            ax.axhline(mm(yy), color=C_WINDOW, lw=1.0, ls=(0, (3, 3)), zorder=5)
        save(fig)

    # 4: fill
    nm = d["n_matrix"]
    fc, fr = d["final_c"], d["final_r"]
    body4 = ("Fillers are placed into the largest remaining\n"
             "voids, biggest first, until the porosity\n"
             "target is met. This also evens out the pore\n"
             "sizes. Fillers need not rest on anything: in\n"
             "3D they would be held out of plane.")
    nfill = d["n_fill"]
    per = max(1, int(round(1.8 * fps / max(nfill, 1))))
    for f in range(per * nfill + fps):
        nshow = min(nfill, f // per + (1 if f % per else 0))
        fig, ax = new_frame(d, "4  Fill the largest voids", body4,
                            f"fillers  {nshow:3d} / {nfill}")
        shift = np.array([0.0, d["margin"]])
        circles(ax, mm(fc[:nm] + shift), mm(fr[:nm]), C_MATRIX)
        circles(ax, mm(fc[nm:nm + nshow] + shift), mm(fr[nm:nm + nshow]), C_FILLER, z=3)
        ax.add_patch(Rectangle((0, mm(y0)), mm(Lx), mm(Ly), fill=False, edgecolor=C_WINDOW,
                               lw=2.5, zorder=6))
        save(fig)

    # 5: periodic tiling
    body5 = ("The finished cell is periodic in both x and y:\n"
             "tiled, it continues seamlessly. This is the\n"
             "cell the phase-field solver and the k_eff\n"
             "homogenization run on.\n\n"
             f"{len(fr)} grains, porosity {d['m']['porosity_achieved']:.3f},\n"
             f"L/R = {Lx / R:.0f}")
    for f in range(int(3.0 * fps)):
        a = min(1.0, f / fps)
        fig = plt.figure(figsize=(16, 9), dpi=120)
        ax = fig.add_axes([0.04, 0.05, 0.50, 0.90])
        ax.set_aspect("equal")
        ax.set_xlim(mm(-1.0 * Lx), mm(2.0 * Lx)); ax.set_ylim(mm(-1.0 * Ly), mm(2.0 * Ly))
        ax.axis("off")
        for sy in (-1, 0, 1):
            for sx in (-1, 0, 1):
                c = fc + [sx * Lx, sy * Ly]
                if sx == 0 and sy == 0:
                    circles(ax, mm(c[:nm]), mm(fr[:nm]), C_MATRIX)
                    circles(ax, mm(c[nm:]), mm(fr[nm:]), C_FILLER)
                else:
                    circles(ax, mm(c), mm(fr), C_IMAGE, alpha=0.55 * a)
        ax.add_patch(Rectangle((0, 0), mm(Lx), mm(Ly), fill=False, edgecolor=C_WINDOW,
                               lw=2.5, zorder=6))
        tx = fig.add_axes([0.58, 0.05, 0.40, 0.90]); tx.axis("off")
        tx.text(0.0, 0.92, "5  A periodic cell", fontsize=26, fontweight="bold",
                color=INK, va="top")
        tx.text(0.0, 0.80, body5, fontsize=15, color=INK, va="top", linespacing=1.5)
        save(fig)
    return fdir


# ----------------------------------------------------------------------------
# Storyboard (one figure for a slide)
# ----------------------------------------------------------------------------

def storyboard(d, out: Path):
    Lx, Ly, H, R = d["Lx"], d["Ly"], d["strip_h"], d["R"]
    trace = d["trace"]
    rads = np.array([rec[0] for rec in trace])
    rest = np.array([rec[1][-1] for rec in trace])
    rest[:, 0] = np.mod(rest[:, 0], Lx)
    y0 = d["margin"]
    fig, axes = plt.subplots(1, 4, figsize=(20, 5.8))
    titles = ["1  Drop and roll", "2  Cut from the interior",
              "3  Close the y seam", "4  Fill voids: periodic cell"]
    for ax, t in zip(axes, titles):
        ax.set_aspect("equal"); ax.axis("off")
        ax.set_title(t, fontsize=16, loc="left", fontweight="bold")

    # 1: zoom on the bed surface where the showcase grains roll into seats
    ax = axes[0]
    g0 = pick_showcase(trace, R, Lx, n_show=3)
    nb = g0
    c, r, base = with_x_images(rest[:nb], rads[:nb], Lx)
    circles(ax, mm(c[:base]), mm(r[:base]), C_MATRIX)
    circles(ax, mm(c[base:]), mm(r[base:]), C_IMAGE, alpha=0.55)
    ends = []
    for g in range(g0, g0 + 3):
        rt = np.array(trace[g][1])
        top = rt[0, 1] + 6 * R
        full = np.vstack([[rt[0, 0], top], rt])
        xs = np.mod(full[:, 0], Lx)
        ax.plot(mm(xs[:2]), mm(full[:2, 1]), "--", color=C_ROUTE, lw=1.2)
        ax.plot(mm(xs[1:]), mm(full[1:, 1]), "-", color=C_WINDOW, lw=2.0, zorder=5)
        circles(ax, mm([[xs[1], full[1, 1]]]), mm([trace[g][0]]), "none", z=4, lw=1.0,
                edge=C_ROUTE)
        circles(ax, mm([[xs[-1], full[-1, 1]]]), mm([trace[g][0]]), C_MATRIX, z=4, lw=1.2)
        ends.append((xs[-1], full[-1, 1]))
    ends = np.array(ends)
    cx, cy = ends[:, 0].mean(), ends[:, 1].mean()
    half = max(8 * R, 0.6 * np.ptp(ends[:, 0]) + 4 * R)
    ax.set_xlim(mm(cx - half), mm(cx + half)); ax.set_ylim(mm(cy - 1.2 * half), mm(cy + 0.8 * half))
    ax.text(0.02, 0.02, "dashed: fall   red: roll to a seat\n(outline: first contact)",
            transform=ax.transAxes, fontsize=11, color=MUTED, va="bottom")

    # 2: full strip + window
    ax = axes[1]
    keep = d["keep"]
    circles(ax, mm(rest[keep]), mm(rads[keep]), C_MATRIX)
    circles(ax, mm(rest[~keep]), mm(rads[~keep]), C_CUT, alpha=0.6)
    ax.add_patch(Rectangle((0, mm(y0)), mm(Lx), mm(Ly), fill=False, edgecolor=C_WINDOW, lw=2.5))
    ax.set_xlim(mm(-0.1 * Lx), mm(1.1 * Lx)); ax.set_ylim(mm(-0.05 * H), mm(1.05 * H))

    # 3: healed window with images and moved grains
    ax = axes[2]
    before, healed, rw = d["before"], d["healed"], d["rad_w"]
    disp = pl.min_image(healed - before, Lx, Ly, True, True)
    moved = np.hypot(*disp.T) > MOVED_TOL * R
    for sy in (-Ly, Ly):
        circles(ax, mm(healed + [0, sy]), mm(rw), C_IMAGE, alpha=0.5)
    circles(ax, mm(healed[~moved]), mm(rw[~moved]), C_MATRIX)
    circles(ax, mm(healed[moved]), mm(rw[moved]), C_MOVED, z=3)
    for (x, y), (dx, dy) in zip(before[moved], disp[moved]):
        if np.hypot(dx, dy) > 0.3 * R:
            ax.annotate("", xy=mm([x + dx, y + dy]), xytext=mm([x, y]),
                        arrowprops=dict(arrowstyle="->", color=INK, lw=0.9), zorder=5)
    ax.add_patch(Rectangle((0, 0), mm(Lx), mm(Ly), fill=False, edgecolor=C_WINDOW, lw=2.5))
    ax.set_xlim(mm(-0.05 * Lx), mm(1.05 * Lx)); ax.set_ylim(mm(-0.25 * Ly), mm(1.25 * Ly))

    # 4: final + tiling hint
    ax = axes[3]
    fc, fr, nm = d["final_c"], d["final_r"], d["n_matrix"]
    for sy in (-1, 0, 1):
        for sx in (-1, 0, 1):
            if sx or sy:
                circles(ax, mm(fc + [sx * Lx, sy * Ly]), mm(fr), C_IMAGE, alpha=0.5)
    circles(ax, mm(fc[:nm]), mm(fr[:nm]), C_MATRIX)
    circles(ax, mm(fc[nm:]), mm(fr[nm:]), C_FILLER, z=3)
    ax.add_patch(Rectangle((0, 0), mm(Lx), mm(Ly), fill=False, edgecolor=C_WINDOW, lw=2.5))
    ax.set_xlim(mm(-0.3 * Lx), mm(1.3 * Lx)); ax.set_ylim(mm(-0.3 * Ly), mm(1.3 * Ly))

    handles = [Circle((0, 0), 1, facecolor=C_MATRIX, edgecolor=C_EDGE),
               Circle((0, 0), 1, facecolor=C_IMAGE, edgecolor=C_EDGE, alpha=0.55),
               Circle((0, 0), 1, facecolor=C_CUT, edgecolor=C_EDGE),
               Circle((0, 0), 1, facecolor=C_MOVED, edgecolor=C_EDGE),
               Circle((0, 0), 1, facecolor=C_FILLER, edgecolor=C_EDGE)]
    fig.legend(handles, ["deposited grain", "periodic image", "discarded (outside domain)",
                         "moved > 0.1 R closing the seam", "filler"],
               loc="lower center", ncol=5, fontsize=12, frameon=False)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    for ext in ("png", "pdf"):
        fig.savefig(out / f"deposition_storyboard.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out / 'deposition_storyboard.png'} (+ .pdf)")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("packing", type=Path)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--fps", type=int, default=30)
    ap.add_argument("--storyboard-only", action="store_true")
    a = ap.parse_args()
    out = a.out or a.packing
    out.mkdir(parents=True, exist_ok=True)
    d = replay(a.packing)
    storyboard(d, out)
    if a.storyboard_only:
        return 0
    fdir = frames(d, out, a.fps)
    return assemble(fdir, out / "deposition.mp4", a.fps)


if __name__ == "__main__":
    sys.exit(main())
