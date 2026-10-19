#!/usr/bin/env python3
"""make_aggregate_movie.py — animate a sintering PACKING of ice grains.

The aggregate counterpart of make_neck_movie.py, with the same three panels,
colours and conventions:

  left     the whole periodic packing at one instant. The vapour field is the
           background, the ice is painted opaque on top, and the zoom window is
           marked.
  top r.   a close-up on ONE fixed grain contact -- the same two grains in every
           frame -- with the waist across the neck drawn on the frame and
           labelled in um.
  bot r.   k_eff / k_eff,r against time [days], the whole curve in grey, the
           elapsed part in colour, a dot at the current frame. k_eff is the
           mean of k_xx and k_yy, divided by its value at theta = t/tau_sub = 30
           (manuscript Figs 6-7). SSA / SSA(opening frame) is read out beside it.

THE BACKGROUND IS SUPERSATURATION, sigma = rho_v/rho_vs(T) - 1, exactly as in
make_neck_movie.py: sigma = 0 is white on a symmetric diverging bar, and it is
always on the bar. Blue is undersaturated (ice sublimates), red is
supersaturated (vapour deposits), white is saturation over flat ice.

THE CONTACT is chosen from the data, not by eye. grains.dat gives every grain
centre and radius; pairs that touch (gap < 5 % of the radii) are scored on the
first and last frames by the width of the ice chord PERPENDICULAR to the line
of centres, minimised along the line (the waist). A pair qualifies when both
grains are still ice at the last frame (neither has been eaten by coarsening),
the waist has at least doubled, ends between 0.4 and 0.85 of the smaller
grain's diameter (a neck, not a coalesced grain), and the window sits inside
the domain. The qualifier closest to 0.65 of the diameter, larger grains
preferred, is used. --contact I,J overrides; --list-contacts prints the ranking and exits.

Snapshots are the log-spaced .vts files, so equal frame spacing is NOT equal
time spacing: the stamp on the frame and the dot on the curve carry the clock.

Usage:
    python make_aggregate_movie.py <run_dir> [--out FILE.mp4] [--fps 3]
        [--frame-png STEP] [--list-contacts] [--contact I,J] [--dpi 120]
        [--tau-sub S] [--sat-clip P] [--keep-ic]
"""

import argparse
import glob
import os
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from matplotlib.colors import AsinhNorm
from scipy.ndimage import map_coordinates
import cmocean

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pplib
from pplib import read_vts, step_of, step_times, drop_ic
from make_neck_movie import (ice_alpha_cmap, centered_cmap, sigma_ticks,
                             ICE_EDGE, ACCENT, TRACE, sigma_field, WANT)

THETA_REF = 30.0
GREY = "#c4c8ce"
INK = "#1e2930"


# ---------------------------------------------------------------- geometry
def read_grains(run_dir, opts):
    """(x, y, r) arrays [m] from the packing file the run was built from."""
    rel = opts.get("-grains_file", "").strip()
    cands = [Path(run_dir) / Path(rel).name,
             Path(__file__).resolve().parents[1] / rel, Path(rel)]
    path = next((p for p in cands if p.is_file()), None)
    if path is None:
        sys.exit(f"packing file {rel} not found (looked in {[str(c) for c in cands]})")
    rows = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        v = [float(t) for t in line.split()]
        if len(v) == 3:
            rows.append(v)
    a = np.asarray(rows)
    return a[:, 0], a[:, 1], a[:, 2]


class Grid:
    """Uniform .vts grid: bilinear sampling of a field at physical points."""
    def __init__(self, X, Y):
        self.x1, self.y1 = X[0, :], Y[:, 0]
        self.x0, self.y0 = self.x1[0], self.y1[0]
        self.dx = (self.x1[-1] - self.x1[0]) / (len(self.x1) - 1)
        self.dy = (self.y1[-1] - self.y1[0]) / (len(self.y1) - 1)

    def sample(self, fld, x, y):
        return map_coordinates(fld, [(np.asarray(y) - self.y0) / self.dy,
                                     (np.asarray(x) - self.x0) / self.dx],
                               order=1, mode="nearest")


def chord_through(g, phi, p, n, level=0.5, half=170e-6, ds=0.5e-6):
    """Width of the ice chord along unit vector n through point p, and its
    end points. 0 if p itself is not ice. Crossings are linear in phi."""
    s = np.arange(-half, half + ds, ds)
    v = g.sample(phi, p[0] + s * n[0], p[1] + s * n[1])
    c = len(s) // 2
    if v[c] < level:
        return 0.0, None
    hi = c
    while hi < len(s) - 1 and v[hi + 1] >= level:
        hi += 1
    lo = c
    while lo > 0 and v[lo - 1] >= level:
        lo -= 1
    def cross(a, b):
        return s[a] + (level - v[a]) * (s[b] - s[a]) / (v[b] - v[a]) if v[b] != v[a] else s[a]
    s_hi = cross(hi, min(hi + 1, len(s) - 1)) if hi < len(s) - 1 else s[hi]
    s_lo = cross(lo, max(lo - 1, 0)) if lo > 0 else s[lo]
    return s_hi - s_lo, (p + s_lo * n, p + s_hi * n)


def waist(g, phi, ci, cj, ri, rj):
    """(width, (end1, end2), centre) of the narrowest perpendicular chord near
    the contact of grains i and j. Scans +-0.35 of the smaller radius."""
    d = np.hypot(*(cj - ci))
    u = (cj - ci) / d
    nrm = np.array([-u[1], u[0]])
    p_c = ci + u * d * ri / (ri + rj)
    best = (np.inf, None, None)
    for off in np.linspace(-0.35, 0.35, 29) * min(ri, rj):
        p = p_c + u * off
        w, ends = chord_through(g, phi, p, nrm)
        if w < best[0]:
            best = (w, ends, p)
    return best if np.isfinite(best[0]) else (0.0, None, p_c)


def candidate_pairs(gx, gy, gr, L, margin):
    """Touching pairs of well-sized grains whose contact is away from the
    domain edge (the window then never wraps)."""
    keep = np.flatnonzero(gr >= 0.7 * gr.mean())
    out = []
    for a in range(len(keep)):
        for b in range(a + 1, len(keep)):
            i, j = keep[a], keep[b]
            d = np.hypot(gx[i] - gx[j], gy[i] - gy[j])
            gap = d - (gr[i] + gr[j])
            if abs(gap) > 0.05 * (gr[i] + gr[j]):
                continue
            c = (np.array([gx[i], gy[i]]) * gr[j] + np.array([gx[j], gy[j]]) * gr[i]) / (gr[i] + gr[j])
            if np.all(c > margin) and np.all(c < L - margin):
                out.append((i, j))
    return out


def rank_contacts(g, phi_first, phi_last, gx, gy, gr, pairs):
    rows = []
    for i, j in pairs:
        ci, cj = np.array([gx[i], gy[i]]), np.array([gx[j], gy[j]])
        # both grains must still be ice after 30 d
        if g.sample(phi_last, [ci[0], cj[0]], [ci[1], cj[1]]).min() < 0.5:
            continue
        w1, _, _ = waist(g, phi_first, ci, cj, gr[i], gr[j])
        w2, _, _ = waist(g, phi_last, ci, cj, gr[i], gr[j])
        rmin = min(gr[i], gr[j])
        ratio = w2 / (2.0 * rmin)
        if not (0.4 <= ratio <= 0.85) or w2 < 2.0 * w1:   # coalesced / still a point contact / no growth
            continue
        # a clear neck: waist about 0.65 of the smaller grain's diameter; larger grains read better
        rows.append((-abs(ratio - 0.65) - 0.002 * (100.0 - rmin * 1e6), i, j, w1, w2))
    rows.sort(reverse=True)
    return rows


# ------------------------------------------------------------------ k_eff
def read_keff(run_dir):
    for name in ("k_eff_tensor.csv", "k_eff.csv"):
        p = Path(run_dir) / name
        if p.is_file():
            a = np.genfromtxt(p, delimiter=",", names=True)
            return name, a["time"], 0.5 * (a["k_00"] + a["k_11"])
    sys.exit(f"no k_eff_tensor.csv / k_eff.csv in {run_dir}")


def tau_sub_of(run_dir, given):
    if given:
        return given
    for p in glob.glob(str(Path(run_dir) / "*.opts")):
        m = re.search(r"tau_sub \((\d+(?:\.\d+)?) s\)", Path(p).read_text())
        if m:
            return float(m.group(1))
    sys.exit("tau_sub not found in the opts comments; pass --tau-sub")


def fmt_time(t_s):
    if t_s < 3600.0:
        return f"{t_s / 60.0:.0f} min"
    if t_s < 86400.0:
        return f"{t_s / 3600.0:.1f} h"
    return f"{t_s / 86400.0:.1f} d"


# -------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--fps", type=int, default=3)
    ap.add_argument("--dpi", type=int, default=120)
    ap.add_argument("--phi", type=float, default=0.5)
    ap.add_argument("--keep-ic", action="store_true")
    ap.add_argument("--frame-png", type=int, default=None,
                    help="render only this step to a PNG and exit")
    ap.add_argument("--list-contacts", action="store_true")
    ap.add_argument("--contact", default=None, help="grain indices I,J (rows of grains.dat, 0-based)")
    ap.add_argument("--tau-sub", type=float, default=None)
    ap.add_argument("--sat-clip", type=float, default=0.1)
    ap.add_argument("--sigma-scale", type=float, default=1e4,
                    help="sigma is shown times this (default 1e4)")
    ap.add_argument("--zoom-um", type=float, default=None,
                    help="zoom window half-width [um] (default 0.8 x the pair's summed radii)")
    args = ap.parse_args()
    run = args.run_dir

    files = drop_ic(sorted(glob.glob(str(run / "vtkOut" / "solV_*.vts")), key=step_of),
                    step_of, args.keep_ic, step_times(run))
    if not files:
        sys.exit(f"no solV_*.vts under {run}/vtkOut")
    tmap = step_times(run)
    t_s = np.array([tmap.get(step_of(f), np.nan) for f in files])
    opts = pplib.read_opts(str(run))
    Lx = float(opts["-Lx"])
    f0, X, Y = read_vts(files[0], want=WANT)
    g = Grid(X, Y)
    fl = read_vts(files[-1], want=("IcePhase",))[0]
    gx, gy, gr = read_grains(run, opts)

    # ---- the contact ------------------------------------------------------
    pairs = candidate_pairs(gx, gy, gr, Lx, margin=0.18e-3)
    ranked = rank_contacts(g, f0["IcePhase"], fl["IcePhase"], gx, gy, gr, pairs)
    if args.list_contacts or not args.contact:
        print(f"{len(pairs)} touching pairs inside the domain, {len(ranked)} qualify; best:")
        for sc, i, j, w1, w2 in ranked[:8]:
            print(f"  grains {i},{j}  R = {gr[i]*1e6:.1f},{gr[j]*1e6:.1f} um  "
                  f"waist {w1*1e6:.1f} -> {w2*1e6:.1f} um  score {sc:.2f}")
    if args.list_contacts:
        return
    if args.contact:
        ci_, cj_ = (int(v) for v in args.contact.split(","))
    elif ranked:
        _, ci_, cj_, _, _ = ranked[0]
    else:
        sys.exit("no contact qualifies; pass --contact I,J")
    ci, cj = np.array([gx[ci_], gy[ci_]]), np.array([gx[cj_], gy[cj_]])
    ri, rj = gr[ci_], gr[cj_]
    d = np.hypot(*(cj - ci)); u = (cj - ci) / d
    centre = ci + u * d * ri / (ri + rj)
    half = (args.zoom_um * 1e-6) if args.zoom_um else 0.8 * (ri + rj)
    print(f"contact: grains {ci_},{cj_} at ({centre[0]*1e6:.0f}, {centre[1]*1e6:.0f}) um, "
          f"window +-{half*1e6:.0f} um")

    # ---- pass 1: waist series and sigma range ------------------------------
    S = args.sigma_scale
    widths, ends, lo_hi = [], [], []
    for k, fn in enumerate(files):
        f, _, _ = read_vts(fn, want=WANT)
        phi = f["IcePhase"]
        w, e, _ = waist(g, phi, ci, cj, ri, rj)
        widths.append(w); ends.append(e)
        pore = phi < args.phi
        sg = S * pplib.supersaturation(f["VaporDensity"], f["Temperature"])[pore]
        lo_hi.append(np.percentile(sg, [args.sat_clip, 100 - args.sat_clip]))
        if k % 20 == 0:
            print(f"  pass 1: {k}/{len(files)}", flush=True)
    widths = np.array(widths) * 1e6
    lo_hi = np.array(lo_hi)
    smin, smax = float(lo_hi[:, 0].min()), float(lo_hi[:, 1].max())
    v = max(abs(smin), abs(smax))
    norm = AsinhNorm(linear_width=v / 300.0, vmin=-v, vmax=v)
    vapcm = centered_cmap(cmocean.cm.balance, norm)
    print(f"waist {widths[0]:.1f} -> {widths[-1]:.1f} um;  sigma p{args.sat_clip:g}: "
          f"{smin:+.3g} .. {smax:+.3g} (x{1/S:g})")

    # ---- k_eff and SSA series ---------------------------------------------
    kname, kt, kv = read_keff(run)
    tau = tau_sub_of(run, args.tau_sub)
    k_ref = float(np.interp(THETA_REF * tau, kt, kv))
    kn = kv / k_ref
    ssa = pplib.load_ssa(str(run))
    ssa_t, ssa_v = ssa[:, 2], ssa[:, 0]
    ssa_open = float(np.interp(t_s[0], ssa_t, ssa_v))
    print(f"k_eff from {kname}; tau_sub = {tau:g} s, k_eff,r at theta=30 (t = "
          f"{THETA_REF*tau/86400:.2f} d) = {k_ref:.4f} W/m/K")

    # ---- figure -------------------------------------------------------------
    fig = plt.figure(figsize=(16, 9))
    W, H = 16.0, 9.0
    def ax_at(l, b, w, h):  # inches -> figure fraction
        return fig.add_axes([l / W, b / H, w / W, h / H])
    axp = ax_at(0.55, 1.45, 7.0, 7.0)                 # packing, square
    axz = ax_at(8.2, 4.55, 3.7, 3.7)                  # zoom, square
    axc = ax_at(8.85, 0.95, 6.75, 3.0)                # curve
    temp_c = f"{float(opts['-temp']):g}".replace("-", "\N{MINUS SIGN}")
    fig.text(0.55 / W, 8.62 / H, "Vapor-transport sintering of a 2-D ice-grain packing  "
             rf"($\varphi$ = {1 - float(f0['IcePhase'].mean()):.3f}, "
             f"T = {temp_c} \N{DEGREE SIGN}C)", fontsize=17, fontweight="bold",
             color=INK, va="bottom")

    icecm = ice_alpha_cmap()
    um = 1e6
    ext_p = (g.x1[0] * um, g.x1[-1] * um, g.y1[0] * um, g.y1[-1] * um)

    def sel(win_lo, win_hi, arr):
        idx = np.flatnonzero((arr >= win_lo) & (arr <= win_hi))
        return idx[0], idx[-1] + 1
    zj = sel(centre[0] - half, centre[0] + half, g.x1)
    zi = sel(centre[1] - half, centre[1] + half, g.y1)
    ext_z = (g.x1[zj[0]] * um, g.x1[zj[1] - 1] * um, g.y1[zi[0]] * um, g.y1[zi[1] - 1] * um)
    zs = np.s_[zi[0]:zi[1], zj[0]:zj[1]]
    XZ, YZ = np.meshgrid(g.x1[zj[0]:zj[1]] * um, g.y1[zi[0]:zi[1]] * um)
    XP, YP = np.meshgrid(g.x1 * um, g.y1 * um)

    s0 = S * pplib.supersaturation(f0["VaporDensity"], f0["Temperature"])
    kw = dict(origin="lower", interpolation="bilinear")
    vP = axp.imshow(s0, extent=ext_p, cmap=vapcm, norm=norm, **kw)
    iP = axp.imshow(f0["IcePhase"], extent=ext_p, cmap=icecm, vmin=0, vmax=1, zorder=2, **kw)
    vZ = axz.imshow(s0[zs], extent=ext_z, cmap=vapcm, norm=norm, **kw)
    iZ = axz.imshow(f0["IcePhase"][zs], extent=ext_z, cmap=icecm, vmin=0, vmax=1, zorder=2, **kw)
    contP = axp.contour(XP, YP, f0["IcePhase"], levels=[args.phi], colors=ICE_EDGE, linewidths=0.5, zorder=4)
    contZ = axz.contour(XZ, YZ, f0["IcePhase"][zs], levels=[args.phi], colors=ICE_EDGE, linewidths=1.1, zorder=4)

    axp.add_patch(plt.Rectangle((ext_z[0], ext_z[2]), ext_z[1] - ext_z[0], ext_z[3] - ext_z[2],
                                fill=False, ec=ACCENT, lw=1.8, ls=(0, (4, 3)), zorder=5))
    axp.text(ext_z[0], ext_z[3] + 12, "zoom", color=ACCENT, fontsize=13, fontweight="bold",
             va="bottom", ha="left", zorder=6, bbox=dict(fc="white", ec="none", alpha=0.8, pad=1.5))
    for ax_, ex in ((axp, ext_p), (axz, ext_z)):
        ax_.set_xlim(ex[0], ex[1]); ax_.set_ylim(ex[2], ex[3]); ax_.set_aspect("equal")
        ax_.set_xticks([]); ax_.set_yticks([])
        for s_ in ax_.spines.values():
            s_.set_edgecolor("#9aa0a8"); s_.set_linewidth(0.8)

    def scalebar(ax, ex, length_um, label):
        x1 = ex[1] - 0.04 * (ex[1] - ex[0]); x0 = x1 - length_um
        y = ex[2] + 0.05 * (ex[3] - ex[2])
        ax.add_patch(plt.Rectangle((x0 - 0.02 * (ex[1] - ex[0]), ex[2] + 0.01 * (ex[3] - ex[2])),
                                   (x1 - x0) + 0.04 * (ex[1] - ex[0]), 0.085 * (ex[3] - ex[2]),
                                   fc="white", ec="none", alpha=0.88, zorder=7))
        ax.plot([x0, x1], [y, y], color=INK, lw=3.2, solid_capstyle="butt", zorder=8)
        ax.text(0.5 * (x0 + x1), y + 0.015 * (ex[3] - ex[2]), label, ha="center", va="bottom",
                fontsize=13, color=INK, zorder=8)
    scalebar(axp, ext_p, 500, r"500 $\mu$m")
    zl = 50 if (ext_z[1] - ext_z[0]) < 400 else 100
    scalebar(axz, ext_z, zl, rf"{zl} $\mu$m")

    (waistZ,) = axz.plot([], [], lw=3.0, color=ACCENT, zorder=7, solid_capstyle="butt",
                         marker="_", markersize=14, markeredgewidth=3.0)
    wlabel = axz.text(0, 0, "", color=ACCENT, fontsize=15, fontweight="bold", ha="left",
                      va="center", zorder=9, bbox=dict(fc="white", ec="none", alpha=0.85, pad=2))
    stamp = axp.text(0.02, 0.975, "", transform=axp.transAxes, fontsize=22, fontweight="bold",
                     color=INK, va="top", ha="left", zorder=9,
                     bbox=dict(fc="white", ec="none", alpha=0.85, pad=3))

    # colour bar under the packing, positioned from the drawn axes
    cax = ax_at(0.55, 0.68, 7.0, 0.22)
    cb = fig.colorbar(vP, cax=cax, orientation="horizontal", ticks=sigma_ticks(norm, min_gap=0.085))
    cb.set_label(r"vapor supersaturation  $\sigma=\rho_v/\rho_{vs}-1$   "
                 rf"[$\times 10^{{{int(round(np.log10(1/S)))}}}$]", fontsize=14)
    cb.ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v_, _p: f"{v_:.3g}"))
    cb.ax.tick_params(labelsize=12)
    cb.ax.text(0.0, 1.9, "undersaturated: ice sublimates", transform=cb.ax.transAxes,
               fontsize=12, color="#555a61", ha="left", va="bottom")
    cb.ax.text(1.0, 1.9, "supersaturated: vapor deposits", transform=cb.ax.transAxes,
               fontsize=12, color="#555a61", ha="right", va="bottom")

    # curve: k_eff / k_eff,r against time, log axis in days
    t_open = t_s[0]
    mk = kt >= t_open
    kd = kt[mk] / 86400.0
    axc.plot(kd, kn[mk], lw=1.8, color=GREY, zorder=1)
    (trace,) = axc.plot([], [], lw=3.0, color=TRACE, zorder=2)
    (dot,) = axc.plot([], [], "o", ms=10, mfc=ACCENT, mec="white", mew=1.5, zorder=3)
    axc.axvline(THETA_REF * tau / 86400.0, color="#9aa0a8", lw=1.0, ls=":", zorder=0)
    axc.text(THETA_REF * tau / 86400.0, 0.02, r" $\theta$ = 30", transform=axc.get_xaxis_transform(),
             fontsize=12, color="#555a61", ha="left", va="bottom")
    axc.set_xscale("log")
    axc.set_xlim(kd[0] * 0.8, 31.0)
    lo, hi = float(kn[mk].min()), float(kn[mk].max())
    axc.set_ylim(lo - 0.06 * (hi - lo), hi + 0.08 * (hi - lo))
    axc.set_xlabel("time  [days]", fontsize=14)
    axc.set_ylabel(r"$k_{\mathrm{eff}}\,/\,k_{\mathrm{eff,r}}$", fontsize=15)
    axc.tick_params(labelsize=12)
    axc.grid(alpha=0.25, lw=0.6, which="both")
    for s_ in ("top", "right"):
        axc.spines[s_].set_visible(False)

    readout = fig.text(12.15 / W, 8.18 / H, "", fontsize=16, color=INK, va="top", ha="left",
                       linespacing=1.55)

    def draw(k):
        nonlocal contP, contZ
        fk = read_vts(files[k], want=WANT)[0]
        phi = fk["IcePhase"]
        sk = S * pplib.supersaturation(fk["VaporDensity"], fk["Temperature"])
        iP.set_data(phi); vP.set_data(sk)
        iZ.set_data(phi[zs]); vZ.set_data(sk[zs])
        contP.remove(); contZ.remove()
        contP = axp.contour(XP, YP, phi, levels=[args.phi], colors=ICE_EDGE, linewidths=0.5, zorder=4)
        contZ = axz.contour(XZ, YZ, phi[zs], levels=[args.phi], colors=ICE_EDGE, linewidths=1.1, zorder=4)

        if ends[k] is not None:
            (a, b) = ends[k]
            waistZ.set_data([a[0] * um, b[0] * um], [a[1] * um, b[1] * um])
            mid = 0.5 * (a + b) * um
            wlabel.set_position((mid[0] + 0.06 * (ext_z[1] - ext_z[0]), mid[1]))
            wlabel.set_text(rf"{widths[k]:.1f} $\mu$m")
        tk = t_s[k]
        stamp.set_text(f"t = {tk / 86400.0:.2f} d" if tk >= 86400 else f"t = {fmt_time(tk)}")
        kcur = float(np.interp(tk, kt, kn)); scur = float(np.interp(tk, ssa_t, ssa_v)) / ssa_open
        sel_ = (kt >= t_open) & (kt <= tk)
        trace.set_data(kt[sel_] / 86400.0, kn[sel_])
        dot.set_data([tk / 86400.0], [kcur])
        readout.set_text(f"neck width   {widths[k]:.1f} $\\mu$m\n"
                         f"SSA / SSA$_0$   {scur:.2f}\n"
                         f"$k_{{\\mathrm{{eff}}}}$ / $k_{{\\mathrm{{eff,r}}}}$   {kcur:.2f}")

    if args.frame_png is not None:
        k = int(np.argmin([abs(step_of(f) - args.frame_png) for f in files]))
        draw(k)
        out = args.out or (run / f"aggregate_frame_{step_of(files[k]):05d}.png")
        fig.savefig(out, dpi=args.dpi)
        print(f"preview -> {out}")
        return

    out = args.out or (run / "plots" / "aggregate_movie.mp4")
    out.parent.mkdir(parents=True, exist_ok=True)
    writer = FFMpegWriter(fps=args.fps, bitrate=-1, codec="libx264",
                          extra_args=["-pix_fmt", "yuv420p", "-crf", "20", "-movflags", "+faststart"],
                          metadata={"title": f"{run.name} aggregate"})
    print(f"rendering {len(files)} frames -> {out}")
    with writer.saving(fig, str(out), dpi=args.dpi):
        for k in range(len(files)):
            draw(k)
            writer.grab_frame()
            if k % 10 == 0:
                print(f"  frame {k}/{len(files)}", flush=True)
    print(f"movie -> {out}")


if __name__ == "__main__":
    main()
