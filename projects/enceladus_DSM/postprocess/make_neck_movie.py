#!/usr/bin/env python3
"""make_neck_movie.py — animate the NECK of a two-grain sintering run.

One frame per snapshot, two synchronised panels:

  top     the ice body around the neck (phi field + the phi = 0.5 contour),
          with the measured waist drawn ON the frame as a capped segment and
          labelled in um. For an axisymmetric r-z run the field is mirrored
          about the axis (r = 0) so the frame shows the physical pair, and the
          drawn segment is the full width 2r, matching what is reported.
  bottom  neck width vs time, the whole curve in grey with the elapsed part
          picked out in colour and a dot at the current frame.

The measurement is neck_width.py's, not a re-invention: same minimum-cross-
section definition, same grain-centre peak split, same sub-grid parabola
refinement (imported from it). So the number burned into the frame is the
number in neck_width.csv, and a frame can be used as evidence for the curve.

Two passes over vtkOut/: one to build the w(t) series (the bottom panel needs
the whole curve from frame 1), one to render. Snapshots are small -- these are
the control-point .vts that plot_fields.py writes, ~0.1 s each -- so the second
read costs less than caching the fields would.

NO EXPERIMENTAL OVERLAY. Molaro's neck widths need the anchoring convention
fit_neck_growth.py applies (their t = 0 is not ours); dropping the raw points
onto a movie axis would compare two different clocks. Use
plot_neck_vs_molaro.py for that comparison.

Usage:
    python make_neck_movie.py <run_dir> [--out FILE.mp4] [--fps 12]
        [--stride N] [--dpi 150] [--frame-png STEP]

Auto-detects -axisym from the run's .opts; force with --axisym/--no-axisym.
"""

import argparse
import glob
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
import cmocean

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pplib
from pplib import read_vts, step_of, step_times
from neck_width import _cross, refine_min

ICE_EDGE = "#1b3a5c"
TRACE = "#3d74d9"
ACCENT = "#d1495b"


def chord_bounds(col, y, level):
    """(y_lo, y_hi) of the outermost `level`-crossings of phi along a column.

    neck_width.chord_width returns only the span; the endpoints themselves are
    needed here to draw the waist on the frame. The crossings come from its
    _cross, i.e. interpolated in LOGIT, not linearly in phi -- on this .vts
    grid (dy = 3.5 eps) a linear crossing carries a sub-cell-phase error that
    shows up as a spurious ripple riding on the neck curve, and a frame
    annotated with it would not match neck_width.csv.
    """
    above = col >= level
    if not above.any():
        return None
    idx = np.flatnonzero(above)
    lo_i, hi_i = idx[0], idx[-1]
    y_lo = (y[lo_i] if lo_i == 0 else
            _cross(y[lo_i - 1], y[lo_i], col[lo_i - 1], col[lo_i], level))
    y_hi = (y[hi_i] if hi_i == len(y) - 1 else
            _cross(y[hi_i], y[hi_i + 1], col[hi_i], col[hi_i + 1], level))
    return y_lo, y_hi


def widths(phi, y, level, axisym):
    """w(x): the ice chord on every column, doubled for an axisymmetric run."""
    w = np.zeros(phi.shape[1])
    for j in range(phi.shape[1]):
        b = chord_bounds(phi[:, j], y, level)
        if b is not None:
            w[j] = b[1] - b[0]
    return 2.0 * w if axisym else w


def grain_centers(w):
    """The two grain-centre columns: outermost prominent local maxima of w(x).

    neck_width.py's rule verbatim -- a column beating its +-5 neighbourhood and
    taller than half the peak. Not the global minimum: the grain TIPS have the
    smallest nonzero w and would hijack the split.
    """
    thr = 0.5 * w.max()
    peaks = [j for j in range(5, len(w) - 5)
             if w[j] >= thr and w[j] == w[j - 5:j + 6].max()]
    if len(peaks) < 2:
        sys.exit("could not locate two grain-center peaks in w(x)")
    return peaks[0], peaks[-1]


def measure(phi, x, y, level, axisym, centers):
    """(neck width, x_neck) for one snapshot, plus the refreshed centers."""
    w = widths(phi, y, level, axisym)
    if centers is None:
        centers = grain_centers(w)
    lo, hi = centers
    interior = np.arange(lo + 1, hi)
    interior = interior[w[interior] > 0]
    if len(interior) == 0:
        return 0.0, np.nan, centers
    jn = interior[np.argmin(w[interior])]
    neck, xneck = refine_min(w, x, jn)
    return neck, xneck, centers


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--out", type=Path, default=None,
                    help="output mp4 (default: <run_dir>/plots/neck_movie.mp4)")
    ap.add_argument("--phi", type=float, default=0.5, help="contour level")
    ap.add_argument("--fps", type=int, default=12)
    ap.add_argument("--stride", type=int, default=1, help="use every Nth snapshot")
    ap.add_argument("--dpi", type=int, default=150)
    ap.add_argument("--zoom", type=float, default=1.6,
                    help="half-window of the zoom panel along x, in units of "
                         "the largest neck width reached (default 1.6)")
    ap.add_argument("--frame-png", type=int, default=None,
                    help="render only this step to a PNG (preview) and exit")
    ax_ = ap.add_mutually_exclusive_group()
    ax_.add_argument("--axisym", dest="axisym", action="store_true", default=None)
    ax_.add_argument("--no-axisym", dest="axisym", action="store_false")
    args = ap.parse_args()

    files = sorted(glob.glob(str(args.run_dir / "vtkOut" / "solV_*.vts")),
                   key=step_of)[:: args.stride]
    if not files:
        sys.exit(f"no solV_*.vts under {args.run_dir}/vtkOut")

    opts = pplib.read_opts(str(args.run_dir))
    if args.axisym is None:
        args.axisym = str(opts.get("-axisym", "0")).strip() in ("1", "true", "True")
    tmap = step_times(args.run_dir)

    # ---- pass 1: the w(t) series -----------------------------------------
    print(f"measuring {len(files)} snapshots "
          f"({'axisymmetric' if args.axisym else 'planar'}) ...")
    centers = None
    ts, ws, xs = [], [], []
    for i, fn in enumerate(files):
        f, X, Y = read_vts(fn, want=("IcePhase",))
        w, xn, centers = measure(f["IcePhase"], X[0, :], Y[:, 0],
                                 args.phi, args.axisym, centers)
        ts.append(tmap.get(step_of(fn), np.nan)); ws.append(w); xs.append(xn)
        if i % 25 == 0:
            print(f"  {i}/{len(files)}", flush=True)
    ts = np.asarray(ts); ws = np.asarray(ws); xs = np.asarray(xs)
    t_min = ts / 60.0
    w_um = ws * 1e6
    print(f"neck width: {w_um[0]:.2f} -> {w_um[-1]:.2f} um "
          f"({100*(w_um[-1]/w_um[0]-1):+.1f}%) over {t_min[-1]:.0f} min")

    # ---- view windows, both fixed over the run --------------------------
    # Fixed, not tracking x_neck: a window centred on the moving neck would
    # make the neck look stationary while the grains slid past it, hiding the
    # plane's migration toward the smaller grain -- part of what this is for.
    _, X, Y = read_vts(files[0], want=("IcePhase",))
    x1d, y1d = X[0, :], Y[:, 0]
    yc = 0.0 if args.axisym else 0.5 * (y1d[0] + y1d[-1])

    def span(lo, hi, arr):
        """Index range of `arr` covering [lo, hi], never narrower than 8."""
        sel = np.flatnonzero((arr >= lo) & (arr <= hi))
        if len(sel) < 8:
            return 0, len(arr)
        return sel[0], sel[-1] + 1

    # Pair panel: the ice bounding box at t = 0, padded. Taken from the field
    # rather than -ice_grain_cx/-ice_grain_R so it holds for any IC.
    f0 = read_vts(files[0], want=("IcePhase",))[0]["IcePhase"]
    ice = f0 >= args.phi
    xi, yi = x1d[ice.any(axis=0)], y1d[ice.any(axis=1)]
    padx = 0.04 * (xi[-1] - xi[0])
    r_max = max(abs(yi[-1] - yc), abs(yi[0] - yc))
    Pj = span(xi[0] - padx, xi[-1] + padx, x1d)
    Pi = span(yc - 1.10 * r_max if not args.axisym else 0.0,
              yc + 1.10 * r_max, y1d)

    # Zoom panel: a few neck widths across, centred on the mean neck plane.
    xc = float(np.nanmean(xs))
    Zj = span(xc - args.zoom * ws.max(), xc + args.zoom * ws.max(), x1d)
    Zi = span(yc - 1.05 * ws.max() if not args.axisym else 0.0,
              yc + 1.05 * ws.max(), y1d)

    def view(phi, win_i, win_j):
        """Crop to a window; mirror about the axis for an axisym run."""
        p = phi[win_i[0]:win_i[1], win_j[0]:win_j[1]]
        yy = y1d[win_i[0]:win_i[1]]
        if args.axisym:
            p = np.vstack([p[:0:-1, :], p])
            yy = np.concatenate([-yy[:0:-1], yy])
        return p, x1d[win_j[0]:win_j[1]] * 1e6, yy * 1e6

    # ---- figure ----------------------------------------------------------
    pP, xP, yP = view(f0, Pi, Pj)
    pZ, xZ, yZ = view(f0, Zi, Zj)
    XP, YP = np.meshgrid(xP, yP)
    XZ, YZ = np.meshgrid(xZ, yZ)

    # Landscape, because the output is a video: the zoom fills the left half,
    # the pair-scale context and the curve stack on the right. Stacking all
    # three vertically instead left both image panels (aspect ~1.5) squeezed
    # by their slot height, wasting a third of the frame width.
    fig = plt.figure(figsize=(12.4, 7.0))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.3, 1.0],
                          height_ratios=[1.0, 1.0],
                          left=0.06, right=0.975, top=0.86, bottom=0.09,
                          wspace=0.22, hspace=0.45)
    axz = fig.add_subplot(gs[:, 0])
    axp = fig.add_subplot(gs[0, 1])
    axc = fig.add_subplot(gs[1, 1])

    def setup(ax, XX, YY, p):
        ax.set_aspect("equal")
        m = ax.pcolormesh(XX, YY, p, cmap=cmocean.cm.ice, vmin=0.0, vmax=1.0,
                          shading="gouraud", rasterized=True)
        ax.set_xlim(XX.min(), XX.max()); ax.set_ylim(YY.min(), YY.max())
        ax.set_xlabel(r"$x$ [$\mu$m]")
        ax.set_ylabel(r"$r$ [$\mu$m]" if args.axisym else r"$y$ [$\mu$m]")
        if args.axisym:
            ax.axhline(0.0, lw=0.7, ls=(0, (6, 4)), color="#ffffff", alpha=0.5)
        return m

    mP = setup(axp, XP, YP, pP)
    mZ = setup(axz, XZ, YZ, pZ)
    contP = axp.contour(XP, YP, pP, levels=[args.phi], colors=ICE_EDGE, linewidths=0.8)
    contZ = axz.contour(XZ, YZ, pZ, levels=[args.phi], colors=ICE_EDGE, linewidths=1.1)

    # the zoom window, drawn on the pair panel
    axp.add_patch(plt.Rectangle((XZ.min(), YZ.min()),
                                XZ.max() - XZ.min(), YZ.max() - YZ.min(),
                                fill=False, ec=ACCENT, lw=1.0, ls=(0, (4, 3)),
                                zorder=5))
    axp.set_title("grain pair", fontsize=10, color="#555a61")

    def waist_artist(ax, lw, ms):
        (ln,) = ax.plot([], [], lw=lw, color=ACCENT, zorder=6,
                        solid_capstyle="butt", marker="_", markersize=ms,
                        markeredgewidth=lw)
        return ln
    waistP = waist_artist(axp, 1.6, 7)
    waistZ = waist_artist(axz, 2.2, 13)
    wlabel = axz.text(0.0, 0.0, "", color=ACCENT, fontsize=12, fontweight="bold",
                      ha="left", va="center", zorder=7,
                      bbox=dict(fc="white", ec="none", alpha=0.78, pad=1.8))

    # bottom: the curve, grey in full, coloured up to the current frame
    good = np.isfinite(t_min) & np.isfinite(w_um)
    axc.plot(t_min[good], w_um[good], lw=1.4, color="#c4c8ce", zorder=1)
    (trace,) = axc.plot([], [], lw=2.2, color=TRACE, zorder=2)
    (dot,) = axc.plot([], [], "o", ms=7, mfc=ACCENT, mec="white", mew=1.2, zorder=3)
    axc.set_xlim(0.0, np.nanmax(t_min) * 1.02)
    pad = 0.08 * (np.nanmax(w_um) - np.nanmin(w_um) + 1e-12)
    axc.set_ylim(np.nanmin(w_um) - pad, np.nanmax(w_um) + pad)
    axc.set_xlabel("time [min]")
    axc.set_ylabel(r"neck width $2r$ [$\mu$m]" if args.axisym
                   else r"neck width [$\mu$m]")
    axc.grid(alpha=0.25, lw=0.6)
    for s_ in ("top", "right"):
        axc.spines[s_].set_visible(False)

    bits = []
    if "-temp" in opts:
        bits.append(f"T = {float(opts['-temp']):g} \u00b0C")
    if "-humidity" in opts:
        bits.append(f"RH = {float(opts['-humidity'])*100:.3f}%")
    subtitle = "    ".join(bits)
    fig.suptitle(f"{args.run_dir.name}\n{subtitle}", fontsize=9, y=0.975)
    title = axz.set_title("", fontsize=13)

    dx_lab = 0.025 * (XZ.max() - XZ.min())

    def draw(k):
        nonlocal contP, contZ
        phi = read_vts(files[k], want=("IcePhase",))[0]["IcePhase"]
        pP_, _, _ = view(phi, Pi, Pj)
        pZ_, _, _ = view(phi, Zi, Zj)
        mP.set_array(pP_.ravel()); mZ.set_array(pZ_.ravel())
        contP.remove(); contZ.remove()
        contP = axp.contour(XP, YP, pP_, levels=[args.phi], colors=ICE_EDGE,
                            linewidths=0.8)
        contZ = axz.contour(XZ, YZ, pZ_, levels=[args.phi], colors=ICE_EDGE,
                            linewidths=1.1)

        xn = xs[k] * 1e6
        if args.axisym:
            lo, hi = -0.5 * w_um[k], 0.5 * w_um[k]
        else:
            b = chord_bounds(phi[:, int(np.argmin(np.abs(x1d - xs[k])))],
                             y1d, args.phi)
            lo, hi = ((b[0] * 1e6, b[1] * 1e6) if b else
                      (yc * 1e6 - 0.5 * w_um[k], yc * 1e6 + 0.5 * w_um[k]))
        waistP.set_data([xn, xn], [lo, hi])
        waistZ.set_data([xn, xn], [lo, hi])
        wlabel.set_position((xn + dx_lab, 0.5 * (lo + hi)))
        wlabel.set_text(f"{w_um[k]:.2f} " + r"$\mu$m")

        trace.set_data(t_min[: k + 1], w_um[: k + 1])
        dot.set_data([t_min[k]], [w_um[k]])
        title.set_text(f"t = {t_min[k]:6.1f} min      "
                       f"neck {w_um[k]:.2f} " + r"$\mu$m"
                       + f"   ({100*(w_um[k]/w_um[0]-1):+.1f}%)")

    if args.frame_png is not None:
        k = int(np.argmin([abs(step_of(f) - args.frame_png) for f in files]))
        draw(k)
        out = args.run_dir / f"neck_frame_{step_of(files[k]):05d}.png"
        fig.savefig(out, dpi=args.dpi)
        print(f"preview -> {out}")
        return

    out = args.out or (args.run_dir / "plots" / "neck_movie.mp4")
    out.parent.mkdir(parents=True, exist_ok=True)
    writer = FFMpegWriter(fps=args.fps, bitrate=-1,
                          metadata={"title": f"{args.run_dir.name} neck"})
    print(f"rendering {len(files)} frames -> {out}")
    with writer.saving(fig, str(out), dpi=args.dpi):
        for k in range(len(files)):
            draw(k)
            writer.grab_frame()
            if k % 25 == 0:
                print(f"  frame {k}/{len(files)}", flush=True)
    print(f"movie -> {out}")


if __name__ == "__main__":
    main()
