#!/usr/bin/env python3
"""make_neck_movie.py — animate the NECK of a two-grain sintering run.

One frame per snapshot, three synchronised panels:

  left    the neck close up. The vapour field is the background, the ice is
          painted opaque on top, and the measured waist is drawn ON the frame
          as a capped segment, labelled in um.
  top r.  the whole grain pair at the same instant, for scale, with the zoom
          window marked.
  bot r.  neck width vs time, the whole curve in grey with the elapsed part
          picked out in colour and a dot at the current frame.

For an axisymmetric r-z run both image panels are mirrored about the axis
(r = 0) so the frame shows the physical pair, and the drawn segment is the
full width 2r, matching what is reported.

THE BACKGROUND IS SUPERSATURATION, sigma = rho_v/rho_vs(T) - 1, not vapour
density: rho_v alone varies in its 4th significant figure and its structure is
invisible without a scale so tight it means nothing. sigma is the quantity the
interface actually responds to -- its sign is the sign of the local growth
rate -- and it makes the movie explain itself: vapour pools in the concave
neck (sigma highest, the ice there is the least-volatile surface in the frame)
and drains toward the undersaturated wall, which is why the waist fills in.

Colours are cmocean: 'ice' for the phase field, 'balance' for sigma. balance
is diverging, so its pale middle is re-sampled onto sigma = 0 (see
centered_cmap) -- blue is undersaturated, red supersaturated, pale is
saturation over flat ice.

Molaro's Fig. 11 points are overlaid on the curve, SHIFTED ONTO THE RUN'S
CLOCK by plot_neck_vs_molaro.py's anchor (t = 0 for them is the moment the
model passes their first measured width, 32.81 um). The experiment moves, not
the model, so the time axis stays the clock in the frame titles. --no-data
drops them.

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
        [--stride N] [--dpi 150] [--frame-png STEP] [--no-vapor]
        [--no-data] [--cmap NAME]

Auto-detects -axisym from the run's .opts; force with --axisym/--no-axisym.
--no-vapor gives the plain ice/air panels; it is also the automatic fallback
for snapshots that carry no VaporDensity/Temperature.
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
from matplotlib.colors import AsinhNorm, ListedColormap
import cmocean

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pplib
from pplib import read_vts, step_of, step_times
from neck_width import _cross, refine_min

ICE_EDGE = "#12263a"
TRACE = "#3d74d9"
ACCENT = "#d1495b"
C_DATA = "#d1495b"         # plot_neck_vs_molaro.py's colour for the experiment
SIGMA_SCALE = 1e4          # sigma is O(1e-4); supersat_probes.py's convention
DEFAULT_ANCHOR_UM = 32.81  # Molaro's first measured neck WIDTH


def ice_alpha_cmap(n=256):
    """cmocean 'ice' over phi, fully TRANSPARENT below phi = 0.5.

    The ice is painted ON TOP of the vapour field as one opaque layer rather
    than the vapour being masked to the pore: every pixel is covered by the
    base, so no pixel can be left unowned at the phi = 0.5 boundary (the
    dropped-pixel failure of thresholding into two disjoint regions).

    phi maps DIRECTLY through the map, so the visible ice (phi in [0.5, 1])
    uses its upper half: grains read light with a mid-tone rim at the
    interface, not the near-black ice(0) rim a [0.5,1]->[0,1] remap would give.
    Gouraud shading interpolates the RGBA, so across the interface only alpha
    ramps -- an anti-aliased edge with no dark fringe.
    """
    phi = np.linspace(0.0, 1.0, n)
    lut = cmocean.cm.ice(phi)
    lut[:, 3] = np.where(phi < 0.5, 0.0, 1.0)
    return ListedColormap(lut)


def centered_cmap(base, norm, n=512):
    """Re-map a DIVERGING colormap so its midpoint lands on sigma = 0.

    balance is diverging, and a diverging map makes a promise: the pale middle
    is the neutral value. A norm maps vmin->0 and vmax->1 regardless, so with
    the asymmetric ranges here (-28.5 .. +0.29) the pale middle would land at
    sigma = -14, i.e. at nothing, and the whole undersaturated field would read
    as "positive". Re-sampling the map so its 0.5 sits at norm(0) keeps the
    promise AND the asinh dynamic range: blue is undersaturated, red is
    supersaturated, pale is saturation over flat ice.

    When the data never reaches zero -- the high-D_v arms are undersaturated
    everywhere -- only the blue half is used, with the pale end as the
    saturation the pore never gets to. That is the honest picture, and it is
    why this is not simply symmetric limits about zero: those would throw away
    most of the range to represent a sign that never occurs.
    """
    if norm.vmin < 0.0 < norm.vmax:
        p0 = float(np.clip(norm(0.0), 0.0, 1.0))
    else:
        p0 = 1.0 if norm.vmax <= 0.0 else 0.0
    p = np.linspace(0.0, 1.0, n)
    if p0 <= 0.0:
        q = 0.5 + 0.5 * p
    elif p0 >= 1.0:
        q = 0.5 * p
    else:
        q = np.where(p <= p0, 0.5 * p / p0, 0.5 + 0.5 * (p - p0) / (1.0 - p0))
    return ListedColormap(base(q))


def read_experiment(path):
    """(t_s, width_m, err+_m, err-_m) from a Molaro Fig. 11 validation CSV."""
    import csv as _csv
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            f = next(_csv.reader([line]))
            rows.append((float(f[0]) * 60.0, float(f[1]) * 1e-6,
                         float(f[2]) * 1e-6, float(f[3]) * 1e-6))
    a = np.asarray(rows)
    return a[:, 0], a[:, 1], a[:, 2], a[:, 3]


def anchor_time(t, w, target):
    """First time the model curve reaches `target` width, interpolated.

    plot_neck_vs_molaro.py's convention, and it has to be: our t = 0 and
    theirs are not the same instant. We start from a chosen initial neck;
    their record starts at 32.81 um at an unknown time after contact. Putting
    both on a common clock means defining t = 0 as the moment each passes the
    SAME physical width.
    """
    for (t0, w0), (t1, w1) in zip(zip(t, w), zip(t[1:], w[1:])):
        if w0 <= target <= w1 and w1 > w0:
            return t0 + (target - w0) * (t1 - t0) / (w1 - w0)
    return None


def sigma_ticks(norm, min_gap=0.055):
    """Decade ticks for an asinh colourbar, thinned in NORM space.

    Two rules the obvious version gets wrong:

    * Thin by BAR POSITION, not by sigma. With the bar running -28.5 .. +0.29,
      a 4 %-of-span guard is +-1.15 in sigma, which swallows every tick between
      -1 and the top.
    * Zero is MANDATORY, not a candidate. It is the one value on the bar that
      means something physical -- saturation over flat ice, the sign change
      between net sublimation and net deposition -- and asinh spacing puts it
      close enough to its neighbours that a plain greedy pass drops it and
      keeps 0.01 instead.

    Decades go in first and 3x-decades fill the gaps, so a narrow bar stays
    readable and a wide one does not go sparse.
    """
    vmin, vmax = float(norm.vmin), float(norm.vmax)
    pos = lambda v: float(norm(v))
    kept = [vmin, vmax] + ([0.0] if vmin < 0.0 < vmax else [])
    decades = [10.0 ** k for k in range(-2, 3)]
    for scale in (1.0, 3.0):
        for d in decades:
            for v in (-scale * d, scale * d):
                if not (vmin < v < vmax):
                    continue
                if all(abs(pos(v) - pos(u)) >= min_gap for u in kept):
                    kept.append(v)
    return sorted(kept)


WANT = ("IcePhase", "VaporDensity", "Temperature")


def sigma_field(f):
    """Supersaturation x 1e4 from a snapshot's fields.

    pplib.supersaturation mirrors the solver's RhoVS_I; the local Temperature
    array is used rather than the -temp option because a run with a gradient
    has neither a single T nor a single rho_vs. Snapshots without the two
    fields are caught in main(), which drops the vapour layer.
    """
    return SIGMA_SCALE * pplib.supersaturation(f["VaporDensity"],
                                               f["Temperature"])


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
    ap.add_argument("--cmap", default="balance",
                    help="cmocean (or matplotlib) colourmap for the vapour "
                         "background (default: cmocean balance)")
    ap.add_argument("--no-center", dest="center", action="store_false",
                    help="do not re-centre the colourmap on sigma = 0; pass "
                         "this when --cmap is sequential rather than diverging")
    ap.add_argument("--no-symmetric", dest="symmetric", action="store_false",
                    help="scale the bar to the data range instead of putting "
                         "sigma = 0 at its midpoint. Recovers contrast on runs "
                         "that are undersaturated everywhere, where half the "
                         "symmetric bar covers a sign the data never reaches")
    ap.add_argument("--data", type=Path, default=None,
                    help="experimental neck-width CSV to overlay on the curve "
                         "(default: the repo's Molaro Fig. 11 T=-20 table)")
    ap.add_argument("--anchor-width", type=float, default=DEFAULT_ANCHOR_UM,
                    help="neck WIDTH [um] at which the experiment's clock is "
                         "matched to the model's (default 32.81)")
    ap.add_argument("--no-data", dest="data_on", action="store_false",
                    help="model curve only, no experimental points")
    ap.add_argument("--sat-vmin", type=float, default=None,
                    help="fix the low end of the supersaturation scale, in "
                         "units of sigma x 1e4 (default: the run's own min)")
    ap.add_argument("--sat-vmax", type=float, default=None)
    ap.add_argument("--no-vapor", dest="vapor", action="store_false",
                    help="plain ice/air panels, no vapour background")
    ax_ = ap.add_mutually_exclusive_group()
    ax_.add_argument("--axisym", dest="axisym", action="store_true", default=None)
    ax_.add_argument("--no-axisym", dest="axisym", action="store_false")
    args = ap.parse_args()

    files = sorted(glob.glob(str(args.run_dir / "vtkOut" / "solV_*.vts")),
                   key=step_of)[:: args.stride]
    if not files:
        sys.exit(f"no solV_*.vts under {args.run_dir}/vtkOut")

    data_path = args.data or (Path(__file__).resolve().parent.parent
                              / "inputs/validation/molaro2019_fig11_T-20.csv")
    if args.data_on and not Path(data_path).is_file():
        print(f"  NOTE: {data_path} not found; model curve only")
        args.data_on = False

    opts = pplib.read_opts(str(args.run_dir))
    if args.axisym is None:
        args.axisym = str(opts.get("-axisym", "0")).strip() in ("1", "true", "True")
    tmap = step_times(args.run_dir)

    # ---- geometry of the pair panel, needed before pass 1 ----------------
    # The vapour colour scale is built over the PAIR WINDOW's pore, not the
    # whole domain: the far field is pinned at the Dirichlet wall value, and a
    # scale stretched to reach it spends most of its range on a boundary layer
    # neither panel shows, flattening both panels to one colour.
    f0f, X, Y = read_vts(files[0], want=WANT)
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
    f0 = f0f["IcePhase"]
    missing = [k for k in WANT[1:] if k not in f0f]
    if args.vapor and missing:
        print(f"  WARNING: {', '.join(missing)} not in the snapshots; "
              f"drawing ice/air only")
        args.vapor = False
    ice0 = f0 >= args.phi
    xi, yi = x1d[ice0.any(axis=0)], y1d[ice0.any(axis=1)]
    padx = 0.04 * (xi[-1] - xi[0])
    r_max = max(abs(yi[-1] - yc), abs(yi[0] - yc))
    Pj = span(xi[0] - padx, xi[-1] + padx, x1d)
    Pi = span(yc - 1.10 * r_max if not args.axisym else 0.0,
              yc + 1.10 * r_max, y1d)
    Pbox = np.s_[Pi[0]:Pi[1], Pj[0]:Pj[1]]

    # ---- pass 1: the w(t) series, and the vapour range -------------------
    print(f"measuring {len(files)} snapshots "
          f"({'axisymmetric' if args.axisym else 'planar'}) ...")
    centers = None
    ts, ws, xs = [], [], []
    smin, smax = np.inf, -np.inf
    for i, fn in enumerate(files):
        f, X, Y = read_vts(fn, want=WANT)
        phi = f["IcePhase"]
        w, xn, centers = measure(phi, X[0, :], Y[:, 0],
                                 args.phi, args.axisym, centers)
        ts.append(tmap.get(step_of(fn), np.nan)); ws.append(w); xs.append(xn)
        if args.vapor:
            pore = phi[Pbox] < args.phi
            if pore.any():
                sig = sigma_field(f)[Pbox][pore]
                smin = min(smin, float(sig.min()))
                smax = max(smax, float(sig.max()))
        if i % 25 == 0:
            print(f"  {i}/{len(files)}", flush=True)
    ts = np.asarray(ts); ws = np.asarray(ws); xs = np.asarray(xs)
    t_min = ts / 60.0
    w_um = ws * 1e6
    print(f"neck width: {w_um[0]:.2f} -> {w_um[-1]:.2f} um "
          f"({100*(w_um[-1]/w_um[0]-1):+.1f}%) over {t_min[-1]:.0f} min")

    # ---- zoom window -----------------------------------------------------
    # Both windows are FIXED over the run, not tracking x_neck: a window
    # centred on the moving neck would hold it still while the grains slid
    # past it, hiding the plane's migration toward the smaller grain -- part
    # of what this is for.
    xc = float(np.nanmean(xs))
    Zj = span(xc - args.zoom * ws.max(), xc + args.zoom * ws.max(), x1d)
    Zi = span(yc - 1.05 * ws.max() if not args.axisym else 0.0,
              yc + 1.05 * ws.max(), y1d)

    def view(fld, win_i, win_j):
        """Crop to a window; mirror about the axis for an axisym run."""
        p = fld[win_i[0]:win_i[1], win_j[0]:win_j[1]]
        yy = y1d[win_i[0]:win_i[1]]
        if args.axisym:
            p = np.vstack([p[:0:-1, :], p])
            yy = np.concatenate([-yy[:0:-1], yy])
        return p, x1d[win_j[0]:win_j[1]] * 1e6, yy * 1e6

    # ---- vapour colour scale --------------------------------------------
    # ONE scale for both panels, asinh-spaced. The panels' ranges differ by
    # roughly a decade (the pair panel reaches the drained field between the
    # grains, the zoom sees only the neck's own boundary layer), so a linear
    # scale covering the pair panel leaves the zoom panel a single flat
    # colour. asinh is log-like over the decades and linear through zero, so
    # it fits both without a false centre -- and unlike a symmetric diverging
    # norm it does not need the data to straddle zero, which for the higher-D_v
    # arms it never does (their pore is undersaturated everywhere).
    if args.vapor:
        vmin = args.sat_vmin if args.sat_vmin is not None else smin
        vmax = args.sat_vmax if args.sat_vmax is not None else smax
        if args.symmetric and args.sat_vmin is None and args.sat_vmax is None:
            # sigma = 0 at the MIDPOINT OF THE BAR, so the colourbar reads as
            # a signed axis: left half sublimation-driving, right half
            # deposition-driving, the tick in the middle the sign change.
            # Re-sampling the colormap (centered_cmap) puts balance's pale
            # middle on sigma = 0 but leaves it at 78 % along the bar; only
            # symmetric limits move the POSITION.
            v = max(abs(smin), abs(smax))
            vmin, vmax = -v, v
        if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
            print("  WARNING: no usable vapour range; drawing ice/air only")
            args.vapor = False
        else:
            # linear_width off the HALF-span: the symmetric bar is twice as
            # wide as the data, and scaling it off the full span would double
            # the linear region and squeeze the decades the neck lives in.
            norm = AsinhNorm(
                linear_width=max(max(abs(vmin), abs(vmax)) / 300.0, 1e-12),
                vmin=vmin, vmax=vmax)
            base = getattr(cmocean.cm, args.cmap, None)
            if base is None:
                base = plt.get_cmap(args.cmap)
            vapcm = centered_cmap(base, norm) if args.center else base
            print(f"  supersaturation: data {smin:+.3g} .. {smax:+.3g}, "
                  f"bar {vmin:+.3g} .. {vmax:+.3g}  (sigma x 1e4, "
                  f"0 at {float(norm(0.0)):.3f} of the bar)")

    # ---- figure ----------------------------------------------------------
    pP, xP, yP = view(f0, Pi, Pj)
    pZ, xZ, yZ = view(f0, Zi, Zj)
    XP, YP = np.meshgrid(xP, yP)
    XZ, YZ = np.meshgrid(xZ, yZ)
    s0 = sigma_field(f0f) if args.vapor else None

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

    icecm = ice_alpha_cmap()

    def setup(ax, XX, YY, p, win):
        """Vapour as a full-window base layer, ice painted opaque on top."""
        ax.set_aspect("equal")
        vap = None
        if args.vapor:
            vap = ax.pcolormesh(XX, YY, view(s0, *win)[0], cmap=vapcm,
                                norm=norm, shading="gouraud", rasterized=True)
        m = ax.pcolormesh(XX, YY, p, cmap=icecm, vmin=0.0, vmax=1.0,
                          shading="gouraud", rasterized=True, zorder=2)
        ax.set_xlim(XX.min(), XX.max()); ax.set_ylim(YY.min(), YY.max())
        ax.set_xlabel(r"$x$ [$\mu$m]")
        ax.set_ylabel(r"$r$ [$\mu$m]" if args.axisym else r"$y$ [$\mu$m]")
        if args.axisym:
            ax.axhline(0.0, lw=0.7, ls=(0, (6, 4)), color=ICE_EDGE,
                       alpha=0.35, zorder=3)
        return m, vap

    mP, vP = setup(axp, XP, YP, pP, (Pi, Pj))
    mZ, vZ = setup(axz, XZ, YZ, pZ, (Zi, Zj))
    contP = axp.contour(XP, YP, pP, levels=[args.phi], colors=ICE_EDGE,
                        linewidths=0.8, zorder=4)
    contZ = axz.contour(XZ, YZ, pZ, levels=[args.phi], colors=ICE_EDGE,
                        linewidths=1.1, zorder=4)

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

    # Molaro's points, SHIFTED ONTO THE RUN'S CLOCK -- the experiment moves,
    # not the model, so the time axis stays the clock in the frame titles and
    # the dot always sits under the t the panel above reports. The shift is
    # plot_neck_vs_molaro.py's anchor: t = 0 for them is the moment the model
    # passes their first measured width, because the two t = 0's are not the
    # same instant (we start from a chosen neck, their record starts at 32.81
    # um an unknown time after contact). Raw, unshifted points would compare
    # two different clocks.
    ylo, yhi = np.nanmin(w_um), np.nanmax(w_um)
    if args.data_on:
        try:
            td, wd, ep, em = read_experiment(data_path)
            t_star = anchor_time(ts, ws, args.anchor_width * 1e-6)
            if t_star is None:
                print(f"  NOTE: the run never crosses the "
                      f"{args.anchor_width:.2f} um anchor "
                      f"({w_um.min():.2f}-{w_um.max():.2f} um); "
                      f"plotting the experiment on its own clock instead")
                t_star = 0.0
            axc.errorbar((t_star + td) / 60.0, wd * 1e6,
                         yerr=np.vstack([em, ep]) * 1e6,
                         fmt="o", ms=4.5, lw=0.0, elinewidth=1.0, capsize=2.5,
                         color=C_DATA, mfc="white", mew=1.3, zorder=4,
                         label=(f"Molaro et al. 2019\n"
                                f"(clock anchored at "
                                f"{args.anchor_width:.2f} $\\mu$m)"))
            axc.legend(loc="lower right", fontsize=7.5, frameon=False,
                       handletextpad=0.6, borderpad=0.2)
            print(f"  experiment anchored at {args.anchor_width:.2f} um: "
                  f"t* = {t_star/60.0:.2f} min")
            ylo = min(ylo, float((wd - em).min() * 1e6))
            yhi = max(yhi, float((wd + ep).max() * 1e6))
        except (OSError, ValueError, IndexError) as e:
            print(f"  WARNING: no experimental overlay ({e})")

    axc.set_xlim(0.0, np.nanmax(t_min) * 1.02)
    pad = 0.08 * (yhi - ylo + 1e-12)
    axc.set_ylim(ylo - pad, yhi + pad)
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

    if args.vapor:
        # Under the zoom panel, positioned from its DRAWN box: the panel is
        # aspect-equal, so its gridspec slot is not where it actually lands.
        fig.canvas.draw()
        bb = axz.get_window_extent().transformed(fig.transFigure.inverted())
        cax = fig.add_axes([bb.x0, max(0.030, bb.y0 - 0.105),
                            bb.width, 0.022])
        cb = fig.colorbar(vZ, cax=cax, orientation="horizontal",
                          ticks=sigma_ticks(norm))
        cb.set_label(r"supersaturation  $\sigma = \rho_v/\rho_{vs}-1$   "
                     r"[$\times 10^{-4}$]", fontsize=9)
        cb.ax.xaxis.set_major_formatter(plt.FuncFormatter(
            lambda v, _p: f"{v:.3g}"))
        cb.ax.tick_params(labelsize=8)

    dx_lab = 0.025 * (XZ.max() - XZ.min())

    def draw(k):
        nonlocal contP, contZ
        fk = read_vts(files[k], want=WANT)[0]
        phi = fk["IcePhase"]
        pP_, _, _ = view(phi, Pi, Pj)
        pZ_, _, _ = view(phi, Zi, Zj)
        mP.set_array(pP_.ravel()); mZ.set_array(pZ_.ravel())
        if args.vapor:
            sk = sigma_field(fk)
            vP.set_array(view(sk, Pi, Pj)[0].ravel())
            vZ.set_array(view(sk, Zi, Zj)[0].ravel())
        contP.remove(); contZ.remove()
        # zorder must match the setup call: the ice layer is at 2, and a
        # contour left at its default would only land on top by draw order.
        contP = axp.contour(XP, YP, pP_, levels=[args.phi], colors=ICE_EDGE,
                            linewidths=0.8, zorder=4)
        contZ = axz.contour(XZ, YZ, pZ_, levels=[args.phi], colors=ICE_EDGE,
                            linewidths=1.1, zorder=4)

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
