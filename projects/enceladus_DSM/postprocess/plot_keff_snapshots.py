#!/usr/bin/env python3
"""plot_keff_snapshots.py — k_eff curve with microstructure snapshots under it.

    python3 plot_keff_snapshots.py --dir <run> [--times 1 10 20 30] [--steps ...]
        [--n-snapshots 4] [--width 7.2] [--absolute] [--iso-only] [--leaders auto|on|off]
        [--title-time STR] [--title-ssa STR] [--no-title]
        [--xlabel-time STR] [--xlabel-ssa STR] [--ylabel STR]
        [--source vts|sol] [--save-dir DIR] [--formats pdf png]

Two page-width manuscript figures, each a row of 3-4 microstructure panels
above one k_eff curve. Each panel is lettered, and the same letter marks its
instant on the curve, joined to it by a leader line:

    keff_time_snapshots.{pdf,png}   k_eff / k_b  vs  time [d]
    keff_ssa_snapshots.{pdf,png}    k_eff / k_b  vs  SSA / SSA_b

The point is to show both things at once: k_eff rises as the packing
sinters, and the microstructure at a few instants shows what that rise
looks like.

PANELS. They follow make_packing_movie.py, and they import its helpers so
the figure and the movie cannot drift apart: supersaturation
sigma = rho_v/rho_vs(T) - 1 on cmocean `balance`, re-centred so the pale
middle is sigma = 0, on an asinh scale; ice painted on top with cmocean
`ice`, transparent below phi = 0.5. The sigma range is taken over the pore
space of the SHOWN snapshots only, and shared by every panel, so the panels
are directly comparable.

CURVE. Same data, normalization and baseline as plot_keff.py, which it
imports: k and SSA are divided by their value at the first sample at
t >= --baseline-days (default 1 d), and the samples before it are the IC
relaxing, drawn grey. --absolute plots k_eff in W m^-1 K^-1 on the time
figure instead; the SSA figure is always normalized on both axes.

SNAPSHOTS. --times (days) picks the snapshot nearest each time; --steps
names them exactly. By default --n-snapshots are spaced evenly in time from
the baseline to the end of the run. The marker sits on the k_eff sample of
the snapshot's own step.

OUTPUT. <dir>/plots/keff/snapshots/ (or --save-dir). PDF is the manuscript
file -- the curve is vector, the fields are embedded rasters -- and PNG is a
preview. Exits 0 without writing anything when the run has no k_eff CSV.
"""
from __future__ import annotations

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
from matplotlib.colors import AsinhNorm, ListedColormap
from matplotlib.cm import ScalarMappable
from matplotlib.patches import ConnectionPatch
import cmocean

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import pplib                                                # noqa: E402
from pplib import read_vts, step_times                      # noqa: E402
from plot_keff import load, DAY, C_XX, C_YY, C_ISO, C_RELAX  # noqa: E402
from make_neck_movie import (SIGMA_SCALE, ice_alpha_cmap,   # noqa: E402
                             centered_cmap, sigma_ticks)

WANT = ("IcePhase", "VaporDensity", "Temperature")
SOL_DOF = {"IcePhase": 0, "Temperature": 1, "VaporDensity": 2}
LETTERS = "abcdefgh"
INK, MUTED = "#1a1a1a", "#5c5c5c"
FS, FS_SMALL = 9, 8                 # page-width figure: 8-9 pt at print size

DEFAULTS = {
    "title_time": "Effective thermal conductivity during dry-snow metamorphism",
    "title_ssa": "Effective thermal conductivity vs specific surface area",
    "xlabel_time": "Time [d]",
    "xlabel_ssa": r"Normalized specific surface area  SSA$\,/\,$SSA$_b$",
    "ylabel_norm": r"Normalized conductivity  $k_\mathrm{eff}\,/\,k_{\mathrm{eff},b}$",
    "ylabel_abs": r"$k_\mathrm{eff}$  [W m$^{-1}$ K$^{-1}$]",
}


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
def snap_step(fn) -> int:
    return int(re.search(r"sol[V]?_(\d+)\.(?:vts|dat)$", str(fn)).group(1))


def make_reader(run: Path, source: str):
    """(files, reader) for the snapshots; see make_packing_movie --source."""
    if source == "sol":
        from igakit.io import PetIGA
        nrb = PetIGA().read(str(run / "igasol.dat"))
        P = nrb.points
        X, Y = P[..., 0].T, P[..., 1].T

        def read(fn, want=WANT):
            sol = PetIGA().read_vec(str(fn), nrb)
            return {k: sol[..., SOL_DOF[k]].T for k in want}, X, Y
        files = glob.glob(str(run / "sol_*.dat"))
    else:
        read = read_vts
        files = glob.glob(str(run / "vtkOut" / "solV_*.vts"))
    return sorted(files, key=snap_step), read


def pick_snapshots(files, tmap, d, times_d, steps, n, t_base):
    """Indices into `files` for the panels, in time order, duplicates removed."""
    fsteps = np.array([snap_step(f) for f in files])
    ftimes = np.array([tmap.get(int(s), np.nan) for s in fsteps])
    if steps:
        want = [int(np.argmin(np.abs(fsteps - s))) for s in steps]
    else:
        if not times_d:
            times_d = np.linspace(t_base / DAY, d["t"][-1] / DAY, n)
        ok = np.isfinite(ftimes)
        cand = np.flatnonzero(ok)
        want = [int(cand[np.argmin(np.abs(ftimes[ok] / DAY - td))]) for td in times_d]
    return sorted(set(want), key=lambda i: fsteps[i]), fsteps, ftimes


def sample_index(d, step, t):
    """Row of the k_eff series for this snapshot: same step, else nearest time."""
    hit = np.flatnonzero(d["step"] == step)
    return int(hit[0]) if hit.size else int(np.argmin(np.abs(d["t"] - t)))


# ---------------------------------------------------------------------------
# Drawing
# ---------------------------------------------------------------------------
def _field(ax, fl, X, Y, norm, vapcm, icecm):
    sig = SIGMA_SCALE * pplib.supersaturation(fl["VaporDensity"], fl["Temperature"])
    XX, YY = X * 1e3, Y * 1e3
    h = 0.5 * (XX[0, 1] - XX[0, 0])
    kw = dict(origin="lower", extent=(XX.min() - h, XX.max() + h,
                                      YY.min() - h, YY.max() + h),
              interpolation="antialiased", interpolation_stage="rgba")
    ax.imshow(sig, cmap=vapcm, norm=norm, **kw)
    ax.imshow(fl["IcePhase"], cmap=icecm, vmin=0.0, vmax=1.0, **kw)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_linewidth(0.6); sp.set_color(MUTED)
    return XX, YY


def _scalebar(ax, XX, YY):
    """A bar ~1/4 of the domain, rounded to 1/2/5 x 10^k mm, bottom left."""
    Lx = XX.max() - XX.min()
    raw = Lx / 4
    e = 10 ** np.floor(np.log10(raw))
    L = max(m * e for m in (1, 2, 5) if m * e <= raw)
    x0 = XX.min() + 0.06 * Lx
    y0 = YY.min() + 0.07 * (YY.max() - YY.min())
    for lw, c in ((3.4, "white"), (1.8, INK)):
        ax.plot([x0, x0 + L], [y0, y0], color=c, lw=lw, solid_capstyle="butt",
                zorder=5)
    txt = f"{L * 1e3:g} µm" if L < 1 else f"{L:g} mm"
    t = ax.text(x0 + L / 2, y0 + 0.025 * Lx, txt, ha="center", va="bottom",
                fontsize=FS_SMALL - 1, color=INK, zorder=5)
    t.set_bbox(dict(facecolor="white", alpha=0.75, lw=0, pad=0.8))


def _curve(ax, x, ys, relax, keys, legend_loc):
    live = ~relax
    lw = {"kxx": 1.0, "kyy": 1.0, "kiso": 1.8}
    col = {"kxx": C_XX, "kyy": C_YY, "kiso": C_ISO}
    lab = {"kxx": r"$k_{xx}$", "kyy": r"$k_{yy}$", "kiso": r"$k_\mathrm{iso}$"}
    for key in keys:
        y = ys[key]
        if relax.any():
            # Grey segment runs through the first live sample so the two join.
            grey = relax.copy()
            if live.any():
                grey[int(np.argmax(live))] = True
            ax.plot(x[grey], y[grey], "-", color=C_RELAX, lw=lw[key], zorder=1,
                    label="IC relaxation (first day)" if key == "kiso" else None)
        ax.plot(x[live], y[live], "-", color=col[key], lw=lw[key], zorder=2,
                label=lab[key])
    # The three curves sit within a few percent of each other, so direct labels
    # at their ends collide; a legend is the readable choice here.
    h, l = ax.get_legend_handles_labels()
    order = sorted(range(len(l)), key=lambda i: "relaxation" in l[i])
    ax.legend([h[i] for i in order], [l[i] for i in order], fontsize=FS_SMALL,
              frameon=False, loc=legend_loc, handlelength=1.6)
    ax.grid(True, alpha=0.25, lw=0.5)
    ax.tick_params(labelsize=FS_SMALL, width=0.6, length=3)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_linewidth(0.6)


def _colorbars(fig, cax_ice, cax_sig, norm, vapcm):
    """Two slim bars, each titled with its symbol so no rotated label is needed."""
    ice_lut = ListedColormap(cmocean.cm.ice(np.linspace(0.5, 1.0, 256)))
    cb = fig.colorbar(ScalarMappable(cmap=ice_lut, norm=plt.Normalize(0.5, 1.0)),
                      cax=cax_ice, ticks=[0.5, 0.75, 1.0], format="%g")
    cax_ice.set_title(r"$\phi_i$", fontsize=FS, pad=4)
    cb.ax.tick_params(labelsize=FS_SMALL - 1, width=0.5, length=2)
    cb.outline.set_linewidth(0.5)

    cb = fig.colorbar(ScalarMappable(cmap=vapcm, norm=norm), cax=cax_sig,
                      ticks=sigma_ticks(norm, min_gap=0.09), extend="both",
                      extendfrac=0.04)
    # Left-aligned: centred, it would run back over the ice bar's ticks.
    cax_sig.set_title(r"$\sigma$ [$\times10^{-4}$]", fontsize=FS_SMALL, pad=4,
                      loc="left")
    cb.ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _p: f"{v:.2g}"))
    cb.ax.tick_params(labelsize=FS_SMALL - 1, width=0.5, length=2)
    cb.outline.set_linewidth(0.5)


def build(kind, snaps, d, x, ys, relax, keys, norm, vapcm, icecm, a):
    """One figure. `snaps` is [(fields, X, Y, t, row)] in time order.

    Laid out in INCHES, not with a gridspec: the panel row, the colour-bar
    strip and the curve all have to line up to the same edges, and the height
    has to follow from the square panels, which a gridspec does not give.
    """
    n = len(snaps)
    W = a.width
    ml, mr = 0.62, 1.00            # y label | colour-bar strip
    gap = 0.08                     # between panels
    pw = (W - ml - mr - (n - 1) * gap) / n
    top = 0.34 if not a.no_title else 0.06
    t_band = 0.24                  # panel titles
    g_mid = 0.42                   # panels -> curve, room for the leaders
    ch = max(1.9, 1.15 * pw)       # curve height
    bot = 0.46                     # x label
    H = top + t_band + pw + g_mid + ch + bot
    fig = plt.figure(figsize=(W, H))
    F = lambda x0, y0, w, h: (x0 / W, y0 / H, w / W, h / H)

    y_pan = bot + ch + g_mid
    axs = [fig.add_axes(F(ml + i * (pw + gap), y_pan, pw, pw)) for i in range(n)]
    xc = ml + n * pw + (n - 1) * gap + 0.14
    cax_ice = fig.add_axes(F(xc, y_pan, 0.08, pw))
    cax_sig = fig.add_axes(F(xc + 0.08 + 0.36, y_pan + 0.05, 0.08, pw - 0.10))
    axc = fig.add_axes(F(ml, bot, n * pw + (n - 1) * gap, ch))
    axc.patch.set_alpha(0.0)       # so the leaders show through to the markers

    _curve(axc, x, ys, relax, keys,
           legend_loc="lower right" if kind == "time" else "lower left")
    if kind == "time" and relax.any():
        axc.axvspan(0, a.baseline_days, color="#f2f2f2", zorder=0, lw=0)
    if not a.absolute or kind == "ssa":
        axc.axhline(1.0, color="#999999", lw=0.6, ls=":", zorder=0)
    if kind == "ssa":
        # SSA falls as the packing sinters, so time runs right to left.
        axc.axvline(1.0, color="#999999", lw=0.6, ls=":", zorder=0)
        axc.annotate("", xy=(0.42, 0.93), xytext=(0.58, 0.93),
                     xycoords="axes fraction",
                     arrowprops=dict(arrowstyle="->", color=MUTED, lw=0.9))
        axc.text(0.59, 0.93, "time", transform=axc.transAxes, va="center",
                 fontsize=FS_SMALL, color=MUTED)

    xlabel = (a.xlabel_time or DEFAULTS["xlabel_time"]) if kind == "time" \
        else (a.xlabel_ssa or DEFAULTS["xlabel_ssa"])
    ylabel = a.ylabel or (DEFAULTS["ylabel_abs"] if (a.absolute and kind == "time")
                          else DEFAULTS["ylabel_norm"])
    axc.set_xlabel(xlabel, fontsize=FS)
    axc.set_ylabel(ylabel, fontsize=FS)
    ylo, yhi = axc.get_ylim()
    axc.set_ylim(ylo, yhi + 0.06 * (yhi - ylo))     # headroom for the letters

    # Leaders only where they cannot cross: against time the markers run left
    # to right like the panels; against SSA they run the other way.
    leaders = a.leaders == "on" or (a.leaders == "auto" and kind == "time")
    ymark = ys["kiso"]
    for i, (fl, X, Y, t, row) in enumerate(snaps):
        ax = axs[i]
        XX, YY = _field(ax, fl, X, Y, norm, vapcm, icecm)
        if i == 0:
            _scalebar(ax, XX, YY)
        L = LETTERS[i]
        ax.set_title(f"({L})  t = {t / DAY:.1f} d", fontsize=FS, color=INK, pad=3)
        px, py = x[row], ymark[row]
        axc.plot([px], [py], "o", ms=6.5, mfc="white", mec=INK, mew=1.1, zorder=6)
        axc.annotate(L, (px, py), xytext=(0, 6), textcoords="offset points",
                     ha="center", va="bottom", fontsize=FS_SMALL, color=INK,
                     fontweight="bold", zorder=7)
        if leaders:
            # Ends just above the letter, so it points at the marker without
            # striking through it.
            con = ConnectionPatch(xyA=(0.5, 0.0), coordsA=ax.transAxes,
                                  xyB=(px, py), coordsB=axc.transData,
                                  color="#a0a0a0", lw=0.6, ls=(0, (3, 2)),
                                  zorder=3, shrinkA=2, shrinkB=18)
            fig.add_artist(con)

    _colorbars(fig, cax_ice, cax_sig, norm, vapcm)

    if not a.no_title:
        title = (a.title_time or DEFAULTS["title_time"]) if kind == "time" \
            else (a.title_ssa or DEFAULTS["title_ssa"])
        fig.suptitle(title, x=(ml + 0.5 * (n * pw + (n - 1) * gap)) / W,
                     y=1.0 - 0.06 / H, ha="center", va="top",
                     fontsize=FS + 1.5, color=INK)
    return fig


# ---------------------------------------------------------------------------
def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dir", default=".", help="run directory (default: cwd)")
    p.add_argument("--save-dir", default=None,
                   help="output directory (default: <dir>/plots/keff/snapshots)")
    p.add_argument("--times", type=float, nargs="+", default=None,
                   help="snapshot times [d]; the nearest snapshot to each is used")
    p.add_argument("--steps", type=int, nargs="+", default=None,
                   help="snapshot steps, overriding --times")
    p.add_argument("--n-snapshots", type=int, default=4,
                   help="default count, evenly spaced baseline..end (default 4)")
    p.add_argument("--baseline-days", type=float, default=1.0,
                   help="normalization baseline, as in plot_keff.py (default 1)")
    p.add_argument("--source", choices=("vts", "sol"), default="vts",
                   help="vts: vtkOut/solV_*.vts (plenty at page width); "
                        "sol: full-resolution sol_*.dat via igakit")
    p.add_argument("--width", type=float, default=7.2,
                   help="figure width [in]; 7.2 is a full two-column page")
    p.add_argument("--absolute", action="store_true",
                   help="time figure: k_eff in W/m/K instead of k/k_b")
    p.add_argument("--iso-only", action="store_true",
                   help="plot k_iso only, not k_xx and k_yy")
    p.add_argument("--leaders", choices=("auto", "on", "off"), default="auto",
                   help="dashed lines from each panel to its marker. auto: on "
                        "the time figure only -- against SSA time runs right "
                        "to left and the lines would cross")
    p.add_argument("--sat-clip", type=float, default=0.1,
                   help="percentile clipped off each end of sigma (default 0.1)")
    p.add_argument("--title-time", default=None)
    p.add_argument("--title-ssa", default=None)
    p.add_argument("--no-title", action="store_true")
    p.add_argument("--xlabel-time", default=None)
    p.add_argument("--xlabel-ssa", default=None)
    p.add_argument("--ylabel", default=None, help="y label for both figures")
    p.add_argument("--formats", nargs="+", default=["pdf", "png"])
    p.add_argument("--dpi", type=int, default=400, help="PNG and raster dpi")
    a = p.parse_args(argv)

    run = Path(a.dir).resolve()
    d = load(run)
    if d is None:
        print(f"  no k_eff CSV (or no SSA_evo.dat / -Lx -Ly) in {run}; nothing to plot")
        return 0
    kf = np.atleast_1d(np.genfromtxt(run / d["csv"], delimiter=",", names=True))
    d["step"] = kf["step"].astype(int)

    files, reader = make_reader(run, a.source)
    if not files:
        print(f"  no {'sol_*.dat' if a.source == 'sol' else 'vtkOut/solV_*.vts'} "
              f"in {run}; nothing to plot")
        return 0
    tmap = step_times(str(run))
    for s, t in zip(d["step"], d["t"]):          # k_eff CSV covers every step
        tmap.setdefault(int(s), float(t))

    relax = d["t"] < a.baseline_days * DAY
    ib = int(np.argmax(~relax)) if not relax.all() else len(relax) - 1
    idx, fsteps, ftimes = pick_snapshots(files, tmap, d, a.times, a.steps,
                                         a.n_snapshots, d["t"][ib])
    if len(idx) > len(LETTERS):
        print(f"  at most {len(LETTERS)} snapshots", file=sys.stderr)
        return 1

    snaps, pore = [], []
    for i in idx:
        fl, X, Y = reader(files[i], want=WANT)
        t = ftimes[i] if np.isfinite(ftimes[i]) else tmap.get(int(fsteps[i]), 0.0)
        snaps.append([fl, X, Y, t, sample_index(d, int(fsteps[i]), t)])
        s = SIGMA_SCALE * pplib.supersaturation(fl["VaporDensity"], fl["Temperature"])
        pore.append(s[fl["IcePhase"] < 0.5])
        print(f"  snapshot step {fsteps[i]}: t = {t / DAY:.2f} d")
    pore = np.concatenate(pore)
    lo, hi = np.percentile(pore, [a.sat_clip, 100.0 - a.sat_clip])
    norm = AsinhNorm(linear_width=max(max(abs(lo), abs(hi)) / 300.0, 1e-12),
                     vmin=lo, vmax=hi)
    vapcm = centered_cmap(cmocean.cm.balance, norm)
    icecm = ice_alpha_cmap()

    keys = ("kiso",) if a.iso_only else ("kxx", "kyy", "kiso")
    kn = {k: d[k] / d[k][ib] for k in ("kxx", "kyy", "kiso")}
    out = Path(a.save_dir) if a.save_dir else run / "plots" / "keff" / "snapshots"
    os.makedirs(out, exist_ok=True)

    plt.rcParams.update({"font.family": "sans-serif", "mathtext.fontset": "dejavusans",
                         "pdf.fonttype": 42, "svg.fonttype": "none"})
    figs = {
        "keff_time_snapshots": build(
            "time", snaps, d, d["t"] / DAY, d if a.absolute else kn, relax, keys,
            norm, vapcm, icecm, a),
        "keff_ssa_snapshots": build(
            "ssa", snaps, d, d["ssa"] / d["ssa"][ib], kn, relax, keys,
            norm, vapcm, icecm, a),
    }
    for stem, fig in figs.items():
        for fmt in a.formats:
            path = out / f"{stem}.{fmt}"
            fig.savefig(path, dpi=a.dpi, bbox_inches="tight", pad_inches=0.03)
            print(f"  wrote {path}")
        plt.close(fig)
    print(f"  sigma x{SIGMA_SCALE:g} range {lo:+.3g} .. {hi:+.3g} over the shown pores")
    return 0


if __name__ == "__main__":
    sys.exit(main())
