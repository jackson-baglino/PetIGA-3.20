#!/usr/bin/env python3
"""poster_jpl_2026.py — TEST figures for the JPL poster day box (11.38 x 12.81 in).

    venv_enceladus/bin/python studies/keff_sintering/figures/poster_jpl_2026.py <campaign dir>
        --out <poster folder> [--seed 1702] [--phi 0.325]

Candidate figures for a space-science audience, not the manuscript: words on
the axes and colour bars, poster type sizes (sans), and two temperatures of
the SAME initial packing picked out in bold colours -- the frame round each
row of microstructures is the colour of that run's curve.

  option1_stacked_<Tw>_<Tc>   two rows x four instants, colour bars at the
                              right; below, k_eff/k_0 vs time (the curves fan
                              out) and vs SSA/SSA_0 (they fall on one curve)
  option1_compact_<Tw>_<Tc>   the same content in less than half the height:
                              two rows x three instants at the left, the two
                              plots stacked at the right
  option2_compact_<Tw>_<Tc>   as option1_compact, but the second plot is the
                              STATE LAW: every well-connected packing
                              (porosity <= 0.375, all seeds, all temperatures)
                              against SSA, with the fitted power law
  grid_only / curves_ssa / curves_law
                              the pieces on their own, for arranging in
                              PowerPoint as last year's poster did: the
                              2 x 4 grid with its colour bars, and the two
                              plots side by side (8.3 in, last year's slot)

DATA. One packing (--seed, --phi) run at five temperatures; the packings are
identical at t = 0, so k_eff,0 and SSA_0 are the same for every curve. Each
quantity is divided by its value at the opening sample (the first with
t >= 1 s). The numbered markers sit on the measured k_eff sample nearest each
snapshot's time. Option 2 follows the manuscript's state law (Figure 7a):
k_eff and SSA relative to their values at sintering age theta = t/tau_sub =
30, one exponent fitted through the origin in log-log per packing, averaged.

Colour bars: ice phase field (cmocean ice) and vapour supersaturation
(cmocean balance, pale = saturated, symmetric about zero, shared by every
panel shown).
"""
from __future__ import annotations

import argparse, re, sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import AsinhNorm
import cmocean

HERE = Path(__file__).resolve().parent
PROJ = HERE.parents[2]
sys.path.insert(0, str(PROJ / "postprocess"))
import pplib  # noqa: E402
from pplib import step_times  # noqa: E402
from plot_keff import load, read_tau_sub  # noqa: E402
from plot_keff_snapshots import make_reader, snap_step, _field, WANT, SIGMA_SCALE  # noqa: E402
from make_neck_movie import centered_cmap, ice_alpha_cmap, sigma_ticks  # noqa: E402
from compare_keff import CMAP  # noqa: E402

DAY = 86400.0
PAT = re.compile(r"phi([\d.]+)_Rave50um_LR40_seed(\d+)_L2mm_eps1000nm_perxy_T(-?\d+)__")
INK, MUTED, GREY = "#1a1a1a", "#5c5c5c", "#a3a3a3"
C_WARM, C_COLD = "#E8532B", "#1E88E5"        # bold: the warm and the cold run
FS, FS_S, FS_T = 15, 13, 12                  # labels | ticks, legends | small (poster captions are 16 pt Arial)
POSTER_RC = {
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold", "pdf.fonttype": 42, "svg.fonttype": "none",
    "axes.unicode_minus": True,
}


def tlabel(T):
    return f"{T} °C".replace("-", "−")


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------
def load_campaign(root: Path):
    """R[(phi, T)][seed] = run dict, opening sample first."""
    R = defaultdict(dict)
    for kf in sorted(root.glob("packing_*/k_eff.csv")):
        m = PAT.search(kf.parent.name)
        if not m or not 1601 <= int(m.group(2)) <= 2005:
            continue
        r = load(kf.parent)
        if r is None:
            continue
        i0 = int(np.argmax(r["t"] >= 1.0))
        for k in ("t", "kiso", "ssa", "step"):
            r[k] = r[k][i0:]
        r["tau"] = read_tau_sub(kf.parent)
        r["th"] = r["t"] / r["tau"]; r["td"] = r["t"] / DAY
        r["kn"] = r["kiso"] / r["kiso"][0]; r["sn"] = r["ssa"] / r["ssa"][0]
        r["dir"] = kf.parent
        R[(float(m.group(1)), int(m.group(3)))][int(m.group(2))] = r
    return R


def snapshots(run, days, stride):
    """[(fields, X, Y, t_s, k_row)] at the snapshots nearest `days`."""
    files, reader = make_reader(run["dir"], "sol")
    tm = step_times(str(run["dir"]))
    st = np.array([snap_step(f) for f in files])
    tt = np.array([tm.get(int(s), np.nan) for s in st])
    out = []
    for d in days:
        i = int(np.nanargmin(np.abs(tt - max(d * DAY, 1.0))))
        fl, X, Y = reader(files[i], want=WANT)
        sub = lambda a_: a_[::stride, ::stride]
        row = int(np.argmin(np.abs(run["t"] - tt[i])))
        out.append(({k: sub(v) for k, v in fl.items()}, sub(X), sub(Y), float(tt[i]), row))
    return out


def sigma_norm(all_snaps):
    pore = np.concatenate([
        (SIGMA_SCALE * pplib.supersaturation(fl["VaporDensity"], fl["Temperature"]))[fl["IcePhase"] < 0.5]
        for snaps in all_snaps for fl, *_ in snaps])
    smin, smax = float(pore.min()), float(pore.max())
    v = min(abs(smin), abs(smax)) or max(abs(smin), abs(smax))
    ext = {(True, True): "both", (True, False): "min", (False, True): "max",
           (False, False): "neither"}[(smin < -v, smax > v)]
    norm = AsinhNorm(linear_width=max(v / 300.0, 1e-12), vmin=-v, vmax=v)
    return norm, ext


def state_law(R, theta_ref=30.0, phi_max=0.375):
    """Per-run (x, y) = (SSA/SSA_r, k/k_r) after theta_ref and the mean exponent."""
    curves, per = [], defaultdict(list)
    for (p, T), d in sorted(R.items()):
        if p > phi_max:
            continue
        for seed, r in d.items():
            th = r["th"]
            if th[-1] < 1.5 * theta_ref:
                continue
            kr = np.interp(theta_ref, th, r["kiso"]); sr = np.interp(theta_ref, th, r["ssa"])
            w = th >= theta_ref
            x, y = r["ssa"][w] / sr, r["kiso"][w] / kr
            curves.append((p, T, seed, x, y))
            per[(p, seed)].append((np.log10(x), np.log10(y)))
    exps = []
    for L in per.values():
        X = np.concatenate([l[0] for l in L]); Y = np.concatenate([l[1] for l in L])
        exps.append(float(np.sum(X * Y) / np.sum(X * X)))
    return curves, float(np.mean(exps)), float(np.std(exps, ddof=1)), len(exps)


# ---------------------------------------------------------------------------
# drawing
# ---------------------------------------------------------------------------
def num_marker(ax, x, y, n, color, ms=15):
    ax.plot([x], [y], "o", ms=ms, mfc="white", mec=color, mew=2.0, zorder=6, clip_on=False)
    ax.text(x, y, str(n), ha="center", va="center", fontsize=FS_T - 1, fontweight="bold",
            color=INK, zorder=7, clip_on=False)


def clean(ax):
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_linewidth(0.9)
    ax.tick_params(labelsize=FS_S, width=0.9, length=4)


def snap_grid(fig, F, x0, y_top, s, gap, rows, norm, vapcm, icecm, with_titles=True):
    """rows = [(label, colour, snaps)]. Top-left corner (x0, y_top) in inches.
    Returns the block's bottom y."""
    n = len(rows[0][2])
    lab_w = 0.34
    for r_i, (label, col, snaps) in enumerate(rows):
        y = y_top - (r_i + 1) * s - r_i * gap
        fig.text(*F(x0 + 0.5 * lab_w - 0.06, y + 0.5 * s, 0, 0)[:2], label, rotation=90,
                 ha="center", va="center", fontsize=FS, fontweight="bold", color=col)
        for c_i, (fl, X, Y, t, _row) in enumerate(snaps):
            ax = fig.add_axes(F(x0 + lab_w + c_i * (s + gap), y, s, s))
            _field(ax, fl, X, Y, norm, vapcm, icecm)
            for sp in ax.spines.values():
                sp.set_linewidth(3.2); sp.set_color(col)
            if r_i == 0 and with_titles:
                d = t / DAY
                ax.set_title(f"{c_i + 1}   " + ("start" if d < 0.01 else f"{d:.0f} days"),
                             fontsize=FS, color=INK, pad=6)
            if r_i == len(rows) - 1 and c_i == 0:           # 500 um scale bar
                x_ = X[0] * 1e3
                L0, x_a = 0.5, x_.min() + 0.07 * (x_.max() - x_.min())
                y_a = (Y[:, 0] * 1e3).min() + 0.08 * (x_.max() - x_.min())
                for lw, c in ((6.0, "white"), (3.2, INK)):
                    ax.plot([x_a, x_a + L0], [y_a, y_a], color=c, lw=lw, solid_capstyle="butt", zorder=5)
                ax.text(x_a + L0 / 2, y_a + 0.05, "0.5 mm", ha="center", va="bottom",
                        fontsize=FS_T, color=INK, zorder=5,
                        bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
    return y_top - len(rows) * s - (len(rows) - 1) * gap, lab_w + n * s + (n - 1) * gap


def colorbars(fig, cax_ice, cax_sig, norm, vapcm, ext, orientation):
    cb = fig.colorbar(ScalarMappable(cmap=cmocean.cm.ice, norm=plt.Normalize(0, 1)), cax=cax_ice,
                      orientation=orientation, ticks=[0, 0.5, 1])
    cb.set_ticklabels(["air", "", "ice"]); cb.ax.tick_params(labelsize=FS_S, width=0.8, length=3)
    cb.set_label("phase", fontsize=FS_S); cb.outline.set_linewidth(0.8)
    cb = fig.colorbar(ScalarMappable(cmap=vapcm, norm=norm), cax=cax_sig, orientation=orientation,
                      extend=ext, extendfrac=0.05, ticks=sigma_ticks(norm, min_gap=0.22))
    fmt = plt.FuncFormatter(lambda v, _p: f"{v:.2g}".replace("-", "−"))
    (cb.ax.xaxis if orientation == "horizontal" else cb.ax.yaxis).set_major_formatter(fmt)
    cb.ax.tick_params(labelsize=FS_S, width=0.8, length=3)
    cb.set_label(r"vapor supersaturation  [$\times10^{-4}$]", fontsize=FS_S)
    cb.outline.set_linewidth(0.8)


K_LABEL = {True: r"$k_\mathrm{eff}\,/\,k_\mathrm{eff,0}$",
           False: r"conductivity  $k_\mathrm{eff}\,/\,k_\mathrm{eff,0}$"}


def panel_time(ax, runs, hot, cold, marks, short=False):
    """k_eff/k_0 against time: the temperatures fan out."""
    def end_label(r, T, col, bold, dx):
        ax.annotate(tlabel(T), (r["td"][-1], r["kn"][-1]), xytext=(dx, 0),
                    textcoords="offset points", color=col, va="center", ha="left",
                    fontsize=FS_S if bold else FS_T, fontweight="bold" if bold else "normal",
                    annotation_clip=False)
    for T, r in sorted(runs.items()):
        if T in (hot, cold):
            continue
        ax.plot(r["td"], r["kn"], color=GREY, lw=1.4, zorder=2)
        end_label(r, T, MUTED, False, 6)
    for T, col in ((hot, C_WARM), (cold, C_COLD)):
        r = runs[T]
        ax.plot(r["td"], r["kn"], color=col, lw=3.2, zorder=4)
        end_label(r, T, col, True, 14)                  # clear of the end marker
        for n, row in enumerate(marks[T], 1):
            num_marker(ax, r["td"][row], r["kn"][row], n, col)
    # the axis stops at the data; the labels sit in the margin to its right
    ax.set_xlim(0, 30.3); ax.set_xticks([0, 10, 20, 30])
    ax.set_xlabel("time  [days]", fontsize=FS)
    ax.set_ylabel(K_LABEL[short], fontsize=FS)
    clean(ax)


def panel_ssa(ax, runs, hot, cold, marks, short=False):
    """The same curves against SSA/SSA_0: one path, travelled at different speeds."""
    for T, r in sorted(runs.items()):
        if T not in (hot, cold):
            ax.plot(r["sn"], r["kn"], color=GREY, lw=1.4, zorder=2)
    ax.plot(runs[hot]["sn"], runs[hot]["kn"], color=C_WARM, lw=5.0, zorder=3, solid_capstyle="round")
    ax.plot(runs[cold]["sn"], runs[cold]["kn"], color=C_COLD, lw=2.2, zorder=4, solid_capstyle="round")
    for T, col, up in ((hot, C_WARM, False), (cold, C_COLD, True)):
        r = runs[T]
        for n, row in enumerate(marks[T], 1):
            if n == 1 and T == cold:
                continue                      # the shared start is marked once
            num_marker(ax, r["sn"][row], r["kn"][row], n, col)
    ax.invert_xaxis()
    ax.set_xlabel(r"specific surface area  SSA$\,/\,$SSA$_0$", fontsize=FS)
    ax.set_ylabel(K_LABEL[short], fontsize=FS)
    ax.annotate("", xy=(0.62, 0.10), xytext=(0.30, 0.10), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.4))
    ax.text(0.46, 0.13, "sintering", transform=ax.transAxes, ha="center", va="bottom",
            fontsize=FS_T, color=MUTED)
    clean(ax)


def panel_law(ax, R, runs_hl, hot, cold, phi_sel, short=False):
    """The state law across packings, with the fitted power law."""
    curves, p0, sd, n = state_law(R)
    phis = sorted({c[0] for c in curves})
    # Greys, light -> dark with porosity: colour is kept for the two runs
    # whose microstructures are shown.
    greys = ["#c4c4c4", "#909090", "#5a5a5a"]
    pc = {p: greys[min(i, 2)] for i, p in enumerate(phis)}
    for p, T, seed, x, y in curves:
        ax.plot(x, y, color=pc[p], lw=0.9, alpha=0.8, zorder=2)
    xx = np.linspace(min(c[3].min() for c in curves), 1.0, 50)
    ax.plot(xx, xx ** p0, color=INK, lw=2.4, ls=(0, (5, 2.5)), zorder=5)
    for T, col in ((hot, C_WARM), (cold, C_COLD)):
        for p, T_, seed, x, y in curves:
            if p == phi_sel and T_ == T and seed == runs_hl:
                ax.plot(x, y, color=col, lw=3.0, zorder=6)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.invert_xaxis()
    ax.set_xticks([1.0, 0.9, 0.8, 0.7]); ax.set_xticklabels(["1.0", "0.9", "0.8", "0.7"])
    yt = [1.0, 1.1, 1.2, 1.3]
    ax.set_yticks(yt); ax.set_yticklabels([f"{v:.1f}" for v in yt])
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel(r"specific surface area  SSA$\,/\,$SSA$_\mathrm{ref}$", fontsize=FS)
    ax.set_ylabel(r"$k_\mathrm{eff}\,/\,k_\mathrm{eff,ref}$" if short else
                  r"conductivity  $k_\mathrm{eff}\,/\,k_\mathrm{eff,ref}$", fontsize=FS)
    h = [plt.Line2D([], [], color=pc[p], lw=2.0, label=f"porosity {p:g}") for p in phis]
    h.append(plt.Line2D([], [], color=INK, lw=2.4, ls=(0, (5, 2.5)),
                        label=rf"fit:  $k_\mathrm{{eff}}\propto\mathrm{{SSA}}^{{{p0:.2f}}}$".replace("-", "−")))
    ax.legend(handles=h, fontsize=FS_T, frameon=False, loc="upper left", handlelength=2.0,
              labelspacing=0.3)
    clean(ax)
    print(f"  state law: exponent {p0:.3f} +- {sd:.3f} over {n} packings, {len(curves)} runs")
    return p0, sd


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------
def build(kind, R, a, hot, cold):
    runs = {T: R[(a.phi, T)][a.seed] for (p, T) in R if p == a.phi and a.seed in R[(p, T)]}
    compact = "compact" in kind
    days = (0.0, 10.0, 30.0) if compact else (0.0, 10.0, 19.4, 30.0)
    S = {T: snapshots(runs[T], days, a.stride) for T in (hot, cold)}
    marks = {T: [s_[4] for s_ in S[T]] for T in S}
    W = a.width_in
    note = (f"identical packings at the start  ·  2 mm × 2 mm periodic cell  ·  "
            f"porosity {a.phi:g}  ·  mean grain radius 50 µm")
    second = panel_law if "option2" in kind or "law" in kind else panel_ssa

    def grid(fig, F, x0, y_top, s, gap):
        norm, ext = sigma_norm(S.values())
        vapcm, icecm = centered_cmap(cmocean.cm.balance, norm), ice_alpha_cmap()
        rows = [(tlabel(hot), C_WARM, S[hot]), (tlabel(cold), C_COLD, S[cold])]
        y_bot, bw = snap_grid(fig, F, x0, y_top, s, gap, rows, norm, vapcm, icecm)
        return y_bot, bw, norm, vapcm, ext

    def second_panel(ax, short):
        if second is panel_law:
            panel_law(ax, R, a.seed, hot, cold, a.phi, short)
        else:
            panel_ssa(ax, runs, hot, cold, marks, short)

    if kind.startswith("curves"):                   # the two plots alone, side by side
        Wc, ph, bot, top = a.curves_width_in, 2.55, 0.68, 0.22
        H = top + ph + bot
        fig = plt.figure(figsize=(Wc, H)); F = lambda x, y, w, h: (x / Wc, y / H, w / Wc, h / H)
        l, mid, r_ = 0.95, 1.75, 0.42              # mid: end labels + next y label
        pw = (Wc - l - mid - r_) / 2
        panel_time(fig.add_axes(F(l, bot, pw, ph)), runs, hot, cold, marks)
        second_panel(fig.add_axes(F(l + pw + mid, bot, pw, ph)), False)
        return fig, (Wc, H)

    if kind == "grid_only" or kind == "option1_stacked":
        gap, cbw = 0.09, 1.55                       # colour bars in a column at the right
        s = (W - 0.10 - 0.34 - cbw - 3 * gap) / 4
        top = 0.42
        if kind == "grid_only":
            lower = 0.42
        else:
            l, mid, r_, ph, bot = 0.95, 1.75, 0.20, 2.75, 0.68
            lower = 0.42 + 0.45 + ph + bot
        H = top + 2 * s + gap + lower
        fig = plt.figure(figsize=(W, H)); F = lambda x, y, w, h: (x / W, y / H, w / W, h / H)
        y_bot, bw, norm, vapcm, ext = grid(fig, F, 0.10, H - top, s, gap)
        xc = 0.10 + bw + 0.42
        colorbars(fig, fig.add_axes(F(xc, y_bot + s + gap + 0.25 * s, 0.16, 0.6 * s)),
                  fig.add_axes(F(xc, y_bot + 0.02 * s, 0.16, 0.95 * s)), norm, vapcm, ext, "vertical")
        fig.text(*F(0.10 + 0.34, y_bot - 0.24, 0, 0)[:2], note, fontsize=FS_T, color=MUTED,
                 ha="left", va="center")
        if kind == "option1_stacked":
            pw = (W - l - mid - r_) / 2
            panel_time(fig.add_axes(F(l, bot, pw, ph)), runs, hot, cold, marks)
            second_panel(fig.add_axes(F(l + pw + mid, bot, pw, ph)), False)
        return fig, (W, H)

    # compact: snapshots at the left, the two plots stacked at the right
    gap, s = 0.08, a.compact_snap_in
    top, cb_strip = 0.40, 1.42                      # bars, their labels, the note
    H = top + 2 * s + gap + cb_strip
    fig = plt.figure(figsize=(W, H)); F = lambda x, y, w, h: (x / W, y / H, w / W, h / H)
    y_bot, bw, norm, vapcm, ext = grid(fig, F, 0.08, H - top, s, gap)
    x_s = 0.08 + 0.34
    colorbars(fig, fig.add_axes(F(x_s + 0.05, y_bot - 0.30, 0.9, 0.13)),
              fig.add_axes(F(x_s + 1.45, y_bot - 0.30, bw - 0.34 - 1.55, 0.13)),
              norm, vapcm, ext, "horizontal")
    fig.text(*F(x_s, 0.22, 0, 0)[:2], note.replace("  ·  porosity", "\nporosity"), fontsize=FS_T - 1,
             color=MUTED, ha="left", va="center", linespacing=1.3)
    xl = 0.08 + bw + 0.92
    pw = W - xl - 0.78                              # 0.78: the curve-end labels
    b0, g_ = 0.66, 0.78
    ph = (H - b0 - g_ - 0.16) / 2
    panel_time(fig.add_axes(F(xl, b0 + ph + g_, pw, ph)), runs, hot, cold, marks, short=True)
    second_panel(fig.add_axes(F(xl, b0, pw, ph)), True)
    return fig, (W, H)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--phi", type=float, default=0.325)
    ap.add_argument("--seed", type=int, default=1702)
    ap.add_argument("--width-in", type=float, default=11.0,
                    help="figure width [in]; the poster box is 11.38 in wide")
    ap.add_argument("--compact-snap-in", type=float, default=1.72)
    ap.add_argument("--stride", type=int, default=2, help="field down-sampling")
    ap.add_argument("--pairs", nargs="+", default=["-5,-20", "-5,-40"],
                    help="warm,cold temperature pairs to build")
    ap.add_argument("--curves-width-in", type=float, default=8.3,
                    help="width of the curves-only pieces [in] (last year's plot slot)")
    ap.add_argument("--kinds", nargs="+",
                    default=["option1_stacked", "option1_compact", "option2_compact",
                             "grid_only", "curves_ssa", "curves_law"])
    ap.add_argument("--dpi", type=int, default=300)
    a = ap.parse_args()
    plt.rcParams.update(POSTER_RC)
    R = load_campaign(a.root)
    a.out.mkdir(parents=True, exist_ok=True)
    for pair in a.pairs:
        hot, cold = (int(v) for v in pair.split(","))
        for kind in a.kinds:
            fig, (Wf, Hf) = build(kind, R, a, hot, cold)
            name = f"{kind}_T{hot}_T{cold}"
            for e in ("png", "pdf"):
                fig.savefig(a.out / f"{name}.{e}", dpi=a.dpi)
            plt.close(fig)
            print(f"  wrote {name}.png/.pdf   {Wf:.2f} x {Hf:.2f} in")


if __name__ == "__main__":
    main()
