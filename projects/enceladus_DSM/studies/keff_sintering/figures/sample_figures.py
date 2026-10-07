#!/usr/bin/env python3
"""Sample manuscript figures for the k_eff results section (pitches, 2026-10-06).

    venv_enceladus/bin/python studies/keff_sintering/figures/sample_figures.py <campaign dir> [--out <dir>]

Four candidate figures, one per claim, built in the manuscript style
(pplib.MANUSCRIPT_RC; 170 mm wide; no titles; (a)(b)(c) panel labels; symbol
labels; >= 8 pt; transparent background):

  figA_collapse      temperature only sets the clock:
                     (a) k/k_0 vs time [d], five temperatures  -> fans out
                     (b) the same vs theta = t/tau_sub          -> one curve
                     (c) k_iso vs SSA, five temperatures        -> one path
  Figure7_state_law       the state law: (a) k_eff vs SSA, each relative to its value
                     at theta = 30, well-connected packings, with the power
                     law; (b) the level vs porosity; (c) the exponent vs porosity
  Figure7_state_law_alt     the same with curves: (a) k_eff(t) by porosity, the gallery's
                     instants marked; (b)-(d) as figB
  Figure8_timescales    the clock extrapolated: time to reach a sintering age as a
                     function of temperature and grain radius (vapour route),
                     with Enceladus conditions and Choukroun's 180 K point
  Figure5_gallery       what happens to the aggregates: one packing of each
                     porosity at four instants, with qualitative colour bars

Seeds are PAIRED across temperature (the packings common to every temperature).
Writes PNG + PDF to <campaign>/compare/figure_samples/ (or --out).
"""
from __future__ import annotations

import argparse, csv, re, sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, AsinhNorm
import matplotlib.patheffects as pe
import cmocean

HERE = Path(__file__).resolve().parent
PROJ = HERE.parents[2]
sys.path.insert(0, str(PROJ / "postprocess"))
import pplib  # noqa: E402
from plot_keff import load, read_tau_sub  # noqa: E402
from compare_keff import CMAP  # noqa: E402

MM = 1 / 25.4
W = 170 * MM
DAY, YEAR = 86400.0, 365.25 * 86400.0
K_ICE = 2.29
PAT = re.compile(r"phi([\d.]+)_Rave50um_LR40_seed(\d+)_L2mm_eps1000nm_perxy_T(-?\d+)__")
INK, MUTED = "#1a1a1a", "#6b6b6b"
FS, FS_S = 9, 8


def panel(ax, s, dx=-0.16, dy=1.02):
    ax.text(dx, dy, f"({s})", transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")


def clean(ax):
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(labelsize=FS_S, width=0.6, length=3)
    ax.xaxis.label.set_size(FS); ax.yaxis.label.set_size(FS)


def save(fig, out, name):
    for ext in ("png", "pdf"):
        fig.savefig(out / f"{name}.{ext}", dpi=400, transparent=True)
    plt.close(fig)
    print(f"  wrote {out / name}.png/.pdf")


def colours(vals, which):
    cm, (a, b) = CMAP[which]
    return {v: cm(a + (b - a) * i / max(1, len(vals) - 1)) for i, v in enumerate(vals)}


def mean_curve(g, xkey, ykey, n=300, logx=True):
    """Seed mean of y on a common x grid, over the span every seed covers."""
    lo = max(r[xkey][r["i0"]] if not logx else max(r[xkey][r["i0"]], 1e-9) for r in g)
    hi = min(r[xkey][-1] for r in g)
    xs = np.geomspace(max(lo, 1e-6), hi, n) if logx else np.linspace(lo, hi, n)
    Y = np.array([np.interp(xs, r[xkey][r["i0"]:], r[ykey][r["i0"]:]) for r in g])
    return xs, Y.mean(0), Y.min(0), Y.max(0)


def fig_collapse(R, out, phi=0.325):
    Ts = sorted({T for (p, T) in R if p == phi})
    common = set.intersection(*(set(R[(phi, T)]) for T in Ts))
    col = colours(Ts, "T")
    fig, ax = plt.subplots(1, 3, figsize=(W, 62 * MM))
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.20, top=0.90, wspace=0.42)
    for T in Ts:
        g = [R[(phi, T)][s] for s in sorted(common)]
        for r in g:
            r["kn"] = r["kiso"] / r["kiso"][r["i0"]]
        x, y, lo, hi = mean_curve(g, "td", "kn", logx=False)
        ax[0].plot(x, y, color=col[T], lw=1.6, label=rf"${T}\,^\circ$C".replace("-", "-"))
        x, y, lo, hi = mean_curve(g, "th", "kn")
        ax[1].plot(x, y, color=col[T], lw=1.6)
        # k vs SSA: mean k and mean SSA at common theta
        xs, ks, _, _ = mean_curve(g, "th", "kiso")
        _, ss, _, _ = mean_curve(g, "th", "ssa")
        ax[2].plot(ss / 1e3, ks, color=col[T], lw=1.6)
    ax[0].set(xlabel="$t$ [d]", ylabel=r"$k_\mathrm{eff}/k_{\mathrm{eff},0}$")
    ax[1].set(xlabel=r"$\theta=t/\tau_\mathrm{sub}$", ylabel=r"$k_\mathrm{eff}/k_{\mathrm{eff},0}$", xscale="log")
    ax[1].set_xlim(1, None)
    ax[2].set(xlabel=r"SSA [mm$^{-1}$]", ylabel=r"$k_\mathrm{eff}$ [W m$^{-1}$ K$^{-1}$]")
    ax[2].invert_xaxis()
    ax[0].legend(frameon=False, fontsize=FS_S, loc="lower right", handlelength=1.4, labelspacing=0.25)
    for a_, s in zip(ax, "abc"):
        clean(a_); panel(a_, s, dx=-0.27)
    save(fig, out, "figA_collapse")


def fig_closure(R, out, theta_ref=30.0, phi_max=0.375):
    """The state law: k_eff against SSA, both taken relative to their values at
    theta_ref (after the width-dependent transient). A power law with one
    exponent for every well-connected packing; porosity sets the level."""
    phis = sorted({p for (p, T) in R})
    col = colours(phis, "phi")
    per = defaultdict(list)
    fig, ax = plt.subplots(1, 3, figsize=(W, 62 * MM))
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.24, top=0.90, wspace=0.50)
    for (p, T), d in sorted(R.items()):
        for seed, r in d.items():
            th = r["th"]
            if th[-1] < 1.5 * theta_ref:
                continue
            kref = np.interp(theta_ref, th, r["kiso"]); sref = np.interp(theta_ref, th, r["ssa"])
            w = th >= theta_ref
            x, y = r["ssa"][w] / sref, r["kiso"][w] / kref
            per[(p, seed)].append((np.log10(x), np.log10(y), kref))
            if p <= phi_max:
                ax[0].plot(x, y, color=col[p], lw=0.8, alpha=0.85)
    p_by, k_by = defaultdict(list), defaultdict(list)
    for (p, seed), L in per.items():
        X = np.concatenate([l[0] for l in L]); Y = np.concatenate([l[1] for l in L])
        p_by[p].append(float(np.sum(X * Y) / np.sum(X * X)))
        k_by[p].append(float(np.mean([l[2] for l in L])))
    good = [p for p in phis if p <= phi_max]
    p0 = np.mean(np.concatenate([p_by[p] for p in good]))
    xx = np.linspace(0.70, 1.0, 30)
    ax[0].plot(xx, xx ** p0, color=INK, lw=1.3, ls="--")
    ax[0].text(0.05, 0.95, rf"$(\mathrm{{SSA}}/\mathrm{{SSA}}_\mathrm{{r}})^{{{p0:.2f}}}$",
               transform=ax[0].transAxes, fontsize=FS_S, va="top", ha="left")
    ax[0].set(xscale="log", yscale="log", xlabel=r"SSA$\,/\,$SSA$_\mathrm{r}$",
              ylabel=r"$k_\mathrm{eff}/k_\mathrm{eff,r}$")
    ax[0].invert_xaxis()
    ax[0].set_xticks([1.0, 0.9, 0.8, 0.7]); ax[0].set_xticklabels(["1.0", "0.9", "0.8", "0.7"])
    ax[0].set_yticks([1.0, 1.1, 1.2, 1.3]); ax[0].set_yticklabels(["1.0", "1.1", "1.2", "1.3"])
    for axis in (ax[0].xaxis, ax[0].yaxis):
        axis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    for p in good:
        ax[0].plot([], [], color=col[p], lw=1.6, label=rf"$\varphi={p:g}$")
    ax[0].legend(frameon=False, fontsize=FS_S, loc="lower right", handlelength=1.2, labelspacing=0.2)
    cf = np.polyfit(np.concatenate([[p] * len(k_by[p]) for p in good]),
                    np.log(np.concatenate([k_by[p] for p in good]) / K_ICE), 1)
    for p in phis:
        ax[1].scatter([p] * len(k_by[p]), np.array(k_by[p]) / K_ICE, s=14, color=col[p],
                      edgecolor="white", linewidth=0.4, zorder=3)
        ax[2].scatter([p] * len(p_by[p]), p_by[p], s=14, color=col[p], edgecolor="white",
                      linewidth=0.4, zorder=3)
    pp = np.linspace(min(phis), phi_max, 20)
    ax[1].plot(pp, np.exp(np.polyval(cf, pp)), color=INK, lw=1.2, ls="--")
    ax[2].axhline(p0, color=INK, lw=1.0, ls="--")
    for a_ in ax[1:]:
        a_.axvspan(phi_max + 0.025, max(phis) + 0.025, color="#efefef", lw=0, zorder=0)
        a_.set_xticks(phis); a_.set_xticklabels([f"{p:g}" for p in phis], rotation=45)
        a_.set_xlim(min(phis) - 0.025, max(phis) + 0.025)
        a_.set_xlabel(r"$\varphi$")
    ax[1].set(ylabel=r"$k_\mathrm{eff,r}/k_\mathrm{ice}$", yscale="log")
    ax[1].set_yticks([0.1, 0.2, 0.3]); ax[1].set_yticklabels(["0.1", "0.2", "0.3"])
    ax[1].yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax[2].set(ylabel=r"exponent $p$")
    for a_, s_ in zip(ax, "abc"):
        clean(a_); panel(a_, s_, dx=-0.32)
    save(fig, out, "Figure7_state_law")
    print(f"  state law: k ~ SSA^{p0:.3f} (phi <= {phi_max}); level k_r/k_ice = "
          f"{np.exp(cf[1]):.3f} exp({cf[0]:.2f} phi); per-phi p: "
          + ", ".join(f"{p:g}: {np.mean(p_by[p]):.2f}±{np.std(p_by[p], ddof=1):.2f}" for p in phis))
    return p0


def fig_porosity(R, out, T=-20, theta_ref=30.0, phi_max=0.375, marks_d=(0.0, 10.0, 19.0, 30.0)):
    """Curves to go with the gallery: (a) k_eff(t) for the five porosities (seed
    mean, band = seed min-max) with the gallery's instants marked; (b) the
    state law k vs SSA; (c) its level; (d) its exponent."""
    phis = sorted({p for (p, TT) in R})
    col = colours(phis, "phi")
    fig, axg = plt.subplots(2, 2, figsize=(W, 120 * MM))
    fig.subplots_adjust(left=0.10, right=0.985, bottom=0.115, top=0.955, wspace=0.34, hspace=0.42)
    ax = axg.ravel()
    for p in phis:
        g = list(R[(p, T)].values())
        x, y, lo, hi = mean_curve(g, "td", "kiso", logx=False)
        ax[0].fill_between(x, lo, hi, color=col[p], alpha=0.15, lw=0)
        ax[0].plot(x, y, color=col[p], lw=1.6)
        ax[0].text(x[-1] + 0.7, y[-1], rf"${p:g}$", fontsize=FS_S, va="center", color=INK)
    for j, tm in enumerate(marks_d):
        ax[0].axvline(tm, color=MUTED, lw=0.6, ls=":")
    ax[0].text(30.7, 1.0, r"$\varphi$", fontsize=FS_S, va="center")
    ax[0].set(xlabel="$t$ [d]", ylabel=r"$k_\mathrm{eff}$ [W m$^{-1}$ K$^{-1}$]", xlim=(-1, 36))
    ax[0].set_xticks([0, 10, 20, 30])
    per = defaultdict(list)
    for (p, TT), d in sorted(R.items()):
        for seed, r in d.items():
            th = r["th"]
            if th[-1] < 1.5 * theta_ref:
                continue
            kref = np.interp(theta_ref, th, r["kiso"]); sref = np.interp(theta_ref, th, r["ssa"])
            w = th >= theta_ref
            x, y = r["ssa"][w] / sref, r["kiso"][w] / kref
            per[(p, seed)].append((np.log10(x), np.log10(y), kref))
            if p <= phi_max:
                ax[1].plot(x, y, color=col[p], lw=0.8, alpha=0.85)
    p_by, k_by = defaultdict(list), defaultdict(list)
    for (p, seed), L in per.items():
        X = np.concatenate([l[0] for l in L]); Y = np.concatenate([l[1] for l in L])
        p_by[p].append(float(np.sum(X * Y) / np.sum(X * X))); k_by[p].append(float(np.mean([l[2] for l in L])))
    good = [p for p in phis if p <= phi_max]
    p0 = np.mean(np.concatenate([p_by[p] for p in good]))
    xx = np.linspace(0.70, 1.0, 30)
    ax[1].plot(xx, xx ** p0, color=INK, lw=1.3, ls="--")
    ax[1].text(0.05, 0.95, rf"$(\mathrm{{SSA}}/\mathrm{{SSA}}_\mathrm{{r}})^{{{p0:.2f}}}$",
               transform=ax[1].transAxes, fontsize=FS_S, va="top")
    ax[1].set(xscale="log", yscale="log", xlabel=r"SSA$\,/\,$SSA$_\mathrm{r}$",
              ylabel=r"$k_\mathrm{eff}/k_\mathrm{eff,r}$")
    ax[1].invert_xaxis()
    ax[1].set_xticks([1.0, 0.9, 0.8, 0.7]); ax[1].set_xticklabels(["1.0", "0.9", "0.8", "0.7"])
    ax[1].set_yticks([1.0, 1.1, 1.2, 1.3]); ax[1].set_yticklabels(["1.0", "1.1", "1.2", "1.3"])
    for axis in (ax[1].xaxis, ax[1].yaxis):
        axis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    cf = np.polyfit(np.concatenate([[p] * len(k_by[p]) for p in good]),
                    np.log(np.concatenate([k_by[p] for p in good]) / K_ICE), 1)
    for p in phis:
        ax[2].scatter([p] * len(k_by[p]), np.array(k_by[p]) / K_ICE, s=16, color=col[p],
                      edgecolor="white", linewidth=0.4, zorder=3)
        ax[3].scatter([p] * len(p_by[p]), p_by[p], s=16, color=col[p], edgecolor="white",
                      linewidth=0.4, zorder=3)
    pp = np.linspace(min(phis), phi_max, 20)
    ax[2].plot(pp, np.exp(np.polyval(cf, pp)), color=INK, lw=1.2, ls="--")
    ax[3].axhline(p0, color=INK, lw=1.0, ls="--")
    for a_ in ax[2:]:
        a_.axvspan(phi_max + 0.025, max(phis) + 0.025, color="#efefef", lw=0, zorder=0)
        a_.set_xticks(phis); a_.set_xticklabels([f"{p:g}" for p in phis])
        a_.set_xlim(min(phis) - 0.025, max(phis) + 0.025); a_.set_xlabel(r"$\varphi$")
    ax[2].set(ylabel=r"$k_\mathrm{eff,r}/k_\mathrm{ice}$", yscale="log")
    ax[2].set_yticks([0.1, 0.2, 0.3]); ax[2].set_yticklabels(["0.1", "0.2", "0.3"])
    ax[2].yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax[3].set(ylabel=r"exponent $p$")
    for a_, s_ in zip(ax, "abcd"):
        clean(a_); panel(a_, s_, dx=-0.20, dy=1.01)
    save(fig, out, "Figure7_state_law_alt")


def psat_ice(T):                                    # Murphy & Koop (2005), Pa
    return np.exp(9.550426 - 5723.265 / T + 3.53068 * np.log(T) - 0.00728332 * T)


def fig_timescales(out, theta_target=331.0, tau_ref=7822.3, T_ref=253.15, R_ref=50e-6):
    """Time to reach theta_target, tau_sub ~ R^2 sqrt(T)/rho_vs(T), anchored on the
    campaign's tau_sub at -20 C and R = 50 um (alpha_c = 1e-3).

    Nothing is written ON the coloured field except two marker labels in white
    boxes: the time classes are read from a discrete colour bar, and the
    Enceladus ranges are brackets OUTSIDE the axes (top: temperature; right:
    plume grain sizes), so every label is dark text on white."""
    from matplotlib.colors import BoundaryNorm, ListedColormap
    rate = lambda T: psat_ice(T) / T / np.sqrt(T)
    T = np.linspace(60, 260, 400); Rg = np.geomspace(0.1e-6, 1e-3, 400)
    TT, RR = np.meshgrid(T, Rg)
    t = theta_target * tau_ref * (rate(T_ref) / rate(TT)) * (RR / R_ref) ** 2
    lev = [3600.0, DAY, YEAR, 1e3 * YEAR, 1e6 * YEAR, 4.5e9 * YEAR]
    lab = ["1 h", "1 d", "1 yr", "1 kyr", "1 Myr", "4.5 Gyr"]
    cols = [cmocean.cm.matter(x) for x in np.linspace(0.04, 0.96, len(lev) + 1)]
    idx = np.digitize(np.log10(t), np.log10(lev))            # class 0..6
    # One AGU column: figures are 85 or 170 mm wide, nothing in between.
    fig = plt.figure(figsize=(85 * MM, 72 * MM))
    ax = fig.add_axes([0.175, 0.17, 0.50, 0.665])
    cmap = ListedColormap(cols)
    im = ax.pcolormesh(TT, RR * 1e6, idx, cmap=cmap, norm=BoundaryNorm(np.arange(-0.5, len(lev) + 1), cmap.N),
                       shading="auto", rasterized=True)
    ax.contour(TT, RR * 1e6, t, levels=lev, colors="white", linewidths=0.6)
    ax.set(yscale="log", xlabel="$T$ [K]", ylabel=r"$R$ [$\mu$m]", xlim=(60, 260), ylim=(0.1, 1e3))
    box = dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.9)
    ax.plot([253.15], [50], "o", ms=5, mfc="white", mec=INK, mew=1.0)
    ax.annotate("this study", (253.15, 50), xytext=(-7, 9), textcoords="offset points",
                ha="right", fontsize=FS_S, color=INK, bbox=box)
    ax.plot([180], [6], "s", ms=5, mfc="white", mec=INK, mew=1.0)
    ax.annotate("Choukroun et al.\n(2020): 15 yr", (180, 6), xytext=(-9, -5), textcoords="offset points",
                ha="right", va="top", fontsize=FS_S, color=INK, bbox=box, multialignment="left")
    # brackets outside the axes
    tr = ax.get_xaxis_transform()
    for a0, a1, s_ in ((60, 80, "surface"), (175, 185, "fractures")):
        ax.plot([a0, a1], [1.035, 1.035], color=INK, lw=2.2, transform=tr, clip_on=False,
                solid_capstyle="butt")
        ax.text((a0 + a1) / 2, 1.06, s_, transform=tr, ha="center", va="bottom", fontsize=FS_S, color=INK)
        for x_ in (a0, a1):
            ax.axvline(x_, color="white", lw=0.5, ls=(0, (2, 2)), alpha=0.7)
    ty = ax.get_yaxis_transform()
    ax.plot([1.03, 1.03], [0.1, 5], color=INK, lw=2.2, transform=ty, clip_on=False, solid_capstyle="butt")
    ax.text(1.055, np.sqrt(0.1 * 5), "plume\ngrains", transform=ty, ha="left", va="center",
            fontsize=FS_S, color=INK)
    ax.axhline(5, color="white", lw=0.5, ls=(0, (2, 2)), alpha=0.7)
    # discrete colour bar: the time classes
    cax = fig.add_axes([0.80, 0.17, 0.035, 0.665])
    cb = fig.colorbar(im, cax=cax, ticks=np.arange(0.5, len(lev)))
    cb.ax.set_yticklabels(lab, fontsize=FS_S); cb.ax.tick_params(length=2, width=0.5)
    cb.outline.set_linewidth(0.5)
    cb.ax.set_title(r"$t$", fontsize=FS, pad=4)
    ax.set_xticks([100, 150, 200, 250])
    clean(ax)
    save(fig, out, "Figure8_timescales")


def fig_gallery(camp, out, T=-20):
    """What happens to the aggregates: five porosities (columns) at four
    instants (rows). A picture for the reader, so the colour bars are
    qualitative: ice, and which way the vapour is driving the surface."""
    from plot_keff_snapshots import make_reader, snap_step, _field, _scalebar, WANT, SIGMA_SCALE, \
        centered_cmap, ice_alpha_cmap
    from pplib import step_times, opening_step
    from matplotlib.cm import ScalarMappable
    seeds = {0.275: 1602, 0.325: 1702, 0.375: 1802, 0.425: 1902, 0.475: 2002}
    frames, times = {}, None
    for p, s_ in seeds.items():
        d = next(camp.glob(f"packing_2D_phi{p}_*seed{s_}_L2mm_eps1000nm_perxy_T{T}__*"))
        files, reader = make_reader(d, "sol")
        tm = step_times(str(d)); st = [snap_step(f) for f in files]
        tt = np.array([tm.get(x, np.nan) for x in st])
        op = opening_step(st, list(tt)); i0 = st.index(op) if op is not None else 0
        te = max(tm.values())
        pick = [i0] + [int(np.nanargmin(np.abs(tt - (tt[i0] + f * (te - tt[i0]))))) for f in (1 / 3, 2 / 3)] \
            + [len(st) - 1]
        frames[p] = [reader(files[i], want=WANT) for i in pick]
        if times is None:
            times = [tt[i] / DAY for i in pick]
    pore = []
    for p in frames:
        for fl, X, Y in frames[p]:
            sg = SIGMA_SCALE * pplib.supersaturation(fl["VaporDensity"], fl["Temperature"])
            pore.append(sg[fl["IcePhase"] < 0.5])
    pore = np.concatenate(pore); v = min(abs(pore.min()), abs(pore.max()))
    norm = AsinhNorm(linear_width=max(v / 300, 1e-12), vmin=-v, vmax=v)
    vapcm, icecm = centered_cmap(cmocean.cm.balance, norm), ice_alpha_cmap()
    n, nr = len(frames), 4
    Wmm, Lm, Rm, gap, top, bot = 170.0, 9.0, 1.0, 1.2, 5.5, 17.0
    cell = (Wmm - Lm - Rm - (n - 1) * gap) / n
    Hmm = top + nr * cell + (nr - 1) * gap + bot
    fig = plt.figure(figsize=(Wmm * MM, Hmm * MM))
    for j, p in enumerate(sorted(frames)):
        for i in range(nr):
            ax = fig.add_axes([(Lm + j * (cell + gap)) / Wmm,
                               (bot + (nr - 1 - i) * (cell + gap)) / Hmm, cell / Wmm, cell / Hmm])
            fl, X, Y = frames[p][i]
            XX, YY = _field(ax, fl, X, Y, norm, vapcm, icecm)
            if i == nr - 1 and j == 0:
                _scalebar(ax, XX, YY)
            if i == 0:
                ax.set_title(rf"$\varphi={p:g}$", fontsize=FS, pad=3)
            if j == 0:
                ax.set_ylabel("0 d" if i == 0 else f"{times[i]:.0f} d", fontsize=FS, labelpad=3)
    # qualitative colour bars
    cax1 = fig.add_axes([(Lm + 6) / Wmm, 10.0 / Hmm, 32 / Wmm, 2.2 / Hmm])
    cb = fig.colorbar(ScalarMappable(cmap=cmocean.cm.ice, norm=plt.Normalize(0, 1)), cax=cax1,
                      orientation="horizontal", ticks=[0, 1])
    cb.ax.set_xticklabels(["air", "ice"], fontsize=FS_S); cb.outline.set_linewidth(0.5)
    cb.ax.tick_params(length=0, pad=2)
    cax2 = fig.add_axes([(Lm + 70) / Wmm, 10.0 / Hmm, 70 / Wmm, 2.2 / Hmm])
    cb = fig.colorbar(ScalarMappable(cmap=vapcm, norm=norm), cax=cax2, orientation="horizontal",
                      ticks=[-v, 0, v])
    cb.ax.set_xticklabels(["undersaturated\n(ice sublimates)", "equilibrium",
                           "supersaturated\n(vapor desublimates)"], fontsize=FS_S)
    cb.minorticks_off()
    cb.outline.set_linewidth(0.5); cb.ax.tick_params(length=0, pad=2)
    fig.text((Lm + 67) / Wmm, 11.1 / Hmm, "vapor", fontsize=FS_S, ha="right", va="center")
    save(fig, out, "Figure5_gallery")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root", type=Path)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--skip-gallery", action="store_true")
    a = ap.parse_args()
    out = a.out or a.root / "compare" / "figure_samples"
    out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(pplib.MANUSCRIPT_RC)
    R = defaultdict(dict)
    for kf in sorted(a.root.glob("packing_*/k_eff.csv")):
        m = PAT.search(kf.parent.name)
        if not m or not 1601 <= int(m.group(2)) <= 2005:
            continue
        r = load(kf.parent)
        r["tau"] = read_tau_sub(kf.parent)
        r["i0"] = int(np.argmax(r["t"] >= 1.0))
        r["th"] = r["t"] / r["tau"]; r["td"] = r["t"] / DAY
        R[(float(m.group(1)), int(m.group(3)))][int(m.group(2))] = r
    fig_collapse(R, out)
    p0 = fig_closure(R, out)
    fig_porosity(R, out)
    fig_timescales(out)
    if not a.skip_gallery:
        fig_gallery(a.root, out)
    print(f"state-law exponent p = {p0:.3f}")


if __name__ == "__main__":
    main()
