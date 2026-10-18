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

Further options, each in its own sub-folder of --out (see FOLDER):

  option3_heat_flow        the local heat flux for a gradient along x at four
                           instants of the -20 C run, each labelled with its
                           conductivity; with and without the ice row
  option4_sintering_clock  five temperatures against time, against sintering
                           age t/tau, and the factor between their clocks
                           (simulated temperatures only -- no extrapolation)
  option5_porosity         one packing per porosity after 30 days, k_eff(t)
                           by porosity (seed mean, band = seed range), and
                           k_eff at 30 days against porosity
  option6_validation       two-grain neck growth against Molaro et al. (2019)

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
import matplotlib.patches, matplotlib.ticker, matplotlib.colors  # noqa: E401
from matplotlib.cm import ScalarMappable
from matplotlib.path import Path as MPath
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
C_SEED_A, C_SEED_B = "#7B3FA0", "#1B9E77"    # two packings at one temperature: not warm/cold hues
FS, FS_S, FS_T = 15, 13, 12                  # labels | ticks, legends | small (poster captions are 16 pt Arial)
# Output sub-folder of each figure kind: one folder per option.
FOLDER = {
    "option1_stacked": "option1_ssa_collapse", "option1_compact": "option1_ssa_collapse",
    "grid_only": "option1_ssa_collapse", "curves_ssa": "option1_ssa_collapse",
    "option2_compact": "option2_state_law", "curves_law": "option2_state_law",
    "heat_flow_row": "option3_heat_flow", "heat_flow_with_microstructure": "option3_heat_flow",
    "sintering_clock": "option4_sintering_clock",
    "porosity": "option5_porosity",
    "validation_molaro": "option6_validation",
}
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
def num_marker(ax, x, y, n, color, ms=15, color2=None):
    """Numbered ring. With color2 the ring is split: left half color, right half color2."""
    ax.plot([x], [y], "o", ms=ms, mfc="white", mec="none" if color2 else color, mew=2.0, zorder=6,
            clip_on=False)
    if color2:
        for (a0, a1), c in (((90, 270), color), ((-90, 90), color2)):
            ax.plot([x], [y], marker=MPath.arc(a0, a1), ms=ms, mfc="none", mec=c, mew=2.0, zorder=6,
                    clip_on=False, ls="none")
    ax.text(x, y, str(n), ha="center", va="center", fontsize=FS_T - 1, fontweight="bold",
            color=INK, zorder=7, clip_on=False)


def pair_markers(ax, runs, hot, cold, marks, xkey, xtol, ytol=0.005, hl=None):
    """Numbered instants of the two highlighted runs. Where instant n of the warm
    run and of the cold run fall on the same point (within xtol, ytol in data
    units) it is drawn once, ringed half in each colour."""
    rh, rc = runs[hot], runs[cold]
    cw, cc = (hl[hot][1], hl[cold][1]) if hl else (C_WARM, C_COLD)
    for n, (ih, ic) in enumerate(zip(marks[hot], marks[cold]), 1):
        ph, pc = (rh[xkey][ih], rh["kn"][ih]), (rc[xkey][ic], rc["kn"][ic])
        if abs(ph[0] - pc[0]) <= xtol and abs(ph[1] - pc[1]) <= ytol:
            num_marker(ax, 0.5 * (ph[0] + pc[0]), 0.5 * (ph[1] + pc[1]), n, cw, color2=cc)
        else:
            num_marker(ax, *ph, n, cw)
            num_marker(ax, *pc, n, cc)


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


def panel_time(ax, runs, hot, cold, marks, short=False, hl=None):
    """k_eff/k_0 against time: the temperatures fan out. hl = {key: (label, colour)}
    replaces the temperature labels and colours (seed mode); the grey curves are
    then left unlabelled."""
    def end_label(r, T, col, bold, dx):
        if hl and T not in hl:
            return
        ax.annotate(hl[T][0] if hl else tlabel(T), (r["td"][-1], r["kn"][-1]), xytext=(dx, 0),
                    textcoords="offset points", color=col, va="center", ha="left",
                    fontsize=FS_S if bold else FS_T, fontweight="bold" if bold else "normal",
                    annotation_clip=False)
    for T, r in sorted(runs.items()):
        if T in (hot, cold):
            continue
        ax.plot(r["td"], r["kn"], color=GREY, lw=1.4, zorder=2)
        end_label(r, T, MUTED, False, 6)
    for T, col in ((hot, hl[hot][1] if hl else C_WARM), (cold, hl[cold][1] if hl else C_COLD)):
        r = runs[T]
        ax.plot(r["td"], r["kn"], color=col, lw=3.2, zorder=4)
        end_label(r, T, col, True, 14)                  # clear of the end marker
    pair_markers(ax, runs, hot, cold, marks, "td", xtol=0.15, hl=hl)
    # the axis stops at the data; the labels sit in the margin to its right
    ax.set_xlim(0, 30.3); ax.set_xticks([0, 10, 20, 30])
    ax.set_xlabel("time  [days]", fontsize=FS)
    ax.set_ylabel(K_LABEL[short], fontsize=FS)
    clean(ax)


def panel_ssa(ax, runs, hot, cold, marks, short=False, hl=None):
    """The same curves against SSA/SSA_0: one path, travelled at different speeds."""
    for T, r in sorted(runs.items()):
        if T not in (hot, cold):
            ax.plot(r["sn"], r["kn"], color=GREY, lw=1.4, zorder=2)
    cw, cc = (hl[hot][1], hl[cold][1]) if hl else (C_WARM, C_COLD)
    ax.plot(runs[hot]["sn"], runs[hot]["kn"], color=cw, lw=5.0, zorder=3, solid_capstyle="round")
    ax.plot(runs[cold]["sn"], runs[cold]["kn"], color=cc, lw=2.2, zorder=4, solid_capstyle="round")
    pair_markers(ax, runs, hot, cold, marks, "sn", xtol=0.005, hl=hl)
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
# further options: heat flow, sintering clock, porosity, validation
# ---------------------------------------------------------------------------
K_ICE, K_AIR = 2.29, 0.02
C_RUN = C_COLD                         # single-run figures: the -20 C master run


def heat_flux(run_dir, step):
    """(|q| per unit mean gradient [W/m/K], phi, x_um, y_um) for a unit
    macroscopic gradient along x, from the corrector the solver wrote
    (-keff_write_corrector). Same tensor law as the solver, as in
    fig_methods_homogenization.py."""
    from igakit.io import PetIGA
    files, reader = make_reader(run_dir, "sol")
    fn = [f for f in files if snap_step(f) == step][0]
    fl, X, Y = reader(fn, want=("IcePhase",))
    phi = fl["IcePhase"]; x = X[0] * 1e6; y = Y[:, 0] * 1e6
    cdir = run_dir / "corrector"
    nrb = PetIGA().read(str(cdir / "igakeff.dat"))
    cor = np.squeeze(PetIGA().read_vec(str(cdir / f"t_vec_{step:05d}_0.dat"), nrb)).T
    gy, gx = np.gradient(cor, y * 1e-6, x * 1e-6)
    pc = np.clip(phi, 0.0, 1.0)
    ka = pc * K_ICE + (1 - pc) * K_AIR
    kh = 1.0 / (pc / K_ICE + (1 - pc) / K_AIR)
    py_, px_ = np.gradient(pc, y * 1e-6, x * 1e-6)
    g2 = px_ ** 2 + py_ ** 2
    Ex, Ey = 1.0 + gx, gy
    with np.errstate(invalid="ignore", divide="ignore"):
        proj = np.where(g2 > 0, (px_ * Ex + py_ * Ey) / g2, 0.0)
    qx = ka * Ex + (kh - ka) * proj * px_
    qy = ka * Ey + (kh - ka) * proj * py_
    return np.hypot(qx, qy), phi, x, y


def build_heat_flow(R, a, with_micro):
    """Where the heat goes: the local flux for a temperature gradient along x,
    at four instants of the master run (-20 C). Each panel carries the
    conductivity of that instant (k_xx of the same solve)."""
    run = R[(a.phi, -20)][a.seed]
    rd = run["dir"]
    rep_ = np.genfromtxt(rd / "corrector" / "k_eff_replay.csv", delimiter=",", names=True)
    steps = [int(v) for v in np.atleast_1d(rep_["step"])]
    st = a.stride
    Q = [heat_flux(rd, s_) for s_ in steps]
    vmax = float(np.percentile(np.concatenate([q[0][::8, ::8].ravel() for q in Q]), 99.9))
    norm = matplotlib.colors.LogNorm(vmin=K_AIR, vmax=vmax)
    n = len(steps); W = a.width_in; gap, cbw = 0.09, 1.15
    s = (W - 0.10 - 0.34 - cbw - (n - 1) * gap) / n
    top, bot = 0.42, 0.50
    H = top + (2 * s + gap if with_micro else s) + bot
    fig = plt.figure(figsize=(W, H)); F = lambda x, y, w, h: (x / W, y / H, w / W, h / H)
    y_flux = bot
    for i, ((q, phi, x, y), step) in enumerate(zip(Q, steps)):
        t_d = float(rep_["time"][i]) / DAY
        ext = (x.min(), x.max(), y.min(), y.max())
        x0 = 0.10 + 0.34 + i * (s + gap)
        if with_micro:
            axm = fig.add_axes(F(x0, y_flux + s + gap, s, s))
            axm.imshow(phi[::st, ::st], origin="lower", extent=ext, cmap=cmocean.cm.ice,
                       vmin=0, vmax=1, interpolation="antialiased")
        ax = fig.add_axes(F(x0, y_flux, s, s))
        im = ax.imshow(q[::st, ::st], origin="lower", extent=ext, cmap=cmocean.cm.thermal,
                       norm=norm, interpolation="antialiased")
        for ax_ in ([axm, ax] if with_micro else [ax]):
            ax_.set_xticks([]); ax_.set_yticks([])
            for sp in ax_.spines.values():
                sp.set_linewidth(3.2); sp.set_color(C_RUN)
        (axm if with_micro else ax).set_title(
            f"{i + 1}   " + ("start" if t_d < 0.01 else f"{t_d:.0f} days"), fontsize=FS, color=INK, pad=6)
        ax.text(0.5, 0.035, rf"$k_\mathrm{{eff}}$ = {float(rep_['k_00'][i]):.2f} W m$^{{-1}}$ K$^{{-1}}$",
                transform=ax.transAxes, ha="center", va="bottom", fontsize=FS_T, color=INK,
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.88))
        if i == 0:
            x_a, y_a = x.min() + 0.07 * np.ptp(x), y.max() - 0.09 * np.ptp(y)
            for lw, c in ((6.0, "white"), (3.2, INK)):
                ax.plot([x_a, x_a + 500], [y_a, y_a], color=c, lw=lw, solid_capstyle="butt")
            ax.text(x_a + 250, y_a - 0.035 * np.ptp(y), "0.5 mm", ha="center", va="top",
                    fontsize=FS_T, color=INK,
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
    if with_micro:
        fig.text(*F(0.10 + 0.11, y_flux + 1.5 * s + gap, 0, 0)[:2], "ice", rotation=90, ha="center",
                 va="center", fontsize=FS, fontweight="bold", color=C_RUN)
    fig.text(*F(0.10 + 0.11, y_flux + 0.5 * s, 0, 0)[:2], "heat flow", rotation=90, ha="center",
             va="center", fontsize=FS, fontweight="bold", color=C_RUN)
    cax = fig.add_axes(F(0.10 + 0.34 + n * s + (n - 1) * gap + 0.22, y_flux + 0.03 * s, 0.16, 0.94 * s))
    cb = fig.colorbar(im, cax=cax); cb.ax.tick_params(labelsize=FS_S, width=0.8, length=3)
    cb.set_label(r"heat flux / gradient  [W m$^{-1}$ K$^{-1}$]", fontsize=FS_T)
    cb.outline.set_linewidth(0.8)
    # the imposed gradient
    xa = 0.10 + 0.34
    fig.patches.append(matplotlib.patches.FancyArrowPatch(
        (xa / W, 0.20 / H), ((xa + 1.5) / W, 0.20 / H), transform=fig.transFigure,
        arrowstyle="-|>", mutation_scale=16, lw=1.6, color=MUTED))
    fig.text((xa + 1.62) / W, 0.20 / H,
             f"imposed temperature gradient   ·   −20 °C   ·   porosity {a.phi:g}   ·   2 mm × 2 mm periodic cell",
             fontsize=FS_T, color=MUTED, ha="left", va="center")
    return fig, (W, H)


def build_clock(R, a):
    """Temperature sets the clock: the same curves against time, against
    sintering age t/tau, and the factor between the clocks."""
    runs = {T: R[(a.phi, T)][a.seed] for (p, T) in R if p == a.phi and a.seed in R[(p, T)]}
    Ts = sorted(runs)
    cm, (c0, c1) = CMAP["T"]
    col = {T: cm(c0 + (c1 - c0) * i / max(1, len(Ts) - 1)) for i, T in enumerate(Ts)}
    W, ph, bot, top = a.width_in, 2.75, 0.70, 0.22
    H = top + ph + bot
    fig = plt.figure(figsize=(W, H)); F = lambda x, y, w, h: (x / W, y / H, w / W, h / H)
    l, mid, r_ = 0.95, 1.05, 0.25
    pw = (W - l - 2 * mid - r_ - 0.75) / 3
    ax1 = fig.add_axes(F(l, bot, pw, ph))
    ax2 = fig.add_axes(F(l + pw + 0.75 + mid, bot, pw, ph))
    ax3 = fig.add_axes(F(l + 2 * pw + 0.75 + 2 * mid, bot, pw, ph))
    for T in Ts:
        r = runs[T]
        ax1.plot(r["td"], r["kn"], color=col[T], lw=2.6)
        ax1.annotate(tlabel(T), (r["td"][-1], r["kn"][-1]), xytext=(6, 0), textcoords="offset points",
                     color=col[T], fontsize=FS_T, fontweight="bold", va="center", annotation_clip=False)
        w = r["th"] >= 1.0
        ax2.plot(r["th"][w], r["kn"][w], color=col[T], lw=2.6)
    ax1.set_xlim(0, 30.3); ax1.set_xticks([0, 10, 20, 30])
    ax1.set_xlabel("time  [days]", fontsize=FS); ax1.set_ylabel(K_LABEL[False], fontsize=FS)
    ax2.set_xscale("log"); ax2.set_xlabel(r"sintering age  $t\,/\,\tau(T)$", fontsize=FS)
    ax2.set_ylabel(K_LABEL[True], fontsize=FS)
    ref = max(Ts)                                  # the warmest run sets the unit
    fac = [runs[T]["tau"] / runs[ref]["tau"] for T in Ts]
    ax3.plot(Ts, fac, color=GREY, lw=1.6, zorder=2)
    for T, f in zip(Ts, fac):
        ax3.plot([T], [f], "o", ms=11, color=col[T], mec="white", mew=1.0, zorder=3)
        ax3.annotate(f"{f:.0f}×" if f >= 10 else f"{f:.1f}×", (T, f), xytext=(9, 4),
                     textcoords="offset points", fontsize=FS_T, color=INK)
    ax3.set_yscale("log"); ax3.set_xticks(Ts); ax3.set_xticklabels([tlabel(T).replace(" °C", "") for T in Ts])
    ax3.set_xlim(min(Ts) - 4, max(Ts) + 7)
    ax3.set_xlabel("temperature  [°C]", fontsize=FS)
    ax3.set_ylabel(f"days to match 1 day at {tlabel(ref)}", fontsize=FS)
    ax3.set_yticks([1, 3, 10, 30]); ax3.set_yticklabels(["1", "3", "10", "30"])
    ax3.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    for ax in (ax1, ax2, ax3):
        clean(ax)
    ax3.tick_params(axis="x", labelsize=FS_T - 1)   # -10 and -5 sit close together
    print("  clock factors: " + ", ".join(f"{tlabel(T)} {f:.2f}x" for T, f in zip(Ts, fac)))
    return fig, (W, H)


def build_porosity(R, a, T=-20):
    """Porosity sets the level: one packing per porosity after 30 days, and
    the conductivity of all five packings per porosity."""
    phis = sorted({p for (p, T_) in R if T_ == T})
    cm, (c0, c1) = CMAP["phi"]
    col = {p: cm(c0 + (c1 - c0) * i / max(1, len(phis) - 1)) for i, p in enumerate(phis)}
    n = len(phis); W = a.width_in; gap = 0.09
    s = (W - 0.20 - (n - 1) * gap) / n
    top, mid, ph, bot = 0.42, 0.80, 2.55, 0.70
    H = top + s + mid + ph + bot
    fig = plt.figure(figsize=(W, H)); F = lambda x, y, w, h: (x / W, y / H, w / W, h / H)
    y_s = bot + ph + mid
    for i, p in enumerate(phis):
        seeds = sorted(R[(p, T)])
        seed = a.seed if a.seed in seeds else seeds[0]
        run = R[(p, T)][seed]
        fl, X, Y, t, row = snapshots(run, (30.0,), a.stride)[0]
        ax = fig.add_axes(F(0.10 + i * (s + gap), y_s, s, s))
        ax.imshow(fl["IcePhase"], origin="lower", cmap=cmocean.cm.ice, vmin=0, vmax=1,
                  extent=(X.min() * 1e3, X.max() * 1e3, Y.min() * 1e3, Y.max() * 1e3),
                  interpolation="antialiased")
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_linewidth(3.2); sp.set_color(col[p])
        ax.set_title(f"porosity {p:g}", fontsize=FS, color=INK, pad=6)
        ax.text(0.5, 0.035, rf"$k_\mathrm{{eff}}$ = {run['kiso'][-1]:.2f} W m$^{{-1}}$ K$^{{-1}}$",
                transform=ax.transAxes, ha="center", va="bottom", fontsize=FS_T, color=INK,
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.88))
        if i == 0:
            for lw, c in ((6.0, "white"), (3.2, INK)):
                ax.plot([0.14, 0.64], [1.83, 1.83], color=c, lw=lw, solid_capstyle="butt")
            ax.text(0.39, 1.77, "0.5 mm", ha="center", va="top", fontsize=FS_T, color=INK,
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
    fig.text(0.10 / W, (y_s - 0.24) / H,
             f"after 30 days at {tlabel(T)}   ·   2 mm × 2 mm periodic cells   ·   ice light, air dark",
             fontsize=FS_T, color=MUTED, ha="left", va="center")
    l, midp, r_ = 0.95, 1.45, 0.30
    pw = (W - l - midp - r_) / 2
    ax1 = fig.add_axes(F(l, bot, pw, ph)); ax2 = fig.add_axes(F(l + pw + midp, bot, pw, ph))
    for p in phis:
        g = list(R[(p, T)].values())
        td = np.linspace(max(r["td"][0] for r in g), min(r["td"][-1] for r in g), 300)
        K = np.array([np.interp(td, r["td"], r["kiso"]) for r in g])
        ax1.fill_between(td, K.min(0), K.max(0), color=col[p], alpha=0.22, lw=0)
        ax1.plot(td, K.mean(0), color=col[p], lw=2.6)
        ax1.annotate(f"{p:g}", (td[-1], K.mean(0)[-1]), xytext=(6, 0), textcoords="offset points",
                     color=col[p], fontsize=FS_T, fontweight="bold", va="center", annotation_clip=False)
        kend = [r["kiso"][-1] for r in g]
        ax2.plot([p] * len(kend), kend, "o", ms=9, color=col[p], mec="white", mew=0.9, zorder=3)
    means = [np.mean([r["kiso"][-1] for r in R[(p, T)].values()]) for p in phis]
    ax2.plot(phis, means, color=GREY, lw=1.6, zorder=2)
    ax1.set_xlim(0, 30.3); ax1.set_xticks([0, 10, 20, 30])
    ax1.set_xlabel("time  [days]", fontsize=FS)
    ax1.set_ylabel(r"conductivity  $k_\mathrm{eff}$  [W m$^{-1}$ K$^{-1}$]", fontsize=FS)
    ax1.set_ylim(0, None)
    ax2.set_xticks(phis); ax2.set_xticklabels([f"{p:g}" for p in phis])
    ax2.set_xlabel("porosity", fontsize=FS)
    ax2.set_ylabel(r"$k_\mathrm{eff}$ after 30 days  [W m$^{-1}$ K$^{-1}$]", fontsize=FS)
    ax2.set_ylim(0, None)
    for ax in (ax1, ax2):
        clean(ax)
    return fig, (W, H)


def build_validation(a):
    """The model against the laboratory: neck growth between two ice grains,
    Molaro et al. (2019), -20 C."""
    from plot_molaro_validation import load_series
    from plot_mechanisms_figure import outline
    run = a.molaro_run.resolve()
    s = load_series("T-20", run)
    m = s["model"]
    W, H = a.validation_width_in, 3.7
    fig = plt.figure(figsize=(W, H)); F = lambda x, y, w, h: (x / W, y / H, w / W, h / H)
    ax = fig.add_axes(F(0.95, 0.70, W - 0.95 - 0.25, H - 0.70 - 0.20))
    ax.plot(m["t"], m["w"], color=C_RUN, lw=3.2, zorder=3, label="model (vapor transport)")
    ax.errorbar(s["td"], s["wd"], yerr=[s["em"], s["ep"]], fmt="o", ms=8, mfc="white", mec=INK,
                mew=1.6, ecolor=INK, elinewidth=1.2, capsize=3, capthick=1.2, zorder=4,
                label="experiment, Molaro et al. (2019)", clip_on=False)
    ax.set_xlim(0, 82); ax.set_ylim(16, 70)       # room below the curve for the inset
    ax.set_xlabel("time  [min]", fontsize=FS); ax.set_ylabel("neck width  [µm]", fontsize=FS)
    ax.legend(fontsize=FS_T, frameon=False, loc="upper left", handlelength=1.8)
    clean(ax)
    # the grain pair, outline only, in the empty lower right
    loops, step = outline(run, a.molaro_step)
    allp = np.vstack(loops)
    ins = ax.inset_axes([0.56, 0.02, 0.42, 0.30])
    for ln in loops:
        ins.add_patch(matplotlib.patches.Polygon(ln, closed=True, fc="#dfeef3", ec=INK, lw=1.2))
    ins.set_xlim(allp[:, 0].min() - 8, allp[:, 0].max() + 8)
    ins.set_ylim(allp[:, 1].min() - 8, allp[:, 1].max() + 8)
    ins.set_aspect("equal"); ins.axis("off")
    ins.text(0.5, 1.02, "two ice grains, −20 °C", transform=ins.transAxes, ha="center", va="bottom",
             fontsize=FS_T, color=MUTED)
    return fig, (W, H)


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------
def build(kind, R, a, hot, cold, seed_T=None):
    """hot, cold: two temperatures of packing --seed; or, with seed_T, two SEEDS
    of porosity --phi at temperature seed_T (the other seeds drawn grey)."""
    if seed_T is None:
        runs = {T: R[(a.phi, T)][a.seed] for (p, T) in R if p == a.phi and a.seed in R[(p, T)]}
        hl = None
    else:
        runs = dict(R[(a.phi, seed_T)])
        hl = {hot: ("A", C_SEED_A), cold: ("B", C_SEED_B)}
    compact = "compact" in kind
    days = (0.0, 10.0, 30.0) if compact else (0.0, 10.0, 19.4, 30.0)
    S = {T: snapshots(runs[T], days, a.stride) for T in (hot, cold)}
    marks = {T: [s_[4] for s_ in S[T]] for T in S}
    W = a.width_in
    note = (f"identical packings at the start  ·  2 mm × 2 mm periodic cell  ·  "
            f"porosity {a.phi:g}  ·  mean grain radius 50 µm")
    if hl:
        note = (f"two different packings at {tlabel(seed_T)} (grey: {len(runs) - 2} more)  ·  "
                f"2 mm × 2 mm periodic cell  ·  porosity {a.phi:g}  ·  mean grain radius 50 µm")
    second = panel_law if "option2" in kind or "law" in kind else panel_ssa

    def grid(fig, F, x0, y_top, s, gap):
        norm, ext = sigma_norm(S.values())
        vapcm, icecm = centered_cmap(cmocean.cm.balance, norm), ice_alpha_cmap()
        rows = ([(f"packing {hl[k][0]}", hl[k][1], S[k]) for k in (hot, cold)] if hl else
                [(tlabel(hot), C_WARM, S[hot]), (tlabel(cold), C_COLD, S[cold])])
        y_bot, bw = snap_grid(fig, F, x0, y_top, s, gap, rows, norm, vapcm, icecm)
        return y_bot, bw, norm, vapcm, ext

    def second_panel(ax, short):
        if second is panel_law:
            panel_law(ax, R, a.seed, hot, cold, a.phi, short)
        else:
            panel_ssa(ax, runs, hot, cold, marks, short, hl=hl)

    if kind.startswith("curves"):                   # the two plots alone, side by side
        Wc, ph, bot, top = a.curves_width_in, 2.55, 0.68, 0.22
        H = top + ph + bot
        fig = plt.figure(figsize=(Wc, H)); F = lambda x, y, w, h: (x / Wc, y / H, w / Wc, h / H)
        l, mid, r_ = 0.95, 1.75, 0.42              # mid: end labels + next y label
        pw = (Wc - l - mid - r_) / 2
        panel_time(fig.add_axes(F(l, bot, pw, ph)), runs, hot, cold, marks, hl=hl)
        second_panel(fig.add_axes(F(l + pw + mid, bot, pw, ph)), False)
        return fig, (Wc, H)

    if kind == "grid_only" or kind == "option1_stacked":
        gap, cbw = 0.09, 1.55                       # colour bars in a column at the right
        s = (W - 0.10 - 0.34 - cbw - 3 * gap) / 4
        top = 0.42
        if kind == "grid_only":
            lower = 0.42
        else:
            l, mid, r_, ph, bot = 1.05, 1.95, 0.20, a.plot_h, 0.80
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
            panel_time(fig.add_axes(F(l, bot, pw, ph)), runs, hot, cold, marks, hl=hl)
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
    pw = W - xl - (0.78 * a.font_scale + 0.2)       # room for the curve-end labels
    b0, g_ = 0.66, 0.78
    ph = (H - b0 - g_ - 0.16) / 2
    panel_time(fig.add_axes(F(xl, b0 + ph + g_, pw, ph)), runs, hot, cold, marks, short=True, hl=hl)
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
    ap.add_argument("--pairs", nargs="*", default=["-5,-20", "-5,-40"],
                    help="warm,cold temperature pairs to build")
    ap.add_argument("--seed-pairs", nargs="*", default=["1702,1701"],
                    help="option 1 only: seedA,seedB pairs of porosity --phi at --seed-T")
    ap.add_argument("--seed-T", type=int, default=-20)
    ap.add_argument("--curves-width-in", type=float, default=8.3,
                    help="width of the curves-only pieces [in] (last year's plot slot)")
    ap.add_argument("--kinds", nargs="+", default=list(FOLDER),
                    help="which figures to build (default: all)")
    ap.add_argument("--molaro-run", type=Path, default=Path.home() / (
        "SimulationResults/HPC_results/enceladus_DSM/GrainPairSintering/"
        "batch_2026-09-08__17.20.46_molaro_T-20_round2/"
        "molaro_2D_L450x225um_eps0.12um_axisym_T-20pair_r14um_dom2__molaro_T-20_h0.99715_2h_a1e-1_dirichlet"),
        help="the -20 C two-grain run, for the validation figure")
    ap.add_argument("--molaro-step", type=int, default=236, help="outline snapshot (t* + 78 min)")
    ap.add_argument("--validation-width-in", type=float, default=5.6)
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--font-scale", type=float, default=1.0,
                    help="multiply every type size (poster body text is 18 pt: use 1.2)")
    ap.add_argument("--plot-h", type=float, default=2.75,
                    help="height of the two plots in the stacked figure [in]")
    a = ap.parse_args()
    global FS, FS_S, FS_T
    FS, FS_S, FS_T = (round(v * a.font_scale, 1) for v in (FS, FS_S, FS_T))
    plt.rcParams.update(POSTER_RC)
    R = load_campaign(a.root)

    def save(fig, size, kind, name):
        d = a.out / FOLDER[kind]
        d.mkdir(parents=True, exist_ok=True)
        for e in ("png", "pdf"):
            fig.savefig(d / f"{name}.{e}", dpi=a.dpi)
        plt.close(fig)
        print(f"  wrote {FOLDER[kind]}/{name}.png/.pdf   {size[0]:.2f} x {size[1]:.2f} in")

    single = {
        "heat_flow_row": lambda: build_heat_flow(R, a, False),
        "heat_flow_with_microstructure": lambda: build_heat_flow(R, a, True),
        "sintering_clock": lambda: build_clock(R, a),
        "porosity": lambda: build_porosity(R, a),
        "validation_molaro": lambda: build_validation(a),
    }
    for kind in a.kinds:
        if kind in single:
            fig, size = single[kind]()
            save(fig, size, kind, kind)
            continue
        for pair in a.pairs:
            hot, cold = (int(v) for v in pair.split(","))
            fig, size = build(kind, R, a, hot, cold)
            save(fig, size, kind, f"{kind}_T{hot}_T{cold}")
        if FOLDER[kind] != "option1_ssa_collapse":
            continue
        for pair in a.seed_pairs:                   # same porosity, two packings
            sa, sb = (int(v) for v in pair.split(","))
            fig, size = build(kind, R, a, sa, sb, seed_T=a.seed_T)
            save(fig, size, kind, f"{kind}_seed{sa}_seed{sb}_T{a.seed_T}")


if __name__ == "__main__":
    main()
