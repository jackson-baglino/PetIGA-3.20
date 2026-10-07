#!/usr/bin/env python3
"""fig_demmenie.py — manuscript figure: neck growth at saturation against t^(1/3).

    venv_enceladus/bin/python studies/molaro_2019/demmenie/fig_demmenie.py <run dir>
        [--relax-tau 11] [--name Figure4_saturated_neck_growth] [--out <dir>] [--copy-to <dir>]

Built to read as the companion of the Molaro figure
(postprocess/plot_molaro_validation.py, molaro_full), and from ITS helpers, so
the two cannot drift apart: the same colour-bar strip, ice over the
supersaturation sigma in the pore, circled instants, panel letters, type sizes.

  (a) the grain pair at the start and at the end of the run, instants 1-2:
      the computed quarter (one grain, half-plane) reflected across the
      mirror plane and the symmetry axis
  (b) neck width against time: the simulation (line), the free fit
      C (t + t0)^a and the one-third fit C (t + t0)^(1/3), both over the
      samples after the relaxation period; instants 1-2 marked
  (c) the free-fit exponent against the start of the fit window, with the
      range Demmenie et al. (2025) measured and the 1/3 law
  (d) the grain diameter D / D_0: the saturation check, on the scale of the
      Molaro figure's shrinkage panels

Fits and numbers come from analyze_demmenie.py, so the figure and the
diagnostic plots cannot disagree.
"""
from __future__ import annotations

import argparse, glob, shutil, sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm, Normalize
from matplotlib.ticker import MaxNLocator
import cmocean

HERE = Path(__file__).resolve().parent
PROJ = HERE.parents[2]
sys.path.insert(0, str(PROJ / "postprocess")); sys.path.insert(0, str(HERE))
import pplib  # noqa: E402
import plot_molaro_validation as pmv  # noqa: E402
from plot_keff_snapshots import (make_reader, snap_step, _field, _scalebar, _snap_title, _mark,  # noqa: E402
                                 WANT, INK, FS, FS_SMALL, MM)
from analyze_demmenie import analyse, _pl, _p3, DEMMENIE, HOUR  # noqa: E402

C_SIM, C_FREE, C_THIRD = pmv.C_T20, INK, pmv.C_T5


def section(run, fn):
    """All fields of the full pair: reflect across r = 0, then across z = 0."""
    _, reader = make_reader(run, "sol")
    fl, X, Y = pmv._mirror(*reader(fn, want=WANT))
    mz = lambda a_: np.hstack([a_[:, :0:-1], a_])
    return {k: mz(v) for k, v in fl.items()}, np.hstack([-X[:, :0:-1], X]), mz(Y)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run", type=Path)
    ap.add_argument("--relax-tau", type=float, default=11.0)
    ap.add_argument("--name", default="Figure4_saturated_neck_growth")
    ap.add_argument("--width-mm", type=float, default=170.0)
    ap.add_argument("--out", type=Path, default=PROJ / "studies/molaro_2019/manuscript")
    ap.add_argument("--copy-to", type=Path, default=None)
    a = ap.parse_args()
    plt.rcParams.update(pplib.MANUSCRIPT_RC)
    A = analyse(a.run, a.relax_tau)
    t, w, tk, f, h, s = A["t"], A["w"], A["tk"], A["free"], A["third"], A["scan"]

    # ---- sections: first and last snapshot, one crop and one stride ----
    files = sorted(glob.glob(str(a.run / "sol_*.dat")), key=snap_step)
    raw = [section(a.run, files[0]), section(a.run, files[-1])]
    box = pmv.auto_crop(*raw[-1])
    box = (-max(abs(box[0]), abs(box[1])), max(abs(box[0]), abs(box[1])), box[2])
    frac = (box[1] - box[0]) * pmv.UM / (raw[0][1].max() - raw[0][1].min())
    stride = max(1, int(np.ceil(frac * raw[0][1].shape[1] / pmv.MAX_PX)))
    sub = lambda a_: a_[::stride, ::stride]
    secs = [pmv.crop({k: sub(v) for k, v in fl.items()}, sub(X), sub(Y), box) for fl, X, Y in raw]
    pore = np.concatenate([pmv.SIGMA_SCALE * pplib.supersaturation(fl["VaporDensity"], fl["Temperature"])[fl["IcePhase"] < 0.5]
                           for fl, _, _ in secs])
    smin, smax = float(pore.min()), float(pore.max())
    # The Molaro figure uses a symmetric asinh bar sized to the SMALLER extreme,
    # which suits a field that changes sign over decades. Here the whole pore
    # is supersaturated against a flat surface (the walls sit at 1 + 2 d0/R)
    # and varies by a factor of four, so that rule would paint every pixel the
    # end colour. Same colours, same zero in the middle, but a LINEAR bar a
    # little past the far-field value, so the depletion at the neck shows.
    if smin * smax > 0:
        v = 1.25 * max(abs(smin), abs(smax))
        norm = Normalize(vmin=-v, vmax=v); ext = "neither"; vapcm = cmocean.cm.balance
    else:
        v = min(abs(smin), abs(smax)) or max(abs(smin), abs(smax))
        ext = {(True, True): "both", (True, False): "min", (False, True): "max", (False, False): "neither"}[(smin < -v, smax > v)]
        norm = AsinhNorm(linear_width=max(v / 300.0, 1e-12), vmin=-v, vmax=v)
        vapcm = pmv.centered_cmap(cmocean.cm.balance, norm)
    icecm = pmv.ice_alpha_cmap()
    print(f"sigma x{pmv.SIGMA_SCALE:g} in the shown pore: {smin:+.3g} .. {smax:+.3g}; bar +-{v:.3g} ({ext})")

    # ---- layout, in inches, as plot_molaro_validation.build_full ----
    W = a.width_mm * MM
    gap = 0.12
    ml, mr = 0.50, 0.08
    axw = W - ml - mr
    s_w = (axw - gap) / 2
    fl0, X0, Y0 = secs[0]
    s_h = s_w * (Y0.max() - Y0.min()) / (X0.max() - X0.min())
    cb_h, cb_lab, cb_gap, t_band, top = 0.07, 0.15, 0.06, 0.19, 0.05
    g_snap, ph, g_row, ph2, bot = 0.34, 1.75, 0.62, 1.55, 0.40
    pgap, ml2 = 0.78, ml + 0.17
    pw = (W - ml2 - mr - pgap) / 2
    H = top + cb_h + cb_lab + cb_gap + t_band + s_h + g_snap + ph + g_row + ph2 + bot
    fig = plt.figure(figsize=(W, H))
    F = lambda x0, y0, ww, hh: (x0 / W, y0 / H, ww / W, hh / H)
    y_neck = bot + ph2 + g_row
    y_snap = y_neck + ph + g_snap

    t_sec = [t[0], t[-1]]
    for i, (fl, X, Y) in enumerate(secs):
        axi = fig.add_axes(F(ml + i * (s_w + gap), y_snap, s_w, s_h))
        XX, YY = _field(axi, fl, X, Y, norm, vapcm, icecm)
        axi.set_aspect("auto")
        if i == 0:
            _scalebar(axi, XX, YY)
        _snap_title(axi, pmv.LETTERS[i], f"{t_sec[i] / HOUR:.0f} h")
    a.sig_extend = ext
    pmv._strip(fig, F, ml, y_snap + s_h + t_band + cb_gap + cb_lab, axw, norm, vapcm, ext)

    def dress(ax):
        ax.patch.set_alpha(0.0)
        ax.tick_params(labelsize=FS_SMALL, width=0.6, length=3, pad=2)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_linewidth(0.6)

    # (b) neck width
    ax = fig.add_axes(F(ml, y_neck, axw, ph)); dress(ax)
    tt = np.linspace(tk[0], tk[-1], 400)
    ax.plot(t / HOUR, w, "-", lw=1.8, color=C_SIM, zorder=2, label="model")
    ax.plot(tt / HOUR, _pl(tt, f["C"], f["t0"], f["a"]), color=C_FREE, lw=1.0, ls=pmv.FIT_LS, dash_capstyle="round",
            zorder=3, label=rf"$C\,(t+t_0)^{{a}}$, $a={f['a']:.2f}$")
    ax.plot(tt / HOUR, _p3(tt, h["C"], h["t0"]), color=C_THIRD, lw=1.5, ls=pmv.ALT_LS, dash_capstyle="round",
            zorder=3, label=r"$C\,(t+t_0)^{1/3}$")
    for i in (0, 1):
        pmv._circled(ax, t_sec[i] / HOUR, float(np.interp(t_sec[i], t, w)), pmv.LETTERS[i], clip_on=False)
    ax.set_xlim(0, t[-1] / HOUR * 1.03)
    ax.xaxis.set_major_locator(MaxNLocator(8, steps=[1, 2, 2.5, 5, 10]))
    ax.yaxis.set_major_locator(MaxNLocator(5, steps=[1, 2, 2.5, 5, 10]))
    ax.set_xlabel("Time [h]", fontsize=FS, labelpad=2); ax.set_ylabel(r"$w$  [$\mu$m]", fontsize=FS, labelpad=3)
    ax.legend(fontsize=FS_SMALL, frameon=False, handlelength=2.2, handletextpad=0.5, loc="lower right")

    # (c) exponent against the start of the fit window
    cx = fig.add_axes(F(ml2, bot, pw, ph2)); dress(cx)
    cx.axhspan(*DEMMENIE, color=C_THIRD, alpha=0.18, lw=0)
    cx.axhline(1 / 3, color=C_THIRD, lw=1.5, ls=pmv.ALT_LS, dash_capstyle="round")
    cx.plot(s[:, 0] / HOUR, s[:, 2], "-", lw=1.8, color=C_SIM)
    cx.text(0.98, np.mean(DEMMENIE), "Demmenie et al. (2025)", transform=cx.get_yaxis_transform(), ha="right",
            va="center", fontsize=FS_SMALL, color=INK)
    cx.set_xlim(0, None); cx.set_ylim(0.19, 0.35)
    cx.yaxis.set_major_locator(MaxNLocator(5, steps=[1, 2, 2.5, 5, 10]))
    cx.set_xlabel("Start of fit window [h]", fontsize=FS, labelpad=2); cx.set_ylabel(r"$a$", fontsize=FS, labelpad=3)

    # (d) grain diameter: the saturation check
    dx = fig.add_axes(F(ml2 + pw + pgap, bot, pw, ph2)); dress(dx)
    gt, gR = A["g"]
    dx.plot(gt / HOUR, gR / gR[0], "-", lw=1.8, color=C_SIM)
    dx.set_xlim(0, t[-1] / HOUR * 1.03); dx.set_ylim(0.95, 1.01)
    dx.yaxis.set_major_locator(MaxNLocator(5, steps=[1, 2, 2.5, 5, 10]))
    dx.set_xlabel("Time [h]", fontsize=FS, labelpad=2); dx.set_ylabel(r"$D\,/\,D_0$", fontsize=FS, labelpad=3)

    for lab, y_top in zip("ab", (y_snap + s_h + 0.5 * t_band, y_neck + ph + 0.14)):
        fig.text(0.02 / W, y_top / H, pplib.bold(f"({lab})"), ha="left", va="center", fontsize=FS, color=INK)
    for i, lab in enumerate("cd"):
        fig.text((0.02 + i * (pw + pgap + ml2 - 0.72)) / W, (bot + ph2 + 0.16) / H, pplib.bold(f"({lab})"),
                 ha="left", va="center", fontsize=FS, color=INK)

    a.out.mkdir(parents=True, exist_ok=True)
    for e in ("pdf", "png"):
        fn = a.out / f"{a.name}.{e}"
        fig.savefig(fn, dpi=600, transparent=True)
        if a.copy_to:
            a.copy_to.mkdir(parents=True, exist_ok=True); shutil.copyfile(fn, a.copy_to / fn.name)
    print(f"wrote {a.out}/{a.name}.pdf/.png ({W / MM:.0f} x {H / MM:.0f} mm); free a = {f['a']:.3f}, one-third rms {h['rms']:.2f} %")


if __name__ == "__main__":
    main()
