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
The exponent's dependence on the fit window and the grain-size check were
panels (c) and (d) until 2026-10-07; they are diagnostic plots now
(analyze_demmenie.py: exponent_window.png, grain_radius.png).

--v4 (2026-10-16) writes <name>_v4: v2's labels with the sigma bar over the
data's own range, as the Molaro figure's full-range version. The pore here is
supersaturated everywhere, so the bar runs over the pore's min .. max on the
red half of the map. Half the pore sits within 0.1 % of the max (the far
field is genuinely flat), so --sig-gamma 3 (the manuscript choice) stretches
the top of the range and the depletion halo round the neck shows; pale is
the pore's minimum, NOT saturation as in the Molaro figure. --sig-amp (the
user's choice, 2026-10-16): cmocean amp from sigma = 0 (white) to the max,
so saturation is always on the bar; with --sig-gamma 3 the halo still shows. _v4_lin is the
linear alternative (--sig-gamma 1 --suffix _lin).
--v2 (advisor meeting 2026-10-15) writes <name>_v2: "simulation" for "model",
the one-third curve labelled as the theory, and an I-bar marking w(t) on each
section whose neck is open (none on the tangent-start 0 h one, where w = 0).

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
    ap.add_argument("--v2", action="store_true", help="the 2026-10-15 labels (see above)")
    ap.add_argument("--v4", action="store_true", help="v2 + sigma bar over 0 .. max (see above)")
    ap.add_argument("--sig-gamma", type=float, default=1.0,
                    help="with --v4: PowerNorm exponent for the sigma bar (1 = linear)")
    ap.add_argument("--sig-amp", action="store_true",
                    help="with --v4: cmocean amp, bar from sigma = 0 (white) to the max")
    ap.add_argument("--suffix", default="")
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
    # Bounds as in the Molaro figure: symmetric about zero, so sigma = 0 is the
    # bar's middle, and sized to THIS run's data. There the size is the smaller
    # of the two extremes; here the pore is supersaturated everywhere against
    # a flat surface (the walls sit at 1 + 2 d0/R), so there is one extreme and
    # the bar ends exactly on it. The mapping is linear, not asinh: the field
    # spans a factor of four, not decades, and asinh would paint the whole
    # pore the end colour and hide the depletion at the neck.
    if a.v4:
        a.v2 = True
    if a.v4 and smin >= 0:
        from matplotlib.colors import PowerNorm
        lo = 0.0 if a.sig_amp else smin          # amp: sigma = 0 always on the bar, white
        norm = (PowerNorm(gamma=a.sig_gamma, vmin=lo, vmax=smax) if a.sig_gamma != 1.0
                else Normalize(vmin=lo, vmax=smax)); ext = "neither"
        vapcm = cmocean.cm.amp if a.sig_amp else pmv.centered_cmap(cmocean.cm.balance, norm)
        v = smax
    elif smin * smax > 0:
        v = max(abs(smin), abs(smax))
        norm = Normalize(vmin=-v, vmax=v); ext = "neither"; vapcm = cmocean.cm.balance
    else:
        v = min(abs(smin), abs(smax)) or max(abs(smin), abs(smax))
        ext = {(True, True): "both", (True, False): "min", (False, True): "max", (False, False): "neither"}[(smin < -v, smax > v)]
        norm = AsinhNorm(linear_width=max(v / 300.0, 1e-12), vmin=-v, vmax=v)
        vapcm = pmv.centered_cmap(cmocean.cm.balance, norm)
    icecm = pmv.ice_alpha_cmap()
    print(f"sigma x{pmv.SIGMA_SCALE:g} in the shown pore: {smin:+.3g} .. {smax:+.3g}; bar +-{v:.3g} ({ext}); "
          f"quantiles 1/10/50/90 %: " + " ".join(f"{q:.4g}" for q in np.percentile(pore, [1, 10, 50, 90])))

    # ---- layout, in inches, as plot_molaro_validation.build_full ----
    W = a.width_mm * MM
    gap = 0.12
    ml, mr = 0.50, 0.08
    axw = W - ml - mr
    s_w = (axw - gap) / 2
    fl0, X0, Y0 = secs[0]
    s_h = s_w * (Y0.max() - Y0.min()) / (X0.max() - X0.min())
    cb_h, cb_lab, cb_gap, t_band, top = 0.07, 0.15, 0.06, 0.19, 0.05
    g_snap, ph, bot = 0.34, 2.30, 0.40          # as plot_molaro_validation.build_combined
    H = top + cb_h + cb_lab + cb_gap + t_band + s_h + g_snap + ph + bot
    fig = plt.figure(figsize=(W, H))
    F = lambda x0, y0, ww, hh: (x0 / W, y0 / H, ww / W, hh / H)
    y_neck = bot
    y_snap = y_neck + ph + g_snap

    t_sec = [t[0], t[-1]]
    for i, (fl, X, Y) in enumerate(secs):
        axi = fig.add_axes(F(ml + i * (s_w + gap), y_snap, s_w, s_h))
        XX, YY = _field(axi, fl, X, Y, norm, vapcm, icecm)
        axi.set_aspect("auto")
        if i == 0:
            _scalebar(axi, XX, YY)
        _snap_title(axi, pmv.LETTERS[i], f"{t_sec[i] / HOUR:.0f} h")
        if a.v2:
            z_mm, w_mm = pmv.neck_at(a.run, t_sec[i])
            if w_mm > 0.02 * (Y0.max() - Y0.min()) * 1e3:      # a visible neck only
                pmv.neck_bar(axi, z_mm, w_mm)
                print(f"w(t) bar on instant {i + 1}: w = {1e3 * w_mm:.2f} um")
    a.sig_extend = ext
    pmv._strip(fig, F, ml, y_snap + s_h + t_band + cb_gap + cb_lab, axw, norm, vapcm, ext)
    if a.v4 and a.sig_amp:
        # sigma_ticks thins by position on the bar, and the gamma stretch
        # crowds 0 .. 0.1 into its first few percent, leaving only the ends.
        cax_sig = fig.axes[-1]                       # _strip adds the sigma bar last
        sig_ticks = [x for x in (0.0, 0.1, 0.15, 0.2) if x < 0.97 * smax] + [smax]   # not tk: that is the fit window
        cax_sig.set_xticks(sig_ticks)
        cax_sig.set_xticklabels([f"${x:.2g}$" for x in sig_ticks])

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
    ax.plot(t / HOUR, w, "-", lw=1.8, color=C_SIM, zorder=2, label="simulation" if a.v2 else "model")
    ax.plot(tt / HOUR, _pl(tt, f["C"], f["t0"], f["a"]), color=C_FREE, lw=1.0, ls=pmv.FIT_LS, dash_capstyle="round",
            zorder=3, label=rf"$C\,(t+t_0)^{{a}}$, $a={f['a']:.2f}$")
    ax.plot(tt / HOUR, _p3(tt, h["C"], h["t0"]), color=C_THIRD, lw=1.5, ls=pmv.ALT_LS, dash_capstyle="round",
            zorder=3, label=(r"theory, $C\,(t+t_0)^{1/3}$" if a.v2 else r"$C\,(t+t_0)^{1/3}$"))
    for i in (0, 1):
        pmv._circled(ax, t_sec[i] / HOUR, float(np.interp(t_sec[i], t, w)), pmv.LETTERS[i], clip_on=False)
    ax.set_xlim(0, t[-1] / HOUR * 1.03)
    ax.xaxis.set_major_locator(MaxNLocator(8, steps=[1, 2, 2.5, 5, 10]))
    ax.yaxis.set_major_locator(MaxNLocator(5, steps=[1, 2, 2.5, 5, 10]))
    ax.set_xlabel("Time [h]", fontsize=FS, labelpad=2); ax.set_ylabel(r"$w$  [$\mu$m]", fontsize=FS, labelpad=3)
    ax.legend(fontsize=FS_SMALL, frameon=False, handlelength=2.2, handletextpad=0.5, loc="lower right")

    for lab, y_top in zip("ab", (y_snap + s_h + 0.5 * t_band, y_neck + ph + 0.14)):
        fig.text(0.02 / W, y_top / H, pplib.bold(f"({lab})"), ha="left", va="center", fontsize=FS, color=INK)
    a.out.mkdir(parents=True, exist_ok=True)
    for e in ("pdf", "png"):
        fn = a.out / f"{a.name}{'_v4' if a.v4 else '_v2' if a.v2 else ''}{a.suffix}.{e}"
        fig.savefig(fn, dpi=600, transparent=True)
        if a.copy_to:
            a.copy_to.mkdir(parents=True, exist_ok=True); shutil.copyfile(fn, a.copy_to / fn.name)
    print(f"wrote {a.out}/{a.name}{'_v4' if a.v4 else '_v2' if a.v2 else ''}{a.suffix}.pdf/.png ({W / MM:.0f} x {H / MM:.0f} mm); free a = {f['a']:.3f}, one-third rms {h['rms']:.2f} %")


if __name__ == "__main__":
    main()
