#!/usr/bin/env python3
"""fig_demmenie_slide.py -- slide figure: neck growth at saturation, fit, theory and experiment.

    venv_enceladus/bin/python studies/molaro_2019/demmenie/fig_demmenie_slide.py <run dir>
        [--relax-tau 11] [--out <dir>] [--name fig4_saturated_slide]

The manuscript's Figure 4(b) draws the free fit C (t+t0)^a and the one-third fit
over the simulation, thin and dashed, on a linear axis. The free fit lies within
0.2 % rms of the simulation and the one-third law within 2.6 % rms, so both hide
under the simulation line and the figure looks as if it had no fit or theory.
This slide version makes each one visible:

  (a) neck width against time. The simulation is drawn as open circles (thinned,
      so the curves show between them); the free fit (a = fitted) is a solid
      line and the one-third law a dashed line, both over the samples after the
      relaxation period (the grey band).
  (b) the free-fit exponent against the start of the fit window, with its 2-sigma
      band, the measured exponents of Demmenie et al. (2025) as a shaded band,
      and the one-third law -- where the experiment sits relative to the model
      and to the theory.

Every number comes from analyze_demmenie.analyse, as in fig_demmenie.py, so the
slide cannot disagree with the manuscript figure or the diagnostic plots.
"""
from __future__ import annotations

import argparse, sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from analyze_demmenie import analyse, _pl, _p3, DEMMENIE, HOUR  # noqa: E402

INK, MUTED = "#1E2930", "#6B7280"
C_SIM, C_FIT, C_THEORY, C_EXP = "#1E2930", "#1E88E5", "#E8532B", "#E8532B"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run", type=Path)
    ap.add_argument("--relax-tau", type=float, default=11.0)
    ap.add_argument("--out", type=Path, default=HERE / "manuscript")
    ap.add_argument("--name", default="fig4_saturated_slide")
    ap.add_argument("--width-in", type=float, default=11.6)
    ap.add_argument("--height-in", type=float, default=4.3)
    a = ap.parse_args()
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                         "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic",
                         "mathtext.bf": "Arial:bold", "pdf.fonttype": 42})
    A = analyse(a.run, a.relax_tau)
    t, w, tk, f, h, s = A["t"], A["w"], A["tk"], A["free"], A["third"], A["scan"]
    tt = np.linspace(tk[0], tk[-1], 400)

    FS, FS_S = 16, 14
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(a.width_in, a.height_in), gridspec_kw={"width_ratios": [1.25, 1]})
    fig.subplots_adjust(left=0.10, right=0.985, bottom=0.17, top=0.96, wspace=0.27)

    # (a) neck width
    ax.axvspan(0, A["t_relax"] / HOUR, color="0.92", lw=0)
    thin = np.unique(np.linspace(0, len(t) - 1, 38).astype(int))
    ax.plot(t[thin] / HOUR, w[thin], "o", ms=6.5, mfc="white", mec=C_SIM, mew=1.4, ls="none", zorder=2,
            label="simulation")
    ax.plot(tt / HOUR, _pl(tt, f["C"], f["t0"], f["a"]), color=C_FIT, lw=3.2, zorder=3,
            label=rf"fit  $C\,(t+t_0)^{{a}}$,  $a={f['a']:.2f}$")
    ax.plot(tt / HOUR, _p3(tt, h["C"], h["t0"]), color=C_THEORY, lw=2.6, ls=(0, (5, 3)), zorder=4,
            label=r"theory  $C\,(t+t_0)^{1/3}$")
    ax.set_xlim(0, t[-1] / HOUR * 1.03); ax.set_ylim(0, None)
    ax.xaxis.set_major_locator(MaxNLocator(6, steps=[1, 2, 5, 10]))
    ax.set_xlabel("time  [h]", fontsize=FS); ax.set_ylabel(r"neck width  $w$  [$\mu$m]", fontsize=FS)
    ax.legend(fontsize=FS_S, frameon=False, loc="lower right", handlelength=2.4)

    # (b) exponent against fit window
    x0 = s[:, 0] / HOUR
    bx.axhspan(*DEMMENIE, color=C_EXP, alpha=0.16, lw=0)
    bx.text(x0[-1], DEMMENIE[1] - 0.004, "Demmenie et al. (2025)\nmeasured, four runs", ha="right", va="top",
            fontsize=FS_S - 1, color=C_THEORY, linespacing=1.15)
    bx.axhline(1 / 3, color=C_THEORY, lw=2.2, ls=(0, (5, 3)))
    bx.text(x0[0], 1 / 3 + 0.004, "theory  1/3", ha="left", va="bottom", fontsize=FS_S - 1, color=C_THEORY)
    bx.fill_between(x0, s[:, 2] - 2 * s[:, 3], s[:, 2] + 2 * s[:, 3], color=C_FIT, alpha=0.2, lw=0)
    bx.plot(x0, s[:, 2], color=C_FIT, lw=3.2, label="simulation: fitted exponent")
    bx.set_ylim(0.15, 0.37)
    bx.set_xlabel("start of the fit window  [h]", fontsize=FS); bx.set_ylabel(r"exponent  $a$", fontsize=FS)
    bx.legend(fontsize=FS_S, frameon=False, loc="lower left", handlelength=2.0)

    for axx, lab in ((ax, "(a)"), (bx, "(b)")):
        axx.patch.set_alpha(0.0)
        axx.tick_params(labelsize=FS_S, width=1.0, length=4)
        for sp in ("top", "right"):
            axx.spines[sp].set_visible(False)
        axx.text(-0.17 if axx is ax else -0.22, 1.02, lab, transform=axx.transAxes, fontsize=FS + 1,
                 fontweight="bold", color=INK, va="top")
    a.out.mkdir(parents=True, exist_ok=True)
    for e in ("png", "pdf"):
        fig.savefig(a.out / f"{a.name}.{e}", dpi=300, transparent=True)
    print(f"  fit a = {f['a']:.3f} +- {f['a_se']:.3f} (rms {f['rms']:.2f} %); one-third rms {h['rms']:.2f} %; "
          f"Demmenie {DEMMENIE}; wrote {a.out}/{a.name}.png/.pdf")


if __name__ == "__main__":
    main()
