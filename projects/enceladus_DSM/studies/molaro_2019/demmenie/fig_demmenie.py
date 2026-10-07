#!/usr/bin/env python3
"""fig_demmenie.py — manuscript figure: neck growth at saturation against t^(1/3).

    venv_enceladus/bin/python studies/molaro_2019/demmenie/fig_demmenie.py <run dir>
        [--relax-tau 11] [--name Figure4_saturated_neck_growth] [--out <dir>] [--copy-to <dir>]

  (a) the grain pair at the start and at the end of the run: the computed
      quarter (one grain, half-plane) reflected across the mirror plane and
      the symmetry axis. Ice only.
  (b) neck width against time: the simulation (line), the free fit
      C (t + t0)^a and the one-third fit C (t + t0)^(1/3), both over the
      samples after the relaxation period (shaded).
  (c) the free-fit exponent against the start of the fit window, with the
      range Demmenie et al. (2025) measured and the 1/3 law.

Fits and numbers come from analyze_demmenie.py, so the figure and the
diagnostic plots cannot disagree. Manuscript style: 170 mm, no titles, symbol
labels, transparent background (pplib.MANUSCRIPT_RC).
"""
from __future__ import annotations

import argparse, glob, shutil, sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cmocean

HERE = Path(__file__).resolve().parent
PROJ = HERE.parents[2]
sys.path.insert(0, str(PROJ / "postprocess")); sys.path.insert(0, str(HERE))
import pplib  # noqa: E402
from plot_keff_snapshots import make_reader, snap_step, INK, MUTED, FS, FS_SMALL, FS_TINY  # noqa: E402
from analyze_demmenie import analyse, _pl, _p3, DEMMENIE, HOUR  # noqa: E402

MM = 1 / 25.4
C_FREE, C_THIRD = "#2a78d6", "#d1495b"


def section(run, fn, stride=3):
    """Ice field of the full pair: reflect across r = 0 and across z = 0."""
    _, reader = make_reader(run, "sol")
    fl, X, Y = reader(fn, want=("IcePhase",))
    p, z, r = fl["IcePhase"][::stride, ::stride], X[0, ::stride] * 1e6, Y[::stride, 0] * 1e6
    p = np.vstack([p[:0:-1], p]); p = np.hstack([p[:, :0:-1], p])
    return p, np.concatenate([-z[:0:-1], z]), np.concatenate([-r[:0:-1], r])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run", type=Path)
    ap.add_argument("--relax-tau", type=float, default=11.0)
    ap.add_argument("--name", default="Figure4_saturated_neck_growth")
    ap.add_argument("--out", type=Path, default=PROJ / "studies/molaro_2019/manuscript")
    ap.add_argument("--copy-to", type=Path, default=None)
    a = ap.parse_args()
    plt.rcParams.update(pplib.MANUSCRIPT_RC)
    A = analyse(a.run, a.relax_tau)
    t, w, tk, f, h, s = A["t"], A["w"], A["tk"], A["free"], A["third"], A["scan"]

    files = sorted(glob.glob(str(a.run / "sol_*.dat")), key=snap_step)
    secs = [section(a.run, files[0]), section(a.run, files[-1])]
    zc = 1.12 * max(np.abs(z[np.any(p >= 0.5, axis=0)]).max() for p, z, r in secs)
    rc = 1.18 * max(np.abs(r[np.any(p >= 0.5, axis=1)]).max() for p, z, r in secs)

    W = 170.0
    L, Rm, gap, top, bot = 13.0, 2.0, 4.0, 7.0, 12.0
    sw = (W - L - Rm - gap) / 2; sh = sw * rc / zc
    ph, row_gap = 52.0, 13.0
    H = top + sh + row_gap + ph + bot
    fig = plt.figure(figsize=(W * MM, H * MM))
    F = lambda x, y, ww, hh: [x / W, y / H, ww / W, hh / H]

    for j, ((p, z, r), lab) in enumerate(zip(secs, ("0 h", f"{t[-1] / HOUR:.0f} h"))):
        ax = fig.add_axes(F(L + j * (sw + gap), bot + ph + row_gap, sw, sh))
        ax.imshow(p, origin="lower", cmap=cmocean.cm.ice, vmin=0, vmax=1, extent=(z[0], z[-1], r[0], r[-1]),
                  interpolation="antialiased", aspect="auto")
        ax.set(xlim=(-zc, zc), ylim=(-rc, rc)); ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_linewidth(0.6); sp.set_color(MUTED)
        ax.set_title(lab, fontsize=FS, pad=3)
        if j == 0:
            x0, y0 = -zc + 0.02 * 2 * zc, rc - 0.115 * 2 * rc   # the dark corner above the grain
            ax.plot([x0, x0 + 50], [y0, y0], color="white", lw=2.2, solid_capstyle="butt")
            ax.text(x0 + 25, y0 + 0.03 * 2 * rc, r"50 $\mu$m", color="white", ha="center", va="bottom", fontsize=FS_TINY)
    fig.text(1.5 / W, (bot + ph + row_gap + sh + 3.0) / H, pplib.bold("(a)"), fontsize=FS, va="center", color=INK)

    pgap = 17.0
    pw = (W - L - Rm - pgap) / 2
    ax = fig.add_axes(F(L, bot, pw, ph))
    tt = np.linspace(tk[0], tk[-1], 400)
    ax.axvspan(0, A["t_relax"] / HOUR, color="0.88", lw=0)
    ax.plot(t / HOUR, w, color=INK, lw=2.0, label="simulation", zorder=3)
    ax.plot(tt / HOUR, _pl(tt, f["C"], f["t0"], f["a"]), color=C_FREE, lw=1.3, ls=(0, (4, 2)),
            label=rf"$C\,(t+t_0)^{{a}}$, $a={f['a']:.2f}$", zorder=4)
    ax.plot(tt / HOUR, _p3(tt, h["C"], h["t0"]), color=C_THIRD, lw=1.3, ls=(0, (1.2, 1.6)),
            label=r"$C\,(t+t_0)^{1/3}$", zorder=4)
    ax.set(xlabel="$t$ [h]", ylabel=r"$w$ [$\mu$m]", xlim=(0, None))
    ax.legend(frameon=False, fontsize=FS_SMALL, loc="lower right", handlelength=2.2, labelspacing=0.3)
    bx = fig.add_axes(F(L + pw + pgap, bot, pw, ph))
    bx.axhspan(*DEMMENIE, color=C_THIRD, alpha=0.16, lw=0)
    bx.axhline(1 / 3, color=C_THIRD, lw=1.0, ls=(0, (1.2, 1.6)))
    bx.plot(s[:, 0] / HOUR, s[:, 2], color=C_FREE, lw=2.0)
    bx.text(0.97, 1 / 3 + 0.004, "1/3", transform=bx.get_yaxis_transform(), ha="right", va="bottom",
            fontsize=FS_SMALL, color=INK)
    bx.text(0.97, np.mean(DEMMENIE), "Demmenie et al. (2025)", transform=bx.get_yaxis_transform(), ha="right",
            va="center", fontsize=FS_SMALL, color=INK)
    bx.set(xlabel="start of fit window [h]", ylabel="$a$", ylim=(0.18, 0.36), xlim=(0, None))
    for x, lab in ((ax, "b"), (bx, "c")):
        for sp in ("top", "right"):
            x.spines[sp].set_visible(False)
        x.tick_params(labelsize=FS_SMALL, width=0.6, length=3)
        x.xaxis.label.set_size(FS); x.yaxis.label.set_size(FS); x.patch.set_alpha(0)
        x.text(-0.17, 1.03, pplib.bold(f"({lab})"), transform=x.transAxes, fontsize=FS, va="bottom", color=INK)

    a.out.mkdir(parents=True, exist_ok=True)
    for e in ("pdf", "png"):
        fn = a.out / f"{a.name}.{e}"
        fig.savefig(fn, dpi=600, transparent=True)
        if a.copy_to:
            a.copy_to.mkdir(parents=True, exist_ok=True); shutil.copyfile(fn, a.copy_to / fn.name)
    print(f"wrote {a.out}/{a.name}.pdf/.png ({W:.0f} x {H:.0f} mm); free a = {f['a']:.3f}, one-third rms {h['rms']:.2f} %")


if __name__ == "__main__":
    main()
