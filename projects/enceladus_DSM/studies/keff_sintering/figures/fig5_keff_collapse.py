#!/usr/bin/env python3
"""Figure 5 (manuscript): k_eff evolution of one packing, and its temperature collapse.

    venv_enceladus/bin/python studies/keff_sintering/figures/fig5_keff_collapse.py <campaign dir>
        [--seed 1702] [--phi 0.325] [--snap-T -20] [--copy-to <dir>] [--out <dir>]

Merges the snapshot figure (plot_keff_snapshots.py) with the temperature
collapse, on ONE packing run at five temperatures:

  (a) four microstructure snapshots of the --snap-T run (opening frame,
      ~t_final/3, ~2 t_final/3, last), numbered 1-4
  (b) k_eff / k_eff,0 against time [d], one curve per temperature; the numbered
      instants sit on the --snap-T curve
  (c) the same curves against sintering age theta = t / tau_sub: one curve
  (d) SSA / SSA_0 against theta: the microstructure itself collapses

One packing, not a seed mean, so the snapshots ARE the curve. k_eff is k_iso.
Manuscript style: 170 mm wide, no titles, symbol labels, transparent
background, >= 8 pt. Writes Figure5_keff_collapse.{pdf,png}; --copy-to also
copies them under that name (never over the Inkscape assemblies there).
"""
from __future__ import annotations

import argparse, re, shutil, sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
import cmocean

HERE = Path(__file__).resolve().parent
PROJ = HERE.parents[2]
sys.path.insert(0, str(PROJ / "postprocess"))
import pplib  # noqa: E402
from pplib import step_times, opening_step  # noqa: E402
from plot_keff import load, read_tau_sub  # noqa: E402
from compare_keff import CMAP  # noqa: E402
from plot_keff_snapshots import (make_reader, snap_step, _field, _scalebar, _colorbars,  # noqa: E402
                                 _mark, _snap_title, WANT, SIGMA_SCALE, centered_cmap,
                                 ice_alpha_cmap, INK, FS, FS_SMALL, FS_TINY)

MM = 1 / 25.4
DAY = 86400.0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root", type=Path)
    ap.add_argument("--phi", default="0.325")
    ap.add_argument("--seed", type=int, default=1702)
    ap.add_argument("--snap-T", dest="snap_T", type=int, default=-20)
    ap.add_argument("--width-mm", type=float, default=170.0)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--copy-to", type=Path, default=None)
    a = ap.parse_args()
    plt.rcParams.update(pplib.MANUSCRIPT_RC)

    runs = {}
    for d in sorted(a.root.glob(f"packing_2D_phi{a.phi}_Rave50um_LR40_seed{a.seed}_L2mm_eps1000nm_perxy_T*__*")):
        T = int(re.search(r"_T(-?\d+)__", d.name).group(1))
        r = load(d)
        if r is None:
            continue
        i0 = int(np.argmax(r["t"] >= 1.0))
        r.update(dir=d, tau=read_tau_sub(d), i0=i0)
        for k in ("t", "kiso", "ssa", "step"):
            r[k] = r[k][i0:]
        r["kn"] = r["kiso"] / r["kiso"][0]; r["sn"] = r["ssa"] / r["ssa"][0]
        r["th"] = r["t"] / r["tau"]
        runs[T] = r
    Ts = sorted(runs)
    cm, (c0, c1) = CMAP["T"]
    col = {T: cm(c0 + (c1 - c0) * i / max(1, len(Ts) - 1)) for i, T in enumerate(Ts)}

    # ---- snapshots of the --snap-T run ----
    S = runs[a.snap_T]
    files, reader = make_reader(S["dir"], "sol")
    tm = step_times(str(S["dir"]))
    st = [snap_step(f) for f in files]; tt = [tm.get(x, np.nan) for x in st]
    op = opening_step(st, tt)
    i_open = st.index(op) if op is not None else 0
    t_end = max(tm.values())
    pick = [i_open] + [int(np.argmin(np.abs(np.array(tt) - (tt[i_open] + f * (t_end - tt[i_open])))))
                       for f in (1 / 3, 2 / 3)] + [len(st) - 1]
    snaps, pore = [], []
    for i in pick:
        fl, X, Y = reader(files[i], want=WANT)
        snaps.append((fl, X, Y, tt[i]))
        s = SIGMA_SCALE * pplib.supersaturation(fl["VaporDensity"], fl["Temperature"])
        pore.append(s[fl["IcePhase"] < 0.5])
    pore = np.concatenate(pore); smin, smax = float(pore.min()), float(pore.max())
    v = min(abs(smin), abs(smax)) or max(abs(smin), abs(smax))
    norm = AsinhNorm(linear_width=max(v / 300.0, 1e-12), vmin=-v, vmax=v)
    ext = {(True, True): "both", (True, False): "min", (False, True): "max",
           (False, False): "neither"}[(smin < -v, smax > v)]
    vapcm, icecm = centered_cmap(cmocean.cm.balance, norm), ice_alpha_cmap()

    # ---- layout (mm) ----
    W = a.width_mm
    L, Rm, gap = 13.0, 2.5, 2.5
    snap_w = (W - L - Rm - 3 * gap) / 4
    cb_h, cb_top, title_h = 2.2, 4.0, 6.0
    cur_h, cur_gap, bot = 44.0, 14.0, 12.0
    H = cb_top + cb_h + 7.0 + title_h + snap_w + cur_gap + cur_h + bot
    fig = plt.figure(figsize=(W * MM, H * MM))
    fx = lambda x: x / W
    fy = lambda y: y / H
    y_snap = bot + cur_h + cur_gap
    y_cb = y_snap + snap_w + title_h + 7.0
    cax_i = fig.add_axes([fx(L + 8), fy(y_cb), fx(38), fy(cb_h)])
    cax_s = fig.add_axes([fx(L + 80), fy(y_cb), fx(W - L - 80 - Rm - 4), fy(cb_h)])
    _colorbars(fig, cax_i, cax_s, norm, vapcm, ext)
    fig.canvas.draw()
    for j, (fl, X, Y, t) in enumerate(snaps):
        ax = fig.add_axes([fx(L + j * (snap_w + gap)), fy(y_snap), fx(snap_w), fy(snap_w)])
        XX, YY = _field(ax, fl, X, Y, norm, vapcm, icecm)
        if j == 0:
            _scalebar(ax, XX, YY)
        _snap_title(ax, str(j + 1), f"{t / DAY:.1f} d" if j else "0 d")
    fig.text(fx(1.5), fy(y_snap + snap_w + title_h - 1), pplib.bold("(a)"), fontsize=FS, va="center")

    cw = (W - L - Rm - 2 * 15.0) / 3
    axs = [fig.add_axes([fx(L + j * (cw + 15.0)), fy(bot), fx(cw), fy(cur_h)]) for j in range(3)]
    for T in Ts:
        r = runs[T]
        lab = rf"${T}\,^\circ$C"
        axs[0].plot(r["t"] / DAY, r["kn"], color=col[T], lw=1.5, label=lab)
        axs[1].plot(r["th"], r["kn"], color=col[T], lw=1.5)
        axs[2].plot(r["th"], r["sn"], color=col[T], lw=1.5)
    for j, (_, _, _, t) in enumerate(snaps):
        _mark(axs[0], t / DAY, float(np.interp(t, S["t"], S["kn"])), str(j + 1))
    axs[0].set(xlabel="$t$ [d]", ylabel=r"$k_\mathrm{eff}\,/\,k_{\mathrm{eff},0}$")
    axs[1].set(xlabel=r"$\theta=t/\tau_\mathrm{sub}$", ylabel=r"$k_\mathrm{eff}\,/\,k_{\mathrm{eff},0}$", xscale="log")
    axs[2].set(xlabel=r"$\theta=t/\tau_\mathrm{sub}$", ylabel=r"SSA$\,/\,$SSA$_0$", xscale="log")
    for x in axs[1:]:
        x.set_xlim(1, None)
    axs[0].legend(frameon=False, fontsize=FS_TINY, loc="lower right", handlelength=1.3,
                  labelspacing=0.2, borderaxespad=0.2)
    for x, s in zip(axs, "bcd"):
        for sp in ("top", "right"):
            x.spines[sp].set_visible(False)
        x.tick_params(labelsize=FS_SMALL, width=0.6, length=3)
        x.xaxis.label.set_size(FS); x.yaxis.label.set_size(FS)
        x.text(-0.30, 1.04, pplib.bold(f"({s})"), transform=x.transAxes, fontsize=FS, va="bottom")

    out = a.out or a.root / "compare" / "figure_samples"
    out.mkdir(parents=True, exist_ok=True)
    for e in ("pdf", "png"):
        f = out / f"Figure5_keff_collapse.{e}"
        fig.savefig(f, dpi=600, transparent=True)
        if a.copy_to:
            a.copy_to.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(f, a.copy_to / f.name)
    print(f"wrote {out}/Figure5_keff_collapse.pdf/.png ({W:.0f} x {H:.0f} mm)"
          + (f"; copied to {a.copy_to}" if a.copy_to else ""))
    print("temperatures:", Ts, "| snapshot times [d]:", [round(s[3] / DAY, 2) for s in snaps])


if __name__ == "__main__":
    main()
