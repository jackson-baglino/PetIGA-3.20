#!/usr/bin/env python3
"""Figure 1, panel (b) candidate: sintering in an aggregate, close up.

    venv_enceladus/bin/python studies/keff_sintering/figures/fig1_aggregate_strip.py <campaign dir>
        [--center-um 700 1450] [--window-um 420] [--copy-to <dir>]

A close-up of a few dozen grains of the MASTER simulation (phi 0.325, seed
1702, -20 C) at four instants: necks form, small pores round off and close.
It carries the two-grain schematic of Fig. 1 over to a packing, without
repeating the full-cell snapshots of the k_eff evolution figure. Ice only (no
vapour field), so it reads as a picture, not as data. Width 170 mm.
Writes aggregate_strip.{pdf,png}.
"""
from __future__ import annotations

import argparse, shutil, sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cmocean

HERE = Path(__file__).resolve().parent
PROJ = HERE.parents[2]
sys.path.insert(0, str(PROJ / "postprocess"))
import pplib  # noqa: E402
from pplib import step_times, opening_step  # noqa: E402
from plot_keff_snapshots import make_reader, snap_step, INK, MUTED, FS, FS_SMALL, FS_TINY  # noqa: E402

MM, DAY = 1 / 25.4, 86400.0
MASTER = "packing_2D_phi0.325_Rave50um_LR40_seed1702_L2mm_eps1000nm_perxy_T-20__snow_T-20_h1.00_30d"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root", type=Path)
    ap.add_argument("--center-um", type=float, nargs=2, default=[700.0, 1450.0])
    ap.add_argument("--window-um", type=float, default=420.0)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--copy-to", type=Path, default=None)
    a = ap.parse_args()
    plt.rcParams.update(pplib.MANUSCRIPT_RC)
    run = a.root / MASTER
    files, reader = make_reader(run, "sol")
    tm = step_times(str(run)); st = [snap_step(f) for f in files]
    tt = np.array([tm.get(s, np.nan) for s in st])
    op = opening_step(st, list(tt)); i0 = st.index(op) if op is not None else 0
    te = max(tm.values())
    pick = [i0] + [int(np.nanargmin(np.abs(tt - (tt[i0] + f * (te - tt[i0]))))) for f in (1 / 3, 2 / 3)] + [len(st) - 1]
    cx, cy = a.center_um; h = a.window_um / 2
    W, gap, L, top = 170.0, 2.0, 1.0, 6.0
    cell = (W - 2 * L - 3 * gap) / 4
    H = cell + top + 1.0
    fig = plt.figure(figsize=(W * MM, H * MM))
    for j, i in enumerate(pick):
        fl, X, Y = reader(files[i], want=("IcePhase",))
        x, y = X[0, :] * 1e6, Y[:, 0] * 1e6
        mx = (x >= cx - h) & (x <= cx + h); my = (y >= cy - h) & (y <= cy + h)
        ax = fig.add_axes([(L + j * (cell + gap)) / W, 1.0 / H, cell / W, cell / H])
        ax.imshow(fl["IcePhase"][np.ix_(my, mx)], origin="lower", cmap=cmocean.cm.ice, vmin=0, vmax=1,
                  extent=(x[mx][0], x[mx][-1], y[my][0], y[my][-1]), interpolation="antialiased")
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_linewidth(0.6); sp.set_color(MUTED)
        ax.set_title("0 d" if j == 0 else f"{tt[i] / DAY:.0f} d", fontsize=FS, pad=3)
        if j == 0:
            ax.plot([cx - h + 20, cx - h + 120], [cy - h + 22, cy - h + 22], color=INK, lw=2.2,
                    solid_capstyle="butt")
            ax.text(cx - h + 70, cy - h + 32, r"100 $\mu$m", ha="center", va="bottom", fontsize=FS_TINY,
                    color=INK, bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.85))
    out = a.out or a.root / "compare" / "figure_samples"
    out.mkdir(parents=True, exist_ok=True)
    for e in ("pdf", "png"):
        f = out / f"aggregate_strip.{e}"
        fig.savefig(f, dpi=500, transparent=True)
        if a.copy_to:
            a.copy_to.mkdir(parents=True, exist_ok=True); shutil.copyfile(f, a.copy_to / f.name)
    print(f"wrote {out}/aggregate_strip.pdf/.png ({W:.0f} x {H:.0f} mm); times [d]: {[round(tt[i] / DAY, 1) for i in pick]}")


if __name__ == "__main__":
    main()
