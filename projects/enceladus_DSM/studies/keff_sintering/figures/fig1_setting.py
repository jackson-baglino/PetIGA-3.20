#!/usr/bin/env python3
"""Figure 1, "setting" panel candidate: what the periodic cell stands for.

    venv_enceladus/bin/python studies/keff_sintering/figures/fig1_setting.py <campaign dir>
        [--copy-to <dir>]

Left: a schematic column of the near subsurface -- plume fallout on top, the
deposit below, a sample about a metre down, heat conducted up through it.
Right: that sample as the model sees it -- the MASTER packing (phi 0.325,
seed 1702) at its opening frame, tiled 3 x 3 to show the periodicity, with
the computed 2 mm cell outlined. The left drawing is a cartoon (random discs,
not data, not to scale); the right one is the simulation's initial condition.
No panel letter: it is meant to be placed in the Fig. 1 assembly.
Writes setting_panel.{pdf,png}. Width 170 mm.
"""
from __future__ import annotations

import argparse, shutil, sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, ConnectionPatch
import cmocean

HERE = Path(__file__).resolve().parent
PROJ = HERE.parents[2]
sys.path.insert(0, str(PROJ / "postprocess"))
import pplib  # noqa: E402
from pplib import step_times, opening_step  # noqa: E402
from plot_keff_snapshots import make_reader, snap_step, INK, MUTED, FS, FS_SMALL  # noqa: E402
from fig1_aggregate_strip import MASTER  # noqa: E402

MM = 1 / 25.4
ICE, PORE = cmocean.cm.ice(1.0), cmocean.cm.ice(0.0)
BOX = dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.9)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root", type=Path)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--copy-to", type=Path, default=None)
    a = ap.parse_args()
    plt.rcParams.update(pplib.MANUSCRIPT_RC)

    W, H = 170.0, 70.0
    fig = plt.figure(figsize=(W * MM, H * MM))

    # ---- left: the column (mm units inside the axes) ----
    cw, ch, surf = 78.0, 64.0, 46.0
    ax = fig.add_axes([12 / W, 3 / H, cw / W, ch / H])
    ax.set(xlim=(0, cw), ylim=(0, ch), aspect="equal"); ax.axis("off")
    clip = Rectangle((0, 0), cw, surf + 2, transform=ax.transData)
    ax.add_patch(Rectangle((0, 0), cw, surf - 1.2, fc=PORE, ec="none"))
    rng = np.random.default_rng(7)
    for _ in range(760):
        r = float(np.clip(rng.lognormal(np.log(0.95), 0.45), 0.35, 2.2))
        c = Circle((rng.uniform(0, cw), rng.uniform(-1, surf - r)), r, fc=ICE, ec="none")
        ax.add_patch(c); c.set_clip_path(clip)
    for _ in range(16):                                   # fallout
        x, y, r = rng.uniform(4, cw - 4), rng.uniform(surf + 3, ch - 9), rng.uniform(0.5, 1.1)
        ax.plot([x, x], [y + r + 0.6, y + r + 3.2], color=MUTED, lw=0.6)
        ax.add_patch(Circle((x, y), r, fc=ICE, ec=INK, lw=0.5))
    ax.text(cw / 2, ch - 0.5, "plume fallout", ha="center", va="top", fontsize=FS, color=INK, bbox=BOX)
    ax.text(2.5, surf - 5, "deposit", ha="left", va="top", fontsize=FS, color=INK, bbox=BOX)
    s, sx, sy = 5.0, 50.0, 14.0                           # the sample
    ax.add_patch(Rectangle((sx, sy), s, s, fc="none", ec="#d1495b", lw=1.4))
    ax.annotate("", (-3.5, sy + s / 2), (-3.5, surf), annotation_clip=False,
                arrowprops=dict(arrowstyle="<->", color=INK, lw=0.9, shrinkA=0, shrinkB=0))
    ax.text(-5, (surf + sy + s / 2) / 2, r"$\sim$1 m", rotation=90, ha="right", va="center",
            fontsize=FS, color=INK)
    ax.annotate("", (14, 30), (14, 6), arrowprops=dict(arrowstyle="-|>", color="#edae49", lw=2.2))
    ax.text(17, 8, "heat", ha="left", va="center", fontsize=FS, color=INK, bbox=BOX)

    # ---- right: the periodic cell, tiled ----
    run = a.root / MASTER
    files, reader = make_reader(run, "sol")
    tm = step_times(str(run)); st = [snap_step(f) for f in files]
    op = opening_step(st, [tm.get(x, np.nan) for x in st])
    fl, X, Y = reader(files[st.index(op) if op is not None else 0], want=("IcePhase",))
    ice = fl["IcePhase"][:-1:3, :-1:3]                     # drop the repeated periodic node
    Lmm = float(X.max() - X.min()) * 1e3
    side = 62.0
    bx = fig.add_axes([(W - side - 3) / W, 2 / H, side / W, side / H])
    bx.imshow(np.tile(ice, (3, 3)), origin="lower", cmap=cmocean.cm.ice, vmin=0, vmax=1,
              extent=(-Lmm, 2 * Lmm, -Lmm, 2 * Lmm), interpolation="antialiased")
    bx.add_patch(Rectangle((0, 0), Lmm, Lmm, fc="none", ec="#d1495b", lw=1.4))
    bx.set_xticks([]); bx.set_yticks([])
    for sp in bx.spines.values():
        sp.set_linewidth(0.6); sp.set_color(MUTED)
    bx.set_title(rf"periodic cell, {Lmm:.0f} mm $\times$ {Lmm:.0f} mm", fontsize=FS, pad=3)
    for ya, yb in ((sy + s, 1.0), (sy, 0.0)):
        fig.add_artist(ConnectionPatch((sx + s, ya), (0, yb), "data", "axes fraction", axesA=ax, axesB=bx,
                                       color="#d1495b", lw=0.7, ls=(0, (3, 2))))

    out = a.out or a.root / "compare" / "figure_samples"
    out.mkdir(parents=True, exist_ok=True)
    for e in ("pdf", "png"):
        f = out / f"setting_panel.{e}"
        fig.savefig(f, dpi=500, transparent=True)
        if a.copy_to:
            a.copy_to.mkdir(parents=True, exist_ok=True); shutil.copyfile(f, a.copy_to / f.name)
    print(f"wrote {out}/setting_panel.pdf/.png ({W:.0f} x {H:.0f} mm)")


if __name__ == "__main__":
    main()
