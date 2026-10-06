#!/usr/bin/env python3
"""Figure 1 (manuscript), assembled: (a) the two sintering mechanisms on a grain
pair, (b) the same process in an aggregate.

    venv_enceladus/bin/python studies/keff_sintering/figures/fig1_mechanisms_aggregate.py
        <campaign dir> --pair-run <grain-pair run> [--schematic-mm 95] [--copy-to <dir>]

(a) is postprocess/plot_mechanisms_figure.py's schematic, drawn by its own
build() so the two cannot drift apart; (b) is fig1_aggregate_strip.py's
close-up of the master run (phi 0.325, seed 1702, -20 C) at four instants.
170 mm wide. Writes figure1_mechanisms_aggregate.{pdf,png}.
"""
from __future__ import annotations

import argparse, shutil, sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
PROJ = HERE.parents[2]
sys.path.insert(0, str(PROJ / "postprocess"))
sys.path.insert(0, str(HERE))
import pplib  # noqa: E402
import plot_mechanisms_figure as mech  # noqa: E402
from fig1_aggregate_strip import add_strip  # noqa: E402
from plot_keff_snapshots import INK, FS  # noqa: E402

MM = 1 / 25.4


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root", type=Path)
    ap.add_argument("--pair-run", type=Path, required=True)
    ap.add_argument("--schematic-mm", type=float, default=95.0)
    ap.add_argument("--center-um", type=float, nargs=2, default=[700.0, 1450.0])
    ap.add_argument("--window-um", type=float, default=420.0)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--copy-to", type=Path, default=None)
    a = ap.parse_args()
    plt.rcParams.update(pplib.MANUSCRIPT_RC)

    # the schematic with its own defaults, at the requested width
    ma = mech.build_parser().parse_args(["--run", str(a.pair_run), "--width-mm", str(a.schematic_mm)])
    fig = mech.build(a.pair_run.resolve(), ma)
    hm = fig.get_figheight() / MM
    W, gap_ab, title, bot = 170.0, 3.0, 6.0, 1.0
    cell = (W - 2 * 1.0 - 3 * 2.0) / 4
    H = hm + gap_ab + title + cell + bot
    fig.set_size_inches(W * MM, H * MM)
    fig.axes[0].set_position([(W - a.schematic_mm) / 2 / W, (H - hm) / H, a.schematic_mm / W, hm / H])
    add_strip(fig, a.root, a.center_um, a.window_um, W, H, y_base=bot)
    fig.text(1.0 / W, (H - 4.0) / H, pplib.bold("(a)"), fontsize=FS, color=INK, ha="left", va="center")
    fig.text(1.0 / W, (bot + cell + title + 0.5) / H, pplib.bold("(b)"), fontsize=FS, color=INK,
             ha="left", va="center")

    out = a.out or a.root / "compare" / "figure_samples"
    out.mkdir(parents=True, exist_ok=True)
    for e in ("pdf", "png"):
        f = out / f"figure1_mechanisms_aggregate.{e}"
        fig.savefig(f, dpi=600, transparent=True)
        if a.copy_to:
            a.copy_to.mkdir(parents=True, exist_ok=True); shutil.copyfile(f, a.copy_to / f.name)
    print(f"wrote {out}/figure1_mechanisms_aggregate.pdf/.png ({W:.0f} x {H:.0f} mm)")


if __name__ == "__main__":
    main()
