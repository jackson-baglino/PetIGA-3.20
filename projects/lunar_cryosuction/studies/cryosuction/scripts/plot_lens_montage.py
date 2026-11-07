"""Initial and final ice of every run in a batch, one panel per run.

The phi = 0.5 contour at the first snapshot (gray, dashed) and the last (black),
over the final phi field. Panels share the x scale; y is to scale in each panel.
Checked with figcheck before saving.

    venv_cryosuction/bin/python studies/cryosuction/scripts/plot_lens_montage.py \
        --batch <batch> --out <png> [--xlim 120 380] [--label-strip channel_2D_ ...]
"""
from __future__ import annotations

import argparse
import glob
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "postprocess"))
import pplib  # noqa: E402


def first_last(run):
    snaps = sorted(glob.glob(os.path.join(run, "sol_*.dat")))
    n = len(snaps)
    got = list(pplib.iter_snapshots(run, every=max(n - 1, 1)))
    return got[0], got[-1]


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--batch", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--xlim", type=float, nargs=2, default=(120.0, 380.0))
    p.add_argument("--label-strip", nargs="*", default=[])
    a = p.parse_args(argv)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from figcheck import check_figure
    runs = sorted(d for d in glob.glob(os.path.join(a.batch, "*__*")) if os.path.isfile(os.path.join(d, "igasol.dat")))
    ncol = 2
    nrow = int(np.ceil(len(runs) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(7.5, 1.2 * nrow + 0.5), constrained_layout=True)
    for ax, run in zip(axs.flat, runs):
        s0, s1 = first_last(run)
        ax.pcolormesh(s1.X * 1e6, s1.Y * 1e6, s1.phi, shading="auto", cmap="Blues", vmin=0, vmax=1.4)
        ax.contour(s0.X * 1e6, s0.Y * 1e6, s0.phi, [0.5], colors="0.45", linewidths=0.8, linestyles="--")
        ax.contour(s1.X * 1e6, s1.Y * 1e6, s1.phi, [0.5], colors="k", linewidths=0.9)
        name = os.path.basename(run)
        for t in a.label_strip:
            name = name.replace(t, "")
        ax.set_title(f"{name}  (t = {s1.t / 86400:.0f} d)", fontsize=8)
        ax.set_xlim(*a.xlim)
        ax.set_aspect("equal")
        ax.tick_params(labelsize=7)
        print(f"{os.path.basename(run)}: t0 {s0.t:.0f} s, t1 {s1.t:.0f} s")
    for ax in list(axs.flat)[len(runs):]:
        ax.axis("off")
    fig.supxlabel("x [µm]   (dashed gray: start; black: end; color: final φ)", fontsize=8)
    fig.supylabel("y [µm]", fontsize=8)
    bad = check_figure(fig)
    fig.savefig(a.out, dpi=200)
    print("figcheck:", bad or "OK", "->", a.out)
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
