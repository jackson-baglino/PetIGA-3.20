#!/usr/bin/env python3
"""render_snapshots.py — high-resolution microstructure PNGs at four instants.

    python3 render_snapshots.py --dir <run> [--fractions 0 0.333 0.667 1] [--dpi-max 900]

One standalone image per instant, written to <run>/plots/snapshots/:

    snap_1_t0.png        the opening frame (pplib's rule: the first snapshot
                         with 1 s <= t <= 1 h, which the movies and the k_eff
                         figures call t = 0)
    snap_2_t1of3.png     the snapshot nearest t_final / 3
    snap_3_t2of3.png     the snapshot nearest 2 t_final / 3
    snap_4_tfinal.png    the last snapshot
    snapshots.txt        step, time and file of each

Each image is the ice phase field over the supersaturation sigma in the pores,
in the colours of plot_keff_snapshots.py, with the same symmetric sigma scale
on all four so they can be compared, a scale bar, and the run's porosity,
temperature, seed and time in the title. The dpi is chosen so the image is at
least as fine as the mesh (one pixel per element, capped by --dpi-max).

Runs on ONE core, after the simulation, as its own SLURM job
(scripts/HPC/render_snapshots_job.sh, submitted by submit_batch.sh with an
afterany dependency), so the simulation's cores are not held for it.
Reads sol_*.dat directly; needs igasol.dat.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pplib  # noqa: E402
from pplib import step_times, opening_step  # noqa: E402
from plot_keff_snapshots import (make_reader, snap_step, _field, _scalebar,  # noqa: E402
                                 _colorbars, WANT, SIGMA_SCALE, centered_cmap,
                                 ice_alpha_cmap, INK)
import cmocean  # noqa: E402

DAY = 86400.0
NAMES = ("snap_1_t0", "snap_2_t1of3", "snap_3_t2of3", "snap_4_tfinal")


def describe(name: str) -> str:
    m = re.search(r"phi([\d.]+)_.*?(?:LR(\d+)_)?seed(\d+)_.*?_T(-?\d+)__", name)
    if not m:
        return name[:60]
    lr = f", L/R {m.group(2)}" if m.group(2) else ""
    return f"φ = {m.group(1)}, T = {m.group(4)} °C, seed {m.group(3)}{lr}"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True, type=Path)
    ap.add_argument("--fractions", type=float, nargs="+", default=[0.0, 1 / 3, 2 / 3, 1.0],
                    help="instants as fractions of t_final (0 = opening frame, 1 = last snapshot)")
    ap.add_argument("--dpi-max", type=int, default=900)
    ap.add_argument("--width-in", type=float, default=6.5)
    a = ap.parse_args(argv)
    run = a.dir.resolve()

    files, reader = make_reader(run, "sol")
    if not files:
        print(f"  no sol_*.dat in {run}; nothing to render")
        return 0
    tmap = step_times(str(run))
    steps = [snap_step(f) for f in files]
    times = [tmap.get(s, np.nan) for s in steps]
    known = [(s, t, f) for s, t, f in zip(steps, times, files) if np.isfinite(t)]
    if not known:
        print("  no snapshot times (SSA_evo.dat / outp.txt missing)", file=sys.stderr)
        return 1
    ks, kt, kf = zip(*known)
    op = opening_step(ks, kt)
    i_open = ks.index(op) if op is not None else 0
    # Thirds of the RUN, not of the snapshot list: if the last snapshot fell
    # short of t_final (a rejected final step ends the run early), the targets
    # still sit at t_final/3 and 2 t_final/3.
    t0, t_end = kt[i_open], max(max(tmap.values()), kt[-1])

    pick = []
    for fr in a.fractions:
        if fr <= 0:
            i = i_open
        elif fr >= 1:
            i = len(ks) - 1
        else:
            target = t0 + fr * (t_end - t0)
            i = int(np.argmin(np.abs(np.array(kt) - target)))
        pick.append(i)

    # One sigma scale for all the frames: symmetric, +-min(|min|, |max|) over
    # the pores, as in the k_eff snapshot figures.
    frames, pore = [], []
    for i in pick:
        fl, X, Y = reader(kf[i], want=WANT)
        frames.append((fl, X, Y, kt[i], ks[i]))
        s = SIGMA_SCALE * pplib.supersaturation(fl["VaporDensity"], fl["Temperature"])
        pore.append(s[fl["IcePhase"] < 0.5])
    pore = np.concatenate(pore)
    smin, smax = float(pore.min()), float(pore.max())
    v = min(abs(smin), abs(smax)) or max(abs(smin), abs(smax))
    norm = AsinhNorm(linear_width=max(v / 300.0, 1e-12), vmin=-v, vmax=v)
    sig_extend = {(True, True): "both", (True, False): "min",
                  (False, True): "max", (False, False): "neither"}[(smin < -v, smax > v)]
    vapcm = centered_cmap(cmocean.cm.balance, norm)
    icecm = ice_alpha_cmap()

    out = run / "plots" / "snapshots"
    out.mkdir(parents=True, exist_ok=True)
    nx = frames[0][0]["IcePhase"].shape[1]
    W = a.width_in
    ax_in = W * 0.92
    dpi = int(min(a.dpi_max, max(300, np.ceil(nx / ax_in))))
    label = describe(run.name)
    lines = [f"# {run.name}", "# file  step  time_s  time_d"]
    names = NAMES if len(pick) == 4 else [f"snap_{j + 1}" for j in range(len(pick))]
    for (fl, X, Y, t, s), name in zip(frames, names):
        H = W * 1.085
        fig = plt.figure(figsize=(W, H))
        ah = 0.92 * W / H
        ax = fig.add_axes([0.04, 0.02, 0.92, ah])
        cax_i = fig.add_axes([0.14, 0.02 + ah + 0.035, 0.30, 0.018])
        cax_s = fig.add_axes([0.60, 0.02 + ah + 0.035, 0.36, 0.018])
        XX, YY = _field(ax, fl, X, Y, norm, vapcm, icecm)
        _scalebar(ax, XX, YY)
        _colorbars(fig, cax_i, cax_s, norm, vapcm, sig_extend)
        tlab = "t = 0 (opening frame, %.3g s)" % t if name.endswith("t0") else f"t = {t / DAY:.2f} d"
        fig.text(0.5, 0.982, f"{label}    {tlab}", ha="center", va="center",
                 fontsize=11, color=INK)
        f = out / f"{name}.png"
        fig.savefig(f, dpi=dpi)
        plt.close(fig)
        lines.append(f"{f.name}  {s}  {t:.6e}  {t / DAY:.4f}")
        print(f"  wrote {f}  (step {s}, t = {t / DAY:.3f} d, {f.stat().st_size / 1e6:.1f} MB, {dpi} dpi)")
    (out / "snapshots.txt").write_text("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
