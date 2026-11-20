#!/usr/bin/env python3
"""Ice in one run at a few times, shown vertically: x = 0 (collector or warm vapor face) at the
top, x = L (dead pore or far end) at the bottom, one panel per time.

Each panel is the ice phase field (air white, ice blue, regolith gray outside the mesh) with the
phi = 0.5 contour. The top and bottom faces are labeled from the run's own options
(read_opts): a vapor-pinned face says what it is pinned to, an unpinned face is a dead end
(zero vapor flux), a periodic run says so; the warm and cold ends come from -grad_temp0.
Checked with figcheck before saving.

    venv_cryosuction/bin/python studies/cryosuction/scripts/plot_snapshots_vertical.py \
        --run <run> --out fig.png [--days 0 10 30 60] [--n 5] [--width 7.0] [--height 4.4]
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "postprocess"))
import pplib  # noqa: E402

DAY = 86400.0


def rho_label(v: float) -> str:
    """rho_v as a multiple of rho_vs(T0); values within 1 % of 1 as 1 +/- offset, so a
    1e-5-level end condition does not print as 1.
    (Offsets of 0.005 or more print as a plain multiple, e.g. 0.991.)"""
    if v == 1.0:
        return "$\\rho_v = \\rho_{vs}(T_0)$"
    if abs(v - 1.0) < 0.005:
        e = math.floor(math.log10(abs(v - 1.0)))
        m = abs(v - 1.0) / 10 ** e
        sign = "+" if v > 1 else "-"
        return f"$\\rho_v = \\rho_{{vs}}(T_0)\\,(1 {sign} {m:.2g}\\times10^{{{e}}})$"
    return f"$\\rho_v = {v:.3g}\\,\\rho_{{vs}}(T_0)$"


def faces(o: dict):
    """(top, bottom) labels for the x = 0 and x = L faces of a run."""
    periodic = int(float(pplib.opt_float(o, "-periodic_x", 0) or 0)) == 1
    g = str(o.get("-grad_temp0", "0,0,0")).split(",")
    G = -float(g[0]) if g and g[0] else 0.0                       # K/m, > 0: warm at x = 0
    warm = ("warm", "cold") if G > 0 else (("cold", "warm") if G < 0 else ("", ""))
    if periodic:
        return "periodic (cell repeats)", "periodic (cell repeats)", G
    pinned = [False, False]
    if pplib.opt_bool(o, "-flag_BC_rhovfix"):
        fl = pplib.face_list(o, "-rhovfix_faces")
        if fl:
            pinned = [(0, 0) in fl, (0, 1) in fl]
        else:
            ax = int(float(pplib.opt_float(o, "-rhovfix_axis", -1) or -1))
            pinned = [ax in (-1, 0)] * 2
    vals = [pplib.opt_float(o, "-rhovfix_lo", pplib.opt_float(o, "-humidity", 1.0)),
            pplib.opt_float(o, "-rhovfix_hi", pplib.opt_float(o, "-humidity", 1.0))]
    collector = pinned[0] and not pinned[1]
    out = []
    for side in (0, 1):
        w = f", {warm[side]}" if warm[side] else ""
        if pinned[side]:
            name = "collector" if (collector and side == 0) else "open end"
            out.append(f"{name}: {rho_label(vals[side])}{w}")
        else:
            out.append(f"dead pore: no vapor flux{w}" if side == 1 else f"closed end: no vapor flux{w}")
    return out[0], out[1], G


def pick_steps(run: str, days, n: int):
    tmap = pplib.snapshot_times(run)
    steps = sorted(s for s in (int(f[4:-4]) for f in pplib._sol_files(run)) if s in tmap)
    t = np.array([tmap[s] for s in steps]) / DAY
    targets = list(days) if days else list(np.linspace(t[0], t[-1], n))
    chosen = []
    for d in targets:
        s = steps[int(np.argmin(np.abs(t - d)))]
        if s not in chosen:
            chosen.append(s)
    return chosen


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--days", type=float, nargs="*", default=None, help="times to show [d] (nearest snapshot)")
    ap.add_argument("--n", type=int, default=5, help="evenly spaced times when --days is not given")
    ap.add_argument("--width", type=float, default=7.0, help="figure width [in]")
    ap.add_argument("--height", type=float, default=4.4, help="figure height [in]")
    ap.add_argument("--n-per-elem", type=int, default=2)
    ap.add_argument("--depth", type=float, nargs=2, default=None, help="x range shown [µm] (default: all)")
    ap.add_argument("--across", type=float, nargs=2, default=None, help="y range shown [µm] (default: all)")
    ap.add_argument("--zoom-ice", action="store_true",
                    help="crop to the box holding ice (phi > 0.5) in any shown snapshot, plus a margin")
    a = ap.parse_args(argv)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    plt.rcParams.update({"font.size": 9})

    run = str(Path(a.run).resolve())
    o = pplib.read_opts(run)
    top, bottom, G = faces(o)
    steps = pick_steps(run, a.days, a.n)
    snaps = list(pplib.iter_snapshots(run, n_per_elem=a.n_per_elem, steps=steps))
    print(f"run {Path(run).name}\n  top (x = 0): {top}\n  bottom (x = L): {bottom}\n  G = {G:g} K/m; "
          f"snapshots at {', '.join(f'{s.t / DAY:.1f}' for s in snaps)} d")

    if a.zoom_ice:
        xs = np.concatenate([s.X[s.phi > 0.5] for s in snaps]) * 1e6
        ys = np.concatenate([s.Y[s.phi > 0.5] for s in snaps]) * 1e6
        X0, X1 = snaps[0].X.min() * 1e6, snaps[0].X.max() * 1e6
        Y0, Y1 = snaps[0].Y.min() * 1e6, snaps[0].Y.max() * 1e6
        mx = max(20.0, 0.25 * (xs.max() - xs.min()))
        my = max(20.0, 0.25 * (ys.max() - ys.min()))
        a.depth = a.depth or [max(X0, xs.min() - mx), min(X1, xs.max() + mx)]
        a.across = a.across or [max(Y0, ys.min() - my), min(Y1, ys.max() + my)]
        print(f"  zoom: depth {a.depth[0]:.0f}-{a.depth[1]:.0f} um, across {a.across[0]:.0f}-{a.across[1]:.0f} um")
    ice = LinearSegmentedColormap.from_list("ice", ["#ffffff", "#2f6fa8"])
    fig, axs = plt.subplots(1, len(snaps), figsize=(a.width, a.height), squeeze=False)
    axs = axs[0]
    for k, (ax, s) in enumerate(zip(axs, snaps)):
        h, v = pplib.to_display(s.X * 1e6, s.Y * 1e6)
        ax.set_facecolor("#9a9a9a")                              # regolith: outside the mesh
        ax.pcolormesh(h, v, np.clip(s.phi, 0, 1), cmap=ice, vmin=0, vmax=1, shading="gouraud", rasterized=True)
        ax.contour(h, v, s.phi, levels=[0.5], colors="k", linewidths=0.7)
        pplib.style_vertical(ax)
        ax.set_xlim(*(a.across or (h.min(), h.max())))
        ax.set_ylim(*((a.depth[1], a.depth[0]) if a.depth else (v.max(), v.min())))
        ax.set_title(f"{s.t / DAY:.0f} d" if s.t / DAY >= 1 else f"{s.t / DAY:.2f} d", fontsize=9)
        if k:
            ax.set_ylabel("")
            ax.tick_params(labelleft=False)
        ax.set_xlabel("")
        ax.tick_params(labelsize=8)
    axs[0].set_ylabel("depth x [µm]")
    fig.text(0.5, 0.045, "across the pore y [µm]", ha="center", va="bottom", fontsize=8)
    fig.text(0.5, 0.995, f"▲ top, x = 0: {top}", ha="center", va="top", fontsize=9)
    fig.text(0.5, 0.005, f"▼ bottom, x = L: {bottom}", ha="center", va="bottom", fontsize=9)
    fig.tight_layout(rect=(0, 0.08, 1, 0.94))
    # x ticks from the drawn panel width: a narrow, tall pore leaves room for one or two labels only
    for ax in axs:
        lo, hi = ax.get_xlim()
        w_in = ax.get_position().width * fig.get_figwidth()
        span = hi - lo
        nice = [v for v in (5, 10, 20, 25, 50, 100, 150, 200, 250, 300, 400, 500) if v <= span]
        step = nice[-1] if nice else span
        first = math.ceil(lo / step) * step
        ticks = [first] if w_in < 0.45 else [first, first + step] if first + step <= hi else [first]
        if w_in < 0.45 and first + step <= hi:
            ticks = [first + step]
        ax.set_xticks(ticks)
    from figcheck import check_figure
    for problem in check_figure(fig):
        print(f"  figcheck: {problem}")
    fig.savefig(a.out, dpi=220)
    plt.close(fig)
    print(f"wrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
