#!/usr/bin/env python3
"""make_heatflow_movie.py — animate where the heat goes as a packing sinters.

    python make_heatflow_movie.py <run_dir> [--corrector-dir DIR] [--out FILE.mp4]
        [--fps 6] [--hold-last 2.0] [--px 3840] [--stride 1] [--frame-png STEP]

Each frame: the ice (left), the local heat flux for a unit temperature
gradient along x (right, log scale, one colour scale for the whole movie),
and underneath the run's conductivity k_xx against time with a marker on the
current frame.

NEEDS THE CORRECTOR FIELDS. The heat flux is q = K (e_x + grad t_x), and t_x
is the corrector of the homogenization solve, which the solver writes only
when asked. One replay over the run's snapshots produces them (k_eff solves
only, about two minutes per snapshot on the workstation):

    R=<run_dir>;  mkdir -p "$R/corrector"
    ./scripts/Studio/run_enceladus.sh <geometry> <experiment> corrector -- \\
        -keff 1 -keff_replay "$R" -keff_write_corrector 1 \\
        -keff_csv "$R/corrector/k_eff_replay.csv" -keff_ksp_type cg -keff_pc_type gamg

which leaves igakeff.dat, t_vec_<step>_{0,1}.dat and k_eff_replay.csv in
<run_dir>/corrector/. A frame is drawn for every snapshot that has both its
sol_<step>.dat and its t_vec_<step>_0.dat, from the opening frame (t >= 1 s)
on; the initial condition is skipped, as in the other movies.

THE FLUX is computed exactly as the solver's tensor law does (src/keff_cell.c,
KeffPointCond): K = k_arith (I - n n) + k_harm n n with n = grad phi/|grad phi|
and phi clamped to [0, 1]. Its cell average is the first row of k_eff, which is
printed per frame against the replay's own k_00 as a check.

The timeline is on a log time axis, because the snapshots are log-spaced: half
of them fall in the first day, where the conductivity rises fastest.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors
import cmocean

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pplib                                              # noqa: E402
from plot_keff import load                                # noqa: E402
from plot_keff_snapshots import make_reader, snap_step    # noqa: E402
from make_packing_movie import assemble                   # noqa: E402

DAY = 86400.0
K_ICE, K_AIR = 2.29, 0.02
INK, MUTED = "#1a1a1a", "#5c5c5c"
C_MARK = "#1E88E5"
RC = {"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
      "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic",
      "mathtext.bf": "Arial:bold", "axes.unicode_minus": True}


def heat_flux(run_dir: Path, step: int, corrector_dir: Path | None = None):
    """(|q| per unit mean gradient [W/m/K], phi, x_um, y_um, (<qx>, <qy>))."""
    from igakit.io import PetIGA
    files, reader = make_reader(run_dir, "sol")
    fn = [f for f in files if snap_step(f) == step][0]
    fl, X, Y = reader(fn, want=("IcePhase",))
    phi = fl["IcePhase"]; x = X[0] * 1e6; y = Y[:, 0] * 1e6
    cdir = corrector_dir or run_dir / "corrector"
    nrb = PetIGA().read(str(cdir / "igakeff.dat"))
    cor = np.squeeze(PetIGA().read_vec(str(cdir / f"t_vec_{step:05d}_0.dat"), nrb)).T
    gy, gx = np.gradient(cor, y * 1e-6, x * 1e-6)
    pc = np.clip(phi, 0.0, 1.0)
    ka = pc * K_ICE + (1 - pc) * K_AIR
    kh = 1.0 / (pc / K_ICE + (1 - pc) / K_AIR)
    py_, px_ = np.gradient(pc, y * 1e-6, x * 1e-6)
    g2 = px_ ** 2 + py_ ** 2
    Ex, Ey = 1.0 + gx, gy
    with np.errstate(invalid="ignore", divide="ignore"):
        proj = np.where(g2 > 0, (px_ * Ex + py_ * Ey) / g2, 0.0)
    qx = ka * Ex + (kh - ka) * proj * px_
    qy = ka * Ey + (kh - ka) * proj * py_
    # drop the repeated periodic row and column before averaging
    return np.hypot(qx, qy), phi, x, y, (float(qx[:-1, :-1].mean()), float(qy[:-1, :-1].mean()))


def time_label(t):
    if t < 3600:
        return f"{t / 60:.0f} min" if t >= 60 else f"{t:.0f} s"
    if t < 2 * DAY:
        return f"{t / 3600:.1f} hours"
    return f"{t / DAY:.1f} days"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--corrector-dir", type=Path, default=None,
                    help="where igakeff.dat and t_vec_*.dat are (default <run_dir>/corrector)")
    ap.add_argument("--out", type=Path, default=None,
                    help="output mp4 (default <run_dir>/plots/heatflow_movie.mp4)")
    ap.add_argument("--fps", type=int, default=6)
    ap.add_argument("--hold-last", type=float, default=2.0, help="seconds to hold the final frame")
    ap.add_argument("--t-min-hours", type=float, default=0.5,
                    help="skip snapshots earlier than this, except the opening frame "
                         "(default 0.5 h: the log-spaced early snapshots show no change)")
    ap.add_argument("--px", type=int, default=3840,
                    help="frame width in pixels (default 3840 = 4K: each panel is then 1416 px, half the "
                         "2831-node field; at 1920 the panels are only 708 px and the grains look soft)")
    ap.add_argument("--stride", type=int, default=1, help="field down-sampling for drawing")
    ap.add_argument("--crf", type=int, default=16)
    ap.add_argument("--frame-png", type=int, default=None,
                    help="render only this step to a PNG (preview) and exit")
    ap.add_argument("--label", default=None,
                    help="text under the panels (default: temperature, porosity, cell size from the opts)")
    a = ap.parse_args()
    plt.rcParams.update(RC)

    run = a.run_dir.resolve()
    cdir = a.corrector_dir or run / "corrector"
    if not (cdir / "igakeff.dat").is_file():
        sys.exit(f"no corrector fields in {cdir} -- see this script's header for the replay that writes them")
    files, _ = make_reader(run, "sol")
    tm = pplib.step_times(str(run))
    steps = [s for s in (snap_step(f) for f in files)
             if (cdir / f"t_vec_{s:05d}_0.dat").is_file() and tm.get(s, 0.0) >= 1.0]
    if not steps:
        sys.exit(f"no snapshot in {run} has a t_vec_<step>_0.dat in {cdir}")
    steps = steps[:1] + [s_ for s_ in steps[1:] if tm[s_] >= a.t_min_hours * 3600.0]
    print(f"  {len(steps)} frame(s): steps {steps[0]}..{steps[-1]}, "
          f"t = {time_label(tm[steps[0]])} .. {time_label(tm[steps[-1]])}")

    d = load(run)                                  # the run's own k_eff(t), for the timeline
    rep = None
    if (cdir / "k_eff_replay.csv").is_file():
        rep = np.atleast_1d(np.genfromtxt(cdir / "k_eff_replay.csv", delimiter=",", names=True))
    opts = pplib.read_opts(str(run))
    T = pplib.opt_float(opts, "-temp", float("nan"))
    Lx = pplib.opt_float(opts, "-Lx", float("nan"))
    label = a.label or (f"{T:g} °C".replace("-", "−") + f"   ·   {Lx * 1e3:g} mm × {Lx * 1e3:g} mm periodic cell"
                        "   ·   temperature gradient imposed left to right")

    # one colour scale for the whole movie, set by the last (most sintered) frame
    q_end = heat_flux(run, steps[-1], cdir)[0]
    norm = matplotlib.colors.LogNorm(vmin=K_AIR, vmax=float(np.percentile(q_end[::4, ::4], 99.9)))

    W, H = 12.8, 7.2                               # 16:9
    dpi = a.px / W
    s = 4.72                                       # panel side [in]
    x1, gap, yb = 0.86, 0.35, 1.50
    out = a.out or run / "plots" / "heatflow_movie.mp4"
    frames = out.parent / (out.stem + "_frames")
    frames.mkdir(parents=True, exist_ok=True)
    todo = [a.frame_png] if a.frame_png is not None else steps

    for step in todo:
        q, phi, x, y, (qx_m, qy_m) = heat_flux(run, step, cdir)
        t = tm[step]
        k_rep = float(rep["k_00"][rep["step"] == step][0]) if rep is not None and (rep["step"] == step).any() else qx_m
        print(f"  step {step:5d}  t = {time_label(t):>10s}   <q_x> = {qx_m:.4f}   k_00 (replay) = {k_rep:.4f} W/m/K")
        fig = plt.figure(figsize=(W, H), facecolor="white")
        F = lambda x_, y_, w_, h_: (x_ / W, y_ / H, w_ / W, h_ / H)
        ext = (x.min() / 1e3, x.max() / 1e3, y.min() / 1e3, y.max() / 1e3)
        st = a.stride
        ax1 = fig.add_axes(F(x1, yb, s, s)); ax2 = fig.add_axes(F(x1 + s + gap, yb, s, s))
        ax1.imshow(phi[::st, ::st], origin="lower", extent=ext, cmap=cmocean.cm.ice, vmin=0, vmax=1,
                   interpolation="antialiased")
        im = ax2.imshow(q[::st, ::st], origin="lower", extent=ext, cmap=cmocean.cm.thermal, norm=norm,
                        interpolation="antialiased")
        for ax, ttl in ((ax1, "ice"), (ax2, "heat flow")):
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_linewidth(1.2); sp.set_color(INK)
            ax.set_title(ttl, fontsize=15, color=INK, pad=5)
        for lw, c in ((6.0, "white"), (3.2, INK)):                # 0.5 mm scale bar
            ax1.plot([ext[0] + 0.12, ext[0] + 0.62], [ext[2] + 0.13] * 2, color=c, lw=lw, solid_capstyle="butt")
        ax1.text(ext[0] + 0.37, ext[2] + 0.18, "0.5 mm", ha="center", va="bottom", fontsize=12, color=INK,
                 bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
        cax = fig.add_axes(F(x1 + 2 * s + gap + 0.18, yb + 0.02 * s, 0.17, 0.96 * s))
        cb = fig.colorbar(im, cax=cax); cb.ax.tick_params(labelsize=12, width=0.8, length=3)
        cb.set_label(r"heat flux / gradient  [W m$^{-1}$ K$^{-1}$]", fontsize=13)
        cb.outline.set_linewidth(0.8)

        # timeline: k_xx(t) on a log time axis, the marker on this frame
        axk = fig.add_axes(F(x1, 0.50, 2 * s + gap, 0.62))
        # the axis starts where the frames after the opening one do
        t_lo = tm[steps[1]] if len(steps) > 1 and a.t_min_hours > 0 else max(tm[steps[0]], 1.0)
        w = d["t"] >= 1.0
        axk.plot(d["t"][w] / DAY, d["kxx"][w], color=MUTED, lw=1.8)
        axk.plot([max(t / DAY, 0.8 * t_lo / DAY)], [k_rep], "o", ms=13, mfc=C_MARK, mec="white",
                 mew=1.6, zorder=5, clip_on=False)
        axk.set_xscale("log"); axk.set_xlim(0.8 * t_lo / DAY, tm[steps[-1]] / DAY * 1.03)
        tk = [(60, "1 min"), (3600, "1 hour"), (6 * 3600, "6 hours"), (DAY, "1 day"), (3 * DAY, "3 days"),
              (10 * DAY, "10 days"), (30 * DAY, "30 days")]
        tk = [(v / DAY, l) for v, l in tk if axk.get_xlim()[0] <= v / DAY <= axk.get_xlim()[1]]
        axk.set_xticks([v for v, _ in tk]); axk.set_xticklabels([l for _, l in tk])
        axk.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        axk.tick_params(labelsize=12, width=0.8, length=4)
        axk.text(0.0, 1.10, r"conductivity  $k_\mathrm{eff}$  [W m$^{-1}$ K$^{-1}$]", transform=axk.transAxes,
                 ha="left", va="bottom", fontsize=12, color=INK)
        for sp in ("top", "right"):
            axk.spines[sp].set_visible(False)
        fig.text((x1 + s + gap / 2) / W, (H - 0.30) / H,
                 f"t = {time_label(t)}        " + rf"$k_\mathrm{{eff}}$ = {k_rep:.2f} W m$^{{-1}}$ K$^{{-1}}$",
                 ha="center", va="center", fontsize=20, color=INK)
        fig.text((x1) / W, 0.13 / H, label, ha="left", va="center", fontsize=11.5, color=MUTED)
        fn = (out.parent / f"{out.stem}_step{step:05d}.png") if a.frame_png is not None \
            else frames / f"frame_{step:05d}.png"
        fig.savefig(fn, dpi=dpi)
        plt.close(fig)
        if a.frame_png is not None:
            print(f"  wrote {fn}")
            return 0

    # hold the last frame: repeat it under later names so the concat keeps the order
    last = frames / f"frame_{steps[-1]:05d}.png"
    for i in range(int(round(a.hold_last * a.fps))):
        (frames / f"frame_{steps[-1]:05d}_hold{i:03d}.png").write_bytes(last.read_bytes())
    return assemble(frames, out, a.fps, a.crf)


if __name__ == "__main__":
    sys.exit(main())
