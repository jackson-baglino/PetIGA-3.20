#!/usr/bin/env python3
"""Methods figure (section 2.3): the homogenization on the master simulation.

    venv_enceladus/bin/python studies/keff_sintering/figures/fig_methods_homogenization.py <campaign dir>
        [--step-near-day 10] [--zoom-center-um 1000 1000] [--zoom-um 45]
        [--corrector-dir <dir>] [--copy-to <dir>]

THE MASTER SIMULATION is phi = 0.325, seed 1702, -20 C: the packing in Fig. 3
and in the porosity gallery. Any single-run k_eff illustration uses it.

  (a) the periodic cell (ice phase field) with the zoom window marked
  (b) the zoom: the diffuse interface on the mesh -- element edges, the
      phi = 0.5 contour, and the phi = 0.01 / 0.99 contours bounding the band
  (c) the corrector t_x for a unit macroscopic gradient along x
  (d) the local heat flux magnitude |K (e_x + grad t_x)|, which shows the
      conduction paths through the necks

(c) and (d) need the corrector fields, which the solver writes with
-keff_write_corrector (igakeff.dat + t_vec_<step>_<m>.dat). Until they exist
those panels are drawn empty. Get them with ONE local replay (about a minute):

    R=<campaign dir>/packing_2D_phi0.325_Rave50um_LR40_seed1702_L2mm_eps1000nm_perxy_T-20__snow_T-20_h1.00_30d
    ./scripts/Studio/run_enceladus.sh packing_2D_phi0.325_Rave50um_LR40_seed1702_L2mm_eps1000nm_perxy_T-20 \\
        snow_T-20_h1.00_30d corrector -- -keff 1 -keff_replay "$R" -keff_write_corrector 1 \\
        -keff_csv "$R/corrector/k_eff_replay.csv" -keff_ksp_type cg -keff_pc_type gamg
    (mkdir "$R/corrector" first)

Writes methods_homogenization.{pdf,png}.
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
from pplib import step_times  # noqa: E402
from plot_keff_snapshots import make_reader, snap_step, INK, MUTED, FS, FS_SMALL, FS_TINY  # noqa: E402

MM, DAY = 1 / 25.4, 86400.0
MASTER = "packing_2D_phi0.325_Rave50um_LR40_seed1702_L2mm_eps1000nm_perxy_T-20__snow_T-20_h1.00_30d"
K_ICE, K_AIR, EPS = 2.29, 0.02, 1.0e-6


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root", type=Path)
    ap.add_argument("--step-near-day", type=float, default=10.0)
    ap.add_argument("--zoom-center-um", type=float, nargs=2, default=None)
    ap.add_argument("--zoom-um", type=float, default=70.0)
    ap.add_argument("--corrector-dir", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--copy-to", type=Path, default=None)
    a = ap.parse_args()
    plt.rcParams.update(pplib.MANUSCRIPT_RC)
    run = a.root / MASTER
    files, reader = make_reader(run, "sol")
    tm = step_times(str(run))
    st = [snap_step(f) for f in files]
    tt = np.array([tm.get(s, np.nan) for s in st])
    i = int(np.nanargmin(np.abs(tt - a.step_near_day * DAY)))
    fl, X, Y = reader(files[i], want=("IcePhase",))
    phi = fl["IcePhase"]; x = X[0, :] * 1e6; y = Y[:, 0] * 1e6
    print(f"master run, step {st[i]}, t = {tt[i] / DAY:.2f} d; grid {phi.shape}")

    # zoom window: by default a NECK -- a saddle of the ice thickness (distance
    # to the pore): a maximum across the bridge, a minimum along it.
    if a.zoom_center_um is None:
        from scipy import ndimage as ndi
        dxu = x[1] - x[0]
        edt = ndi.gaussian_filter(ndi.distance_transform_edt(phi > 0.5) * dxu, 8)
        gy, gx = np.gradient(edt)
        hyy, hyx = np.gradient(gy); hxy, hxx = np.gradient(gx)
        det = hxx * hyy - hxy * hyx
        sad = (det < 0) & (np.hypot(gx, gy) < 0.02) & (edt > 9) & (edt < 22)
        yy, xx = np.nonzero(sad)
        j = int(np.argmin(np.hypot(x[xx] - 1000, y[yy] - 1000) - 300 * det[yy, xx] / np.abs(det[yy, xx]).max()))
        cx, cy = float(x[xx[j]]), float(y[yy[j]])
        print(f"neck found: half-width {edt[yy[j], xx[j]]:.1f} um")
    else:
        cx, cy = a.zoom_center_um
    h = a.zoom_um / 2
    print(f"zoom centre ({cx:.0f}, {cy:.0f}) um, window {a.zoom_um:g} um")

    cor = None
    cdir = a.corrector_dir or run / "corrector"
    tf = cdir / f"t_vec_{st[i]:05d}_0.dat"
    if (cdir / "igakeff.dat").is_file() and tf.is_file():
        from igakit.io import PetIGA
        nrb = PetIGA().read(str(cdir / "igakeff.dat"))
        cor = np.squeeze(PetIGA().read_vec(str(tf), nrb)).T
        print(f"corrector read: {tf.name}, shape {cor.shape}")

    W, H = 170.0, 131.0
    fig = plt.figure(figsize=(W * MM, H * MM))
    s_ = 58.0
    pos = {"a": (10, 70), "b": (95, 70), "c": (10, 5), "d": (95, 5)}
    axs = {k: fig.add_axes([px / W, py / H, s_ / W, s_ / H]) for k, (px, py) in pos.items()}
    ext = (x.min(), x.max(), y.min(), y.max())
    icecm = cmocean.cm.ice

    ax = axs["a"]
    ax.imshow(phi, origin="lower", extent=ext, cmap=icecm, vmin=0, vmax=1, interpolation="antialiased")
    ax.add_patch(plt.Rectangle((cx - h, cy - h), 2 * h, 2 * h, fill=False, ec="#D55E00", lw=1.4, zorder=5))
    from matplotlib.patches import ConnectionPatch
    for (ya, yb) in ((cy + h, 1.0), (cy - h, 0.0)):
        fig.add_artist(ConnectionPatch(xyA=(cx + h, ya), coordsA=ax.transData, xyB=(0.0, yb),
                                       coordsB=axs["b"].transAxes, color="#D55E00", lw=0.7))
    ax.plot([120, 620], [110, 110], color=INK, lw=2.2, solid_capstyle="butt", zorder=6)
    ax.text(370, 140, r"500 $\mu$m", color=INK, ha="center", va="bottom", fontsize=FS_TINY, zorder=6,
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.85))

    ax = axs["b"]
    mx = (x >= cx - h) & (x <= cx + h); my = (y >= cy - h) & (y <= cy + h)
    xz, yz, pz = x[mx], y[my], phi[np.ix_(my, mx)]
    ax.pcolormesh(xz, yz, pz, cmap=icecm, vmin=0, vmax=1, shading="nearest", rasterized=True)
    dx = xz[1] - xz[0]
    for g in xz[::1]:
        ax.axvline(g - dx / 2, color="white", lw=0.12, alpha=0.5)
    for g in yz[::1]:
        ax.axhline(g - dx / 2, color="white", lw=0.12, alpha=0.5)
    ax.contour(xz, yz, pz, levels=[0.01, 0.99], colors=["#D55E00"], linewidths=0.8, linestyles="--")
    ax.contour(xz, yz, pz, levels=[0.5], colors=["#D55E00"], linewidths=1.2)
    ax.set_xlim(cx - h, cx + h); ax.set_ylim(cy - h, cy + h)
    ax.plot([cx - h + 4, cx - h + 14], [cy - h + 4, cy - h + 4], color=INK, lw=2.2, solid_capstyle="butt", zorder=6)
    ax.text(cx - h + 9, cy - h + 5.6, r"10 $\mu$m", color=INK, ha="center", va="bottom", fontsize=FS_TINY, zorder=6,
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.85))

    if cor is not None:
        tx = cor * 1e6                                    # m -> um, for a unit gradient
        v = np.percentile(np.abs(tx), 99.5)
        im = axs["c"].imshow(tx, origin="lower", extent=ext, cmap=cmocean.cm.balance, vmin=-v, vmax=v,
                             interpolation="antialiased")
        axs["c"].contour(x, y, phi, levels=[0.5], colors=[INK], linewidths=0.25)
        gy, gx = np.gradient(cor, y * 1e-6, x * 1e-6)
        # The solver's tensor law (src/keff_cell.c, KeffPointCond):
        # K = k_arith (I - n n) + k_harm n n, n = grad phi / |grad phi|,
        # phi clamped to [0, 1]; q = K (e_x + grad t_x).
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
        q = np.hypot(qx, qy)
        # Check: the cell average of q is the first row of k_eff. Drop the
        # repeated periodic row and column before averaging.
        print(f"cell average of q: ({qx[:-1, :-1].mean():.4f}, {qy[:-1, :-1].mean():.4f}) W/m/K"
              "  -- compare k_00, k_01 of this step in corrector/k_eff_replay.csv")
        im2 = axs["d"].imshow(q, origin="lower", extent=ext, cmap=cmocean.cm.thermal,
                              norm=matplotlib.colors.LogNorm(vmin=K_AIR, vmax=np.percentile(q, 99.9)),
                              interpolation="antialiased")
        for k_, im_, lab in (("c", im, r"$t_x$ [$\mu$m]"), ("d", im2, r"$|\mathbf{q}|\,/\,|\nabla\bar{T}|$ [W m$^{-1}$ K$^{-1}$]")):
            p_ = axs[k_].get_position()
            cax = fig.add_axes([p_.x1 + 0.008, p_.y0, 0.012, p_.height])
            cb = fig.colorbar(im_, cax=cax); cb.ax.tick_params(labelsize=FS_TINY, width=0.5, length=2)
            cb.set_label(lab, fontsize=FS_SMALL); cb.outline.set_linewidth(0.5)
    else:
        for k_, txt in (("c", r"corrector $t_x$"), ("d", r"heat flux $|\mathbf{q}_x|$")):
            axs[k_].text(0.5, 0.5, txt + "\n(needs the corrector replay;\nsee this script's header)",
                         transform=axs[k_].transAxes, ha="center", va="center", fontsize=FS_SMALL, color=MUTED)
    for k_, ax in axs.items():
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_linewidth(0.6); sp.set_color(MUTED)
        # Figure-level text above the zoom lines, on a small white pad: the
        # upper zoom line runs through the (b) label's position.
        fig.text(*fig.transFigure.inverted().transform(ax.transAxes.transform((-0.03, 1.0))),
                 pplib.bold(f"({k_})"), ha="right", va="top", fontsize=FS, zorder=20,
                 bbox=dict(boxstyle="circle,pad=0.12", fc="white", ec="none", alpha=0.5) if k_ == "b" else None)

    out = a.out or a.root / "compare" / "figure_samples"
    out.mkdir(parents=True, exist_ok=True)
    for e in ("pdf", "png"):
        f = out / f"methods_homogenization.{e}"
        fig.savefig(f, dpi=500, transparent=True)
        if a.copy_to:
            a.copy_to.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(f, a.copy_to / f.name)
    print(f"wrote {out}/methods_homogenization.pdf/.png")


if __name__ == "__main__":
    main()
