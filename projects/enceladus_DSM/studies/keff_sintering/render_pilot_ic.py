#!/usr/bin/env python3
"""Render the pilot's initial condition as the SOLVER sees it.

The generator's own preview.png draws sharp discs. That is the geometry as
specified, not the field the run starts from: the solver initialises phi by
the additive Molaro convention (-ic_grain_union 0),

    phi = clamp( sum_k [ 0.5 - 0.5 tanh(0.5 (r_k - R_k)/eps) ], 0, 1 )

with eps = 1 um here. At the 2 mm domain scale the 9.2*eps = 9.2 um band is
0.46% of the frame and invisible, so the bottom row zooms in far enough to
show it.

Writes pilot_ic.png next to this file.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

HERE = Path(__file__).parent
PROJ = HERE.parents[1]
sys.path.insert(0, str(PROJ / "preprocess"))
sys.path.insert(0, str(PROJ / "studies/packing_design"))
import figstyle as fs                                       # noqa: E402
from measure_bias import load, phi_field                    # noqa: E402

EPS = 1.0e-6
BAND = 9.2 * EPS
RAVE = 50e-6
PACKDIR = PROJ / "inputs/packings/pilot_LR40"
N_FULL = 1100
N_ZOOM = 420


def main() -> int:
    import packing_lib as pl
    seeds = [1, 2, 3, 4]
    fig = plt.figure(figsize=(10.0, 7.4))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 0.95], hspace=0.42,
                          wspace=0.28, left=0.07, right=0.985,
                          top=0.90, bottom=0.135)
    gtop = gs[0, :].subgridspec(1, 4, wspace=0.10)

    # ---- top row: the four initial conditions, full domain ----------------
    for col, s_ in enumerate(seeds):
        d = PACKDIR / f"pilot_phi0.325_Rave50um_LR40_seed{s_}"
        m, cen, rad = load(d)
        Lx, Ly = m["Lx"], m["Ly"]
        phi = phi_field(cen, rad, Lx, Ly, EPS, N_FULL)
        ax = fig.add_subplot(gtop[0, col])
        ax.imshow(phi, origin="lower", extent=[0, Lx * 1e3, 0, Ly * 1e3],
                  cmap="bone_r", vmin=0, vmax=1, interpolation="nearest")
        ax.set_aspect("equal")
        ax.set_title(f"seed {s_}   N={m['n_grains']}", fontsize=fs.FS_C_LABEL,
                     color=fs.INK, pad=4)
        ax.set_xlabel(r"x  [mm]" + "\n"
                      + f"$\\phi$={m['porosity_achieved']:.3f}  "
                        f"$Z_{{band}}$={m['coordination_at_band']:.2f}",
                      fontsize=fs.FS_NOTE)
        ax.set_xticks([0, 1, 2])
        if col == 0:
            ax.set_ylabel("y  [mm]", fontsize=fs.FS_C_LABEL)
            ax.set_yticks([0, 1, 2])
        else:
            ax.set_yticks([])
        ax.tick_params(labelsize=fs.FS_C_TICK)
        if col == 0:
            keep = (cen, rad, Lx, Ly)

    # ---- bottom row: zoom in on a NEAR-TANGENT CONTACT --------------------
    # Centring on a grain would show an isolated disc and none of the physics.
    # The interesting object is a throat the band can weld shut, so find a pair
    # whose surface gap sits inside the band and zoom on the gap itself.
    cen, rad, Lx, Ly = keep
    bonds = pl.delaunay_bonds(cen, Lx, Ly, True, True)
    i, j = bonds[:, 0], bonds[:, 1]
    dd = cen[j] + bonds[:, 2:4] * np.array([Lx, Ly]) - cen[i]
    gap = np.hypot(dd[:, 0], dd[:, 1]) - (rad[i] + rad[j])
    cand = np.flatnonzero((gap > 0.2 * BAND) & (gap < 0.7 * BAND)
                          & (rad[i] > 0.7 * RAVE) & (rad[j] > 0.7 * RAVE))
    k = int(cand[len(cand) // 2])
    pa = cen[i[k]]
    pb = cen[j[k]] + bonds[k, 2:4] * np.array([Lx, Ly])
    cx, cy = 0.5 * (pa + pb)
    the_gap = gap[k]

    for col, half in enumerate((300e-6, 90e-6, 22e-6)):
        x0, y0 = cx - half, cy - half
        w = h = 2 * half
        phi = phi_field_window(cen, rad, Lx, Ly, EPS, x0, y0, w, h, N_ZOOM)
        ax = fig.add_subplot(gs[1, col])
        ext = [-half * 1e6, half * 1e6, -half * 1e6, half * 1e6]
        if col < 2:
            ax.imshow(phi, origin="lower", extent=ext, cmap="bone_r",
                      vmin=0, vmax=1, interpolation="nearest")
            ax.plot(0, 0, "o", mfc="none", mec=fs.C[1], mew=2.0, ms=16)
        else:
            cmap = ListedColormap([fs.C[1], "#dcdcdc", fs.C[0]])
            ax.imshow(np.digitize(phi, [0.01, 0.99]), origin="lower",
                      extent=ext, cmap=cmap, vmin=0, vmax=2,
                      interpolation="nearest")
            x1 = ext[0] + 0.08 * w * 1e6
            yb = ext[2] + 0.10 * h * 1e6
            ax.plot([x1, x1 + BAND * 1e6], [yb, yb], "-", color="k", lw=3.5,
                    solid_capstyle="butt")
            ax.text(x1 + BAND * 1e6 / 2, yb + 0.035 * h * 1e6,
                    f"band = {BAND*1e6:.1f} " + r"$\mu$m", ha="center",
                    va="bottom", fontsize=fs.FS_NOTE, color="k")
            ax.text(0.04, 0.96,
                    "grey crosses the\nthroat unbroken:\n"
                    r"$\phi$ never reaches 0",
                    transform=ax.transAxes, ha="left", va="top",
                    fontsize=fs.FS_NOTE, color=fs.INK,
                    bbox=dict(boxstyle="round,pad=0.28", fc="white",
                              ec="none", alpha=0.9))
        ax.set_aspect("equal")
        ax.set_xlabel(r"x  [$\mu$m]", fontsize=fs.FS_C_LABEL)
        if col == 0:
            ax.set_ylabel(r"y  [$\mu$m]", fontsize=fs.FS_C_LABEL)
        ax.tick_params(labelsize=fs.FS_C_TICK)
        titles = [r"zoom 600 $\mu$m  ($R_{ave}$ = 50 $\mu$m)",
                  r"zoom 180 $\mu$m  — on the ringed throat",
                  f"zoom 44 " + r"$\mu$m  — gap "
                  + f"{the_gap/BAND:.2f}" + r"$\times$band"]
        ax.set_title(titles[col], fontsize=fs.FS_NOTE + 0.5, color=fs.INK,
                     loc="left", pad=5)

    handles = [plt.Rectangle((0, 0), 1, 1, fc=fs.C[0], label=r"ice  $\phi>0.99$"),
               plt.Rectangle((0, 0), 1, 1, fc="#dcdcdc", label="diffuse band"),
               plt.Rectangle((0, 0), 1, 1, fc=fs.C[1], label=r"pore  $\phi<0.01$")]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               fontsize=fs.FS_LEG, bbox_to_anchor=(0.5, 0.005))
    fig.suptitle("Pilot initial condition — what the solver starts from "
                 r"($\varepsilon$ = 1 $\mu$m, $L$ = 2 mm, $L/R_{ave}$ = 40)",
                 fontsize=fs.FS_TITLE, x=0.02, ha="left", y=0.975)
    fs.save(fig, HERE, "pilot_ic", dpi=185)
    return 0


def phi_field_window(cen, rad, Lx, Ly, eps, x0, y0, w, h, n):
    """Additive phi on an arbitrary window, with periodic images."""
    X, Y = np.meshgrid(np.linspace(x0, x0 + w, n), np.linspace(y0, y0 + h, n))
    phi = np.zeros_like(X)
    reach = 12.0 * eps
    for (gx, gy), R in zip(cen, rad):
        for ox in (-Lx, 0.0, Lx):
            for oy in (-Ly, 0.0, Ly):
                px, py = gx + ox, gy + oy
                if (px + R + reach < x0 or px - R - reach > x0 + w
                        or py + R + reach < y0 or py - R - reach > y0 + h):
                    continue
                phi += 0.5 - 0.5 * np.tanh(0.5 * (np.hypot(X - px, Y - py) - R) / eps)
    return np.clip(phi, 0.0, 1.0)


if __name__ == "__main__":
    raise SystemExit(main())
