#!/usr/bin/env python3
"""How much do the three packing 'problems' actually bias the results?

The design goal is not a perfect microstructure -- 2D forbids one -- but the
LEAST BIASED one that still satisfies the homogenization assumptions. That
turns the three problems from judgement calls into numbers, all measurable at
t = 0 without a solver run.

The lever is eps. Every problem is a diffuse-band artifact, so sweeping eps on
a FIXED packing separates artifact from material: whatever changes with eps is
the artifact, and the eps -> 0 limit is the material.

  k_eff(eps)   spurious solid bridges add conduction paths, so k_eff reads
               high. The slope says by how much at the eps in use.
  SSA(eps)     a bridged throat has no surface, so SSA reads low. SSA is
               exactly int|grad phi| by the co-area formula, which for a fixed
               sharp geometry is eps-INDEPENDENT -- so any drift is bridging,
               with nothing else it could be.
  D_eff(eps)   effective vapour diffusivity with coefficient phi_air: the
               continuous replacement for a percolation test whose answer is
               always 'disconnected'. It measures how much vapour still moves
               through the thin films, which is the quantity that actually
               matters for whether pore fragmentation changes the physics.

CONVERGENCE, WHICH IS NOT THE SAME FOR ALL THREE. Checked on the LR20 packing
at eps = 45 nm over rasters 384 / 512 / 768 / 1024:

    k_eff   0.6330  0.6337  0.6344  0.6347    converged, 0.3% spread
    SSA     394544  399689  404192  406017    converged to ~1-2%
    D_eff   7.6e-3  4.2e-3  3.0e-3  1.4e-3    NOT CONVERGED, still falling

So k_eff and SSA are quantitative here and D_eff IS AN UPPER BOUND. D_eff is
set by the narrowest constrictions, and coarse pixels bridge them, so every
refinement removes transport; 1024 is simply the finest raster a direct solve
fits in memory. It is not a floor artifact -- transport is film-carried, not
through the ice -- which was checked separately: at N = 1024 the floor ladder
1e-4 ... 1e-8 gives 5.1e-3, 1.9e-3, 1.40e-3, 1.33e-3, 1.325e-3, i.e. converged
in the floor by 1e-7. Quote D_eff as "at most", never as a value.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

PROJ = Path(__file__).resolve().parents[2]      # projects/enceladus_DSM
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(PROJ / "preprocess"))
from cell_solve import k_eff                                # noqa: E402

K_ICE, K_AIR = 2.29, 0.02


def load(d: Path):
    m = json.loads((d / "metadata.json").read_text())
    rows = [l.split() for l in (d / "grains.dat").read_text().splitlines()
            if l and not l.startswith("#")]
    raw = np.array([[float(v) for v in r] for r in rows[1:]])
    c, r = raw[:, :2], raw[:, 2]
    k = ((c[:, 0] >= 0) & (c[:, 0] < m["Lx"])
         & (c[:, 1] >= 0) & (c[:, 1] < m["Ly"]))
    return m, c[k].copy(), r[k].copy()


def phi_field(cen, rad, Lx, Ly, eps, n):
    """Additive phi on an n x n periodic raster, as -ic_grain_union 0 builds it."""
    xs = (np.arange(n) + 0.5) * Lx / n
    ys = (np.arange(n) + 0.5) * Ly / n
    X, Y = np.meshgrid(xs, ys)
    phi = np.zeros_like(X)
    reach = 12.0 * eps
    for (cx, cy), R in zip(cen, rad):
        for ox in (-Lx, 0.0, Lx):
            for oy in (-Ly, 0.0, Ly):
                px, py = cx + ox, cy + oy
                if (px + R + reach < 0 or px - R - reach > Lx
                        or py + R + reach < 0 or py - R - reach > Ly):
                    continue
                r = np.hypot(X - px, Y - py)
                phi += 0.5 - 0.5 * np.tanh(0.5 * (r - R) / eps)
    return np.clip(phi, 0.0, 1.0)


def ssa(phi, h):
    """Interface length per unit area = int |grad phi|, by the co-area formula.

    Exact for a fixed sharp geometry at ANY eps -- the profile integrates to
    the jump in phi, which is 1, regardless of how wide the band is. So this
    number CANNOT drift with eps unless interface is being destroyed, which is
    what a bridged throat does.
    """
    gx = (np.roll(phi, -1, axis=1) - np.roll(phi, 1, axis=1)) / (2 * h)
    gy = (np.roll(phi, -1, axis=0) - np.roll(phi, 1, axis=0)) / (2 * h)
    return float(np.mean(np.hypot(gx, gy)))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--packing", type=Path,
                    default=HERE / "packings/LR20_phi0.325_seed1")
    # Ladder floor set by RESOLUTION, not by taste: the raster must resolve the
    # band, or the sweep measures pixellation rather than interface width. The
    # rule of thumb below is band/h >= 4, and the run prints band/h per row so
    # it can be checked rather than trusted.
    ap.add_argument("--eps", type=float, nargs="+",
                    default=[90e-9, 63e-9, 45e-9, 32e-9, 22.5e-9])
    ap.add_argument("--raster", type=int, default=1024)
    ap.add_argument("--floor", type=float, default=1e-6,
                    help="phi_air floor for the vapour solve; the ice interior "
                         "is genuinely ~0 there and CG needs a positive "
                         "coefficient. Sensitivity is reported.")
    ap.add_argument("--out", type=Path, default=HERE / "bias.csv")
    a = ap.parse_args()

    m, cen, rad = load(a.packing)
    Lx, Ly = m["Lx"], m["Ly"]
    h = Lx / a.raster
    rows = []
    print(f"packing {a.packing.name}   raster {a.raster}   "
          f"R_ave {m['mean_r_m_requested']*1e6:.2f} um")
    print(f"{'eps[nm]':>8} {'band/h':>7} {'band/R':>7} {'k_xx':>8} {'k_yy':>8} "
          f"{'k_iso':>8} {'SSA[1/m]':>11} {'D/D0 x':>9} {'D/D0 y':>9}")
    for eps in a.eps:
        phi = phi_field(cen, rad, Lx, Ly, eps, a.raster)
        K = K_AIR + (K_ICE - K_AIR) * phi
        ke = k_eff(K, h)
        air = np.maximum(1.0 - phi, a.floor)
        de = k_eff(air, h)                      # coefficient IS phi_air
        s = ssa(phi, h)
        rows.append({
            "eps_m": eps, "band_per_Rave": 9.2 * eps / m["mean_r_m_requested"],
            "porosity": float(1.0 - phi.mean()),
            "k_xx": ke[0, 0], "k_yy": ke[1, 1], "k_iso": 0.5 * (ke[0, 0] + ke[1, 1]),
            "k_anis": ke[1, 1] / ke[0, 0],
            "ssa_per_m": s,
            "D_xx_over_D0": de[0, 0], "D_yy_over_D0": de[1, 1],
        })
        r = rows[-1]
        r["band_per_h"] = 9.2 * eps / h
        flag = "" if r["band_per_h"] >= 4.0 else "  <-- UNDER-RESOLVED"
        print(f"{eps*1e9:8.2f} {r['band_per_h']:7.2f} {r['band_per_Rave']:7.3f} "
              f"{ke[0,0]:8.4f} {ke[1,1]:8.4f} {r['k_iso']:8.4f} {s:11.0f} "
              f"{de[0,0]:9.2e} {de[1,1]:9.2e}{flag}")

    with a.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
