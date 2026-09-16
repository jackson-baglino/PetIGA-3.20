#!/usr/bin/env python3
"""Is the +33% k_eff bias real, or is it a sub-REV sample-size artifact?

THE QUESTION. Spurious solid bridges across sub-band gaps raise k_eff. If that
+33% comes from a handful of bridges in a domain too small to average them out,
it is a sample-size artifact and would shrink as the domain grows. If instead
it is a bulk property of the microstructure at this eps, it will be the same at
every domain size, with only the seed-to-seed SCATTER shrinking as 1/sqrt(N).

Those two predictions differ and the test separates them:

    mean bias vs domain size      flat  -> bulk property
                                  falls -> sample-size artifact
    scatter vs domain size        falls -> the domain is averaging, as it must

Three domain sizes at L/R_ave = 10, 20, 40 (26, 83, 312 grains), three seeds
each. The bias is taken from a two-point eps ladder: with k(eps) linear in eps
-- which the five-point ladder in bias.csv establishes -- the eps -> 0 intercept
is k0 = 2*k(eps1) - k(eps2) for eps2 = 2*eps1, and bias = k(eps1)/k0 - 1.

Raster is chosen per size to hold h fixed in grain radii where affordable, and
the band/h achieved is reported. k_eff is insensitive to this (0.3% over a 2.7x
raster change, measured), which is why it can be compared across sizes at all.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
PROJ = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(PROJ / "preprocess"))
from cell_solve import k_eff                                # noqa: E402
from measure_bias import load, phi_field, K_ICE, K_AIR      # noqa: E402
import packing_lib as pl                                    # noqa: E402

EPS1, EPS2 = 45e-9, 90e-9
RASTER = {10: 512, 20: 1024, 40: 1024}


def bias_of(d: Path, raster: int) -> dict:
    m, cen, rad = load(d)
    Lx, Ly = m["Lx"], m["Ly"]
    h = Lx / raster
    k = {}
    for e in (EPS1, EPS2):
        phi = phi_field(cen, rad, Lx, Ly, e, raster)
        ke = k_eff(K_AIR + (K_ICE - K_AIR) * phi, h)
        k[e] = 0.5 * (ke[0, 0] + ke[1, 1])
        if e == EPS1:
            anis = ke[1, 1] / ke[0, 0]
    k0 = 2 * k[EPS1] - k[EPS2]                 # linear extrapolation to eps = 0
    # how much of the contact count is spurious, for context
    b = pl.delaunay_bonds(cen, Lx, Ly, True, True)
    i, j = b[:, 0], b[:, 1]
    dd = cen[j] + b[:, 2:4] * np.array([Lx, Ly]) - cen[i]
    gap = np.hypot(dd[:, 0], dd[:, 1]) - (rad[i] + rad[j])
    tol = 0.02 * rad.min()
    return {
        "name": d.name, "n_grains": len(cen), "raster": raster,
        "band_per_h": 9.2 * EPS1 / h,
        "k_at_eps1": k[EPS1], "k_sharp": k0,
        "bias": k[EPS1] / k0 - 1.0, "k_anis": anis,
        "frac_contacts_spurious": float(np.sum((gap > tol) & (gap < 9.2 * EPS1))
                                        / max(np.sum(gap < 9.2 * EPS1), 1)),
    }


def main() -> int:
    rows = []
    print(f"{'packing':22s} {'n':>5} {'band/h':>7} {'k(45nm)':>9} {'k(0)':>8} "
          f"{'bias':>7} {'anis':>6} {'spurious':>9}")
    for lr in (10, 20, 40):
        for s in (1, 2, 3):
            d = HERE / f"packings/LR{lr}_phi0.325_seed{s}"
            r = bias_of(d, RASTER[lr]); r["L_over_R"] = lr
            rows.append(r)
            print(f"{r['name']:22s} {r['n_grains']:5d} {r['band_per_h']:7.2f} "
                  f"{r['k_at_eps1']:9.4f} {r['k_sharp']:8.4f} {r['bias']:+7.1%} "
                  f"{r['k_anis']:6.3f} {r['frac_contacts_spurious']:9.1%}", flush=True)

    with (HERE / "rev_bias.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    print(f"\n{'L/R_ave':>8} {'grains':>7} {'mean bias':>11} {'sd(bias)':>10} "
          f"{'sd(k at 45nm)/mean':>19}")
    for lr in (10, 20, 40):
        g = [r for r in rows if r["L_over_R"] == lr]
        b = np.array([r["bias"] for r in g]); kk = np.array([r["k_at_eps1"] for r in g])
        print(f"{lr:8d} {g[0]['n_grains']:7d} {b.mean():+11.1%} {b.std(ddof=1):10.1%} "
              f"{kk.std(ddof=1)/kk.mean():19.1%}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
