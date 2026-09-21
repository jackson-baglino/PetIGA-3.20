#!/usr/bin/env python3
"""Correct the k_eff eps bias from saved snapshots -- no new simulations.

    venv_enceladus/bin/python studies/keff_sintering/eps_correction.py <run_dir>

WHY A CORRECTION IS NEEDED AT ALL. eps is a numerical parameter; the physical
answer is the eps -> 0 limit. k_eff is biased HIGH and linearly in eps (see
studies/packing_design/), which would be harmless for a study reporting only
ratios IF the bias were constant in time. It is not: it scales with eps*SSA,
and SSA falls ~37% over a run, so the bias shrinks and does NOT cancel in
k(t)/k(0) -- the one number this campaign reports.

WHY THIS NEEDS NO NEW RUNS. The bias is a property of the FIELD at each
instant, not of how the field got there. So it can be measured on snapshots
that already exist, by rebuilding the same geometry at other interface widths
and extrapolating.

    phi = sigma(s/eps)  =>  s = eps * logit(phi)

recovers the signed distance, and re-evaluating sigma(s/eps') re-widens the
same interface to eps'. Exact for the equilibrium logistic profile, which is
what the solver relaxes to.

WIDEN, NEVER SHARPEN. The .vts output grid is coarser than the solve grid
(945 vs 2829 here), so at the production eps the band is already only ~4 cells
across. Sharpening would push it below the grid and measure pixellation.
Widening is always resolvable, so the ladder runs UPWARD from the production
eps and the intercept is an extrapolation -- stated as such.

WHAT THIS DOES NOT COVER. Two different eps effects exist:

  (a) the k_eff MEASUREMENT bias on a given geometry -- this, removable here;
  (b) whether the EVOLUTION itself would differ at smaller eps -- not
      removable here, since the geometry was evolved at one eps. That needs a
      rerun, and (a) should be settled first because it is free and may be the
      whole story.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
PROJ = HERE.parents[1]
sys.path.insert(0, str(PROJ / "postprocess"))
sys.path.insert(0, str(PROJ / "studies/packing_design"))
import pplib                                                # noqa: E402
from cell_solve import k_eff                                # noqa: E402

K_ICE, K_AIR = 2.29, 0.02
DAY = 86400.0


def rewiden(phi: np.ndarray, eps0: float, eps1: float, clip: float = 1e-6):
    """Same interface, width eps1 instead of eps0, via the logistic profile."""
    p = np.clip(phi, clip, 1.0 - clip)
    s = eps0 * np.log(p / (1.0 - p))          # signed distance
    return 1.0 / (1.0 + np.exp(-s / eps1))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run", type=Path)
    ap.add_argument("--eps0", type=float, default=1.0e-6)
    ap.add_argument("--factors", type=float, nargs="+", default=[1.0, 1.5, 2.0, 3.0])
    ap.add_argument("--L", type=float, default=2.0e-3)
    ap.add_argument("--times", type=float, nargs="+", default=[0.0, 30.0],
                    help="days at which to correct")
    a = ap.parse_args()

    vts = sorted((a.run / "vtkOut").glob("solV_*.vts"),
                 key=lambda p: int(re.search(r"solV_(\d+)", p.name).group(1)))
    if not vts:
        raise SystemExit(f"no vtkOut/*.vts under {a.run}")
    ssa = pplib.load_ssa(str(a.run))
    step2t = dict(zip(ssa[:, 3].astype(int), ssa[:, 2]))
    steps = np.array([int(re.search(r"solV_(\d+)", p.name).group(1)) for p in vts])
    times = np.array([step2t.get(s, np.nan) for s in steps])

    out = {}
    for tday in a.times:
        i = int(np.nanargmin(np.abs(times - tday * DAY)))
        fl, _, _ = pplib.read_vts(vts[i], want=["IcePhase"])
        phi0 = fl["IcePhase"]
        h = a.L / phi0.shape[1]
        print(f"\nt = {times[i]/DAY:.2f} d   (step {steps[i]}, grid {phi0.shape[0]})")
        print(f"  {'eps[um]':>8} {'band/h':>7} {'k_iso':>9}")
        es, ks = [], []
        for f in a.factors:
            e = a.eps0 * f
            phi = phi0 if f == 1.0 else rewiden(phi0, a.eps0, e)
            ke = k_eff(K_AIR + (K_ICE - K_AIR) * phi, h)
            kiso = 0.5 * (ke[0, 0] + ke[1, 1])
            es.append(e); ks.append(kiso)
            print(f"  {e*1e6:8.2f} {9.2*e/h:7.2f} {kiso:9.4f}")
        es, ks = np.array(es), np.array(ks)
        sl, ic = np.polyfit(es, ks, 1)
        out[tday] = (ks[0], ic)
        print(f"  linear fit -> eps=0 intercept {ic:.4f}   "
              f"(measured at eps0: {ks[0]:.4f}, bias {ks[0]/ic-1:+.1%})")

    if len(out) >= 2:
        ta, tb = sorted(out)[0], sorted(out)[-1]
        m0, c0 = out[ta]
        m1, c1 = out[tb]
        print(f"\n--- the number the campaign reports ---")
        print(f"  measured  k({tb:.0f}d)/k({ta:.0f}d) = {m1/m0:.4f}  ({m1/m0-1:+.1%})")
        print(f"  corrected k({tb:.0f}d)/k({ta:.0f}d) = {c1/c0:.4f}  ({c1/c0-1:+.1%})")
        print(f"  the eps bias changes the headline by "
              f"{(c1/c0-1)/(m1/m0-1)-1:+.0%} of itself")
        print("\n  Intercepts are extrapolations from a ladder that only widens")
        print("  (the output grid cannot resolve a sharper band), so quote them")
        print("  with the ladder, not alone.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
