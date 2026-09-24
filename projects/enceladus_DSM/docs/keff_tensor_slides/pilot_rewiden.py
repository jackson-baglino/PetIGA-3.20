#!/usr/bin/env python3
"""k_eff of one pilot packing snapshot re-widened to larger eps, three ways.

Data for the slides' packing-level figures (slides 1 and 7). No simulation:
the stored phase field of pilot seed 1 is re-widened analytically,
phi = sigma(s/eps) -> s = eps logit(phi) -> sigma(s/eps'), exactly as in
studies/keff_sintering/eps_correction.py, and k_eff is computed with the numpy
finite-volume cell solver (studies/packing_design/cell_solve.py) under

  arithmetic  K = K_a + (K_i - K_a) phi                  the pre-09-23 law
  tensor      arithmetic along, harmonic across          (diagonal-only FV
              approximation of the tensor; see _faces_tensor)
  sharp       K on the thresholded geometry phi >= 1/2   independent of eps
              by construction; carries an O(h) staircase error instead

This is the 945^2 OUTPUT grid, not the 2829^2 solve grid, so the band is only
~4 cells across at the production eps, which is why the ladder only widens.

    venv_enceladus/bin/python docs/keff_tensor_slides/pilot_rewiden.py
writes pilot_rewiden.csv next to this script.
"""
from __future__ import annotations

import csv
import re
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PROJ = HERE.parents[1]
sys.path.insert(0, str(PROJ / "postprocess"))
sys.path.insert(0, str(PROJ / "studies/packing_design"))
sys.path.insert(0, str(PROJ / "studies/keff_sintering"))
import pplib                                                   # noqa: E402
from cell_solve import k_eff, k_eff_tensor, k_eff_sharp         # noqa: E402
from eps_correction import rewiden                              # noqa: E402

RUN = Path.home() / ("SimulationResults/HPC_results/enceladus_DSM/GrainPackingSintering/"
                     "batch_2026-09-16__13.16.11_pilot_keff/"
                     "packing_2D_pilot_phi0.325_Rave50um_LR40_seed1_L2mm_eps1000nm_perxy_T-20")
K_ICE, K_AIR = 2.29, 0.02
EPS0, L = 1.0e-6, 2.0e-3
FACTORS = (1.0, 1.5, 2.0, 3.0)
TIMES_D = (0.0, 1.0, 30.0)
DAY = 86400.0


def main() -> int:
    vts = sorted((RUN / "vtkOut").glob("solV_*.vts"),
                 key=lambda p: int(re.search(r"solV_(\d+)", p.name).group(1)))
    ssa = pplib.load_ssa(str(RUN))
    step2t = dict(zip(ssa[:, 3].astype(int), ssa[:, 2]))
    steps = [int(re.search(r"solV_(\d+)", p.name).group(1)) for p in vts]
    times = np.array([step2t.get(s, np.nan) for s in steps])

    rows = []
    for tday in TIMES_D:
        i = int(np.nanargmin(np.abs(times - tday * DAY)))
        fl, _, _ = pplib.read_vts(vts[i], want=["IcePhase"])
        phi0 = fl["IcePhase"]
        h = L / phi0.shape[1]
        t_act = times[i] / DAY
        ks = k_eff_sharp(phi0, h, K_AIR, K_ICE)
        k_sharp = 0.5 * (ks[0, 0] + ks[1, 1])
        for f in FACTORS:
            e = EPS0 * f
            phi = phi0 if f == 1.0 else rewiden(phi0, EPS0, e)
            t0 = time.time()
            ka = k_eff(K_AIR + (K_ICE - K_AIR) * phi, h)
            kt = k_eff_tensor(phi, h, K_AIR, K_ICE)
            row = dict(t_days=t_act, step=steps[i], eps=e, phi_bar=float(phi.mean()),
                       k_arithmetic=0.5 * (ka[0, 0] + ka[1, 1]),
                       k_tensor=0.5 * (kt[0, 0] + kt[1, 1]), k_sharp=k_sharp)
            rows.append(row)
            print(f"t={t_act:6.2f} d eps={e*1e6:.1f} um  arithmetic {row['k_arithmetic']:.4f}"
                  f"  tensor {row['k_tensor']:.4f}  sharp {k_sharp:.4f}"
                  f"  ({time.time()-t0:.0f} s)", flush=True)

    with (HERE / "pilot_rewiden.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
