#!/usr/bin/env python3
"""Write the HPC test spec for the disk-array ice-fraction sweep.

The disk ladder (../disk/) tests both conductivity laws at ONE ice fraction,
f = 0.196. This sweep repeats it over six disk radii, so each interface width
gives a curve k_eff(f) to put against Rayleigh's sharp curve. It shows two
things a single point cannot: the arithmetic law's bias at the production
resolution, and that the bias is not a constant fraction but grows with f.
That second point is why the bias survives the ratio k(t)/k(0) in a sintering
run, where the microstructure changes.

The width is held at a fixed fraction of the RADIUS, eps/R, because that is how
the production runs are resolved: eps = 1 um on R_ave = 50 um is eps/R = 0.02.
Mesh resolution scales with eps (eps/h = 5.12, as in both ladders), so the
interface error is not confounded with mesh convergence.

    venv_enceladus/bin/python studies/keff_sharp_limit/disk_sweep/make_sweep_spec.py

writes sweep_tests.txt (for submit_batch.sh --tests-file) and
sweep_manifest.csv (what each label means, read by collect_sweep.py).
"""

from __future__ import annotations

import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent

GEOM = "singleice_2D_L1mm_R250um_keff"   # -RCice is overridden per job
EXP = "snow_T-20_h1.00_1d"               # -keff_only exits before integrating
L = 1.0e-3
EPS_PER_ELEM = 5.12

RADII_UM = (125, 175, 250, 300, 350, 400)   # f = 0.049 ... 0.503
EPS_OVER_R = (0.04, 0.02, 0.01)             # 0.02 is the production resolution
# The solver's flag value for the old law is "arith"; everything we write
# says "arithmetic".
LAWS = {"arithmetic": "arith", "tensor": "tensor"}


def mesh(eps: float) -> int:
    n = round(EPS_PER_ELEM * L / eps)
    return n + (n % 2)                      # even keeps the centre on a mesh line


def main() -> None:
    rows = []
    for R_um in RADII_UM:
        R = R_um * 1e-6
        for ratio in EPS_OVER_R:
            eps = ratio * R
            N = mesh(eps)
            for law, flag in LAWS.items():
                label = f"{law}_R{R_um}_e{ratio:g}"
                opts = (f"--label {label} -keff 1 -keff_only 1 -keff_interp {flag} "
                        f"-eps {eps:.8e} -Nx {N} -Ny {N} -RCice {R:.6e} "
                        f"-keff_ksp_type cg -keff_pc_type gamg")
                rows.append(dict(label=label, law=law, R=R, eps_over_R=ratio,
                                 eps=eps, N=N, spec=f"{GEOM}:{EXP}:{opts}"))

    (HERE / "sweep_tests.txt").write_text("".join(r["spec"] + "\n" for r in rows))
    with (HERE / "sweep_manifest.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["label", "law", "R", "eps_over_R", "eps", "N"],
                           lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in w.fieldnames})

    print(f"{len(rows)} jobs -> sweep_tests.txt, sweep_manifest.csv")
    for r in rows[::2]:
        print(f"  R = {r['R']*1e6:5.0f} um  eps/R = {r['eps_over_R']:<5g} "
              f"eps = {r['eps']*1e6:6.3f} um  N = {r['N']}")


if __name__ == "__main__":
    main()
