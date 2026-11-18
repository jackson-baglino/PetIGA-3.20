"""Compare two runs of the same case: A1 against the derived model (or any two runs on the same mesh and ranks).

    venv_lunar/bin/python studies/derived_model/verification/compare_runs.py <run_A1> <run_derived>

Prints, from SSA_evo.dat (monitoring.c), the largest relative difference over time of the interface measure
(sub_interf/eps), the ice volume, the vapor content and the total water mass, with the derived run
interpolated onto the A1 times; and, from the last sol_*.dat of each run (PETSc binary Vec, 3 DOFs
interleaved: phi, T, rho_v), the time of that snapshot and the largest pointwise difference of each field.
The two runs must use the same mesh and number of ranks, so the Vec orderings agree. Nothing is judged here.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "postprocess"))
import pplib  # noqa: E402

VEC_CLASSID = 1211214


def read_vec(path):
    with open(path, "rb") as f:
        head = np.fromfile(f, dtype=">i4", count=2)
        if head[0] != VEC_CLASSID:
            raise SystemExit(f"{path}: not a PETSc Vec (class id {head[0]})")
        return np.fromfile(f, dtype=">f8", count=int(head[1]))


def last_snapshot(run):
    files = sorted(Path(run).glob("sol_*.dat"))
    if not files:
        raise SystemExit(f"{run}: no sol_*.dat")
    f = files[-1]
    step = int(f.stem.split("_")[1])
    ssa = pplib.load_ssa(str(run))
    rows = ssa[ssa[:, pplib.STEP] == step]
    t = float(rows[0, pplib.TIME]) if len(rows) else float("nan")
    return f, step, t


def main(run_a, run_d):
    a, d = pplib.load_ssa(run_a), pplib.load_ssa(run_d)
    print(f"A1      : {run_a}\n          {len(a)} steps, t_end = {a[-1, pplib.TIME]:.6e} s")
    print(f"derived : {run_d}\n          {len(d)} steps, t_end = {d[-1, pplib.TIME]:.6e} s")
    t = a[:, pplib.TIME]
    keep = t <= d[-1, pplib.TIME]
    for name, col in (("interface measure", pplib.SSA), ("ice volume", pplib.ICE),
                      ("vapor content", pplib.RHOV), ("water mass", pplib.MASS)):
        ya = a[keep, col]
        yd = np.interp(t[keep], d[:, pplib.TIME], d[:, col])
        scale = np.max(np.abs(ya))
        rel = np.max(np.abs(yd - ya)) / scale if scale > 0 else float("nan")
        print(f"  {name:18s} max |derived - A1| / max |A1| = {rel:.3e}   "
              f"(A1 start {ya[0]:.6e}, end {ya[-1]:.6e}; derived end {yd[-1]:.6e})")
    fa, sa, ta = last_snapshot(run_a)
    fd, sd, td = last_snapshot(run_d)
    va, vd = read_vec(fa), read_vec(fd)
    print(f"  last snapshot: A1 {fa.name} (t = {ta:.6e} s), derived {fd.name} (t = {td:.6e} s)")
    if va.size != vd.size:
        raise SystemExit("snapshots differ in size: not the same mesh")
    for k, name in enumerate(("phi", "T", "rho_v")):
        xa, xd = va[k::3], vd[k::3]
        print(f"  {name:6s} max |derived - A1| = {np.max(np.abs(xd - xa)):.3e}   (A1 range {xa.min():.6e} .. {xa.max():.6e})")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    main(sys.argv[1], sys.argv[2])
