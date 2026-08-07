#!/usr/bin/env python3
"""
estimate_memory.py — peak RAM for a run, from its geometry .opts.

Written after a local smoke test of calib_2D_L500um_eps0.50um_pack_s025 (24 M DOF) exhausted a
64 GB machine. Mesh size alone is a bad proxy: the Jacobian and especially the
ILU(3) fill dominate, and solver.opts uses -sub_pc_factor_levels 3.

Model, for 2D tensor-product B-splines of degree p with dof unknowns per node:
  nonzeros/row = (2p+1)^2 * dof        each basis fn overlaps (2p+1)^2 others
  Jacobian     = DOF * nnz/row * (sizeof(PetscScalar) + sizeof(PetscInt))
  ILU(k) fill  ~ k * Jacobian          (empirical, level-dependent)
  vectors      ~ 12 * DOF * 8 B        TS/SNES/KSP work vectors

Deliberately an over-estimate: being wrong high costs a queued job, being wrong
low costs the workstation.

  ./preprocess/estimate_memory.py inputs/geometry/calib_pack_s050.opts
  ./preprocess/estimate_memory.py --budget-gb 8 inputs/geometry/calib_*.opts
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


def opt(text: str, key: str, default=None):
    m = re.search(rf"^{re.escape(key)}\s+(\S+)", text, re.MULTILINE)
    return m.group(1) if m else default


def estimate(geom: Path, solver: Path):
    g = geom.read_text()
    s = solver.read_text() if solver.is_file() else ""
    p = int(opt(s, "-p", 2))
    dof = int(opt(s, "-dof", 3))
    lev = int(opt(s, "-sub_pc_factor_levels", 3))
    dim = int(opt(g, "-dim", 2))
    nx = int(opt(g, "-Nx", 1))
    ny = int(opt(g, "-Ny", 1))
    nz = int(opt(g, "-Nz", 1)) if dim == 3 else 1

    ndof = dof * nx * ny * nz
    nnz_row = ((2 * p + 1) ** dim) * dof
    jac = ndof * nnz_row * (8 + 4)
    ilu = jac * lev
    vecs = 12 * ndof * 8
    return dict(nx=nx, ny=ny, nz=nz, dof=ndof, p=p, lev=lev,
                jac=jac, ilu=ilu, vecs=vecs, total=jac + ilu + vecs)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("geometry", nargs="+", type=Path)
    ap.add_argument("--solver", type=Path,
                    default=Path(__file__).resolve().parent.parent / "inputs/solver.opts")
    ap.add_argument("--budget-gb", type=float, default=None,
                    help="flag anything above this as too big for the machine")
    args = ap.parse_args(argv)

    G = 2 ** 30
    print(f"{'geometry':>30} {'Nx':>6} {'Ny':>6} {'DOF':>12} "
          f"{'Jac':>8} {'ILU':>8} {'TOTAL':>9}"
          + ("  verdict" if args.budget_gb else ""))
    over = 0
    for f in args.geometry:
        e = estimate(f, args.solver)
        line = (f"{f.stem:>30} {e['nx']:>6} {e['ny']:>6} {e['dof']:>12,} "
                f"{e['jac']/G:>7.1f}G {e['ilu']/G:>7.1f}G {e['total']/G:>8.1f}G")
        if args.budget_gb:
            ok = e["total"] / G <= args.budget_gb
            line += "  local" if ok else "  HPC"
            over += 0 if ok else 1
        print(line)
    if args.budget_gb:
        print(f"\n  budget {args.budget_gb:g} GB -> {over} of {len(args.geometry)} "
              f"need the cluster")
    return 0


if __name__ == "__main__":
    sys.exit(main())
