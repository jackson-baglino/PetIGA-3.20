#!/usr/bin/env python3
"""
plotpermafrost.py  —  Convert PetIGA sol_*.dat output to VTK files.

Exports four DOFs plus one derived field:
  IcePhase     — φ_i          (DOF 0)
  Temperature  — T            (DOF 1)
  VaporDensity — ρ_v          (DOF 2)
  SedPhase     — φ_s          (DOF 3)
  AirPhase     — 1 − φ_i − φ_s  (derived, clipped to [0, 1])

VTK files are written to ./vtkOut/ and skipped if already present.

Usage
-----
  python plotpermafrost.py                 # run from the output directory
  python plotpermafrost.py --dir /path/to/run
"""

import argparse
import glob
import os

import numpy as np
from igakit.io import PetIGA, VTK


def convert(run_dir: str = ".", iga_file: str = "igasol.dat",
            force: bool = False):
    iga_path = os.path.join(run_dir, iga_file)
    if not os.path.isfile(iga_path):
        raise FileNotFoundError(f"IGA geometry file not found: {iga_path}")

    out_dir = os.path.join(run_dir, "vtkOut")
    os.makedirs(out_dir, exist_ok=True)

    nrb = PetIGA().read(iga_path)

    sol_files = sorted(glob.glob(os.path.join(run_dir, "sol*.dat")))
    if not sol_files:
        print(f"No sol*.dat files found in '{run_dir}'")
        return

    for infile in sol_files:
        name   = os.path.splitext(os.path.basename(infile))[0]
        number = name.split("l")[1]
        outfile = os.path.join(out_dir, f"solV{number}.vtk")

        if not force and os.path.isfile(outfile):
            print(f"  Skipping (exists): {outfile}")
            continue

        sol = PetIGA().read_vec(infile, nrb)

        scalars = {
            "IcePhase":    0,
            "Temperature": 1,
            "VaporDensity": 2,
            "SedPhase":    3,
        }
        # Only export DOFs present in the solution vector
        if sol.ndim < 2 or sol.shape[-1] < 4:
            scalars.pop("SedPhase")

        # Append AirPhase = 1 - φ_i - φ_s as a derived field
        if sol.ndim >= 2 and sol.shape[-1] >= 4:
            air = np.clip(1.0 - sol[..., 0:1] - sol[..., 3:4], 0.0, 1.0)
            sol = np.concatenate([sol, air], axis=-1)
            scalars["AirPhase"] = sol.shape[-1] - 1

        VTK().write(outfile, nrb, fields=sol, scalars=scalars)
        print(f"  Written: {outfile}")


def parse_args():
    p = argparse.ArgumentParser(description="Convert PetIGA output to VTK.")
    p.add_argument("--dir",   default=".",          help="Run directory (default: .)")
    p.add_argument("--iga",   default="igasol.dat", help="IGA geometry file")
    p.add_argument("--force", action="store_true",  help="Overwrite existing VTK files")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    convert(run_dir=args.dir, iga_file=args.iga, force=args.force)
