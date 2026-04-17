#!/usr/bin/env python3
"""
makevtk.py

Convert PetIGA solution .dat files to VTK format.

Sediment is now DOF 3 in the solution vector (ice=0, temp=1, rhov=2, sed=3).
AirPhase = 1 - IcePhase - Sediment is computed and embedded automatically.

Usage:
    python makevtk.py [--sol-pattern SOL_PATTERN] [--vtk-dir VTK_DIR] [--force]

    --sol-pattern : glob pattern for solution files  (default: sol*.dat)
    --vtk-dir     : output directory                 (default: vtkOut)
    --force       : overwrite existing VTK files     (default: skip)
"""

import os
import glob
import argparse
import numpy as np
from igakit.io import PetIGA, VTK


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def read_nurbs(path: str):
    """Read an IGA geometry file and return the NURBS object."""
    try:
        nrb = PetIGA().read(path)
        print(f"  ✅ Read geometry: {path}")
        return nrb
    except Exception as e:
        raise RuntimeError(f"Failed to read NURBS from '{path}': {e}") from e


def read_solution_vec(dat_path: str, nrb):
    """Read a PetIGA solution vector and return the raw array."""
    try:
        sol = PetIGA().read_vec(dat_path, nrb)
        return sol
    except Exception as e:
        raise RuntimeError(f"Failed to read solution vector '{dat_path}': {e}") from e


def extract_number(filename: str) -> str:
    """
    Extract the numeric suffix from a filename of the form sol<N>.dat.
    e.g.  sol_00042.dat  ->  '_00042'
          sol1.dat       ->  '1'
    """
    name = os.path.splitext(os.path.basename(filename))[0]  # e.g. 'sol_00042'
    digits = name.lstrip("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ")
    if not digits:
        raise ValueError(f"Cannot extract numeric suffix from '{filename}'")
    return digits


# ---------------------------------------------------------------------------
# Core conversion
# ---------------------------------------------------------------------------

def convert_solution_files(
    nrb_path:   str,
    pattern:    str,
    scalar_map: dict,
    vtk_prefix: str,
    vtk_dir:    str  = "vtkOut",
    force:      bool = False,
):
    """
    Convert solution .dat files matched by *pattern* to VTK.

    AirPhase = 1 - IcePhase - Sediment is appended automatically when both
    IcePhase (comp 0) and Sediment (comp 3) are present in scalar_map.
    """
    print(f"\n--- Converting '{pattern}' ---")

    infiles = sorted(glob.glob(pattern))
    if not infiles:
        print(f"  ⚠️  No files matched pattern '{pattern}' — nothing to do.")
        return

    nrb = read_nurbs(nrb_path)
    os.makedirs(vtk_dir, exist_ok=True)

    add_air = ("IcePhase" in scalar_map and "Sediment" in scalar_map)

    for infile in infiles:
        try:
            number = extract_number(infile)
        except ValueError as e:
            print(f"  ⚠️  Skipping '{infile}': {e}")
            continue

        outfile = os.path.join(vtk_dir, f"{vtk_prefix}{number}.vtk")

        if os.path.isfile(outfile) and not force:
            print(f"  ⏭️  Skipping existing: {outfile}")
            continue

        try:
            sol = read_solution_vec(infile, nrb)

            combined_scalars = dict(scalar_map)

            if add_air:
                if sol.ndim == 2:
                    sol = sol[..., np.newaxis]
                ice = sol[..., 0:1]
                sed = sol[..., 3:4]
                air = 1.0 - ice - sed
                ncomp = sol.shape[-1]
                combined_scalars["AirPhase"] = ncomp
                sol = np.concatenate([sol, air], axis=-1)

            VTK().write(outfile, nrb, fields=sol, scalars=combined_scalars)
            print(f"  ✅ Wrote: {outfile}")

        except Exception as e:
            print(f"  ❌ Error processing '{infile}': {e}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert PetIGA .dat files to VTK, embedding AirPhase."
    )
    parser.add_argument(
        "--sol-pattern", default="sol*.dat",
        help="Glob pattern for time-series solution files (default: sol*.dat)"
    )
    parser.add_argument(
        "--vtk-dir", default="vtkOut",
        help="Output directory for VTK files (default: vtkOut)"
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Overwrite existing VTK files (default: skip)"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    convert_solution_files(
        nrb_path   = "igasol.dat",
        pattern    = args.sol_pattern,
        scalar_map = {
            "IcePhase":    0,
            "Temperature": 1,
            "VaporDensity": 2,
            "Sediment":    3,
        },
        vtk_prefix = "solV",
        vtk_dir    = args.vtk_dir,
        force      = args.force,
    )


if __name__ == "__main__":
    main()
