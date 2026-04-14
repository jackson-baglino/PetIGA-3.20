#!/usr/bin/env python3
"""
makevtk.py

Convert PetIGA solution .dat files to VTK format.

The sediment field (soil.dat / igasoil.dat) is read once and merged into
every sol*.vtk output.  AirPhase = 1 - IcePhase - Sediment is also computed
and embedded, so all four fields are available in ParaView on the time series
with no additional filters needed.

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
# Sediment loader
# ---------------------------------------------------------------------------

def load_sediment(soil_dat: str, igasoil_dat: str):
    """
    Read the static sediment field once.
    Returns (nrb_soil, sol_soil), or (None, None) if either file is missing.
    """
    if not os.path.isfile(igasoil_dat):
        print(f"  ⚠️  Sediment geometry not found: {igasoil_dat}  (Sediment will NOT be embedded)")
        return None, None
    if not os.path.isfile(soil_dat):
        print(f"  ⚠️  Sediment solution not found: {soil_dat}  (Sediment will NOT be embedded)")
        return None, None

    print(f"\n--- Loading sediment field ---")
    nrb_soil = read_nurbs(igasoil_dat)
    sol_soil  = read_solution_vec(soil_dat, nrb_soil)
    print(f"  ✅ Sediment field loaded from '{soil_dat}'")
    print(f"  ℹ️  sol_soil shape: {sol_soil.shape}")
    return nrb_soil, sol_soil


# ---------------------------------------------------------------------------
# Core conversion
# ---------------------------------------------------------------------------

def convert_solution_files(
    nrb_path:      str,
    pattern:       str,
    scalar_map:    dict,
    vtk_prefix:    str,
    vtk_dir:       str        = "vtkOut",
    force:         bool       = False,
    sediment_sol:  np.ndarray = None,
):
    """
    Convert solution .dat files matched by *pattern* to VTK.

    If sediment_sol is provided it is appended to every output file as the
    'Sediment' field, and AirPhase = 1 - IcePhase - Sediment is also written.
    IcePhase is assumed to be component 0 of the primary solution (scalar_map).

    Parameters
    ----------
    nrb_path     : path to the IGA geometry file for this solution set
    pattern      : glob pattern for input .dat files
    scalar_map   : {field_name: component_index} for the primary solution
    vtk_prefix   : prefix for output VTK filenames
    vtk_dir      : output directory
    force        : if True, overwrite existing VTK files
    sediment_sol : solution array from soil.dat (same mesh, may be 2D or 3D)
    """
    print(f"\n--- Converting '{pattern}' ---")

    infiles = sorted(glob.glob(pattern))
    if not infiles:
        print(f"  ⚠️  No files matched pattern '{pattern}' — nothing to do.")
        return

    nrb = read_nurbs(nrb_path)
    os.makedirs(vtk_dir, exist_ok=True)

    # Normalise sediment to (nx, ny, 1) once, outside the loop
    sed_norm = None
    if sediment_sol is not None:
        if sediment_sol.ndim == 2:
            # Single-field solution returned as (nx, ny) — add component axis
            sed_norm = sediment_sol[..., np.newaxis]
        elif sediment_sol.ndim == 3:
            # Multi-field solution — Sediment is component 0
            sed_norm = sediment_sol[..., 0:1]
        else:
            print(f"  ⚠️  Unexpected sediment_sol shape {sediment_sol.shape} — skipping embed.")

        if sed_norm is not None:
            print(f"  ℹ️  Normalised sediment shape: {sed_norm.shape}")

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

            if sed_norm is not None:
                # Guarantee primary sol is 3D: (nx, ny, ncomp)
                if sol.ndim == 2:
                    sol = sol[..., np.newaxis]

                # print(f"    sol shape: {sol.shape}  |  sed shape: {sed_norm.shape}")

                # IcePhase is component 0 of the primary solution
                ice = sol[..., 0:1]         # (nx, ny, 1)
                air = 1.0 - ice - sed_norm  # (nx, ny, 1)

                ncomp = sol.shape[-1]
                combined_scalars["Sediment"] = ncomp
                combined_scalars["AirPhase"] = ncomp + 1

                sol = np.concatenate([sol, sed_norm, air], axis=-1)
                # print(f"    combined sol shape: {sol.shape}")

            VTK().write(outfile, nrb, fields=sol, scalars=combined_scalars)
            print(f"  ✅ Wrote: {outfile}")

        except Exception as e:
            print(f"  ❌ Error processing '{infile}': {e}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert PetIGA .dat files to VTK, embedding Sediment and AirPhase."
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

    # ------------------------------------------------------------------
    # 1. Load the static sediment field (once)
    # ------------------------------------------------------------------
    _nrb_soil, sol_soil = load_sediment(
        soil_dat    = "soil.dat",
        igasoil_dat = "igasoil.dat",
    )

    # ------------------------------------------------------------------
    # 2. Convert time-series solution files, embedding Sediment + AirPhase
    # ------------------------------------------------------------------
    convert_solution_files(
        nrb_path     = "igasol.dat",
        pattern      = args.sol_pattern,
        scalar_map   = {"IcePhase": 0, "Temperature": 1, "VaporDensity": 2},
        vtk_prefix   = "solV",
        vtk_dir      = args.vtk_dir,
        force        = args.force,
        sediment_sol = sol_soil,
    )

    # ------------------------------------------------------------------
    # 3. Write the standalone soil VTK for inspection
    # ------------------------------------------------------------------
    convert_solution_files(
        nrb_path   = "igasoil.dat",
        pattern    = "soil.dat",
        scalar_map = {"Sediment": 0},
        vtk_prefix = "soil",
        vtk_dir    = args.vtk_dir,
        force      = args.force,
    )


if __name__ == "__main__":
    main()