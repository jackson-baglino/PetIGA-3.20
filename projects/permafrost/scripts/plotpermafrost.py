#!/usr/bin/env python3
import os
import glob
from igakit.io import PetIGA, VTK
from numpy import linspace

def convert_solution_files(nrb_path, pattern, scalar_map, vtk_prefix, vtk_dir="vtkOut"):
    """
    Convert solution .dat files to VTK format if not already present.

    Parameters:
        nrb_path (str): Path to the IGA geometry file.
        pattern (str): Glob pattern to match input solution files.
        scalar_map (dict): Mapping of scalar field names to their indices.
        vtk_prefix (str): Prefix for VTK output filenames.
        vtk_dir (str): Directory to store output VTK files.
    """
    print(f"Processing files matching pattern: {pattern}")
    
    # Read IGA geometry
    try:
        nrb = PetIGA().read(nrb_path)
    except Exception as e:
        print(f"❌ Failed to read NURBS from {nrb_path}: {e}")
        return

    # Ensure output directory exists
    os.makedirs(vtk_dir, exist_ok=True)

    for infile in glob.glob(pattern):
        name = os.path.splitext(os.path.basename(infile))[0]
        try:
            number = name.split("l")[1]
        except IndexError:
            print(f"⚠️ Skipping file {infile}: unable to extract number from name")
            continue

        outfile = os.path.join(vtk_dir, f"{vtk_prefix}{number}.vtk")
        if not os.path.isfile(outfile):
            try:
                sol = PetIGA().read_vec(infile, nrb)
                VTK().write(outfile, nrb, fields=sol, scalars=scalar_map)
                print(f"✅ Wrote: {outfile}")
            except Exception as e:
                print(f"❌ Error processing {infile}: {e}")
        else:
            print(f"⏭️ Skipping existing file: {outfile}")

def main():
    convert_solution_files(
        nrb_path="igasol.dat",
        pattern="sol*.dat",
        scalar_map={'IcePhase': 0, 'Temperature': 1, 'VaporDensity': 2},
        vtk_prefix="solV"
    )

    convert_solution_files(
        nrb_path="igasoil.dat",
        pattern="soil.dat",
        scalar_map={'Sediment': 0},
        vtk_prefix="soil"
    )

if __name__ == "__main__":
    main()