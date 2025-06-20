from igakit.io import PetIGA, VTK
import numpy as np
import os
import glob

def ensure_directory(path):
    if not os.path.exists(path):
        os.makedirs(path)
        print(f"[INFO] Created directory: {path}")

def convert_field_data(file_pattern, iga_file, output_dir, field_name, scalar_map):
    # Read the geometry only once
    try:
        nrb = PetIGA().read(iga_file)
    except Exception as e:
        print(f"[ERROR] Failed to read {iga_file}: {e}")
        return

    files = sorted(glob.glob(file_pattern))
    if not files:
        print(f"[WARNING] No files found matching pattern: {file_pattern}")
        return

    for infile in files:
        basename = os.path.splitext(os.path.basename(infile))[0]
        index = basename.lstrip("sol")  # Extract numerical index from "sol00001"
        vtk_filename = f"{field_name}{index}.vtk"
        vtk_path = os.path.join(output_dir, vtk_filename)

        if os.path.isfile(vtk_path):
            print(f"[SKIP] {vtk_path} already exists.")
            continue

        try:
            sol = PetIGA().read_vec(infile, nrb)
            VTK().write(vtk_path, nrb, fields=sol, scalars=scalar_map)
            print(f"[WRITE] {vtk_path}")
        except Exception as e:
            print(f"[ERROR] Failed to process {infile}: {e}")

def main():
    output_dir = "vtkOut"
    ensure_directory(output_dir)

    # Convert volumetric simulation results
    convert_field_data(
        file_pattern="sol*.dat",
        iga_file="igasol.dat",
        output_dir=output_dir,
        field_name="solV",
        scalar_map={'IcePhase': 0, 'Temperature': 1, 'VaporDensity': 2}
    )

    # Convert static soil data
    if os.path.exists("soil.dat") and os.path.exists("igasoil.dat"):
        convert_field_data(
            file_pattern="soil.dat",
            iga_file="igasoil.dat",
            output_dir=output_dir,
            field_name="soil",
            scalar_map={'Sediment': 0}
        )
    else:
        print("[INFO] Skipping soil data: files not found")

if __name__ == "__main__":
    main()