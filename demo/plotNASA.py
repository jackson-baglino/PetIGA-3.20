import numpy as np
import glob
import os
from igakit.io import PetIGA, VTK
from numpy import linspace
from skimage import measure
import meshio
import pyvista as pv

from rich.traceback import install
install()

def create_vtk_files():
    """
    Creates VTK files for ice and sediment grains and converts them to STL files.

    Reads input data files and generates VTK files for visualization.
    The VTK files contain information about the ice and sediment grains,
    including ice phase, temperature, vapor density, and sediment phase.
    """

    # Get the values of Nx, Ny, Nz from environment variables
    maxNum = 300
    numX = int(os.getenv("Nx", maxNum))
    numY = int(os.getenv("Ny", maxNum))
    numZ = int(os.getenv("Nz", maxNum))

    num_points = max(numX, numY, numZ)

    # Print the number of points
    print(f"Number of points: {num_points}.\n")

    # Read the input data file for ice grains
    nrb = PetIGA().read("igasol.dat")

    # Define a uniform sampling function
    uniform = lambda U: linspace(U[0], U[-1], num_points)

    # Import ice grains:
    for infile in glob.glob("sol*.dat"):
        name = infile.split(".")[0]
        number = name.split("l")[1]
        vtk_root = f'./vtkOut/solV_{number}.vtk'
        stl_root = f'./stlOut/IcePhase_{number}.stl'

        # Read the solution vector from the input file
        sol = PetIGA().read_vec(infile, nrb)

        # Remove the existing VTK file if it exists
        if os.path.exists(vtk_root):
            os.remove(vtk_root)

        # Write the VTK file with ice grain information
        VTK().write(vtk_root,
                    nrb, fields=sol, 
                    sampler=uniform, 
                    scalars={'IcePhase': 0, 'Temperature': 1, 'VaporDensity': 2})
        print(f"Created: {vtk_root}.\n")

        # Convert the VTK file to STL
        convert_vtk_to_stl(vtk_root, stl_root)
        print(f"Converted: {vtk_root} to {stl_root}.\n")

def convert_vtk_to_stl(vtk_file, stl_file):
    """
    Converts a VTK file to STL file.

    Parameters:
    vtk_file (str): Path to the input VTK file.
    stl_file (str): Path to the output STL file.
    """
    # Read the VTK file
    mesh = pv.read(vtk_file)

    # Extract the 'IcePhase' field
    if 'IcePhase' in mesh.array_names:
        ice_phase = mesh['IcePhase']

        # Add 'IcePhase' as the active scalars
        mesh.point_data['scalars'] = ice_phase
        mesh.set_active_scalars('scalars')

        # Apply threshold to the 'IcePhase' field
        threshold = 0.5
        thresholded_mesh = mesh.threshold(threshold, scalars='scalars')

        # Convert UnstructuredGrid to PolyData for smoothing
        poly_data = thresholded_mesh.extract_surface()

        # Smooth the mesh
        smoothed_mesh = poly_data.smooth(n_iter=50, relaxation_factor=0.1)

        # Save the smoothed mesh to STL
        print(f"Writing: {stl_file}....................", end="")
        smoothed_mesh.save(stl_file)
        print(" Complete!\n")

    else:
        print(f"Error: 'IcePhase' field not found in {vtk_file}")

# Call the function to create VTK files
def main():

    create_vtk_files()

if __name__ == "__main__":
    main()