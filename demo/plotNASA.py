# import numpy as np
# import glob
# import os
# from igakit.io import PetIGA, VTK
# from numpy import linspace
# from skimage import measure
# import pyvista as pv
# import vtk

# from rich.traceback import install
# install()

# def create_vtk_files():
#     maxNum = 300
#     numX = int(os.getenv("Nx", maxNum))
#     numY = int(os.getenv("Ny", maxNum))
#     numZ = int(os.getenv("Nz", maxNum))

#     num_points = max(numX, numY, numZ)
#     print(f"Number of points: {num_points}.\n")

#     nrb = PetIGA().read("igasol.dat")
#     uniform = lambda U: linspace(U[0], U[-1], num_points)

#     for infile in glob.glob("sol*.dat"):
#         name = infile.split(".")[0]
#         number = name.split("l")[1]
#         vtk_root = f'./vtkOut/solV_{number}.vtk'
#         # stl_root = f'./stlOut/IcePhase_{number}.stl'

#         if not os.path.exists(vtk_root):
#             sol = PetIGA().read_vec(infile, nrb)
            
#             VTK().write(vtk_root,
#                         nrb, fields=sol, 
#                         sampler=uniform, 
#                         scalars={'IcePhase': 0, 'Temperature': 1, 'VaporDensity': 2})
#             print(f"Created: {vtk_root}.\n")

#             verify_vtk_fields(vtk_root)
#             # convert_vtk_to_stl(vtk_root, stl_root)
#             # print(f"Converted: {vtk_root} to {stl_root}.\n")
#         else:
#             print(f"File already exists: {vtk_root}")

# def verify_vtk_fields(vtk_file):
#     reader = vtk.vtkStructuredGridReader()
#     reader.SetFileName(vtk_file)
#     reader.Update()
#     vtk_data = reader.GetOutput()

#     cell_data = vtk_data.GetCellData()
#     point_data = vtk_data.GetPointData()
    
#     print(f"Fields in {vtk_file}:")
#     print("Cell Data Fields:")
#     for i in range(cell_data.GetNumberOfArrays()):
#         print(f"  {cell_data.GetArrayName(i)}")
        
#     print("Point Data Fields:")
#     for i in range(point_data.GetNumberOfArrays()):
#         print(f"  {point_data.GetArrayName(i)}")

# def apply_contour_filter(vtk_data, field_name, threshold_value, vtk_file):
#     if field_name not in [vtk_data.GetPointData().GetArrayName(i) for i in range(vtk_data.GetPointData().GetNumberOfArrays())]:
#         print(f"Field '{field_name}' not found in the VTK file.")
#         verify_vtk_fields(vtk_file)
#         return None

#     print("Applying contour filter...")
#     contour_filter = vtk.vtkContourFilter()
#     contour_filter.SetInputData(vtk_data)
#     contour_filter.SetValue(0, threshold_value)
#     contour_filter.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, field_name)
#     contour_filter.Update()
#     return contour_filter.GetOutput()

# def save_vtk_file(vtk_data, output_filename):
#     writer = vtk.vtkXMLPolyDataWriter()
#     writer.SetFileName(output_filename)
#     writer.SetInputData(vtk_data)
#     writer.Write()

# def process_files(input_folder, output_folder, field_name, threshold_value):
#     if not os.path.exists(output_folder):
#         os.makedirs(output_folder)
    
#     vtk_files = glob.glob(os.path.join(input_folder, '*.vtk'))
    
#     for vtk_file in vtk_files:
#         print(f'Processing {vtk_file}...')
#         reader = vtk.vtkStructuredGridReader()
#         reader.SetFileName(vtk_file)
#         reader.Update()
#         vtk_data = reader.GetOutput()

#         contour_data = apply_contour_filter(vtk_data, field_name, threshold_value, vtk_file)
        
#         if contour_data is not None:
#             output_filename = os.path.join(output_folder, os.path.basename(vtk_file).replace('.vtk', '_contour.vtp'))
#             save_vtk_file(contour_data, output_filename)
#             print(f'Saved contour to {output_filename}')
#         else:
#             print(f"Error processing {vtk_file}")

# if __name__ == '__main__':
#     input_folder = './vtkOut/'
#     output_folder = './stlOut/'
#     field_name = 'IcePhase'
#     threshold_value = 0.5
    
#     create_vtk_files()
#     process_files(input_folder, output_folder, field_name, threshold_value)


import numpy as np
import glob
import os
from igakit.io import PetIGA
from numpy import linspace
from rich.traceback import install

install()

# Define user-defined variables
tempFlag = 0
vaporFlag = 0

def create_dat_files():
    maxNum = 300
    numX = int(os.getenv("Nx", maxNum))
    numY = int(os.getenv("Ny", maxNum))
    numZ = int(os.getenv("Nz", maxNum))

    num_points = (numX, numY, numZ)
    print(f"Number of points: {num_points}.\n")

    nrb = PetIGA().read("igasol.dat")
    uniform = lambda U: linspace(U[0], U[-1], num_points[0])

    # Ensure directories exist
    if not os.path.exists('output/icePhase/'):
        os.makedirs('output/icePhase/')
    
    if tempFlag:
        if not os.path.exists('output/Temperature/'):
            os.makedirs('output/Temperature/')
    
    if vaporFlag:
        if not os.path.exists('output/VaporDensity/'):
            os.makedirs('output/VaporDensity/')

    for infile in glob.glob("sol*.dat"):
        name = infile.split(".")[0]
        number = name.split("l")[1]

        # Read the solution (assuming this is a 3D field)
        sol = PetIGA().read_vec(infile, nrb)

        # Extract fields (assume 3D array)
        icePhase = sol[..., 0]  # Assuming IcePhase is the first field

        if tempFlag:
            temperature = sol[..., 1]  # Assuming Temperature is the second field

        if vaporFlag:
            vaporDensity = sol[..., 2]  # Assuming VaporDensity is the third field

        # Save 3D fields to .dat files with coordinates
        save_field_to_dat_3D(icePhase, num_points, f'output/icePhase/icePhase_{number}.dat')

        if tempFlag:
            save_field_to_dat_3D(temperature, num_points, f'output/Temperature/Temperature_{number}.dat')
        if vaporFlag:
             save_field_to_dat_3D(vaporDensity, num_points, f'output/VaporDensity/VaporDensity_{number}.dat')

        print(f"Created .dat files for sol{number}.\n")

def save_field_to_dat_3D(field_data, num_points, output_filename):
    """ Save the 3D field data to a .dat file with coordinates (x, y, z, value). """
    
    numX, numY, numZ = num_points
    with open(output_filename, 'w') as f:
        for i in range(numX):
            for j in range(numY):
                for k in range(numZ):
                    # Save coordinates and corresponding value from the field
                    f.write(f"{i}, {j}, {k}, {field_data[i, j, k]}\n")
    
    print(f"Saved: {output_filename}.\n")

if __name__ == '__main__':
    create_dat_files()