from igakit.io import PetIGA,VTK
from numpy import linspace
import glob
import os

nrb = PetIGA().read("igasol.dat")
#uniform = lambda U: linspace(U[0], U[-1], 400)
for infile in glob.glob("sol*.dat"):
	name = infile.split(".")[0]
	number = name.split("l")[1]
	root = './solV'+number+'.vtk'
	if not os.path.isfile(root):
		sol = PetIGA().read_vec(infile,nrb)
		outfile = './vtkOut/solV'+number+'.vtk' 
		VTK().write(outfile,  
	            nrb,             
	            fields=sol,     
    		    scalars = {'IcePhase':0,'Temperature':1,'VaporDensity':2})







# ------------------------------------------------------------------------------
# ------------------------- CODE HERE STILL NEEDS WORK -------------------------
# ----------------------------- ONLY WORKS FOR 3D ------------------------------
# ------------------------------------------------------------------------------

# from igakit.io import PetIGA
# import glob
# import os
# import numpy as np
# from scipy.interpolate import RegularGridInterpolator

# def read_nrb(file_path):
#     """Reads the IGA solution file."""
#     return PetIGA().read(file_path)

# def interpolate_data_3d(data, x_dim, y_dim, z_dim):
#     """Interpolate 3D data using RegularGridInterpolator."""
#     # Original grid coordinates (assuming regular spacing)
#     x = np.linspace(0, 1, data.shape[0])
#     y = np.linspace(0, 1, data.shape[1])
#     z = np.linspace(0, 1, data.shape[2])

#     # Create an interpolator for the regular grid
#     interpolator = RegularGridInterpolator((x, y, z), data)

#     # New grid coordinates for interpolation
#     xnew = np.linspace(0, 1, x_dim)
#     ynew = np.linspace(0, 1, y_dim)
#     znew = np.linspace(0, 1, z_dim)
#     Xnew, Ynew, Znew = np.meshgrid(xnew, ynew, znew)

#     # Perform the interpolation
#     interpolated_data = interpolator((Xnew, Ynew, Znew))

#     return interpolated_data

# def write_dat_file(output_dir, field_name, sequence_number, data):
#     """Writes a .dat file for a specific field."""
#     os.makedirs(output_dir, exist_ok=True)  # Ensure directory exists
#     file_name = f'{output_dir}/{field_name}_{sequence_number:04d}.dat'
#     print(f"Writing {file_name}")
    
#     np.savetxt(file_name, data.reshape(-1), fmt='%.6e')

# def process_solution_files(nrb, input_pattern="sol*.dat", output_dir="./output", interpolate=False, x_dim=None, y_dim=None, z_dim=None):
#     """Processes all solution files and saves .dat files with optional interpolation."""
    
#     # Create output directories for each field
#     ice_phase_dir = os.path.join(output_dir, "IcePhase")
#     temp_dir = os.path.join(output_dir, "Temperature")
#     vapor_dir = os.path.join(output_dir, "VaporDensity")
    
#     os.makedirs(ice_phase_dir, exist_ok=True)
#     os.makedirs(temp_dir, exist_ok=True)
#     os.makedirs(vapor_dir, exist_ok=True)
    
#     # Sort input files to maintain order
#     input_files = sorted(glob.glob(input_pattern), key=lambda x: int(x.split("l")[1].split(".")[0]))

#     # Sequential numbering for output (0001, 0002, 0003, ...)
#     for sequence_number, infile in enumerate(input_files, start=1):
#         print(f"Processing {infile}...")

#         # Read the solution vector from the .dat file
#         sol = PetIGA().read_vec(infile, nrb)

#         # Debug: Print the shape of the data being processed
#         print(f"Shape of data in {infile}: {sol.shape}")

#         # Check if the data is 4D (3D spatial dimensions and 3 fields)
#         if sol.ndim == 4 and sol.shape[-1] == 3:
#             print("Processing 4D data (3D space + fields)...")

#             # Extract the fields (IcePhase, Temperature, VaporDensity)
#             ice_phase_data = sol[:, :, :, 0]  # Assuming IcePhase is field 0
#             temperature_data = sol[:, :, :, 1]  # Assuming Temperature is field 1
#             vapor_density_data = sol[:, :, :, 2]  # Assuming VaporDensity is field 2

#             # Optionally interpolate the data for 3D
#             if interpolate and x_dim is not None and y_dim is not None and z_dim is not None:
#                 ice_phase_data = interpolate_data_3d(ice_phase_data, x_dim, y_dim, z_dim)
#                 temperature_data = interpolate_data_3d(temperature_data, x_dim, y_dim, z_dim)
#                 vapor_density_data = interpolate_data_3d(vapor_density_data, x_dim, y_dim, z_dim)
        
#         else:
#             raise ValueError(f"Unsupported data dimensionality: {sol.ndim}D for {infile}. Expected 4D data (3D space + fields).")

#         # Write the .dat files for each field
#         write_dat_file(ice_phase_dir, "IcePhase", sequence_number, ice_phase_data)
#         write_dat_file(temp_dir, "Temperature", sequence_number, temperature_data)
#         write_dat_file(vapor_dir, "VaporDensity", sequence_number, vapor_density_data)

# if __name__ == "__main__":
#     nrb_file = "igasol.dat"
#     nrb = read_nrb(nrb_file)
    
#     # Example of running the script with interpolation and dimensions (optional)
#     process_solution_files(nrb, interpolate=True, x_dim=500, y_dim=500, z_dim=500)