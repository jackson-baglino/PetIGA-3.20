from igakit.io import PetIGA, VTK
import glob
import os

# Read the IGA solution data structure from the input file
nrb = PetIGA().read("igasol.dat")

# Loop through all solution files matching the pattern "sol*.dat"
for infile in glob.glob("sol*.dat"):
    # Extract the base name and numeric identifier from the file name
    name = infile.split(".")[0]  # Get file name without extension
    number = name.split("l")[1].zfill(4)  # Add leading zeros to make it 4 digits

    # Construct the output file path
    root = f'./solV{number}.vtk'
    
    # Check if the output file already exists to avoid reprocessing
    if not os.path.isfile(root):
        # Read the solution vector for the current file using the PetIGA library
        sol = PetIGA().read_vec(infile, nrb)
        
        # Define the output file path within the vtkOut directory
        outfile = f'./vtkOut/solV{number}.vtk'
        
        # Write the solution to a VTK file, including scalar fields for visualization
        VTK().write(
            outfile,      # Output file path
            nrb,          # IGA structure
            fields=sol,   # Solution vector fields
            scalars={     # Map scalar field names to their indices
                'IcePhase': 0,
                'Temperature': 1,
                'VaporDensity': 2
            }
        )

# Notes:
# - The igakit.io.PetIGA class reads IGA data structures and solution vectors.
# - The igakit.io.VTK class converts IGA data into VTK format for visualization.
# - The scalars dictionary maps descriptive field names to their indices.
# - Ensure the "vtkOut" directory exists before running the script.