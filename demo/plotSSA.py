from igakit.io import PetIGA, VTK
import glob
import os

def read_nrb(file_path):
    """Reads the IGA solution file."""
    return PetIGA().read(file_path)

def generate_vtk(nrb, infile, output_dir, sequence_number):
    """Generates the VTK file with sequential numbering (4 digits) and returns the path."""
    file_number = f'{sequence_number:04d}'  # Ensures 4-digit numbering (0001, 0002, etc.)
    output_file = f'{output_dir}/g_{file_number}.vtk'

    if not os.path.isfile(output_file):
        print(f"Generating VTK: {infile} -> {output_file}")
        sol = PetIGA().read_vec(infile, nrb)
        VTK().write(output_file,  
                    nrb,             
                    fields=sol,     
                    scalars={'IcePhase': 0, 'Temperature': 1, 'VaporDensity': 2})
    return output_file

def process_solution_files(nrb, input_pattern="sol*.dat", output_dir="./vtkOut", contour_value=0.5):
    """Processes all solution files, generates VTK with sequential numbering, and saves them."""
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Sort input files to maintain order
    input_files = sorted(glob.glob(input_pattern), key=lambda x: int(x.split("l")[1].split(".")[0]))

    # Sequential numbering for VTK output (g_0001, g_0002, g_0003, ...)
    for sequence_number, infile in enumerate(input_files, start=1):
        generate_vtk(nrb, infile, output_dir, sequence_number)

if __name__ == "__main__":
    nrb_file = "igasol.dat"
    nrb = read_nrb(nrb_file)
    process_solution_files(nrb)