import numpy as np
import sys
import matplotlib.pyplot as plt

# input_file="/Users/jacksonbaglino/SimulationResults/DrySed_Metamorphism/NASAv2/NASAv2_10G_2D_T-20.0_hum0.70_2025-03-13__14.20.59/sol_00430.dat"  
# output_file="/Users/jacksonbaglino/SimulationResults/ThermalConductivity/ice_field.txt"  
# Nx=276  # Grid size in X
# Ny=276  # Grid size in Y

def extract_ice_field(input_file, output_file, dof=3):
    """
    Extracts the ice phase field (DOF 0) from a PETSc binary solution file.
    
    Args:
        input_file (str): Path to the PETSc binary file (sol_xxxxx.dat).
        output_file (str): Path to the output text file for the ice field.
        dof (int): Degrees of freedom per node in the PETSc vector (default is 3).
    """
    try:
        # Read binary data (Big Endian 64-bit float format)
        data = np.fromfile(input_file, dtype='>f8')
        data = data[1:]  # Skip header if needed

        # Ensure the file size is a multiple of the number of DOFs
        if len(data) % dof != 0:
            print(f"Error: File size {len(data)} is not a multiple of DOF {dof}.")
            sys.exit(1)

        # Reshape and extract the first DOF (ice field)
        num_points = len(data) // dof
        ice_field = data.reshape(num_points, dof)[:, 0]  # Extract DOF 0

        # Save the extracted data to a text file
        np.savetxt(output_file, ice_field, fmt="%.8e")
        print(f"Successfully extracted ice field from {input_file} and saved to {output_file}")

        return ice_field

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

def plot_ice_field(ice_field, Nx, Ny):
    """
    Plots the ice field and keeps the plot window open until manually closed.

    Args:
        ice_field (numpy array): 1D array containing ice field values.
        Nx (int): Number of grid points in the x-direction.
        Ny (int): Number of grid points in the y-direction.
    """
    try:
        # Reshape into 2D array
        ice_field_2d = ice_field.reshape((Ny, Nx))

        # Create the figure
        plt.figure(figsize=(6, 5))
        plt.imshow(ice_field_2d, cmap="Blues", origin="lower", interpolation="nearest")
        plt.colorbar(label="Ice Phase")
        plt.title("Extracted Ice Field")
        plt.xlabel("X Index")
        plt.ylabel("Y Index")

        # Keep plot open indefinitely
        print("Close the plot window to continue...")
        plt.show(block=True)  # Keeps the window open until manually closed

    except Exception as e:
        print(f"Plotting error: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python extract_ice_field.py <input_file> <output_file> <Nx> <Ny>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    Nx = int(sys.argv[3])  # Grid size in X
    Ny = int(sys.argv[4])  # Grid size in Y

    # Extract ice field data
    ice_field = extract_ice_field(input_file, output_file)

    # Show the plot (this will now stay open until you manually close it)
    plot_ice_field(ice_field, Nx, Ny)