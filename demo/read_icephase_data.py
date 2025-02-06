import numpy as np
import re
import matplotlib.pyplot as plt
import os

def read_binary_simulation_data(filename, Nx, Ny, Nz, disp_flag=0):
    """
    Reads a binary simulation output file and extracts ice phase, temperature, and vapor density.

    Parameters:
        filename (str): Path to the binary data file.
        Nx (int): Number of grid points in x-direction.
        Ny (int): Number of grid points in y-direction.
        Nz (int): Number of grid points in z-direction.
        disp_flag (int, optional): Display flag (default: 0).

    Returns:
        tuple: rho_v (ice phase), phi_i (temperature), T (vapor density) as 3D NumPy arrays.
    """
    dof = 3  # Degrees of freedom (3 fields per grid point)
    data_size = dof * (Nx * Ny * Nz)  # Expected number of elements

    # Read binary data
    try:
        with open(filename, "rb") as f:
            file_data = np.fromfile(f, dtype=np.dtype('>f8'), count=data_size)
    except FileNotFoundError:
        raise FileNotFoundError(f"File {filename} could not be opened.")

    if file_data.size != data_size:
        raise ValueError(f"File {filename} does not contain the expected amount of data. "
                         f"Expected {data_size}, but got {file_data.size} elements.")

    # Extract fields
    T = file_data[0::dof]  # Ice phase (every 3rd element, starting at 0)
    rho_v = file_data[1::dof]  # Temperature (every 3rd element, starting at 1)
    phi_i = file_data[2::dof]      # Vapor density (every 3rd element, starting at 2)

    # Reshape to 3D arrays
    T = T.reshape((Nx, Ny, Nz))
    rho_v = rho_v.reshape((Nx, Ny, Nz))
    phi_i = phi_i.reshape((Nx, Ny, Nz))

    if disp_flag:
        print(f"Loaded simulation data from {filename}:")
        print(f"rho_v shape: {rho_v.shape}, phi_i shape: {phi_i.shape}, T shape: {T.shape}")

    return rho_v, phi_i, T

def plot_midplane(rho_v, phi_i, T, Nz):
    """
    Plots the mid-plane (Nz//2) of each field with colorbars and closes on any key press.
    """
    mid_z = Nz // 2

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    im1 = axes[0].imshow(phi_i[:, :, mid_z], cmap='Blues', origin='lower')
    axes[0].set_title("Ice Phase (phi_i)")
    plt.colorbar(im1, ax=axes[0])
    
    im2 = axes[1].imshow(T[:, :, mid_z], cmap='hot', origin='lower')
    axes[1].set_title("Temperature (T)")
    plt.colorbar(im2, ax=axes[1])
    
    im3 = axes[2].imshow(rho_v[:, :, mid_z], cmap='coolwarm', origin='lower')
    axes[2].set_title("Vapor Density (rho_v)")
    plt.colorbar(im3, ax=axes[2])
    
    plt.tight_layout()
    
    # Close plot on any key press
    def on_key(event):
        plt.close(fig)
    
    fig.canvas.mpl_connect("key_press_event", on_key)
    plt.show()

def main():
    """
    Main function to loop over all solution files in the directory and process the data.
    """
    directory = "/Users/jacksonbaglino/SimulationResults/DrySed_Metamorphism/NASAv2/NASAv2_2G-Molaro_3D_T-188.0_hum0.70_2025-02-05__14.35.43"
    sol_files = sorted([f for f in os.listdir(directory) if f.startswith("sol_")])
    
    Nx = 135
    Ny = 215
    Nz = 135

    for sol_file in sol_files:
        input_file = os.path.join(directory, sol_file)
        print(f"Processing file: {input_file}")
        
        try:
            results = read_binary_simulation_data(input_file, Nx, Ny, Nz, disp_flag=0)
            
            if len(results) == 3:  # 3D case
                phi_i, T, rho_v = results
                print(f"Data Loaded!")

                print("Max values in IcePhase:", phi_i.max())
                print("Max values in vapor density:", rho_v.max())
                print("Max values in Temp:", T.max())

                print("\nMin values in IcePhase:", phi_i.min())
                print("Min values in vapor density:", rho_v.min()) 
                print("Min values in Temp:", T.min())

                plot_midplane(rho_v, phi_i, T, Nz)
            else:  # 2D case
                print(f"Data was not loaded correctly.")

            # Optional: Display array shapes
            print(f"IcePhase shape: {phi_i.shape}")
            print(f"AirPhase shape: {rho_v.shape}")
            print(f"Temp shape: {T.shape}")

        except Exception as e:
            print(f"Error processing {input_file}: {e}")

if __name__ == "__main__":
    main()
