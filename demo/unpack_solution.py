import numpy as np
import sys
import os
import argparse
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

def unpack_solution(output_path, Nx, Ny, Nz, dof=3, save_plots=False):
    """
    Unpacks the ice field and, if applicable, the temperature and vapor density fields from a binary PETSc Vec solution file.

    Args:
        output_path (str): Path to the binary PETSc Vec file.
        Nx, Ny, Nz (int): Number of elements in each spatial dimension.
        dof (int): Degrees of freedom (1 for Ice only, 3 for Ice, Temperature, Vapor Density).
        save_plots (bool): If True, saves plots of the extracted fields.

    Outputs:
        Saves the extracted fields as `.dat` files in the same directory.
        Saves images of the fields if save_plots=True.
    """
    # Ensure file exists
    if not os.path.exists(output_path):
        print(f"❌ Error: File {output_path} does not exist.")
        sys.exit(1)

    # Read binary file (Big Endian 64-bit float)
    data = np.fromfile(output_path, dtype=">f8")

    # Remove the first data point corresponding to the IGA mesh
    data = data[1:]

    num_points = Nx * Ny * (Nz if Nz > 1 else 1)

    # Ensure expected data size matches
    expected_size = num_points * dof
    if len(data) != expected_size:
        print(f"❌ Error: Expected {expected_size} values but found {len(data)}.")
        print(f"   Check your Nx, Ny, Nz values: {Nx}, {Ny}, {Nz}")
        sys.exit(1)

    # Reshape data into (num_points, dof)
    data = data.reshape((num_points, dof))

    # Extract fields based on DOF
    ice_field = data[:, 0]  # Ice Phase
    temperature_field = data[:, 1] if dof == 3 else None  # Temperature
    vapor_density_field = data[:, 2] if dof == 3 else None  # Vapor Density

    # Reshape into structured 3D or 2D arrays
    def reshape_field(field):
        return field.reshape((Nz, Ny, Nx)) if Nz > 1 else field.reshape((Ny, Nx))

    ice_field = reshape_field(ice_field)
    temperature_field = reshape_field(temperature_field) if dof == 3 else None
    vapor_density_field = reshape_field(vapor_density_field) if dof == 3 else None

    # Define output directory
    base_dir = os.path.dirname(output_path)

    # Define output file paths
    ice_file = os.path.join(base_dir, "ice_field.dat")
    np.savetxt(ice_file, ice_field.flatten(), fmt="%.8e", delimiter=" ")
    print(f"✅ Ice Field saved: {ice_file}")

    if dof == 3:
        temp_file = os.path.join(base_dir, "temperature_field.dat")
        vapor_file = os.path.join(base_dir, "vapor_density_field.dat")
        np.savetxt(temp_file, temperature_field.flatten(), fmt="%.8e", delimiter=" ")
        np.savetxt(vapor_file, vapor_density_field.flatten(), fmt="%.8e", delimiter=" ")
        print(f"✅ Temperature Field saved: {temp_file}")
        print(f"✅ Vapor Density Field saved: {vapor_file}")

    # Save plots if enabled
    if save_plots:
        save_plot_images(ice_field, temperature_field, vapor_density_field, Nx, Ny, Nz, base_dir, dof)


def save_plot_images(ice, temp, vapor, Nx, Ny, Nz, save_dir, dof):
    """
    Saves the plots of the Ice Phase, and optionally Temperature and Vapor Density fields.

    Args:
        ice, temp, vapor (np.ndarray): Field data arrays.
        Nx, Ny, Nz (int): Grid dimensions.
        save_dir (str): Directory to save the images.
        dof (int): Degrees of freedom (1 for Ice only, 3 for all).
    """
    mid_slice = Nz // 2 if Nz > 1 else 0
    slice_idx = mid_slice if Nz > 1 else slice(None)

    def set_standard_format(cb):
        cb.formatter = ScalarFormatter(useMathText=False)
        cb.formatter.set_scientific(False)
        cb.update_ticks()

    def set_scientific_format(cb):
        cb.formatter = ScalarFormatter(useMathText=True)
        cb.formatter.set_scientific(True)
        cb.update_ticks()

    # Plot Ice Field
    plt.figure(figsize=(6, 5))
    im = plt.imshow(ice[slice_idx], cmap="Blues", origin="lower", interpolation="nearest", aspect="auto")
    cb = plt.colorbar(im, label="Ice Phase")
    set_standard_format(cb)
    plt.xlabel("X Index")
    plt.ylabel("Y Index")
    plt.title("Ice Phase Field")
    plt.savefig(os.path.join(save_dir, "ice_field.png"), dpi=150)
    plt.close()

    if dof == 3:
        # Plot Temperature Field
        plt.figure(figsize=(6, 5))
        im = plt.imshow(temp[slice_idx], cmap="magma", origin="lower", interpolation="nearest", aspect="auto")
        cb = plt.colorbar(im, label="Temperature (K)")
        set_standard_format(cb)
        plt.xlabel("X Index")
        plt.ylabel("Y Index")
        plt.title("Temperature Field")
        plt.savefig(os.path.join(save_dir, "temperature_field.png"), dpi=150)
        plt.close()

        # Plot Vapor Density Field (Scientific Notation)
        plt.figure(figsize=(6, 5))
        im = plt.imshow(vapor[slice_idx], cmap="coolwarm", origin="lower", interpolation="nearest", aspect="auto")
        cb = plt.colorbar(im, label="Vapor Density")
        set_scientific_format(cb)
        plt.xlabel("X Index")
        plt.ylabel("Y Index")
        plt.title("Vapor Density Field")
        plt.savefig(os.path.join(save_dir, "vapor_density_field.png"), dpi=150)
        plt.close()

    print(f"📸 Saved plots in {save_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Unpack and save PETSc solution fields from a binary file.")
    parser.add_argument("output_file", type=str, help="Path to the binary PETSc Vec file.")
    parser.add_argument("Nx", type=int, help="Number of elements in X direction.")
    parser.add_argument("Ny", type=int, help="Number of elements in Y direction.")
    parser.add_argument("Nz", type=int, help="Number of elements in Z direction (use 1 for 2D).")
    parser.add_argument("--dof", type=int, choices=[1, 3], default=3, help="Degrees of freedom (1 for Ice only, 3 for Ice, Temp, and Vapor).")
    parser.add_argument("--save_plots", action="store_true", help="Save images of the fields.")

    args = parser.parse_args()

    unpack_solution(args.output_file, args.Nx, args.Ny, args.Nz, dof=args.dof, save_plots=args.save_plots)
