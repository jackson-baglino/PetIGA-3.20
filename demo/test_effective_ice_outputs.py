import numpy as np
import matplotlib.pyplot as plt
import os

# ============================
# 🔧 User-Defined Paths (Hardcoded)
# ============================
folder = "ThermalSim_2025-04-01__11.39.27"
base_dir = "/Users/jacksonbaglino/SimulationResults/ThermalConductivity/" + folder
save_dir = "/Users/jacksonbaglino/PetIGA-3.20/demo/input/Thermal_IO"
ice_field_file = "/Users/jacksonbaglino/SimulationResults/ThermalConductivity/"+ folder + "/ice_field.dat"  # Ice field output from Python script

# Ensure save directory exists
os.makedirs(save_dir, exist_ok=True)

# ============================
# 📂 File Paths
# ============================
temp_file = os.path.join(base_dir, "temperature.bin")

# ============================
# 📥 Load Temperature Data
# ============================
temperature = np.fromfile(temp_file, dtype=">f8")

# Grid size (Ensure these match the simulation settings)
Nx, Ny = 128, 128
Nx_temp, Ny_temp = Nx+1, Ny+1  # Structured grid (Temperature)

# Debugging: Print number of data points
print(f"Temperature data points: {temperature.size}")

# Ensure data sizes match expected grid dimensions
while temperature.size > Nx_temp * Ny_temp:
    temperature = temperature[1:]  # Remove extra points if needed
    print("Removed extra data point from temperature")

# Reshape temperature data into 2D array
temperature = temperature.reshape((Ny_temp, Nx_temp))

# ============================
# 📥 Load Ice Field Data from .dat File
# ============================
try:
    ice_field = np.loadtxt(ice_field_file)  # Read from space-separated .dat file
    print(f"The size of the ice field data: {ice_field.size}")
except Exception as e:
    print(f"❌ Error reading ice field file: {e}")
    exit(1)

# Debugging: Print shape of ice field
print(f"Ice field raw shape: {ice_field.shape}")

# Ensure the ice field is reshaped into a 2D array (Ny, Nx)
ice_field = ice_field.reshape((Ny, Nx))

# Debugging: Print the reshaped ice field
print(f"Ice field reshaped to: {ice_field.shape}")

# ============================
# 📊 Generate and Save Plots
# ============================
fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=150)

# 🌊 Ice Phase Field
ax1 = axes[0]
im1 = ax1.imshow(ice_field, cmap="Blues", origin="lower", interpolation="nearest", aspect="auto", extent=[0, Nx, 0, Ny])
fig.colorbar(im1, ax=ax1, label="Ice Phase")
ax1.set_xlabel("X Index")
ax1.set_ylabel("Y Index")
ax1.set_title("Ice Phase Field")

# Save Ice Phase plot
ice_plot_path = os.path.join(save_dir, "ice_phase_TS.png")
fig.savefig(ice_plot_path, dpi=150, bbox_inches="tight")
print(f"✅ Ice phase plot saved: {ice_plot_path}")

# 🔥 Temperature Field with Ice Contour
ax2 = axes[1]
im2 = ax2.imshow(temperature, cmap="magma", origin="lower", interpolation="nearest", aspect="auto", extent=[0, Nx, 0, Ny])
fig.colorbar(im2, ax=ax2, label="Temperature (K)")
ax2.set_xlabel("X Index")
ax2.set_ylabel("Y Index")
ax2.set_title("Temperature Field with Ice Contour")

# Overlay the ice phase contour at 0.5
contour_levels = [0.5]  # Contour at IcePhase = 0.5
ax2.contour(ice_field, levels=contour_levels, colors="red", linewidths=1.5, extent=[0, Nx, 0, Ny])

# Save Temperature plot
temp_plot_path = os.path.join(save_dir, "temperature_field_TS.png")
fig.savefig(temp_plot_path, dpi=150, bbox_inches="tight")
print(f"✅ Temperature plot saved: {temp_plot_path}")

plt.close(fig)  # Close figure to free memory
