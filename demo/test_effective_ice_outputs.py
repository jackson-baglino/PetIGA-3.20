import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import zoom

# Define the base directory
base_dir = "/Users/jacksonbaglino/SimulationResults/ThermalConductivity/"\
           "ThermalSim_2025-03-17__15.45.21"

# File paths for temperature and ice phase fields
temp_file = f"{base_dir}/temperature.bin"
# ice_file = f"{base_dir}/ice_field.bin"

# Load binary data (Big Endian 64-bit float)
temperature = np.fromfile(temp_file, dtype='>f8')
# ice_field = np.fromfile(ice_file, dtype='>f8')

# Grid size (Ensure these match the simulation settings)
Nx, Ny = 256, 256
Nx_temp, Ny_temp = Nx + 1, Ny + 1  # Structured grid (Temperature)
Nx_ice,  Ny_ice  = Nx + 1, Ny + 1     # Staggered grid (Ice field)

# Print the number of data points for debugging
print(f"Temperature data points: {temperature.size}")
# print(f"Ice field data points: {ice_field.size}")

# Ensure data sizes match expected grid dimensions
while temperature.size > Nx_temp * Ny_temp:
    temperature = temperature[1:]  # Remove extra points if needed
    print("Removed extra data point from temperature")

# while ice_field.size > Nx_ice * Ny_ice:
#     ice_field = ice_field[1:]  # Remove extra points if needed
#     print("Removed extra data point from ice field")

# Reshape to 2D arrays
temperature = temperature.reshape((Ny_temp, Nx_temp))
# ice_field = ice_field.reshape((Ny_ice, Nx_ice))

# Interpolate ice field to match temperature grid resolution
scale_x = Nx_temp / Nx_ice
scale_y = Ny_temp / Ny_ice
# ice_interp = zoom(ice_field, (scale_y, scale_x), order=1)  # Linear interpolation

# Create figure with two subplots
fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=150)

# Plot Ice Field
ax1 = axes[0]
# im1 = ax1.imshow(ice_field, cmap='Blues', origin='lower', interpolation="nearest")
# fig.colorbar(im1, ax=ax1, label="Ice Phase")
ax1.set_xlabel("X Index")
ax1.set_ylabel("Y Index")
ax1.set_title("Ice Phase Field")

# Plot Temperature Field with Ice Contour
ax2 = axes[1]
im2 = ax2.imshow(temperature, cmap='magma', origin='lower', interpolation="nearest")
fig.colorbar(im2, ax=ax2, label="Temperature (K)")
ax2.set_xlabel("X Index")
ax2.set_ylabel("Y Index")
ax2.set_title("Temperature Field with Ice Contour")

# Overlay the ice phase contour at 0.5
contour_levels = [0.5]  # Contour at IcePhase = 0.5
# ax2.contour(ice_interp, levels=contour_levels, colors='white', linewidths=1.5)

# Show both plots
plt.tight_layout()
plt.show()