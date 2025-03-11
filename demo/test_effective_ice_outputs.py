import numpy as np
import matplotlib.pyplot as plt

# Define the file path
filename = "/Users/jacksonbaglino/SimulationResults/ThermalConductivity/"\
    "ThermalSim_2025-03-11__12.51.25" \
    "/ice_field.bin"

# Load binary data (Big Endian 64-bit float)
data = np.fromfile(filename, dtype='>f8')

# Grid size (adjust as necessary)
Nx, Ny = 1000, 1000  # Ensure these match your simulation settings

# Reshape data to 2D array
while data.size > Nx * Ny:
    data = data[1:]
    print("Removed extra data point")

for i in range(Nx):
    for j in range(Ny):
        print(f"{data[i * Nx + j]:0.2f}", end=" ")
    print()

data = data.reshape((Ny, Nx))

# Plot the data without smoothing/interpolation
plt.figure(figsize=(6, 5), dpi=150)
plt.imshow(data, cmap='magma', origin='lower', interpolation="nearest")
plt.colorbar(label="Thermal Conductivity (W/mK)")
plt.xlabel("X Index")
plt.ylabel("Y Index")
plt.title("Thermal Conductivity Field")
plt.show()