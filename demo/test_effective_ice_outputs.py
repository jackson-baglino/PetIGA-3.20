import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# Define file directory and filename
filedir = (
    "/Users/jacksonbaglino/SimulationResults/ThermalConductivity/"
    "ThermalSim_2025-03-11__09.53.05/"
)
filename = filedir + "ice_field.bin"  # Change to "thermal_conductivity.bin" if needed

# Load binary data (Big Endian 64-bit float)
data = np.fromfile(filename, dtype='>f8')  

# Expected grid size (adjust as necessary)
Nx, Ny = 32, 32  # Ensure these match your simulation settings

# Remove the first element (assumed header)
data = data[1:]

# Check for invalid values
if np.isnan(data).any() or np.isinf(data).any():
    raise ValueError("Data contains NaN or Inf values. Check the simulation output.")

# Attempt to reshape the data, accounting for PETSc padding
try:
    reshaped = data.reshape((Ny, Nx))
    print(f"Successfully reshaped to ({Ny}, {Nx})")
except ValueError:
    raise ValueError(f"Error reshaping: Expected {(Nx) * (Ny)}, got {data.size}")

# Strip ghost cells (if needed, uncomment the next line)
# trimmed = reshaped[1:-1, 1:-1]  # Removes first and last row/column
trimmed = reshaped  # Keep full grid for now

# Create the figure with high-quality settings
fig, ax = plt.subplots(figsize=(6, 5), dpi=150)

# Display the data with a smooth, high-resolution colormap
cax = ax.imshow(trimmed, cmap='magma', origin='lower', aspect='auto', interpolation='bilinear')
cbar = fig.colorbar(cax, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label("Thermal Conductivity (W/mK)", fontsize=14, fontweight='bold')
cbar.ax.tick_params(labelsize=12)

# Add fine grid overlay (mesh)
ax.set_xticks(np.linspace(0, Nx, 5))
ax.set_yticks(np.linspace(0, Ny, 5))
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
ax.tick_params(axis='both', which='major', labelsize=12, width=1.5, length=6)
ax.tick_params(axis='both', which='minor', width=1.0, length=4)
ax.grid(which="major", color="w", linestyle='-', linewidth=0.5, alpha=0.7)

# Axis labels and title with enhanced aesthetics
ax.set_xlabel("X Index", fontsize=14, fontweight='bold', labelpad=10)
ax.set_ylabel("Y Index", fontsize=14, fontweight='bold', labelpad=10)
ax.set_title("Thermal Conductivity Field", fontsize=16, fontweight='bold', pad=15)

# Remove excessive borders for a sleek look
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(1.5)
ax.spines['bottom'].set_linewidth(1.5)

# Show the plot with tight layout for perfect spacing
plt.tight_layout()
plt.show()
