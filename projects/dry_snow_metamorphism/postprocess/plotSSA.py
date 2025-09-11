import os
import re
import numpy as np
import matplotlib.pyplot as plt

# ----------------------------
# Configuration from environment
# ----------------------------
# Required: inputFile (path to grains.dat)
inputFile = os.getenv("inputFile")
# Optional: dim ("2" or "3"), title for figure name, output folder
dim = os.getenv("dim", "2").strip()
title = os.getenv("title", "").strip() or "run"
outputfolder = os.getenv("outputfolder")

print(f"inputFile:  {inputFile}")

# ----------------------------
# Resolve output path & safe filename
# ----------------------------
safe_title = re.sub(r"[\\/]+", "_", title).strip() or "run"
output_dir = outputfolder if outputfolder else os.getcwd()
os.makedirs(output_dir, exist_ok=True)
output_path = os.path.join(output_dir, f"ssa_evolution_plot_{safe_title}.png")

# ----------------------------
# Load SSA evolution data
# Expected columns: area, volume, time(sec)
# ----------------------------
ssa_path = os.path.join(os.getcwd(), "SSA_evo.dat")
try:
    ssa = np.loadtxt(ssa_path)
except Exception as e:
    raise RuntimeError(f"Failed to read SSA data from {ssa_path}: {e}")

if ssa.ndim != 2 or ssa.shape[1] < 3:
    raise ValueError(f"SSA_evo.dat must have at least 3 columns (area, volume, time[sec]). Got shape {ssa.shape}.")

area_data = ssa[:, 0]
volume_data = ssa[:, 1]
time_data = ssa[:, 2] / 3600.0  # seconds -> hours
print(f"The number of time_data points is: {len(time_data)}")

# ----------------------------
# Load grain radii from grains.dat
# The generator writes a header line (Lx Ly Lz) then rows of centers/radius.
# We read the LAST column as radius, and skip the first header row.
# Works for either 3-col (x y r) or 4-col (x y z r) formats.
# ----------------------------
if not inputFile or not os.path.isfile(inputFile):
    raise FileNotFoundError(f"grains.dat not found at: {inputFile}")

try:
    radii = np.loadtxt(inputFile, comments="#", usecols=(-1,), ndmin=1, skiprows=1)
except Exception:
    # Fallback without skipping header in case file lacks header
    radii = np.loadtxt(inputFile, comments="#", usecols=(-1,), ndmin=1)

if radii.size == 0:
    raise ValueError("No radii parsed from grains.dat; check file format.")

# ----------------------------
# Compute initial area/volume for normalization
# ----------------------------
if dim == "2":
    area0 = np.sum(2.0 * np.pi * radii)
    volume0 = np.sum(np.pi * radii ** 2)
elif dim == "3":
    area0 = np.sum(4.0 * np.pi * radii ** 2)
    volume0 = np.sum((4.0 / 3.0) * np.pi * radii ** 3)
else:
    raise ValueError(f"Unsupported dim={dim}. Expected '2' or '3'.")

if area0 <= 0 or volume0 <= 0:
    raise ValueError(f"Non-positive area0/volume0 (area0={area0}, volume0={volume0}). Check radii.")

# ----------------------------
# Normalize and compute SSA
# ----------------------------
area_norm = area_data / area0
volume_norm = volume_data / volume0
ssa_data = area_norm / volume_norm
ssa_data = ssa_data / ssa_data[0]

# ----------------------------
# Plot
# ----------------------------
plt.figure(figsize=(10, 6))
plt.plot(time_data, ssa_data, label="Specific Surface Area (normalized)")
plt.xlabel("Time [hours]", fontsize=14)
plt.ylabel("SSA / SSA$_0$", fontsize=14)
plt.title("Surface Area Evolution", fontsize=18)
plt.grid(True, alpha=0.25)
plt.legend()
plt.tight_layout()

# Save figure
plt.savefig(output_path)
print(f"[ok] Wrote: {output_path}")

# Non-blocking close (safe in headless runs)
plt.close()