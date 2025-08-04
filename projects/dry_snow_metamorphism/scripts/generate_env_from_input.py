import os
import sys
import select
import numpy as np
from time import time, sleep

# -----------------------------
# Physical Constants
# -----------------------------
alpha = 1e-2
m = 2.99e-26  # kg
k = 1.38e-23  # J/K
K_i = 2.29    # W/m/K
K_a = 0.02
C_i = 1.8e6   # J/m^3/K
C_a = 1.044e3
rhoI = 919.0  # kg/m^3
Dv0 = 2.178e-5  # m^2/s
Kj = [-0.5865e4, 0.2224e2, 0.1375e-1, -0.3403e-4, 0.2697e-7, 0.6918]
Patm = 1.013250
rhoatm = 1.341
Rave = 8.2883e-05

# -----------------------------
# Input Args
# -----------------------------
if len(sys.argv) < 5:
    print("Usage: python generate_env_from_input.py <input_file> <output_env_path> <Lx> <Ly> [Lz]")
    sys.exit(1)

input_path = sys.argv[1]
output_path = sys.argv[2]
Lx = float(sys.argv[3])
Ly = float(sys.argv[4])
Lz = None
if len(sys.argv) >= 6:
    Lz = float(sys.argv[5])

print(f"[INFO] Reading grain centers from {input_path}")
print(f"[INFO] Writing .env file to {output_path}")

# -----------------------------
# Read grain centers
# -----------------------------
def read_grain_centers(filepath):
    centers = []
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.replace(',', ' ').split()
            if len(parts) >= 3:
                try:
                    x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
                    centers.append((abs(x), abs(y), abs(z)))
                except ValueError:
                    continue
    return np.array(centers)

centers = read_grain_centers(input_path)

# -----------------------------
# Determine Lz
# -----------------------------
if Lz is None:
    max_coords = centers.max(axis=0)
    Lz = max_coords[2]
    print(f"[INFO] Auto-detected Lz: {Lz:.2e}")
else:
    print(f"[INFO] Using user-specified Lz: {Lz:.2e}")

# -----------------------------
# Compute epsilon
# -----------------------------
T0 = 253.15  # K (default for T = -20 C)

PvsT = sum(Kj[i] * T0**(i - 1) for i in range(5))
PvsT = np.exp(PvsT + Kj[-1] * np.log(T0))
xsT = 0.62 * PvsT / (Patm - PvsT)
rhoVS = rhoatm * xsT
rhoVS_rhoI = rhoVS / rhoI
betaP0 = (1 / alpha) * np.sqrt((2 * np.pi * m) / (k * T0))
beta0 = betaP0 / rhoVS_rhoI

eps1 = (K_i / C_i) * rhoVS_rhoI * beta0
eps2 = Dv0 * rhoVS_rhoI * beta0
eps3 = 1 / (1 / Rave)

eps = min(eps1, eps2, eps3) * 0.975

# -----------------------------
# Compute Nx, Ny, Nz
# -----------------------------
factor = 2  # Factor to ensure domain size is a multiple of eps
Nx = int(np.ceil(Lx / eps / factor))
Ny = int(np.ceil(Ly / eps / factor))
Nz = int(np.ceil(Lz / eps / factor))

# -----------------------------
# Write to .env file
# -----------------------------
with open(output_path, 'w') as f:
    f.write(f"Lx={Lx:.6e}\n")
    f.write(f"Ly={Ly:.6e}\n")
    f.write(f"Lz={Lz:.6e}\n")
    f.write(f"Nx={Nx}\n")
    f.write(f"Ny={Ny}\n")
    f.write(f"Nz={Nz}\n")
    f.write(f"eps={eps:.12e}\n")

print(f"[INFO] .env file written to {output_path}")