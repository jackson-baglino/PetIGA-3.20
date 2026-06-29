#!/usr/bin/env python3
"""
generate_opts_from_input.py — compute mesh parameters from a grains.dat file
and write them to a PETSc .opts file.

Usage:
    python3 preprocess/generate_opts_from_input.py <grains.dat> <output.opts> <Lx> <Ly> [Lz]

Replaces the legacy generate_env_from_input.py (which wrote shell .env format).
"""
import sys
import numpy as np

# ---- Physical constants ----
alpha  = 1e-2
m      = 2.99e-26   # kg, water molecule mass
k      = 1.38e-23   # J/K, Boltzmann
K_i    = 2.29       # W/m/K, ice thermal conductivity
C_i    = 1.8e6      # J/m³/K, ice volumetric heat capacity
rhoI   = 919.0      # kg/m³, ice density
Dv0    = 2.178e-5   # m²/s, vapor diffusivity at reference T
Kj     = [-0.5865e4, 0.2224e2, 0.1375e-1, -0.3403e-4, 0.2697e-7, 0.6918]
Patm   = 1.013250   # atm
rhoatm = 1.341      # kg/m³, air density

if len(sys.argv) < 5:
    print("Usage: python3 generate_opts_from_input.py <grains.dat> <output.opts> <Lx> <Ly> [Lz]")
    sys.exit(1)

input_path  = sys.argv[1]
output_path = sys.argv[2]
Lx = float(sys.argv[3])
Ly = float(sys.argv[4])
Lz = float(sys.argv[5]) if len(sys.argv) >= 6 else None

print(f"[INFO] Reading grain centers from: {input_path}")

def read_grain_centers(filepath):
    centers, radii = [], []
    with open(filepath) as f:
        for line in f:
            parts = line.replace(',', ' ').split()
            if len(parts) >= 4:
                try:
                    x, y, z, r = float(parts[0]), float(parts[1]), float(parts[2]), float(parts[3])
                    centers.append((abs(x), abs(y), abs(z)))
                    radii.append(abs(r))
                except ValueError:
                    continue
    return np.array(centers), np.array(radii)

centers, radii = read_grain_centers(input_path)
if len(radii) == 0:
    print("[ERROR] No grain data found in input file.")
    sys.exit(1)

Rave = np.mean(radii)
print(f"[INFO] Number of grains: {len(radii)},  mean radius: {Rave:.3e} m")

if Lz is None:
    max_z = centers[:, 2].max()
    Lz = 2 * max_z if max_z > 0 else 1.0e-4
    print(f"[INFO] Auto-detected Lz: {Lz:.3e} m")

# ---- Compute eps from sharp-interface capillary length ----
T0 = 273.15  # K (0°C); lower T gives larger eps — conservative choice
PvsT = sum(Kj[i] * T0**(i - 1) for i in range(5))
PvsT = np.exp(PvsT + Kj[-1] * np.log(T0))
xsT  = 0.62 * PvsT / (Patm - PvsT)
rhoVS = rhoatm * xsT
rhoVS_rhoI = rhoVS / rhoI
betaP0 = (1 / alpha) * np.sqrt((2 * np.pi * m) / (k * T0))
beta0  = betaP0 / rhoVS_rhoI

eps1 = (K_i / C_i) * rhoVS_rhoI * beta0
eps2 = Dv0 * rhoVS_rhoI * beta0
eps3 = Rave  # cannot exceed grain radius
eps  = min(eps1, eps2, eps3) * 0.975

print(f"[INFO] eps = {eps:.6e} m  (eps1={eps1:.2e}, eps2={eps2:.2e}, eps3(Rave)={eps3:.2e})")

# ---- Mesh resolution: ~2 elements per eps ----
factor = 2
Nx = int(np.ceil(Lx / eps / factor))
Ny = int(np.ceil(Ly / eps / factor))
Nz = int(np.ceil(Lz / eps / factor))
dim = 2 if Nz <= 1 else 3

print(f"[INFO] Mesh: Nx={Nx}  Ny={Ny}  Nz={Nz}  dim={dim}")

with open(output_path, 'w') as f:
    f.write(f"# Auto-generated geometry opts for {input_path}\n\n")
    f.write(f"-dim {dim}\n")
    f.write(f"-Lx {Lx:.6e}\n")
    f.write(f"-Ly {Ly:.6e}\n")
    f.write(f"-Lz {Lz:.6e}\n")
    f.write(f"-Nx {Nx}\n")
    f.write(f"-Ny {Ny}\n")
    f.write(f"-Nz {Nz}\n")
    f.write(f"-eps {eps:.12e}\n")
    f.write(f"-readFlag 1\n")
    f.write(f"-grains_file {input_path}\n")

print(f"[INFO] .opts file written to: {output_path}")
