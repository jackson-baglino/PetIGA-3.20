"""
CompEps.py — Compute interface width epsilon and required mesh resolution
based on Kaempfer & Plapp (2009).

Key distinction: eps (model parameter) vs actual interface width w = 2*sqrt(2)*eps.
Mesh must resolve the actual interface width, not eps itself.
"""

import numpy as np

# =========================================================================
# Physical constants
# =========================================================================
alpha   = 1.0e-2        # Condensation coefficient (1e-3 to 1e-1)
m       = 2.99e-26      # Mass of water molecule (kg)
k_B     = 1.38e-23      # Boltzmann constant (J/K)
K_i     = 2.29          # Thermal conductivity of ice (W/m/K)
C_i     = 1.8e6         # Volumetric heat capacity of ice (J/m^3/K)
rho_ice = 919.0         # Density of ice (kg/m^3)
Dv0     = 2.178e-5      # Vapor diffusion coefficient (m^2/s)

# =========================================================================
# Tunable parameters
# =========================================================================
T0   = 273.15 - 20.0   # Mean temperature (K)
Rave = 3.0e-5          # Representative grain radius (m)
Lx   = 1.25e-4         # Domain size x (m)
Ly   = 2.5e-4          # Domain size y (m)
Lz   = 0.0             # Domain size z (m), 0 for 2D

# =========================================================================
# Saturation vapor density
# =========================================================================
Kj      = [-0.5865e4, 0.2224e2, 0.1375e-1, -0.3403e-4, 0.2697e-7, 0.6918]
Patm    = 101325.0
rho_air = 1.341
bb      = 0.62

exponent = sum(Kj[j] * T0**(j - 1) for j in range(5)) + Kj[5] * np.log(T0)
Pvs      = np.exp(exponent)
rho_vs   = rho_air * bb * Pvs / (Patm - Pvs)
rho_ratio = rho_vs / rho_ice

# =========================================================================
# Kinetic coefficient
# =========================================================================
beta_prime = (1.0 / alpha) * np.sqrt(2.0 * np.pi * m / (k_B * T0))
beta0      = beta_prime / rho_ratio

# =========================================================================
# Upper bounds on eps from Kaempfer & Plapp (2009) eqs. 44-46
# These are UPPER bounds — eps must be SMALLER than these values
# =========================================================================
eps_heat  = (K_i / C_i) * rho_ratio * beta0
eps_vapor = Dv0          * rho_ratio * beta0
eps_geom  = Rave          # eps << Rave for sharp-interface limit

eps_max = min(eps_heat, eps_vapor, eps_geom)

# =========================================================================
# Actual interface width in the simulation
# For equilibrium tanh profile: w_actual = 2*sqrt(2)*eps ~ 2.83*eps
# This is what must be resolved by the mesh, not eps itself
# =========================================================================
w_actual = 2.0 * np.sqrt(2.0) * eps_max

# =========================================================================
# Mesh sizing — matching your MATLAB script logic
# Your MATLAB script uses: Nx = ceil(Lx / eps_choose)
# This gives one element per eps, meaning the actual interface
# spans ~2.83 elements — the bare minimum for any resolution at all.
#
# For reliable results you want n_per_interface elements across w_actual:
#   dx = w_actual / n_per_interface
#   Nx = ceil(Lx / dx)
#
# With p=2 C=1 IGA, 4-5 elements across w_actual is generally sufficient.
# With p=3 C=2, 3-4 elements may suffice due to higher accuracy per element.
# =========================================================================

# Reproduce your MATLAB result exactly:
Nx_matlab = int(np.ceil(Lx / eps_max))
Ny_matlab = int(np.ceil(Ly / eps_max))
Nz_matlab = int(np.ceil(Lz / eps_max)) if Lz > 0 else 0

# Better estimate: resolve actual interface width with n elements
for n_per_interface in [3, 4, 5, 8]:
    dx  = w_actual / n_per_interface
    Nx  = int(np.ceil(Lx / dx))
    Ny  = int(np.ceil(Ly / dx))

# =========================================================================
# Output
# =========================================================================
print("=" * 65)
print("  Interface width and mesh sizing — Kaempfer & Plapp (2009)")
print("=" * 65)

print(f"\n--- Thermodynamic state at T0 = {T0:.2f} K ({T0-273.15:.1f} C) ---")
print(f"  rho_vs                 = {rho_vs:.4e} kg/m^3")
print(f"  rho_vs / rho_ice       = {rho_ratio:.4e}")
print(f"  beta0 (normalized)     = {beta0:.4e} s/m")

print(f"\n--- Epsilon upper bounds (eps must be LESS than these) ---")
print(f"  From heat diffusion    = {eps_heat:.4e} m")
print(f"  From vapor diffusion   = {eps_vapor:.4e} m")
print(f"  From grain geometry    = {eps_geom:.4e} m  (Rave)")
print(f"  Maximum allowable eps  = {eps_max:.4e} m")
print(f"  Actual interface width = {w_actual:.4e} m  (= 2*sqrt(2)*eps)")
print(f"  eps / Rave             = {eps_max/Rave:.4f}  (should be << 1)")

print(f"\n--- Mesh: MATLAB-equivalent (1 element per eps) ---")
print(f"  dx = eps             = {eps_max:.4e} m")
print(f"  Nx = {Nx_matlab},  Ny = {Ny_matlab}")
print(f"  Elements across interface (w_actual/dx) = {w_actual/eps_max:.1f}")
print(f"  NOTE: ~2.8 elements across interface is very coarse.")
print(f"  May be adequate for p=2 or p=3 IGA but worth checking.")

print(f"\n--- Mesh: resolving actual interface width ---")
print(f"  {'n/interface':<14} {'dx (m)':<14} {'Nx':<8} {'Ny':<8}")
print(f"  {'-'*50}")
for n in [3, 4, 5, 8]:
    dx = w_actual / n
    Nx = int(np.ceil(Lx / dx))
    Ny = int(np.ceil(Ly / dx))
    print(f"  {n:<14} {dx:<14.4e} {Nx:<8} {Ny:<8}")

print(f"\n--- Validity checks ---")
vn_typical = 1.0e-9
print(f"  eps * vn / Dv          = {eps_max * vn_typical / Dv0:.4e}  (should be << 1)")
if eps_max / Rave < 0.01:
    print(f"  eps/Rave check: OK (< 0.01)")
elif eps_max / Rave < 0.1:
    print(f"  eps/Rave check: CAUTION (0.01 to 0.1)")
else:
    print(f"  eps/Rave check: WARNING (> 0.1, sharp-interface limit questionable)")

print("=" * 65)