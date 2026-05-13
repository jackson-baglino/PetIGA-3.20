"""
comp_eps.py — Compute interface width epsilon and required mesh resolution
based on Kaempfer & Plapp (2009).

Key distinction: eps (model parameter) vs actual interface width w = 2*sqrt(2)*eps.
Mesh must resolve the actual interface width, not eps itself.

Usage (CLI)
-----------
  python preprocess/comp_eps.py --Lx 1e-4 --dim 1 [--Rave 3e-5] [--n 4]
  python preprocess/comp_eps.py --Lx 1e-4 --Ly 2.5e-4 --dim 2
  python preprocess/comp_eps.py --Lx 1e-4 --dim 1 --patch inputs/tests/mytest.opts

  --n        elements across the actual interface width (default 4)
  --patch    opts file to update: replaces -eps and -Nx/-Ny/-Nz in-place
  --quiet    print only eps and Nx (one per line) for scripting
"""

import argparse
import math
import re
import sys
from pathlib import Path

import numpy as np

# =========================================================================
# Physical constants (fixed)
# =========================================================================
_ALPHA   = 1.0e-2        # Condensation coefficient
_M_H2O   = 2.99e-26      # Mass of water molecule (kg)
_K_B     = 1.38e-23      # Boltzmann constant (J/K)
_K_I     = 2.29          # Thermal conductivity of ice (W/m/K)
_C_I     = 1.8e6         # Volumetric heat capacity of ice (J/m³/K)
_RHO_ICE = 919.0         # Density of ice (kg/m³)
_DV0     = 2.178e-5      # Vapour diffusion coefficient (m²/s)
_PATM    = 101325.0      # Atmospheric pressure (Pa)
_RHO_AIR = 1.341         # Air density (kg/m³)
_BB      = 0.62

# ASHRAE saturation-pressure polynomial coefficients (same as RhoVS_I in C)
_KJ = [-0.5865e4, 0.2224e2, 0.1375e-1, -0.3403e-4, 0.2697e-7, 0.6918]


def rho_vs_sat(T_C: float) -> float:
    """Saturation vapour density [kg/m³] — matches RhoVS_I() in material_properties.c."""
    T  = T_C + 273.15
    kj = _KJ
    exp_arg = (kj[0]/T + kj[1] + kj[2]*T + kj[3]*T**2
               + kj[4]*T**3 + kj[5]*math.log(T))
    Pvs = math.exp(exp_arg)
    return _RHO_AIR * _BB * Pvs / (_PATM - Pvs)


def compute_params(
    Lx: float,
    Ly: float = 0.0,
    Lz: float = 0.0,
    Rave: float = 3.0e-5,
    T0_C: float = -20.0,
    n_per_interface: int = 4,
) -> dict:
    """
    Compute eps and element counts for a given domain and grain size.

    Parameters
    ----------
    Lx, Ly, Lz       : domain side lengths (m); set to 0 for unused dimensions
    Rave              : representative grain / feature radius (m)
    T0_C              : mean temperature (°C)
    n_per_interface   : target number of elements across the diffuse interface
                        (actual width w = 2*sqrt(2)*eps).  Default 4.

    Returns
    -------
    dict with keys:
      eps       – recommended interface width parameter (m)
      eps_max   – upper bound on eps from Kaempfer & Plapp constraints
      w_actual  – actual diffuse interface width 2*sqrt(2)*eps (m)
      Nx, Ny, Nz – element counts (Ny/Nz = 0 when Ly/Lz = 0)
      dx        – element size (m)
    """
    T0 = T0_C + 273.15

    # Saturation vapour density
    rho_vs  = rho_vs_sat(T0_C)
    rho_rat = rho_vs / _RHO_ICE

    # Kinetic coefficient β₀  [s/m]
    beta_prime = (1.0 / _ALPHA) * math.sqrt(2.0 * math.pi * _M_H2O / (_K_B * T0))
    beta0      = beta_prime / rho_rat

    # Upper bounds on eps (Kaempfer & Plapp 2009, eqs. 44-46)
    eps_heat  = (_K_I / _C_I) * rho_rat * beta0
    eps_vapor = _DV0           * rho_rat * beta0
    eps_geom  = Rave
    eps_max   = min(eps_heat, eps_vapor, eps_geom) / 2.0  # add safety factor to ensure all constraints are satisfied

    # Actual interface width and element size
    w_actual = 2.0 * math.sqrt(2.0) * eps_max
    dx       = w_actual / n_per_interface

    Nx = math.ceil(Lx / dx)          if Lx > 0 else 0
    Ny = math.ceil(Ly / dx)          if Ly > 0 else 0
    Nz = math.ceil(Lz / dx)          if Lz > 0 else 0

    return dict(
        eps=eps_max, eps_max=eps_max,
        w_actual=w_actual, dx=dx,
        Nx=Nx, Ny=Ny, Nz=Nz,
        rho_vs=rho_vs, rho_rat=rho_rat,
        eps_heat=eps_heat, eps_vapor=eps_vapor, eps_geom=eps_geom,
    )


def patch_opts(opts_path: Path, params: dict, dim: int) -> None:
    """
    Update -eps, -Nx (and -Ny, -Nz for dim ≥ 2) in an opts file in-place.
    Adds the option if it is not already present.
    """
    text = opts_path.read_text()

    def _replace_or_append(text: str, key: str, value: str) -> str:
        pattern = rf"^({re.escape(key)}\s+)\S+(\s*)$"
        if re.search(pattern, text, flags=re.MULTILINE):
            return re.sub(pattern, rf"{key} {value}\2", text, flags=re.MULTILINE)
        # append before the first blank line or at the end
        return text.rstrip() + f"\n{key} {value}\n"

    text = _replace_or_append(text, "-eps", f"{params['eps']:.4e}")
    text = _replace_or_append(text, "-Nx",  str(params["Nx"]))
    if dim >= 2 and params["Ny"] > 0:
        text = _replace_or_append(text, "-Ny", str(params["Ny"]))
    if dim == 3 and params["Nz"] > 0:
        text = _replace_or_append(text, "-Nz", str(params["Nz"]))

    opts_path.write_text(text)
    print(f"  Patched {opts_path}")


# =========================================================================
# CLI
# =========================================================================
def _cli():
    ap = argparse.ArgumentParser(
        description="Compute eps and Nx/Ny/Nz from domain geometry.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--Lx",   type=float, required=True,  help="Domain length x (m)")
    ap.add_argument("--Ly",   type=float, default=0.0,    help="Domain length y (m); 0 = unused")
    ap.add_argument("--Lz",   type=float, default=0.0,    help="Domain length z (m); 0 = unused")
    ap.add_argument("--Rave", type=float, default=3.0e-5, help="Representative grain radius (m)")
    ap.add_argument("--T0",   type=float, default=-20.0,  help="Mean temperature (°C)")
    ap.add_argument("--n",    type=int,   default=4,
                    help="Elements across interface width (higher = finer)")
    ap.add_argument("--dim",  type=int,   default=None,
                    help="Problem dimension (1/2/3); inferred from Ly/Lz if omitted")
    ap.add_argument("--patch", type=Path, default=None,
                    help="Opts file to update with computed eps and Nx/Ny/Nz")
    ap.add_argument("--quiet", action="store_true",
                    help="Print only eps and Nx (one per line) for scripting")
    args = ap.parse_args()

    # infer dimension
    dim = args.dim
    if dim is None:
        if args.Lz > 0:   dim = 3
        elif args.Ly > 0: dim = 2
        else:             dim = 1

    p = compute_params(
        Lx=args.Lx, Ly=args.Ly, Lz=args.Lz,
        Rave=args.Rave, T0_C=args.T0,
        n_per_interface=args.n,
    )

    if args.quiet:
        print(f"{p['eps']:.4e}")
        print(p["Nx"])
        if dim >= 2: print(p["Ny"])
        if dim == 3: print(p["Nz"])
        if args.patch:
            patch_opts(args.patch, p, dim)
        return

    print("=" * 65)
    print("  Interface width and mesh sizing — Kaempfer & Plapp (2009)")
    print("=" * 65)
    T0 = args.T0
    print(f"\n--- Thermodynamic state at T0 = {T0+273.15:.2f} K ({T0:.1f} °C) ---")
    print(f"  rho_vs                 = {p['rho_vs']:.4e} kg/m³")
    print(f"  rho_vs / rho_ice       = {p['rho_rat']:.4e}")

    print(f"\n--- Epsilon upper bounds (eps must be LESS than these) ---")
    print(f"  From heat diffusion    = {p['eps_heat']:.4e} m")
    print(f"  From vapour diffusion  = {p['eps_vapor']:.4e} m")
    print(f"  From grain geometry    = {p['eps_geom']:.4e} m  (Rave)")
    print(f"  Maximum allowable eps  = {p['eps_max']:.4e} m")
    print(f"  Actual interface width = {p['w_actual']:.4e} m  (2√2·eps)")
    print(f"  eps / Rave             = {p['eps_max']/args.Rave:.4f}")

    print(f"\n--- Mesh: resolving interface with --n = {args.n} elements ---")
    print(f"  dx  = {p['dx']:.4e} m")
    print(f"  Nx  = {p['Nx']}", end="")
    if dim >= 2: print(f",  Ny = {p['Ny']}", end="")
    if dim == 3: print(f",  Nz = {p['Nz']}", end="")
    print()

    print(f"\n--- Recommended opts entries ---")
    print(f"  -eps  {p['eps']:.4e}")
    print(f"  -Nx   {p['Nx']}")
    if dim >= 2: print(f"  -Ny   {p['Ny']}")
    if dim == 3: print(f"  -Nz   {p['Nz']}")

    if args.patch:
        print()
        patch_opts(args.patch, p, dim)

    print("=" * 65)


if __name__ == "__main__":
    _cli()
