"""
comp_eps.py — Compute the diffuse-interface parameter `eps` and the required
mesh resolution, following Kaempfer & Plapp, Phys. Rev. E 79, 031502 (2009).

The bounds implemented are the paper's equations (42), (43), and (45) — the
sharp-interface limit applicability conditions — plus a geometric constraint
`W ≪ 1/K ≈ R_ave` from the same section.

  Eq (42), using the scaled kinetic coefficient β₀′ = (ρ_vs/ρ_i) · β₀ :
      W / (λ_i / C_i)  ≪  β₀′
      W / (λ_a / C_a)  ≪  β₀′
      W /  D_v         ≪  β₀′
  Eq (43) = Eq (42) rewritten with the unscaled β₀ pulled out.
  Eq (45)                        :  W  ≪  d₀ / (β₀ · v_n)
  Geometric (eq. text near 44)   :  W  ≪  1 / K  ≈  R_ave

The script's variable `eps` is identified with the paper's `W` (Karma-style
symmetric well, φ ∈ [-1,+1]). Our solver uses the φ ∈ [0,1] well
`Ψ = η φ²(1-φ)²`, which is the same well under the affine map ψ = 2φ-1, so
the interface-width parameter has the same meaning.

The "actual interface width" `w = 2√2·eps` follows the Karma-Plapp convention
(integrated surface-tension definition). The visible 5%-95% band of our φ ∈
[0,1] equilibrium tanh profile is ≈ 6·eps — wider than 2√2·eps ≈ 2.83·eps.
For mesh sizing we use the 2√2·eps convention with `--n_per_interface`
elements across it; n=4 puts about 10 elements across the 1%-99% band, which
matches the visual rule of thumb of "~7.5 elements across phi=0.01-0.99"
once a safety factor is applied.

The paper uses `≪` everywhere, not `≤`. Use `--safety` to apply a margin
below the binding bound (default 0.5; lower = stricter).

Usage (CLI)
-----------
  python preprocess/comp_eps.py --Lx 1e-4 --dim 1 [--Rave 3e-5] [--n 4]
  python preprocess/comp_eps.py --Lx 1e-4 --Ly 2.5e-4 --dim 2
  python preprocess/comp_eps.py --Lx 1e-4 --dim 1 --patch inputs/geometry/foo.opts

  --safety       multiplicative safety margin on eps_max (default 0.5)
  --vn           assumed front velocity for Eq 45 (default 1e-9 m/s)
  --n            elements across w_actual = 2√2·eps (default 4)
  --patch        opts file to update in place (replaces -eps, -Nx/-Ny/-Nz)
  --quiet        print only eps and Nx/Ny/Nz for scripting
"""

import argparse
import math
import re
import sys
from pathlib import Path

# =========================================================================
# Physical constants (fixed; consistent with material_properties.c)
# =========================================================================
_ALPHA    = 1.0e-2        # Condensation coefficient α  [dimensionless]
_M_H2O    = 2.99e-26      # Mass of a water molecule    [kg]
_K_B      = 1.38e-23      # Boltzmann constant          [J/K]
_K_I      = 2.29          # Thermal conductivity of ice [W/m/K]
_K_A      = 0.02          # Thermal conductivity of (dry) air [W/m/K]
_C_I      = 1.8e6         # Volumetric heat capacity of ice [J/m³/K]
_C_A      = 1.341 * 1004  # Volumetric heat capacity of air ≈ rho_air * cp_air [J/m³/K]
_RHO_ICE  = 919.0         # Ice density [kg/m³]
_DV0      = 2.178e-5      # Vapor diffusivity in air @ 273.15 K [m²/s]
_PATM     = 101325.0      # Atmospheric pressure [Pa]
_RHO_AIR  = 1.341         # Air density [kg/m³]
_BB       = 0.62          # ASHRAE coefficient (mass ratio)
_GAMMA_IV = 0.109         # Ice-vapor surface energy [J/m²]
_VM_ICE   = 1.96e-5       # Molar volume of ice [m³/mol] (= M / ρ)
_R_GAS    = 8.314         # Gas constant [J/mol/K]

# ASHRAE saturation-pressure polynomial coefficients (matches RhoVS_I() in C)
_KJ = [-0.5865e4, 0.2224e2, 0.1375e-1, -0.3403e-4, 0.2697e-7, 0.6918]


def rho_vs_sat(T_C: float) -> float:
    """Saturation vapor density [kg/m³] — matches RhoVS_I() in material_properties.c."""
    T  = T_C + 273.15
    exp_arg = (_KJ[0]/T + _KJ[1] + _KJ[2]*T + _KJ[3]*T**2
               + _KJ[4]*T**3 + _KJ[5]*math.log(T))
    Pvs = math.exp(exp_arg)
    return _RHO_AIR * _BB * Pvs / (_PATM - Pvs)


def hertz_knudsen_beta(T_C: float, alpha: float = _ALPHA) -> float:
    """
    Hertz-Knudsen kinetic coefficient β_HK = (1/α)·√(2π m / k_B T) [s/m].
    In the paper's notation this equals β₀' (the scaled kinetic coefficient).
    The unscaled β₀ = β_HK / (ρ_vs/ρ_i).
    """
    T_K = T_C + 273.15
    return (1.0 / alpha) * math.sqrt(2.0 * math.pi * _M_H2O / (_K_B * T_K))


def capillary_length_d0(T_C: float, gamma: float = _GAMMA_IV) -> float:
    """
    Capillary length d₀ = γ · v_m / (R_gas · T)  [m], Kelvin/Gibbs-Thomson.
    """
    T_K = T_C + 273.15
    return gamma * _VM_ICE / (_R_GAS * T_K)


def compute_params(
    Lx: float,
    Ly: float = 0.0,
    Lz: float = 0.0,
    Rave: float = 3.0e-5,
    T0_C: float = -20.0,
    n_per_interface: int = 4,
    safety: float = 0.5,
    v_n: float = 1.0e-9,
) -> dict:
    """
    Compute the interface-width parameter eps and element counts using the
    Kaempfer & Plapp (2009) bounds (eqs 42, 43, 45) plus the geometric
    bound R_ave.

    Parameters
    ----------
    Lx, Ly, Lz       : domain side lengths (m); set to 0 for unused dimensions
    Rave             : representative grain / feature radius (m)  — geometric bound
    T0_C             : mean temperature (°C)
    n_per_interface  : target elements across the Karma-convention interface
                       width w_actual = 2√2·eps (default 4)
    safety           : multiplicative margin (eps = safety · min(bounds));
                       paper uses `≪`, so a value < 1 is appropriate
    v_n              : assumed normal front velocity for Eq 45 (default 1e-9 m/s)

    Returns
    -------
    dict with the five eps bounds, the binding bound, eps (post-safety),
    w_actual, dx, Nx/Ny/Nz, and the intermediate thermodynamic quantities.
    """
    T0_K = T0_C + 273.15

    # Thermodynamic state
    rho_vs  = rho_vs_sat(T0_C)
    rho_rat = rho_vs / _RHO_ICE                       # paper: ρ_vs/ρ_i

    # Kinetic coefficients (paper notation)
    beta0p  = hertz_knudsen_beta(T0_C)                # β₀′ (scaled) = β_HK
    beta0   = beta0p / rho_rat                        # β₀ (unscaled)

    # Capillary length
    d0      = capillary_length_d0(T0_C)

    # --- Eq (42) bounds — three thermal/vapor constraints using β₀′ -----
    alpha_T_ice = _K_I / _C_I                         # λ_i / C_i  [m²/s]
    alpha_T_air = _K_A / _C_A                         # λ_a / C_a  [m²/s]
    Dv          = _DV0                                # vapor diffusivity (T₀-corrected below)

    # Optional T-correction for D_v ~ (T/273.15)^1.81 (matches VaporDiffus()):
    Dv_T = _DV0 * (T0_K / 273.15) ** 1.81

    eps_heat_ice  = alpha_T_ice * beta0p              # Eq 42a / 43a
    eps_heat_air  = alpha_T_air * beta0p              # Eq 42b / 43b
    eps_vapor     = Dv_T        * beta0p              # Eq 42c / 43c

    # --- Eq (45) bound — kinetic / front-velocity ---------------------
    # W ≪ d₀ / (β₀ · v_n)
    eps_kinetic   = d0 / (beta0 * v_n)

    # --- Geometric bound — W ≪ 1/K  ≈  R_ave -------------------------
    eps_geom      = Rave

    # Binding bound and safety-applied eps
    bounds = {
        "eps_heat_ice": eps_heat_ice,
        "eps_heat_air": eps_heat_air,
        "eps_vapor":    eps_vapor,
        "eps_kinetic":  eps_kinetic,
        "eps_geom":     eps_geom,
    }
    eps_max  = min(bounds.values())
    eps      = safety * eps_max

    # Karma-convention "interface width" + mesh sizing
    w_actual = 2.0 * math.sqrt(2.0) * eps
    dx       = w_actual / n_per_interface

    Nx = math.ceil(Lx / dx) if Lx > 0 else 0
    Ny = math.ceil(Ly / dx) if Ly > 0 else 0
    Nz = math.ceil(Lz / dx) if Lz > 0 else 0

    return dict(
        eps=eps, eps_max=eps_max, safety=safety,
        w_actual=w_actual, dx=dx,
        Nx=Nx, Ny=Ny, Nz=Nz,
        rho_vs=rho_vs, rho_rat=rho_rat,
        beta0=beta0, beta0p=beta0p, d0=d0, v_n=v_n,
        Dv_T=Dv_T,
        binding=min(bounds, key=bounds.get),
        **bounds,
    )


def patch_opts(opts_path: Path, params: dict, dim: int) -> None:
    """Replace -eps, -Nx, and (if dim≥2/3) -Ny/-Nz in `opts_path` in place."""
    text = opts_path.read_text()

    def _replace_or_append(text: str, key: str, value: str) -> str:
        pattern = rf"^({re.escape(key)}\s+)\S+(\s*)$"
        if re.search(pattern, text, flags=re.MULTILINE):
            return re.sub(pattern, rf"{key} {value}\2", text, flags=re.MULTILINE)
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
        description="Compute eps and Nx/Ny/Nz per Kaempfer & Plapp 2009 (eqs 42-46).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--Lx",     type=float, required=True,  help="Domain length x [m]")
    ap.add_argument("--Ly",     type=float, default=0.0,    help="Domain length y [m]; 0=unused")
    ap.add_argument("--Lz",     type=float, default=0.0,    help="Domain length z [m]; 0=unused")
    ap.add_argument("--Rave",   type=float, default=3.0e-5, help="Representative grain radius [m]")
    ap.add_argument("--T0",     type=float, default=-20.0,  help="Mean temperature [°C]")
    ap.add_argument("--n",      type=int,   default=4,
                    help="Elements across w_actual = 2√2·eps")
    ap.add_argument("--safety", type=float, default=0.5,
                    help="Safety factor on eps (eps = safety · min(bounds)). "
                         "Paper uses ≪, so safety < 1 honours the inequality.")
    ap.add_argument("--vn",     type=float, default=1.0e-9,
                    help="Assumed normal front velocity for Eq 45 [m/s]")
    ap.add_argument("--dim",    type=int,   default=None,
                    help="Problem dimension (1/2/3); inferred from Ly/Lz if omitted")
    ap.add_argument("--patch",  type=Path,  default=None,
                    help="Opts file to update with computed eps and Nx/Ny/Nz")
    ap.add_argument("--quiet",  action="store_true",
                    help="Machine-readable output: print eps and Nx/Ny/Nz, one per line")
    args = ap.parse_args()

    dim = args.dim
    if dim is None:
        if args.Lz > 0:   dim = 3
        elif args.Ly > 0: dim = 2
        else:             dim = 1

    p = compute_params(
        Lx=args.Lx, Ly=args.Ly, Lz=args.Lz,
        Rave=args.Rave, T0_C=args.T0,
        n_per_interface=args.n,
        safety=args.safety, v_n=args.vn,
    )

    if args.quiet:
        print(f"{p['eps']:.4e}")
        print(p["Nx"])
        if dim >= 2: print(p["Ny"])
        if dim == 3: print(p["Nz"])
        if args.patch:
            patch_opts(args.patch, p, dim)
        return

    print("=" * 70)
    print("  Interface width and mesh sizing — Kaempfer & Plapp (2009) eqs 42-46")
    print("=" * 70)

    print(f"\n--- Thermodynamic state at T₀ = {args.T0:.1f} °C ({args.T0+273.15:.2f} K) ---")
    print(f"  ρ_vs(T₀)                = {p['rho_vs']:.4e} kg/m³")
    print(f"  ρ_vs / ρ_i              = {p['rho_rat']:.4e}")
    print(f"  β₀′ = β_HK              = {p['beta0p']:.4e} s/m   (paper's scaled β₀′)")
    print(f"  β₀  = β_HK / (ρ_vs/ρ_i) = {p['beta0']:.4e} s/m   (paper's unscaled β₀)")
    print(f"  d₀  = γ·vₘ/(R·T)        = {p['d0']:.4e} m       (Kelvin capillary length)")
    print(f"  D_v(T₀) = D_v0·(T/273.15)^1.81 = {p['Dv_T']:.4e} m²/s")

    print(f"\n--- Upper bounds on eps  (Kaempfer & Plapp 2009) ---")
    print(f"  Eq 42a — heat in ice:    (λ_i/C_i)·β₀′           = {p['eps_heat_ice']:.4e} m")
    print(f"  Eq 42b — heat in air:    (λ_a/C_a)·β₀′           = {p['eps_heat_air']:.4e} m")
    print(f"  Eq 42c — vapor:           D_v·β₀′                = {p['eps_vapor']:.4e} m")
    print(f"  Eq 45  — kinetic vₙ={args.vn:.1e}:  d₀/(β₀·vₙ)   = {p['eps_kinetic']:.4e} m")
    print(f"  Geometric — R_ave:                                  {p['eps_geom']:.4e} m")
    print(f"  → binding bound is {p['binding']:<14s}        = {p['eps_max']:.4e} m")

    print(f"\n--- Chosen eps after safety factor {p['safety']:.3g} ---")
    print(f"  eps                     = {p['eps']:.4e} m")
    print(f"  w_actual = 2√2·eps      = {p['w_actual']:.4e} m   (Karma convention)")
    print(f"  eps / R_ave             = {p['eps']/args.Rave:.4f}")

    print(f"\n--- Mesh: --n = {args.n} elements across w_actual ---")
    print(f"  dx                      = {p['dx']:.4e} m")
    print(f"  Nx                      = {p['Nx']}")
    if dim >= 2: print(f"  Ny                      = {p['Ny']}")
    if dim == 3: print(f"  Nz                      = {p['Nz']}")

    print(f"\n--- Recommended opts entries ---")
    print(f"  -eps  {p['eps']:.4e}")
    print(f"  -Nx   {p['Nx']}")
    if dim >= 2: print(f"  -Ny   {p['Ny']}")
    if dim == 3: print(f"  -Nz   {p['Nz']}")

    if args.patch:
        print()
        patch_opts(args.patch, p, dim)

    print("=" * 70)


if __name__ == "__main__":
    _cli()
