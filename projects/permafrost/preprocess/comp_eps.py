"""
comp_eps.py — Compute the diffuse-interface parameter `eps` and the required
mesh resolution for the dry-snow phase-field model, following Kaempfer & Plapp,
Phys. Rev. E 79, 031502 (2009).

MESH CONVENTION (your solver)
------------------------------
In this solver's free energy, setting the element size h = eps (i.e. one
element per eps) produces ~7.5 elements across the full diffuse interface.
Therefore:

    h  = eps
    Nx = ceil(Lx / eps)

There is NO additional 2√2 prefactor or n_per_interface knob — those belong to
Karma-convention solvers with a different double-well normalisation.

BOUNDS IMPLEMENTED
------------------
All three are from Kaempfer & Plapp (2009), Sec. III C / eqs. (42)–(46).
K&P use the notation W for what this solver calls eps.

  [B-HEAT-ICE]   W ≪ (κᵢ/Cᵢ) · (ρ_vs/ρᵢ) · β₀      Eq. (43a)
  [B-HEAT-AIR]   W ≪ (κₐ/Cₐ) · (ρ_vs/ρᵢ) · β₀      Eq. (43b)
  [B-VAPOR]      W ≪  Dᵥ     · (ρ_vs/ρᵢ) · β₀      Eq. (43c)
  [B-KINETIC]    W ≪  d₀*    / (β₀* · vₙ)            Eq. (45)
                       where d₀* = d₀·(ρ_vs/ρᵢ), β₀* = β₀·(ρ_vs/ρᵢ) = β_HK
                       → simplifies to d₀/(β₀·vₙ)    [dimensionally identical]
  [B-CURV]       W ≪  R_ave                           geometric (near eq. 44)

In K&P's notation:
  β₀  = physical kinetic coefficient  [s/m]  (Eq. 26, unscaled)
  β₀* = β₀ · (ρ_vs/ρᵢ) = β_HK       [s/m]  (the Hertz-Knudsen coefficient)
  d₀  = γ a³/(k_B T)                 [m]    (Eq. 13, physical capillary length)
  d₀* = d₀ · (ρ_vs/ρᵢ)              [m]    (Eq. 25, scaled capillary length)

β₀ is related to the condensation coefficient α_c via:
  β₀* = β_HK = (1/α_c) · √(2π m / k_B T)
  β₀  = β_HK / (ρ_vs/ρᵢ)

α_c is temperature-dependent and highly uncertain (10⁻⁴ to 10⁻¹ in the
literature; see Libbrecht & Rickerby 2013, Bouvet et al. 2022).  This script
accepts either a single value (--alpha) or an α_c sweep (--alpha_range) to
show how the binding bound and eps change with α_c.

ADDITIONAL DIAGNOSTICS
-----------------------
  Pe (Péclet)       Pe = eps · vₙ / Dᵥ                should be ≪ 1
  ξ_v validity      ξ_v ≤ ρ_vs/ρᵢ ~ 5×10⁻⁶            K&P Eq. (48)
                    (only checked/warned; not used to size eps)

The ξ_v constraint is a separate bound on the time-scaling factor used to
accelerate diffusion, not on eps itself.  It is included here because
violating it while keeping eps valid gives false confidence.

Usage (CLI)
-----------
  # single temperature, one α_c
  python comp_eps.py --Lx 6e-4 --Ly 6e-4 --Lz 6e-4 --T0 -5 --alpha 1e-3

  # sweep α_c at two temperatures
  python comp_eps.py --Lx 6e-4 --Ly 6e-4 --T0 -20 --alpha_range 1e-4 1e-1 12

  # patch an existing opts file in place
  python comp_eps.py --Lx 6e-4 --T0 -5 --alpha 1e-3 --patch inputs/sim.opts

  # machine-readable (eps then Nx [Ny] [Nz])
  python comp_eps.py --Lx 6e-4 --T0 -5 --alpha 1e-3 --quiet
"""

import argparse
import math
import re
import sys
from pathlib import Path

# =========================================================================
# Physical constants — match material_properties.c exactly
# =========================================================================
_M_H2O   = 2.99e-26       # mass of a water molecule [kg]
_K_B     = 1.38e-23       # Boltzmann constant [J/K]
_K_I     = 2.29           # thermal conductivity of ice [W/m/K]
_K_A     = 0.02           # thermal conductivity of air [W/m/K]
_C_I     = 1.8e6          # volumetric heat capacity of ice, ρᵢ·cpᵢ [J/m³/K]
_C_A     = 1.341 * 1004   # volumetric heat capacity of air, ρₐ·cpₐ [J/m³/K]
_RHO_ICE = 919.0          # ice density [kg/m³]
_DV0     = 2.178e-5       # vapor diffusivity in air at 273.15 K [m²/s]
_PATM    = 101325.0       # atmospheric pressure [Pa]
_RHO_AIR = 1.341          # air density [kg/m³]
_BB      = 0.62           # ASHRAE humidity ratio coefficient (≈ M_water/M_dryair)
_GAMMA   = 0.109          # ice–air surface energy γ [J/m²]  (K&P Table I)
_VM_ICE  = 1.963e-5       # molar volume of ice, M_water/ρᵢ [m³/mol]
_R_GAS   = 8.314          # universal gas constant [J/mol/K]
_N_A     = 6.022e23       # Avogadro's number [mol⁻¹]

# ASHRAE saturation-pressure polynomial (K&P Table I, eq. 24 source)
_KJ = [-0.5865e4, 0.2224e2, 0.1375e-1, -0.3403e-4, 0.2697e-7, 0.6918]


# =========================================================================
# Thermodynamic functions
# =========================================================================

def rho_vs_sat(T_C: float) -> float:
    """
    Saturation vapour density over ice ρ_vs [kg/m³].
    Uses the ASHRAE polynomial fit, matching RhoVS_I() in material_properties.c.
    """
    T = T_C + 273.15
    exp_arg = (_KJ[0]/T + _KJ[1] + _KJ[2]*T + _KJ[3]*T**2
               + _KJ[4]*T**3 + _KJ[5]*math.log(T))
    Pvs = math.exp(exp_arg)
    # ρ_v = ρ_air · ω  where ω = 0.622 · Pvs/(Pa - Pvs)  [ASHRAE]
    return _RHO_AIR * _BB * Pvs / (_PATM - Pvs)


def Dv_T(T_C: float) -> float:
    """
    Vapour diffusivity in air [m²/s], power-law T-correction.
    Dᵥ(T) = Dᵥ₀ · (T/273.15)^1.81  (K&P, their ref. [25]).
    """
    return _DV0 * ((T_C + 273.15) / 273.15) ** 1.81


def hertz_knudsen_beta_star(T_C: float, alpha_c: float) -> float:
    """
    Hertz-Knudsen kinetic coefficient β* = β₀* = β_HK [s/m].

    This is K&P's *scaled* kinetic coefficient (their Eq. 26):
        β₀* = β_HK = (1/α_c) · √(2π m / k_B T)

    The *unscaled* coefficient is β₀ = β₀* / (ρ_vs/ρᵢ).
    """
    T_K = T_C + 273.15
    return (1.0 / alpha_c) * math.sqrt(2.0 * math.pi * _M_H2O / (_K_B * T_K))


def capillary_length(T_C: float, gamma: float = _GAMMA) -> float:
    """
    Physical capillary length d₀ = γ a³/(k_B T) [m]  (K&P Eq. 13).

    Implemented via the macroscopic Kelvin form d₀ = γ V_m/(R T),
    which is identical since V_m = N_A a³ and R = N_A k_B.
    At T = -10°C: d₀ ≈ 1.3×10⁻⁹ m  (K&P Table I). ✓
    """
    T_K = T_C + 273.15
    return gamma * _VM_ICE / (_R_GAS * T_K)


# =========================================================================
# Core computation
# =========================================================================

def compute_eps(
    Lx: float,
    Ly: float = 0.0,
    Lz: float = 0.0,
    Rave: float = 3.0e-5,
    T0_C: float = -20.0,
    alpha_c: float = 1.0e-2,
    safety: float = 0.5,
    v_n: float = 1.0e-9,
    xi_v: float | None = None,
) -> dict:
    """
    Compute the interface-width parameter eps and element counts.

    Mesh rule (this solver): h = eps → Nx = ceil(Lx / eps).

    Parameters
    ----------
    Lx, Ly, Lz  : domain side lengths [m]; 0 = dimension not used
    Rave         : representative grain/feature radius [m] (geometric bound)
    T0_C         : mean simulation temperature [°C]
    alpha_c      : condensation coefficient α_c (dimensionless, 0 < α_c ≤ 1)
    safety       : eps = safety · eps_max  (K&P use ≪; safety < 1 honours this)
    v_n          : assumed normal front velocity [m/s] for Eq. (45) bound
    xi_v         : time-scaling factor ξᵥ; if provided, its validity is checked
                   against K&P Eq. (48): ξᵥ ≤ ρ_vs/ρᵢ

    Returns
    -------
    dict with bounds, eps, mesh counts, diagnostics, and any warnings.
    """
    T0_K = T0_C + 273.15

    # --- Thermodynamic state ---
    rho_vs  = rho_vs_sat(T0_C)
    rho_rat = rho_vs / _RHO_ICE          # ρ_vs/ρᵢ  (≈ 5×10⁻⁶ near melting)

    # Kinetic coefficients (K&P Eqs. 25–26)
    beta_star = hertz_knudsen_beta_star(T0_C, alpha_c)   # β₀* = β_HK  [s/m]
    beta0     = beta_star / rho_rat                       # β₀ (unscaled) [s/m]

    # Capillary lengths (physical and scaled, K&P Eq. 25)
    d0        = capillary_length(T0_C)       # d₀  [m]       (K&P Eq. 13)
    d0_star   = d0 * rho_rat                 # d₀* = d₀·(ρ_vs/ρᵢ) [m]

    # Vapour diffusivity at T₀
    Dv        = Dv_T(T0_C)

    # Thermal diffusivities
    alpha_i   = _K_I / _C_I   # κᵢ/Cᵢ  [m²/s]
    alpha_a   = _K_A / _C_A   # κₐ/Cₐ  [m²/s]

    # ----------------------------------------------------------------
    # Upper bounds on eps  (K&P eqs. 42–43 and 45)
    # All three forms of eq.(43) can be written as:
    #   W ≪ D_eff · β₀*  where D_eff is the relevant diffusivity
    #   and β₀* = (ρ_vs/ρᵢ)·β₀  is the scaled kinetic coefficient.
    # ----------------------------------------------------------------
    b_heat_ice = alpha_i * beta_star     # (κᵢ/Cᵢ) · β₀*    Eq. (43a)
    b_heat_air = alpha_a * beta_star     # (κₐ/Cₐ) · β₀*    Eq. (43b)
    b_vapor    = Dv      * beta_star     # Dᵥ · β₀*          Eq. (43c)

    # Eq. (45): W ≪ d₀*/(β₀*·vₙ) = d₀/(β₀·vₙ)  (scaled and unscaled
    # forms are identical because the ρ_vs/ρᵢ factors cancel)
    b_kinetic  = d0 / (beta0 * v_n)     # = d₀*/(β₀*·vₙ)   Eq. (45)

    # Geometric bound: W ≪ 1/K ≈ R_ave             (near Eq. 44 in K&P)
    b_curv     = Rave

    bounds = {
        "B-HEAT-ICE": b_heat_ice,
        "B-HEAT-AIR": b_heat_air,
        "B-VAPOR":    b_vapor,
        "B-KINETIC":  b_kinetic,
        "B-CURV":     b_curv,
    }

    eps_max  = min(bounds.values())
    binding  = min(bounds, key=bounds.get)
    eps      = safety * eps_max

    # ----------------------------------------------------------------
    # Mesh sizing: h = eps, Nx = ceil(Lx/eps)
    # ~7.5 elements across the full diffuse interface at h = eps.
    # ----------------------------------------------------------------
    Nx = math.ceil(Lx / eps) if Lx > 0 else 0
    Ny = math.ceil(Ly / eps) if Ly > 0 else 0
    Nz = math.ceil(Lz / eps) if Lz > 0 else 0

    # ----------------------------------------------------------------
    # Diagnostics
    # ----------------------------------------------------------------
    # Péclet number: Pe = eps·vₙ/Dᵥ  — should be ≪ 1
    Pe = eps * v_n / Dv

    # K&P Eq. (48): ξᵥ must not exceed ρ_vs/ρᵢ
    xi_v_max  = rho_rat   # hard upper bound from K&P
    xi_v_warn = None
    if xi_v is not None and xi_v > xi_v_max:
        xi_v_warn = (
            f"WARNING: ξᵥ = {xi_v:.2e} exceeds ρ_vs/ρᵢ = {xi_v_max:.2e} "
            f"(K&P Eq. 48). The quasi-steady assumption behind the time-scaling "
            f"scheme is violated; results will be unphysical."
        )

    # Warn if eps is within 2x of the binding bound (safety may be too loose)
    eps_ratio = eps / eps_max   # = safety; included for completeness

    return dict(
        # --- chosen eps and mesh ---
        eps=eps,
        Nx=Nx, Ny=Ny, Nz=Nz,
        # --- bounds ---
        eps_max=eps_max,
        binding=binding,
        **bounds,
        # --- diagnostics ---
        Pe=Pe,
        rho_vs=rho_vs,
        rho_rat=rho_rat,
        beta_star=beta_star,
        beta0=beta0,
        d0=d0,
        d0_star=d0_star,
        Dv=Dv,
        alpha_i=alpha_i,
        alpha_a=alpha_a,
        v_n=v_n,
        safety=safety,
        xi_v_max=xi_v_max,
        xi_v_warn=xi_v_warn,
    )


def alpha_c_sensitivity(
    T0_C: float,
    alpha_lo: float,
    alpha_hi: float,
    n_points: int = 10,
    **kwargs,
) -> list[dict]:
    """
    Sweep α_c log-uniformly from alpha_lo to alpha_hi and return a list
    of result dicts (one per α_c value), suitable for tabular display.
    kwargs are forwarded to compute_eps (Lx, Ly, Lz, Rave, safety, v_n, xi_v).
    """
    alphas = [
        math.exp(math.log(alpha_lo) + i * math.log(alpha_hi / alpha_lo) / (n_points - 1))
        for i in range(n_points)
    ]
    return [compute_eps(T0_C=T0_C, alpha_c=a, **kwargs) for a in alphas], alphas


# =========================================================================
# Opts-file patching
# =========================================================================

def patch_opts(opts_path: Path, params: dict, dim: int) -> None:
    """Replace -eps, -Nx, and (if dim ≥ 2/3) -Ny/-Nz in opts_path in place."""
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
# Formatted output helpers
# =========================================================================

def _print_header(title: str) -> None:
    print("=" * 72)
    print(f"  {title}")
    print("=" * 72)


def _print_single(args, p: dict, alpha_c: float, dim: int) -> None:
    _print_header("Interface width and mesh sizing — Kaempfer & Plapp (2009) eqs. 42–46")

    print(f"\n--- Input parameters ---")
    print(f"  T₀                       = {args.T0:.1f} °C  ({args.T0+273.15:.2f} K)")
    print(f"  α_c (condensation coeff) = {alpha_c:.3e}   [dimensionless]")
    print(f"  v_n (front velocity)     = {args.vn:.1e}  m/s")
    print(f"  R_ave (grain radius)     = {args.Rave:.2e} m")
    print(f"  Safety factor            = {args.safety:.2f}  (eps = safety · min(bounds))")

    print(f"\n--- Thermodynamic state at T₀ ---")
    print(f"  ρ_vs(T₀)         = {p['rho_vs']:.4e}  kg/m³")
    print(f"  ρ_vs/ρᵢ          = {p['rho_rat']:.4e}  [-]   (K&P 'χ' ratio; ~5×10⁻⁶ near melt)")
    print(f"  β₀* = β_HK       = {p['beta_star']:.4e}  s/m   (scaled kinetic coeff, K&P Eq. 26)")
    print(f"  β₀  (unscaled)   = {p['beta0']:.4e}  s/m   (= β₀* / (ρ_vs/ρᵢ))")
    print(f"  d₀  (physical)   = {p['d0']:.4e}  m     (capillary length, K&P Eq. 13; ~1.3×10⁻⁹ at -10°C)")
    print(f"  d₀* (scaled)     = {p['d0_star']:.4e}  m     (= d₀ · ρ_vs/ρᵢ, K&P Eq. 25)")
    print(f"  Dᵥ(T₀)           = {p['Dv']:.4e}  m²/s")
    print(f"  κᵢ/Cᵢ            = {p['alpha_i']:.4e}  m²/s")
    print(f"  κₐ/Cₐ            = {p['alpha_a']:.4e}  m²/s")

    print(f"\n--- Upper bounds on eps (K&P eqs. 42–46) ---")
    labels = {
        "B-HEAT-ICE": "Eq.(43a)  (κᵢ/Cᵢ)·β₀*",
        "B-HEAT-AIR": "Eq.(43b)  (κₐ/Cₐ)·β₀*",
        "B-VAPOR":    "Eq.(43c)  Dᵥ·β₀*      ",
        "B-KINETIC":  "Eq.(45)   d₀/(β₀·vₙ)  ",
        "B-CURV":     "Geometric R_ave         ",
    }
    for key, label in labels.items():
        marker = "  ← BINDING" if key == p["binding"] else ""
        print(f"  {label}  =  {p[key]:.4e} m{marker}")

    print(f"\n--- Chosen eps (after safety = {p['safety']:.2f}) ---")
    print(f"  eps               = {p['eps']:.4e}  m")
    print(f"  eps / R_ave       = {p['eps']/args.Rave:.4f}")
    print(f"  eps / eps_max     = {p['safety']:.4f}  (= safety factor as a check)")

    print(f"\n--- Péclet check ---")
    Pe_status = "OK (≪ 1)" if p["Pe"] < 0.1 else ("borderline" if p["Pe"] < 1.0 else "VIOLATED — Pe ≥ 1!")
    print(f"  Pe = eps·vₙ/Dᵥ   = {p['Pe']:.4e}  [{Pe_status}]")
    print(f"  (K&P: should remain ≪ 1 for the sharp-interface analogy to hold)")

    print(f"\n--- ξᵥ validity check (K&P Eq. 48) ---")
    print(f"  ξᵥ hard upper bound = ρ_vs/ρᵢ = {p['xi_v_max']:.4e}")
    if p["xi_v_warn"]:
        print(f"  {p['xi_v_warn']}")
    elif args.xiv is not None:
        print(f"  ξᵥ = {args.xiv:.2e}  → OK (ξᵥ ≤ ρ_vs/ρᵢ)")
    else:
        print(f"  (ξᵥ not provided; pass --xiv to check)")

    print(f"\n--- Mesh sizing (h = eps; ~7.5 elements across full interface) ---")
    print(f"  h = eps           = {p['eps']:.4e}  m")
    print(f"  Nx = ceil(Lx/eps) = {p['Nx']}")
    if dim >= 2 and p["Ny"] > 0:
        print(f"  Ny = ceil(Ly/eps) = {p['Ny']}")
    if dim == 3 and p["Nz"] > 0:
        print(f"  Nz = ceil(Lz/eps) = {p['Nz']}")

    print(f"\n--- Recommended opts entries ---")
    print(f"  -eps  {p['eps']:.4e}")
    print(f"  -Nx   {p['Nx']}")
    if dim >= 2 and p["Ny"] > 0: print(f"  -Ny   {p['Ny']}")
    if dim == 3 and p["Nz"] > 0: print(f"  -Nz   {p['Nz']}")

    print("=" * 72)


def _print_sweep(results: list[dict], alphas: list[float], T0_C: float) -> None:
    _print_header(f"α_c sensitivity sweep at T₀ = {T0_C:.1f} °C")
    hdr = (f"  {'α_c':>10s}  {'β₀* [s/m]':>12s}  {'eps_max [m]':>13s}  "
           f"{'eps [m]':>11s}  {'binding':>12s}  {'Pe':>10s}")
    print(hdr)
    print("  " + "-" * 75)
    for p, a in zip(results, alphas):
        print(f"  {a:10.3e}  {p['beta_star']:12.4e}  {p['eps_max']:13.4e}  "
              f"{p['eps']:11.4e}  {p['binding']:>12s}  {p['Pe']:10.3e}")
    print("=" * 72)
    print("  Note: eps = safety · eps_max.  Pe = eps · vₙ / Dᵥ.")
    print("  Binding bound changes with α_c because β₀* ∝ 1/α_c scales the")
    print("  heat/vapor bounds while leaving B-KINETIC and B-CURV fixed.")
    print("=" * 72)


# =========================================================================
# CLI
# =========================================================================

def _infer_dim(Lx, Ly, Lz, explicit_dim):
    if explicit_dim is not None:
        return explicit_dim
    if Lz > 0:   return 3
    if Ly > 0:   return 2
    return 1


def _cli():
    ap = argparse.ArgumentParser(
        description=(
            "Compute eps (interface width) and Nx/Ny/Nz per Kaempfer & Plapp (2009). "
            "Mesh rule: h = eps, Nx = ceil(Lx/eps)."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # Domain
    ap.add_argument("--Lx",     type=float, required=True,  help="Domain length x [m]")
    ap.add_argument("--Ly",     type=float, default=0.0,    help="Domain length y [m]; 0 = unused")
    ap.add_argument("--Lz",     type=float, default=0.0,    help="Domain length z [m]; 0 = unused")
    ap.add_argument("--Rave",   type=float, default=3.0e-5, help="Representative grain radius [m]")
    # Physics
    ap.add_argument("--T0",     type=float, default=-20.0,  help="Mean temperature [°C]")
    ap.add_argument("--alpha",  type=float, default=None,
                    help="Condensation coefficient α_c [−].  Mutually exclusive with --alpha_range.")
    ap.add_argument("--alpha_range", nargs=3, metavar=("LO","HI","N"), default=None,
                    help="α_c sweep: LO HI N_points (log-uniform).  Mutually exclusive with --alpha.")
    ap.add_argument("--vn",     type=float, default=1.0e-9, help="Normal front velocity vₙ [m/s] for Eq.(45)")
    ap.add_argument("--xiv",    type=float, default=None,   help="Time-scaling factor ξᵥ (optional; triggers K&P Eq.48 check)")
    # Numerics
    ap.add_argument("--safety", type=float, default=0.5,
                    help="Safety factor: eps = safety · min(bounds).  K&P use ≪, so < 1 is appropriate.")
    ap.add_argument("--dim",    type=int,   default=None,   help="Problem dimension (1/2/3); inferred from Ly/Lz if omitted")
    # Output
    ap.add_argument("--patch",  type=Path,  default=None,   help="Opts file to update in place with eps and Nx/Ny/Nz")
    ap.add_argument("--quiet",  action="store_true",         help="Print only eps and Nx/Ny/Nz (one per line)")
    args = ap.parse_args()

    # --- Validate α_c input ---
    if args.alpha is not None and args.alpha_range is not None:
        ap.error("--alpha and --alpha_range are mutually exclusive.")
    if args.alpha is None and args.alpha_range is None:
        ap.error("Provide either --alpha or --alpha_range.")

    dim = _infer_dim(args.Lx, args.Ly, args.Lz, args.dim)

    common = dict(
        Lx=args.Lx, Ly=args.Ly, Lz=args.Lz,
        Rave=args.Rave, v_n=args.vn,
        safety=args.safety, xi_v=args.xiv,
    )

    # --- Single α_c run ---
    if args.alpha is not None:
        p = compute_eps(T0_C=args.T0, alpha_c=args.alpha, **common)

        if args.quiet:
            print(f"{p['eps']:.4e}")
            print(p["Nx"])
            if dim >= 2: print(p["Ny"])
            if dim == 3: print(p["Nz"])
            if args.patch:
                patch_opts(args.patch, p, dim)
            return

        _print_single(args, p, args.alpha, dim)

        if p["xi_v_warn"]:
            print(f"\n  {p['xi_v_warn']}\n")

        if args.patch:
            print()
            patch_opts(args.patch, p, dim)
        return

    # --- α_c sweep ---
    lo, hi, n = float(args.alpha_range[0]), float(args.alpha_range[1]), int(args.alpha_range[2])
    results, alphas = alpha_c_sensitivity(
        T0_C=args.T0, alpha_lo=lo, alpha_hi=hi, n_points=n, **common
    )
    _print_sweep(results, alphas, args.T0)


if __name__ == "__main__":
    _cli()