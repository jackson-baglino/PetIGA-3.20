"""
comp_eps.py — Compute the diffuse-interface parameter `eps`, mesh resolution,
and all derived phase-field model parameters for the dry-snow phase-field model,
following Kaempfer & Plapp, Phys. Rev. E 79, 031502 (2009), and Moure & Fu,
Cryst. Growth Des. 24, 5687 (2024) — sublimation reduction.

MESH CONVENTION (this solver)
------------------------------
Setting h = ε/√2 gives ~7.5 elements across the phi=0.05–0.95 visible band:
    h  = eps / sqrt(2)
    Nx = ceil(Lx / (eps / sqrt(2)))  =  ceil(sqrt(2) * Lx / eps)

This is sqrt(2) MORE elements per direction than the old h=eps rule.

The logistic equilibrium profile phi(x) = 1/(1+exp(-x/eps)) has these
band widths expressed in units of eps:
    phi=0.050–0.950  (5%–95%):      5.89·eps  →  8.33 elements at h=ε/√2  ≈ ~7.5
    phi=0.010–0.990  (1%–99%):      9.19·eps  → 13.0  elements at h=ε/√2
    phi=0.005–0.995  (0.5%–99.5%): 10.59·eps  → 14.97 elements at h=ε/√2
No 2√2 prefactor in the physics — that belongs to Karma-convention solvers.
The 1/√2 here is purely a mesh-resolution choice targeting ~7.5 visible
(phi=0.05–0.95) interface elements.

BOUNDS ON ε (K&P eqs. 42–46, equivalent to M&F SI Cond. 1–3)
------------------------------------------------------------
  [B-HEAT-ICE]   W ≪ (κᵢ/Cᵢ) · (ρ_vs/ρᵢ) · β₀     Eq. (43a)
  [B-HEAT-AIR]   W ≪ (κₐ/Cₐ) · (ρ_vs/ρᵢ) · β₀     Eq. (43b)
  [B-VAPOR]      W ≪  Dᵥ     · (ρ_vs/ρᵢ) · β₀     Eq. (43c)
  [B-KINETIC]    W ≪  d₀ / (β₀ · vₙ)               Eq. (45)
  [B-CURV]       W ≪  R_ave                         geometric (near Eq. 44)

Here β₀ is the PHYSICAL (unscaled) kinetic coefficient (M&F's β_sub, K&P's β₀).

DERIVED PHASE-FIELD MODEL PARAMETERS (M&F SI Eq. 9, verbatim)
-------------------------------------------------------------
Once ε, β_sub, d_sub are fixed, ALL of these are algebraically determined:

   d_sub · (ρ_vs/ρᵢ)  =  a₁ · ε / λ_sub
   β_sub · (ρ_vs/ρᵢ)  =  a₁ · [τ_sub/(ε·λ_sub) − a₂·ε/D*_ia − a₂·ε/Dᵥ]
   M_sub = ε / (3·τ_sub)     α_sub = λ_sub / τ_sub

where D*_ia = (Dᵢ + Dₐ)/2, and the M&F-specific asymptotic constants
(from their SI footnote 1, derived from Karma-Rappel applied to their
free energy):  a₁ ≈ 5,  a₂ ≈ 0.1581.
(These differ from Karma-Rappel's a₁=5√2/8, a₂=0.6267, which are specific
to Karma's ±1 double-well; M&F's φ∈[0,1] well gives different values.)

β_sub INTERPRETATION (important!)
---------------------------------
Both K&P and M&F use β_sub / β₀ to mean the *unscaled* physical kinetic
coefficient in the Gibbs-Thomson condition ρᵥ − ρvs = β_sub · ρvs · vₙ + …
Hertz-Knudsen gives β_sub = β_HK = (1/α_c)·√(2πm/kT).
(K&P sometimes uses β' or β₀' for the *scaled* form β₀·(ρ_vs/ρᵢ) — we
avoid that notation and use β_sub throughout.)

α_c IS THE ONLY GENUINE FREE PARAMETER
--------------------------------------
Uncertain by orders of magnitude (10⁻⁴ to 10⁻¹; Libbrecht & Rickerby 2013,
Bouvet et al. 2022, Libbrecht 2017).

DIAGNOSTICS
-----------
  β_eff / β_target     Thin-interface correction. Should be > 0.9 per D channel.
  Kinetic/curvature    β_sub·vₙ / (d₀·κ). Should be ≫ 1 for kinetic regime.
  Pe                   ε·vₙ/Dᵥ. Should be ≪ 1.
  ξᵥ ≤ ρ_vs/ρᵢ         K&P Eq. (48).
  τ_sub bracket        Kinetic vs correction dominance (M&F regime check).

Usage (CLI)
-----------
  python comp_eps.py --Lx 6e-4 --Ly 6e-4 --Lz 6e-4 --T0 -10 --alpha 1e-3
  python comp_eps.py --Lx 6e-4 --T0 -20 --alpha_range 1e-4 1e-1 12
  python comp_eps.py --Lx 6e-4 --T0 -5 --alpha 1e-3 --patch inputs/sim.opts
  python comp_eps.py --Lx 6e-4 --T0 -5 --alpha 1e-3 --quiet
"""

import argparse
import math
import re
from pathlib import Path

# =========================================================================
# Physical constants — match material_properties.c
# =========================================================================
_M_H2O   = 2.99e-26       # water molecule mass [kg]
_K_B     = 1.38e-23       # Boltzmann [J/K]
_K_I     = 2.29           # thermal conductivity of ice [W/m/K]
_K_A     = 0.02           # thermal conductivity of air [W/m/K]
_C_I     = 1.8e6          # volumetric heat capacity of ice [J/m³/K]
_C_A     = 1.341 * 1004   # volumetric heat capacity of air [J/m³/K]
_RHO_ICE = 919.0          # ice density [kg/m³]
_DV0     = 2.178e-5       # vapor diffusivity at 273.15 K [m²/s]
_PATM    = 101325.0       # atmospheric pressure [Pa]
_RHO_AIR = 1.341          # air density [kg/m³]
_BB      = 0.62           # ASHRAE humidity ratio coeff
_GAMMA   = 0.109          # ice–vapor surface energy γ [J/m²]  (K&P Table I)
_VM_ICE  = 1.963e-5       # molar volume of ice [m³/mol]
_R_GAS   = 8.314          # gas constant [J/mol/K]

# M&F asymptotic constants (SI footnote 1; from Karma-Rappel procedure
# applied to M&F's φ∈[0,1] well, NOT Karma's ±1 well).
_A1 = 5.0            # (M&F: "a₁ ≈ 5")
_A2 = 0.1581         # (M&F: "a₂ ≈ 0.1581")

# ASHRAE saturation-pressure polynomial
_KJ = [-0.5865e4, 0.2224e2, 0.1375e-1, -0.3403e-4, 0.2697e-7, 0.6918]


# =========================================================================
# Thermodynamic functions
# =========================================================================

def rho_vs_sat(T_C: float) -> float:
    """Saturation vapor density over ice ρ_vs [kg/m³] (ASHRAE polynomial)."""
    T = T_C + 273.15
    exp_arg = (_KJ[0]/T + _KJ[1] + _KJ[2]*T + _KJ[3]*T**2
               + _KJ[4]*T**3 + _KJ[5]*math.log(T))
    Pvs = math.exp(exp_arg)
    return _RHO_AIR * _BB * Pvs / (_PATM - Pvs)


def Dv_T(T_C: float) -> float:
    """Vapor diffusivity Dᵥ(T) = Dᵥ₀·(T/273.15)^1.81 [m²/s]."""
    return _DV0 * ((T_C + 273.15) / 273.15) ** 1.81


def beta_HK(T_C: float, alpha_c: float) -> float:
    """
    Hertz-Knudsen kinetic coefficient β_HK [s/m]:
        β_HK = (1/α_c) · √(2π m / k_B T)

    IMPORTANT convention warning:
    * K&P denote this as β' or β₀' (the *scaled* form).
    * M&F write β_sub·(ρ_vs/ρᵢ) — this product also equals β_HK.
    * K&P's β₀ (unscaled) = β_HK / (ρ_vs/ρᵢ) — this is what M&F call β_sub
      in Table S1 (~10⁴–10⁶ s/m).

    Both conventions give the same numerical result for K&P Eq. 43 bounds
    and M&F SI Eq. 9 model parameters, as long as scaling is consistent.
    In this script:
       beta_HK        = Hertz-Knudsen (SCALED,  K&P β',   M&F β_sub·(ρ_vs/ρᵢ))
       beta_unscaled  = β_HK/(ρ_vs/ρᵢ)  (K&P β₀, M&F β_sub)
    """
    T_K = T_C + 273.15
    return (1.0 / alpha_c) * math.sqrt(2.0 * math.pi * _M_H2O / (_K_B * T_K))


def capillary_length(T_C: float, gamma: float = _GAMMA) -> float:
    """Physical capillary length d₀ = γ V_m/(R T) [m] (K&P Eq. 13).

    NOTE: the C code (monitoring.c flag_Tdep branch) hardcodes
    d0 = 2.548e-7/T_K, i.e. γ·V_m/R = 2.548e-7, ~1% below this function's
    0.109*1.963e-5/8.314 = 2.574e-7. Reconcile before quantitative
    validation runs; comp_eps prints both when run verbosely.
    """
    T_K = T_C + 273.15
    return gamma * _VM_ICE / (_R_GAS * T_K)


def capillary_length_code(T_C: float) -> float:
    """d₀ exactly as hardcoded in monitoring.c: 2.548e-7 / T_K [m]."""
    return 2.548e-7 / (T_C + 273.15)


# =========================================================================
# Libbrecht temperature-dependent attachment kinetics
# (exact ports of Sigma0() in material_properties.c and the alpha model in
#  the flag_Tdep branch of monitoring.c)
# =========================================================================

_SIG0_S = [3.0e-3, 4.1e-3, 5.5e-3, 8.0e-3, 4.0e-3,
           6.0e-3, 3.5e-2, 7.0e-2, 1.1e-1, 0.75]
_SIG0_T = [-0.0001, -2.0, -4.0, -6.0, -7.0,
           -10.0, -20.0, -30.0, -40.0, -100.0]


def sigma0(T_C: float) -> float:
    """Libbrecht critical supersaturation σ₀(T) [-].

    Exact port of Sigma0() in material_properties.c: log10-log10
    interpolation of the lookup table in |T| (°C), flat beyond both table
    ends. NOTE the table is non-monotonic (bump at -6/-7 °C), so σ₀(T) has a
    kink there — this is inherited from the C code by design.
    """
    if T_C > _SIG0_T[0]:
        return _SIG0_S[0]
    interv = 0
    for ii in range(10):
        if T_C <= _SIG0_T[ii]:
            interv = ii
    if interv == 9:
        return _SIG0_S[9]
    t0, t1 = abs(_SIG0_T[interv]), abs(_SIG0_T[interv + 1])
    s0, s1 = _SIG0_S[interv], _SIG0_S[interv + 1]
    return 10.0 ** (math.log10(s0)
                    + (math.log10(s1) - math.log10(s0))
                    / (math.log10(t1) - math.log10(t0))
                    * (math.log10(abs(T_C)) - math.log10(t0)))


def alpha_libbrecht(T_C: float, sigma_surf: float) -> float:
    """Libbrecht attachment coefficient α(T, σ) = exp(-σ₀(T)/σ) [-].

    Exact port of the model in monitoring.c's flag_Tdep branch, including
    its floor: σ_surf < σ₀/ln(1e30) = σ₀/69.0775  ->  α = 1e-30.

    σ_surf is the LOCAL supersaturation |ρᵥ - ρ_vs|/ρ_vs. In a saturated
    sintering problem the driving σ is Gibbs-Thomson-scale (d₀·κ ~ 1e-5 to
    1e-3), which puts α in the exp(-30)..exp(-3000) regime — enormously
    smaller than the constant α_c ~ 1e-3..1e-2 typically assumed. The α you
    feed the eps bounds therefore depends strongly on which σ you consider
    characteristic; sweep it before trusting a mesh.
    """
    s0 = sigma0(T_C)
    if sigma_surf < s0 / 69.0775:
        return 1.0e-30
    return math.exp(-s0 / sigma_surf)


# =========================================================================
# Derived phase-field model parameters (M&F SI Eq. 9)
# =========================================================================

def derived_pf_params(eps: float, T0_C: float, d0: float, beta_hk: float,
                      alpha_i: float, alpha_a: float, Dv: float) -> dict:
    """
    Solve M&F SI Eq. (9) for the phase-field model parameters (λ, τ, M, α_source)
    given ε, physical target (β_sub, d_sub) and diffusivities.

    beta_hk is the Hertz-Knudsen value = β_sub·(ρ_vs/ρᵢ) in M&F notation.
    (Equivalently, if you have M&F's unscaled β_sub, then β_HK = β_sub·rho_rat.)
    """
    rho_vs  = rho_vs_sat(T0_C)
    rho_rat = rho_vs / _RHO_ICE
    D_ia    = 0.5 * (alpha_i + alpha_a)          # D*_ia arithmetic mean

    # M&F SI Eq. 9, first line:  d_sub · (ρ_vs/ρᵢ) = a₁ · ε / λ_sub
    lam_sub = _A1 * eps / (d0 * rho_rat)

    # M&F SI Eq. 9, second line: β_sub·(ρ_vs/ρᵢ) = β_HK
    #   β_HK = a₁ · [τ/(ε·λ) − a₂·ε/D*_ia − a₂·ε/Dᵥ]
    # ⇒ τ_sub = ε·λ_sub · [β_HK/a₁ + a₂·ε/D*_ia + a₂·ε/Dᵥ]
    kinetic_term = beta_hk / _A1
    corr_thermal = _A2 * eps / D_ia
    corr_vapor   = _A2 * eps / Dv
    tau_sub = eps * lam_sub * (kinetic_term + corr_thermal + corr_vapor)

    M_sub     = eps / (3.0 * tau_sub)
    alpha_src = lam_sub / tau_sub

    # Thin-interface correction per K&P Eq. 40 form:
    #   β_eff / β_target = 1 − a₁·a₂·ε / (D · β_HK)
    def _beta_ratio(D_channel):
        return 1.0 - _A1 * _A2 * eps / (D_channel * beta_hk)
    beta_ratio_heat_ice = _beta_ratio(alpha_i)
    beta_ratio_heat_air = _beta_ratio(alpha_a)
    beta_ratio_vapor    = _beta_ratio(Dv)

    return dict(
        lam_sub=lam_sub, tau_sub=tau_sub,
        M_sub=M_sub, alpha_src=alpha_src,
        D_ia=D_ia,
        kinetic_term=kinetic_term,
        corr_thermal=corr_thermal,
        corr_vapor=corr_vapor,
        beta_ratio_heat_ice=beta_ratio_heat_ice,
        beta_ratio_heat_air=beta_ratio_heat_air,
        beta_ratio_vapor=beta_ratio_vapor,
    )


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
    safety: float = 1.0,
    v_n: float = 1.0e-9,
    xi_v: float | None = None,
) -> dict:
    """Compute ε, mesh, and derived phase-field model parameters."""
    rho_vs   = rho_vs_sat(T0_C)
    rho_rat  = rho_vs / _RHO_ICE
    beta_hk  = beta_HK(T0_C, alpha_c)              # β_HK = β_scaled (K&P β')
    beta_uns = beta_hk / rho_rat                    # β_unscaled (K&P β₀ = M&F β_sub)
    d0       = capillary_length(T0_C)
    Dv       = Dv_T(T0_C)
    alpha_i  = _K_I / _C_I
    alpha_a  = _K_A / _C_A

    # K&P Eq. 43 bounds: W ≪ D · β' (uses the SCALED coefficient β' = β_HK)
    b_heat_ice = alpha_i * beta_hk
    b_heat_air = alpha_a * beta_hk
    b_vapor    = Dv      * beta_hk
    # K&P Eq. 45: W ≪ d₀ / (β₀·vₙ) (uses the UNSCALED coefficient β₀)
    b_kinetic  = d0 / (beta_uns * v_n)
    b_curv     = Rave

    bounds = {
        "B-HEAT-ICE": b_heat_ice,
        "B-HEAT-AIR": b_heat_air,
        "B-VAPOR":    b_vapor,
        "B-KINETIC":  b_kinetic,
        "B-CURV":     b_curv,
    }
    eps_max = min(bounds.values())
    binding = min(bounds, key=bounds.get)
    eps     = safety * eps_max

    # Mesh: h = eps/sqrt(2) → ~7.5 elements across phi=0.05–0.95 visible band
    h  = eps / math.sqrt(2.0)
    Nx = math.ceil(Lx / h) if Lx > 0 else 0
    Ny = math.ceil(Ly / h) if Ly > 0 else 0
    Nz = math.ceil(Lz / h) if Lz > 0 else 0

    # Derived PF parameters (M&F SI Eq. 9)
    dpf = derived_pf_params(eps, T0_C, d0, beta_hk, alpha_i, alpha_a, Dv)

    # Diagnostics
    Pe = eps * v_n / Dv

    xi_v_max  = rho_rat
    xi_v_warn = None
    if xi_v is not None and xi_v > xi_v_max:
        xi_v_warn = (f"WARNING: ξᵥ = {xi_v:.2e} exceeds ρ_vs/ρᵢ = "
                     f"{xi_v_max:.2e} (K&P Eq. 48). Quasi-steady assumption "
                     f"violated.")

    eps_over_Rave = eps / Rave if Rave > 0 else 0.0
    geom_warn = None
    if eps_over_Rave > 0.10:
        geom_warn = (f"WARNING: ε/R_ave = {eps_over_Rave:.1%} > 10%. "
                     f"Curvature-driven dynamics will be substantially wrong.")
    elif eps_over_Rave > 0.05:
        geom_warn = (f"NOTE: ε/R_ave = {eps_over_Rave:.1%} (5–10%). "
                     f"Borderline Gibbs-Thomson accuracy.")

    kappa = 1.0 / Rave if Rave > 0 else 0.0
    # d₀ inflation regime check: β·vₙ (Gibbs-Thomson kinetic term) vs d₀·κ (capillary term)
    # Use unscaled β for consistency with the physical Gibbs-Thomson condition.
    kin_over_curv = ((beta_uns * v_n) / (d0 * kappa)
                     if (d0 * kappa) > 0 else float("inf"))

    tau_warn = None
    if dpf["tau_sub"] <= 0:
        tau_warn = "ERROR: τ_sub ≤ 0. Asymptotic expansion is ill-posed."

    # τ_sub bracket dominance (kinetic fraction — M&F sharp-interface regime check)
    ttot = dpf["kinetic_term"] + dpf["corr_thermal"] + dpf["corr_vapor"]
    kinetic_frac = dpf["kinetic_term"] / ttot if ttot > 0 else 0.0

    bracket_warn = None
    if kinetic_frac < 0.10:
        bracket_warn = (f"WARNING: kinetic term is only {kinetic_frac:.1%} of τ_sub "
                        f"bracket. τ is dominated by thin-interface corrections; "
                        f"the sharp-interface analogy is severely stretched.")
    elif kinetic_frac < 0.50:
        bracket_warn = (f"NOTE: kinetic term is {kinetic_frac:.1%} of τ_sub bracket. "
                        f"Sharp-interface regime is marginal.")

    min_beta_ratio = min(dpf["beta_ratio_heat_ice"],
                         dpf["beta_ratio_heat_air"],
                         dpf["beta_ratio_vapor"])

    return dict(
        eps=eps, Nx=Nx, Ny=Ny, Nz=Nz,
        eps_max=eps_max, binding=binding, **bounds,
        Pe=Pe, rho_vs=rho_vs, rho_rat=rho_rat,
        beta_hk=beta_hk, beta_uns=beta_uns,
        d0=d0, Dv=Dv,
        alpha_i=alpha_i, alpha_a=alpha_a,
        v_n=v_n, safety=safety, Rave=Rave, alpha_c=alpha_c,
        lam_sub=dpf["lam_sub"], tau_sub=dpf["tau_sub"],
        M_sub=dpf["M_sub"], alpha_src=dpf["alpha_src"],
        D_ia=dpf["D_ia"],
        kinetic_term=dpf["kinetic_term"],
        corr_thermal=dpf["corr_thermal"],
        corr_vapor=dpf["corr_vapor"],
        kinetic_frac=kinetic_frac,
        beta_ratio_heat_ice=dpf["beta_ratio_heat_ice"],
        beta_ratio_heat_air=dpf["beta_ratio_heat_air"],
        beta_ratio_vapor=dpf["beta_ratio_vapor"],
        min_beta_ratio=min_beta_ratio,
        kin_over_curv=kin_over_curv,
        xi_v_max=xi_v_max, xi_v_warn=xi_v_warn,
        eps_over_Rave=eps_over_Rave, geom_warn=geom_warn,
        bracket_warn=bracket_warn, tau_warn=tau_warn,
    )


def alpha_c_sensitivity(T0_C, alpha_lo, alpha_hi, n_points=10, **kwargs):
    """Sweep α_c log-uniformly. Returns (results, alphas)."""
    alphas = [
        math.exp(math.log(alpha_lo) + i * math.log(alpha_hi/alpha_lo)/(n_points-1))
        for i in range(n_points)
    ]
    return [compute_eps(T0_C=T0_C, alpha_c=a, **kwargs) for a in alphas], alphas


# =========================================================================
# Opts-file patching
# =========================================================================

def patch_opts(opts_path: Path, params: dict, dim: int) -> None:
    """Replace ε, mesh, and derived model parameters in opts_path."""
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
    text = _replace_or_append(text, "-M_sub",     f"{params['M_sub']:.4e}")
    text = _replace_or_append(text, "-alpha_sub", f"{params['alpha_src']:.4e}")
    text = _replace_or_append(text, "-tau_sub",   f"{params['tau_sub']:.4e}")
    text = _replace_or_append(text, "-lam_sub",   f"{params['lam_sub']:.4e}")
    text = _replace_or_append(text, "-beta_sub",  f"{params['beta_sub']:.4e}")
    text = _replace_or_append(text, "-d0_sub",    f"{params['d0']:.4e}")

    opts_path.write_text(text)
    print(f"  Patched {opts_path}")


# =========================================================================
# Output
# =========================================================================

def _print_header(title: str) -> None:
    print("=" * 76)
    print(f"  {title}")
    print("=" * 76)


def _print_single(args, p: dict, alpha_c: float, dim: int) -> None:
    _print_header("Phase-field parameters — K&P (2009) + M&F (2024) SI Eq. 9 [sublimation]")

    print(f"\n--- Input parameters ---")
    print(f"  T₀                       = {args.T0:.1f} °C  ({args.T0+273.15:.2f} K)")
    print(f"  α_c (condensation coeff) = {alpha_c:.3e}   [only genuine free parameter]")
    print(f"  vₙ (front velocity)      = {args.vn:.1e}  m/s")
    print(f"  R_ave (grain radius)     = {args.Rave:.2e} m")
    print(f"  Safety factor            = {args.safety:.2f}  (ε = safety · min(bounds))")
    print(f"  Asymptotic constants     : a₁ = {_A1}, a₂ = {_A2}  (M&F SI footnote 1)")

    print(f"\n--- Thermodynamic state at T₀ ---")
    print(f"  ρ_vs(T₀)         = {p['rho_vs']:.4e}  kg/m³")
    print(f"  ρ_vs/ρᵢ  (χ)     = {p['rho_rat']:.4e}  [-]")
    print(f"  Dᵥ(T₀)           = {p['Dv']:.4e}  m²/s")
    print(f"  κᵢ/Cᵢ  (D_ice)   = {p['alpha_i']:.4e}  m²/s")
    print(f"  κₐ/Cₐ  (D_air)   = {p['alpha_a']:.4e}  m²/s")
    print(f"  D*_ia (mean)     = {p['D_ia']:.4e}  m²/s")

    print(f"\n--- Physical targets (Gibbs-Thomson parameters) ---")
    print(f"  β_HK (Hertz-Knudsen)   = {p['beta_hk']:.4e}  s/m   [K&P β', SCALED form]")
    print(f"                           [= β_sub·(ρ_vs/ρᵢ) in M&F notation]")
    print(f"  β_sub (physical)       = {p['beta_uns']:.4e}  s/m   [K&P β₀, M&F Table S1 form]")
    print(f"                           [= β_HK/(ρ_vs/ρᵢ); matches K&P range 3e4–3e6]")
    print(f"  d₀_sub  (physical)     = {p['d0']:.4e}  m     [γ·V_m/RT; K&P Eq. 13]")

    print(f"\n--- Upper bounds on ε (K&P eqs. 42–46) ---")
    labels = {
        "B-HEAT-ICE": "Eq.(43a)  (κᵢ/Cᵢ)·β_HK  ",
        "B-HEAT-AIR": "Eq.(43b)  (κₐ/Cₐ)·β_HK  ",
        "B-VAPOR":    "Eq.(43c)  Dᵥ·β_HK        ",
        "B-KINETIC":  "Eq.(45)   d₀/(β_sub·vₙ) ",
        "B-CURV":     "Geometric R_ave           ",
    }
    for key, label in labels.items():
        marker = "  ← BINDING" if key == p["binding"] else ""
        print(f"  {label}  =  {p[key]:.4e} m{marker}")

    h = p['eps'] / math.sqrt(2.0)
    print(f"\n--- Chosen ε and mesh ---")
    print(f"  ε                        = {p['eps']:.4e}  m")
    print(f"  h = ε/√2                 = {h:.4e}  m    (mesh rule: ~7.5 elements across phi=0.05-0.95)")
    print(f"  Nx = ceil(Lx·√2/ε)       = {p['Nx']}")
    if dim >= 2 and p["Ny"] > 0: print(f"  Ny = ceil(Ly·√2/ε)       = {p['Ny']}")
    if dim == 3 and p["Nz"] > 0: print(f"  Nz = ceil(Lz·√2/ε)       = {p['Nz']}")

    print(f"\n--- Derived phase-field model parameters (M&F SI Eq. 9) ---")
    print(f"  These are your solver inputs — recompute at every ε.")
    print(f"")
    print(f"  λ_sub  (coupling)    = {p['lam_sub']:.4e}  [-]     [= a₁·ε/(d₀·(ρ_vs/ρᵢ))]")
    print(f"  τ_sub  (relaxation)  = {p['tau_sub']:.4e}  s       [from M&F SI Eq. 9]")
    print(f"  M_sub  (mobility)    = {p['M_sub']:.4e}  m/s     [= ε/(3·τ_sub)]")
    print(f"  α_sub  (source rate) = {p['alpha_src']:.4e}  1/s     [= λ_sub/τ_sub]")

    ttot = p['kinetic_term'] + p['corr_thermal'] + p['corr_vapor']
    print(f"\n  τ_sub bracket [β_HK/a₁ + a₂·ε/D*_ia + a₂·ε/Dᵥ] decomposition:")
    print(f"    β_HK/a₁       kinetic  = {p['kinetic_term']:.4e} s/m  ({p['kinetic_term']/ttot:5.1%})")
    print(f"    a₂·ε/D*_ia    thermal  = {p['corr_thermal']:.4e} s/m  ({p['corr_thermal']/ttot:5.1%})")
    print(f"    a₂·ε/Dᵥ       vapor    = {p['corr_vapor']:.4e} s/m  ({p['corr_vapor']/ttot:5.1%})")
    if p["bracket_warn"]:
        print(f"  {p['bracket_warn']}")
    else:
        print(f"  → Kinetic-dominated τ_sub ⇒ we are in the correct sharp-interface regime.")

    if p["tau_warn"]:
        print(f"\n  {p['tau_warn']}")

    print(f"\n--- Thin-interface correction: β_eff/β_target (K&P Eq. 40 form) ---")
    print(f"  = 1 − a₁·a₂·ε / (D · β_HK). Should be > 0.9 per D channel.")
    print(f"    heat-in-ice channel   = {p['beta_ratio_heat_ice']:.4f}")
    print(f"    heat-in-air channel   = {p['beta_ratio_heat_air']:.4f}")
    print(f"    vapor channel         = {p['beta_ratio_vapor']:.4f}")

    print(f"\n--- Kinetic vs. curvature dominance (d₀ inflation regime check) ---")
    print(f"  β_sub·vₙ / (d₀·κ)  with κ ≈ 1/R_ave  = {p['kin_over_curv']:.2e}")
    if p['kin_over_curv'] > 100:
        print(f"  → Kinetic-dominated regime. Physical d₀ is safe; inflation OK too.")
    elif p['kin_over_curv'] > 10:
        print(f"  → Mostly kinetic-driven. Physical d₀ recommended; inflation borderline.")
    elif p['kin_over_curv'] > 1:
        print(f"  → Comparable regimes. Use physical d₀; do NOT inflate.")
    else:
        print(f"  → Curvature-dominated (ETM-like). Physical d₀ mandatory; ξ_v suspect.")

    print(f"\n--- Péclet check ---")
    Pe_status = ("OK (≪ 1)" if p["Pe"] < 0.1
                 else ("borderline" if p["Pe"] < 1.0
                       else "VIOLATED — Pe ≥ 1!"))
    print(f"  Pe = ε·vₙ/Dᵥ    = {p['Pe']:.4e}  [{Pe_status}]")

    print(f"\n--- Geometric accuracy (ε ≪ R_grain) ---")
    geom_status = ("OK (< 5%)" if p["eps_over_Rave"] < 0.05
                   else ("BORDERLINE (5–10%)" if p["eps_over_Rave"] < 0.10
                         else "VIOLATED (> 10%)"))
    print(f"  ε / R_ave         = {p['eps_over_Rave']:.1%}  [{geom_status}]")
    if p["geom_warn"]:
        print(f"  {p['geom_warn']}")

    print(f"\n--- ξᵥ validity check (K&P Eq. 48) ---")
    print(f"  ξᵥ hard upper bound = ρ_vs/ρᵢ = {p['xi_v_max']:.4e}")
    if p["xi_v_warn"]:
        print(f"  {p['xi_v_warn']}")
    elif args.xiv is not None:
        print(f"  ξᵥ = {args.xiv:.2e}  → OK (ξᵥ ≤ ρ_vs/ρᵢ)")
    else:
        print(f"  (ξᵥ not provided; pass --xiv to check)")

    print(f"\n--- Recommended opts entries ---")
    print(f"  -eps        {p['eps']:.4e}")
    print(f"  -Nx         {p['Nx']}")
    if dim >= 2 and p["Ny"] > 0: print(f"  -Ny         {p['Ny']}")
    if dim == 3 and p["Nz"] > 0: print(f"  -Nz         {p['Nz']}")
    print(f"  -M_sub      {p['M_sub']:.4e}")
    print(f"  -alpha_sub  {p['alpha_src']:.4e}")
    print(f"  -tau_sub    {p['tau_sub']:.4e}")
    print(f"  -lam_sub    {p['lam_sub']:.4e}")
    print(f"  -beta_sub   {p['beta_uns']:.4e}")
    print(f"  -d0_sub     {p['d0']:.4e}")

    print("=" * 76)


def _print_sweep(results, alphas, T0_C: float) -> None:
    _print_header(f"α_c sensitivity sweep at T₀ = {T0_C:.1f} °C")
    hdr = (f"  {'α_c':>10s} {'β_sub[s/m]':>11s} {'eps_max[m]':>11s} "
           f"{'eps[m]':>10s} {'binding':>11s} {'M_sub[m/s]':>11s} "
           f"{'α_sub[1/s]':>11s} {'kin.frac':>8s}")
    print(hdr)
    print("  " + "-" * 87)
    for p, a in zip(results, alphas):
        print(f"  {a:10.3e} {p['beta_uns']:11.3e} {p['eps_max']:11.3e} "
              f"{p['eps']:10.3e} {p['binding']:>11s} {p['M_sub']:11.3e} "
              f"{p['alpha_src']:11.3e} {p['kinetic_frac']:8.2%}")
    print("=" * 76)
    print("  Note: ε = safety · eps_max.  M_sub, α_sub from M&F SI Eq. 9.")
    print("  kin.frac = kinetic term / (bracket total). Should be ≥ 50% for the")
    print("  sharp-interface regime; smaller values ⇒ τ is dominated by corrections.")
    print("=" * 76)


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
        description=("Compute ε, mesh, and derived phase-field model parameters "
                     "per K&P (2009) and M&F (2024) SI Eq. 9. Mesh rule: h = ε/√2."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--Lx",    type=float, required=True,  help="Domain length x [m]")
    ap.add_argument("--Ly",    type=float, default=0.0,    help="Domain length y [m]; 0 = unused")
    ap.add_argument("--Lz",    type=float, default=0.0,    help="Domain length z [m]; 0 = unused")
    ap.add_argument("--Rave",  type=float, default=3.0e-5, help="Representative grain radius [m]")
    ap.add_argument("--T0",    type=float, default=-20.0,  help="Mean temperature [°C]")
    ap.add_argument("--alpha", type=float, default=None,
                    help="Condensation coefficient α_c [−]. Mutually exclusive with --alpha_range.")
    ap.add_argument("--alpha_range", nargs=3, metavar=("LO","HI","N"), default=None,
                    help="α_c sweep: LO HI N_points (log-uniform).")
    ap.add_argument("--sigma_surf", type=float, default=None,
                    help="Characteristic supersaturation σ [-]: derive α_c from the "
                         "code's Libbrecht model α = exp(-σ₀(T0)/σ) (Sigma0 lookup, "
                         "monitoring.c flag_Tdep). Mutually exclusive with "
                         "--alpha/--alpha_range.")
    ap.add_argument("--vn",     type=float, default=1.0e-9,
                    help="Normal front velocity vₙ [m/s] for Eq.(45)")
    ap.add_argument("--xiv",    type=float, default=None,
                    help="Time-scaling factor ξᵥ (optional; triggers K&P Eq. 48 check)")
    ap.add_argument("--safety", type=float, default=0.5,
                    help="Safety factor: ε = safety · min(bounds).")
    ap.add_argument("--dim",    type=int,   default=None,
                    help="Problem dimension (1/2/3); inferred from Ly/Lz if omitted")
    ap.add_argument("--patch",  type=Path,  default=None,
                    help="Opts file to update in place")
    ap.add_argument("--quiet",  action="store_true",
                    help="Print only ε and Nx/Ny/Nz (one per line)")
    args = ap.parse_args()

    n_modes = sum(x is not None for x in (args.alpha, args.alpha_range, args.sigma_surf))
    if n_modes > 1:
        ap.error("--alpha, --alpha_range, and --sigma_surf are mutually exclusive.")
    if n_modes == 0:
        ap.error("Provide one of --alpha, --alpha_range, or --sigma_surf.")

    if args.sigma_surf is not None:
        s0 = sigma0(args.T0)
        args.alpha = alpha_libbrecht(args.T0, args.sigma_surf)
        print(f"\n  Libbrecht kinetics mode: sigma0({args.T0:g} C) = {s0:.4e}, "
              f"sigma_surf = {args.sigma_surf:.3e}")
        print(f"  -> alpha_c = exp(-sigma0/sigma_surf) = {args.alpha:.4e}"
              + ("   [FLOORED at 1e-30: sigma_surf < sigma0/69.08 — kinetics "
                 "effectively frozen at this supersaturation]"
                 if args.alpha <= 1.0e-30 else ""))

    dim = _infer_dim(args.Lx, args.Ly, args.Lz, args.dim)

    common = dict(
        Lx=args.Lx, Ly=args.Ly, Lz=args.Lz,
        Rave=args.Rave, v_n=args.vn,
        safety=args.safety, xi_v=args.xiv,
    )

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

        d0_code = capillary_length_code(args.T0)
        if abs(d0_code - p["d0"]) / p["d0"] > 1e-3:
            print(f"\n  NOTE: monitoring.c (flag_Tdep) hardcodes d0 = 2.548e-7/T_K "
                  f"= {d0_code:.4e} m,\n        {100*(d0_code/p['d0']-1):+.1f}% vs "
                  f"this script's gamma*V_m/(R*T) = {p['d0']:.4e} m — reconcile "
                  f"before quantitative validation.")

        if p["xi_v_warn"]:
            print(f"\n  {p['xi_v_warn']}\n")

        if args.patch:
            print()
            patch_opts(args.patch, p, dim)
        return

    lo, hi, n = float(args.alpha_range[0]), float(args.alpha_range[1]), int(args.alpha_range[2])
    results, alphas = alpha_c_sensitivity(
        T0_C=args.T0, alpha_lo=lo, alpha_hi=hi, n_points=n, **common
    )
    _print_sweep(results, alphas, args.T0)


if __name__ == "__main__":
    _cli()