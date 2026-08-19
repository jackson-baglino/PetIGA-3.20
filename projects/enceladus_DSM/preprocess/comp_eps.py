"""
comp_eps.py — Compute the diffuse-interface parameter `eps`, the mesh, and the
derived phase-field parameters for the dry-snow sublimation model.

Sources: Kaempfer & Plapp, Phys. Rev. E 79, 031502 (2009) ["K&P"], and
Moure & Fu, Cryst. Growth Des. 24, 7808 (2024) ["M&F"], sublimation reduction.

WIDTH CONVENTION
----------------
This solver uses M&F's well, (1/2)*phi^2*(1-phi)^2 with an eps^2 gradient
term, so the 1D equilibrium profile is logistic:

    phi(x) = 1/(1 + exp(-x/eps)) = 0.5*(1 + tanh(x/(2*eps)))

`eps` here is exactly M&F's eps, and a1 = 5, a2 = 0.1581 below are the
asymptotic constants for THAT well. K&P's Eq.(33) uses a +-1 well whose
profile is tanh(x/(sqrt(2)*W)), so their W = sqrt(2)*eps and their constants
are a1 = 5*sqrt(2)/8, a2 = 0.6267. The two sets agree: 5 == (5*sqrt(2)/8)*4,
and a1*a2 == 0.79 in eps-units vs 0.78 in W-units. Nothing here needs a
sqrt(2): every width-carrying quantity goes through a1/a2, and the bounds
below are written in eps the way M&F write them.

MESH
----
    h = eps/sqrt(2)      Nx = ceil(Lx/h) = ceil(sqrt(2)*Lx/eps)

That is 2 elements per asymptotic width W, and (since the logistic profile
has no edge) 8.3 elements across phi=0.05-0.95 and 13.0 across 0.01-0.99.
See docs/interface_width_conventions.md before comparing against a ParaView
measurement.

BOUNDS ON eps
-------------
  [B-HEAT]     eps <= safety * D_heat * beta_HK    K&P Eq.(43a/b)
  [B-VAPOR]    eps <= safety * D_v    * beta_HK    K&P Eq.(43c)
  [B-KINETIC]  eps <= safety * d0/(beta_sub*vn)    K&P Eq.(45)
  [B-CURV]     eps <= eps_over_R * R_ave           geometric, K&P "W << 1/chi"

B-HEAT/B-VAPOR control the thin-interface expansion parameter

    delta = a1*a2*eps/(D*beta_HK) = a1*a2*safety     (on the binding channel)

which is what "<<" means quantitatively. tau_sub already COMPENSATES the
first order in delta (it carries the a2 terms), so delta is the size of the
correction being applied, not the residual error; the residual is O(delta^2).
safety=0.5 -> delta=0.40, safety=0.25 -> delta=0.20.

D_heat is set by --Dchannel. "mean" (default) = D*_ia = (D_i+D_a)/2, the
diffusivity tau_sub compensates against and the one the solver assembles
(enceladus_main.c diff_sub). "ice" = kappa_i/C_i, ~6x tighter. The ice channel
is knowingly violated: enforcing it drove the Demmenie mesh to ~1e9 nodes
(tried 2026-08-11, reverted). M&F violate it by 100x in print and say so
("we consider a value eps < 10^-6 m, which allows us to use a coarser spatial
discretization ... less accurate in capturing interface kinetics along the
air-ice ... interfaces"). delta_ice is reported as a diagnostic instead.

B-CURV is deliberately crude: it is 1-2 decades looser than the binding bound
in every case we run, and the true smallest curvature (a fresh neck fillet)
starts near-singular and then relaxes, so sizing a mesh from it is both
unstable and needlessly expensive. eps <= 0.05*R_ave, no safety factor.

DERIVED PARAMETERS (M&F SI Eq. 9)
---------------------------------
   d_sub*(rho_vs/rho_i) = a1*eps/lam_sub
   beta_sub*(rho_vs/rho_i) = a1*[tau_sub/(eps*lam_sub) - a2*eps/D*_ia - a2*eps/Dv]
   M_sub = eps/(3*tau_sub)      alpha_sub = lam_sub/tau_sub

beta_sub / beta_0 always mean the UNSCALED coefficient in the Gibbs-Thomson
condition; beta_HK = beta_sub*(rho_vs/rho_i) = (1/alpha_c)*sqrt(2*pi*m/kT) is
the scaled one (K&P's beta'). alpha_c is the only genuine free parameter
(10^-4 to 10^-1; Libbrecht 2017, Braun et al. 2024).

Usage
-----
  python comp_eps.py --Lx 6e-4 --Ly 6e-4 --T0 -10 --alpha 1e-3
  python comp_eps.py --Lx 6e-4 --T0 -20 --alpha_range 1e-4 1e-1 12
  python comp_eps.py --Lx 6e-4 --T0 -5 --alpha 1e-3 --patch inputs/sim.opts
"""

from __future__ import annotations   # PEP 604 (`float | None`) on Python 3.9

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
_GAMMA   = 0.109          # ice-vapor surface energy γ [J/m²]  (K&P Table I)
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
      in Table S1 (~10⁴-10⁶ s/m).

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
# Arrhenius attachment kinetics — the default alpha_c(T) model
# =========================================================================
# alpha_libbrecht() above is not usable for mesh sizing: its sigma0 table is
# non-monotonic (kink at -6/-7 C, so eps jumps with a 1 C change in T0), and
# alpha = exp(-sigma0/sigma) underflows at sintering supersaturations
# (alpha(-20 C) = 1e-30 at sigma = 4.5e-4, i.e. no sintering at all). Because
# sigma0 spans a factor of 27 over -2..-40 C, no single reference sigma fits
# both ends -- the incompatibility is structural.
#
# Instead: a smooth alpha_c(T) = A*exp(-Ea/(R*T)) anchored to the literature
# band (Libbrecht 2017; Braun, Fourteau & Lowe, The Cryosphere 18, 1653, 2024).
# Libbrecht's sigma0(T) sets only the SIGN of the T-dependence. Defaults imply
# Ea = 31.8 kJ/mol.
#
# alpha_c must be resolved BEFORE eps: alpha_c -> beta_sub -> binding eps bound
# -> mesh, and it also sets the transport regime via L* = beta_HK*D_v (below
# L* attachment-limited, above it vapour-diffusion-limited).

_R_GAS = 8.314462618            # J/(mol K)

# Literature band (default anchors and default clamp). Braun et al. 2024 find
# 1e-3 < alpha < 1e-1, but all their values sit at T ~ -8 C and drift within a
# single sample as the snow coarsens -- the spread is SUPERSATURATION
# dependence, not temperature. Treat the band as a consistency check on
# alpha_c(T, sigma), not a fit target.
ALPHA_LIT_LO, ALPHA_LIT_HI = 1.0e-3, 1.0e-1
ALPHA_ANCHOR_WARM = (-2.0, 1.0e-1)
ALPHA_ANCHOR_COLD = (-40.0, 1.0e-3)


def arrhenius_params(anchor_warm=ALPHA_ANCHOR_WARM, anchor_cold=ALPHA_ANCHOR_COLD):
    """(A, Ea) of alpha = A*exp(-Ea/(R*T)) through two (T_C, alpha) anchors.

    Returns (A [-], Ea [J/mol]). Ea > 0 means alpha falls as T falls, which
    is the direction Libbrecht's sigma0(T) implies.
    """
    (T1, a1), (T2, a2) = anchor_warm, anchor_cold
    T1K, T2K = T1 + 273.15, T2 + 273.15
    if a1 <= 0 or a2 <= 0 or abs(T1K - T2K) < 1e-9:
        raise ValueError("anchors must have positive alpha and distinct T")
    Ea = _R_GAS * math.log(a1 / a2) / (1.0 / T2K - 1.0 / T1K)
    A = a1 * math.exp(Ea / (_R_GAS * T1K))
    return A, Ea


def alpha_arrhenius(T_C: float,
                    anchor_warm=ALPHA_ANCHOR_WARM,
                    anchor_cold=ALPHA_ANCHOR_COLD,
                    clamp=(ALPHA_LIT_LO, ALPHA_LIT_HI)) -> float:
    """Condensation coefficient alpha_c(T) [-] from the anchored Arrhenius.

    Clamped to `clamp` (default the literature band) so that EXTRAPOLATING
    past the anchors cannot silently hand the mesh sizer an alpha_c the
    literature does not support -- which is exactly the failure mode of the
    Libbrecht form at sintering supersaturations.
    """
    A, Ea = arrhenius_params(anchor_warm, anchor_cold)
    a = A * math.exp(-Ea / (_R_GAS * (T_C + 273.15)))
    if clamp:
        a = min(max(a, clamp[0]), clamp[1])
    return a


# -------------------------------------------------------------------------
# Two-parameter Libbrecht fit: alpha_c(T, rho_v)
# -------------------------------------------------------------------------
#
# Preferred over the T-only Arrhenius. Libbrecht's form has TWO free
# parameters, A and sigma0; monitoring.c hardwires A = 1, which is what makes
# it underflow. Freeing A and rescaling the sigma0 table by a single factor f
# (keeping its shape) gives two unknowns for two anchors:
#
#     alpha_c(T, sigma) = A * exp(-f * sigma0(T) / sigma)
#
# A is the rough-surface ceiling on alpha -- exactly the quantity the
# literature bounds -- and f ~ 1/100 turns the -6/-7 C table kink into an 8%
# step while keeping Libbrecht's T-dependence rather than an invented Ea.
#
# The clamp floor is still required: alpha collapses as sigma -> 0 and
# beta_sub ~ 1/alpha diverges. The mesh does NOT blow up (Eq.45 becomes
# binding and is beta-independent under --vn_feature), but physical t_final
# scales as 1/alpha_c -- a 2.5 h experiment needs ~65 h simulated at
# alpha_c = 1e-4 and ~26 days at 1e-5. Runtime, not resolution, is the cost.

LIBBRECHT_SIGMA_CHAR = 4.5e-4     # characteristic sintering supersaturation


def libbrecht2_params(sigma_ref: float = LIBBRECHT_SIGMA_CHAR,
                      anchor_warm=ALPHA_ANCHOR_WARM,
                      anchor_cold=ALPHA_ANCHOR_COLD):
    """(A, f) for alpha = A*exp(-f*sigma0(T)/sigma), fitted at sigma_ref.

    f rescales the sigma0 table without changing its shape. Returns
    (A, f, sigma_ref) so callers can record what the fit was pinned to --
    A and f are only meaningful alongside it.
    """
    (Tw, aw), (Tc, ac) = anchor_warm, anchor_cold
    s0w, s0c = sigma0(Tw), sigma0(Tc)
    if abs(s0c - s0w) < 1e-15:
        raise ValueError("anchors must straddle a change in sigma0(T)")
    f = sigma_ref * math.log(aw / ac) / (s0c - s0w)
    A = aw * math.exp(f * s0w / sigma_ref)
    return A, f, sigma_ref


def alpha_libbrecht2(T_C: float, sigma: float,
                     sigma_ref: float = LIBBRECHT_SIGMA_CHAR,
                     anchor_warm=ALPHA_ANCHOR_WARM,
                     anchor_cold=ALPHA_ANCHOR_COLD,
                     clamp=(ALPHA_LIT_LO, ALPHA_LIT_HI)) -> float:
    """alpha_c(T, sigma) from the two-parameter fit, clamped.

    `sigma` is the local supersaturation |rho_v - rho_vs(T)|/rho_vs(T), so
    this is the alpha_c(T, rho_v) the solver would evaluate pointwise.
    """
    A, f, _ = libbrecht2_params(sigma_ref, anchor_warm, anchor_cold)
    if sigma <= 0.0:
        a = clamp[0] if clamp else 0.0
    else:
        expo = f * sigma0(T_C) / sigma
        a = 0.0 if expo > 700.0 else A * math.exp(-expo)
    if clamp:
        a = min(max(a, clamp[0]), clamp[1])
    return a


def sigma_char_of(T_C: float, r_feat: float) -> float:
    """Gibbs-Thomson supersaturation at the smallest resolved feature.

    sigma = d0*kappa = d0/r_feat. This is the driving supersaturation a neck
    fillet of radius r_feat actually sees, and it is what the alpha_c fit
    should be pinned to -- not a snow-crystal-growth sigma, which is orders
    of magnitude larger.
    """
    return capillary_length(T_C) / r_feat


def libbrecht_sigma_ref_for(alpha_target: float, T_C: float):
    """sigma that makes alpha_libbrecht(T_C, sigma) equal alpha_target.

    Diagnostic only. Inverting alpha = exp(-sigma0/sigma) gives
    sigma = sigma0(T) / ln(1/alpha). Evaluating this at one anchor and then
    applying it at the other is how the structural incompatibility in note 2
    above is demonstrated numerically rather than asserted.
    """
    if not (0.0 < alpha_target < 1.0):
        return float("nan")
    return sigma0(T_C) / math.log(1.0 / alpha_target)


# =========================================================================
# Heat-channel ε ceiling — single source of truth
# =========================================================================

def D_heat_of(thermal_channel: str = "mean") -> float:
    """Diffusivity for the Eq.(43) heat channel. See --Dchannel."""
    D_i = _K_I / _C_I
    D_a = _K_A / _C_A
    if thermal_channel == "mean":
        return 0.5 * (D_i + D_a)          # D*_ia; matches src/<project>_main.c:558
    if thermal_channel == "ice":
        return D_i                        # conservative single-sided K&P 43a
    raise ValueError(f"thermal_channel must be 'mean' or 'ice', got "
                     f"{thermal_channel!r}")


def eps_ceiling(T_C: float, alpha_c: float,
                thermal_channel: str = "mean",
                safety: float = 0.5,
                delta_target: float | None = None) -> float:
    """ε ceiling from the K&P Eq.(43) heat channel.

    safety (default 0.5):  ε = safety · D_heat · β_HK.
    delta_target (overrides safety): size ε from the thin-interface expansion
        parameter directly, ε = delta_target · D_heat · β_HK / (a₁·a₂).
    The two are the same knob: delta = a₁·a₂·safety = 0.79·safety.
    """
    bound = D_heat_of(thermal_channel) * beta_HK(T_C, alpha_c)
    if delta_target is not None:
        return delta_target * bound / (_A1 * _A2)
    return safety * bound


# =========================================================================
# Derived phase-field model parameters (M&F SI Eq. 9)
# =========================================================================

def derived_pf_params(eps: float, T0_C: float, d0: float, beta_hk: float,
                      alpha_i: float, alpha_a: float, Dv: float,
                      thermal_channel: str = "mean") -> dict:
    """
    Solve M&F SI Eq. (9) for (λ, τ, M, α_source) given ε, the physical targets
    (β_sub, d_sub) and the diffusivities.

    beta_hk is the scaled coefficient β_sub·(ρ_vs/ρᵢ). thermal_channel selects
    the D reported in the delta diagnostic; τ_sub always compensates against
    D*_ia, which is what the solver assembles (enceladus_main.c diff_sub).
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

    # Thin-interface expansion parameter (K&P Eq. 40 bracket):
    #   delta = a₁·a₂·ε/(D·β_HK)
    # tau_sub above compensates delta to first order, so this is the SIZE of
    # the correction, not the residual error; the residual is O(delta^2).
    def _delta(D_channel):
        return _A1 * _A2 * eps / (D_channel * beta_hk)
    D_heat = D_heat_of(thermal_channel)

    return dict(
        lam_sub=lam_sub, tau_sub=tau_sub,
        M_sub=M_sub, alpha_src=alpha_src,
        D_ia=D_ia, D_heat=D_heat, thermal_channel=thermal_channel,
        kinetic_term=kinetic_term,
        corr_thermal=corr_thermal,
        corr_vapor=corr_vapor,
        delta_heat=_delta(D_heat),
        delta_ia=_delta(D_ia),          # what tau_sub actually compensates
        delta_heat_ice=_delta(alpha_i),
        delta_heat_air=_delta(alpha_a),
        delta_vapor=_delta(Dv),
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
    safety: float = 0.5,   # match the CLI default and eps_ceiling()
    v_n: float = 1.0e-9,
    xi_v: float | None = None,
    thermal_channel: str = "mean",
    eps_over_R: float = 0.05,
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

    # K&P Eq. 43 bounds, in eps (M&F's own form). Evaluated against beta_hk,
    # the scaled coefficient beta'. --Dchannel picks the heat channel; the
    # single-channel readings are reported via delta_* but never enforced.
    D_heat     = D_heat_of(thermal_channel)
    b_heat     = D_heat * beta_hk             # Eq. 43a/b
    b_vapor    = Dv     * beta_hk             # Eq. 43c
    # K&P Eq. 45 uses the UNSCALED coefficient beta_0.
    b_kinetic  = d0 / (beta_uns * v_n)

    # Bounds that the safety factor applies to (they are all "<<" statements).
    bounds = {
        "B-HEAT":    b_heat,
        "B-VAPOR":   b_vapor,
        "B-KINETIC": b_kinetic,
    }
    eps_max = min(bounds.values())
    eps     = safety * eps_max
    binding = min(bounds, key=bounds.get)

    # B-CURV is already a fraction, so it does NOT get the safety factor on
    # top. It is 1-2 decades looser than the others in practice.
    b_curv = eps_over_R * Rave
    if b_curv < eps:
        eps, binding = b_curv, "B-CURV"
    bounds["B-CURV"] = b_curv

    _n_est = 1
    for _L in (Lx, Ly, Lz):
        if _L > 0:
            _n_est *= math.ceil(_L * math.sqrt(2.0) / eps)
    intractable = _n_est > 5.0e7

    # Mesh: h = eps/sqrt(2) → ~7.5 elements across phi=0.05-0.95 visible band
    h  = eps / math.sqrt(2.0)
    Nx = math.ceil(Lx / h) if Lx > 0 else 0
    Ny = math.ceil(Ly / h) if Ly > 0 else 0
    Nz = math.ceil(Lz / h) if Lz > 0 else 0

    # Derived PF parameters (M&F SI Eq. 9)
    dpf = derived_pf_params(eps, T0_C, d0, beta_hk, alpha_i, alpha_a, Dv,
                            thermal_channel=thermal_channel)

    # Diagnostics
    Pe = eps * v_n / Dv

    # xi_v: the M&F temporal scaling. It multiplies vapour DIFFUSION and the
    # rho_ice source in R_vap but NOT the storage term (src/assembly.c), so the
    # residual storage error is rho_vs/(xi_v*rho_i) -- meaning LARGER xi_v is
    # MORE quasi-steady, not less. This check used to warn when xi_v exceeded
    # rho_vs/rho_i "(K&P Eq. 48) -- quasi-steady assumption violated", which has
    # the sense backwards: Eq. 48 is what you need to RESOLVE the true vapour
    # transient, and exceeding it deliberately is the whole point of the
    # scaling. See studies/scratch/snow_thermal/CALIBRATION.md:92-96.
    xi_v_min       = 10.0 * rho_rat          # keep the storage error under ~10 %
    xi_v_store_err = (rho_rat / xi_v) if xi_v else None
    xi_v_tau       = (Lx * Lx / (xi_v * Dv)) if (xi_v and Lx > 0) else None
    xi_v_warn = None
    if xi_v is not None and xi_v < xi_v_min:
        xi_v_warn = (f"WARNING: ξᵥ = {xi_v:.2e} is below 10·ρ_vs/ρᵢ = "
                     f"{xi_v_min:.2e}. The unscaled vapour STORAGE term is then "
                     f"{100*xi_v_store_err:.1f}% of the source and the "
                     f"quasi-steady approximation degrades. Raise ξᵥ.")

    eps_over_Rave = eps / Rave if Rave > 0 else 0.0
    geom_warn = None
    if eps_over_Rave > 0.10:
        geom_warn = (f"WARNING: ε/R_ave = {eps_over_Rave:.1%} > 10%. "
                     f"Curvature-driven dynamics will be substantially wrong.")
    elif eps_over_Rave > 0.05:
        geom_warn = (f"NOTE: ε/R_ave = {eps_over_Rave:.1%} (5-10%). "
                     f"Borderline Gibbs-Thomson accuracy.")

    # d₀ inflation regime check: β·vₙ (kinetic term) vs d₀·χ (capillary term)
    # in the Gibbs-Thomson condition; χ is the curvature. Unscaled β.
    chi = 1.0 / Rave if Rave > 0 else 0.0
    kin_over_curv = ((beta_uns * v_n) / (d0 * chi)
                     if (d0 * chi) > 0 else float("inf"))

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

    max_delta = max(dpf["delta_heat"], dpf["delta_vapor"])

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
        D_ia=dpf["D_ia"], D_heat=dpf["D_heat"],
        thermal_channel=thermal_channel,
        kinetic_term=dpf["kinetic_term"],
        corr_thermal=dpf["corr_thermal"],
        corr_vapor=dpf["corr_vapor"],
        kinetic_frac=kinetic_frac,
        intractable=intractable,
        delta_heat=dpf["delta_heat"],
        delta_ia=dpf["delta_ia"],
        delta_heat_ice=dpf["delta_heat_ice"],
        delta_heat_air=dpf["delta_heat_air"],
        delta_vapor=dpf["delta_vapor"],
        max_delta=max_delta,
        kin_over_curv=kin_over_curv,
        xi_v_min=xi_v_min, xi_v_store_err=xi_v_store_err,
        xi_v_tau=xi_v_tau, xi_v_warn=xi_v_warn,
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

def patch_opts(opts_path: Path, params: dict, dim: int, T0_C: float) -> None:
    """Write ε, mesh, and the kinetics into opts_path.

    Emits only flags the solver reads. Deliberately does NOT emit -mob_sub or
    -alph_sub: the solver derives both from (ε, β_sub0, d0_sub0) so that
    α_sub/M_sub = 3λ/ε exactly, and overriding either one alone breaks that.

    CONVENTION WARNING -- the two kinetics flags are scaled differently:
      -beta_sub0 is the UNSCALED β₀ = β_HK/(ρ_vs/ρᵢ)
      -d0_sub0   is the RAW physical capillary length γ·V_m/(R·T)
    The solver multiplies BOTH by ρ_vs/ρᵢ, recovering β_HK and d₀·(ρ_vs/ρᵢ).
    Writing a -d0_sub0 taken from M&F's already-scaled d_sub is wrong by ~10⁶.
    """
    text = opts_path.read_text()

    def _replace_or_append(text: str, key: str, value: str, comment: str = "") -> str:
        # Overwrite any existing trailing comment too: for the keys this
        # function owns, the comment is provenance (which temperature and
        # safety factor produced the number), so a stale one is worse than
        # none -- it would describe the previous patch's parameters.
        tail = f"    # {comment}" if comment else ""
        pattern = rf"^{re.escape(key)}\s+\S+[^\S\n]*(?:#.*)?$"
        if re.search(pattern, text, flags=re.MULTILINE):
            return re.sub(pattern, f"{key} {value}{tail}", text, flags=re.MULTILINE)
        return text.rstrip() + f"\n{key} {value}{tail}\n"

    text = _replace_or_append(text, "-eps", f"{params['eps']:.6e}",
                              f"comp_eps.py @ T={T0_C:g}C, safety={params['safety']:g}, "
                              f"binding bound: {params['binding']}")
    text = _replace_or_append(text, "-Nx", str(params["Nx"]),
                              "ceil(sqrt(2)*Lx/eps)")
    if dim >= 2 and params["Ny"] > 0:
        text = _replace_or_append(text, "-Ny", str(params["Ny"]))
    if dim == 3 and params["Nz"] > 0:
        text = _replace_or_append(text, "-Nz", str(params["Nz"]))

    text = _replace_or_append(
        text, "-beta_sub0", f"{params['beta_uns']:.6e}",
        f"K&P beta0 = beta_HK/(rho_vs/rho_i) at {T0_C:g}C, alpha_c={params['alpha_c']:.3e}")
    text = _replace_or_append(
        text, "-d0_sub0", f"{params['d0']:.6e}",
        f"physical gamma*V_m/(R*T) at {T0_C:g}C (K&P Eq. 13)")
    text = _replace_or_append(
        text, "-eps_valid_temp", f"{T0_C:g}",
        "solver ABORTS if -temp differs by >1C")

    opts_path.write_text(text)
    print(f"  Patched {opts_path}")


# =========================================================================
# Output
# =========================================================================

def _print_header(title: str) -> None:
    print("=" * 76)
    print(f"  {title}")
    print("=" * 76)


def eps_over_vn(p: dict) -> float:
    """Time for the interface to advance one ε — the timescale τ_vap must beat."""
    return p["eps"] / p["v_n"] if p["v_n"] else float("inf")


def _print_single(args, p: dict, alpha_c: float, dim: int) -> None:
    _print_header("Phase-field parameters — K&P (2009) + M&F (2024) SI Eq. 9 [sublimation]")

    print(f"\n--- Input parameters ---")
    print(f"  T₀                       = {args.T0:.1f} °C  ({args.T0+273.15:.2f} K)")
    print(f"  α_c (condensation coeff) = {alpha_c:.3e}   [only genuine free parameter]")
    print(f"  vₙ (front velocity)      = {args.vn:.1e}  m/s")
    print(f"  R_ave (grain radius)     = {args.Rave:.2e} m")
    print(f"  Safety factor            = {args.safety:.2f}  (ε = safety · min(bounds); δ = a₁a₂·safety = {_A1*_A2*args.safety:.2f})")
    print(f"  Asymptotic constants     : a₁ = {_A1}, a₂ = {_A2}  (M&F SI footnote 1)")

    print(f"\n--- Thermodynamic state at T₀ ---")
    print(f"  ρ_vs(T₀)         = {p['rho_vs']:.4e}  kg/m³")
    print(f"  ρ_vs/ρᵢ          = {p['rho_rat']:.4e}  [-]")
    print(f"  Dᵥ(T₀)           = {p['Dv']:.4e}  m²/s")
    print(f"  κᵢ/Cᵢ  (D_ice)   = {p['alpha_i']:.4e}  m²/s")
    print(f"  κₐ/Cₐ  (D_air)   = {p['alpha_a']:.4e}  m²/s")
    print(f"  D*_ia (mean)     = {p['D_ia']:.4e}  m²/s")

    print(f"\n--- Physical targets (Gibbs-Thomson parameters) ---")
    print(f"  β_HK (Hertz-Knudsen)   = {p['beta_hk']:.4e}  s/m   [K&P β', SCALED form]")
    print(f"                           [= β_sub·(ρ_vs/ρᵢ) in M&F notation]")
    print(f"  β_sub (physical)       = {p['beta_uns']:.4e}  s/m   [K&P β₀, M&F Table S1 form]")
    print(f"                           [= β_HK/(ρ_vs/ρᵢ); matches K&P range 3e4-3e6]")
    print(f"  d₀_sub  (physical)     = {p['d0']:.4e}  m     [γ·V_m/RT; K&P Eq. 13]")

    print(f"\n--- Upper bounds on ε (K&P Eqs. 43, 45) ---")
    chan = "D*_ia" if p["thermal_channel"] == "mean" else "κᵢ/Cᵢ"
    labels = {
        "B-HEAT":     f"Eq.(43a/b) {chan}·β_HK  ",
        "B-VAPOR":    "Eq.(43c)   Dᵥ·β_HK       ",
        "B-KINETIC":  "Eq.(45)    d₀/(β_sub·vₙ) ",
        "B-CURV":     "Geometric  0.05·R_ave    ",
    }
    for key, label in labels.items():
        marker = "  ← BINDING" if key == p["binding"] else ""
        print(f"  {label}  =  {p[key]:.4e} m{marker}")
    print(f"  (safety applies to the first three; B-CURV is already a fraction)")

    h = p['eps'] / math.sqrt(2.0)
    if p.get("intractable"):
        print(f"\n  *** MESH IS INTRACTABLE: {p['Nx']} x {p.get('Ny',1)} "
              f"= {p['Nx']*max(p.get('Ny',1),1)/1e6:.0f}M nodes ***")
        print(f"      Driven by {p['binding']}. The Eq.43 ceilings scale as β_HK ∝ 1/α_c,")
        print(f"      so this usually means α_c is too LARGE, not that the mesh must be")
        print(f"      this fine. Check --alpha_arrhenius / --alpha_libbrecht2 first.")

    print(f"\n--- Chosen ε and mesh ---")
    print(f"  ε                        = {p['eps']:.4e}  m")
    print(f"  h = ε/√2                 = {h:.4e}  m    (2 elements per W; 8.3 across phi=0.05-0.95)")
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

    if max(p["delta_heat_ice"], p["delta_vapor"]) > 1.0:
        # delta = a1*a2*eps/(D*beta_HK) and beta_HK ~ 1/alpha_c, so every K&P
        # ceiling scales as 1/alpha_c. The neck-fillet criterion does NOT --
        # it is alpha_c-independent -- so overriding eps with the fillet can
        # silently blow past the thin-interface regime when alpha_c is raised.
        DELTA_REF = 0.273          # delta_ice of the validated alpha_c = 1.34e-2 runs
        eps_ref = DELTA_REF * p["alpha_i"] * p["beta_hk"] / (_A1 * _A2)
        print(f"\n  *** delta > 1 ON A CHANNEL ***")
        print(f"      The thin-interface CORRECTION exceeds the kinetic term it")
        print(f"      corrects (delta_ice = {p['delta_heat_ice']:.2f}). beta_HK ~ 1/alpha_c, so every")
        print(f"      K&P ceiling tightens in proportion to alpha_c -- but the")
        print(f"      neck-fillet criterion does not, so overriding eps with the")
        print(f"      fillet can walk straight out of the asymptotic regime.")
        print(f"      At alpha_c = {p['alpha_c']:.1e}, matching the delta_ice = {DELTA_REF} of the")
        print(f"      validated alpha_c = 1.34e-2 runs needs eps <= {eps_ref:.3e} m.")
        print(f"      This run reports eps = {p['eps']:.3e} m ({p['eps']/eps_ref:.1f}x that).")
        print(f"      Decide deliberately; do not let the fillet override hide it.")

    print(f"\n--- Thin-interface expansion parameter δ = a₁a₂·ε/(D·β_HK) ---")
    print(f"  τ_sub compensates δ to first order, so the residual error is O(δ²).")
    print(f"    heat channel ({chan}, active)    δ = {p['delta_heat']:.4f}")
    print(f"    heat-in-ice κᵢ/Cᵢ  [violated]   δ = {p['delta_heat_ice']:.4f}")
    print(f"    heat-in-air κₐ/Cₐ               δ = {p['delta_heat_air']:.4f}")
    print(f"    vapor                           δ = {p['delta_vapor']:.4f}")
    print(f"  τ_sub always compensates at D*_ia (δ = {p['delta_ia']:.4f}), so the ice-side")
    print(f"  correction is under-applied by δ_ice − δ(D*_ia) = {p['delta_heat_ice'] - p['delta_ia']:.3f} — that is the")
    print(f"  honest uncertainty on β. Enforcing the ice channel is intractable;")
    print(f"  M&F violate it by 100x in print for the same reason.")

    print(f"\n--- Kinetic vs. curvature dominance (d₀ inflation regime check) ---")
    print(f"  β_sub·vₙ / (d₀·χ)  with χ ≈ 1/R_ave  = {p['kin_over_curv']:.2e}")
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
                   else ("BORDERLINE (5-10%)" if p["eps_over_Rave"] < 0.10
                         else "VIOLATED (> 10%)"))
    print(f"  ε / R_ave         = {p['eps_over_Rave']:.1%}  [{geom_status}]")
    if p["geom_warn"]:
        print(f"  {p['geom_warn']}")

    print(f"\n--- ξᵥ (M&F temporal scaling) ---")
    print(f"  ξᵥ scales vapour DIFFUSION and the ρᵢ source, NOT the storage term,")
    print(f"  so the residual storage error is ρ_vs/(ξᵥ·ρᵢ): LARGER ξᵥ is MORE")
    print(f"  quasi-steady. Recommended floor 10·ρ_vs/ρᵢ = {p['xi_v_min']:.4e}.")
    if args.xiv is not None:
        print(f"  ξᵥ = {args.xiv:.2e}   storage error = ρ_vs/(ξᵥρᵢ) = "
              f"{100*p['xi_v_store_err']:.3f} %")
        if p["xi_v_tau"]:
            print(f"  τ_vap = Lx²/(ξᵥ·Dᵥ) = {p['xi_v_tau']:.2f} s   "
                  f"(must stay well below ε/vₙ = {eps_over_vn(p):.0f} s)")
        print(f"  {p['xi_v_warn']}" if p["xi_v_warn"] else "  → OK")
    else:
        print(f"  (ξᵥ not provided; pass --xiv to check. Solver default is 1e-3.)")

    print(f"\n--- Recommended opts entries ---")
    print(f"  -eps        {p['eps']:.4e}")
    print(f"  -Nx         {p['Nx']}")
    if dim >= 2 and p["Ny"] > 0: print(f"  -Ny         {p['Ny']}")
    if dim == 3 and p["Nz"] > 0: print(f"  -Nz         {p['Nz']}")
    print(f"  -beta_sub0  {p['beta_uns']:.4e}   # K&P beta0 (unscaled)")
    print(f"  -d0_sub0    {p['d0']:.4e}   # physical gamma*V_m/(R*T)")
    print(f"  -eps_valid_temp {args.T0:g}")
    print(f"  (M_sub/alpha_sub/tau_sub/lam_sub are DERIVED by the solver from")
    print(f"   eps + beta_sub0 + d0_sub0 -- do not set them separately.)")

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
    print("  Note: ε = min(safety · eps_max, eps_over_R · R_ave).  M_sub, α_sub")
    print("  from M&F SI Eq. 9; eps_max is the min of B-HEAT/B-VAPOR/B-KINETIC.")
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
                         "--alpha/--alpha_range. NOTE: at sintering σ (~5e-4) this "
                         "underflows — see --alpha_arrhenius.")
    ap.add_argument("--alpha_arrhenius", action="store_true",
                    help="Derive α_c(T0) from the anchored Arrhenius "
                         "α = A·exp(-Ea/RT) instead of a constant. This is the "
                         "recommended mode: smooth, monotonic, and pinned to the "
                         "α_c ∈ [1e-4, 1e-3] band that Libbrecht (2017) and Braun "
                         "et al. (2024) support. Mutually exclusive with the above.")
    ap.add_argument("--alpha_libbrecht2", action="store_true",
                    help="Derive α_c(T0, σ) from the TWO-parameter Libbrecht fit "
                         "α = A·exp(-f·σ₀(T)/σ), with A and the σ₀ rescaling f "
                         "fitted to the literature anchors. Keeps Libbrecht's "
                         "σ₀(T) shape, and tames the −6/−7 °C table kink to ~8%. "
                         "Pair with --sigma_char. Mutually exclusive with the above.")
    ap.add_argument("--sigma_char", type=float, default=None,
                    help="Characteristic supersaturation σ the α_c fit is pinned "
                         "to [-]. Default: d₀/R_feat (the Gibbs-Thomson σ at the "
                         "smallest resolved feature), using --vn_feature if given "
                         f"else Rave/50. Reference value {LIBBRECHT_SIGMA_CHAR:g}.")
    ap.add_argument("--alpha_anchors", default=None,
                    metavar="T1,a1,T2,a2",
                    help="Arrhenius anchors as 'T_warm,alpha_warm,T_cold,alpha_cold' "
                         f"(default {ALPHA_ANCHOR_WARM[0]:g},{ALPHA_ANCHOR_WARM[1]:g},"
                         f"{ALPHA_ANCHOR_COLD[0]:g},{ALPHA_ANCHOR_COLD[1]:g} → "
                         "Ea = 31.8 kJ/mol).")
    ap.add_argument("--alpha_clamp", default=None, metavar="LO,HI",
                    help=f"Clamp α_c(T) to this band (default "
                         f"{ALPHA_LIT_LO:g},{ALPHA_LIT_HI:g}, the literature range). "
                         "Pass '0,1' to disable and allow free extrapolation.")
    ap.add_argument("--vn",     type=float, default=1.0e-9,
                    help="Normal front velocity vₙ [m/s] for Eq.(45)")
    ap.add_argument("--vn_feature", type=float, default=None, metavar="R_FEAT",
                    help="Derive vₙ self-consistently from capillary-driven "
                         "kinetics: vₙ = d₀/(β_sub·R_FEAT) for the smallest "
                         "curvature feature you must resolve [m]. This makes "
                         "the Eq.(45) bound β-independent (≈ R_FEAT), fixing "
                         "the artifact where colder T (larger β_sub) demanded "
                         "a finer mesh while the physics actually slowed "
                         "down. Overrides --vn.")
    ap.add_argument("--Dchannel", choices=("mean", "ice"), default="mean",
                    help="Diffusivity for the Eq.(43) heat bound and its "
                         "thin-interface check. 'mean' (default) = D*_ia, "
                         "consistent with tau_sub and src/<project>_main.c:558. "
                         "'ice' = kappa_i/C_i, ~6x more restrictive on eps.")
    ap.add_argument("--xiv",    type=float, default=None,
                    help="Time-scaling factor ξᵥ (optional; triggers K&P Eq. 48 check)")
    ap.add_argument("--safety", type=float, default=0.5,
                    help="ε = safety · min(B-HEAT, B-VAPOR, B-KINETIC). Equivalent "
                         "to a thin-interface expansion parameter δ = a₁a₂·safety.")
    ap.add_argument("--eps_over_R", type=float, default=0.05,
                    help="Geometric bound ε <= eps_over_R · R_ave. Crude on "
                         "purpose: 1-2 decades looser than the binding bound.")
    ap.add_argument("--dim",    type=int,   default=None,
                    help="Problem dimension (1/2/3); inferred from Ly/Lz if omitted")
    ap.add_argument("--patch",  type=Path,  default=None,
                    help="Opts file to update in place")
    ap.add_argument("--quiet",  action="store_true",
                    help="Print only ε and Nx/Ny/Nz (one per line)")
    args = ap.parse_args()

    n_modes = sum(x is not None for x in (args.alpha, args.alpha_range, args.sigma_surf))
    n_modes += 1 if args.alpha_arrhenius else 0
    n_modes += 1 if args.alpha_libbrecht2 else 0
    if n_modes > 1:
        ap.error("--alpha, --alpha_range, --sigma_surf, --alpha_arrhenius and "
                 "--alpha_libbrecht2 are mutually exclusive.")
    if n_modes == 0:
        ap.error("Provide one of --alpha, --alpha_range, --sigma_surf, "
                 "--alpha_arrhenius or --alpha_libbrecht2.")

    # ---- alpha_c(T, sigma) from the two-parameter Libbrecht fit. Resolved
    # HERE for the same reason as the Arrhenius below: alpha_c -> beta_sub ->
    # the binding eps ceiling -> the mesh.
    if args.alpha_libbrecht2:
        _say = (lambda *a, **k: None) if args.quiet else print
        aw, ac = ALPHA_ANCHOR_WARM, ALPHA_ANCHOR_COLD
        if args.alpha_anchors:
            v = [float(x) for x in args.alpha_anchors.split(",")]
            if len(v) != 4:
                ap.error("--alpha_anchors needs 'T1,a1,T2,a2'")
            aw, ac = (v[0], v[1]), (v[2], v[3])
        clamp = (ALPHA_LIT_LO, ALPHA_LIT_HI)
        if args.alpha_clamp:
            c = [float(x) for x in args.alpha_clamp.split(",")]
            if len(c) != 2:
                ap.error("--alpha_clamp needs 'LO,HI'")
            clamp = (min(c), max(c))

        rfeat = args.vn_feature if args.vn_feature else args.Rave / 50.0
        sig = args.sigma_char if args.sigma_char else sigma_char_of(args.T0, rfeat)
        A, f, _ = libbrecht2_params(sig, aw, ac)
        raw = alpha_libbrecht2(args.T0, sig, sig, aw, ac, clamp=None)
        args.alpha = alpha_libbrecht2(args.T0, sig, sig, aw, ac, clamp)

        _say(f"\n  Libbrecht 2-parameter mode: alpha_c(T, sigma) = A*exp(-f*sigma0(T)/sigma)")
        _say(f"    anchors  {aw[1]:.3e} @ {aw[0]:g} C  and  {ac[1]:.3e} @ {ac[0]:g} C")
        _say(f"    sigma_char = {sig:.3e}"
              + ("  (given)" if args.sigma_char else f"  (= d0/R_feat, R_feat = {rfeat:.3e} m)"))
        _say(f"    fit: A = {A:.4e}  (rough-surface ceiling),  f = {f:.4e} "
              f"(sigma0 rescaled x1/{1/f:.0f})")
        _say(f"    -> alpha_c({args.T0:g} C) = {args.alpha:.4e}"
              + (f"   [CLAMPED from {raw:.3e} into [{clamp[0]:g}, {clamp[1]:g}]]"
                 if abs(raw - args.alpha) > 1e-12 * max(abs(raw), 1e-300) else ""))
        k6, k7 = alpha_libbrecht2(-6.0, sig, sig, aw, ac, clamp), \
                 alpha_libbrecht2(-7.0, sig, sig, aw, ac, clamp)
        _say(f"    [kink check] alpha(-6 C)/alpha(-7 C) = {k6/k7:.3f} "
              f"(unscaled f=1 this ratio is {math.exp((sigma0(-7)-sigma0(-6))/sig):.2e})")
        _say(f"    NOTE: the clamp floor caps beta_sub at "
              f"{beta_HK(args.T0, clamp[0])/(rho_vs_sat(args.T0)/_RHO_ICE):.3e} s/m. "
              f"Without it beta_sub diverges as sigma -> 0")
        _say(f"      and physical t_final scales as 1/alpha_c (26x at 1e-4 vs 1e-2) "
              f"— the mesh is unaffected, the runtime is not.")

    # ---- Arrhenius alpha_c(T). Resolved HERE, before --vn_feature evaluates
    # beta and before compute_eps sizes the mesh: alpha_c -> beta_sub -> the
    # binding eps bound -> Nx/Ny, so it cannot be applied afterwards.
    if args.alpha_arrhenius:
        _say = (lambda *a, **k: None) if args.quiet else print
        aw, ac = ALPHA_ANCHOR_WARM, ALPHA_ANCHOR_COLD
        if args.alpha_anchors:
            v = [float(x) for x in args.alpha_anchors.split(",")]
            if len(v) != 4:
                ap.error("--alpha_anchors needs 'T1,a1,T2,a2'")
            aw, ac = (v[0], v[1]), (v[2], v[3])
        clamp = (ALPHA_LIT_LO, ALPHA_LIT_HI)
        if args.alpha_clamp:
            c = [float(x) for x in args.alpha_clamp.split(",")]
            if len(c) != 2:
                ap.error("--alpha_clamp needs 'LO,HI'")
            clamp = (c[0], c[1]) if c[1] > c[0] else (c[1], c[0])

        A, Ea = arrhenius_params(aw, ac)
        raw = A * math.exp(-Ea / (_R_GAS * (args.T0 + 273.15)))
        args.alpha = alpha_arrhenius(args.T0, aw, ac, clamp)

        _say(f"\n  Arrhenius kinetics mode: alpha_c(T) = A*exp(-Ea/RT)")
        _say(f"    anchors  {aw[1]:.3e} @ {aw[0]:g} C   and   "
              f"{ac[1]:.3e} @ {ac[0]:g} C")
        _say(f"    Ea = {Ea/1000.0:.2f} kJ/mol,  A = {A:.4e}")
        _say(f"    -> alpha_c({args.T0:g} C) = {args.alpha:.4e}"
              + (f"   [CLAMPED from {raw:.3e} to the "
                 f"[{clamp[0]:g}, {clamp[1]:g}] band]"
                 if abs(raw - args.alpha) > 1e-12 * max(raw, 1e-300) else ""))

        # Show, numerically, why the Libbrecht form is not used directly.
        s_eq = libbrecht_sigma_ref_for(aw[1], aw[0])
        a_cold_lib = alpha_libbrecht(ac[0], s_eq)
        _say(f"    [Libbrecht cross-check] the sigma that reproduces the WARM "
              f"anchor is {s_eq:.3e};")
        _say(f"      at that same sigma the form gives alpha({ac[0]:g} C) = "
              f"{a_cold_lib:.3e}, vs the cold anchor {ac[1]:.3e}.")
        _say(f"      sigma0 spans {sigma0(aw[0]):.2e}->{sigma0(ac[0]):.2e} "
              f"(x{sigma0(ac[0])/sigma0(aw[0]):.0f}), and ln(alpha) scales with "
              f"it, so NO single sigma")
        _say(f"      fits both anchors — the incompatibility is structural, "
              f"not a matter of tuning sigma.")

    if args.sigma_surf is not None:
        s0 = sigma0(args.T0)
        args.alpha = alpha_libbrecht(args.T0, args.sigma_surf)
        print(f"\n  Libbrecht kinetics mode: sigma0({args.T0:g} C) = {s0:.4e}, "
              f"sigma_surf = {args.sigma_surf:.3e}")
        print(f"  -> alpha_c = exp(-sigma0/sigma_surf) = {args.alpha:.4e}"
              + ("   [FLOORED at 1e-30: sigma_surf < sigma0/69.08 — kinetics "
                 "effectively frozen at this supersaturation]"
                 if args.alpha <= 1.0e-30 else ""))

    if args.vn_feature is not None:
        rho_rat_ = rho_vs_sat(args.T0) / _RHO_ICE
        beta_uns_ = beta_HK(args.T0, args.alpha) / rho_rat_ if args.alpha else None
        if beta_uns_ is None:
            ap.error("--vn_feature needs --alpha (or --sigma_surf) to evaluate beta.")
        d0_ = capillary_length(args.T0)
        args.vn = d0_ / (beta_uns_ * args.vn_feature)
        print(f"\n  Self-consistent vn: d0/(beta_sub*R_feat) = {d0_:.3e}/"
              f"({beta_uns_:.3e}*{args.vn_feature:.2e}) = {args.vn:.3e} m/s")

    dim = _infer_dim(args.Lx, args.Ly, args.Lz, args.dim)

    common = dict(
        Lx=args.Lx, Ly=args.Ly, Lz=args.Lz,
        Rave=args.Rave, v_n=args.vn,
        safety=args.safety, xi_v=args.xiv,
        thermal_channel=args.Dchannel,
        eps_over_R=args.eps_over_R,
    )

    if args.alpha is not None:
        p = compute_eps(T0_C=args.T0, alpha_c=args.alpha, **common)

        if args.quiet:
            print(f"{p['eps']:.4e}")
            print(p["Nx"])
            if dim >= 2: print(p["Ny"])
            if dim == 3: print(p["Nz"])
            if args.patch:
                patch_opts(args.patch, p, dim, args.T0)
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
            patch_opts(args.patch, p, dim, args.T0)
        return

    lo, hi, n = float(args.alpha_range[0]), float(args.alpha_range[1]), int(args.alpha_range[2])
    results, alphas = alpha_c_sensitivity(
        T0_C=args.T0, alpha_lo=lo, alpha_hi=hi, n_points=n, **common
    )
    _print_sweep(results, alphas, args.T0)


if __name__ == "__main__":
    _cli()