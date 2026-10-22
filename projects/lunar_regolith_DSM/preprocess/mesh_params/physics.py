"""Pure formulas (no I/O, no printing). Ported unchanged from preprocess/comp_eps.py so numbers match exactly.

References: Kaempfer and Plapp, Phys. Rev. E 79, 031502 (2009) [KP]; Moure and Fu, Cryst. Growth Des. 24 (2024) [MF].
Temperatures are in degrees C unless a name ends in _k.
"""
from __future__ import annotations

import math

from .constants import CONSTANTS as C


def rho_vs_sat_kg_per_m3(temperature_c: float) -> float:
    """Saturation vapor density over ice (ASHRAE polynomial), as in material_properties.c."""
    t = temperature_c + 273.15
    k = C.ashrae_kj
    p_vs = math.exp(k[0] / t + k[1] + k[2] * t + k[3] * t**2 + k[4] * t**3 + k[5] * math.log(t))
    return C.rho_air_kg_per_m3 * C.bb * p_vs / (C.p_atm_pa - p_vs)


def dv_m2_per_s(temperature_c: float) -> float:
    """Vapor diffusivity in AIR at 1 atm, Dv0 (T/273.15)^1.81. Not valid for pores without background gas."""
    return C.dv0_m2_per_s * ((temperature_c + 273.15) / 273.15) ** 1.81


def beta_hk_s_per_m(temperature_c: float, alpha_c: float) -> float:
    """Hertz-Knudsen kinetic coefficient, SCALED form: KP beta', equal to MF beta_sub * (rho_vs/rho_i)."""
    t = temperature_c + 273.15
    return (1.0 / alpha_c) * math.sqrt(2.0 * math.pi * C.m_h2o_kg / (C.k_b_j_per_k * t))


def capillary_length_m(temperature_c: float) -> float:
    """d0 = gamma Vm / (R T) (KP Eq. 13). This tool's d0; sizes eps."""
    return C.gamma_iv_j_per_m2 * C.vm_ice_m3_per_mol / (C.r_gas_j_per_mol_k * (temperature_c + 273.15))


def capillary_length_solver_m(temperature_c: float) -> float:
    """d0 as the solver computes it once at temp0: D0_CAPILLARY = gamma a^3/(k_B T), a^3 = m_H2O/rho_ice."""
    a3 = C.m_h2o_d0_kg / C.rho_ice_kg_per_m3
    return C.gamma_iv_j_per_m2 * a3 / (C.k_b_d0_j_per_k * (temperature_c + 273.15))


def sigma0(temperature_c: float) -> float:
    """Libbrecht critical supersaturation, exact port of Sigma0() in material_properties.c (log-log table)."""
    s, tt = C.sigma0_table_s, C.sigma0_table_t_c
    if temperature_c > tt[0]:
        return s[0]
    interv = 0
    for ii in range(10):
        if temperature_c <= tt[ii]:
            interv = ii
    if interv == 9:
        return s[9]
    t0, t1 = abs(tt[interv]), abs(tt[interv + 1])
    s0, s1 = s[interv], s[interv + 1]
    return 10.0 ** (math.log10(s0) + (math.log10(s1) - math.log10(s0)) / (math.log10(t1) - math.log10(t0))
                    * (math.log10(abs(temperature_c)) - math.log10(t0)))


def alpha_libbrecht(temperature_c: float, sigma_surf: float) -> float:
    """alpha(T, sigma) = exp(-sigma0(T)/sigma), floored at 1e-30 like monitoring.c."""
    s0 = sigma0(temperature_c)
    if sigma_surf < s0 / 69.0775:
        return 1.0e-30
    return math.exp(-s0 / sigma_surf)


def d_heat_m2_per_s(thermal_channel: str = "mean") -> float:
    """Diffusivity for the KP Eq. 43 heat channel: 'mean' = D*_ia = (k_i/C_i + k_a/C_a)/2, 'ice' = k_i/C_i."""
    d_i = C.k_ice_w_per_m_k / C.c_ice_j_per_m3_k
    d_a = C.k_air_w_per_m_k / C.c_air_j_per_m3_k
    if thermal_channel == "mean":
        return 0.5 * (d_i + d_a)
    if thermal_channel == "ice":
        return d_i
    raise ValueError(f"thermal_channel must be 'mean' or 'ice', got {thermal_channel!r}")


def derived_pf_params(eps_m, temperature_c, d0_m, beta_hk, alpha_i, alpha_a, dv, thermal_channel="mean") -> dict:
    """Solve MF SI Eq. 9 for lambda_sub, tau_sub, M_sub, alpha_sub given eps and the physical targets."""
    rho_rat = rho_vs_sat_kg_per_m3(temperature_c) / C.rho_ice_kg_per_m3
    d_ia = 0.5 * (alpha_i + alpha_a)
    lam = C.a1 * eps_m / (d0_m * rho_rat)
    kinetic = beta_hk / C.a1
    corr_th = C.a2 * eps_m / d_ia
    corr_v = C.a2 * eps_m / dv
    tau = eps_m * lam * (kinetic + corr_th + corr_v)

    def ratio(d_channel):
        return 1.0 - C.a1 * C.a2 * eps_m / (d_channel * beta_hk)

    d_heat = d_heat_m2_per_s(thermal_channel)
    return dict(lam_sub=lam, tau_sub=tau, M_sub=eps_m / (3.0 * tau), alpha_src=lam / tau, D_ia=d_ia, D_heat=d_heat,
                kinetic_term=kinetic, corr_thermal=corr_th, corr_vapor=corr_v,
                beta_ratio_heat=ratio(d_heat), beta_ratio_heat_ice=ratio(alpha_i),
                beta_ratio_heat_air=ratio(alpha_a), beta_ratio_vapor=ratio(dv))
