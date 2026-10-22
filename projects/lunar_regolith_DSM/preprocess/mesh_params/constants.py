"""Physical constants and model constants used by the mesh-parameter tool.

Single source of truth for this tool. Every value that is duplicated in the C solver is listed in
C_SOURCE_LINKS so a test can detect drift (preprocess/tests/test_mesh_params.py::test_constants_match_solver_source).
Values are unchanged from the original preprocess/comp_eps.py.
"""
from __future__ import annotations

from dataclasses import dataclass

VERSION = "1.0.0"


@dataclass(frozen=True)
class Constants:
    m_h2o_kg: float = 2.99e-26          # water molecule mass
    k_b_j_per_k: float = 1.38e-23       # Boltzmann constant
    k_ice_w_per_m_k: float = 2.29       # thermal conductivity of ice
    k_air_w_per_m_k: float = 0.02       # thermal conductivity of air
    c_ice_j_per_m3_k: float = 1.8e6     # volumetric heat capacity of ice
    c_air_j_per_m3_k: float = 1.341 * 1004.0
    rho_ice_kg_per_m3: float = 919.0
    dv0_m2_per_s: float = 2.178e-5      # vapor diffusivity in air at 273.15 K
    p_atm_pa: float = 101325.0
    rho_air_kg_per_m3: float = 1.341
    bb: float = 0.62                    # ASHRAE humidity-ratio coefficient
    gamma_iv_j_per_m2: float = 0.109    # ice-vapor surface energy (Kaempfer and Plapp 2009, Table I)
    vm_ice_m3_per_mol: float = 1.963e-5
    r_gas_j_per_mol_k: float = 8.314
    a1: float = 5.0                     # Moure and Fu SI footnote 1 (phi in [0,1] double well)
    a2: float = 0.1581                  # idem; a Karma-Rappel-type constant for the TWO-SIDED problem
    ashrae_kj: tuple = (-0.5865e4, 0.2224e2, 0.1375e-1, -0.3403e-4, 0.2697e-7, 0.6918)
    sigma0_table_s: tuple = (3.0e-3, 4.1e-3, 5.5e-3, 8.0e-3, 4.0e-3, 6.0e-3, 3.5e-2, 7.0e-2, 1.1e-1, 0.75)
    sigma0_table_t_c: tuple = (-0.0001, -2.0, -4.0, -6.0, -7.0, -10.0, -20.0, -30.0, -40.0, -100.0)


CONSTANTS = Constants()

# (constant name, value, file, text that must appear in that file). Checked by the drift test.
# The molecule mass differs slightly in monitoring.c (3.0e-26 vs 2.99e-26) and is listed as a known difference.
C_SOURCE_LINKS = [
    ("k_ice_w_per_m_k", 2.29, "src/lunar_main.c", "thcond_ice = 2.29"),
    ("k_air_w_per_m_k", 0.02, "src/lunar_main.c", "thcond_air = 0.02"),
    ("rho_ice_kg_per_m3", 919.0, "src/lunar_main.c", "rho_ice    = 919.0"),
    ("dv0_m2_per_s", 2.178e-5, "src/lunar_main.c", "dif_vap    = 2.178e-5"),
    ("rho_air_kg_per_m3", 1.341, "src/lunar_main.c", "rho_air    = 1.341"),
    ("p_atm_pa", 101325.0, "src/material_properties.c", "Patm = 101325.0"),
    ("bb", 0.62, "src/material_properties.c", "bb = 0.62"),
    ("a1", 5.0, "src/lunar_main.c", "a1 = 5.0"),
    ("a2", 0.1581, "src/lunar_main.c", "a2 = 0.1581"),
]
KNOWN_DIFFERENCES = {
    "m_h2o_kg": "monitoring.c uses 3.0e-26 kg in the Hertz-Knudsen kinetic term (0.3 % above this value)",
    "gamma_iv_j_per_m2": "monitoring.c hard-codes d0 = 2.548e-7/T_K, about 1 % below gamma*Vm/(R*T) with this gamma",
}
