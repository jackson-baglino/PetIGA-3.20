"""mesh_params: choose the interface parameter eps, the mesh, and the derived phase-field parameters.

This is the agent-facing tool. One job: given the simulation parameters (domain, temperature, condensation
coefficient, ...), return eps, the element counts, the derived model parameters and a list of structured checks.
It is pure and deterministic (no files, no printing); writing input files is the caller's job.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, asdict, field
from typing import Optional

from .constants import CONSTANTS as C, VERSION
from . import physics as P

# Width of the logistic profile phi = 1/(1+exp(-x/eps)) between two levels, in units of eps: 2 ln(p/(1-p)).
BAND_5_95_IN_EPS = 2.0 * math.log(0.95 / 0.05)     # 5.889
BAND_1_99_IN_EPS = 2.0 * math.log(0.99 / 0.01)     # 9.190
H_OVER_EPS = 1.0 / math.sqrt(2.0)                  # mesh rule of this solver: h = eps / sqrt(2)


class ToolInputError(ValueError):
    """Raised with ALL input problems at once, each stating the parameter, its value and its unit."""


@dataclass(frozen=True)
class MeshParamsInput:
    """All lengths in metres, temperature in degrees C, velocity in m/s."""
    domain_x_m: float
    domain_y_m: float = 0.0              # 0 means unused (1D)
    domain_z_m: float = 0.0              # 0 means unused (1D/2D)
    grain_radius_m: float = 3.0e-5       # smallest curvature radius that must be resolved (KP "R_ave")
    temperature_c: float = -20.0
    alpha_c: Optional[float] = None      # condensation coefficient, exactly one of alpha_c / sigma_surf
    sigma_surf: Optional[float] = None   # characteristic supersaturation; alpha_c = exp(-sigma0(T)/sigma_surf)
    front_velocity_m_per_s: float = 1.0e-9
    feature_radius_m: Optional[float] = None   # if set, v_n = d0 / (beta_sub * feature_radius) overrides the velocity
    safety: float = 0.5                  # eps = safety * min(bounds)
    thermal_channel: str = "mean"        # 'mean' (D*_ia, matches the solver) or 'ice' (k_i/C_i, ~6x tighter)
    xi_v: Optional[float] = None         # time-scaling factor; triggers the KP Eq. 48 check
    dimension: Optional[int] = None      # 1, 2 or 3; inferred from the domain if omitted


@dataclass
class Check:
    name: str
    status: str            # "ok" | "note" | "warn" | "fail"
    value: Optional[float]
    limit: Optional[float]
    message: str


def _validate(i: MeshParamsInput) -> None:
    errs = []

    def need(cond, msg):
        if not cond:
            errs.append(msg)

    def finite(x):
        return x is not None and isinstance(x, (int, float)) and math.isfinite(x)

    need(finite(i.domain_x_m) and i.domain_x_m > 0, f"domain_x_m must be > 0 m (got {i.domain_x_m!r})")
    need(finite(i.domain_y_m) and i.domain_y_m >= 0, f"domain_y_m must be >= 0 m (got {i.domain_y_m!r})")
    need(finite(i.domain_z_m) and i.domain_z_m >= 0, f"domain_z_m must be >= 0 m (got {i.domain_z_m!r})")
    need(not (i.domain_z_m > 0 and i.domain_y_m == 0), "domain_z_m > 0 requires domain_y_m > 0")
    need(finite(i.grain_radius_m) and i.grain_radius_m > 0, f"grain_radius_m must be > 0 m (got {i.grain_radius_m!r})")
    need(finite(i.temperature_c) and -150.0 < i.temperature_c < -0.001,
         f"temperature_c must be in (-150, -0.001) degrees C, the range where the ice formulas apply (got {i.temperature_c!r})")
    need((i.alpha_c is None) != (i.sigma_surf is None), "give exactly one of alpha_c and sigma_surf")
    need(i.alpha_c is None or (finite(i.alpha_c) and 0 < i.alpha_c <= 1), f"alpha_c must be in (0, 1] (got {i.alpha_c!r})")
    need(i.sigma_surf is None or (finite(i.sigma_surf) and i.sigma_surf > 0), f"sigma_surf must be > 0 (got {i.sigma_surf!r})")
    need(finite(i.front_velocity_m_per_s) and i.front_velocity_m_per_s > 0,
         f"front_velocity_m_per_s must be > 0 m/s (got {i.front_velocity_m_per_s!r})")
    need(i.feature_radius_m is None or (finite(i.feature_radius_m) and i.feature_radius_m > 0),
         f"feature_radius_m must be > 0 m (got {i.feature_radius_m!r})")
    need(finite(i.safety) and 0 < i.safety <= 1, f"safety must be in (0, 1] (got {i.safety!r})")
    need(i.thermal_channel in ("mean", "ice"), f"thermal_channel must be 'mean' or 'ice' (got {i.thermal_channel!r})")
    need(i.xi_v is None or (finite(i.xi_v) and i.xi_v > 0), f"xi_v must be > 0 (got {i.xi_v!r})")
    need(i.dimension in (None, 1, 2, 3), f"dimension must be 1, 2, 3 or omitted (got {i.dimension!r})")
    if errs:
        raise ToolInputError("; ".join(errs))


def compute_mesh_params(inp: MeshParamsInput) -> dict:
    """Return a JSON-serialisable dict: eps_m, mesh, resolution, eps_bounds_m, derived, thermo, checks, solver_options.

    Raises ToolInputError for invalid input. Does not read or write files.
    """
    _validate(inp)
    T = inp.temperature_c
    alpha_c = inp.alpha_c if inp.alpha_c is not None else P.alpha_libbrecht(T, inp.sigma_surf)

    rho_vs = P.rho_vs_sat_kg_per_m3(T)
    rho_rat = rho_vs / C.rho_ice_kg_per_m3
    beta_hk = P.beta_hk_s_per_m(T, alpha_c)
    beta_uns = beta_hk / rho_rat
    d0 = P.capillary_length_m(T)
    dv = P.dv_m2_per_s(T)
    alpha_i = C.k_ice_w_per_m_k / C.c_ice_j_per_m3_k
    alpha_a = C.k_air_w_per_m_k / C.c_air_j_per_m3_k
    v_n = d0 / (beta_uns * inp.feature_radius_m) if inp.feature_radius_m else inp.front_velocity_m_per_s

    d_heat = P.d_heat_m2_per_s(inp.thermal_channel)
    bounds = {"B-HEAT": d_heat * beta_hk, "B-VAPOR": dv * beta_hk,
              "B-KINETIC": d0 / (beta_uns * v_n), "B-CURV": inp.grain_radius_m}
    eps_max = min(bounds.values())
    binding = min(bounds, key=bounds.get)
    eps = inp.safety * eps_max

    h = eps * H_OVER_EPS
    nx = math.ceil(inp.domain_x_m / h)
    ny = math.ceil(inp.domain_y_m / h) if inp.domain_y_m > 0 else 0
    nz = math.ceil(inp.domain_z_m / h) if inp.domain_z_m > 0 else 0
    dim = inp.dimension or (3 if nz else 2 if ny else 1)

    d = P.derived_pf_params(eps, T, d0, beta_hk, alpha_i, alpha_a, dv, inp.thermal_channel)
    pe = eps * v_n / dv
    kin_over_curv = (beta_uns * v_n) / (d0 / inp.grain_radius_m)
    tot = d["kinetic_term"] + d["corr_thermal"] + d["corr_vapor"]
    kin_frac = d["kinetic_term"] / tot if tot > 0 else 0.0
    min_ratio = min(d["beta_ratio_heat"], d["beta_ratio_vapor"])
    eps_over_r = eps / inp.grain_radius_m

    checks: list[Check] = []
    checks.append(Check("tau_sub_positive", "ok" if d["tau_sub"] > 0 else "fail", d["tau_sub"], 0.0,
                        "tau_sub > 0" if d["tau_sub"] > 0 else "tau_sub <= 0: the asymptotic expansion is ill-posed"))
    st = "warn" if eps_over_r > 0.10 else "note" if eps_over_r > 0.05 else "ok"
    checks.append(Check("eps_over_grain_radius", st, eps_over_r, 0.05,
                        {"warn": "eps/R > 10 %: curvature-driven dynamics will be substantially wrong",
                         "note": "eps/R is 5-10 %: borderline Gibbs-Thomson accuracy", "ok": "eps/R <= 5 %"}[st]))
    st = "warn" if kin_frac < 0.10 else "note" if kin_frac < 0.50 else "ok"
    checks.append(Check("kinetic_fraction_of_tau_sub", st, kin_frac, 0.5,
                        {"warn": "kinetic term is under 10 % of the tau_sub bracket: the sharp-interface analogy is severely stretched",
                         "note": "kinetic term is 10-50 % of the tau_sub bracket: marginal sharp-interface regime",
                         "ok": "kinetic term dominates tau_sub"}[st]))
    if inp.xi_v is not None:
        ok = inp.xi_v <= rho_rat
        checks.append(Check("xi_v_quasi_steady", "ok" if ok else "warn", inp.xi_v, rho_rat,
                            "xi_v <= rho_vs/rho_i (KP Eq. 48)" if ok else "xi_v exceeds rho_vs/rho_i: quasi-steady vapor assumption violated (KP Eq. 48)"))
    checks.append(Check("thin_interface_ratio", "ok" if min_ratio >= 0.9 else "note", min_ratio, 0.9,
                        "beta_eff/beta_target >= 0.9" if min_ratio >= 0.9 else
                        "beta_eff/beta_target < 0.9. This is expected with safety = 0.5 (about 0.605); the 0.9 reading is about 4x tighter in eps"))
    d0_code = P.capillary_length_code_m(T)
    mism = abs(d0_code - d0) / d0
    checks.append(Check("d0_matches_solver", "ok" if mism <= 1e-3 else "note", mism, 1e-3,
                        "d0 agrees with monitoring.c" if mism <= 1e-3 else
                        f"monitoring.c hard-codes d0 = 2.548e-7/T_K, {100*(d0_code/d0-1):+.1f} % vs this tool; reconcile before quantitative validation"))
    checks.append(Check("thin_interface_two_sided_constants", "note", C.a2, None,
                        "a1, a2 are the Moure-Fu constants derived by analogy to the TWO-SIDED Karma-Rappel analysis; the solver's vapor "
                        "diffusivity is one-sided and has no anti-trapping current (see NOTES.md). tau_sub and beta_ratio are therefore "
                        "estimates, not controlled results"))
    checks.append(Check("vapor_diffusivity_is_air", "note", dv, None,
                        "Dv is the diffusivity of water vapor in AIR at 1 atm. For lunar pores without background gas the regime is "
                        "Knudsen/free-molecular and eps_bounds B-VAPOR/Pe are not meaningful; see the transport-law review"))

    cp = [nx + 2] + ([ny + 2] if ny else []) + ([nz + 2] if nz else [])
    ncp = math.prod(cp)
    return {
        "tool": "mesh_params", "version": VERSION,
        "inputs": {**asdict(inp), "alpha_c_used": alpha_c, "front_velocity_used_m_per_s": v_n, "dimension_used": dim},
        "eps_m": eps,
        "mesh": {"dimension": dim, "h_m": h, "nx": nx, "ny": ny, "nz": nz,
                 "elements_total": math.prod(n for n in (nx, ny, nz) if n),
                 "control_points_p2_c1": ncp, "dofs_two_phase_3dof": 3 * ncp},
        "resolution": {"band_5_95_m": BAND_5_95_IN_EPS * eps, "band_1_99_m": BAND_1_99_IN_EPS * eps,
                       "elements_across_5_95": BAND_5_95_IN_EPS * eps / h, "elements_across_1_99": BAND_1_99_IN_EPS * eps / h,
                       "note": "eps is a decay-length scale; judge resolution by the visible 1-99 % band, not by eps"},
        "eps_bounds_m": {**bounds, "eps_max_m": eps_max, "binding": binding},
        "derived": {"tau_sub_s": d["tau_sub"], "lam_sub": d["lam_sub"], "M_sub": d["M_sub"], "alpha_src_per_s": d["alpha_src"],
                    "kinetic_term": d["kinetic_term"], "corr_thermal": d["corr_thermal"], "corr_vapor": d["corr_vapor"],
                    "kinetic_fraction": kin_frac, "min_beta_ratio": min_ratio,
                    "beta_ratio_heat": d["beta_ratio_heat"], "beta_ratio_vapor": d["beta_ratio_vapor"]},
        "thermo": {"rho_vs_kg_per_m3": rho_vs, "rho_vs_over_rho_ice": rho_rat, "dv_m2_per_s": dv, "d_heat_m2_per_s": d_heat,
                   "beta_hk_s_per_m": beta_hk, "beta_sub_s_per_m": beta_uns, "d0_m": d0, "peclet": pe, "kinetic_over_curvature": kin_over_curv},
        "checks": [asdict(c) for c in checks],
        "solver_options": solver_options(eps, nx, ny, nz, dim),
        "provenance": {"formulas": "Kaempfer and Plapp 2009 Eqs. 13, 43-48; Moure and Fu 2024 SI Eq. 9",
                       "mesh_rule": "h = eps/sqrt(2)", "constants_version": VERSION},
    }


def solver_options(eps: float, nx: int, ny: int, nz: int, dim: int) -> dict:
    """The only options this tool is allowed to set in a solver opts file: -eps and the element counts."""
    o = {"-eps": f"{eps:.4e}", "-Nx": str(nx)}
    if dim >= 2 and ny:
        o["-Ny"] = str(ny)
    if dim == 3 and nz:
        o["-Nz"] = str(nz)
    return o


def sweep_alpha_c(base: MeshParamsInput, alpha_lo: float, alpha_hi: float, n_points: int = 10) -> list:
    """Log-uniform sweep of alpha_c; returns [(alpha_c, result_dict), ...]. base.alpha_c / sigma_surf are ignored."""
    if not (0 < alpha_lo < alpha_hi <= 1) or n_points < 2:
        raise ToolInputError("need 0 < alpha_lo < alpha_hi <= 1 and n_points >= 2")
    out = []
    for k in range(n_points):
        a = math.exp(math.log(alpha_lo) + k * math.log(alpha_hi / alpha_lo) / (n_points - 1))
        out.append((a, compute_mesh_params(MeshParamsInput(**{**asdict(base), "alpha_c": a, "sigma_surf": None}))))
    return out


TOOL_SPEC = {
    "name": "mesh_params",
    "description": (
        "Compute the diffuse-interface parameter eps, the element counts and the derived phase-field parameters for the "
        "lunar_regolith_DSM solver from the simulation parameters. Units are in the argument names (m, degrees C, m/s). "
        "Give exactly one of alpha_c (condensation coefficient) or sigma_surf. Returns eps_m, mesh (h, nx, ny, nz), "
        "resolution (elements across the 5-95 % and 1-99 % bands), eps bounds with the binding one, derived tau_sub/lam_sub/M_sub, "
        "and a list of checks with status ok/note/warn/fail. Never hand-tune eps: use this tool. Does not write files."),
    "input_schema": {
        "type": "object",
        "properties": {
            "domain_x_m": {"type": "number", "exclusiveMinimum": 0, "description": "Domain length in x [m]"},
            "domain_y_m": {"type": "number", "minimum": 0, "description": "Domain length in y [m]; 0 = unused"},
            "domain_z_m": {"type": "number", "minimum": 0, "description": "Domain length in z [m]; 0 = unused"},
            "grain_radius_m": {"type": "number", "exclusiveMinimum": 0, "description": "Smallest curvature radius to resolve [m]"},
            "temperature_c": {"type": "number", "description": "Mean temperature [degrees C], between -150 and 0"},
            "alpha_c": {"type": "number", "exclusiveMinimum": 0, "maximum": 1, "description": "Condensation coefficient [-]"},
            "sigma_surf": {"type": "number", "exclusiveMinimum": 0, "description": "Characteristic supersaturation [-]; alternative to alpha_c"},
            "front_velocity_m_per_s": {"type": "number", "exclusiveMinimum": 0, "description": "Normal front velocity v_n [m/s] for KP Eq. 45"},
            "feature_radius_m": {"type": "number", "exclusiveMinimum": 0, "description": "If given, v_n = d0/(beta_sub*feature_radius_m) overrides the velocity"},
            "safety": {"type": "number", "exclusiveMinimum": 0, "maximum": 1, "description": "eps = safety * min(bounds)"},
            "thermal_channel": {"type": "string", "enum": ["mean", "ice"], "description": "Heat-channel diffusivity: mean = D*_ia (solver), ice = k_i/C_i"},
            "xi_v": {"type": "number", "exclusiveMinimum": 0, "description": "Time-scaling factor; enables the KP Eq. 48 check"},
            "dimension": {"type": "integer", "enum": [1, 2, 3], "description": "Problem dimension; inferred if omitted"},
        },
        "required": ["domain_x_m"],
    },
}
