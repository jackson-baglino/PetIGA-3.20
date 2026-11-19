"""M1 check (Pint): every unit in the ## Notation table of spec_problem.md is
consistent with the formulas that use it, and every term in each equation row
has the same units.

Usage: python m1_units.py
"""
import json
import sys

import pint

ureg = pint.UnitRegistry()
U = ureg.Unit
# Notation table, transcribed from spec_problem.md ## Notation
NOT = {
    "x": "m", "t": "s", "L": "m", "T": "K", "T0": "K",
    "k": "W/(m*K)", "rho": "kg/m**3", "cp": "J/(kg*K)", "kappa": "m**2/s",
    "q": "W/m**3", "j": "W/m**2", "dTq": "K",
    "xi": "dimensionless", "tau": "dimensionless", "theta": "dimensionless",
    "Q": "dimensionless", "n": "dimensionless",
}
print(json.dumps({"notation_units": NOT}))
u = {s: (1 * U(v)) for s, v in NOT.items()}
results = []


def same(name, claim, a, b):
    ok = a.dimensionality == b.dimensionality
    results.append(ok)
    print(json.dumps({"check": name, "claim": claim, "pass": ok,
                      "details": {"lhs": str(a.dimensionality), "rhs": str(b.dimensionality)}}))


same("kappa_def", "kappa = k/(rho c_p) has units m^2/s", u["k"] / (u["rho"] * u["cp"]), u["kappa"])
same("flux_def", "j = -k dT/dx has units W/m^2", u["k"] * u["T"] / u["x"], u["j"])
same("dTq_def", "Delta T_q = q L^2/k has units K", u["q"] * u["L"]**2 / u["k"], u["dTq"])
same("Q_def", "Q = q L^2/(k T0) is dimensionless", u["q"] * u["L"]**2 / (u["k"] * u["T0"]), u["Q"])
same("xi_def", "xi = x/L is dimensionless", u["x"] / u["L"], u["xi"])
same("tau_def", "tau = kappa t/L^2 is dimensionless", u["kappa"] * u["t"] / u["L"]**2, u["tau"])
same("theta_def", "theta = (T-T0)/Delta T_q is dimensionless", (u["T"] - u["T0"]) / u["dTq"], u["theta"])
# heat equation row: rho cp dT/dt, k d2T/dx2, q all W/m^3
storage = u["rho"] * u["cp"] * u["T"] / u["t"]
cond = u["k"] * u["T"] / u["x"]**2
same("pde_storage_vs_source", "rho c_p dT/dt has units of q (W/m^3)", storage, u["q"])
same("pde_conduction_vs_source", "k d2T/dx2 has units of q (W/m^3)", cond, u["q"])
same("flux_div_vs_source", "dj/dx has units of q", u["j"] / u["x"], u["q"])
# caloric law e = rho cp (T - Tref): energy per volume
same("caloric", "e = rho c_p (T - T_ref) has units J/m^3", u["rho"] * u["cp"] * u["T"], 1 * U("J/m**3"))
same("caloric_rate", "de/dt has units of q", u["rho"] * u["cp"] * u["T"] / u["t"], u["q"])
# scales
same("t_c", "t_c = L^2/kappa has units s", u["L"]**2 / u["kappa"], u["t"])
same("t_slow", "t_slow = L^2/(pi^2 kappa) has units s", u["L"]**2 / u["kappa"], u["t"])
same("dTmax", "q L^2/(8k) has units K", u["q"] * u["L"]**2 / u["k"], u["T"])
same("parabola", "q x (L-x)/(2k) has units K", u["q"] * u["x"] * u["L"] / u["k"], u["T"])
same("end_flux", "q L/2 has units of flux j (W/m^2)", u["q"] * u["L"], u["j"])
same("balance", "k DeltaT/L^2 has units of q", u["k"] * u["dTq"] / u["L"]**2, u["q"])
same("beta_dT", "beta Delta T_max is dimensionless with beta = d ln k/dT in 1/K",
     (1 / U("K")) * u["dTq"], 1 * U("dimensionless"))
same("Q_form_source", "in theta'=T/T0 form, source term Q matches theta'_tau units",
     u["Q"], u["theta"] / u["tau"])

sys.exit(0 if all(results) else 1)
