"""M1 check: recompute every scaling number and dimensionless group in
spec_problem.md (## Mechanisms, ## Scales) from the reference inputs, with Pint
units, and compare to the claimed values within rel_tol (1 %).

Usage: python m1_scales.py [inputs.json]
"""
import json
import math
import os
import sys

import pint

HERE = os.path.dirname(os.path.abspath(__file__))
inp_path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "m1_inputs.json")
with open(inp_path) as f:
    INP = json.load(f)

ureg = pint.UnitRegistry()
Qty = ureg.Quantity
print(json.dumps({"inputs_file": inp_path, "reference_case": INP["reference_case"],
                  "rel_tol": INP["rel_tol"]}))

R = {k: Qty(v[0], v[1]) for k, v in INP["reference_case"].items()}
C = {k: Qty(v[0], v[1] if v[1] else "dimensionless") for k, v in INP["claimed"].items()}
TOL = INP["rel_tol"]
L, k, rho, cp, q, T0 = R["L"], R["k"], R["rho"], R["cp"], R["q"], R["T0"]

results = []


def check(name, claim, computed, claimed, tol=TOL):
    details = {}
    try:
        comp = computed.to(claimed.units)
        details["units_ok"] = True
    except pint.DimensionalityError as e:
        results.append(False)
        print(json.dumps({"check": name, "claim": claim, "pass": False,
                          "details": {"units_ok": False, "error": str(e),
                                      "computed": str(computed), "claimed": str(claimed)}}))
        return
    c, s = float(comp.magnitude), float(claimed.magnitude)
    rel = abs(c - s) / abs(c) if c != 0 else abs(s)
    ok = rel <= tol
    details.update({"computed": c, "claimed": s, "units": str(claimed.units),
                    "rel_diff": rel, "tol": tol})
    results.append(ok)
    print(json.dumps({"check": name, "claim": claim, "pass": ok, "details": details}))


kappa = k / (rho * cp)
check("kappa", "kappa = k/(rho c_p) = 1 m^2/s", kappa, C["kappa"])
check("x_c", "x_c = L = 1 m", L, C["x_c"])
t_c = L**2 / kappa
check("t_c", "t_c = L^2/kappa = 1 s", t_c, C["t_c"])
t_slow = L**2 / (math.pi**2 * kappa)
check("t_slow", "t_slow = L^2/(pi^2 kappa) ~= 0.1013 s", t_slow, C["t_slow"])
check("ln100", "ln 100 ~= 4.605", Qty(math.log(100.0), "dimensionless"), C["ln100"])
t_99 = t_slow * math.log(100.0)
check("t_99", "t_99 = t_slow ln 100 ~= 0.4666 s", t_99, C["t_99"])
# independent: time at which exp(-pi^2 kappa t / L^2) = 0.01
t_99_direct = -math.log(0.01) * L**2 / (math.pi**2 * kappa)
check("t_99_direct", "mode sin(pi x/L) at 1 % of initial amplitude at t ~= 0.4666 s",
      t_99_direct, C["t_99"])
check("tau_steady", "steady state at tau >~ 0.47 (= t_99/t_c)",
      (t_99 / t_c).to("dimensionless"), C["tau_steady"])
dT_q = q * L**2 / k
check("dT_q", "Delta T_q = q L^2/k = 1 K", dT_q, C["dT_q"])
dT_max = q * L**2 / (8 * k)
check("dT_max", "Delta T_max = q L^2/(8k) = 0.125 K", dT_max, C["dT_max"])
check("dTmax_over_T0", "Delta T_max / T_0 = 0.125", (dT_max / T0).to("dimensionless"),
      C["dTmax_over_T0"])
Qg = q * L**2 / (k * T0)
check("Q", "Q = q L^2/(k T_0) = 1 (dimensionless)", Qg, C["Q"])
check("Q_equals_dTq_over_T0", "Q = Delta T_q / T_0", (dT_q / T0).to("dimensionless"), Qg)
check("Q_over_8", "Delta T_max / T_0 = Q/8 = 0.125", Qg / 8, C["Q_over_8"])
check("end_flux", "steady end flux -k dT/dx at x=L equals q L/2 (0.5 W/m^2 ref case)",
      q * L / 2, C["end_flux"])
# Mechanism row 4: "error of order 1e-1 times that power" -> 0.125 is order 1e-1
oom = math.floor(math.log10(float((dT_max / T0).to("dimensionless").magnitude)))
ok = oom == -1
results.append(ok)
print(json.dumps({"check": "order_of_magnitude", "claim":
                  "Delta T_max/T_0 is of order 1e-1 (Mechanisms, T-dependence row)",
                  "pass": ok, "details": {"log10_floor": oom}}))

sys.exit(0 if all(results) else 1)
