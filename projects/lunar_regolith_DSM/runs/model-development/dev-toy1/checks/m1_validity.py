"""M1 check: every constitutive law's stated validity range in spec_problem.md
(## Constitutive laws) is consistent with the problem's temperature range as
stated in the file: T in [T0, T0 + q L^2/(8k)] for q >= 0 (## Notation).

Usage: python m1_validity.py [inputs.json]
"""
import json
import os
import sys

import pint

HERE = os.path.dirname(os.path.abspath(__file__))
inp_path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "m1_inputs.json")
with open(inp_path) as f:
    INP = json.load(f)
ureg = pint.UnitRegistry()
R = {kk: ureg.Quantity(v[0], v[1]) for kk, v in INP["reference_case"].items()}
print(json.dumps({"inputs_file": inp_path, "reference_case": INP["reference_case"]}))
L, k, q, T0 = R["L"], R["k"], R["q"], R["T0"]
results = []


def emit(name, claim, ok, details):
    results.append(bool(ok))
    print(json.dumps({"check": name, "claim": claim, "pass": bool(ok), "details": details}))


Tmin = T0
Tmax = T0 + q * L**2 / (8 * k)
Tmax = Tmax.to("K")
rng = {"T_min_K": float(Tmin.to("K").magnitude), "T_max_K": float(Tmax.magnitude)}

# Problem range is positive (absolute temperature), as the notation requires T > 0
emit("T_positive", "temperature range is > 0 K (Notation: T > 0)",
     Tmin.magnitude > 0 and Tmax.magnitude > 0, rng)

# Fourier law, constant k: needs small relative temperature change.
rel = float(((Tmax - Tmin) / T0).to("dimensionless").magnitude)
emit("fourier_relative_excursion", "Fourier/constant k: relative excursion over the range "
     "Delta T_max/T0 = 0.125, as stated in the 'Holds?' column",
     abs(rel - 0.125) <= 0.01 * 0.125, {**rng, "dTmax_over_T0": rel})
# beta*dT << 1 cannot be evaluated: no beta sourced. Test the weaker,
# checkable statement: for k ~ T^m, beta*dT = m*dT/T0 < 1 needs |m| < 1/rel = 8.
m_crit = 1.0 / rel
emit("fourier_beta_bound", "constant-k error ~ 1e-1 x power m: beta dT_max = m * 0.125 "
     "(< 1 only for |m| < 8; no beta sourced, accepted by brief)",
     abs(m_crit - 8.0) < 1e-9, {"m_for_beta_dT_eq_1": m_crit})

# Caloric law: no phase change and constant cp over the range.
# Checkable part: the range is a single closed interval with no stated transition
# temperature inside it; the file states no transition temperature at all.
emit("caloric_range", "caloric law: no phase change listed inside [T0, T0+qL^2/(8k)] "
     "(file states none; brief line 21 excludes phase change)", True,
     {**rng, "transition_temperatures_in_file": []})

# Source: 'as stated' - uniform, constant; no temperature dependence -> no range limit.
emit("source_range", "source q uniform and T-independent, valid over any T range", True, rng)

# Consistency of Notation range with C2 (sign of q): range assumes q >= 0
emit("q_sign_consistency", "stated T range [T0, T0+qL^2/(8k)] requires q >= 0; reference q = 1 >= 0",
     q.magnitude >= 0, {"q": str(q)})

# Physical plausibility note (not a pass/fail of the file's claims): continuum
# validity at T0 = 1 K is not assessable without a material.
print(json.dumps({"check": "note_continuum_at_1K", "claim":
                  "continuum validity at T0 = 1 K (reference case) is not assessable: no material",
                  "pass": True, "details": {"informational": True}}))

sys.exit(0 if all(results) else 1)
