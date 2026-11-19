"""Tests of the T07 grader: right phase-field equations pass in other forms; each planted error fails a check
that needs no key.

Run: venv_workflows/bin/python -m unittest discover -s workflows/model-development/agents/beth-mathematician/tests
"""
import importlib.util
import unittest
from pathlib import Path

GRADER = Path(__file__).resolve().parents[1] / "tasks" / "T07_phase_field" / "grade.py"
spec = importlib.util.spec_from_file_location("grade_T07", GRADER)
g = importlib.util.module_from_spec(spec)
spec.loader.exec_module(g)

P, TT, RV = "phi(x, y, z, t)", "T(x, y, z, t)", "rho_v(x, y, z, t)"
LAP = "(" + " + ".join(f"Derivative({P}, ({c}, 2))" for c in "xyz") + ")"
MU = f"(6*gamma/eps*{P}*(1 - {P})*(1 - 2*{P}) - 6*gamma*eps*{LAP})"
SYMS = [{"name": "Lphi(T)", "meaning": "Allen-Cahn mobility", "units": "m^3 J^-1 s^-1", "typical": 1e-3},
        {"name": "beta(T)", "meaning": "kinetic coefficient", "units": "m s^-1", "typical": 1e-2},
        {"name": "rho_ice", "meaning": "ice density", "units": "kg m^-3", "typical": 920}]
DRIVE = f"beta({TT})/eps*30*{P}**2*(1 - {P})**2*({RV} - rho_vs({TT}))/rho_ice"
RIGHT = f"Derivative({P}, t) + Lphi({TT})*{MU} - {DRIVE}"


def failed(expr, symbols=SYMS):
    ok, details = g.grade({"status": "done", "result_expr": expr, "symbols": list(symbols), "neglected": []})
    return ok, details if isinstance(details, str) else [c["check"] for c in details["checks"] if not c["pass"]]


class Passes(unittest.TestCase):
    def test_key_check(self):
        ok, why = g.key_check()
        self.assertTrue(ok, why)

    def test_mobility_times_variational_derivative(self):
        self.assertEqual(failed(RIGHT), (True, []))

    def test_relative_supersaturation_driving_force(self):
        drive = f"beta({TT})/eps*{P}**2*(1 - {P})**2*({RV}/rho_vs({TT}) - 1)"
        self.assertEqual(failed(f"Derivative({P}, t) + Lphi({TT})*{MU} - {drive}"), (True, []))


    def test_unused_symbol_with_typical_zero(self):
        """Found in run 20261118-162313: the agent declared the residual R with typical value 0. A symbol the
        result does not use needs no size."""
        syms = SYMS + [{"name": "R", "meaning": "residual", "units": "s^-1", "typical": 0}]
        self.assertEqual(failed(RIGHT, syms), (True, []))

    def test_used_symbol_with_typical_zero_refused(self):
        syms = [dict(s, typical=0) if s["name"] == "rho_ice" else s for s in SYMS]
        ok, why = failed(RIGHT, syms)
        self.assertFalse(ok)


class PlantedErrorsFail(unittest.TestCase):
    def test_nucleation_in_gas(self):
        drive = f"beta({TT})/eps*({RV} - rho_vs({TT}))/rho_ice"
        self.assertIn("bulk", failed(f"Derivative({P}, t) + Lphi({TT})*{MU} - {drive}")[1])

    def test_driving_force_sign(self):
        self.assertIn("growth", failed(f"Derivative({P}, t) + Lphi({TT})*{MU} + {DRIVE}")[1])

    def test_curvature_not_from_F(self):
        wrong_mu = f"(6*gamma/eps*{P}*(1 - {P})*(1 - 2*{P}) - gamma*eps*{LAP})"
        self.assertIn("relaxation", failed(f"Derivative({P}, t) + Lphi({TT})*{wrong_mu} - {DRIVE}")[1])

    def test_uphill(self):
        self.assertIn("relaxation", failed(f"Derivative({P}, t) - Lphi({TT})*{MU} - {DRIVE}")[1])


if __name__ == "__main__":
    unittest.main()
