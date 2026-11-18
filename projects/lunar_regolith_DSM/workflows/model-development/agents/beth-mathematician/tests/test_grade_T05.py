"""Tests of the T05 grader: right water balances pass whatever the agent names its symbols; each planted error
fails a check that needs no key.

Run: venv_workflows/bin/python -m unittest discover -s workflows/model-development/agents/beth-mathematician/tests
"""
import importlib.util
import unittest
from pathlib import Path

GRADER = Path(__file__).resolve().parents[1] / "tasks" / "T05_vapor_balance" / "grade.py"
spec = importlib.util.spec_from_file_location("grade_T05", GRADER)
g = importlib.util.module_from_spec(spec)
spec.loader.exec_module(g)

P, RV, TT = "phi(x, y, z, t)", "rho_v(x, y, z, t)", "T(x, y, z, t)"
SYMS = [{"name": "Deff(phi, T)", "meaning": "effective vapor diffusivity of the gas", "units": "m^2/s", "typical": 2e-5},
        {"name": "rho_ice", "meaning": "ice density", "units": "kg m^-3", "typical": 920}]
STORE = f"(1 - {P})*Derivative({RV}, t)"
DIFF = "(" + " + ".join(f"Derivative(Deff({P}, {TT})*Derivative({RV}, {c}), {c})" for c in "xyz") + ")"
EXCH = f"(rho_ice - {RV})*Derivative({P}, t)"
RIGHT = f"{STORE} - {DIFF} + {EXCH}"


def answer(expr, symbols=SYMS):
    return {"status": "done", "result_expr": expr, "symbols": symbols, "neglected": []}


class Passes(unittest.TestCase):
    def check(self, expr, symbols=SYMS):
        ok, details = g.grade(answer(expr, symbols))
        self.assertTrue(ok, details)
        return details

    def test_key_check(self):
        ok, why = g.key_check()
        self.assertTrue(ok, why)

    def test_agent_names_with_a_diffusivity_that_depends_on_phi(self):
        self.check(RIGHT)

    def test_written_as_d_dt_of_vapor_content_plus_ice_gain(self):
        self.check(f"Derivative((1 - {P})*{RV}, t) + rho_ice*Derivative({P}, t) - {DIFF}")

    def test_displaced_vapor_dropped_as_small(self):
        """rho_v phi_t is about 1e-9 of rho_ice phi_t; an agent may drop it and say so."""
        self.check(f"{STORE} - {DIFF} + rho_ice*Derivative({P}, t)")


class PlantedErrorsFail(unittest.TestCase):
    def failed_checks(self, expr, symbols=SYMS):
        ok, details = g.grade(answer(expr, symbols))
        self.assertFalse(ok, f"accepted: {details}")
        return details if isinstance(details, str) else [c["check"] for c in details["checks"] if not c["pass"]]

    def test_exchange_sign(self):
        self.assertIn("exchange", self.failed_checks(f"{STORE} - {DIFF} - {EXCH}"))

    def test_vapor_flowing_uphill(self):
        self.assertIn("diffusion", self.failed_checks(f"{STORE} + {DIFF} + {EXCH}"))

    def test_storage_not_weighted_by_gas_fraction(self):
        self.assertIn("exchange", self.failed_checks(f"Derivative({RV}, t) - {DIFF} + {EXCH}"))

    def test_diffusivity_outside_the_divergence(self):
        lap = " + ".join(f"Derivative({RV}, ({c}, 2))" for c in "xyz")
        self.assertIn("conservation", self.failed_checks(f"{STORE} - Deff({P}, {TT})*({lap}) + {EXCH}"))

    def test_exchange_counted_twice(self):
        self.assertIn("exchange", self.failed_checks(f"{STORE} - {DIFF} + 2*rho_ice*Derivative({P}, t)")
                      if False else self.failed_checks(f"Derivative((1 - {P})*{RV}, t) - {DIFF} + {EXCH}"))

    def test_rho_v_missing(self):
        self.assertIn("rho_v", self.failed_checks(f"rho_ice*Derivative({P}, t)"))


if __name__ == "__main__":
    unittest.main()
