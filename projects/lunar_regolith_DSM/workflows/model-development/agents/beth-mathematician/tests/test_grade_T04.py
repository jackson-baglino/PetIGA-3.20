"""Tests of the T04 grader: right energy equations pass whatever the agent names its symbols; each planted
physical error fails a check that needs no key; terms the current model lacks are listed, not failed.

Run: venv_workflows/bin/python -m unittest discover -s workflows/model-development/agents/beth-mathematician/tests
"""
import importlib.util
import unittest
from pathlib import Path

GRADER = Path(__file__).resolve().parents[1] / "tasks" / "T04_heat_latent" / "grade.py"
spec = importlib.util.spec_from_file_location("grade_T04", GRADER)
g = importlib.util.module_from_spec(spec)
spec.loader.exec_module(g)

P, TT = "phi(x, y, z, t)", "T(x, y, z, t)"
# The agent's own names: a volumetric heat capacity C(phi, T), a conductivity lam(phi, T), rho_ice and Ls.
SYMS = [{"name": "C(phi, T)", "meaning": "volumetric heat capacity", "units": "J m^-3 K^-1", "typical": 1.5e6},
        {"name": "lam(phi, T)", "meaning": "thermal conductivity", "units": "W m^-1 K^-1", "typical": 2},
        {"name": "rho_ice", "meaning": "ice density", "units": "kg m^-3", "typical": 920},
        {"name": "Ls", "meaning": "latent heat of sublimation", "units": "J/kg", "typical": 2.84e6}]
STORE = f"C({P}, {TT})*Derivative({TT}, t)"
COND = "(" + " + ".join(f"Derivative(lam({P}, {TT})*Derivative({TT}, {c}), {c})" for c in "xyz") + ")"
LAT = f"rho_ice*Ls*Derivative({P}, t)"
RIGHT = f"{STORE} - {COND} - {LAT}"


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

    def test_agent_names_and_T_dependent_properties(self):
        d = self.check(RIGHT)
        self.assertEqual(d["comparison (for the user)"]["agent's terms by role"]["other"], [])

    def test_times_minus_one(self):
        self.check(f"-({RIGHT})")

    def test_conduction_written_out_by_the_product_rule(self):
        expanded = " + ".join(f"lam({P}, {TT})*Derivative({TT}, ({c}, 2)) + Derivative(lam({P}, {TT}), {c})"
                              f"*Derivative({TT}, {c})" for c in "xyz")
        self.check(f"{STORE} - ({expanded}) - {LAT}")

    def test_extra_term_is_listed_not_failed(self):
        """A vapor-density correction to the latent term: not in the current model, for the user to judge."""
        syms = SYMS + [{"name": "rho_v", "meaning": "vapor density", "units": "kg m^-3", "typical": 1e-3}]
        d = self.check(f"{STORE} - {COND} - (rho_ice - rho_v)*Ls*Derivative({P}, t)", syms)
        self.assertIn("phase change", d["comparison (for the user)"]["agent's terms by role"])


class PlantedErrorsFail(unittest.TestCase):
    def failed_checks(self, expr, symbols=SYMS):
        ok, details = g.grade(answer(expr, symbols))
        self.assertFalse(ok, f"accepted: {details}")
        return details if isinstance(details, str) else [c["check"] for c in details["checks"] if not c["pass"]]

    def test_latent_sign(self):
        self.assertIn("latent sign", self.failed_checks(f"{STORE} - {COND} + {LAT}"))

    def test_heat_flowing_uphill(self):
        self.assertIn("conduction", self.failed_checks(f"{STORE} + {COND} - {LAT}"))

    def test_conductivity_outside_the_divergence(self):
        lap = " + ".join(f"Derivative({TT}, ({c}, 2))" for c in "xyz")
        self.assertIn("conservation", self.failed_checks(f"{STORE} - lam({P}, {TT})*({lap}) - {LAT}"))

    def test_latent_heat_per_volume_missing_density(self):
        self.assertIn("units", self.failed_checks(f"{STORE} - {COND} - Ls*Derivative({P}, t)"))

    def test_wrong_declared_units(self):
        syms = [dict(s, units="J kg^-1 K^-1") if s["name"].startswith("C(") else s for s in SYMS]
        self.assertIn("units", self.failed_checks(RIGHT, syms))

    def test_source_out_of_nothing(self):
        syms = SYMS + [{"name": "q", "meaning": "heat source", "units": "W m^-3", "typical": 100}]
        self.assertIn("equilibrium", self.failed_checks(f"{RIGHT} - q", syms))

    def test_undeclared_name(self):
        self.assertIn("neither in the work order nor declared", self.failed_checks(f"{RIGHT} - k0"))

    def test_code_refused(self):
        self.assertIn("characters", self.failed_checks("__import__('os').system('true')"))

    def test_latent_term_missing(self):
        self.assertIn("latent sign", self.failed_checks(f"{STORE} - {COND}"))

    def test_missing_typical_value(self):
        syms = [{k: v for k, v in s.items() if k != "typical"} for s in SYMS]
        self.assertIn("typical", self.failed_checks(RIGHT, syms))


if __name__ == "__main__":
    unittest.main()
