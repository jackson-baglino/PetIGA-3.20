"""Tests of the T02 grader: a correct answer passes in several forms; each planted error fails.

Run: venv_workflows/bin/python -m unittest discover -s workflows/model-development/agents/beth-mathematician/tests
"""
import importlib.util
import unittest
from pathlib import Path

GRADER = Path(__file__).resolve().parents[1] / "tasks" / "T02_heat_nondimensionalize" / "grade.py"
spec = importlib.util.spec_from_file_location("grade_T02", GRADER)
g = importlib.util.module_from_spec(spec)
spec.loader.exec_module(g)

TS = "Derivative(theta(xi, s), s)"
TXX = "Derivative(theta(xi, s), (xi, 2))"
GROUP = "Q0*L**2/(k0*DT)"
CORRECT = f"{TS} - {TXX} - {GROUP}"


class Passes(unittest.TestCase):
    def check(self, text):
        ok, details = g.grade(text)
        self.assertTrue(ok, details)

    def test_reference_form(self):
        self.check(CORRECT)

    def test_times_minus_one(self):
        self.check(f"-({CORRECT})")

    def test_other_spellings(self):
        self.check("Derivative(theta(xi, s), s) - Derivative(theta(xi, s), xi, xi) - L**2*Q0/(DT*k0)")


class PlantedErrorsFail(unittest.TestCase):
    """Each answer below has one deliberate mistake; the grader must reject it."""

    def check(self, text, why):
        ok, details = g.grade(text)
        self.assertFalse(ok, f"{why} was accepted: {details}")

    def test_diffusion_sign_flipped(self):
        self.check(f"{TS} + {TXX} - {GROUP}", "a flipped diffusion sign")

    def test_source_missing(self):
        self.check(f"{TS} - {TXX}", "a missing source")

    def test_wrong_group_length_power(self):
        self.check(f"{TS} - {TXX} - Q0*L/(k0*DT)", "L in place of L**2")

    def test_group_upside_down(self):
        self.check(f"{TS} - {TXX} - k0*DT/(Q0*L**2)", "the inverted group")

    def test_dimensional_prefactor_left(self):
        self.check(f"k0*DT/L**2*({CORRECT})", "an equation still carrying k0 DT / L^2")

    def test_rho_c_left_in_group(self):
        self.check(f"{TS} - {TXX} - Q0*L**2/(rho*c*DT)", "rho c in place of k0")

    def test_dimensional_names_refused(self):
        self.check("Derivative(T(x, t), t)", "dimensional names")

    def test_code_refused(self):
        self.check("__import__('os').system('true')", "code")


if __name__ == "__main__":
    unittest.main()
