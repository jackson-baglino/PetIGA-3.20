"""Tests of the T03 grader: a correct answer passes in several forms; each planted error fails.

Run: venv_workflows/bin/python -m unittest discover -s workflows/model-development/agents/beth-mathematician/tests
"""
import importlib.util
import unittest
from pathlib import Path

GRADER = Path(__file__).resolve().parents[1] / "tasks" / "T03_allen_cahn" / "grade.py"
spec = importlib.util.spec_from_file_location("grade_T03", GRADER)
g = importlib.util.module_from_spec(spec)
spec.loader.exec_module(g)

P = "phi(x, y, t)"
PT = f"Derivative({P}, t)"
LAP = f"(Derivative({P}, (x, 2)) + Derivative({P}, (y, 2)))"
WELL = f"2*W*{P}*(1 - {P})*(1 - 2*{P})"
CORRECT = f"{PT} - M*kappa*{LAP} + M*{WELL}"


class Passes(unittest.TestCase):
    def check(self, text):
        ok, details = g.grade(text)
        self.assertTrue(ok, details)

    def test_reference_form(self):
        self.check(CORRECT)

    def test_times_minus_one(self):
        self.check(f"-({CORRECT})")

    def test_expanded_well(self):
        self.check(f"{PT} - M*kappa*{LAP} + M*W*(4*{P}**3 - 6*{P}**2 + 2*{P})")

    def test_reference_matches_euler_lagrange(self):
        self.assertEqual(g.sp.simplify(g.mu_reference() - g.mu_from_euler_lagrange()), 0)


class PlantedErrorsFail(unittest.TestCase):
    """Each answer below has one deliberate mistake; the grader must reject it."""

    def check(self, text, why):
        ok, details = g.grade(text)
        self.assertFalse(ok, f"{why} was accepted: {details}")

    def test_laplacian_sign_flipped(self):
        self.check(f"{PT} + M*kappa*{LAP} + M*{WELL}", "a flipped Laplacian sign")

    def test_half_kept_on_kappa(self):
        self.check(f"{PT} - M*kappa/2*{LAP} + M*{WELL}", "kappa/2 in place of kappa")

    def test_factor_two_lost_in_well(self):
        self.check(f"{PT} - M*kappa*{LAP} + M*W*{P}*(1 - {P})*(1 - 2*{P})", "a lost factor 2")

    def test_wrong_well_derivative(self):
        self.check(f"{PT} - M*kappa*{LAP} + 2*M*W*{P}*(1 - {P})", "phi(1 - phi) in place of its derivative")

    def test_y_derivative_missing(self):
        self.check(f"{PT} - M*kappa*Derivative({P}, (x, 2)) + M*{WELL}", "a missing y-derivative")

    def test_uphill_flow(self):
        self.check(f"{PT} + M*kappa*{LAP} - M*{WELL}", "phi_t = +M mu")

    def test_mobility_missing(self):
        self.check(f"{PT} - kappa*{LAP} + {WELL}", "a missing mobility")

    def test_energy_density_instead_of_mu(self):
        self.check(f"{PT} + M*(kappa/2*(Derivative({P}, x)**2 + Derivative({P}, y)**2) + W*{P}**2*(1 - {P})**2)",
                   "the energy density in place of its variational derivative")


if __name__ == "__main__":
    unittest.main()
