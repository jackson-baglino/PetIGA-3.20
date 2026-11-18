"""Tests of the T06 grader: right free energies pass in other forms; each planted error fails a check that
needs no key.

Run: venv_workflows/bin/python -m unittest discover -s workflows/model-development/agents/beth-mathematician/tests
"""
import importlib.util
import unittest
from pathlib import Path

GRADER = Path(__file__).resolve().parents[1] / "tasks" / "T06_free_energy" / "grade.py"
spec = importlib.util.spec_from_file_location("grade_T06", GRADER)
g = importlib.util.module_from_spec(spec)
spec.loader.exec_module(g)

P = "phi(x, y, z, t)"
G2 = " + ".join(f"Derivative({P}, {c})**2" for c in "xyz")


def failed(expr, symbols=()):
    ok, details = g.grade({"status": "done", "result_expr": expr, "symbols": list(symbols), "neglected": []})
    return ok, details if isinstance(details, str) else [c["check"] for c in details["checks"] if not c["pass"]]


class Passes(unittest.TestCase):
    def test_key_check(self):
        ok, why = g.key_check()
        self.assertTrue(ok, why)

    def test_kappa_and_W_written_out(self):
        """kappa = 6 gamma eps, W = 3 gamma / eps, written as kappa/2 |grad|^2 + W phi^2 (1 - phi)^2."""
        self.assertEqual(failed(f"(6*gamma*eps)/2*({G2}) + (3*gamma/eps)*{P}**2*(1 - {P})**2"), (True, []))

    def test_obstacle_free_quartic_with_another_width_convention(self):
        """The same energy with eps defined as twice the decay length: still energy gamma."""
        self.assertEqual(failed(f"3*gamma*(eps/2)*({G2}) + 3*gamma*{P}**2*(1 - {P})**2/(eps/2)"), (True, []))


class PlantedErrorsFail(unittest.TestCase):
    def test_unequal_wells(self):
        self.assertIn("wells", failed(f"3*gamma*eps*({G2}) + 3*gamma*{P}**2*(1 - {P})**2/eps + gamma*{P}/eps")[1])

    def test_interface_energy_off_by_a_factor(self):
        self.assertIn("interface", failed(f"gamma*eps*({G2}) + gamma*{P}**2*(1 - {P})**2/eps")[1])

    def test_second_derivatives_refused(self):
        ok, why = failed(f"-3*gamma*eps*{P}*Derivative({P}, (x, 2)) + 3*gamma*{P}**2*(1 - {P})**2/eps")
        self.assertFalse(ok)
        self.assertIn("form", why)

    def test_declared_function_refused(self):
        syms = [{"name": "a(phi)", "meaning": "a coefficient", "units": "J m^-3", "typical": 1}]
        ok, why = failed(f"3*gamma*eps*({G2}) + a({P})*{P}**2*(1 - {P})**2", syms)
        self.assertFalse(ok)


if __name__ == "__main__":
    unittest.main()
