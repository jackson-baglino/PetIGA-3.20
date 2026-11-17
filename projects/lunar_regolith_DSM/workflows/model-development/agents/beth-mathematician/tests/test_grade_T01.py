"""Tests of the T01 grader: a correct answer passes in several forms; each planted error fails.

Run: venv_workflows/bin/python -m unittest discover -s workflows/model-development/agents/beth-mathematician/tests
"""
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

GRADER = Path(__file__).resolve().parents[1] / "tasks" / "T01_heat_weak_form" / "grade.py"
spec = importlib.util.spec_from_file_location("grade_T01", GRADER)
g = importlib.util.module_from_spec(spec)
spec.loader.exec_module(g)

TX = "Derivative(T(x, t), x)"
TT = "Derivative(T(x, t), t)"
WX = "Derivative(w(x), x)"
BND = f"Subs(w(x)*k(x)*{TX}, x, L) - Subs(w(x)*k(x)*{TX}, x, 0)"

CORRECT = (f"Integral(rho*c*w(x)*{TT}, (x, 0, L)) + Integral({WX}*k(x)*{TX}, (x, 0, L))"
           f" - Integral(w(x)*Q(x), (x, 0, L)) - ({BND})")


class Passes(unittest.TestCase):
    def check(self, text):
        ok, details = g.grade(text)
        self.assertTrue(ok, details)

    def test_reference_form(self):
        self.check(CORRECT)

    def test_whole_equation_times_minus_one(self):
        self.check(f"-({CORRECT})")

    def test_one_integral_and_point_values(self):
        self.check(f"Integral(rho*c*w(x)*{TT} + {WX}*k(x)*{TX} - w(x)*Q(x), (x, 0, L))"
                   f" - w(L)*k(L)*Subs({TX}, x, L) + w(0)*k(0)*Subs({TX}, x, 0)")


class PlantedErrorsFail(unittest.TestCase):
    """Each answer below has one deliberate mistake; the grader must reject it."""

    def check(self, text, why):
        ok, details = g.grade(text)
        self.assertFalse(ok, f"{why} was accepted: {details}")

    def test_flux_sign_flipped(self):
        self.check(CORRECT.replace(f"+ Integral({WX}", f"- Integral({WX}"), "flipped conduction sign")

    def test_boundary_term_dropped(self):
        self.check(CORRECT.replace(f" - ({BND})", ""), "missing boundary term")

    def test_conductivity_dropped(self):
        self.check(CORRECT.replace(f"Integral({WX}*k(x)*", f"Integral({WX}*"), "missing k in the flux")

    def test_source_sign_flipped(self):
        self.check(CORRECT.replace("- Integral(w(x)*Q(x)", "+ Integral(w(x)*Q(x)"), "flipped source sign")

    def test_not_integrated_by_parts(self):
        # Equal in value to the weak form, but the second derivative is still there.
        self.check(f"Integral(w(x)*(rho*c*{TT} - Derivative(k(x)*{TX}, x) - Q(x)), (x, 0, L))",
                   "strong form under an integral")

    def test_blocked_or_empty(self):
        self.check("", "empty answer")


class UntrustedTextRejected(unittest.TestCase):
    def check(self, text):
        ok, details = g.grade(text)
        self.assertFalse(ok)
        return details

    def test_import(self):
        self.assertIn("characters outside", self.check("__import__('os').system('true')"))

    def test_import_without_quotes(self):
        self.assertIn("not allowed", self.check("__import__(os)"))

    def test_attribute_access(self):
        self.assertIn("'.'", self.check("x.func"))

    def test_unknown_function(self):
        self.assertIn("not allowed", self.check("sin(x)"))


class CommandLine(unittest.TestCase):
    def run_grader(self, result):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "result.json"
            p.write_text(json.dumps(result))
            r = subprocess.run([sys.executable, str(GRADER), str(p)], capture_output=True, text=True)
        return r.returncode, json.loads(r.stdout)

    def test_pass_exit_0(self):
        code, out = self.run_grader({"status": "done", "result_expr": CORRECT, "steps": [],
                                     "assumptions_used": [], "questions": []})
        self.assertEqual((code, out["pass"]), (0, True))

    def test_blocked_exit_1(self):
        code, out = self.run_grader({"status": "blocked", "result_expr": "", "steps": [],
                                     "assumptions_used": [], "questions": ["?"]})
        self.assertEqual((code, out["pass"]), (1, False))


if __name__ == "__main__":
    unittest.main()
